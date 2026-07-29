//
//  Grid3Drcfs.h
//  ttcr
//
//  Created by Bernard Giroux on 16-01-07.
//  Copyright © 2016 Bernard Giroux. All rights reserved.
//

/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/**
 * @file Grid3Drcfs.h
 * @brief Fast sweeping solver on a 3-D rectilinear grid, from a cell slowness
 *        model.
 *
 * Declares ttcr::Grid3Drcfs, the 3-D counterpart of ttcr::Grid2Drcfs. As there,
 * the class accepts a cell slowness model but derives from the node-slowness
 * base and averages the values onto the nodes; see @ref g3drcfs_hybrid.
 *
 * @warning This solver supports **cubic cells only**, and that is not
 *          validated anywhere. @sa @ref g3drcfs_cubic
 *
 * @sa Grid3Drn.h, Grid3Drc.h, Grid2Drcfs.h, Grid3Drcfs_OpenCL.h
 */

#ifndef ttcr_Grid3Drcfs_h
#define ttcr_Grid3Drcfs_h

#include <cmath>
#include <stdexcept>

#include "Grid3Drn.h"
#include "Node3Dn.h"

// TODO: change base class for Grid3Drc (?) and use tt_from_rp to improve accuracy
namespace ttcr {

    /**
     * @brief Fast sweeping eikonal solver taking a cell slowness model.
     *
     * @tparam T1 floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2 integer type of node and cell indices.
     *
     * @section g3drcfs_hybrid Why it derives from Grid3Drn, not Grid3Drc
     * As with ttcr::Grid2Drcfs, the @c rc in the name describes the model the
     * caller supplies, not the internal representation. Fast sweeping needs a
     * slowness value at every node, so @ref setSlowness accepts one value per
     * cell and averages them onto the nodes — which is why the base is
     * ttcr::Grid3Drn. A side effect is that sharp contrasts are smoothed over
     * one cell before the solve begins.
     *
     * @section g3drcfs_cubic Cubic cells only
     * The constructor takes a **single** cell size @p ddx and passes it as all
     * three spacings, so this solver only represents cubic cells.
     *
     * @warning Nothing checks that. `grids.h` reads three independent spacings
     *          from the model file but hands only @c d[0] to this constructor,
     *          so a model with @f$d_x \neq d_y@f$ or @f$d_x \neq d_z@f$ is
     *          silently solved on a grid whose y and z spacings have been
     *          replaced by @f$d_x@f$ — wrong traveltimes, no diagnostic. The
     *          shortest-path and dynamic shortest-path builders pass all three
     *          spacings and are unaffected.
     *
     * @note Not templated on a @c CELL policy, so isotropic only.
     *
     * @sa Grid3Drn.h, Grid3Drc.h, Grid2Drcfs.h, Grid3Drcfs_OpenCL.h
     */
    template<typename T1, typename T2>
    class Grid3Drcfs : public Grid3Drn<T1,T2,Node3Dn<T1,T2>> {
    public:
        /**
         * @brief Build the grid and its nodes.
         *
         * @param nx    number of cells along x.
         * @param ny    number of cells along y.
         * @param nz    number of cells along z.
         * @param ddx   cell size, used for **all three** axes.
         *              @sa @ref g3drcfs_cubic
         * @param minx  x coordinate of the grid origin.
         * @param miny  y coordinate of the grid origin.
         * @param minz  z coordinate of the grid origin.
         * @param eps   convergence tolerance, as a mean per-node traveltime
         *              change; scaled internally by the node count.
         * @param maxit maximum number of sweep iterations.
         * @param w     use the 3rd-order WENO stencil.
         * @param ttrp  recompute receiver traveltimes along the raypath.
         * @param intVel interpolate velocity rather than slowness
         *              (ttcr::input_parameters::processVel).
         * @param nt    number of threads.
         * @param _translateOrigin shift the grid origin to (0,0,0).
         *
         * @post Nodes and neighbour lists are built. Slowness is **not** set —
         *       call @ref setSlowness before raytracing.
         */
        Grid3Drcfs(const T2 nx, const T2 ny, const T2 nz, const T1 ddx,
                   const T1 minx, const T1 miny, const T1 minz,
                   const T1 eps, const int maxit, const bool w,
                   const bool ttrp=true, const bool intVel=false,
                   const size_t nt=1,
                   const bool _translateOrigin=false) :
        Grid3Drn<T1,T2,Node3Dn<T1,T2>>(nx, ny, nz, ddx, ddx, ddx, minx, miny, minz, ttrp, intVel, nt, _translateOrigin),
        epsilon(eps), nitermax(maxit), niter_final(0), niterw_final(0), weno3(w)
        {
            this->buildGridNodes();
            this->template buildGridNeighbors<Node3Dn<T1,T2>>(this->nodes);
            epsilon *= static_cast<T1>(this->nodes.size());  // per-node tol -> L1-sum threshold (nodes built)
        }

        /// Destructor.
        virtual ~Grid3Drcfs() {
        }

        /**
         * @brief Set the cell slowness model and average it onto the nodes.
         * @param s one slowness per cell, @f$n_{cx}n_{cy}n_{cz}@f$ values in the
         *          x-fastest order of @ref g3drc_numbering.
         * @throws std::length_error if @p s has the wrong size.
         * @note A corner node takes its single adjacent cell, an edge node the
         *       mean of two, a face node of four and an interior node of eight.
         *       @sa @ref g3drcfs_hybrid
         * @note Unlike ttcr::Grid2Drcfs this class keeps no copy of the cell
         *       values, so there is no @c getCellSlowness here.
         */
        void setSlowness(const std::vector<T1>& s);

        /// @return Number of sweep iterations the last solve took.
        const int get_niter() const { return niter_final; }
        /// @return Number of WENO sweep iterations the last solve took.
        const int get_niterw() const { return niterw_final; }

    protected:
        /// Convergence threshold, holding the **scaled** value: the constructor
        /// multiplies the supplied tolerance by the node count, so the loop can
        /// test an L1 sum while the user supplies a mean per-node change.
        T1 epsilon;
        int nitermax;             ///< Iteration cap for the sweeps.
        mutable int niter_final;  ///< Iterations used by the last solve; @c mutable so the @c const raytrace can record it.
        mutable int niterw_final; ///< WENO iterations used by the last solve.
        bool weno3;               ///< Use the 3rd-order WENO stencil.

    private:
        /// @name Non-copyable
        /// @{
        Grid3Drcfs() {}
        Grid3Drcfs(const Grid3Drcfs<T1,T2>& g) {}
        Grid3Drcfs<T1,T2>& operator=(const Grid3Drcfs<T1,T2>& g) = delete;
        /// @}

        /**
         * @brief Propagate the traveltime field and evaluate it at the receivers.
         * @param Tx       source positions.
         * @param t0       origin time of each source.
         * @param Rx       receiver positions.
         * @param threadNo thread to compute on.
         * @note Sweeps until the change falls below @ref epsilon or
         *       @ref nitermax is reached, recording the count in
         *       @ref niter_final.
         */
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<sxyz<T1>>& Rx,
                      const size_t threadNo=0) const;

        /**
         * @brief Propagate once and evaluate at several receiver sets.
         * @param Tx       source positions.
         * @param t0       origin time of each source.
         * @param Rx       the receiver sets.
         * @param threadNo thread to compute on.
         */
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      const size_t threadNo=0) const;

    };


    template<typename T1, typename T2>
    void Grid3Drcfs<T1,T2>::setSlowness(const std::vector<T1>& s) {

        if ( static_cast<size_t>(this->ncx)*this->ncy*this->ncz != s.size() ) {
            throw std::length_error("Error: slowness vectors of incompatible size.");
        }

        // interpolate slowness at grid nodes

        const size_t nx = this->ncx;
        const size_t ny = this->ncy;
        const size_t nz = this->ncz;

        // corners
        this->nodes[                       0].setNodeSlowness( s[                       0] );
        this->nodes[                      nx].setNodeSlowness( s[                    nx-1] );
        this->nodes[           ny *(nx+1)   ].setNodeSlowness( s[          (ny-1)*nx     ] );
        this->nodes[           ny *(nx+1)+nx].setNodeSlowness( s[          (ny-1)*nx+nx-1] );
        this->nodes[(nz*(ny+1)   )*(nx+1)   ].setNodeSlowness( s[((nz-1)*ny     )*nx     ] );
        this->nodes[(nz*(ny+1)   )*(nx+1)+nx].setNodeSlowness( s[((nz-1)*ny     )*nx+nx-1] );
        this->nodes[(nz*(ny+1)+ny)*(nx+1)   ].setNodeSlowness( s[((nz-1)*ny+ny-1)*nx     ] );
        this->nodes[(nz*(ny+1)+ny)*(nx+1)+nx].setNodeSlowness( s[((nz-1)*ny+ny-1)*nx+nx-1] );

        // edges
        for ( size_t i=1; i<nx; ++i ) {
            this->nodes[                      i].setNodeSlowness( 0.5*(s[                    i]+s[                    i-1]) );
            this->nodes[           ny *(nx+1)+i].setNodeSlowness( 0.5*(s[          (ny-1)*nx+i]+s[          (ny-1)*nx+i-1]) );
            this->nodes[(nz*(ny+1)   )*(nx+1)+i].setNodeSlowness( 0.5*(s[((nz-1)*ny     )*nx+i]+s[((nz-1)*ny     )*nx+i-1]) );
            this->nodes[(nz*(ny+1)+ny)*(nx+1)+i].setNodeSlowness( 0.5*(s[((nz-1)*ny+ny-1)*nx+i]+s[((nz-1)*ny+ny-1)*nx+i-1]) );
        }
        for ( size_t j=1; j<ny; ++j ) {
            this->nodes[           j *(nx+1)   ].setNodeSlowness( 0.5*(s[            j*nx     ]+s[          (j-1)*nx     ]) );
            this->nodes[           j *(nx+1)+nx].setNodeSlowness( 0.5*(s[            j*nx+nx-1]+s[          (j-1)*nx+nx-1]) );
            this->nodes[(nz*(ny+1)+j)*(nx+1)   ].setNodeSlowness( 0.5*(s[((nz-1)*ny+j)*nx     ]+s[((nz-1)*ny+j-1)*nx     ]) );
            this->nodes[(nz*(ny+1)+j)*(nx+1)+nx].setNodeSlowness( 0.5*(s[((nz-1)*ny+j)*nx+nx-1]+s[((nz-1)*ny+j-1)*nx+nx-1]) );
        }
        for ( size_t k=1; k<nz; ++k ) {
            this->nodes[(k*(ny+1)   )*(nx+1)   ].setNodeSlowness( 0.5*(s[(k*ny     )*nx     ]+s[((k-1)*ny     )*nx     ]) );
            this->nodes[(k*(ny+1)   )*(nx+1)+nx].setNodeSlowness( 0.5*(s[(k*ny     )*nx+nx-1]+s[((k-1)*ny     )*nx+nx-1]) );
            this->nodes[(k*(ny+1)+ny)*(nx+1)   ].setNodeSlowness( 0.5*(s[(k*ny+ny-1)*nx     ]+s[((k-1)*ny+ny-1)*nx     ]) );
            this->nodes[(k*(ny+1)+ny)*(nx+1)+nx].setNodeSlowness( 0.5*(s[(k*ny+ny-1)*nx+nx-1]+s[((k-1)*ny+ny-1)*nx+nx-1]) );
        }

        // faces
        for ( size_t i=1; i<nx; ++i ) {
            for ( size_t j=1; j<ny; ++j ) {
                this->nodes[           j *(nx+1)+i].setNodeSlowness( 0.25*(s[             j *nx+i]+s[             j *nx+i-1]+
                                                                           s[          (j-1)*nx+i]+s[          (j-1)*nx+i-1]) );
                this->nodes[(nz*(ny+1)+j)*(nx+1)+i].setNodeSlowness( 0.25*(s[((nz-1)*ny+  j)*nx+i]+s[((nz-1)*ny+  j)*nx+i-1]+
                                                                           s[((nz-1)*ny+j-1)*nx+i]+s[((nz-1)*ny+j-1)*nx+i-1]) );
            }
        }
        for ( size_t i=1; i<nx; ++i ) {
            for ( size_t k=1; k<nz; ++k ) {
                this->nodes[(k*(ny+1)   )*(nx+1)+i].setNodeSlowness( 0.25*(s[(   k *ny     )*nx+i]+s[(   k *ny     )*nx+i-1]+
                                                                           s[((k-1)*ny     )*nx+i]+s[((k-1)*ny     )*nx+i-1]) );
                this->nodes[(k*(ny+1)+ny)*(nx+1)+i].setNodeSlowness( 0.25*(s[(   k *ny+ny-1)*nx+i]+s[(   k *ny+ny-1)*nx+i-1]+
                                                                           s[((k-1)*ny+ny-1)*nx+i]+s[((k-1)*ny+ny-1)*nx+i-1]) );
            }
        }
        for ( size_t j=1; j<ny; ++j ) {
            for ( size_t k=1; k<nz; ++k ) {
                this->nodes[(k*(ny+1)+j)*(nx+1)   ].setNodeSlowness( 0.25*(s[(k*ny+  j)*nx     ]+s[((k-1)*ny+  j)*nx     ]+
                                                                           s[(k*ny+j-1)*nx     ]+s[((k-1)*ny+j-1)*nx     ]) );
                this->nodes[(k*(ny+1)+j)*(nx+1)+nx].setNodeSlowness( 0.25*(s[(k*ny+  j)*nx+nx-1]+s[((k-1)*ny+  j)*nx+nx-1]+
                                                                           s[(k*ny+j-1)*nx+nx-1]+s[((k-1)*ny+j-1)*nx+nx-1]) );
            }
        }

        // interior
        for ( size_t i=1; i<nx; ++i ) {
            for ( size_t j=1; j<ny; ++j ) {
                for ( size_t k=1; k<nz; ++k ) {
                    this->nodes[(k*(ny+1)+j)*(nx+1)+i].setNodeSlowness( 0.125*(s[(    k*ny+j  )*nx+i  ]+
                                                                               s[(    k*ny+j  )*nx+i-1]+
                                                                               s[(    k*ny+j-1)*nx+i  ]+
                                                                               s[(    k*ny+j-1)*nx+i-1]+
                                                                               s[((k-1)*ny+j  )*nx+i  ]+
                                                                               s[((k-1)*ny+j  )*nx+i-1]+
                                                                               s[((k-1)*ny+j-1)*nx+i  ]+
                                                                               s[((k-1)*ny+j-1)*nx+i-1]) );
                }
            }
        }
    }


    template<typename T1, typename T2>
    void Grid3Drcfs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<sxyz<T1>>& Rx,
                                     const size_t threadNo) const {
        if ( verbose > 2 ) {
            std::cout << "\nIn Grid3Drcfs::raytrace(Tx, t0, Rx, threadNo)\n" << std::endl;
        }
        this->checkPts(Tx, true);
        this->checkPts(Rx, true);

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        // Set Tx pts
        std::vector<bool> frozen( this->nodes.size(), false );
        int npts = weno3 ? 2 : 1;
        this->initFSM(Tx, t0, frozen, npts, threadNo);

        std::vector<T1> times( this->nodes.size() );
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            times[n] = this->nodes[n].getTT( threadNo );
        }

        T1 change = std::numeric_limits<T1>::max();
        if ( weno3 == true ) {
            int niter = 0;
            int niterw = 0;
            if ( this->dx != this->dz || this->dx != this->dy ) {
                throw std::logic_error("Error: WENO stencil needs dx equal to dz");
            }
            while ( change >= epsilon && niter<nitermax ) {
                this->sweep(frozen, threadNo);
                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = std::abs( times[n] - this->nodes[n].getTT(threadNo) );

                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                niter++;
            }
            change = std::numeric_limits<T1>::max();
            while ( change >= epsilon && niterw<nitermax ) {
                this->sweep_weno3(frozen, threadNo);
                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = std::abs( times[n] - this->nodes[n].getTT(threadNo) );

                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                niterw++;
            }
            niter_final = niter;
            niterw_final = niterw;
        } else {
            int niter = 0;
            while ( change >= epsilon && niter<nitermax ) {
                this->sweep(frozen, threadNo);

                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = std::abs( times[n] - this->nodes[n].getTT(threadNo) );

                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                niter++;
            }
            niter_final = niter;
        }
    }

    template<typename T1, typename T2>
    void Grid3Drcfs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<std::vector<sxyz<T1>>>& Rx,
                                     const size_t threadNo) const {

        this->checkPts(Tx, true);
        for ( size_t n=0; n<Rx.size(); ++n ) {
            this->checkPts(Rx[n], true);
        }

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        // Set Tx pts
        std::vector<bool> frozen( this->nodes.size(), false );
        int npts = 1;
        if ( weno3 == true ) npts = 2;
        this->initFSM(Tx, t0, frozen, npts, threadNo);

        std::vector<T1> times( this->nodes.size() );
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            times[n] = this->nodes[n].getTT( threadNo );
        }

        T1 change = std::numeric_limits<T1>::max();
        if ( weno3 == true ) {
            int niter = 0;
            int niterw = 0;
            if ( this->dx != this->dz || this->dx != this->dy ) {
                throw std::logic_error("Error: WENO stencil needs dx equal to dz");
            }
            while ( change >= epsilon && niter<nitermax ) {
                this->sweep(frozen, threadNo);
                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = std::abs( times[n] - this->nodes[n].getTT(threadNo) );

                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                niter++;
            }
            change = std::numeric_limits<T1>::max();
            while ( change >= epsilon && niterw<nitermax ) {
                this->sweep_weno3(frozen, threadNo);
                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = std::abs( times[n] - this->nodes[n].getTT(threadNo) );

                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                niterw++;
            }
            niter_final = niter;
            niterw_final = niterw;
        } else {
            int niter = 0;
            while ( change >= epsilon && niter<nitermax ) {
                this->sweep(frozen, threadNo);

                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = std::abs( times[n] - this->nodes[n].getTT(threadNo) );

                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                niter++;
            }
            niter_final = niter;
        }
    }

}

#endif /* Grid3Drcfs_h */

//
//  Grid3Drnfs.h
//  ttcr
//
//  Created by Bernard Giroux on 15-12-27.ncz
//  Copyright © 2015 Bernard Giroux. All rights reserved.
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
 * @file Grid3Drnfs.h
 * @brief Fast sweeping solver on a 3-D rectilinear grid with node-based slowness.
 *
 * Declares ttcr::Grid3Drnfs. Slowness is taken at the nodes directly, unlike
 * ttcr::Grid3Drcfs which accepts a cell model and averages it onto them.
 *
 * @warning Cubic cells only, unvalidated — see the class warning.
 *
 * @sa Grid3Drn.h, Grid3Drcfs.h, Grid2Drnfs.h, Grid3Drnfs_OpenCL.h
 */

#ifndef ttcr_Grid3Drnfs_h
#define ttcr_Grid3Drnfs_h

#include <cmath>
#include <utility>

#include "Grid3Drn.h"
#include "Node3Dn.h"

namespace ttcr {

    template<typename T1, typename T2>
    /**
     * @brief Fast sweeping eikonal solver with node-based slowness.
     *
     * @tparam T1 floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2 integer type of node and cell indices.
     *
     * The 3-D counterpart of ttcr::Grid2Drnfs and the node-slowness sibling of
     * ttcr::Grid3Drcfs. Slowness is taken at the nodes directly, with none of
     * the cell-to-node averaging ttcr::Grid3Drcfs performs. The sweep and
     * single-node update stencils live in ttcr::Grid3Drn
     * (@ref g3drn_sweeps); this class supplies the iteration driver and the
     * convergence test.
     *
     * @warning **Cubic cells only.** The constructor takes a single cell size
     *          @p ddx and passes it as all three spacings, and nothing validates
     *          that: `grids.h` reads three independent spacings from the model
     *          file but hands only @c d[0] to this constructor. A model with
     *          @f$d_x \neq d_y@f$ or @f$d_x \neq d_z@f$ is silently solved on a
     *          grid whose y and z spacings have been replaced by @f$d_x@f$.
     *          All four 3-D fast sweeping classes share this restriction; the
     *          shortest-path and dynamic shortest-path builders pass all three
     *          spacings and are unaffected.
     *
     * @section g3drnfs_conv Convergence
     * Sweeping stops when the summed change in nodal traveltime falls below a
     * threshold. The constructor multiplies the supplied @p eps by the node
     * count, so the user-facing tolerance is a **mean** per-node change while
     * the loop tests an L1 sum; @p maxit caps the iterations.
     *
     * @sa Grid3Drn.h, Grid3Drcfs.h, Grid2Drnfs.h, Grid3Drnfs_OpenCL.h
     */
    class Grid3Drnfs : public Grid3Drn<T1,T2,Node3Dn<T1,T2>> {
    public:
        /**
         * @brief Build the grid and its nodes.
         *
         * @param nx    number of cells along x.
         * @param ny    number of cells along y.
         * @param nz    number of cells along z.
         * @param ddx   cell size, used for **all three** axes — see the class
         *              warning.
         * @param minx  x coordinate of the grid origin.
         * @param miny  y coordinate of the grid origin.
         * @param minz  z coordinate of the grid origin.
         * @param eps   convergence tolerance, as a mean per-node traveltime
         *              change; scaled internally by the node count.
         *              @sa @ref g3drnfs_conv
         * @param maxit maximum number of sweep iterations.
         * @param w     use the 3rd-order WENO stencil.
         * @param ttrp  recompute receiver traveltimes along the raypath.
         * @param intVel interpolate velocity rather than slowness.
         *              @sa @ref g3drn_procvel
         * @param nt    number of threads.
         * @param _translateOrigin shift the grid origin to (0,0,0).
         *
         * @post Nodes and neighbour lists are built. Slowness is **not** set.
         */
        Grid3Drnfs(const T2 nx, const T2 ny, const T2 nz, const T1 ddx,
                   const T1 minx, const T1 miny, const T1 minz,
                   const T1 eps, const int maxit, const bool w,
                   const bool ttrp=true, const bool intVel=false,
                   const size_t nt=1, const bool _translateOrigin=false) :
        Grid3Drn<T1,T2,Node3Dn<T1,T2>>(nx, ny, nz, ddx, ddx, ddx, minx, miny, minz, ttrp, intVel, nt, _translateOrigin),
        epsilon(eps), nitermax(maxit), niter_final(0), niterw_final(0), weno3(w)
        {
            this->buildGridNodes();
            this->template buildGridNeighbors<Node3Dn<T1,T2>>(this->nodes);
            epsilon *= static_cast<T1>(this->nodes.size());  // per-node tol -> L1-sum threshold (nodes built)
        }

        /// Destructor.
        ~Grid3Drnfs() {

        }

        /// @return Number of sweep iterations the last solve took.
        const int get_niter() const { return niter_final; }
        /// @return Number of WENO sweep iterations the last solve took.
        const int get_niterw() const { return niterw_final; }

    protected:
        /// Convergence threshold, holding the **scaled** value: the constructor
        /// multiplies the supplied tolerance by the node count.
        /// @sa @ref g3drnfs_conv
        T1 epsilon;
        int nitermax;             ///< Iteration cap for the sweeps.
        mutable int niter_final;  ///< Iterations used by the last solve; @c mutable so the @c const raytrace can record it.
        mutable int niterw_final; ///< WENO iterations used by the last solve.
        bool weno3;               ///< Use the 3rd-order WENO stencil.

    private:
        /// @name Non-copyable
        /// @{
        Grid3Drnfs() {}
        Grid3Drnfs(const Grid3Drnfs<T1,T2>& g) {}
        Grid3Drnfs<T1,T2>& operator=(const Grid3Drnfs<T1,T2>& g) = delete;
        /// @}

        /**
         * @brief Propagate the traveltime field and evaluate it at the receivers.
         * @param Tx       source positions.
         * @param t0       origin time of each source.
         * @param Rx       receiver positions.
         * @param threadNo thread to compute on.
         * @note Freezes a halo around each source via
         *       ttcr::Grid3Drn::initFSM, then sweeps until the change falls
         *       below @ref epsilon or @ref nitermax is reached.
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
    void Grid3Drnfs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<sxyz<T1>>& Rx,
                                     const size_t threadNo) const {

        this->checkPts(Tx, true);
        this->checkPts(Rx, true);

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

    template<typename T1, typename T2>
    void Grid3Drnfs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
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

#endif /* Grid3Drnfs_h */

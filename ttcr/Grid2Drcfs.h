//
//  Grid2Drcfs.h
//  ttcr
//
//  Created by Bernard Giroux on 15-12-23.
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
 * @file Grid2Drcfs.h
 * @brief Fast sweeping solver on a 2-D rectilinear grid, from a cell slowness model.
 *
 * Declares ttcr::Grid2Drcfs, which solves the eikonal equation by Gauss-Seidel
 * sweeps in alternating directions rather than by graph search.
 *
 * @section g2drcfs_hybrid Why it derives from Grid2Drn, not Grid2Drc
 * Despite the @c rc in its name this class derives from ttcr::Grid2Drn — the
 * **node**-slowness base — not from ttcr::Grid2Drc. The reason is that the
 * fast sweeping update needs a slowness value at each node, whereas the @c rc
 * model supplies one per cell. The class bridges the two: it accepts a cell
 * slowness vector, keeps a copy, and averages it onto the nodes.
 *
 * ttcr::Grid2Drcfs::setSlowness does that averaging — a corner node takes the
 * value of its single adjacent cell, an edge node the mean of the two cells
 * sharing it, and an interior node the mean of the four. So the @c rc in the
 * name describes the **model the caller supplies**, not the representation used
 * internally.
 *
 * A consequence worth knowing: the model is smoothed by that averaging, so a
 * sharp slowness contrast is blurred over one cell before the solve begins.
 * ttcr::Grid2Drcfs::getSlowness returns the original cell values, not the
 * averaged nodal ones.
 *
 * @section g2drcfs_conv Convergence
 * The sweeps stop when the mean change in nodal traveltime falls below
 * @p eps times the **range** of the traveltime field, so the tolerance
 * (ttcr::input_parameters::epsilon) is dimensionless and the same value
 * means the same thing whatever units the model is expressed in. The loop
 * itself tests an L1 sum, and ttcr::fsmTolerance folds in the node count
 * and the range each sweep. @p maxit caps the iterations.
 *
 * The WENO3 update is non-monotone early on: the change per sweep rises
 * before it falls, from a first value that shrinks as @f$O(h^2)@f$ under
 * refinement. A loop that simply tested the tolerance would therefore stop
 * on the way up on a fine enough grid, so a pass never terminates while
 * its change is still at or above the previous sweep's.
 *
 * @sa Grid2Drn.h, Grid2Drc.h, Grid2Drcsp.h, Grid2Drcfs_OpenCL.h
 */

#ifndef ttcr_Grid2Drcfs_h
#define ttcr_Grid2Drcfs_h

#include <cmath>
#include <stdexcept>

#include "Grid2Drn.h"
#include "Node2Dn.h"

namespace ttcr {

    /**
     * @brief Fast sweeping eikonal solver taking a cell slowness model.
     *
     * @tparam T1 floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2 integer type of node and cell indices.
     * @tparam S  point type, @ref sxz or @ref sxyz.
     *
     * @note Not templated on a @c CELL policy — unlike ttcr::Grid2Drcsp, this
     *       class is isotropic only. Anisotropy would have to survive the
     *       cell-to-node averaging, which it does not.
     */
    template<typename T1, typename T2, typename S>
    class Grid2Drcfs : public Grid2Drn<T1,T2,S,Node2Dn<T1,T2>> {
    public:
        /**
         * @brief Build the grid and its nodes.
         *
         * @param nx    number of cells along x.
         * @param nz    number of cells along z.
         * @param ddx   cell size along x.
         * @param ddz   cell size along z.
         * @param minx  x coordinate of the grid origin.
         * @param minz  z coordinate of the grid origin.
         * @param eps   convergence tolerance, **relative**: the sweeps stop
         *              once the mean per-node change in traveltime falls
         *              below this fraction of the traveltime range.
         *              @sa @ref g2drcfs_conv
         * @param maxit maximum number of sweep iterations.
         * @param w     use the 3rd-order WENO stencil instead of the first-order
         *              one (ttcr::input_parameters::weno3).
         * @param rt    use rotated stencils as well as axis-aligned ones
         *              (ttcr::input_parameters::rotated_template).
         * @param ttrp  recompute receiver traveltimes along the raypath.
         * @param nt    number of threads.
         *
         * @post Nodes and their neighbour lists are built. Slowness is **not**
         *       set — call @ref setSlowness before raytracing, or
         *       @ref getCellSlowness will throw.
         */
        Grid2Drcfs(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
                   const T1 minx, const T1 minz, const T1 eps, const int maxit,
                   const bool w, const bool rt, const bool ttrp,
                   const size_t nt=1) :
        Grid2Drn<T1,T2,S,Node2Dn<T1,T2>>(nx,nz,ddx,ddz,minx,minz,ttrp,nt),
        epsilon(eps), nitermax(maxit), niter_final(0), niterw_final(0),
        weno3(w), rotated_template(rt), hasCellSlown(false), slowness()
        {
            buildGridNodes();
            this->template buildGridNeighbors<Node2Dn<T1,T2>>(this->nodes);
        }

        /// Destructor.
        virtual ~Grid2Drcfs() {
        }

        /**
         * @brief Set the cell slowness model and average it onto the nodes.
         * @param s one slowness per cell, @f$n_{cx}n_{cz}@f$ values in the
         *          column-wise order of @ref g2drc_numbering.
         * @throws std::length_error if @p s has the wrong size.
         * @post A copy is kept for @ref getSlowness and @ref getCellSlowness,
         *       and every node receives the mean of its adjacent cells.
         *       @sa @ref g2drcfs_hybrid
         */
        void setSlowness(const std::vector<T1>& s);
        /**
         * @brief Retrieve the cell slowness model.
         * @param[out] s the values passed to @ref setSlowness.
         * @note Returns the original **cell** values, not the averaged nodal
         *       ones the solver actually used.
         */
        void getSlowness(std::vector<T1>& s) const {
            s = slowness;
        }

        /// @return Number of sweep iterations the last solve took.
        const int get_niter() const { return niter_final; }
        /// @return Number of WENO sweep iterations the last solve took.
        const int get_niterw() const { return niterw_final; }
        /// @return True once @ref setSlowness has been called.
        const bool hasCellSlowness() const { return hasCellSlown; }
        /**
         * @brief Slowness of one cell.
         * @param cell_no cell number, per @ref g2drc_numbering.
         * @return That cell's slowness.
         * @throws std::runtime_error if no slowness model has been set.
         * @warning @p cell_no is not range-checked.
         */
        const T1 getCellSlowness(const size_t cell_no) const {
            if ( hasCellSlown ) {
                return slowness[cell_no];
            } else {
                throw std::runtime_error("slowness data not assigned");
            }
        }
    private:
        /// Convergence threshold. Holds the **scaled** value: the constructor
        /// multiplies the supplied tolerance by the node count.
        /// @sa @ref g2drcfs_conv
        T1 epsilon;
        int nitermax;             ///< Iteration cap for the sweeps.
        mutable int niter_final;  ///< Iterations used by the last solve; @c mutable so the @c const raytrace can record it.
        mutable int niterw_final; ///< WENO iterations used by the last solve.
        bool weno3;               ///< Use the 3rd-order WENO stencil.
        bool rotated_template;    ///< Use rotated stencils in addition to axis-aligned ones.
        bool hasCellSlown;        ///< Whether @ref slowness has been populated.
        std::vector<T1> slowness; ///< Copy of the cell slowness model as supplied.

        /// @name Non-copyable
        /// Copy assignment is deleted, so assigning one grid to another is a
        /// compile error. The default and copy constructors still use the older
        /// private-and-defined idiom and would leave members uninitialised if
        /// ever called from within the class.
        /// @{
        Grid2Drcfs() {}
        Grid2Drcfs(const Grid2Drcfs<T1,T2,S>& g) {}
        Grid2Drcfs<T1,T2,S>& operator=(const Grid2Drcfs<T1,T2,S>& g) = delete;
        /// @}

        /// @brief Create the grid nodes and set their positions.
        /// @note Only primary nodes; fast sweeping needs no secondary ones.
        void buildGridNodes();

        /**
         * @brief Propagate the traveltime field and evaluate it at the receivers.
         * @param Tx       source positions.
         * @param t0       origin time of each source.
         * @param Rx       receiver positions.
         * @param threadNo thread to compute on.
         * @note Private here: callers reach it through the ttcr::Grid2D
         *       interface. Sweeps until the change falls below @ref epsilon or
         *       @ref nitermax is reached, then records the count in
         *       @ref niter_final.
         */
        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<S>& Rx,
                      const size_t threadNo=0) const;

        /**
         * @brief Propagate once and evaluate at several receiver sets.
         * @param Tx       source positions.
         * @param t0       origin time of each source.
         * @param Rx       pointers to the receiver sets.
         * @param threadNo thread to compute on.
         */
        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<const std::vector<S>*>& Rx,
                      const size_t threadNo=0) const;

    };

    template<typename T1, typename T2, typename S>
    void Grid2Drcfs<T1,T2,S>::setSlowness(const std::vector<T1>& s) {

        if ( static_cast<size_t>(this->ncx)*this->ncz != s.size() ) {
            throw std::length_error("Error: slowness vectors of incompatible size.");
        }

        // keep a copy of slowness values for cells
        slowness = s;
        hasCellSlown = true;

        // interpolate slowness at grid nodes

        const size_t nx = this->ncx;
        const size_t nz = this->ncz;

        // four corners
        this->nodes[0].setNodeSlowness( s[0] );
        this->nodes[nz].setNodeSlowness( s[nz-1] );
        this->nodes[nx*(nz+1)].setNodeSlowness( s[nz*(nx-1)] );
        this->nodes[(nx+1)*(nz+1)-1].setNodeSlowness( s[nx*nz-1] );

        // sides
        for ( size_t j=1; j<nz; ++j ) {
            this->nodes[j].setNodeSlowness( 0.5*(s[j]+s[j-1]) );
            this->nodes[nx*(nz+1)+j].setNodeSlowness( 0.5*(s[nz*(nx-1)+j]+s[nz*(nx-1)+j-1]) );
        }
        // top & bottom
        for ( size_t i=1; i<nx; ++i ) {
            this->nodes[i*(nz+1)].setNodeSlowness( 0.5*(s[i*nz]+s[(i-1)*nz]) );
            this->nodes[i*(nz+1)+nz].setNodeSlowness( 0.5*(s[(i+1)*nz-1]+s[i*nz-1]) );
        }
        for ( size_t i=1; i<nx; ++i ) {
            for ( size_t j=1; j<nz; ++j ) {
                this->nodes[i*(nz+1)+j].setNodeSlowness(0.25*(s[i*nz+j]+
                                                              s[i*nz+j-1]+
                                                              s[(i-1)*nz+j]+
                                                              s[(i-1)*nz+j-1]));
            }
        }
    }


    template<typename T1, typename T2, typename S>
    void Grid2Drcfs<T1,T2,S>::buildGridNodes() {

        T2 cell_upLeft = std::numeric_limits<T2>::max();
        T2 cell_upRight = std::numeric_limits<T2>::max();
        T2 cell_downLeft = 0;
        T2 cell_downRight = 0;

        for ( T2 n=0, nc=0; nc<=this->ncx; ++nc ) {

            T1 x = this->xmin + nc*this->dx;

            for ( T2 nr=0; nr<=this->ncz; ++nr ) {

                T1 z = this->zmin + nr*this->dz;

                if ( nr < this->ncz && nc < this->ncx ) {
                    cell_downRight = nc*this->ncz + nr;
                }
                else {
                    cell_downRight = std::numeric_limits<T2>::max();
                }

                if ( nr > 0 && nc < this->ncx ) {
                    cell_upRight = nc*this->ncz + nr - 1;
                }
                else {
                    cell_upRight = std::numeric_limits<T2>::max();
                }

                if ( nr < this->ncz && nc > 0 ) {
                    cell_downLeft = (nc-1)*this->ncz + nr;
                }
                else {
                    cell_downLeft = std::numeric_limits<T2>::max();
                }

                if ( nr > 0 && nc > 0 ) {
                    cell_upLeft = (nc-1)*this->ncz + nr - 1;
                }
                else {
                    cell_upLeft = std::numeric_limits<T2>::max();
                }

                if ( cell_upLeft != std::numeric_limits<T2>::max() ) {
                    this->nodes[n].pushOwner( cell_upLeft );
                }
                if ( cell_downLeft != std::numeric_limits<T2>::max() ) {
                    this->nodes[n].pushOwner( cell_downLeft );
                }
                if ( cell_upRight != std::numeric_limits<T2>::max() ) {
                    this->nodes[n].pushOwner( cell_upRight );
                }
                if ( cell_downRight != std::numeric_limits<T2>::max() ) {
                    this->nodes[n].pushOwner( cell_downRight );
                }

                this->nodes[n].setX( x );
                this->nodes[n].setZ( z );
                this->nodes[n].setGridIndex( n );
                this->nodes[n].setPrimary(true);

                ++n;
            }
        }
    }


    template<typename T1, typename T2, typename S>
    void Grid2Drcfs<T1,T2,S>::raytrace(const std::vector<S>& Tx,
                                       const std::vector<T1>& t0,
                                       const std::vector<S>& Rx,
                                       const size_t threadNo) const {

        this->checkPts(Tx);
        this->checkPts(Rx);

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        // Set Tx pts
        std::vector<bool> frozen( this->nodes.size(), false );
        int npts = 1;
        if (weno3 == true) npts = 2;
        this->initFSM(Tx, t0, frozen, npts, threadNo);

        std::vector<T1> times( this->nodes.size() );
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            times[n] = this->nodes[n].getTT( threadNo );
        }

        T1 change = std::numeric_limits<T1>::max();
        T1 prev = 0.0;   // previous sweep's change, for the non-monotone guard
        T1 tref = 0.0;   // traveltime range, refreshed by fsmChange
        T1 tol = 0.0;    // absolute stop threshold, from epsilon and tref
        if ( weno3 == true ) {
            int niter = 0;
            int niterw = 0;
            if ( this->dx != this->dz ) {
                while ( niter<nitermax && ( niter<2 || change >= tol || change > prev ) ) {
                    this->sweep_xz(frozen, threadNo);
                    prev = change;
                    change = fsmChange(this->nodes, times, threadNo, tref);
                    tol = fsmTolerance(epsilon, tref, this->nodes.size());
                    niter++;
                }
                if ( niter == nitermax && change >= tol ) {
                    warnFSMnotConverged("first-order", niter, change, epsilon,
                                        tref, this->nodes.size());
                }
                change = std::numeric_limits<T1>::max();
                prev = 0.0;
                while ( niterw<nitermax && ( niterw<2 || change >= tol || change >= prev ) ) {
                    this->sweep_weno3_xz(frozen, threadNo);
                    prev = change;
                    change = fsmChange(this->nodes, times, threadNo, tref);
                    tol = fsmTolerance(epsilon, tref, this->nodes.size());
                    niterw++;
                }
                if ( niterw == nitermax && change >= tol ) {
                    warnFSMnotConverged("WENO3", niterw, change, epsilon,
                                        tref, this->nodes.size());
                }
            } else {
                while ( niter<nitermax && ( niter<2 || change >= tol || change >= prev ) ) {
                    this->sweep(frozen, threadNo);
                    prev = change;
                    change = fsmChange(this->nodes, times, threadNo, tref);
                    tol = fsmTolerance(epsilon, tref, this->nodes.size());
                    niter++;
                }
                if ( niter == nitermax && change >= tol ) {
                    warnFSMnotConverged("first-order", niter, change, epsilon,
                                        tref, this->nodes.size());
                }
                change = std::numeric_limits<T1>::max();
                prev = 0.0;
                while ( niterw<nitermax && ( niterw<2 || change >= tol || change >= prev ) ) {
                    this->sweep_weno3(frozen, threadNo);
                    prev = change;
                    change = fsmChange(this->nodes, times, threadNo, tref);
                    tol = fsmTolerance(epsilon, tref, this->nodes.size());
                    niterw++;
                }
                if ( niterw == nitermax && change >= tol ) {
                    warnFSMnotConverged("WENO3", niterw, change, epsilon,
                                        tref, this->nodes.size());
                }
            }
            niter_final = niter;
            niterw_final = niterw;
        } else {
            int niter = 0;
            while ( niter<nitermax && ( niter<2 || change >= tol || change >= prev ) ) {
                if ( this->dx == this->dz ) {
                    this->sweep(frozen, threadNo);
                    if ( rotated_template == true ) {
                        this->sweep45(frozen, threadNo);
                    }
                } else {
                    this->sweep_xz(frozen, threadNo);
                }

                prev = change;
                change = fsmChange(this->nodes, times, threadNo, tref);
                tol = fsmTolerance(epsilon, tref, this->nodes.size());
                niter++;
            }
            if ( niter == nitermax && change >= tol ) {
                warnFSMnotConverged("first-order", niter, change, epsilon,
                                    tref, this->nodes.size());
            }
            niter_final = niter;
        }
    }


    template<typename T1, typename T2, typename S>
    void Grid2Drcfs<T1,T2,S>::raytrace(const std::vector<S>& Tx,
                                       const std::vector<T1>& t0,
                                       const std::vector<const std::vector<S>*>& Rx,
                                       const size_t threadNo) const {

        this->checkPts(Tx);
        for ( size_t n=0; n<Rx.size(); ++n ) {
            this->checkPts(*Rx[n]);
        }

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        // Set Tx pts
        std::vector<bool> frozen( this->nodes.size(), false );
        int npts = 1;
        if ( weno3 == true) npts = 2;
        this->initFSM(Tx, t0, frozen, npts, threadNo);

        std::vector<T1> times( this->nodes.size() );
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            times[n] = this->nodes[n].getTT( threadNo );
        }

        T1 change = std::numeric_limits<T1>::max();
        T1 prev = 0.0;   // previous sweep's change, for the non-monotone guard
        T1 tref = 0.0;   // traveltime range, refreshed by fsmChange
        T1 tol = 0.0;    // absolute stop threshold, from epsilon and tref
        if ( weno3 == true ) {
            int niter = 0;
            int niterw = 0;
            if ( this->dx != this->dz ) {
                while ( niter<nitermax && ( niter<2 || change >= tol || change >= prev ) ) {
                    this->sweep_xz(frozen, threadNo);
                    prev = change;
                    change = fsmChange(this->nodes, times, threadNo, tref);
                    tol = fsmTolerance(epsilon, tref, this->nodes.size());
                    niter++;
                }
                if ( niter == nitermax && change >= tol ) {
                    warnFSMnotConverged("first-order", niter, change, epsilon,
                                        tref, this->nodes.size());
                }
                change = std::numeric_limits<T1>::max();
                prev = 0.0;
                while ( niterw<nitermax && ( niterw<2 || change >= tol || change >= prev ) ) {
                    this->sweep_weno3_xz(frozen, threadNo);
                    prev = change;
                    change = fsmChange(this->nodes, times, threadNo, tref);
                    tol = fsmTolerance(epsilon, tref, this->nodes.size());
                    niterw++;
                }
                if ( niterw == nitermax && change >= tol ) {
                    warnFSMnotConverged("WENO3", niterw, change, epsilon,
                                        tref, this->nodes.size());
                }
            } else {
                while ( niter<nitermax && ( niter<2 || change >= tol || change >= prev ) ) {
                    this->sweep(frozen, threadNo);
                    prev = change;
                    change = fsmChange(this->nodes, times, threadNo, tref);
                    tol = fsmTolerance(epsilon, tref, this->nodes.size());
                    niter++;
                }
                if ( niter == nitermax && change >= tol ) {
                    warnFSMnotConverged("first-order", niter, change, epsilon,
                                        tref, this->nodes.size());
                }
                change = std::numeric_limits<T1>::max();
                prev = 0.0;
                while ( niterw<nitermax && ( niterw<2 || change >= tol || change >= prev ) ) {
                    this->sweep_weno3(frozen, threadNo);
                    prev = change;
                    change = fsmChange(this->nodes, times, threadNo, tref);
                    tol = fsmTolerance(epsilon, tref, this->nodes.size());
                    niterw++;
                }
                if ( niterw == nitermax && change >= tol ) {
                    warnFSMnotConverged("WENO3", niterw, change, epsilon,
                                        tref, this->nodes.size());
                }
            }
            niter_final = niter;
            niterw_final = niterw;
        } else {
            int niter = 0;
            while ( niter<nitermax && ( niter<2 || change >= tol || change >= prev ) ) {
                if ( this->dx == this->dz ) {
                    this->sweep(frozen, threadNo);
                    if ( rotated_template == true ) {
                        this->sweep45(frozen, threadNo);
                    }
                } else {
                    this->sweep_xz(frozen, threadNo);
                }

                prev = change;
                change = fsmChange(this->nodes, times, threadNo, tref);
                tol = fsmTolerance(epsilon, tref, this->nodes.size());
                niter++;
            }
            if ( niter == nitermax && change >= tol ) {
                warnFSMnotConverged("first-order", niter, change, epsilon,
                                    tref, this->nodes.size());
            }
            niter_final = niter;
        }
    }

}

#endif /* Grid2Drcfs_h */

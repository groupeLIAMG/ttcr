//
//  Grid2Drnfs.h
//  ttcr
//
//  Created by Bernard Giroux on 2015-09-22.
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

/*

 Fast sweeping method (serial) of

 ARTICLE{zhao05,
 author = {Zhao, Hongkai},
 title = {A Fast Sweeping Method for Eikonal Equations},
 journal = {Mathematics of Computation},
 year = {2005},
 volume = {74},
 pages = {603--627},
 number = {250},
 month = apr,
 abstract = {In this paper a fast sweeping method for computing the numerical solution
 of Eikonal equations on a rectangular grid is presented. The method
 is an iterative method which uses upwind difference for discretization
 and uses Gauss-Seidel iterations with alternating sweeping ordering
 to solve the discretized system, The crucial idea is that each sweeping
 ordering follows a family of characteristics of the corresponding
 Eikonal equation in a certain direction simultaneously. The method
 has an optimal complexity of O(N) for N grid points and is extremely
 simple to implement in any number of dimensions. Monotonicity and
 stability properties of the fast sweeping algorithm are proven. Convergence
 and error estimates of the algorithm for computing the distance function
 is studied in detail. It is shown that 2n Gauss-Seidel iterations
 is enough for the distance function in n dimensions. An estimation
 of the number of iterations for general Eikonal equations is also
 studied. Numerical examples are used to verify the analysis.},
 issn = {00255718},
 publisher = {American Mathematical Society},
 url = {http://www.jstor.org/stable/4100081}
 }

 A posteriori ray tracing done following

 @Article{aldridge93,
 Title                    = {Two-dimensional tomographic inversion with finite-difference traveltimes},
 Author                   = {Aldridge, D.F. and Oldenburg, D.W.},
 Journal                  = {J. Seism. Explor.},
 Year                     = {1993},
 Pages                    = {257--274},
 Volume                   = {2}
 }

 */

/**
 * @file Grid2Drnfs.h
 * @brief Fast sweeping solver on a 2-D rectilinear grid with node-based slowness.
 *
 * Declares ttcr::Grid2Drnfs, the serial fast sweeping solver of Zhao (2005),
 * cited in full above. It takes nodal slowness directly, in contrast to
 * ttcr::Grid2Drcfs which accepts a cell model and averages it onto the nodes.
 *
 * The sweep stencils live in Grid2Drn.h; this file supplies the iteration
 * driver and convergence test.
 *
 * @sa Grid2Drn.h, Grid2Drcfs.h, Grid2Drnfs_OpenCL.h
 */

#ifndef ttcr_Grid2Drnfs_h
#define ttcr_Grid2Drnfs_h

#include <cmath>

#include "Grid2Drn.h"
#include "Node2Dn.h"

namespace ttcr {

    /**
     * @brief Fast sweeping eikonal solver with node-based slowness.
     *
     * @tparam T1 floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2 integer type of node and cell indices.
     * @tparam S  point type, @ref sxz or @ref sxyz.
     *
     * Implements Zhao's fast sweeping method (see the citation above): upwind
     * differences solved by Gauss-Seidel iterations with alternating sweep
     * orderings, each ordering following one family of characteristics, giving
     * @f$O(N)@f$ complexity for @f$N@f$ nodes.
     *
     * The sweep and single-node update stencils themselves live in
     * ttcr::Grid2Drn (@ref g2drn_sweeps); this class supplies the iteration
     * driver, the convergence test and the node construction.
     *
     * Unlike ttcr::Grid2Drcfs, which accepts a cell model and averages it onto
     * the nodes, this class takes nodal slowness directly through the inherited
     * ttcr::Grid2Drn::setSlowness — one value per node, no smoothing.
     *
     * @section g2drnfs_conv Convergence
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
     * @sa Grid2Drn.h, Grid2Drcfs.h, Grid2Drnfs_OpenCL.h
     */
    template<typename T1, typename T2, typename S>
    class Grid2Drnfs : public Grid2Drn<T1,T2,S,Node2Dn<T1,T2>> {
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
         *              @sa @ref g2drnfs_conv
         * @param maxit maximum number of sweep iterations.
         * @param w     use the 3rd-order WENO stencil
         *              (ttcr::input_parameters::weno3).
         * @param rt    use rotated stencils as well as axis-aligned ones
         *              (ttcr::input_parameters::rotated_template).
         * @param ttrp  recompute receiver traveltimes along the raypath.
         * @param nt    number of threads.
         *
         * @post Nodes and neighbour lists are built. Slowness is **not** set.
         * @note Both cell sizes are honoured here, unlike the 3-D
         *       ttcr::Grid3Drcfs which is cubic-only.
         */
        Grid2Drnfs(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
                   const T1 minx, const T1 minz, const T1 eps, const int maxit,
                   const bool w, const bool rt,
                   const bool ttrp, const size_t nt=1) :
        Grid2Drn<T1,T2,S,Node2Dn<T1,T2>>(nx,nz,ddx,ddz,minx,minz,ttrp,nt),
        epsilon(eps), nitermax(maxit), niter_final(0), niterw_final(0),
        weno3(w), rotated_template(rt)
        {
            buildGridNodes();
            this->template buildGridNeighbors<Node2Dn<T1,T2>>(this->nodes);
        }

        /// Destructor.
        virtual ~Grid2Drnfs() {
        }

        /// @return Number of sweep iterations the last solve took.
        const int get_niter() const { return niter_final; }
        /// @return Number of WENO sweep iterations the last solve took.
        const int get_niterw() const { return niterw_final; }

    private:
        /// Convergence tolerance, dimensionless: the fraction of the traveltime
        /// range the mean per-node change must fall below.  The loops compare
        /// against an L1 sum, so ttcr::fsmTolerance folds in the node count and
        /// the range each sweep rather than scaling this at construction.
        /// @sa @ref g2drnfs_conv
        T1 epsilon;
        int nitermax;             ///< Iteration cap for the sweeps.
        mutable int niter_final;  ///< Iterations used by the last solve; @c mutable so the @c const raytrace can record it.
        mutable int niterw_final; ///< WENO iterations used by the last solve.
        bool weno3;               ///< Use the 3rd-order WENO stencil.
        bool rotated_template;    ///< Use rotated stencils in addition to axis-aligned ones.

        /// @name Non-copyable
        /// @{
        Grid2Drnfs() {}
        Grid2Drnfs(const Grid2Drnfs<T1,T2,S>& g) {}
        Grid2Drnfs<T1,T2,S>& operator=(const Grid2Drnfs<T1,T2,S>& g) = delete;
        /// @}

        /// @brief Create the grid nodes and set their positions.
        /// @note Primary nodes only; sweeping needs no secondary ones.
        void buildGridNodes();

        /**
         * @brief Propagate the traveltime field and evaluate it at the receivers.
         * @param Tx       source positions.
         * @param t0       origin time of each source.
         * @param Rx       receiver positions.
         * @param threadNo thread to compute on.
         * @note Freezes a halo around each source via
         *       ttcr::Grid2Drn::initFSM, then sweeps until the change falls
         *       below @ref epsilon or @ref nitermax is reached.
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
    void Grid2Drnfs<T1,T2,S>::buildGridNodes() {

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
    void Grid2Drnfs<T1,T2,S>::raytrace(const std::vector<S>& Tx,
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
    void Grid2Drnfs<T1,T2,S>::raytrace(const std::vector<S>& Tx,
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

#endif /* Grid2Drnfs_h */

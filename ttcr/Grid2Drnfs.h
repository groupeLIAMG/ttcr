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

#ifndef ttcr_Grid2Drnfs_h
#define ttcr_Grid2Drnfs_h

#include <cmath>

#include "Grid2Drn.h"
#include "Node2Dn.h"

namespace ttcr {

    template<typename T1, typename T2, typename S>
    class Grid2Drnfs : public Grid2Drn<T1,T2,S,Node2Dn<T1,T2>> {
    public:
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

        virtual ~Grid2Drnfs() {
        }

        const int get_niter() const { return niter_final; }
        const int get_niterw() const { return niterw_final; }

    private:
        T1 epsilon;
        int nitermax;
        mutable int niter_final;
        mutable int niterw_final;
        bool weno3;
        bool rotated_template;

        Grid2Drnfs() {}
        Grid2Drnfs(const Grid2Drnfs<T1,T2,S>& g) {}
        Grid2Drnfs<T1,T2,S>& operator=(const Grid2Drnfs<T1,T2,S>& g) {}

        void buildGridNodes();

        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<S>& Rx,
                      const size_t threadNo=0) const;

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
        if ( weno3 == true ) {
            int niter = 0;
            int niterw = 0;
            if ( this->dx != this->dz ) {
                while ( change >= epsilon && niter<nitermax ) {
                    this->sweep_xz(frozen, threadNo);
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
                    this->sweep_weno3_xz(frozen, threadNo);
                    change = 0.0;
                    for ( size_t n=0; n<this->nodes.size(); ++n ) {
                        T1 dt = std::abs( times[n] - this->nodes[n].getTT(threadNo) );

                        change += dt;
                        times[n] = this->nodes[n].getTT(threadNo);
                    }
                    niterw++;
                }
            } else {
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
            }
            niter_final = niter;
            niterw_final = niterw;
        } else {
            int niter = 0;
            while ( change >= epsilon && niter<nitermax ) {
                if ( this->dx == this->dz ) {
                    this->sweep(frozen, threadNo);
                    if ( rotated_template == true ) {
                        this->sweep45(frozen, threadNo);
                    }
                } else {
                    this->sweep_xz(frozen, threadNo);
                }

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
        if ( weno3 == true ) {
            int niter = 0;
            int niterw = 0;
            if ( this->dx != this->dz ) {
                while ( change >= epsilon && niter<nitermax ) {
                    this->sweep_xz(frozen, threadNo);
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
                    this->sweep_weno3_xz(frozen, threadNo);
                    change = 0.0;
                    for ( size_t n=0; n<this->nodes.size(); ++n ) {
                        T1 dt = std::abs( times[n] - this->nodes[n].getTT(threadNo) );

                        change += dt;
                        times[n] = this->nodes[n].getTT(threadNo);
                    }
                    niterw++;
                }
            } else {
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
            }
            niter_final = niter;
            niterw_final = niterw;
        } else {
            int niter = 0;
            while ( change >= epsilon && niter<nitermax ) {
                if ( this->dx == this->dz ) {
                    this->sweep(frozen, threadNo);
                    if ( rotated_template == true ) {
                        this->sweep45(frozen, threadNo);
                    }
                } else {
                    this->sweep_xz(frozen, threadNo);
                }

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

#endif /* Grid2Drnfs_h */

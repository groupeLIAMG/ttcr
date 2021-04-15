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

#ifndef ttcr_Grid3Drnfs_h
#define ttcr_Grid3Drnfs_h

#include <cmath>
#include <utility>

#include "Grid3Drn.h"
#include "Node3Dn.h"

namespace ttcr {

    template<typename T1, typename T2>
    class Grid3Drnfs : public Grid3Drn<T1,T2,Node3Dn<T1,T2>> {
    public:
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
        }

        ~Grid3Drnfs() {

        }

        const int get_niter() const { return niter_final; }
        const int get_niterw() const { return niterw_final; }

    protected:
        T1 epsilon;
        int nitermax;
        mutable int niter_final;
        mutable int niterw_final;
        bool weno3;

    private:
        Grid3Drnfs() {}
        Grid3Drnfs(const Grid3Drnfs<T1,T2>& g) {}
        Grid3Drnfs<T1,T2>& operator=(const Grid3Drnfs<T1,T2>& g) {}

        void raytrace(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<sxyz<T1>>& Rx,
                      const size_t threadNo=0) const;
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

        this->checkPts(Tx);
        this->checkPts(Rx);

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
            //std::cout << Tx[0] << "    times " << times[0] << '\t' << this->nodes[0].getNodeSlowness() << '\n';
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
        //        for ( size_t n=0; n<this->nodes.size(); ++n ) {
        //            std::cout << this->nodes[n].getTT(threadNo) << '\n';
        //        }
    }

    template<typename T1, typename T2>
    void Grid3Drnfs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<std::vector<sxyz<T1>>>& Rx,
                                     const size_t threadNo) const {

        this->checkPts(Tx);
        for ( size_t n=0; n<Rx.size(); ++n ) {
            this->checkPts(Rx[n]);
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

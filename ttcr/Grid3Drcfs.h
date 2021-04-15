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

#ifndef ttcr_Grid3Drcfs_h
#define ttcr_Grid3Drcfs_h

#include <cmath>
#include <stdexcept>

#include "Grid3Drn.h"
#include "Node3Dn.h"

namespace ttcr {

    template<typename T1, typename T2>
    class Grid3Drcfs : public Grid3Drn<T1,T2,Node3Dn<T1,T2>> {
    public:
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
        }

        virtual ~Grid3Drcfs() {
        }

        void setSlowness(const std::vector<T1>& s);

        const int get_niter() const { return niter_final; }
        const int get_niterw() const { return niterw_final; }

    protected:
        T1 epsilon;
        int nitermax;
        mutable int niter_final;
        mutable int niterw_final;
        bool weno3;

    private:
        Grid3Drcfs() {}
        Grid3Drcfs(const Grid3Drcfs<T1,T2>& g) {}
        Grid3Drcfs<T1,T2>& operator=(const Grid3Drcfs<T1,T2>& g) {}

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
    void Grid3Drcfs<T1,T2>::setSlowness(const std::vector<T1>& s) {

        if ( this->ncx*this->ncy*this->ncz != s.size() ) {
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

#endif /* Grid3Drcfs_h */

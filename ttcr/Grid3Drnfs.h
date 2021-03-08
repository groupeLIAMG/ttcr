//
//  Grid3Drnfs.h
//  ttcr
//
//  Created by Bernard Giroux on 15-12-27.ncz
//  Copyright © 2015 Bernard Giroux. All rights reserved.
//

#ifndef Grid3Drnfs_h
#define Grid3Drnfs_h

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
                   const size_t nt=1) :
        Grid3Drn<T1,T2,Node3Dn<T1,T2>>(nx, ny, nz, ddx, ddx, ddx, minx, miny, minz, ttrp, intVel, nt),
        epsilon(eps), nitermax(maxit), niter_final(0), niterw_final(0), weno3(w)
        {
            buildGridNodes();
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
        
        void buildGridNodes();
        
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
                      const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                      const size_t threadNo=0) const;

    };
    
    template<typename T1, typename T2>
    void Grid3Drnfs<T1,T2>::buildGridNodes() {
        
        T2 cell_XmYmZm; 	// cell in the (x-,y-,z-) direction from the node
        T2 cell_XpYmZm; 	// cell in the (x+,y-,z-) direction from the node
        T2 cell_XmYpZm;
        T2 cell_XpYpZm;
        T2 cell_XmYmZp;
        T2 cell_XpYmZp;
        T2 cell_XmYpZp;
        T2 cell_XpYpZp;
        
        T2 n=0;
        for ( T2 nk=0; nk<=this->ncz; ++nk ) {
            
            T1 z = this->zmin + nk*this->dz;
            
            for ( T2 nj=0; nj<=this->ncy; ++nj ) {
                
                T1 y = this->ymin + nj*this->dy;
                
                for (T2 ni=0; ni<=this->ncx; ++ni){
                    
                    T1 x = this->xmin + ni*this->dx;
                    
                    // Find the adjacent cells for each primary node
                    
                    if (ni < this->ncx && nj < this->ncy && nk < this->ncz){
                        cell_XpYpZp = nj*this->ncx + nk*(this->ncx*this->ncy) + ni;
                    }
                    else {
                        cell_XpYpZp = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni > 0 && nj < this->ncy && nk < this->ncz){
                        cell_XmYpZp = nj*this->ncx + nk*(this->ncx*this->ncy) + ni - 1;
                    }
                    else {
                        cell_XmYpZp = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni < this->ncx && nj > 0 && nk < this->ncz){
                        cell_XpYmZp = (nj-1)*this->ncx + nk*(this->ncx*this->ncy) + ni;
                    }
                    else {
                        cell_XpYmZp = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni > 0 && nj > 0 && nk < this->ncz){
                        cell_XmYmZp = (nj-1)*this->ncx + nk*(this->ncx * this->ncy) + ni - 1;
                    }
                    else {
                        cell_XmYmZp = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni < this->ncx && nj < this->ncy && nk > 0){
                        cell_XpYpZm = nj*this->ncx + (nk-1)*(this->ncx*this->ncy) + ni;
                    }
                    else {
                        cell_XpYpZm = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni > 0 && nj < this->ncy && nk > 0){
                        cell_XmYpZm = nj*this->ncx + (nk-1)*(this->ncx*this->ncy) + ni - 1;
                    }
                    else {
                        cell_XmYpZm = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni < this->ncx && nj > 0 && nk > 0){
                        cell_XpYmZm = (nj-1)*this->ncx + (nk-1)*(this->ncx*this->ncy) + ni;
                    }
                    else {
                        cell_XpYmZm = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni > 0 && nj > 0 && nk > 0){
                        cell_XmYmZm = (nj-1)*this->ncx + (nk-1)*(this->ncx*this->ncy) + ni - 1;
                    }
                    else {
                        cell_XmYmZm = std::numeric_limits<T2>::max();
                    }
                    
                    
                    // Index the primary nodes owners
                    
                    if (cell_XmYmZm != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cell_XmYmZm );
                    }
                    if (cell_XpYmZm != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cell_XpYmZm );
                    }
                    if (cell_XmYpZm != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cell_XmYpZm );
                    }
                    if (cell_XpYpZm != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cell_XpYpZm );
                    }
                    if (cell_XmYmZp != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cell_XmYmZp );
                    }
                    if (cell_XpYmZp != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cell_XpYmZp );
                    }
                    if (cell_XmYpZp != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cell_XmYpZp );
                    }
                    if (cell_XpYpZp != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cell_XpYpZp );
                    }
                    
                    this->nodes[n].setXYZindex( x, y, z, n );
                    this->nodes[n].setPrimary(true);
                    ++n;
                }
            }
        }
    }
    
    
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
//        for ( size_t n=0; n<this->nodes.size(); ++n ) {
//            if ( frozen[n] ) {
//                AtomicWriter aw;
//                aw << Tx[0] << " frozen at " << n << " tt = " << std::setprecision(9) << this->nodes[n].getTT(threadNo) << '\n';
//            }
//        }
        
        std::vector<T1> times( this->nodes.size() );
        for ( size_t n=0; n<this->nodes.size(); ++n )
            times[n] = this->nodes[n].getTT( threadNo );
        
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
                                     const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                                     const size_t threadNo) const {
        
        this->checkPts(Tx);
        for ( size_t n=0; n<Rx.size(); ++n )
        this->checkPts(*Rx[n]);
        
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }
        
        // Set Tx pts
        std::vector<bool> frozen( this->nodes.size(), false );
        int npts = 1;
        if ( weno3 == true ) npts = 2;
        this->initFSM(Tx, t0, frozen, npts, threadNo);
        
        std::vector<T1> times( this->nodes.size() );
        for ( size_t n=0; n<this->nodes.size(); ++n )
        times[n] = this->nodes[n].getTT( threadNo );
        
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

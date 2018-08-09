//
//  Grid3Drnfs.h
//  ttcr
//
//  Created by Bernard Giroux on 15-12-27.ncz
//  Copyright © 2015 Bernard Giroux. All rights reserved.
//

#ifndef Grid3Drnfs_h
#define Grid3Drnfs_h

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
                   const size_t nt=1) :
        Grid3Drn<T1,T2,Node3Dn<T1,T2>>(nx, ny, nz, ddx, ddx, ddx, minx, miny, minz, nt),
        epsilon(eps), nitermax(maxit), niter(0), niterw(0), weno3(w)
        {
            buildGridNodes();
            this->buildGridNeighbors();
        }
        
        ~Grid3Drnfs() {
            
        }
        
        const int get_niter() const { return niter; }
        const int get_niterw() const { return niterw; }
        
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                     const std::vector<T1>& t0,
                     const std::vector<sxyz<T1>>& Rx,
                     std::vector<T1>& traveltimes,
                     const size_t threadNo=0) const;
        
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                     const std::vector<T1>& t0,
                     const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                     std::vector<std::vector<T1>*>& traveltimes,
                     const size_t threadNo=0) const;
        
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                     const std::vector<T1>& t0,
                     const std::vector<sxyz<T1>>& Rx,
                     std::vector<T1>& traveltimes,
                     std::vector<std::vector<sxyz<T1>>>& r_data,
                     const size_t threadNo=0) const;
        
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                     const std::vector<T1>& t0,
                     const std::vector<sxyz<T1>>& Rx,
                     std::vector<T1>& traveltimes,
                     std::vector<std::vector<sxyz<T1>>>& r_data,
                     T1& v0,
                     const size_t threadNo=0) const;
        
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                     const std::vector<T1>& t0,
                     const std::vector<sxyz<T1>>& Rx,
                     std::vector<T1>& traveltimes,
                     std::vector<std::vector<sxyz<T1>>>& r_data,
                     T1& v0,
                     std::vector<std::vector<sijv<T1>>>& m_data,
                     const size_t threadNo=0) const;

        void raytrace(const std::vector<sxyz<T1>>& Tx,
                     const std::vector<T1>& t0,
                     const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                     std::vector<std::vector<T1>*>& traveltimes,
                     std::vector<std::vector<std::vector<sxyz<T1>>>*>& r_data,
                     const size_t threadNo=0) const;
        
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                     const std::vector<T1>& t0,
                     const std::vector<sxyz<T1>>& Rx,
                     std::vector<T1>& traveltimes,
                     std::vector<std::vector<sxyz<T1>>>& r_data,
                     std::vector<std::vector<siv<T1>>>& l_data,
                     const size_t threadNo=0) const;
        
    protected:
        T1 epsilon;
        int nitermax;
        mutable int niter;
        mutable int niterw;
        bool weno3;
        
        void buildGridNodes();
        
    private:
        Grid3Drnfs() {}
        Grid3Drnfs(const Grid3Drnfs<T1,T2>& g) {}
        Grid3Drnfs<T1,T2>& operator=(const Grid3Drnfs<T1,T2>& g) {}
        
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
                    
                    ++n;
                }
            }
        }
    }
    
    
    template<typename T1, typename T2>
    void Grid3Drnfs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<sxyz<T1>>& Rx,
                                     std::vector<T1>& traveltimes,
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
        for ( size_t n=0; n<this->nodes.size(); ++n )
            times[n] = this->nodes[n].getTT( threadNo );
        
        T1 change = std::numeric_limits<T1>::max();
        if ( weno3 == true ) {
            niter=0;
            niterw=0;
            if ( this->dx != this->dz || this->dx != this->dy ) {
                throw std::logic_error("Error: WENO stencil needs dx equal to dz");
            }
            while ( change >= epsilon && niter<nitermax ) {
                this->sweep(frozen, threadNo);
                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
                    
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
                    T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
                    
                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                niterw++;
            }
        } else {
            niter=0;
            while ( change >= epsilon && niter<nitermax ) {
                this->sweep(frozen, threadNo);
                
                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
                    
                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                niter++;
            }
        }
        
        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        
        for (size_t n=0; n<Rx.size(); ++n) {
            traveltimes[n] = this->getTraveltime(Rx[n], threadNo);
        }
    }
    
    template<typename T1, typename T2>
    void Grid3Drnfs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                                     std::vector<std::vector<T1>*>& traveltimes,
                                     const size_t threadNo) const {
        
        this->checkPts(Tx);
        for ( size_t n=0; n<Rx.size(); ++n )
            this->checkPts(*Rx[n]);
        
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }
        
        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        
        // Set Tx pts
        std::vector<bool> frozen( this->nodes.size(), false );
        int npts = 1;
        if ( weno3 == true) npts = 2;
        this->initFSM(Tx, t0, frozen, npts, threadNo);
        
        std::vector<T1> times( this->nodes.size() );
        for ( size_t n=0; n<this->nodes.size(); ++n )
            times[n] = this->nodes[n].getTT( threadNo );
        
        T1 change = std::numeric_limits<T1>::max();
        if ( weno3 == true ) {
            niter=0;
            niterw=0;
            if ( this->dx != this->dz || this->dx != this->dy ) {
                throw std::logic_error("Error: WENO stencil needs dx equal to dz");
            }
            while ( change >= epsilon && niter<nitermax ) {
                this->sweep(frozen, threadNo);
                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
                    
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
                    T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
                    
                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                niterw++;
            }
        } else {
            niter=0;
            while ( change >= epsilon && niter<nitermax ) {
                this->sweep(frozen, threadNo);
                
                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
                    
                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                niter++;
            }
        }
        
        for (size_t nr=0; nr<Rx.size(); ++nr) {
            traveltimes[nr]->resize( Rx[nr]->size() );
            for (size_t n=0; n<Rx[nr]->size(); ++n)
                (*traveltimes[nr])[n] = this->getTraveltime((*Rx[nr])[n], threadNo);
        }
    }
    
    template<typename T1, typename T2>
    void Grid3Drnfs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<sxyz<T1>>& Rx,
                                     std::vector<T1>& traveltimes,
                                     std::vector<std::vector<sxyz<T1>>>& r_data,
                                     const size_t threadNo) const {

        raytrace(Tx, t0, Rx, traveltimes, threadNo);
        
        if ( r_data.size() != Rx.size() ) {
            r_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<r_data.size(); ++ni ) {
            r_data[ni].resize( 0 );
        }
        
        for (size_t n=0; n<Rx.size(); ++n) {
            this->getRaypath(Tx, Rx[n], r_data[n], threadNo);
        }
    }
    
    template<typename T1, typename T2>
    void Grid3Drnfs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<sxyz<T1>>& Rx,
                                     std::vector<T1>& traveltimes,
                                     std::vector<std::vector<sxyz<T1>>>& r_data,
                                     T1& v0,
                                     const size_t threadNo) const {
        
        raytrace(Tx, t0, Rx, traveltimes, threadNo);
        
        if ( r_data.size() != Rx.size() ) {
            r_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<r_data.size(); ++ni ) {
            r_data[ni].resize( 0 );
        }
        
        v0 = 0.0;
        for ( size_t n=0; n<Tx.size(); ++n ) {
            v0 += this->computeSlowness( Tx[n] );
        }
        v0 = Tx.size() / v0;

        for (size_t n=0; n<Rx.size(); ++n) {
            this->getRaypath(Tx, Rx[n], r_data[n], threadNo);
        }
    }
    
    template<typename T1, typename T2>
    void Grid3Drnfs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<sxyz<T1>>& Rx,
                                     std::vector<T1>& traveltimes,
                                     std::vector<std::vector<sxyz<T1>>>& r_data,
                                     T1& v0,
                                     std::vector<std::vector<sijv<T1>>>& m_data,
                                     const size_t threadNo) const {
        
        raytrace(Tx, t0, Rx, traveltimes, threadNo);
        
        if ( r_data.size() != Rx.size() ) {
            r_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<r_data.size(); ++ni ) {
            r_data[ni].resize( 0 );
        }
        if ( m_data.size() != Rx.size() ) {
            m_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<m_data.size(); ++ni ) {
            m_data[ni].resize( 0 );
        }
        
        v0 = 0.0;
        for ( size_t n=0; n<Tx.size(); ++n ) {
            v0 += this->computeSlowness( Tx[n] );
        }
        v0 = Tx.size() / v0;

        for (size_t n=0; n<Rx.size(); ++n) {
            this->getRaypath(Tx, Rx[n], r_data[n], m_data[n], n, threadNo);
        }
    }
    

    template<typename T1, typename T2>
    void Grid3Drnfs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                                     std::vector<std::vector<T1>*>& traveltimes,
                                     std::vector<std::vector<std::vector<sxyz<T1>>>*>& r_data,
                                     const size_t threadNo) const {
        
        raytrace(Tx, t0, Rx, traveltimes, threadNo);
        
        if ( r_data.size() != Rx.size() ) {
            r_data.resize( Rx.size() );
        }
        
        for (size_t nr=0; nr<Rx.size(); ++nr) {
            r_data[nr]->resize( Rx[nr]->size() );
            for ( size_t ni=0; ni<r_data[nr]->size(); ++ni ) {
                (*r_data[nr])[ni].resize( 0 );
            }
            
            for (size_t n=0; n<Rx[nr]->size(); ++n) {
                this->getRaypath(Tx, (*Rx[nr])[n], (*r_data[nr])[n], threadNo);
            }
        }
    }
    
    template<typename T1, typename T2>
    void Grid3Drnfs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<sxyz<T1>>& Rx,
                                     std::vector<T1>& traveltimes,
                                     std::vector<std::vector<sxyz<T1>>>& r_data,
                                     std::vector<std::vector<siv<T1>>>& l_data,
                                     const size_t threadNo) const {
        
        raytrace(Tx, t0, Rx, traveltimes, threadNo);
        
        if ( r_data.size() != Rx.size() ) {
            r_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<r_data.size(); ++ni ) {
            r_data[ni].resize( 0 );
        }
        if ( l_data.size() != Rx.size() ) {
            l_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<l_data.size(); ++ni ) {
            l_data[ni].resize( 0 );
        }
        
        siv<T1> cell;
        for (size_t n=0; n<Rx.size(); ++n) {
            this->getRaypath(Tx, Rx[n], r_data[n], threadNo);
            
            for (size_t ns=0; ns<r_data[n].size()-1; ++ns) {
                sxyz<T1> m = static_cast<T1>(0.5)*(r_data[n][ns]+r_data[n][ns+1]);  // ps @ middle of segment
                cell.i = this->getCellNo( m );
                cell.v = r_data[n][ns].getDistance( r_data[n][ns+1] );
                
                bool found=false;
                for (size_t nc=0; nc<l_data[n].size(); ++nc) {
                    if ( l_data[n][nc].i == cell.i ) {
                        l_data[n][nc].v += cell.v;  // must add in case we pass through secondary nodes along edge
                        found = true;
                        break;
                    }
                }
                if ( found == false ) {
                    l_data[n].push_back( cell );
                }
            }
            //  must be sorted to build matrix L
            sort(l_data[n].begin(), l_data[n].end(), CompareSiv_i<T1>());
        }
    }
}

#endif /* Grid3Drnfs_h */

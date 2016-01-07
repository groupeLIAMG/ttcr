//
//  Grid3Drifs.h
//  ttcr
//
//  Created by Bernard Giroux on 15-12-27.ncz
//  Copyright © 2015 Bernard Giroux. All rights reserved.
//

#ifndef Grid3Drifs_h
#define Grid3Drifs_h

#include <utility>

#include "Grid3Dri.h"
#include "Node3Di.h"

template<typename T1, typename T2>
class Grid3Drifs : public Grid3Dri<T1,T2,Node3Disp<T1,T2>> {
public:
    Grid3Drifs(const T2 nx, const T2 ny, const T2 nz,
               const T1 ddx, const T1 ddy, const T1 ddz,
               const T1 minx, const T1 miny, const T1 minz,
               const T1 eps, const int maxit,
               const size_t nt=1) :
    Grid3Dri<T1,T2,Node3Disp<T1,T2>>(nx, ny, nz, ddx, ddy, ddz, minx, miny, minz, nt),
    epsilon(eps), nitermax(maxit)
    {
        buildGridNodes();
        this->buildGridNeighbors();
    }
    
    ~Grid3Drifs() {
        
    }
    
    int raytrace(const std::vector<sxyz<T1>>& Tx,
                 const std::vector<T1>& t0,
                 const std::vector<sxyz<T1>>& Rx,
                 std::vector<T1>& traveltimes,
                 const size_t threadNo=0) const;
    
    int raytrace(const std::vector<sxyz<T1>>& Tx,
                 const std::vector<T1>& t0,
                 const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                 std::vector<std::vector<T1>*>& traveltimes,
                 const size_t threadNo=0) const;
    
    int raytrace(const std::vector<sxyz<T1>>& Tx,
                 const std::vector<T1>& t0,
                 const std::vector<sxyz<T1>>& Rx,
                 std::vector<T1>& traveltimes,
                 std::vector<std::vector<sxyz<T1>>>& r_data,
                 const size_t threadNo=0) const;
    
    int raytrace(const std::vector<sxyz<T1>>& Tx,
                 const std::vector<T1>& t0,
                 const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                 std::vector<std::vector<T1>*>& traveltimes,
                 std::vector<std::vector<std::vector<sxyz<T1>>>*>& r_data,
                 const size_t threadNo=0) const;
    
    
protected:
    T1 epsilon;
    int nitermax;
    
    void buildGridNodes();

    void updateNode(const size_t, const size_t, const size_t, const size_t=0) const;
    void initTx(const std::vector<sxyz<T1>>& Tx,
                const std::vector<T1>& t0,
                std::vector<bool>& frozen,
                const size_t threadNo) const;
    void sweep(const std::vector<bool>& frozen,
               const size_t threadNo) const;
    
private:
    Grid3Drifs() {}
    Grid3Drifs(const Grid3Drifs<T1,T2>& g) {}
    Grid3Drifs<T1,T2>& operator=(const Grid3Drifs<T1,T2>& g) {}
    
};

template<typename T1, typename T2>
void Grid3Drifs<T1,T2>::buildGridNodes() {
    
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
void Grid3Drifs<T1,T2>::updateNode(const size_t i, const size_t j, const size_t k,
                                   const size_t threadNo) const {
    T1 a1, a2, a3, t;
    
    if (k==0)
        a1 = this->nodes[ ((k+1)*(this->ncy+1)+j)*(this->ncx+1)+i ].getTT(threadNo);
    else if (i==this->ncz)
        a1 = this->nodes[ ((k-1)*(this->ncy+1)+j)*(this->ncx+1)+i ].getTT(threadNo);
    else {
        a1 = this->nodes[ ((k-1)*(this->ncy+1)+j)*(this->ncx+1)+i ].getTT(threadNo);
        t  = this->nodes[ ((k+1)*(this->ncy+1)+j)*(this->ncx+1)+i ].getTT(threadNo);
        a1 = a1<t ? a1 : t;
    }

    if (j==0)
        a2 = this->nodes[ (k*(this->ncy+1)+j+1)*(this->ncx+1)+i ].getTT(threadNo);
    else if (j==this->ncy)
        a2 = this->nodes[ (k*(this->ncy+1)+j-1)*(this->ncx+1)+i ].getTT(threadNo);
    else {
        a2 = this->nodes[ (k*(this->ncy+1)+j-1)*(this->ncx+1)+i ].getTT(threadNo);
        t  = this->nodes[ (k*(this->ncy+1)+j+1)*(this->ncx+1)+i ].getTT(threadNo);
        a2 = a2<t ? a2 : t;
    }
    
    if (i==0)
        a3 = this->nodes[ (k*(this->ncy+1)+j)*(this->ncx+1)+i+1 ].getTT(threadNo);
    else if (k==this->ncx)
        a3 = this->nodes[ (k*(this->ncy+1)+j)*(this->ncx+1)+i-1 ].getTT(threadNo);
    else {
        a3 = this->nodes[ (k*(this->ncy+1)+j)*(this->ncx+1)+i-1 ].getTT(threadNo);
        t  = this->nodes[ (k*(this->ncy+1)+j)*(this->ncx+1)+i+1 ].getTT(threadNo);
        a3 = a3<t ? a3 : t;
    }
    
    if ( a1>a2 ) std::swap(a1, a2);
    if ( a1>a3 ) std::swap(a1, a3);
    if ( a2>a3 ) std::swap(a2, a3);
    
    T1 fh = this->nodes[(k*(this->ncy+1)+j)*(this->ncx+1)+i].getNodeSlowness() * this->dx;

    t = a1 + fh;
    if ( t > a2 ) {
        if (fabs(a1-a2) >= fh) {
            t = (a1<a2 ? a1 : a2) + fh;
        } else {
            t = 0.5*(a1+a2+sqrt(2.*fh*fh - (a1-a2)*(a1-a2)));
        }
    }
    if ( t > a3 ) {
        t = a2 + fh;
        if ( t > a3 ) {
            if (fabs(a2-a3) >= fh) {
                t = (a2<a3 ? a2 : a3) + fh;
            } else {
                t = 0.5*(a2+a3+sqrt(2.*fh*fh - (a2-a3)*(a2-a3)));
            }
        }
    }
    
    if ( t<this->nodes[(k*(this->ncy+1)+j)*(this->ncx+1)+i].getTT(threadNo) )
        this->nodes[(k*(this->ncy+1)+j)*(this->ncx+1)+i].setTT(t,threadNo);

}


template<typename T1, typename T2>
void Grid3Drifs<T1,T2>::initTx(const std::vector<sxyz<T1>>& Tx,
                               const std::vector<T1>& t0,
                               std::vector<bool>& frozen,
                               const size_t threadNo) const {
    
    for (size_t n=0; n<Tx.size(); ++n) {
        bool found = false;
        for ( size_t nn=0; nn<this->nodes.size(); ++nn ) {
            if ( this->nodes[nn] == Tx[n] ) {
                found = true;
                this->nodes[nn].setTT( t0[n], threadNo );
                frozen[nn] = true;
                break;
            }
        }
        if ( found==false ) {

            // find cell where Tx resides
            T2 cellNo = this->getCellNo(Tx[n]);

            T2 ck = cellNo/(this->ncy*this->ncx);
            T2 cj = (cellNo-ck*this->ncy*this->ncx)/this->ncx;
            T2 ci = cellNo - (ck*this->ncy+cj)*this->ncx;
            
            T2 nn = (ck*(this->ncy+1)+cj)*(this->ncx+1)+ci;
            T1 tt = this->nodes[nn].getDistance(Tx[n]) * this->nodes[nn].getNodeSlowness();
            this->nodes[nn].setTT( tt, threadNo );
            frozen[nn] = true;
            
            nn++;
            tt = this->nodes[nn].getDistance(Tx[n]) * this->nodes[nn].getNodeSlowness();
            this->nodes[nn].setTT( tt, threadNo );
            frozen[nn] = true;
            
            nn = (ck*(this->ncy+1)+cj+1)*(this->ncx+1)+ci;
            tt = this->nodes[nn].getDistance(Tx[n]) * this->nodes[nn].getNodeSlowness();
            this->nodes[nn].setTT( tt, threadNo );
            frozen[nn] = true;
            
            nn++;
            tt = this->nodes[nn].getDistance(Tx[n]) * this->nodes[nn].getNodeSlowness();
            this->nodes[nn].setTT( tt, threadNo );
            frozen[nn] = true;
            
            nn = ((ck+1)*(this->ncy+1)+cj)*(this->ncx+1)+ci;
            tt = this->nodes[nn].getDistance(Tx[n]) * this->nodes[nn].getNodeSlowness();
            this->nodes[nn].setTT( tt, threadNo );
            frozen[nn] = true;
            
            nn++;
            tt = this->nodes[nn].getDistance(Tx[n]) * this->nodes[nn].getNodeSlowness();
            this->nodes[nn].setTT( tt, threadNo );
            frozen[nn] = true;
            
            nn = ((ck+1)*(this->ncy+1)+cj+1)*(this->ncx+1)+ci;
            tt = this->nodes[nn].getDistance(Tx[n]) * this->nodes[nn].getNodeSlowness();
            this->nodes[nn].setTT( tt, threadNo );
            frozen[nn] = true;
            
            nn++;
            tt = this->nodes[nn].getDistance(Tx[n]) * this->nodes[nn].getNodeSlowness();
            this->nodes[nn].setTT( tt, threadNo );
            frozen[nn] = true;

        }
    }
}


template<typename T1, typename T2>
int Grid3Drifs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                const std::vector<T1>& t0,
                                const std::vector<sxyz<T1>>& Rx,
                                std::vector<T1>& traveltimes,
                                const size_t threadNo) const {
    
    if ( this->check_pts(Tx) == 1 ) return 1;
    if ( this->check_pts(Rx) == 1 ) return 1;
    
    for ( size_t n=0; n<this->nodes.size(); ++n ) {
        this->nodes[n].reinit( threadNo );
    }
    
    // Set Tx pts
    std::vector<bool> frozen( this->nodes.size(), false );
    initTx(Tx, t0, frozen, threadNo);
    
    T1 change = std::numeric_limits<T1>::max();
    std::vector<T1> times( this->nodes.size() );
    for ( size_t n=0; n<this->nodes.size(); ++n )
        times[n] = this->nodes[n].getTT( threadNo );
    
    int niter=0;
    while ( change >= epsilon && niter<nitermax ) {
        sweep(frozen, threadNo);
        
        change = 0.0;
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
            
            change += dt;
            times[n] = this->nodes[n].getTT(threadNo);
        }
    }
    std::cout << niter << " iterations were needed with epsilon = " << epsilon << '\n';
    
    if ( traveltimes.size() != Rx.size() ) {
        traveltimes.resize( Rx.size() );
    }
    
    for (size_t n=0; n<Rx.size(); ++n) {
        traveltimes[n] = this->getTraveltime(Rx[n], this->nodes, threadNo);
    }
    return 0;
}

template<typename T1, typename T2>
int Grid3Drifs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                const std::vector<T1>& t0,
                                const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                                std::vector<std::vector<T1>*>& traveltimes,
                                const size_t threadNo) const {
    
    if ( this->check_pts(Tx) == 1 ) return 1;
    for ( size_t n=0; n<Rx.size(); ++n )
        if ( this->check_pts(*Rx[n]) == 1 ) return 1;
    
    for ( size_t n=0; n<this->nodes.size(); ++n ) {
        this->nodes[n].reinit( threadNo );
    }
    
    if ( traveltimes.size() != Rx.size() ) {
        traveltimes.resize( Rx.size() );
    }
    
    // Set Tx pts
    std::vector<bool> frozen( this->nodes.size(), false );
    initTx(Tx, t0, frozen, threadNo);
    
    T1 change = std::numeric_limits<T1>::max();
    std::vector<T1> times( this->nodes.size() );
    for ( size_t n=0; n<this->nodes.size(); ++n )
        times[n] = this->nodes[n].getTT( threadNo );
    
    int niter=0;
    while ( change >= epsilon && niter<nitermax ) {
        sweep(frozen, threadNo);
        
        change = 0.0;
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
            
            change += dt;
            times[n] = this->nodes[n].getTT(threadNo);
        }
        niter++;
    }
    std::cout << niter << " iterations were needed with epsilon = " << epsilon << '\n';
    
    for (size_t nr=0; nr<Rx.size(); ++nr) {
        traveltimes[nr]->resize( Rx[nr]->size() );
        for (size_t n=0; n<Rx[nr]->size(); ++n)
            (*traveltimes[nr])[n] = this->getTraveltime((*Rx[nr])[n], this->nodes, threadNo);
    }
    return 0;
}

template<typename T1, typename T2>
int Grid3Drifs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                const std::vector<T1>& t0,
                                const std::vector<sxyz<T1>>& Rx,
                                std::vector<T1>& traveltimes,
                                std::vector<std::vector<sxyz<T1>>>& r_data,
                                const size_t threadNo) const {
    
    if ( this->check_pts(Tx) == 1 ) return 1;
    if ( this->check_pts(Rx) == 1 ) return 1;
    
    for ( size_t n=0; n<this->nodes.size(); ++n ) {
        this->nodes[n].reinit( threadNo );
    }
    
    // Set Tx pts
    std::vector<bool> frozen( this->nodes.size(), false );
    initTx(Tx, t0, frozen, threadNo);
    
    T1 change = std::numeric_limits<T1>::max();
    std::vector<T1> times( this->nodes.size() );
    for ( size_t n=0; n<this->nodes.size(); ++n )
        times[n] = this->nodes[n].getTT( threadNo );
    
    int niter=0;
    while ( change >= epsilon && niter<nitermax ) {
        sweep(frozen, threadNo);
        
        change = 0.0;
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
            
            change += dt;
            times[n] = this->nodes[n].getTT(threadNo);
        }
    }
    std::cout << niter << " iterations were needed with epsilon = " << epsilon << '\n';
    
    if ( traveltimes.size() != Rx.size() ) {
        traveltimes.resize( Rx.size() );
    }
    if ( r_data.size() != Rx.size() ) {
        r_data.resize( Rx.size() );
    }
    
    for (size_t n=0; n<Rx.size(); ++n) {
        traveltimes[n] = this->getTraveltime(Rx[n], this->nodes, threadNo);
        this->getRaypath(Tx, Rx[n], r_data[n], threadNo);
    }
    return 0;
}

template<typename T1, typename T2>
int Grid3Drifs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                const std::vector<T1>& t0,
                                const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                                std::vector<std::vector<T1>*>& traveltimes,
                                std::vector<std::vector<std::vector<sxyz<T1>>>*>& r_data,
                                const size_t threadNo) const {
    
    if ( this->check_pts(Tx) == 1 ) return 1;
    for ( size_t n=0; n<Rx.size(); ++n )
        if ( this->check_pts(*Rx[n]) == 1 ) return 1;
    
    for ( size_t n=0; n<this->nodes.size(); ++n ) {
        this->nodes[n].reinit( threadNo );
    }
    
    if ( traveltimes.size() != Rx.size() ) {
        traveltimes.resize( Rx.size() );
    }
    
    // Set Tx pts
    std::vector<bool> frozen( this->nodes.size(), false );
    initTx(Tx, t0, frozen, threadNo);
    
    T1 change = std::numeric_limits<T1>::max();
    std::vector<T1> times( this->nodes.size() );
    for ( size_t n=0; n<this->nodes.size(); ++n )
        times[n] = this->nodes[n].getTT( threadNo );
    
    int niter=0;
    while ( change >= epsilon && niter<nitermax ) {
        sweep(frozen, threadNo);
        
        change = 0.0;
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
            
            change += dt;
            times[n] = this->nodes[n].getTT(threadNo);
        }
        niter++;
    }
    std::cout << niter << " iterations were needed with epsilon = " << epsilon << '\n';
    
    if ( traveltimes.size() != Rx.size() ) {
        traveltimes.resize( Rx.size() );
    }
    if ( r_data.size() != Rx.size() ) {
        r_data.resize( Rx.size() );
    }
    
    for (size_t nr=0; nr<Rx.size(); ++nr) {
        traveltimes[nr]->resize( Rx[nr]->size() );
        r_data[nr]->resize( Rx[nr]->size() );
        for ( size_t ni=0; ni<r_data[nr]->size(); ++ni ) {
            (*r_data[nr])[ni].resize( 0 );
        }
        
        for (size_t n=0; n<Rx[nr]->size(); ++n) {
            
            (*traveltimes[nr])[n] = this->getTraveltime((*Rx[nr])[n], this->nodes, threadNo);
            
            this->getRaypath(Tx, (*Rx[nr])[n], (*r_data[nr])[n], threadNo);
        }
    }
    return 0;
}

template<typename T1, typename T2>
void Grid3Drifs<T1,T2>::sweep(const std::vector<bool>& frozen,
                              const size_t threadNo) const {
    
    // sweep first direction
    for ( size_t k=0; k<=this->ncz; ++k ) {
        for ( size_t j=0; j<=this->ncy; ++j ) {
            for ( size_t i=0; i<=this->ncx; ++i ) {
                if ( !frozen[ (k*(this->ncy+1)+j)*(this->ncx+1)+i ] ) {
                    updateNode(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep second direction
    for ( size_t k=0; k<=this->ncz; ++k ) {
        for ( size_t j=0; j<=this->ncy; ++j ) {
            for ( long int i=this->ncx; i>=0; --i ) {
                if ( !frozen[ (k*(this->ncy+1)+j)*(this->ncx+1)+i ] ) {
                    updateNode(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep third direction
    for ( size_t k=0; k<=this->ncz; ++k ) {
        for ( long int j=this->ncy; j>=0; --j ) {
            for ( size_t i=0; i<=this->ncx; ++i ) {
                if ( !frozen[ (k*(this->ncy+1)+j)*(this->ncx+1)+i ] ) {
                    updateNode(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep fourth direction
    for ( size_t k=0; k<=this->ncz; ++k ) {
        for ( long int j=this->ncy; j>=0; --j ) {
            for ( long int i=this->ncx; i>=0; --i ) {
                if ( !frozen[ (k*(this->ncy+1)+j)*(this->ncx+1)+i ] ) {
                    updateNode(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep fifth direction
    for ( long int k=this->ncz; k>=0; --k ) {
        for ( size_t j=0; j<=this->ncy; ++j ) {
            for ( size_t i=0; i<=this->ncx; ++i ) {
                if ( !frozen[ (k*(this->ncy+1)+j)*(this->ncx+1)+i ] ) {
                    updateNode(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep sixth direction
    for ( long int k=this->ncz; k>=0; --k ) {
        for ( size_t j=0; j<=this->ncy; ++j ) {
            for ( long int i=this->ncx; i>=0; --i ) {
                if ( !frozen[ (k*(this->ncy+1)+j)*(this->ncx+1)+i ] ) {
                    updateNode(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep seventh direction
    for ( long int k=this->ncz; k>=0; --k ) {
        for ( long int j=this->ncy; j>=0; --j ) {
            for ( size_t i=0; i<=this->ncx; ++i ) {
                if ( !frozen[ (k*(this->ncy+1)+j)*(this->ncx+1)+i ] ) {
                    updateNode(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep eighth direction
    for ( long int k=this->ncz; k>=0; --k ) {
        for ( long int j=this->ncy; j>=0; --j ) {
            for ( long int i=this->ncx; i>=0; --i ) {
                if ( !frozen[ (k*(this->ncy+1)+j)*(this->ncx+1)+i ] ) {
                    updateNode(i, j, k, threadNo);
                }
            }
        }
    }
}




#endif /* Grid3Drifs_h */

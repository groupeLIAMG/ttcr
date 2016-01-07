//
//  Grid2Drcfs.h
//  ttcr
//
//  Created by Bernard Giroux on 15-12-23.
//  Copyright © 2015 Bernard Giroux. All rights reserved.
//

#ifndef Grid2Drcfs_h
#define Grid2Drcfs_h

#include "Grid2Dri.h"
#include "Node2Di.h"

template<typename T1, typename T2>
class Grid2Drcfs : public Grid2Dri<T1,T2,Node2Di<T1,T2>> {
public:
    Grid2Drcfs(const T2 nx, const T2 nz, const T1 ddx,
               const T1 minx, const T1 minz, const T1 eps, const int maxit,
               const size_t nt=1);
    
    virtual ~Grid2Drcfs() {
    }
    
    int setSlowness(const std::vector<T1>& s);
    
    int raytrace(const std::vector<sxz<T1>>& Tx,
                 const std::vector<T1>& t0,
                 const std::vector<sxz<T1>>& Rx,
                 std::vector<T1>& traveltimes,
                 const size_t threadNo=0) const;
    
    int raytrace(const std::vector<sxz<T1>>& Tx,
                 const std::vector<T1>& t0,
                 const std::vector<const std::vector<sxz<T1>>*>& Rx,
                 std::vector<std::vector<T1>*>& traveltimes,
                 const size_t threadNo=0) const;
    
    int raytrace(const std::vector<sxz<T1>>& Tx,
                 const std::vector<T1>& t0,
                 const std::vector<sxz<T1>>& Rx,
                 std::vector<T1>& traveltimes,
                 std::vector<std::vector<sxz<T1>>>& r_data,
                 const size_t threadNo=0) const;
    
    int raytrace(const std::vector<sxz<T1>>& Tx,
                 const std::vector<T1>& t0,
                 const std::vector<const std::vector<sxz<T1>>*>& Rx,
                 std::vector<std::vector<T1>*>& traveltimes,
                 std::vector<std::vector<std::vector<sxz<T1>>>*>& r_data,
                 const size_t threadNo=0) const;
    
    
protected:
    T1 epsilon;
    int nitermax;
    
    void buildGridNodes();
    
    void updateNode(const size_t, const size_t, const size_t=0) const;
    void initTx(const std::vector<sxz<T1>>& Tx,
                const std::vector<T1>& t0,
                std::vector<bool>& frozen,
                const size_t threadNo) const;
    void sweep(const std::vector<bool>& frozen,
               const size_t threadNo) const;
    
private:
    Grid2Drcfs() {}
    Grid2Drcfs(const Grid2Drcfs<T1,T2>& g) {}
    Grid2Drcfs<T1,T2>& operator=(const Grid2Drcfs<T1,T2>& g) {}
    
};


template<typename T1, typename T2>
Grid2Drcfs<T1,T2>::Grid2Drcfs(const T2 nx, const T2 nz, const T1 ddx,
                              const T1 minx, const T1 minz,
                              const T1 eps, const int maxit,
                              const size_t nt) :
Grid2Dri<T1,T2,Node2Di<T1,T2>>(nx,nz,ddx,ddx,minx,minz,nt),
epsilon(eps), nitermax(maxit)
{
    buildGridNodes();
    this->buildGridNeighbors();
}


template<typename T1, typename T2>
int Grid2Drcfs<T1,T2>::setSlowness(const std::vector<T1>& s) {
    
    if ( this->nCellx*this->nCellz != s.size() ) {
        std::cerr << "Error: slowness vectors of incompatible size.";
        return 1;
    }
    
    // interpolate slowness at grid nodes
    
    const size_t nx = this->nCellx;
    const size_t nz = this->nCellz;
    
    // four courners
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
            this->nodes[i*(nz+1)+j].setNodeSlowness(0.25*(
                                                          s[i*nz+j]+
                                                          s[i*nz+j-1]+
                                                          s[(i-1)*nz+j]+
                                                          s[(i-1)*nz+j-1]));
        }
    }
    
//    if ( false ) {
//        for ( size_t i=0; i<nx; ++i ) {
//            for ( size_t j=0; j<nz; ++j ) {
//                std::cout << (this->xmin+(i+0.5)*this->dx) << '\t' << (this->zmin+(j+0.5)*this->dz) << '\t' << s[i*nz+j] << '\n';
//            }
//        }
//        std::cout << '\n';
//        std::cout << '\n';
//        
//        for ( size_t i=0; i<nx+1; ++i ) {
//            for ( size_t j=0; j<nz+1; ++j ) {
//                std::cout << this->nodes[i*(nz+1)+j].getX() << '\t' << this->nodes[i*(nz+1)+j].getZ() << '\t' << this->nodes[i*(nz+1)+j].getNodeSlowness() << '\n';
//            }
//        }
//    }
    
    return 0;
}


template<typename T1, typename T2>
void Grid2Drcfs<T1,T2>::buildGridNodes() {
    
    T2 cell_upLeft = std::numeric_limits<T2>::max();
    T2 cell_upRight = std::numeric_limits<T2>::max();
    T2 cell_downLeft = 0;
    T2 cell_downRight = 0;
    
    for ( T2 n=0, nc=0; nc<=this->nCellx; ++nc ) {
        
        double x = this->xmin + nc*this->dx;
        
        for ( T2 nr=0; nr<=this->nCellz; ++nr ) {
            
            double z = this->zmin + nr*this->dz;
            
            if ( nr < this->nCellz && nc < this->nCellx ) {
                cell_downRight = nc*this->nCellz + nr;
            }
            else {
                cell_downRight = std::numeric_limits<T2>::max();
            }
            
            if ( nr > 0 && nc < this->nCellx ) {
                cell_upRight = nc*this->nCellz + nr - 1;
            }
            else {
                cell_upRight = std::numeric_limits<T2>::max();
            }
            
            if ( nr < this->nCellz && nc > 0 ) {
                cell_downLeft = (nc-1)*this->nCellz + nr;
            }
            else {
                cell_downLeft = std::numeric_limits<T2>::max();
            }
            
            if ( nr > 0 && nc > 0 ) {
                cell_upLeft = (nc-1)*this->nCellz + nr - 1;
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
            
            ++n;
        }
    }
}


template<typename T1, typename T2>
void Grid2Drcfs<T1,T2>::updateNode(const size_t i, const size_t j,
                                   const size_t threadNo) const {
    T1 a, b, t;
    if (i==0)
        a = this->nodes[ (i+1)*(this->nCellz+1)+j ].getTT(threadNo);
    else if (i==this->nCellx)
        a = this->nodes[ (i-1)*(this->nCellz+1)+j ].getTT(threadNo);
    else {
        a = this->nodes[ (i-1)*(this->nCellz+1)+j ].getTT(threadNo);
        t = this->nodes[ (i+1)*(this->nCellz+1)+j ].getTT(threadNo);
        a = a<t ? a : t;
    }
    
    if (j==0)
        b = this->nodes[ i*(this->nCellz+1)+j+1 ].getTT(threadNo);
    else if (j==this->nCellz)
        b = this->nodes[ i*(this->nCellz+1)+j-1 ].getTT(threadNo);
    else {
        b = this->nodes[ i*(this->nCellz+1)+j-1 ].getTT(threadNo);
        t = this->nodes[ i*(this->nCellz+1)+j+1 ].getTT(threadNo);
        b = b<t ? b : t;
    }
    
    T1 fh = this->nodes[i*(this->nCellz+1)+j].getNodeSlowness() *
    this->dx;
    
    if ( fabs(a-b) >= fh )
        t = (a<b ? a : b) + fh;
    else
        t = 0.5*( a+b + sqrt(2.*fh*fh - (a-b)*(a-b) ));
    
    if ( t<this->nodes[i*(this->nCellz+1)+j].getTT(threadNo) )
        this->nodes[i*(this->nCellz+1)+j].setTT(t,threadNo);
    
}


template<typename T1, typename T2>
void Grid2Drcfs<T1,T2>::initTx(const std::vector<sxz<T1>>& Tx,
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
            
            T2 ci = cellNo/this->nCellz;
            T2 ck = cellNo - ci*this->nCellz;
            
            // upper left node
            T2 nn = ci*(this->nCellz+1) + ck;
            T1 tt = this->nodes[nn].getDistance(Tx[n]) * this->nodes[nn].getNodeSlowness();
            this->nodes[nn].setTT( tt, threadNo );
            frozen[nn] = true;
            
            // lower left node
            nn++;
            tt = this->nodes[nn].getDistance(Tx[n]) * this->nodes[nn].getNodeSlowness();
            this->nodes[nn].setTT( tt, threadNo );
            frozen[nn] = true;
            
            // lower right node
            nn += (this->nCellz+1);
            tt = this->nodes[nn].getDistance(Tx[n]) * this->nodes[nn].getNodeSlowness();
            this->nodes[nn].setTT( tt, threadNo );
            frozen[nn] = true;
            
            // upper right node
            nn--;
            tt = this->nodes[nn].getDistance(Tx[n]) * this->nodes[nn].getNodeSlowness();
            this->nodes[nn].setTT( tt, threadNo );
            frozen[nn] = true;
        }
    }
}

template<typename T1, typename T2>
int Grid2Drcfs<T1,T2>::raytrace(const std::vector<sxz<T1>>& Tx,
                                const std::vector<T1>& t0,
                                const std::vector<sxz<T1>>& Rx,
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
int Grid2Drcfs<T1,T2>::raytrace(const std::vector<sxz<T1>>& Tx,
                                const std::vector<T1>& t0,
                                const std::vector<const std::vector<sxz<T1>>*>& Rx,
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
int Grid2Drcfs<T1,T2>::raytrace(const std::vector<sxz<T1>>& Tx,
                                const std::vector<T1>& t0,
                                const std::vector<sxz<T1>>& Rx,
                                std::vector<T1>& traveltimes,
                                std::vector<std::vector<sxz<T1>>>& r_data,
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
int Grid2Drcfs<T1,T2>::raytrace(const std::vector<sxz<T1>>& Tx,
                                const std::vector<T1>& t0,
                                const std::vector<const std::vector<sxz<T1>>*>& Rx,
                                std::vector<std::vector<T1>*>& traveltimes,
                                std::vector<std::vector<std::vector<sxz<T1>>>*>& r_data,
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
void Grid2Drcfs<T1,T2>::sweep(const std::vector<bool>& frozen,
                              const size_t threadNo) const {
    
    // sweep first direction
    for ( size_t i=0; i<=this->nCellx; ++i ) {
        for ( size_t j=0; j<=this->nCellz; ++j ) {
            if ( !frozen[ i*(this->nCellz+1)+j ] ) {
                updateNode(i, j, threadNo);
            }
        }
    }
    // sweep second direction
    for ( long int i=this->nCellx; i>=0; --i ) {
        for ( size_t j=0; j<=this->nCellz; ++j ) {
            if ( !frozen[ i*(this->nCellz+1)+j ] ) {
                updateNode(i, j, threadNo);
            }
        }
    }
    // sweep third direction
    for ( long int i=this->nCellx; i>=0; --i ) {
        for ( long int j=this->nCellz; j>=0; --j ) {
            if ( !frozen[ i*(this->nCellz+1)+j ] ) {
                updateNode(i, j, threadNo);
            }
        }
    }
    // sweep fourth direction
    for ( size_t i=0; i<=this->nCellx; ++i ) {
        for ( long int j=this->nCellz; j>=0; --j ) {
            if ( !frozen[ i*(this->nCellz+1)+j ] ) {
                updateNode(i, j, threadNo);
            }
        }
    }
}



#endif /* Grid2Drcfs_h */


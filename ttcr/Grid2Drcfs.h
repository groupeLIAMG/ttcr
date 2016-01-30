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
    Grid2Drcfs(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
               const T1 minx, const T1 minz, const T1 eps, const int maxit,
               const bool w, const bool rt, const size_t nt=1);
    
    virtual ~Grid2Drcfs() {
    }
    
    int setSlowness(const std::vector<T1>& s);

    const int get_niter() const { return niter; }
    const int get_niterw() const { return niterw; }

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

    int raytrace(const std::vector<sxz<T1>>& Tx,
                 const std::vector<T1>& t0,
                 const std::vector<sxz<T1>>& Rx,
                 std::vector<T1>& traveltimes,
                 std::vector<std::vector<sxz<T1>>>& r_data,
                 std::vector<std::vector<siv<T1>>>& l_data,
                 const size_t threadNo=0) const;
    
protected:
    T1 epsilon;
    int nitermax;
    mutable int niter;
    mutable int niterw;
    bool weno3;
    bool rotated_template;
    
    void buildGridNodes();
        
private:
    Grid2Drcfs() {}
    Grid2Drcfs(const Grid2Drcfs<T1,T2>& g) {}
    Grid2Drcfs<T1,T2>& operator=(const Grid2Drcfs<T1,T2>& g) {}
    
};


template<typename T1, typename T2>
Grid2Drcfs<T1,T2>::Grid2Drcfs(const T2 nx, const T2 nz,
                              const T1 ddx, const T1 ddz,
                              const T1 minx, const T1 minz,
                              const T1 eps, const int maxit, const bool w,
                              const bool rt, const size_t nt) :
Grid2Dri<T1,T2,Node2Di<T1,T2>>(nx,nz,ddx,ddz,minx,minz,nt),
epsilon(eps), nitermax(maxit), niter(0), niterw(0), weno3(w), rotated_template(rt)
{
    buildGridNodes();
    this->buildGridNeighbors();
}


template<typename T1, typename T2>
int Grid2Drcfs<T1,T2>::setSlowness(const std::vector<T1>& s) {
    
    if ( this->ncx*this->ncz != s.size() ) {
        std::cerr << "Error: slowness vectors of incompatible size.";
        return 1;
    }
    
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
    
    for ( T2 n=0, nc=0; nc<=this->ncx; ++nc ) {
        
        double x = this->xmin + nc*this->dx;
        
        for ( T2 nr=0; nr<=this->ncz; ++nr ) {
            
            double z = this->zmin + nr*this->dz;
            
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
            
            ++n;
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
        if ( this->dx != this->dz ) {
            while ( change >= epsilon && niter<nitermax ) {
                this->sweep_xz(frozen, threadNo);
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
                this->sweep_weno3_xz(frozen, threadNo);
                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
                    
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
        }
    } else {
        niter=0;
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
        if ( this->dx != this->dz ) {
            while ( change >= epsilon && niter<nitermax ) {
                this->sweep_xz(frozen, threadNo);
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
                this->sweep_weno3_xz(frozen, threadNo);
                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
                    
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
        }
    } else {
        niter=0;
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
    return 0;
}

template<typename T1, typename T2>
int Grid2Drcfs<T1,T2>::raytrace(const std::vector<sxz<T1>>& Tx,
                                const std::vector<T1>& t0,
                                const std::vector<sxz<T1>>& Rx,
                                std::vector<T1>& traveltimes,
                                std::vector<std::vector<sxz<T1>>>& r_data,
                                const size_t threadNo) const {
    
    int check = raytrace(Tx, t0, Rx, traveltimes, threadNo);
    if ( check == 1 ) return 1;

    if ( r_data.size() != Rx.size() ) {
        r_data.resize( Rx.size() );
    }
    for ( size_t ni=0; ni<r_data.size(); ++ni ) {
        r_data[ni].resize( 0 );
    }

    for (size_t n=0; n<Rx.size(); ++n) {
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
    
    int check = raytrace(Tx, t0, Rx, traveltimes, threadNo);
    if ( check == 1 ) return 1;

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
    return 0;
}

template<typename T1, typename T2>
int Grid2Drcfs<T1,T2>::raytrace(const std::vector<sxz<T1>>& Tx,
                                const std::vector<T1>& t0,
                                const std::vector<sxz<T1>>& Rx,
                                std::vector<T1>& traveltimes,
                                std::vector<std::vector<sxz<T1>>>& r_data,
                                std::vector<std::vector<siv<T1>>>& l_data,
                                const size_t threadNo) const {
    
    int check = raytrace(Tx, t0, Rx, traveltimes, threadNo);
    if ( check == 1 ) return 1;
    
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
            sxz<T1> m = 0.5*(r_data[n][ns]+r_data[n][ns+1]);  // ps @ middle of segment
            cell.i = this->getCellNo( m );
            cell.v = r_data[n][ns].getDistance( r_data[n][ns+1] );
            
            bool found=false;
            for (size_t nc=0; nc<l_data[n].size(); ++nc) {
                if ( l_data[n][nc].i == cell.i ) {
                    l_data[n][nc].v += cell.v;
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
    return 0;
}



#endif /* Grid2Drcfs_h */


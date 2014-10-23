//
//  Grid3Duifs.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-04-21.
//  Copyright (c) 2014 Bernard Giroux. All rights reserved.
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

#ifndef ttcr_Grid3Duifs_h
#define ttcr_Grid3Duifs_h

#include <cmath>
#include <fstream>
#include <queue>
#include <vector>

#include "Grid3Dui.h"
#include "Node3Di.h"
#include "Metric.h"

template<typename T1, typename T2>
class Grid3Duifs : public Grid3Dui<T1,T2,Node3Di<T1,T2>> {
public:
	Grid3Duifs(const std::vector<sxyz<T1>>& no,
			   const std::vector<tetrahedronElem<T2>>& tet,
			   const T1 eps, const int maxit, const bool rp=false,
               const size_t nt=1) :
    Grid3Dui<T1,T2,Node3Di<T1,T2>>(no, tet, nt),
	rp_ho(rp), epsilon(eps), nitermax(maxit), S()
	{
		buildGridNodes(no, nt);
		this->buildGridNeighbors();
	}
	
	~Grid3Duifs() {
	}
	
	void initOrdering(const std::vector<sxyz<T1>>& refPts, const int order);
	
	int raytrace(const std::vector<sxyz<T1>>& Tx,
				 const std::vector<T1>& t0,
				 const std::vector<sxyz<T1>>& Rx,
				 std::vector<T1>& traveltimes,
				 const size_t threadNo=0) const;
	
	int raytrace(const std::vector<sxyz<T1>>&,
				 const std::vector<T1>&,
				 const std::vector<const std::vector<sxyz<T1>>*>&,
				 std::vector<std::vector<T1>*>&,
				 const size_t=0) const;
    
	int raytrace(const std::vector<sxyz<T1>>&,
                 const std::vector<T1>& ,
                 const std::vector<sxyz<T1>>&,
                 std::vector<T1>&,
                 std::vector<std::vector<sxyz<T1>>>&,
				 const size_t=0) const;
    
    int raytrace(const std::vector<sxyz<T1>>&,
                 const std::vector<T1>&,
                 const std::vector<const std::vector<sxyz<T1>>*>&,
                 std::vector<std::vector<T1>*>&,
                 std::vector<std::vector<std::vector<sxyz<T1>>>*>&,
                 const size_t=0) const;
	
private:
    bool rp_ho;
	T1 epsilon;
	int nitermax;
	std::vector<std::vector<Node3Di<T1,T2>*>> S;
	
	void buildGridNodes(const std::vector<sxyz<T1>>&, const size_t);
	
	void initTx(const std::vector<sxyz<T1>>& Tx, const std::vector<T1>& t0,
				std::vector<bool>& frozen, const size_t threadNo) const;
	
	void initBand(const std::vector<sxyz<T1>>& Tx,
				  const std::vector<T1>& t0,
				  std::priority_queue<Node3Di<T1,T2>*,
				  std::vector<Node3Di<T1,T2>*>,
				  CompareNodePtr<T1>>&,
				  std::vector<Node3Di<T1,T2>>&,
				  std::vector<bool>&,
				  std::vector<bool>&,
				  const size_t) const;
	
	void propagate(std::priority_queue<Node3Di<T1,T2>*,
				   std::vector<Node3Di<T1,T2>*>,
				   CompareNodePtr<T1>>&,
				   std::vector<bool>&,
				   std::vector<bool>&,
				   const size_t) const;
	
};

template<typename T1, typename T2>
void Grid3Duifs<T1,T2>::buildGridNodes(const std::vector<sxyz<T1>>& no,
									   const size_t nt) {
	
	// primary nodes
	for ( T2 n=0; n<no.size(); ++n ) {
		this->nodes[n].setXYZindex( no[n].x, no[n].y, no[n].z, n );
	}
	
	for ( T2 ntet=0; ntet<this->tetrahedra.size(); ++ntet ) {
		
		// for each triangle
		for ( T2 ntri=0; ntri<4; ++ntri ) {
			
			// push owner for primary nodes
			this->nodes[ this->tetrahedra[ntet].i[ntri] ].pushOwner( ntet );
			
		}
	}
}

template<typename T1, typename T2>
void Grid3Duifs<T1,T2>::initOrdering(const std::vector<sxyz<T1>>& refPts,
									 const int order) {
	S.resize( refPts.size() );
	
	Metric<T1> *m;
	if ( order == 1 )
		m = new Metric1<T1>();
	else
		m = new Metric2<T1>();
	
	std::priority_queue<siv<T1>,std::vector<siv<T1>>,CompareSiv_vr<T1>> queue;
	
	for ( size_t np=0; np<refPts.size(); ++np ) {
		
		for ( size_t n=0; n<this->nodes.size(); ++n ) {
			queue.push( {n, m->l(this->nodes[n], refPts[np])} );
		}
		
		while ( !queue.empty() ) {
			siv<T1> s = queue.top();
			queue.pop();
			S[np].push_back( &(this->nodes[s.i]) );
		}
	}
	
	delete m;
}

template<typename T1, typename T2>
int Grid3Duifs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
								const std::vector<T1>& t0,
								const std::vector<sxyz<T1>>& Rx,
								std::vector<T1>& traveltimes,
								const size_t threadNo) const {
    
    if ( this->check_pts(Tx) == 1 ) return 1;
    if ( this->check_pts(Rx) == 1 ) return 1;
    
    for ( size_t n=0; n<this->nodes.size(); ++n ) {
        this->nodes[n].reinit( threadNo );
    }
    
	std::vector<bool> frozen( this->nodes.size(), false );
	initTx(Tx, t0, frozen, threadNo);
	
	T1 error = std::numeric_limits<T1>::max();
	std::vector<T1> times( this->nodes.size() );
	for ( size_t n=0; n<this->nodes.size(); ++n )
		times[n] = this->nodes[n].getTT( threadNo );
	
	int niter=0;
	while ( error >= epsilon && niter<nitermax ) {
		
		for ( size_t i=0; i<S.size(); ++i ) {
			
			// ascending
			for ( auto vertexC=S[i].begin(); vertexC!=S[i].end(); ++vertexC ) {
                if ( !frozen[(*vertexC)->getGridIndex()] )
					//                    this->local_3Dsolver(*vertexC, threadNo);
                    this->local_update3D(*vertexC, threadNo);
			}
			
			error = 0.0;
			for ( size_t n=0; n<this->nodes.size(); ++n ) {
				T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
				
				error += dt;
				times[n] = this->nodes[n].getTT(threadNo);
			}
			if ( error < epsilon ) {
                break;
            }
			
			// descending
			for ( auto vertexC=S[i].rbegin(); vertexC!=S[i].rend(); ++vertexC ) {
                if ( !frozen[(*vertexC)->getGridIndex()] )
					//                    this->local_3Dsolver(*vertexC, threadNo);
                    this->local_update3D(*vertexC, threadNo);
			}
			
			error = 0.0;
			for ( size_t n=0; n<this->nodes.size(); ++n ) {
				T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
				
				error += dt;
				times[n] = this->nodes[n].getTT(threadNo);
			}
			if ( error < epsilon ) {
				break;
            }
		}
        niter++;
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
int Grid3Duifs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
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
    
	std::vector<bool> frozen( this->nodes.size(), false );
	initTx(Tx, t0, frozen, threadNo);
    
	T1 error = std::numeric_limits<T1>::max();
	std::vector<T1> times( this->nodes.size() );
	for ( size_t n=0; n<this->nodes.size(); ++n )
		times[n] = this->nodes[n].getTT( threadNo );
	
    int niter=0;
	while ( error >= epsilon && niter<nitermax ) {
		
		for ( size_t i=0; i<S.size(); ++i ) {
			
			// ascending
			for ( auto vertexC=S[i].begin(); vertexC!=S[i].end(); ++vertexC ) {
                if ( !frozen[(*vertexC)->getGridIndex()] )
					//                    this->local_3Dsolver(*vertexC, threadNo);
                    this->local_update3D(*vertexC, threadNo);
			}
			
			//			char fname[200];
			//			sprintf(fname, "fsm%06d_%zd_a.dat",niter+1,i+1);
			//			saveTT(fname, threadNo);
			
			error = 0.0;
			for ( size_t n=0; n<this->nodes.size(); ++n ) {
				T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
				
				error += dt;
				times[n] = this->nodes[n].getTT(threadNo);
			}
			if ( error < epsilon ) {
                break;
            }
			
			// descending
			for ( auto vertexC=S[i].rbegin(); vertexC!=S[i].rend(); ++vertexC ) {
                if ( !frozen[(*vertexC)->getGridIndex()] )
					//                    this->local_3Dsolver(*vertexC, threadNo);
                    this->local_update3D(*vertexC, threadNo);
			}
			//			sprintf(fname, "fsm%06d_%zd_d.dat",niter+1,i+1);
			//			saveTT(fname, threadNo);
			
			error = 0.0;
			for ( size_t n=0; n<this->nodes.size(); ++n ) {
				T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
				
				error += dt;
				times[n] = this->nodes[n].getTT(threadNo);
			}
			if ( error < epsilon ) {
				break;
            }
            
		}
        niter++;
		std::cout << niter << std::endl;
	}
    std::cout << niter << " iterations were needed with epsilon = " << epsilon << '\n';
	
	
    if ( traveltimes.size() != Rx.size() ) {
        traveltimes.resize( Rx.size() );
    }
    
    for (size_t nr=0; nr<Rx.size(); ++nr) {
        traveltimes[nr]->resize( Rx[nr]->size() );
        for (size_t n=0; n<Rx[nr]->size(); ++n)
            (*traveltimes[nr])[n] = this->getTraveltime((*Rx[nr])[n], this->nodes, threadNo);
    }
	return 0;
}


template<typename T1, typename T2>
int Grid3Duifs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
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
    
	std::vector<bool> frozen( this->nodes.size(), false );
	initTx(Tx, t0, frozen, threadNo);
	
	T1 error = std::numeric_limits<T1>::max();
	std::vector<T1> times( this->nodes.size() );
	for ( size_t n=0; n<this->nodes.size(); ++n )
		times[n] = this->nodes[n].getTT( threadNo );
	
	int niter=0;
	while ( error >= epsilon && niter<nitermax ) {
		
		for ( size_t i=0; i<S.size(); ++i ) {
			
			// ascending
			for ( auto vertexC=S[i].begin(); vertexC!=S[i].end(); ++vertexC ) {
                if ( !frozen[(*vertexC)->getGridIndex()] )
					//                    this->local_3Dsolver(*vertexC, threadNo);
                    this->local_update3D(*vertexC, threadNo);
			}
			
			error = 0.0;
			for ( size_t n=0; n<this->nodes.size(); ++n ) {
				T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
				
				error += dt;
				times[n] = this->nodes[n].getTT(threadNo);
			}
			if ( error < epsilon ) {
                break;
            }
			
			// descending
			for ( auto vertexC=S[i].rbegin(); vertexC!=S[i].rend(); ++vertexC ) {
                if ( !frozen[(*vertexC)->getGridIndex()] )
					//                    this->local_3Dsolver(*vertexC, threadNo);
                    this->local_update3D(*vertexC, threadNo);
			}
			
			error = 0.0;
			for ( size_t n=0; n<this->nodes.size(); ++n ) {
				T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
				
				error += dt;
				times[n] = this->nodes[n].getTT(threadNo);
			}
			if ( error < epsilon ) {
				break;
            }
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
    for ( size_t ni=0; ni<r_data.size(); ++ni ) {
        r_data[ni].resize( 0 );
    }
    
    if ( rp_ho ) {
        for (size_t n=0; n<Rx.size(); ++n) {
            traveltimes[n] = this->getTraveltime(Rx[n], this->nodes, threadNo);
            this->getRaypath_ho(Tx, Rx[n], traveltimes[n], r_data[n], threadNo);
        }
    } else {
        for (size_t n=0; n<Rx.size(); ++n) {
            traveltimes[n] = this->getTraveltime(Rx[n], this->nodes, threadNo);
            this->getRaypath(Tx, Rx[n], traveltimes[n], r_data[n], threadNo);
        }
    }
	
	return 0;
}

template<typename T1, typename T2>
int Grid3Duifs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
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
    
	std::vector<bool> frozen( this->nodes.size(), false );
	initTx(Tx, t0, frozen, threadNo);
    
	T1 error = std::numeric_limits<T1>::max();
	std::vector<T1> times( this->nodes.size() );
	for ( size_t n=0; n<this->nodes.size(); ++n )
		times[n] = this->nodes[n].getTT( threadNo );
	
    int niter=0;
	while ( error >= epsilon && niter<nitermax ) {
		
		for ( size_t i=0; i<S.size(); ++i ) {
			
			// ascending
			for ( auto vertexC=S[i].begin(); vertexC!=S[i].end(); ++vertexC ) {
                if ( !frozen[(*vertexC)->getGridIndex()] )
					//                    this->local_3Dsolver(*vertexC, threadNo);
                    this->local_update3D(*vertexC, threadNo);
			}
			
			//			char fname[200];
			//			sprintf(fname, "fsm%06d_%zd_a.dat",niter+1,i+1);
			//			saveTT(fname, threadNo);
			
			error = 0.0;
			for ( size_t n=0; n<this->nodes.size(); ++n ) {
				T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
				
				error += dt;
				times[n] = this->nodes[n].getTT(threadNo);
			}
			if ( error < epsilon ) {
                break;
            }
			
			// descending
			for ( auto vertexC=S[i].rbegin(); vertexC!=S[i].rend(); ++vertexC ) {
                if ( !frozen[(*vertexC)->getGridIndex()] )
					//                    this->local_3Dsolver(*vertexC, threadNo);
                    this->local_update3D(*vertexC, threadNo);
			}
			//			sprintf(fname, "fsm%06d_%zd_d.dat",niter+1,i+1);
			//			saveTT(fname, threadNo);
			
			error = 0.0;
			for ( size_t n=0; n<this->nodes.size(); ++n ) {
				T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
				
				error += dt;
				times[n] = this->nodes[n].getTT(threadNo);
			}
			if ( error < epsilon ) {
				break;
            }
            
		}
        niter++;
		std::cout << niter << std::endl;
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
		
        if ( rp_ho ) {
            for (size_t n=0; n<Rx[nr]->size(); ++n) {
                (*traveltimes[nr])[n] = this->getTraveltime((*Rx[nr])[n], this->nodes, threadNo);
                this->getRaypath_ho(Tx, (*Rx[nr])[n], (*traveltimes[nr])[n], (*r_data[nr])[n], threadNo);
            }
        } else {
            for (size_t n=0; n<Rx[nr]->size(); ++n) {
                (*traveltimes[nr])[n] = this->getTraveltime((*Rx[nr])[n], this->nodes, threadNo);
                this->getRaypath(Tx, (*Rx[nr])[n], (*traveltimes[nr])[n], (*r_data[nr])[n], threadNo);
            }
        }
    }
	return 0;
}


template<typename T1, typename T2>
void Grid3Duifs<T1,T2>::initTx(const std::vector<sxyz<T1>>& Tx,
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
				
                if ( Grid3Dui<T1,T2,Node3Di<T1,T2>>::source_radius == 0.0 ) {
                    // populate around Tx
                    for ( size_t no=0; no<this->nodes[nn].getOwners().size(); ++no ) {
                    
                        T2 cellNo = this->nodes[nn].getOwners()[no];
                        for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                            T2 neibNo = this->neighbors[cellNo][k];
                            if ( neibNo == nn ) continue;
                            T1 dt = this->computeDt(this->nodes[nn], this->nodes[neibNo]);
                            
                            if ( t0[n]+dt < this->nodes[neibNo].getTT(threadNo) ) {
                                this->nodes[neibNo].setTT( t0[n]+dt, threadNo );
                                frozen[neibNo] = true;
                            }
                        }
                    }
                } else {
                    // find nodes within source radius
                    size_t nodes_added = 0;
                    for ( size_t no=0; no<this->nodes.size(); ++no ) {
                        
                        if ( no == nn ) continue;
                        
                        T1 d = this->nodes[nn].getDistance( this->nodes[no] );
                        if ( d <= Grid3Dui<T1,T2,Node3Di<T1,T2>>::source_radius ) {
                            
                            T1 dt = this->computeDt(this->nodes[nn], this->nodes[no] );
                            
                            if ( t0[n]+dt < this->nodes[no].getTT(threadNo) ) {
                                if ( this->nodes[no].getTT(threadNo) == std::numeric_limits<T1>::max() ) nodes_added++;
                                this->nodes[no].setTT( t0[n]+dt, threadNo );
                            }
                        }
                    }
                    if ( nodes_added == 0 ) {
                        std::cerr << "Error: no nodes found within source radius, aborting" << std::endl;
                        abort();
                    } else {
                        std::cout << "(found " << nodes_added << " nodes arounf Tx point)\n";
                    }
                }
				
                break;
            }
        }
        if ( found==false ) {
			
			T2 cellNo = this->getCellNo(Tx[n]);
            if ( Grid3Dui<T1,T2,Node3Di<T1,T2>>::source_radius == 0.0 ) {
                for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                    T2 neibNo = this->neighbors[cellNo][k];
				
                    // compute dt
                    T1 dt = this->nodes[neibNo].getDistance(Tx[n])*this->nodes[neibNo].getNodeSlowness();
                    
                    this->nodes[neibNo].setTT( t0[n]+dt, threadNo );
                    frozen[neibNo] = true;
                }
            } else if ( Tx.size()==1 ) { // look into source radius only for point sources
                // find nodes within source radius
                size_t nodes_added = 0;
                for ( size_t no=0; no<this->nodes.size(); ++no ) {
                    
                    T1 d = this->nodes[no].getDistance( Tx[n] );
                    if ( d <= Grid3Dui<T1,T2,Node3Di<T1,T2>>::source_radius ) {
                        
                        T1 dt = this->nodes[no].getDistance(Tx[n])*this->nodes[no].getNodeSlowness();
                        
                        if ( t0[n]+dt < this->nodes[no].getTT(threadNo) ) {
                            if ( this->nodes[no].getTT(threadNo) == std::numeric_limits<T1>::max() ) nodes_added++;
                            this->nodes[no].setTT( t0[n]+dt, threadNo );
                        }
                    }
                }
                if ( nodes_added == 0 ) {
                    std::cerr << "Error: no nodes found within source radius, aborting" << std::endl;
                    abort();
                } else {
                    std::cout << "(found " << nodes_added << " nodes around Tx point)\n";
                }
            }
		}
    }
}


#endif

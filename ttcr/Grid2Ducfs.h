//
//  Grid2Ducfs.h
//  ttcr_u
//
//  Created by Bernard Giroux on 2014-01-30.
//  Copyright (c) 2014 Bernard Giroux. All rights reserved.
//


//
//@ARTICLE{qian07,
//	author = {Qian, Jianliang and Zhang, Yong-Tao and Zhao, Hong-Kai},
//	title = {Fast Sweeping Methods for Eikonal Equations on Triangular Meshes},
//	journal = {SIAM Journal on Numerical Analysis},
//	year = {2007},
//	volume = {45},
//	pages = {83--107},
//	number = {1},
//	doi = {10.1137/050627083},
//	issn = {00361429},
//	publisher = {Society for Industrial and Applied Mathematics},
//	url = {http://www.jstor.org/stable/40232919}
//	}
//




#ifndef ttcr_u_Grid2Ducfs_h
#define ttcr_u_Grid2Ducfs_h

#include <fstream>
#include <queue>

#include "Grid2Duc.h"
#include "Node2Dc.h"
#include "Metric.h"

template<typename T1, typename T2>
class Grid2Ducfs : public Grid2Duc<T1,T2,Node2Dc<T1,T2>,sxz<T1>> {
public:
	Grid2Ducfs(const std::vector<sxz<T1>>& no,
				const std::vector<triangleElem<T2>>& tri,
				const T1 eps, const int maxit, const size_t nt=1) :
	Grid2Duc<T1,T2,Node2Dc<T1,T2>,sxz<T1>>(no, tri, nt),
	epsilon(eps), nitermax(maxit), S()
	{
		buildGridNodes(no, nt);
		this->buildGridNeighbors();
		this->processObtuse();
	}
	
	~Grid2Ducfs() {
	}
	
	void initOrdering(const std::vector<sxz<T1>>& refPts, const int order);
	
	int raytrace(const std::vector<sxz<T1>>&,
				 const std::vector<T1>&,
				 const std::vector<sxz<T1>>&,
				 std::vector<T1>&,
				 const size_t=0) const;
    
    int raytrace(const std::vector<sxz<T1>>&,
                 const std::vector<T1>&,
                 const std::vector<const std::vector<sxz<T1>>*>&,
                 std::vector<std::vector<T1>*>&,
                 const size_t=0) const;
	
	int raytrace(const std::vector<sxz<T1>>&,
                 const std::vector<T1>& ,
                 const std::vector<sxz<T1>>&,
                 std::vector<T1>&,
                 std::vector<std::vector<sxz<T1>>>&,
				 const size_t=0) const;
    
    int raytrace(const std::vector<sxz<T1>>&,
                 const std::vector<T1>&,
                 const std::vector<const std::vector<sxz<T1>>*>&,
                 std::vector<std::vector<T1>*>&,
                 std::vector<std::vector<std::vector<sxz<T1>>>*>&,
                 const size_t=0) const;
	
private:
	T1 epsilon;
	int nitermax;
	std::vector<std::vector<Node2Dc<T1,T2>*>> S;
	
	void buildGridNodes(const std::vector<sxz<T1>>&,
						const size_t);
	
	void initTx(const std::vector<sxz<T1>>& Tx, const std::vector<T1>& t0,
				std::vector<bool>& frozen, const size_t threadNo) const;
	
    
};

template<typename T1, typename T2>
void Grid2Ducfs<T1,T2>::buildGridNodes(const std::vector<sxz<T1>>& no,
									   const size_t nt) {
	
	// primary nodes
	for ( T2 n=0; n<no.size(); ++n ) {
		this->nodes[n].setXZindex( no[n].x, no[n].z, n );
	}
	for ( T2 ntri=0; ntri<this->triangles.size(); ++ntri ) {
		
		for ( size_t nl=0; nl<3; ++nl ) {
			
			// push owner for primary nodes
			this->nodes[ this->triangles[ntri].i[nl] ].pushOwner( ntri );
		}
		
		// distance between node 1 & 2 (opposite of node 0)
		T1 a = this->nodes[ this->triangles[ntri].i[1] ].getDistance( this->nodes[ this->triangles[ntri].i[2] ] );
		
		// distance between node 0 & 2 (opposite of node 1)
		T1 b = this->nodes[ this->triangles[ntri].i[0] ].getDistance( this->nodes[ this->triangles[ntri].i[2] ] );
		
		// distance between node 0 & 1 (opposite of node 2]
		T1 c = this->nodes[ this->triangles[ntri].i[0] ].getDistance( this->nodes[ this->triangles[ntri].i[1] ] );
		
		this->triangles[ntri].l[0] = a;
		this->triangles[ntri].l[1] = b;
		this->triangles[ntri].l[2] = c;
		
		// angle at node 0
		this->triangles[ntri].a[0] = acos((b*b + c*c - a*a)/(2.*b*c));
		
		// angle at node 1
		this->triangles[ntri].a[1] = acos((c*c + a*a - b*b)/(2.*a*c));
		
		// angle at node 2
		this->triangles[ntri].a[2] = acos((a*a + b*b - c*c)/(2.*a*b));
		
	}
}


template<typename T1, typename T2>
void Grid2Ducfs<T1,T2>::initOrdering(const std::vector<sxz<T1>>& refPts,
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
int Grid2Ducfs<T1,T2>::raytrace(const std::vector<sxz<T1>>& Tx,
							  const std::vector<T1>& t0,
							  const std::vector<sxz<T1>>& Rx,
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
                    this->local_solver(*vertexC, threadNo);
			}
			
			error = 0.0;
			for ( size_t n=0; n<this->nodes.size(); ++n ) {
				T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
				
				error += dt;
				times[n] = this->nodes[n].getTT(threadNo);
			}
			if ( error < epsilon ) {
				niter++;
                break;
            }
			
			// descending
			for ( auto vertexC=S[i].rbegin(); vertexC!=S[i].rend(); ++vertexC ) {
                if ( !frozen[(*vertexC)->getGridIndex()] )
                    this->local_solver(*vertexC, threadNo);
			}
			
			error = 0.0;
			for ( size_t n=0; n<this->nodes.size(); ++n ) {
				T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
				
				error += dt;
				times[n] = this->nodes[n].getTT(threadNo);
			}
			if ( error < epsilon ) {
                niter++;
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
int Grid2Ducfs<T1,T2>::raytrace(const std::vector<sxz<T1>>& Tx,
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
                    this->local_solver(*vertexC, threadNo);
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
				niter++;
                break;
            }
			
			// descending
			for ( auto vertexC=S[i].rbegin(); vertexC!=S[i].rend(); ++vertexC ) {
                if ( !frozen[(*vertexC)->getGridIndex()] )
                    this->local_solver(*vertexC, threadNo);
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
                niter++;
				break;
            }
            
		}
        niter++;
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
int Grid2Ducfs<T1,T2>::raytrace(const std::vector<sxz<T1>>& Tx,
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
                    this->local_solver(*vertexC, threadNo);
			}
			
			error = 0.0;
			for ( size_t n=0; n<this->nodes.size(); ++n ) {
				T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
				
				error += dt;
				times[n] = this->nodes[n].getTT(threadNo);
			}
			if ( error < epsilon ) {
				niter++;
                break;
            }
			
			// descending
			for ( auto vertexC=S[i].rbegin(); vertexC!=S[i].rend(); ++vertexC ) {
                if ( !frozen[(*vertexC)->getGridIndex()] )
                    this->local_solver(*vertexC, threadNo);
			}
			
			error = 0.0;
			for ( size_t n=0; n<this->nodes.size(); ++n ) {
				T1 dt = fabs( times[n] - this->nodes[n].getTT(threadNo) );
				
				error += dt;
				times[n] = this->nodes[n].getTT(threadNo);
			}
			if ( error < epsilon ) {
                niter++;
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

    for (size_t n=0; n<Rx.size(); ++n) {
        traveltimes[n] = this->getTraveltime(Rx[n], this->nodes, threadNo);
        this->getRaypath(Tx, Rx[n], traveltimes[n], r_data[n], threadNo);
    }
	
	return 0;
}


template<typename T1, typename T2>
int Grid2Ducfs<T1,T2>::raytrace(const std::vector<sxz<T1>>& Tx,
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
                    this->local_solver(*vertexC, threadNo);
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
				niter++;
                break;
            }
			
			// descending
			for ( auto vertexC=S[i].rbegin(); vertexC!=S[i].rend(); ++vertexC ) {
                if ( !frozen[(*vertexC)->getGridIndex()] )
                    this->local_solver(*vertexC, threadNo);
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
                niter++;
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

    for (size_t nr=0; nr<Rx.size(); ++nr) {
        traveltimes[nr]->resize( Rx[nr]->size() );
        r_data[nr]->resize( Rx[nr]->size() );
        for ( size_t ni=0; ni<r_data[nr]->size(); ++ni ) {
            (*r_data[nr])[ni].resize( 0 );
        }

        for (size_t n=0; n<Rx[nr]->size(); ++n) {
            (*traveltimes[nr])[n] = this->getTraveltime((*Rx[nr])[n], this->nodes, threadNo);
			
            this->getRaypath(Tx, (*Rx[nr])[n], (*traveltimes[nr])[n], (*r_data[nr])[n], threadNo);
		}
    }
	return 0;
}


template<typename T1, typename T2>
void Grid2Ducfs<T1,T2>::initTx(const std::vector<sxz<T1>>& Tx,
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
				
                // populate around Tx
                for ( size_t no=0; no<this->nodes[nn].getOwners().size(); ++no ) {
                    
                    T2 cellNo = this->nodes[nn].getOwners()[no];
                    for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                        T2 neibNo = this->neighbors[cellNo][k];
                        if ( neibNo == nn ) continue;
                        T1 dt = this->computeDt(this->nodes[nn], this->nodes[neibNo], cellNo);
							
                        if ( t0[n]+dt < this->nodes[neibNo].getTT(threadNo) ) {
                            this->nodes[neibNo].setTT( t0[n]+dt, threadNo );
//                            this->nodes[neibNo].setnodeParent(this->nodes[nn].getGridIndex(),threadNo);
//                            this->nodes[neibNo].setCellParent(cellNo, threadNo );
                            //frozen[neibNo] = true;
                        }
					}
				}
				
                break;
            }
        }
        if ( found==false ) {
			
			T2 cellNo = this->getCellNo(Tx[n]);
			for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                T2 neibNo = this->neighbors[cellNo][k];
				
				// compute dt
                T1 dt = this->computeDt(this->nodes[neibNo], Tx[n], cellNo);
				
				this->nodes[neibNo].setTT( t0[n]+dt, threadNo );
                frozen[neibNo] = true;
				
			}
		}
    }
}



#endif

//
//  Grid2Ducsp.h
//  ttcr_u
//
//  Created by Bernard Giroux on 2012-09-20.
//  Copyright (c) 2012 Bernard Giroux. All rights reserved.
//




#ifndef __GRID2DUCSP_H__
#define __GRID2DUCSP_H__

#include <array>
#include <fstream>
#include <queue>

#include "Grid2Duc.h"
#include "Node2Dcsp.h"

template<typename T1, typename T2>
class Grid2Ducsp : public Grid2Duc<T1,T2,Node2Dcsp<T1,T2>> {
public:
	Grid2Ducsp(const std::vector<sxz<T1>>& no,
			   const std::vector<triangleElem<T2> >& tri,
			   const T2 ns, const size_t nt=1) :
	Grid2Duc<T1,T2,Node2Dcsp<T1,T2>>(no, tri, nt)
	{
		buildGridNodes(no, ns, nt);
		this->buildGridNeighbors();
	}
	
	~Grid2Ducsp() {
	}
	
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
	
    int raytrace(const std::vector<sxz<T1> >&,
                 const std::vector<T1>& ,
                 const std::vector<sxz<T1> >&,
                 std::vector<T1>&,
                 std::vector<std::vector<sxz<T1> > >&,
                 std::vector<std::vector<siv<T1> > >&,
				 const size_t=0) const;
	
    
//	int getTraveltimes(const std::vector<sxz<T1> >&,
//					   std::vector<T1>&,
//					   const size_t=0) const;
//	
//	int getTraveltimes(const std::vector<sxz<T1>>&,
//					   const std::vector<sxz<T1>>&,
//					   std::vector<T1>&,
//					   std::vector<std::vector<sxz<T1>>>&,
//					   const size_t=0) const;
	
private:
	void buildGridNodes(const std::vector<sxz<T1>>&,
						const int,
						const size_t);
	
	void initQueue(const std::vector<sxz<T1> >& Tx,
				   const std::vector<T1>& t0,
				   std::priority_queue<Node2Dcsp<T1,T2>*,
				   std::vector<Node2Dcsp<T1,T2>*>,
				   CompareNodePtr<T1> >& queue,
				   std::vector<Node2Dcsp<T1,T2> >& txNodes,
				   std::vector<bool>& inQueue,
				   std::vector<bool>& frozen,
				   const size_t threadNo) const;
	
	void propagate(std::priority_queue<Node2Dcsp<T1,T2>*,
				   std::vector<Node2Dcsp<T1,T2>*>,
				   CompareNodePtr<T1> >& queue,
				   std::vector<bool>& inQueue,
				   std::vector<bool>& frozen,
				   const size_t threadNo) const;
	
};

template<typename T1, typename T2>
void Grid2Ducsp<T1,T2>::buildGridNodes(const std::vector<sxz<T1>>& no,
									   const int nsecondary,
									   const size_t nt) {
	
	// primary nodes
	for ( T2 n=0; n<no.size(); ++n ) {
		this->nodes[n].setXZindex( no[n].x, no[n].z, n );
	}
	T2 nNodes = static_cast<T2>(this->nodes.size());
	
	std::map<std::array<T2,2>,std::vector<T2>> lineMap;
	std::array<T2,2> lineKey;
	typename std::map<std::array<T2,2>,std::vector<T2>>::iterator lineIt;
    
    size_t estLineNo = (this->triangles.size()+this->triangles.size()/10) * 3/2;
	this->nodes.reserve( nNodes + estLineNo*nsecondary );
	
	// edge nodes
	Node2Dcsp<T1,T2> tmpNode(nt);
	for ( T2 ntri=0; ntri<this->triangles.size(); ++ntri ) {
				
		for ( size_t nl=0; nl<3; ++nl ) {
		
			// push owner for primary nodes
			this->nodes[ this->triangles[ntri].i[nl] ].pushOwner( ntri );
			
			if ( nsecondary>0 ) {
				
				lineKey = { this->triangles[ntri].i[nl],
					this->triangles[ntri].i[(nl+1)%3] };
				std::sort(lineKey.begin(), lineKey.end());
				
				lineIt = lineMap.find( lineKey );
				if ( lineIt == lineMap.end() ) {
					// not found, insert new pair
					lineMap[ lineKey ] = std::vector<T2>(nsecondary);
				} else {
					for ( size_t n=0; n<lineIt->second.size(); ++n ) {
						// setting owners
						this->nodes[ lineIt->second[n] ].pushOwner( ntri );
					}
					continue;
				}
				
				sxz<T1> d = { (no[lineKey[1]].x-no[lineKey[0]].x)/(nsecondary+1),
					          (no[lineKey[1]].z-no[lineKey[0]].z)/(nsecondary+1) };
			
				for ( size_t n2=0; n2<nsecondary; ++n2 ) {
					tmpNode.setXZindex(no[lineKey[0]].x+(1+n2)*d.x,
									   no[lineKey[0]].z+(1+n2)*d.z,
									   nNodes );
					lineMap[lineKey][n2] = nNodes++;
					this->nodes.push_back( tmpNode );
					this->nodes.back().pushOwner( ntri );
				}
			}
		}
	}
	
	this->nodes.shrink_to_fit();
}

template<typename T1, typename T2>
int Grid2Ducsp<T1,T2>::raytrace(const std::vector<sxz<T1> >& Tx,
								const std::vector<T1>& t0,
								const std::vector<sxz<T1> >& Rx,
								std::vector<T1>& traveltimes,
								const size_t threadNo) const {
    
    if ( this->check_pts(Tx) == 1 ) return 1;
    if ( this->check_pts(Rx) == 1 ) return 1;
    
    for ( size_t n=0; n<this->nodes.size(); ++n ) {
        this->nodes[n].reinit( threadNo );
    }
    
    CompareNodePtr<T1> cmp(threadNo);
    std::priority_queue< Node2Dcsp<T1,T2>*, std::vector<Node2Dcsp<T1,T2>*>,
    CompareNodePtr<T1> > queue( cmp );
    
    std::vector<Node2Dcsp<T1,T2> > txNodes;
    std::vector<bool> inQueue( this->nodes.size(), false );
    std::vector<bool> frozen( this->nodes.size(), false );
    
    initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
    
    propagate(queue, inQueue, frozen, threadNo);
    
    if ( traveltimes.size() != Rx.size() ) {
        traveltimes.resize( Rx.size() );
    }
    
    for (size_t n=0; n<Rx.size(); ++n) {
        traveltimes[n] = this->getTraveltime(Rx[n], this->nodes, threadNo);
    }
    return 0;
}

template<typename T1, typename T2>
int Grid2Ducsp<T1,T2>::raytrace(const std::vector<sxz<T1> >& Tx,
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
    
    CompareNodePtr<T1> cmp(threadNo);
    std::priority_queue< Node2Dcsp<T1,T2>*, std::vector<Node2Dcsp<T1,T2>*>,
    CompareNodePtr<T1> > queue( cmp );
    
    std::vector<Node2Dcsp<T1,T2> > txNodes;
    std::vector<bool> inQueue( this->nodes.size(), false );
    std::vector<bool> frozen( this->nodes.size(), false );
    
    initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
    
    propagate(queue, inQueue, frozen, threadNo);
    
    if ( traveltimes.size() != Rx.size() ) {
        traveltimes.resize( Rx.size() );
    }
    
    for (size_t nr=0; nr<Rx.size(); ++nr) {
        traveltimes[nr]->resize( Rx[nr]->size() );
        for (size_t n=0; n<Rx[nr]->size(); ++n)
            (*traveltimes[nr])[n] = this->getTraveltime((*Rx[nr])[n],
														this->nodes, threadNo);
    }
	return 0;
}

template<typename T1, typename T2>
int Grid2Ducsp<T1,T2>::raytrace(const std::vector<sxz<T1>>& Tx,
								const std::vector<T1>& t0,
								const std::vector<sxz<T1>>& Rx,
								std::vector<T1>& traveltimes,
								std::vector<std::vector<sxz<T1>>>& r_data,
								const size_t threadNo) const {
    
	//	std::cout << "   running in thread no " << threadNo << std::endl;
    if ( this->check_pts(Tx) == 1 ) return 1;
    if ( this->check_pts(Rx) == 1 ) return 1;
    
    for ( size_t n=0; n<this->nodes.size(); ++n ) {
        this->nodes[n].reinit( threadNo );
    }
    
    CompareNodePtr<T1> cmp(threadNo);
    std::priority_queue< Node2Dcsp<T1,T2>*, std::vector<Node2Dcsp<T1,T2>*>,
    CompareNodePtr<T1> > queue( cmp );
    
    std::vector<Node2Dcsp<T1,T2> > txNodes;
    std::vector<bool> inQueue( this->nodes.size(), false );
    std::vector<bool> frozen( this->nodes.size(), false );
    
    initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
    
    propagate(queue, inQueue, frozen, threadNo);
    
    if ( traveltimes.size() != Rx.size() ) {
        traveltimes.resize( Rx.size() );
    }
    if ( r_data.size() != Rx.size() ) {
        r_data.resize( Rx.size() );
    }
    for ( size_t ni=0; ni<r_data.size(); ++ni ) {
        r_data[ni].resize( 0 );
    }
    T2 nodeParentRx;
    T2 cellParentRx;
    
    for (size_t n=0; n<Rx.size(); ++n) {
        
        traveltimes[n] = this->getTraveltime(Rx[n], this->nodes, nodeParentRx,
											 cellParentRx, threadNo);
        
        bool flag=false;
        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx[n] == Tx[ns] ) {
                
                r_data[n].resize( 1 );
                r_data[n][0].x = Rx[n].x;
                r_data[n][0].z = Rx[n].z;
                
                flag = true;
                break;
            }
        }
        if ( flag ) continue;
        
        // Rx are in nodes (not txNodes)
        std::vector<Node2Dcsp<T1,T2> > *node_p;
        node_p = &(this->nodes);
        
        std::vector<sxz<T1> > r_tmp;
        T2 iChild, iParent = nodeParentRx;
        sxz<T1> child;
		
        // store the son's coord
        child.x = Rx[n].x;
        child.z = Rx[n].z;
        while ( (*node_p)[iParent].getNodeParent(threadNo) !=
               std::numeric_limits<T2>::max() ) {
 			
			r_tmp.push_back( child );
						
			// we now go up in time - parent becomes the child of grand'pa
			iChild = iParent;
			child.x = (*node_p)[iChild].getX();
			child.z = (*node_p)[iChild].getZ();
            
			// grand'pa is now papa
			iParent = (*node_p)[iChild].getNodeParent(threadNo);
            if ( iParent >= this->nodes.size() ) {
                node_p = &txNodes;
                iParent -= this->nodes.size();
            }
            else {
                node_p = &(this->nodes);
            }
		}
		
		// parent is now at Tx
		r_tmp.push_back( child );
				
		// finally, store Tx position
		child.x = (*node_p)[iParent].getX();
		child.z = (*node_p)[iParent].getZ();
		r_tmp.push_back( child );
		
        // the order should be from Tx to Rx, so we reorder...
        iParent = static_cast<T2>(r_tmp.size());
        r_data[n].resize( r_tmp.size() );
        for ( size_t nn=0; nn<r_data[n].size(); ++nn ) {
            r_data[n][nn].x = r_tmp[ iParent-1-nn ].x;
            r_data[n][nn].z = r_tmp[ iParent-1-nn ].z;
        }
    }
	return 0;
}

template<typename T1, typename T2>
int Grid2Ducsp<T1,T2>::raytrace(const std::vector<sxz<T1> >& Tx,
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
    
    CompareNodePtr<T1> cmp(threadNo);
    std::priority_queue< Node2Dcsp<T1,T2>*, std::vector<Node2Dcsp<T1,T2>*>,
    CompareNodePtr<T1> > queue( cmp );
    
    std::vector<Node2Dcsp<T1,T2> > txNodes;
    std::vector<bool> inQueue( this->nodes.size(), false );
    std::vector<bool> frozen( this->nodes.size(), false );
    
    initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
    
    propagate(queue, inQueue, frozen, threadNo);
    
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
        
        T2 nodeParentRx;
        T2 cellParentRx;
        
        for (size_t n=0; n<Rx[nr]->size(); ++n) {
            
            (*traveltimes[nr])[n] = this->getTraveltime((*Rx[nr])[n], this->nodes,
                                                  nodeParentRx, cellParentRx,
                                                  threadNo);
            
            bool flag=false;
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                if ( (*Rx[nr])[n] == Tx[ns] ) {
                    
                    (*r_data[nr])[n].resize( 1 );
                    (*r_data[nr])[n][0].x = (*Rx[nr])[n].x;
                    (*r_data[nr])[n][0].z = (*Rx[nr])[n].z;
                    
                    flag = true;
                    break;
                }
            }
            if ( flag ) continue;
            
            // Rx are in nodes (not txNodes)
            std::vector<Node2Dcsp<T1,T2> > *node_p;
            node_p = &(this->nodes);
            
            std::vector<sxz<T1> > r_tmp;
            T2 iChild, iParent = nodeParentRx;
            sxz<T1> child;
            
            // store the son's coord
            child.x = (*Rx[nr])[n].x;
            child.z = (*Rx[nr])[n].z;
            while ( (*node_p)[iParent].getNodeParent(threadNo) !=
                   std::numeric_limits<T2>::max() ) {
                
                r_tmp.push_back( child );
                
                // we now go up in time - parent becomes the child of grand'pa
                iChild = iParent;
                child.x = (*node_p)[iChild].getX();
                child.z = (*node_p)[iChild].getZ();
                
                // grand'pa is now papa
                iParent = (*node_p)[iChild].getNodeParent(threadNo);
                if ( iParent >= this->nodes.size() ) {
                    node_p = &txNodes;
                    iParent -= this->nodes.size();
                }
                else {
                    node_p = &(this->nodes);
                }
            }
            
            // parent is now at Tx
            r_tmp.push_back( child );
            
            // finally, store Tx position
            child.x = (*node_p)[iParent].getX();
            child.z = (*node_p)[iParent].getZ();
            r_tmp.push_back( child );
            
            // the order should be from Tx to Rx, so we reorder...
            iParent = static_cast<T2>(r_tmp.size());
            (*r_data[nr])[n].resize( r_tmp.size() );
            for ( size_t nn=0; nn<(*r_data[nr])[n].size(); ++nn ) {
                (*r_data[nr])[n][nn].x = r_tmp[ iParent-1-nn ].x;
                (*r_data[nr])[n][nn].z = r_tmp[ iParent-1-nn ].z;
            }
            
        }
    }
	return 0;
}

template<typename T1, typename T2>
int Grid2Ducsp<T1,T2>::raytrace(const std::vector<sxz<T1> >& Tx,
								const std::vector<T1>& t0,
								const std::vector<sxz<T1> >& Rx,
								std::vector<T1>& traveltimes,
								std::vector<std::vector<sxz<T1> > >& r_data,
								std::vector<std::vector<siv<T1> > >& l_data,
								const size_t threadNo) const {
    
	//	std::cout << "   running in thread no " << threadNo << std::endl;
    if ( this->check_pts(Tx) == 1 ) return 1;
    if ( this->check_pts(Rx) == 1 ) return 1;
    
    for ( size_t n=0; n<this->nodes.size(); ++n ) {
        this->nodes[n].reinit( threadNo );
    }
    
	CompareNodePtr<T1> cmp(threadNo);
    std::priority_queue< Node2Dcsp<T1,T2>*, std::vector<Node2Dcsp<T1,T2>*>,
    CompareNodePtr<T1> > queue( cmp );
    
    std::vector<Node2Dcsp<T1,T2> > txNodes;
    std::vector<bool> inQueue( this->nodes.size(), false );
    std::vector<bool> frozen( this->nodes.size(), false );
    
    initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
    
    propagate(queue, inQueue, frozen, threadNo);
    
    if ( traveltimes.size() != Rx.size() ) {
        traveltimes.resize( Rx.size() );
    }
    if ( l_data.size() != Rx.size() ) {
        l_data.resize( Rx.size() );
    }
    for ( size_t ni=0; ni<l_data.size(); ++ni ) {
        l_data[ni].resize( 0 );
    }
    if ( r_data.size() != Rx.size() ) {
        r_data.resize( Rx.size() );
    }
    for ( size_t ni=0; ni<r_data.size(); ++ni ) {
        r_data[ni].resize( 0 );
    }
    T2 nodeParentRx;
    T2 cellParentRx;
    
    for (size_t n=0; n<Rx.size(); ++n) {
        
        traveltimes[n] = this->getTraveltime(Rx[n], this->nodes, nodeParentRx, cellParentRx,
									   threadNo);
        
        bool flag=false;
        for ( size_t ns=0; ns<Tx.size(); ++n ) {
            if ( Rx[n] == Tx[ns] ) {
                
                r_data[n].resize( 1 );
                r_data[n][0].x = Rx[n].x;
                r_data[n][0].z = Rx[n].z;
                
                // no need to update l_data: ray length is zero
                
                flag = true;
            }
        }
        if ( flag ) continue;

        // Rx are in nodes (not txNodes)
        std::vector<Node2Dcsp<T1,T2> > *node_p;
        node_p = &(this->nodes);
        
        std::vector<sxz<T1> > r_tmp;
        T2 iChild, iParent = nodeParentRx;
        sxz<T1> child;
        siv<T1> cell;
		
        // store the son's coord
        child.x = Rx[n].x;
        child.z = Rx[n].z;
        cell.i = cellParentRx;
        while ( (*node_p)[iParent].getNodeParent(threadNo) !=
               std::numeric_limits<T2>::max() ) {
 			
			r_tmp.push_back( child );
			
			cell.v = (*node_p)[iParent].getDistance( child );
			bool found=false;
			for (size_t nc=0; nc<l_data[n].size(); ++nc) {
				if ( l_data[n][nc].i == cell.i ) {
					l_data[n][nc].v += cell.v;
					found = true;
				}
			}
			if ( found == false ) {
				l_data[n].push_back( cell );
			}
			
			// we now go up in time - parent becomes the child of grand'pa
			iChild = iParent;
			child.x = (*node_p)[iChild].getX();
			child.z = (*node_p)[iChild].getZ();
            cell.i = (*node_p)[iChild].getCellParent(threadNo);
			
			// grand'pa is now papa
			iParent = (*node_p)[iChild].getNodeParent(threadNo);
            if ( iParent >= this->nodes.size() ) {
                node_p = &txNodes;
                iParent -= this->nodes.size();
            }
            else {
                node_p = &(this->nodes);
            }
		}
		
		// parent is now at Tx
		r_tmp.push_back( child );
		
		cell.v = (*node_p)[iParent].getDistance( child );
		bool found=false;
		for (size_t nc=0; nc<l_data[n].size(); ++nc) {
			if ( l_data[n][nc].i == cell.i ) {
				l_data[n][nc].v += cell.v;
				found = true;
			}
		}
		if ( found == false ) {
			l_data[n].push_back( cell );
		}
		
		// finally, store Tx position
		child.x = (*node_p)[iParent].getX();
		child.z = (*node_p)[iParent].getZ();
		r_tmp.push_back( child );
		
        //  must be sorted to build matrix L
        sort(l_data[n].begin(), l_data[n].end(), CompareSiv_i<T1>());
        
        // the order should be from Tx to Rx, so we reorder...
        iParent = static_cast<T2>(r_tmp.size());
        r_data[n].resize( r_tmp.size() );
        for ( size_t nn=0; nn<r_data[n].size(); ++nn ) {
            r_data[n][nn].x = r_tmp[ iParent-1-nn ].x;
            r_data[n][nn].z = r_tmp[ iParent-1-nn ].z;
        }
    }
	return 0;
}


//template<typename T1, typename T2>
//int Grid2Ducsp<T1,T2>::getTraveltimes(const std::vector<sxz<T1> >& Rx,
//									std::vector<T1>& traveltimes,
//									const size_t threadNo) const {
//    
//    if ( check_pts(Rx) == 1 ) return 1;
//    
//    if ( traveltimes.size() != Rx.size() ) {
//        traveltimes.resize( Rx.size() );
//    }
//    
//    for (size_t n=0; n<Rx.size(); ++n) {
//        traveltimes[n] = getTraveltime(Rx[n], this->nodes, threadNo);
//    }
//    return 0;
//}

//template<typename T1, typename T2>
//int Grid2Ducsp<T1,T2>::getTraveltimes(const std::vector<sxz<T1> >& Tx,
//									  const std::vector<sxz<T1> >& Rx,
//									  std::vector<T1>& traveltimes,
//									  std::vector<std::vector<sxz<T1> > >& r_data,
//									  const size_t threadNo) const {
//    
//    if ( check_pts(Rx) == 1 ) return 1;
//	
//    std::vector<Node2Dcsp<T1,T2> > txNodes;
//	for (size_t n=0; n<Tx.size(); ++n) {
//        bool found = false;
//        for ( size_t nn=0; nn<this->nodes.size(); ++nn ) {
//            if ( this->nodes[nn] == Tx[n] ) {
//                found = true;
//                break;
//            }
//        }
//        if ( found==false ) {
//            txNodes.push_back( Node2Dcsp<T1,T2>(0.0, Tx[n].x, Tx[n].z,
//											  this->nThreads, threadNo) );
//            // we belong to cell index no
//            txNodes.back().pushOwner( getCellNo(Tx[n]) );
//            txNodes.back().setGridIndex( static_cast<T2>(this->nodes.size()+
//                                                         txNodes.size())-1 );
//        }
//    }
//    
//    if ( traveltimes.size() != Rx.size() ) {
//        traveltimes.resize( Rx.size() );
//    }
//    if ( r_data.size() != Rx.size() ) {
//        r_data.resize( Rx.size() );
//    }
//    for ( size_t ni=0; ni<r_data.size(); ++ni ) {
//        r_data[ni].resize( 0 );
//    }
//    T2 nodeParentRx;
//    T2 cellParentRx;
//    
//    for (size_t n=0; n<Rx.size(); ++n) {
//        traveltimes[n] = this->getTraveltime(Rx[n], this->nodes, nodeParentRx, cellParentRx,
//									   threadNo);
//        
//        // Rx are in nodes (not txNodes)
//        std::vector<Node2Dcsp<T1,T2> > *node_p;
//        node_p = &(this->nodes);
//        
//        std::vector<sxz<T1> > r_tmp;
//        T2 iChild, iParent = nodeParentRx;
//        sxz<T1> child;
//		
//        // store the son's coord
//        child.x = Rx[n].x;
//        child.z = Rx[n].z;
//        while ( (*node_p)[iParent].getNodeParent(threadNo) !=
//               std::numeric_limits<T2>::max() ) {
// 			
//			r_tmp.push_back( child );
//			
//			// we now go up in time - parent becomes the child of grand'pa
//			iChild = iParent;
//			child.x = (*node_p)[iChild].getX();
//			child.z = (*node_p)[iChild].getZ();
//            
//			// grand'pa is now papa
//			iParent = (*node_p)[iChild].getNodeParent(threadNo);
//            if ( iParent >= this->nodes.size() ) {
//                node_p = &txNodes;
//                iParent -= this->nodes.size();
//            }
//            else {
//                node_p = &(this->nodes);
//            }
//		}
//		
//		// parent is now at Tx
//		r_tmp.push_back( child );
//		
//		// finally, store Tx position
//		child.x = (*node_p)[iParent].getX();
//		child.z = (*node_p)[iParent].getZ();
//		r_tmp.push_back( child );
//		
//        // the order should be from Tx to Rx, so we reorder...
//        iParent = static_cast<T2>(r_tmp.size());
//        r_data[n].resize( r_tmp.size() );
//        for ( size_t nn=0; nn<r_data[n].size(); ++nn ) {
//            r_data[n][nn].x = r_tmp[ iParent-1-nn ].x;
//            r_data[n][nn].z = r_tmp[ iParent-1-nn ].z;
//        }
//    }
//    return 0;
//}


template<typename T1, typename T2>
void Grid2Ducsp<T1,T2>::initQueue(const std::vector<sxz<T1> >& Tx,
								  const std::vector<T1>& t0,
								  std::priority_queue<Node2Dcsp<T1,T2>*,
								  std::vector<Node2Dcsp<T1,T2>*>,
								  CompareNodePtr<T1> >& queue,
								  std::vector<Node2Dcsp<T1,T2> >& txNodes,
								  std::vector<bool>& inQueue,
								  std::vector<bool>& frozen,
								  const size_t threadNo) const {
    
    for (size_t n=0; n<Tx.size(); ++n) {
        bool found = false;
        for ( size_t nn=0; nn<this->nodes.size(); ++nn ) {
            if ( this->nodes[nn] == Tx[n] ) {
                found = true;
                this->nodes[nn].setTT( t0[n], threadNo );
                queue.push( &(this->nodes[nn]) );
                inQueue[nn] = true;
                frozen[nn] = true;
                break;
            }
        }
        if ( found==false ) {
            txNodes.push_back( Node2Dcsp<T1,T2>(t0[n], Tx[n].x, Tx[n].z, this->nThreads,
												threadNo) );
            // we belong to cell index no
            txNodes.back().pushOwner( this->getCellNo(Tx[n]) );
            txNodes.back().setGridIndex( static_cast<T2>(this->nodes.size()+
														 txNodes.size()-1) );
            
            queue.push( &(txNodes.back()) );
            inQueue.push_back( true );
            frozen.push_back( true );
        }
    }
}


template<typename T1, typename T2>
void Grid2Ducsp<T1,T2>::propagate(std::priority_queue<Node2Dcsp<T1,T2>*,
								  std::vector<Node2Dcsp<T1,T2>*>,
								  CompareNodePtr<T1> >& queue,
								  std::vector<bool>& inQueue,
								  std::vector<bool>& frozen,
								  const size_t threadNo) const {
//    size_t n=1;
    while ( !queue.empty() ) {
        
        
//        char filename[200];
//        sprintf(filename, "/tmp/spm_tt_%06zd.xyz", n);
//        saveTT(filename,1,threadNo);
//        n++;
        
        
        
        const Node2Dcsp<T1,T2>* source = queue.top();
        queue.pop();
        inQueue[ source->getGridIndex() ] = false;
        
        for ( size_t no=0; no<source->getOwners().size(); ++no ) {
            
            T2 cellNo = source->getOwners()[no];
            
            for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                T2 neibNo = this->neighbors[cellNo][k];
                if ( neibNo == source->getGridIndex() || frozen[neibNo] ) {
                    continue;
                }
                
                // compute dt
                T1 dt = this->computeDt(*source, this->nodes[neibNo], cellNo);
				
                if (source->getTT(threadNo)+dt < this->nodes[neibNo].getTT(threadNo)) {
                    this->nodes[neibNo].setTT( source->getTT(threadNo)+dt, threadNo );
                    this->nodes[neibNo].setnodeParent(source->getGridIndex(),threadNo);
                    this->nodes[neibNo].setCellParent(cellNo, threadNo );
                    
                    if ( !inQueue[neibNo] ) {
                        queue.push( &(this->nodes[neibNo]) );
                        inQueue[neibNo] = true;
                    }
                }
            }
        }
    }
}


#endif

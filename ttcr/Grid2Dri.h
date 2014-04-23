//
//  Grid2Dri.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-04-22.
//  Copyright (c) 2014 Bernard Giroux. All rights reserved.
//

#ifndef ttcr_Grid2Dri_h
#define ttcr_Grid2Dri_h



#include <algorithm>
#include <cstring>
#include <iostream>
#include <fstream>
#include <queue>
#include <vector>

#include "Grid2D.h"
#include "Node2Disp.h"

template<typename T1, typename T2>
class Grid2Dri : public Grid2D<T1,T2,sxz<T1>> {
public:
    Grid2Dri(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
			 const T1 minx, const T1 minz, const T2 nnx, const T2 nnz,
			 const size_t nt=1);
    
    virtual ~Grid2Dri() {
    }
    
    T1 getDx() const { return dx; }
    T1 getDz() const { return dz; }
    T1 getXmin() const { return xmin; }
    T1 getXmax() const { return xmax; }
    T1 getZmin() const { return zmin; }
    T1 getZmax() const { return zmax; }
    T2 getNcellx() const { return nCellx; }
    T2 getNcellz() const { return nCellz; }
    T2 getNsnx() const { return nsnx; }
    T2 getNsnz() const { return nsnz; }
    
	T1 getSlowness(const size_t n) const {
		return nodes[n].getNodeSlowness();
	}
    
    void setSlowness(const T1 s) {
        for ( size_t n=0; n<nodes.size(); ++n ) {
        	nodes[n].setNodeSlowness( s );
        }
    }
    
    int setSlowness(const std::vector<T1>& s) {
        if ( nPrimary != s.size() ) {
            std::cerr << "Error: slowness vectors of incompatible size.";
            return 1;
        }
		size_t nn=0;
        for ( size_t n=0; n<nodes.size(); ++n ) {
			if ( nodes[n].getPrimary() == 5 )
				nodes[n].setNodeSlowness( s[nn++] );
        }
		interpSlownessSecondary();
        return 0;
    }
    
    int raytrace(const std::vector<sxz<T1>>& Tx,
                 const std::vector<T1>& t0,
                 const std::vector<sxz<T1>>& Rx,
                 std::vector<T1>& traveltimes,
				 const size_t threadNo=0) const;
    
    int raytrace(const std::vector<sxz<T1>>&,
                 const std::vector<T1>&,
                 const std::vector<const std::vector<sxz<T1>>*>&,
                 std::vector<std::vector<T1>*>&,
                 const size_t=0) const;
	
    int raytrace(const std::vector<sxz<T1>>& Tx,
                 const std::vector<T1>& t0,
                 const std::vector<sxz<T1>>& Rx,
                 std::vector<T1>& traveltimes,
                 std::vector<std::vector<sxz<double>>>& r_data,
				 const size_t threadNo=0) const;
    
    int raytrace(const std::vector<sxz<T1>>&,
                 const std::vector<T1>&,
                 const std::vector<const std::vector<sxz<T1>>*>&,
                 std::vector<std::vector<T1>*>&,
                 std::vector<std::vector<std::vector<sxz<T1>>>*>&,
                 const size_t=0) const;
	
    int raytrace(const std::vector<sxz<T1>>& Tx,
                 const std::vector<T1>& t0,
                 const std::vector<sxz<T1>>& Rx,
                 std::vector<T1>& traveltimes,
                 std::vector<std::vector<sxz<double>>>& r_data,
                 std::vector<std::vector<siv<double>>>& l_data,
				 const size_t threadNo=0) const;
		
    size_t getNumberOfNodes() const { return nodes.size(); }
    
    void saveSlownessXYZ(const char filename[]) const {
        std::ofstream fout( filename );
        
		for ( size_t n=0; n<nodes.size(); ++n ) {
			fout << nodes[n].getX() << "   " << nodes[n].getZ() << "   " << nodes[n].getNodeSlowness() << '\n';
        }
        
        fout.close();
    }
    
    T2 getCellNo(const sxz<T1>& pt) const {
        T1 x = xmax-pt.x < small ? xmax-.5*dx : pt.x;
        T1 z = zmax-pt.z < small ? zmax-.5*dz : pt.z;
        T2 nx = static_cast<T2>( small + (x-xmin)/dx );
        T2 nz = static_cast<T2>( small + (z-zmin)/dz );
        return nx*nCellz + nz;
    }
    
    
protected:
	size_t nThreads;
    T1 dx;           // cell size in x
    T1 dz;           // cell size in z
    T1 xmin;         // x origin of the grid
    T1 zmin;         // z origin of the grid
    T1 xmax;         // x end of the grid
    T1 zmax;         // z end of the grid
    T2 nCellx;  // number of cells in x
    T2 nCellz;  // number of cells in x
    T2 nsnx;    // number of secondary nodes in x
    T2 nsnz;    // number of secondary nodes in z
    T2 nsgx;    // number of subgrid cells in x
    T2 nsgz;    // number of subgrid cells in z
	T2 nPrimary;
    
    mutable std::vector<Node2Disp<T1,T2>> nodes;
    
    std::vector<std::vector<T2>> neighbors;  // nodes common to a cell
    
    void buildGridNodes();
    void buildGridNeighbors();
	
	void interpSlownessSecondary();

    void propagate(std::priority_queue<Node2Disp<T1,T2>*,
				   std::vector<Node2Disp<T1,T2>*>,
                   CompareNodePtr<T1>>& queue,
                   std::vector<bool>& inQueue,
                   std::vector<bool>& frozen,
				   const size_t threadNo) const;

	T1 computeDt(const Node2Disp<T1,T2>& source, const Node2Disp<T1,T2>& node) const {
		return (node.getNodeSlowness()+source.getNodeSlowness())/2 * source.getDistance( node );
	}
    
	T1 computeDt(const Node2Disp<T1,T2>& source, const sxz<T1>& node, T1 slo) const {
		return (slo+source.getNodeSlowness())/2 * source.getDistance( node );
	}

    
    int check_pts(const std::vector<sxz<T1>>&) const;
    
    void initQueue(const std::vector<sxz<T1>>& Tx,
                   const std::vector<T1>& t0,
                   std::priority_queue<Node2Disp<T1,T2>*,
				   std::vector<Node2Disp<T1,T2>*>,
                   CompareNodePtr<T1>>& queue,
                   std::vector<Node2Disp<T1,T2>>& txNodes,
                   std::vector<bool>& inQueue,
                   std::vector<bool>& frozen,
				   const size_t threadNo) const;
    
    bool inPolygon(const sxz<T1>& p, const sxz<T1> poly[], const size_t N) const;
    
    T1 getTraveltime(const sxz<T1>& Rx, const std::vector<Node2Disp<T1,T2>>& nodes,
					 const size_t threadNo) const;
    
    T1 getTraveltime(const sxz<T1>& Rx, const std::vector<Node2Disp<T1,T2>>& nodes,
					 T2& nodeParentRx, T2& cellParentRx,
					 const size_t threadNo) const;
	
//    void save(const char filename[]) const;
    
	
private:
    Grid2Dri() {}
    Grid2Dri(const Grid2Dri<T1,T2>& g) {}
    Grid2Dri<T1,T2>& operator=(const Grid2Dri<T1,T2>& g) {}
    
};

template<typename T1, typename T2>
Grid2Dri<T1,T2>::Grid2Dri(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
						  const T1 minx, const T1 minz, const T2 nnx, const T2 nnz,
						  const size_t nt) : nThreads(nt),
dx(ddx), dz(ddz), xmin(minx), zmin(minz), xmax(minx+nx*ddx), zmax(minz+nz*ddz),
nCellx(nx), nCellz(nz), nsnx(nnx), nsnz(nnz), nsgx(0), nsgz(0),
nPrimary((nCellx+1) * (nCellz+1)),
nodes(std::vector<Node2Disp<T1,T2>>( // noeuds secondaires
									nCellx*nsnx*(nCellz+1) + nCellz*nsnz*(nCellx+1) +
									// noeuds primaires
									(nCellx+1) * (nCellz+1), Node2Disp<T1,T2>(nt) )),
neighbors(std::vector<std::vector<T2>>(nCellx*nCellz))
{
    buildGridNodes();
    buildGridNeighbors();
}

template<typename T1, typename T2>
void Grid2Dri<T1,T2>::buildGridNodes() {
    T1 dxs = dx/(nsnx+1);
    T1 dzs = dz/(nsnz+1);
    
    T2 cell_upLeft = std::numeric_limits<T2>::max();
    T2 cell_upRight = std::numeric_limits<T2>::max();
    T2 cell_downLeft = 0;
    T2 cell_downRight = 0;
    
    for ( T2 n=0, nc=0; nc<=nCellx; ++nc ) {
        
        double x = xmin + nc*dx;
        
        for ( T2 nr=0; nr<=nCellz; ++nr ) {
            
            double z = zmin + nr*dz;
            
            if ( nr < nCellz && nc < nCellx ) {
                cell_downRight = nc*nCellz + nr;
            }
            else {
                cell_downRight = std::numeric_limits<T2>::max();
            }
            
            if ( nr > 0 && nc < nCellx ) {
                cell_upRight = nc*nCellz + nr - 1;
            }
            else {
                cell_upRight = std::numeric_limits<T2>::max();
            }
            
            if ( nr < nCellz && nc > 0 ) {
                cell_downLeft = (nc-1)*nCellz + nr;
            }
            else {
                cell_downLeft = std::numeric_limits<T2>::max();
            }
            
            if ( nr > 0 && nc > 0 ) {
                cell_upLeft = (nc-1)*nCellz + nr - 1;
            }
            else {
                cell_upLeft = std::numeric_limits<T2>::max();
            }
            
            // primary nodes
            //            std::cout << n << "\t p\t-\t" << x << '\t' << z
            //            << "\t-\t" << cell_upLeft
            //            << '\t' << cell_downLeft
            //            << '\t' << cell_upRight
            //            << '\t' << cell_downRight << '\n';
            
            if ( cell_upLeft != std::numeric_limits<T2>::max() ) {
                nodes[n].pushOwner( cell_upLeft );
            }
            if ( cell_downLeft != std::numeric_limits<T2>::max() ) {
                nodes[n].pushOwner( cell_downLeft );
            }
            if ( cell_upRight != std::numeric_limits<T2>::max() ) {
                nodes[n].pushOwner( cell_upRight );
            }
            if ( cell_downRight != std::numeric_limits<T2>::max() ) {
                nodes[n].pushOwner( cell_downRight );
            }
            
            nodes[n].setX( x );
            nodes[n].setZ( z );
            nodes[n].setGridIndex( n );
			nodes[n].setPrimary(5);
            
            ++n;
            
            // secondary nodes on the vertical
            if ( nr < nCellz ) {
                for (T2 ns=0; ns<nsnz; ++ns, ++n ) {
                    
                    double zsv = zmin + nr*dz + (ns+1)*dzs;
                    
                    //                    std::cout << n << "\tsv\t-\t" << x << '\t' << zsv << "\t-\t"
                    //                    << cell_downLeft << '\t' << cell_downRight << '\n';
                    
                    if ( cell_downLeft != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cell_downLeft );
                    }
                    if ( cell_downRight != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cell_downRight );
                    }
                    
                    nodes[n].setX( x );
                    nodes[n].setZ( zsv );
                    nodes[n].setGridIndex( n );
                }
            }
            
            // secondary nodes on the horizontal
            if ( nc < nCellx ) {
                for ( T2 ns=0; ns<nsnx; ++ns, ++n ) {
                    
                    double xsh = xmin + nc*dx + (ns+1)*dxs;
                    
                    //                    std::cout << n << "\tsh\t-\t" << xsh << '\t' << z << "\t-\t"
                    //                    << cell_upRight << '\t' << cell_downRight << '\n';
                    
                    if ( cell_upRight != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cell_upRight );
                    }
                    if ( cell_downRight != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cell_downRight );
                    }
                    
                    nodes[n].setX( xsh );
                    nodes[n].setZ( z );
                    nodes[n].setGridIndex( n );
                }
            }
            //            std::cout << '\n';
        }
        //        std::cout << '\n';
    }
}

template<typename T1, typename T2>
void Grid2Dri<T1,T2>::buildGridNeighbors() {
    for ( T2 n=0; n<nodes.size(); ++n ) {
        for ( size_t n2=0; n2<nodes[n].getOwners().size(); ++n2) {
            neighbors[ nodes[n].getOwners()[n2] ].push_back(n);
        }
    }
}

template<typename T1, typename T2>
void Grid2Dri<T1,T2>::interpSlownessSecondary() {  // TODO : test this
	
	T1 dxs = dx/(nsnx+1);
    T1 dzs = dz/(nsnz+1);
    
    for ( T2 n=0, nc=0; nc<=nCellx; ++nc ) {
        
        for ( T2 nr=0; nr<=nCellz; ++nr ) {
            
			T2 np1 = n;
			T2 np2v = np1 + nsnx+nsnz;
			T2 np2h = np1 + nCellz*nsnx + (nCellz-1)*nsnz;
            
            ++n;
			
            // secondary nodes on the vertical
            if ( nr < nCellz ) {
				T1 slope = (nodes[np2v].getNodeSlowness()-nodes[np1].getNodeSlowness())/dz;
                for (T2 ns=0; ns<nsnz; ++ns, ++n ) {
                    
                    T1 s = nodes[np1].getNodeSlowness() + slope * (ns+1)*dzs;
                    
                    nodes[n].setNodeSlowness( s );
                }
            }
			
            // secondary nodes on the horizontal
            if ( nc < nCellx ) {
				T1 slope = (nodes[np2h].getNodeSlowness()-nodes[np1].getNodeSlowness())/dx;
                for ( T2 ns=0; ns<nsnx; ++ns, ++n ) {
                    
                    T1 s = nodes[np1].getNodeSlowness() + slope * (ns+1)*dxs;
                    
                    nodes[n].setNodeSlowness( s );
                }
            }
        }
    }
}

template<typename T1, typename T2>
int Grid2Dri<T1,T2>::check_pts(const std::vector<sxz<T1>>& pts) const {
    for (size_t n=0; n<pts.size(); ++n) {
        if ( pts[n].x < xmin || pts[n].x > xmax ||
            pts[n].z < zmin || pts[n].z > zmax ) {
            std::cerr << "Error: point no " << (n+1)
            << " outside the grid.\n";
            return 1;
        }
    }
    return 0;
}


template<typename T1, typename T2>
void Grid2Dri<T1,T2>::initQueue(const std::vector<sxz<T1>>& Tx,
								const std::vector<T1>& t0,
								std::priority_queue<Node2Disp<T1,T2>*,
								std::vector<Node2Disp<T1,T2>*>,
								CompareNodePtr<T1>>& queue,
								std::vector<Node2Disp<T1,T2>>& txNodes,
								std::vector<bool>& inQueue,
								std::vector<bool>& frozen,
								const size_t threadNo) const {
    
    for (size_t n=0; n<Tx.size(); ++n) {
        bool found = false;
        for ( size_t nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == Tx[n] ) {
                found = true;
                nodes[nn].setTT( t0[n], threadNo );
                queue.push( &(nodes[nn]) );
                inQueue[nn] = true;
                frozen[nn] = true;
                break;
            }
        }
        if ( found==false ) {
            txNodes.push_back( Node2Disp<T1,T2>(t0[n], Tx[n], nThreads,
												threadNo) );
            // we belong to cell index no
            txNodes.back().pushOwner( getCellNo(Tx[n]) );
            txNodes.back().setGridIndex( static_cast<T2>(nodes.size()+
														 txNodes.size()-1) );
            
            queue.push( &(txNodes.back()) );
            inQueue.push_back( true );
            frozen.push_back( true );
        }
    }
}


template<typename T1, typename T2>
int Grid2Dri<T1,T2>::raytrace(const std::vector<sxz<T1>>& Tx,
							  const std::vector<T1>& t0,
							  const std::vector<sxz<T1>>& Rx,
							  std::vector<T1>& traveltimes,
							  const size_t threadNo) const {
    
    if ( check_pts(Tx) == 1 ) return 1;
    if ( check_pts(Rx) == 1 ) return 1;
    
    for ( size_t n=0; n<nodes.size(); ++n ) {
        nodes[n].reinit( threadNo );
    }
    
    CompareNodePtr<T1> cmp(threadNo);
    std::priority_queue< Node2Disp<T1,T2>*, std::vector<Node2Disp<T1,T2>*>,
    CompareNodePtr<T1>> queue( cmp );
    
    std::vector<Node2Disp<T1,T2>> txNodes;
    std::vector<bool> inQueue( nodes.size(), false );
    std::vector<bool> frozen( nodes.size(), false );
    
    initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
    
    propagate(queue, inQueue, frozen, threadNo);
    
    if ( traveltimes.size() != Rx.size() ) {
        traveltimes.resize( Rx.size() );
    }
    
    for (size_t n=0; n<Rx.size(); ++n) {
        traveltimes[n] = getTraveltime(Rx[n], nodes, threadNo);
    }
    return 0;
}

template<typename T1, typename T2>
int Grid2Dri<T1,T2>::raytrace(const std::vector<sxz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<const std::vector<sxz<T1>>*>& Rx,
                              std::vector<std::vector<T1>*>& traveltimes,
                              const size_t threadNo) const {
    if ( check_pts(Tx) == 1 ) return 1;
    for ( size_t n=0; n<Rx.size(); ++n )
        if ( check_pts(*Rx[n]) == 1 ) return 1;
    
    for ( size_t n=0; n<nodes.size(); ++n ) {
        nodes[n].reinit( threadNo );
    }
    
    CompareNodePtr<T1> cmp(threadNo);
    std::priority_queue< Node2Disp<T1,T2>*, std::vector<Node2Disp<T1,T2>*>,
    CompareNodePtr<T1>> queue( cmp );
    
    std::vector<Node2Disp<T1,T2>> txNodes;
    std::vector<bool> inQueue( nodes.size(), false );
    std::vector<bool> frozen( nodes.size(), false );
    
    initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
    
    propagate(queue, inQueue, frozen, threadNo);
    
    if ( traveltimes.size() != Rx.size() ) {
        traveltimes.resize( Rx.size() );
    }
    
    for (size_t nr=0; nr<Rx.size(); ++nr) {
        traveltimes[nr]->resize( Rx[nr]->size() );
        for (size_t n=0; n<Rx[nr]->size(); ++n)
            (*traveltimes[nr])[n] = getTraveltime((*Rx[nr])[n], nodes, threadNo);
    }
    return 0;
	
}

template<typename T1, typename T2>
int Grid2Dri<T1,T2>::raytrace(const std::vector<sxz<T1>>& Tx,
							  const std::vector<T1>& t0,
							  const std::vector<sxz<T1>>& Rx,
							  std::vector<T1>& traveltimes,
							  std::vector<std::vector<sxz<double>>>& r_data,
							  const size_t threadNo) const {
    
    if ( check_pts(Tx) == 1 ) return 1;
    if ( check_pts(Rx) == 1 ) return 1;
    
    for ( size_t n=0; n<nodes.size(); ++n ) {
        nodes[n].reinit( threadNo );
    }
    
    CompareNodePtr<T1> cmp(threadNo);
    std::priority_queue< Node2Disp<T1,T2>*, std::vector<Node2Disp<T1,T2>*>,
    CompareNodePtr<T1>> queue( cmp );
    std::vector<Node2Disp<T1,T2>> txNodes;
    std::vector<bool> inQueue( nodes.size(), false );
    std::vector<bool> frozen( nodes.size(), false );
    
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
        
        traveltimes[n] = getTraveltime(Rx[n], nodes, nodeParentRx, cellParentRx,
									   threadNo);
        
        // Rx are in nodes (not txNodes)
        std::vector<Node2Disp<T1,T2>> *node_p;
        node_p = &nodes;
        
        std::vector<sxz<double>> r_tmp;
        T2 iChild, iParent = nodeParentRx;
        sxz<double> child;
		
        // store the son's coord
        child.x = Rx[n].x;
        child.z = Rx[n].z;
        while ( (*node_p)[iParent].getNodeParent(threadNo) != std::numeric_limits<T2>::max() ) {
 			
			r_tmp.push_back( child );
			
			// we now go up in time - parent becomes the child of grand'pa
			iChild = iParent;
			child.x = (*node_p)[iChild].getX();
			child.z = (*node_p)[iChild].getZ();
			
			// grand'pa is now papa
			iParent = (*node_p)[iChild].getNodeParent(threadNo);
            if ( iParent >= nodes.size() ) {
                node_p = &txNodes;
                iParent -= nodes.size();
            }
            else {
                node_p = &nodes;
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
int Grid2Dri<T1,T2>::raytrace(const std::vector<sxz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<const std::vector<sxz<T1>>*>& Rx,
                              std::vector<std::vector<T1>*>& traveltimes,
                              std::vector<std::vector<std::vector<sxz<T1>>>*>& r_data,
                              const size_t threadNo) const {
    
    if ( check_pts(Tx) == 1 ) return 1;
    for ( size_t n=0; n<Rx.size(); ++n )
        if ( check_pts(*Rx[n]) == 1 ) return 1;
    
    for ( size_t n=0; n<nodes.size(); ++n ) {
        nodes[n].reinit( threadNo );
    }
    
    CompareNodePtr<T1> cmp(threadNo);
    std::priority_queue< Node2Disp<T1,T2>*, std::vector<Node2Disp<T1,T2>*>,
    CompareNodePtr<T1>> queue( cmp );
    
    std::vector<Node2Disp<T1,T2>> txNodes;
    std::vector<bool> inQueue( nodes.size(), false );
    std::vector<bool> frozen( nodes.size(), false );
    
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
            
            (*traveltimes[nr])[n] = getTraveltime((*Rx[nr])[n], nodes,
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
            std::vector<Node2Disp<T1,T2>> *node_p;
            node_p = &nodes;
            
            std::vector<sxz<T1>> r_tmp;
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
                if ( iParent >= nodes.size() ) {
                    node_p = &txNodes;
                    iParent -= nodes.size();
                }
                else {
                    node_p = &nodes;
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
int Grid2Dri<T1,T2>::raytrace(const std::vector<sxz<T1>>& Tx,
							  const std::vector<T1>& t0,
							  const std::vector<sxz<T1>>& Rx,
							  std::vector<T1>& traveltimes,
							  std::vector<std::vector<sxz<double>>>& r_data,
							  std::vector<std::vector<siv<double>>>& l_data,
							  const size_t threadNo) const {
    
    if ( check_pts(Tx) == 1 ) return 1;
    if ( check_pts(Rx) == 1 ) return 1;
    
    for ( size_t n=0; n<nodes.size(); ++n ) {
        nodes[n].reinit( threadNo );
    }
    
    CompareNodePtr<T1> cmp(threadNo);
    std::priority_queue< Node2Disp<T1,T2>*, std::vector<Node2Disp<T1,T2>*>,
    CompareNodePtr<T1>> queue( cmp );
    std::vector<Node2Disp<T1,T2>> txNodes;
    std::vector<bool> inQueue( nodes.size(), false );
    std::vector<bool> frozen( nodes.size(), false );
    
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
        
        traveltimes[n] = getTraveltime(Rx[n], nodes, nodeParentRx, cellParentRx,
									   threadNo);
        
        // Rx are in nodes (not txNodes)
        std::vector<Node2Disp<T1,T2>> *node_p;
        node_p = &nodes;
        
        std::vector<sxz<double>> r_tmp;
        T2 iChild, iParent = nodeParentRx;
        sxz<double> child;
        siv<double> cell;
		
        // store the son's coord
        child.x = Rx[n].x;
        child.z = Rx[n].z;
        cell.i = cellParentRx;
        while ( (*node_p)[iParent].getNodeParent(threadNo) != std::numeric_limits<T2>::max() ) {
 			
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
            if ( iParent >= nodes.size() ) {
                node_p = &txNodes;
                iParent -= nodes.size();
            }
            else {
                node_p = &nodes;
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
        iParent = r_tmp.size();
        r_data[n].resize( r_tmp.size() );
        for ( size_t nn=0; nn<r_data[n].size(); ++nn ) {
            r_data[n][nn].x = r_tmp[ iParent-1-nn ].x;
            r_data[n][nn].z = r_tmp[ iParent-1-nn ].z;
        }
    }
    return 0;
}


template<typename T1, typename T2>
void Grid2Dri<T1,T2>::propagate( std::priority_queue<Node2Disp<T1,T2>*, std::vector<Node2Disp<T1,T2>*>,
								CompareNodePtr<T1>>& queue,
								std::vector<bool>& inQueue,
								std::vector<bool>& frozen,
								const size_t threadNo) const {
    
    while ( !queue.empty() ) {
        const Node2Disp<T1,T2>* source = queue.top();
        queue.pop();
        inQueue[ source->getGridIndex() ] = false;
        
        for ( size_t no=0; no<source->getOwners().size(); ++no ) {
            
            T2 cellNo = source->getOwners()[no];
            
            for ( size_t k=0; k< neighbors[cellNo].size(); ++k ) {
                T2 neibNo = neighbors[cellNo][k];
                if ( neibNo == source->getGridIndex() || frozen[neibNo] ) {
                    continue;
                }
                
                // compute dt
                T1 dt = computeDt(*source, nodes[neibNo]);
				
                if ( source->getTT(threadNo)+dt < nodes[neibNo].getTT(threadNo) ) {
                    nodes[neibNo].setTT( source->getTT(threadNo)+dt, threadNo );
                    nodes[neibNo].setnodeParent( source->getGridIndex(), threadNo );
                    nodes[neibNo].setCellParent( cellNo, threadNo );
                    
                    if ( !inQueue[neibNo] ) {
                        queue.push( &(nodes[neibNo]) );
                        inQueue[neibNo] = true;
                    }
                }
            }
        }
    }
}





template<typename T1, typename T2>
bool Grid2Dri<T1,T2>::inPolygon(const sxz<T1>& p, const sxz<T1> poly[], const size_t N) const {
    bool c = false;
    for (size_t i = 0, j = N-1; i < N; j = i++) {
        if ((((poly[i].z <= p.z) && (p.z < poly[j].z)) ||
             ((poly[j].z <= p.z) && (p.z < poly[i].z))) &&
            (p.x < (poly[j].x - poly[i].x) * (p.z - poly[i].z) / (poly[j].z - poly[i].z) + poly[i].x))
            c = !c;
    }
    return c;
}


template<typename T1, typename T2>
T1 Grid2Dri<T1,T2>::getTraveltime(const sxz<T1>& Rx,
								  const std::vector<Node2Disp<T1,T2>>& nodes,
								  const size_t threadNo) const {
    
    for ( size_t nn=0; nn<nodes.size(); ++nn ) {
        if ( nodes[nn] == Rx ) {
            return nodes[nn].getTT(threadNo);
        }
    }
    
    T2 cellNo = getCellNo( Rx );
    T2 neibNo = neighbors[cellNo][0];
    T1 dt = computeDt(nodes[neibNo], Rx, cellNo);
    
    T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
    for ( size_t k=1; k< neighbors[cellNo].size(); ++k ) {
        neibNo = neighbors[cellNo][k];
        dt = computeDt(nodes[neibNo], Rx, cellNo);
        if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
            traveltime =  nodes[neibNo].getTT(threadNo)+dt;
        }
    }
    return traveltime;
}


template<typename T1, typename T2>
T1 Grid2Dri<T1,T2>::getTraveltime(const sxz<T1>& Rx,
								  const std::vector<Node2Disp<T1,T2>>& nodes,
								  T2& nodeParentRx, T2& cellParentRx,
								  const size_t threadNo) const {
    
    for ( size_t nn=0; nn<nodes.size(); ++nn ) {
        if ( nodes[nn] == Rx ) {
            nodeParentRx = nodes[nn].getNodeParent(threadNo);
            cellParentRx = nodes[nn].getCellParent(threadNo);
            return nodes[nn].getTT(threadNo);
        }
    }
    
    T2 cellNo = getCellNo( Rx );
    T2 neibNo = neighbors[cellNo][0];
    T1 dt = computeDt(nodes[neibNo], Rx, cellNo);
    
    T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
    nodeParentRx = neibNo;
    cellParentRx = cellNo;
    for ( size_t k=1; k< neighbors[cellNo].size(); ++k ) {
        neibNo = neighbors[cellNo][k];
        dt = computeDt(nodes[neibNo], Rx, cellNo);
        if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
            traveltime =  nodes[neibNo].getTT(threadNo)+dt;
            nodeParentRx = neibNo;
        }
    }
    return traveltime;
}

//template<typename T1, typename T2>
//void Grid2Dri<T1,T2>::save(const char filename[]) const {
//	std::ofstream fout( filename );
//	
//	fout << dx << ' ' << dz << ' ' << xmin << ' ' << zmin << ' ' << xmax << ' ' << zmax << '\n';
//	fout << nCellx << ' ' << nCellz << ' ' << nsnx << ' ' << nsnz << ' ' << nsgx << ' ' << nsgz << '\n';
//	fout << nodes.size() << '\n';
//	for ( size_t n=0; n < nodes.size(); ++n ) {
//		fout << nodes[n].getNodeSlowness() << '\n';
//	}
//	fout.close();
//}


#endif

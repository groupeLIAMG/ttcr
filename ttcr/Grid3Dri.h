//
//  Grid3Dri.h
//  ttcr.v2
//
//  Created by Giroux Bernard on 12-08-15.
//  Copyright (c) 2012 INRS-ETE. All rights reserved.
//

#ifndef __GRID3DRI_H__
#define __GRID3DRI_H__


#include <algorithm>
#include <cstring>
#include <iostream>
#include <fstream>
#include <queue>
#include <vector>
#include <ctime>

#ifdef VTK
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkXMLRectilinearGridWriter.h"
#endif

#include "Interpolator.h"


template<typename T1, typename T2>
class Grid3Dri : public Grid3Dr<T1,T2> {
public:
    Grid3Dri(const T2 nx, const T2 ny, const T2 nz,
			 const T1 ddx, const T1 ddy, const T1 ddz,
			 const T1 minx, const T1 miny, const T1 minz,
			 const T2 nnx, const T2 nny, const T2 nnz,
			 const size_t nt=1, const bool invDist=false);
    
    ~Grid3Dri() {}
    
    T1 getSlowness(const size_t n) const {
		return nodes[n].getNodeSlowness();
	}
    
    void setSlowness(const T1 s) {
        for ( size_t n=0; n<nodes.size(); ++n ) {
        	nodes[n].setNodeSlowness( s );
        }
    }
	
    int setSlowness(const std::vector<T1>& s);
    
    int linearInterpolation() const;
    
    int invDistInterpolation() const;
    
    size_t getNumberOfNodes() const { return nodes.size(); }
	
    int raytrace(const std::vector<sxyz<T1> >& Tx,
                 const std::vector<T1>& t0,
                 const std::vector<sxyz<T1> >& Rx,
                 std::vector<T1>& traveltimes,
                 const size_t threadNo=0) const;
	
	int raytrace(const std::vector<sxyz<T1> >& Tx,
                 const std::vector<T1>& t0,
                 const std::vector<sxyz<T1> >& Rx,
                 std::vector<T1>& traveltimes,
				 std::vector<std::vector<sxyz<T1> > >& r_data,
                 const size_t threadNo=0) const;
	
    void saveSlownessXYZ(const char filename[]) const {
    	//Saves the Slowness of the primary nodes
        std::ofstream fout( filename );
        for ( size_t n=0; n< Grid3Dr<T1,T2>::nodes.size(); ++n ) {
        	if (Grid3Dr<T1,T2>::nodes[n].isPrimary() == 5 ){
				fout << Grid3Dr<T1,T2>::nodes[n].getX() << "   "
				<< Grid3Dr<T1,T2>::nodes[n].getY() << "   "
				<< Grid3Dr<T1,T2>::nodes[n].getZ() << "   "
				<< Grid3Dr<T1,T2>::nodes[n].getNodeSlowness() << '\n';
        	}
        }
        fout.close();
    }
        
    void save(const char filename[]) const;
    void dsave(const char filename[]) const;
    void savefast(const char filename[]) const;
    void savePrimary(const char filename[], const size_t nt=0,
					 const bool vtkFormat=0) const;
	
    size_t getNeighborsSize() const {
        size_t n_elem = 0;
        for ( size_t n=0; n<neighbors.size(); ++n ) {
            n_elem += neighbors[n].size();
        }
        return n_elem*sizeof(size_t);
    }
    size_t getNodesSize() const {
        size_t size = 0;
        for ( size_t n=0; n<nodes.size(); ++n ) {
            size += nodes[n].getSize();
        }
        return size;
    }
	
private:
    
    bool inverseDistance;
    mutable std::vector<Node3Di<T1,T2> > nodes;
    std::vector<std::vector<T2> > neighbors;  // nodes common to a cell
    
    void buildGridNodes();
    void buildGridNeighbors();
	
    void propagate(std::priority_queue<Node3Di<T1,T2>*, std::vector<Node3Di<T1,T2>*>,
                   CompareNodePtr<T1> >& queue,
                   std::vector<bool>& inQueue,
                   std::vector<bool>& frozen,
                   size_t threadNo) const;
    
    void prepropagate(const Node3Di<T1,T2>& node,
                      std::priority_queue<Node3Di<T1,T2>*, std::vector<Node3Di<T1,T2>*>,
                      CompareNodePtr<T1> >& queue,
                      std::vector<bool>& inQueue,
                      std::vector<bool>& frozen,
                      size_t threadNo) const;
	
	T1 computeDt(const Node3Di<T1,T2>& source, const Node3Di<T1,T2>& node) const {
		return (node.getNodeSlowness()+source.getNodeSlowness())/2 * source.getDistance( node );
	}
    
	T1 computeDt(const Node3Di<T1,T2>& source, const sxyz<T1>& node, T1 slo) const {
		return (slo+source.getNodeSlowness())/2 * source.getDistance( node );
	}
    
    void initQueue(const std::vector<sxyz<T1> >& Tx,
                   const std::vector<T1>& t0,
                   std::priority_queue<Node3Di<T1,T2>*,
                   std::vector<Node3Di<T1,T2>*>,
                   CompareNodePtr<T1> >& queue,
                   std::vector<Node3Di<T1,T2> >& txNodes,
                   std::vector<bool>& inQueue,
                   std::vector<bool>& frozen,
                   const size_t threadNo) const;
    
    T1 getTraveltime(const sxyz<T1>& Rx,
					 const std::vector<Node3Di<T1,T2> >& nodes,
					 const size_t threadNo) const;
    
    T1 getTraveltime(const sxyz<T1>& Rx,
					 const std::vector<Node3Di<T1,T2> >& nodes,
					 T2&, T2&, const size_t threadNo) const;
	
    // Linear interpolation functions:
    T1 linX(const Node3Di<T1,T2>& node, const size_t C1, const size_t C2) const;
    T1 linX(const Node3Di<T1,T2>& node, const size_t C1, const size_t C2, const T1 S1, const T1 S2) const;
    T1 linX(const sxyz<T1>& node, const size_t C1, const size_t C2) const;
    T1 linX(const sxyz<T1>& node, const size_t C1, const size_t C2, const T1 S1, const T1 S2) const;
	
    T1 linY(const Node3Di<T1,T2>& node, const size_t C1, const size_t C2) const;
    T1 linY(const Node3Di<T1,T2>& node, const size_t C1, const size_t C2, const T1 S1, const T1 S2) const;
    T1 linY(const sxyz<T1>& node, const size_t C1, const size_t C2) const;
    T1 linY(const sxyz<T1>& node, const size_t C1, const size_t C2, const T1 S1, const T1 S2) const;
	
    T1 linZ(const Node3Di<T1,T2>& node, const size_t C1, const size_t C2) const;
    T1 linZ(const Node3Di<T1,T2>& node, const size_t C1, const size_t C2, const T1 S1, const T1 S2) const;
    T1 linZ(const sxyz<T1>& node, const size_t C1, const size_t C2) const;
    T1 linZ(const sxyz<T1>& node, const size_t C1, const size_t C2, const T1 S1, const T1 S2) const;
    
    T1 linXY(const Node3Di<T1,T2>& node, const size_t N1, const size_t N2, const size_t N3, const size_t N4) const ;
    T1 linXY(const sxyz<T1>& node, const size_t N1, const size_t N2, const size_t N3, const size_t N4) const ;
	
    T1 linXZ(const Node3Di<T1,T2>& node, const size_t N1, const size_t N2, const size_t N3, const size_t N4) const ;
    T1 linXZ(const sxyz<T1>& node, const size_t N1, const size_t N2, const size_t N3, const size_t N4) const ;
	
    T1 linYZ(const Node3Di<T1,T2>& node, const size_t N1, const size_t N2, const size_t N3, const size_t N4) const ;
    T1 linYZ(const sxyz<T1>& node, const size_t N1, const size_t N2, const size_t N3, const size_t N4) const ;
	
    T1 computeSlowness(const sxyz<T1>& Rx ) const;
	
    bool isNearInt( double value ) const {
		return ( remainder(value, 1.)  <= small );
	}
	
    Grid3Dri() {}
    Grid3Dri(const Grid3Dri<T1,T2>& g) {}
    Grid3Dri<T1,T2>& operator=(const Grid3Dri<T1,T2>& g) { return *this; }
    
	};
	
	
	/* Constructor Format:
	 Grid3Dri<T1,T2>::Grid3Dri(nb cells in x, nb cells in y, nb cells in z,
	 x cells size, y cells size, z cells size,
	 x origin, y origin, z origin,
	 nb sec. cells in x, nb sec. cells in y, nb sec. cells in z,
	 index of the thread)
	 */
	template<typename T1, typename T2>
	Grid3Dri<T1,T2>::Grid3Dri(const T2 nx, const T2 ny, const T2 nz,
							  const T1 ddx, const T1 ddy, const T1 ddz,
							  const T1 minx, const T1 miny, const T1 minz,
							  const T2 nnx, const T2 nny, const T2 nnz,
							  const size_t nt, const bool invDist) :
	Grid3Dr<T1,T2>(nx, ny, nz, ddx, ddy, ddz, minx, miny, minz, nnx, nny, nnz, nt),
    inverseDistance(invDist),
	nodes(std::vector<Node3Di<T1,T2> >(// secondary nodes on the edges
									   nx*nnx*((ny+1)*(nz+1)) +
									   ny*nny*((nx+1)*(nz+1)) +
									   nz*nnz*((nx+1)*(ny+1)) +
									   // secondary nodes on the faces
									   (nnx*nny)*(nx*ny*(nz+1))+
									   (nnx*nnz)*(nx*nz*(ny+1))+
									   (nny*nnz)*(ny*nz*(nx+1))+
									   // primary nodes
									   (nx+1) * (ny+1) * (nz+1),
									   Node3Di<T1,T2>(nt) )),
	neighbors(std::vector<std::vector<T2> >(nx*ny*nz))
	{
		buildGridNodes();
		buildGridNeighbors();
	}
	
	template<typename T1, typename T2>
	void Grid3Dri<T1,T2>::buildGridNodes() {
		
		// Create the grid, assign a number for each node, determine the type of the node and find the owners
		// Nodes and cells are first indexed in z, then y, and x.
		// Secondary nodes are placed on the faces and edges of every cells.
		// Ex: the node in "node[A]=(i,j,k)" is followed by the node in "node[A+1]=(i+dx,j,k)"
		
		T1 dxs = Grid3Dr<T1,T2>::dx/(Grid3Dr<T1,T2>::nsnx+1); 	// distance between secondary nodes in x
		T1 dys = Grid3Dr<T1,T2>::dy/(Grid3Dr<T1,T2>::nsny+1);
		T1 dzs = Grid3Dr<T1,T2>::dz/(Grid3Dr<T1,T2>::nsnz+1);
		
		T2 cell_XmYmZm; 	// cell in the (x-,y-,z-) direction from the node
		T2 cell_XpYmZm; 	// cell in the (x+,y-,z-) direction from the node
		T2 cell_XmYpZm;
		T2 cell_XpYpZm;
		T2 cell_XmYmZp;
		T2 cell_XpYmZp;
		T2 cell_XmYpZp;
		T2 cell_XpYpZp;
		
		T2 n=0;
		for ( T2 nk=0; nk<=Grid3Dr<T1,T2>::ncz; ++nk ) {
			
			T1 z = Grid3Dr<T1,T2>::zmin + nk*Grid3Dr<T1,T2>::dz;
			
			for ( T2 nj=0; nj<=Grid3Dr<T1,T2>::ncy; ++nj ) {
				
				T1 y = Grid3Dr<T1,T2>::ymin + nj*Grid3Dr<T1,T2>::dy;
				
				for (T2 ni=0; ni<=Grid3Dr<T1,T2>::ncx; ++ni){
					
					T1 x = Grid3Dr<T1,T2>::xmin + ni*Grid3Dr<T1,T2>::dx;
					
					// Find the adjacent cells for each primary node
					
					if (ni < Grid3Dr<T1,T2>::ncx && nj < Grid3Dr<T1,T2>::ncy && nk < Grid3Dr<T1,T2>::ncz){
						cell_XpYpZp = nj*Grid3Dr<T1,T2>::ncx + nk*(Grid3Dr<T1,T2>::ncx*Grid3Dr<T1,T2>::ncy) + ni;
					}
					else {
						cell_XpYpZp = std::numeric_limits<T2>::max();
					}
					
					if (ni > 0 && nj < Grid3Dr<T1,T2>::ncy && nk < Grid3Dr<T1,T2>::ncz){
						cell_XmYpZp = nj*Grid3Dr<T1,T2>::ncx + nk*(Grid3Dr<T1,T2>::ncx*Grid3Dr<T1,T2>::ncy) + ni - 1;
					}
					else {
						cell_XmYpZp = std::numeric_limits<T2>::max();
					}
					
					if (ni < Grid3Dr<T1,T2>::ncx && nj > 0 && nk < Grid3Dr<T1,T2>::ncz){
						cell_XpYmZp = (nj-1)*Grid3Dr<T1,T2>::ncx + nk*(Grid3Dr<T1,T2>::ncx*Grid3Dr<T1,T2>::ncy) + ni;
					}
					else {
						cell_XpYmZp = std::numeric_limits<T2>::max();
					}
					
					if (ni > 0 && nj > 0 && nk < Grid3Dr<T1,T2>::ncz){
						cell_XmYmZp = (nj-1)*Grid3Dr<T1,T2>::ncx + nk*(Grid3Dr<T1,T2>::ncx * Grid3Dr<T1,T2>::ncy) + ni - 1;
					}
					else {
						cell_XmYmZp = std::numeric_limits<T2>::max();
					}
					
					if (ni < Grid3Dr<T1,T2>::ncx && nj < Grid3Dr<T1,T2>::ncy && nk > 0){
						cell_XpYpZm = nj*Grid3Dr<T1,T2>::ncx + (nk-1)*(Grid3Dr<T1,T2>::ncx*Grid3Dr<T1,T2>::ncy) + ni;
					}
					else {
						cell_XpYpZm = std::numeric_limits<T2>::max();
					}
					
					if (ni > 0 && nj < Grid3Dr<T1,T2>::ncy && nk > 0){
						cell_XmYpZm = nj*Grid3Dr<T1,T2>::ncx + (nk-1)*(Grid3Dr<T1,T2>::ncx*Grid3Dr<T1,T2>::ncy) + ni - 1;
					}
					else {
						cell_XmYpZm = std::numeric_limits<T2>::max();
					}
					
					if (ni < Grid3Dr<T1,T2>::ncx && nj > 0 && nk > 0){
						cell_XpYmZm = (nj-1)*Grid3Dr<T1,T2>::ncx + (nk-1)*(Grid3Dr<T1,T2>::ncx*Grid3Dr<T1,T2>::ncy) + ni;
					}
					else {
						cell_XpYmZm = std::numeric_limits<T2>::max();
					}
					
					if (ni > 0 && nj > 0 && nk > 0){
						cell_XmYmZm = (nj-1)*Grid3Dr<T1,T2>::ncx + (nk-1)*(Grid3Dr<T1,T2>::ncx*Grid3Dr<T1,T2>::ncy) + ni - 1;
					}
					else {
						cell_XmYmZm = std::numeric_limits<T2>::max();
					}
					
					
					// Index the primary nodes owners
					
					if (cell_XmYmZm != std::numeric_limits<T2>::max() ) {
						nodes[n].pushOwner( cell_XmYmZm );
					}
					if (cell_XpYmZm != std::numeric_limits<T2>::max() ) {
						nodes[n].pushOwner( cell_XpYmZm );
					}
					if (cell_XmYpZm != std::numeric_limits<T2>::max() ) {
						nodes[n].pushOwner( cell_XmYpZm );
					}
					if (cell_XpYpZm != std::numeric_limits<T2>::max() ) {
						nodes[n].pushOwner( cell_XpYpZm );
					}
					if (cell_XmYmZp != std::numeric_limits<T2>::max() ) {
						nodes[n].pushOwner( cell_XmYmZp );
					}
					if (cell_XpYmZp != std::numeric_limits<T2>::max() ) {
						nodes[n].pushOwner( cell_XpYmZp );
					}
					if (cell_XmYpZp != std::numeric_limits<T2>::max() ) {
						nodes[n].pushOwner( cell_XmYpZp );
					}
					if (cell_XpYpZp != std::numeric_limits<T2>::max() ) {
						nodes[n].pushOwner( cell_XpYpZp );
					}
					
					nodes[n].setXYZindex( x, y, z, n );
					nodes[n].setPrimary(5);
					
					++n;
					
					// Secondary nodes on x edge
					if ( ni < Grid3Dr<T1,T2>::ncx ) {
						for (T2 ns=0; ns< Grid3Dr<T1,T2>::nsnx; ++ns, ++n ) {
							
							T1 xsv = Grid3Dr<T1,T2>::xmin + ni*Grid3Dr<T1,T2>::dx + (ns+1)*dxs;
							
							if ( cell_XpYmZm != std::numeric_limits<T2>::max() ) {
								nodes[n].pushOwner( cell_XpYmZm );
							}
							if ( cell_XpYpZm != std::numeric_limits<T2>::max() ) {
								nodes[n].pushOwner( cell_XpYpZm );
							}
							if ( cell_XpYmZp != std::numeric_limits<T2>::max() ) {
								nodes[n].pushOwner( cell_XpYmZp );
							}
							if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
								nodes[n].pushOwner( cell_XpYpZp );
							}
							nodes[n].setXYZindex( xsv, y, z, n );
							
							if (nj >0 && nk>0){
								nodes[n].setPrimary(28);
							}
							else if (nj==0 && nk>0){
								nodes[n].setPrimary(27);
							}
							else if (nj>0 && nk==0){
								nodes[n].setPrimary(26);
							}
							else if (nj==0 && nk==0){
								nodes[n].setPrimary(25);
							}
						}
					}
					
					// Secondary nodes on y edge
					if ( nj < Grid3Dr<T1,T2>::ncy ) {
						for (T2 ns=0; ns< Grid3Dr<T1,T2>::nsny; ++ns, ++n ) {
							
							T1 ysv = Grid3Dr<T1,T2>::ymin + nj* Grid3Dr<T1,T2>::dy + (ns+1)*dys;
							
							if ( cell_XmYpZm != std::numeric_limits<T2>::max() ) {
								nodes[n].pushOwner( cell_XmYpZm );
							}
							if ( cell_XpYpZm != std::numeric_limits<T2>::max() ) {
								nodes[n].pushOwner( cell_XpYpZm );
							}
							if ( cell_XmYpZp != std::numeric_limits<T2>::max() ) {
								nodes[n].pushOwner( cell_XmYpZp );
							}
							if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
								nodes[n].pushOwner( cell_XpYpZp );
							}
							nodes[n].setXYZindex( x, ysv, z, n );
							
							if (ni >0 && nk>0){
								nodes[n].setPrimary(38);
							}
							else if (ni>0 && nk==0){
								nodes[n].setPrimary(37);
							}
							else if (ni==0 && nk>0){
								nodes[n].setPrimary(36);
							}
							else if (ni==0 && nk==0){
								nodes[n].setPrimary(35);
							}
						}
					}
					
					// Secondary nodes on z edge
					if ( nk < Grid3Dr<T1,T2>::ncz ) {
						for (T2 ns=0; ns< Grid3Dr<T1,T2>::nsnz; ++ns, ++n ) {
							
							T1 zsv = Grid3Dr<T1,T2>::zmin + nk* Grid3Dr<T1,T2>::dz + (ns+1)*dzs;
							
							if ( cell_XmYmZp != std::numeric_limits<T2>::max() ) {
								nodes[n].pushOwner( cell_XmYmZp );
							}
							if ( cell_XpYmZp != std::numeric_limits<T2>::max() ) {
								nodes[n].pushOwner( cell_XpYmZp );
							}
							if ( cell_XmYpZp != std::numeric_limits<T2>::max() ) {
								nodes[n].pushOwner( cell_XmYpZp );
							}
							if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
								nodes[n].pushOwner( cell_XpYpZp );
							}
							nodes[n].setXYZindex( x, y, zsv, n );
							
							if (ni >0 && nj>0){
								nodes[n].setPrimary(48);
							}
							else if (ni>0 && nj==0){
								nodes[n].setPrimary(47);
							}
							else if (ni==0 && nj>0){
								nodes[n].setPrimary(46);
							}
							else if (ni==0 && nj==0){
								nodes[n].setPrimary(45);
							}
						}
					}
					
					// Secondary nodes on the xy0 planes
					if ( ni < Grid3Dr<T1,T2>::ncx && nj < Grid3Dr<T1,T2>::ncy ) {
						for(T2 sy=0; sy < Grid3Dr<T1,T2>::nsny; ++sy){
							for(T2 sx=0; sx < Grid3Dr<T1,T2>::nsnx; ++sx, n++){
								
								T1 ysv= Grid3Dr<T1,T2>::ymin+ nj* Grid3Dr<T1,T2>::dy+ (sy+1)*dys;
								T1 xsv= Grid3Dr<T1,T2>::xmin+ ni* Grid3Dr<T1,T2>::dx+ (sx+1)*dxs;
								
								if ( cell_XpYpZm != std::numeric_limits<T2>::max() ) {
									nodes[n].pushOwner( cell_XpYpZm );
								}
								if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
									nodes[n].pushOwner( cell_XpYpZp );
								}
								nodes[n].setXYZindex( xsv, ysv, z, n );
								
								if (nk>0){
									nodes[n].setPrimary(51);
								}
								else if (nk==0){
									nodes[n].setPrimary(50);
								}
							}
						}
					}
					
					// Secondary nodes on the x0z planes
					if ( ni < Grid3Dr<T1,T2>::ncx && nk < Grid3Dr<T1,T2>::ncz ) {
						for(T2 sz=0; sz < Grid3Dr<T1,T2>::nsnz; ++sz){
							for(T2 sx=0; sx < Grid3Dr<T1,T2>::nsnx; ++sx, n++){
								
								T1 zsv= Grid3Dr<T1,T2>::zmin+ nk* Grid3Dr<T1,T2>::dz+ (sz+1)*dzs;
								T1 xsv= Grid3Dr<T1,T2>::xmin+ ni* Grid3Dr<T1,T2>::dx+ (sx+1)*dxs;
								
								if ( cell_XpYmZp != std::numeric_limits<T2>::max() ) {
									nodes[n].pushOwner( cell_XpYmZp );
								}
								if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
									nodes[n].pushOwner( cell_XpYpZp );
								}
								nodes[n].setXYZindex( xsv, y, zsv, n );
								
								if (nj>0){
									nodes[n].setPrimary(61);
								}
								else if (nj==0){
									nodes[n].setPrimary(60);
								}
							}
						}
					}
					
					// Secondary nodes on the 0yz planes
					if ( nj < Grid3Dr<T1,T2>::ncy && nk < Grid3Dr<T1,T2>::ncz ) {
						for(T2 sz=0; sz < Grid3Dr<T1,T2>::nsnz; ++sz){
							for(T2 sy=0; sy < Grid3Dr<T1,T2>::nsny; ++sy, n++){
								
								T1 zsv= Grid3Dr<T1,T2>::zmin + nk*Grid3Dr<T1,T2>::dz + (sz+1)*dzs;
								T1 ysv= Grid3Dr<T1,T2>::ymin + nj*Grid3Dr<T1,T2>::dy + (sy+1)*dys;
								
								if ( cell_XmYpZp != std::numeric_limits<T2>::max() ) {
									nodes[n].pushOwner( cell_XmYpZp );
								}
								if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
									nodes[n].pushOwner( cell_XpYpZp );
								}
								nodes[n].setXYZindex( x, ysv, zsv, n );
								
								if (ni>0){
									nodes[n].setPrimary(71);
								}
								else if (ni==0){
									nodes[n].setPrimary(70);
								}
							}
						}
					}
				}
			}
		}
		// sanity check
		if ( n != nodes.size() ) {
			std::cerr << "Error building grid, wrong number of nodes\n";
			abort();
		}
	}
	
	template<typename T1, typename T2>
	void Grid3Dri<T1,T2>::buildGridNeighbors() {
		
		//Index the neighbors nodes of each cell
		for ( T2 n=0; n<nodes.size(); ++n ) {
			for ( size_t n2=0; n2<nodes[n].getOwners().size(); ++n2) {
				neighbors[ nodes[n].getOwners()[n2] ].push_back(n);
			}
		}
	}
	
	template<typename T1, typename T2>
	int Grid3Dri<T1,T2>::setSlowness(const std::vector<T1>& s) {
		
		if ( ((Grid3Dr<T1,T2>::ncx+1)*(Grid3Dr<T1,T2>::ncy+1)*(Grid3Dr<T1,T2>::ncz+1)) != s.size() ) {
			std::cerr << "Error: slowness vector of incompatible size.";
			return 1;
		}
		//Set the slowness for primary nodes
		size_t i=0;
		for ( size_t n=0; n<nodes.size(); ++n ) {
			if (nodes[n].getPrimary() == 5){
				nodes[n].setNodeSlowness( s[i] );
				i++;
			}
		}
		
        if ( inverseDistance ) {
            return invDistInterpolation();
        } else {
            return linearInterpolation();
        }
    }
    
	template<typename T1, typename T2>
	int Grid3Dri<T1,T2>::linearInterpolation() const {
        
        std::vector<size_t> list;
        list.reserve(8);
        T1 x[3], y[3], s[4];
        
		//Interpolation for the secondary nodes
		for (size_t n=0; n<nodes.size(); ++n ){
			T2 ind = nodes[n].getPrimary();
			if( ind != 5 ){
				
				//Lists the primary nodes around the secondary node
				list.resize(0);
				T2 cellno = nodes[n].getOwners()[0];
				for (size_t n3=0; n3 < neighbors[ cellno ].size(); n3++){
					if( nodes[neighbors[ cellno ][n3] ].getPrimary() == 5 ){
						list.push_back(neighbors[ cellno ][n3]);
					}
				}
				
				// list elements are as following:
				//
				// list[0] = x_min, y_min, z_min
				// list[1] = x_max, y_min, z_min
				// list[2] = x_min, y_max, z_min
				// list[3] = x_max, y_max, z_min
				// list[4] = x_min, y_min, z_max
				// list[5] = x_max, y_min, z_max
				// list[6] = x_min, y_max, z_max
				// list[7] = x_max, y_max, z_max
				
                switch (ind){
					case 0:
						std::cerr << "Error: the node " << n
						<< " is not set, check the grid construction.\n";
						return 1;
					case 25:
					case 27:
					case 45:
					case 47:
					case 60:
                        x[0] = nodes[n].getX();
                        y[0] = nodes[n].getZ();
                        
                        x[1] = nodes[ list[0] ].getX();
                        x[2] = nodes[ list[1] ].getX();
                        
                        y[1] = nodes[ list[0] ].getZ();
                        y[2] = nodes[ list[4] ].getZ();
                        
                        s[0] = nodes[ list[0] ].getNodeSlowness();
                        s[1] = nodes[ list[4] ].getNodeSlowness();
                        s[2] = nodes[ list[1] ].getNodeSlowness();
                        s[3] = nodes[ list[5] ].getNodeSlowness();
                        
                        break;
						
                    case 35:
                    case 36:
                    case 46:
                    case 70:
                        x[0] = nodes[n].getY();
                        y[0] = nodes[n].getZ();
                        
						x[1] = nodes[ list[0] ].getY();
                        x[2] = nodes[ list[2] ].getY();
                        
                        y[1] = nodes[ list[0] ].getZ();
                        y[2] = nodes[ list[4] ].getZ();
						
						s[0] = nodes[ list[0] ].getNodeSlowness();
                        s[1] = nodes[ list[4] ].getNodeSlowness();
                        s[2] = nodes[ list[2] ].getNodeSlowness();
                        s[3] = nodes[ list[6] ].getNodeSlowness();
                        
						break;
                        
                    case 37:
                    case 38:
                    case 48:
					case 71:
						x[0] = nodes[n].getY();
                        y[0] = nodes[n].getZ();
                        
						x[1] = nodes[ list[1] ].getY();
                        x[2] = nodes[ list[3] ].getY();
                        
                        y[1] = nodes[ list[1] ].getZ();
                        y[2] = nodes[ list[5] ].getZ();
						
						s[0] = nodes[ list[1] ].getNodeSlowness();
                        s[1] = nodes[ list[5] ].getNodeSlowness();
                        s[2] = nodes[ list[3] ].getNodeSlowness();
                        s[3] = nodes[ list[7] ].getNodeSlowness();
                        
                        break;
                        
					case 26:
					case 50:
						x[0] = nodes[n].getX();
                        y[0] = nodes[n].getY();
						
						x[1] = nodes[ list[0] ].getX();
                        x[2] = nodes[ list[1] ].getX();
                        
						y[1] = nodes[ list[0] ].getY();
                        y[2] = nodes[ list[2] ].getY();
						
						s[0] = nodes[ list[0] ].getNodeSlowness();
                        s[1] = nodes[ list[2] ].getNodeSlowness();
                        s[2] = nodes[ list[1] ].getNodeSlowness();
                        s[3] = nodes[ list[3] ].getNodeSlowness();
                        
						break;
						
					case 28:
					case 51:
						x[0] = nodes[n].getX();
                        y[0] = nodes[n].getY();
						
						x[1] = nodes[ list[4] ].getX();
                        x[2] = nodes[ list[5] ].getX();
						
						y[1] = nodes[ list[4] ].getY();
                        y[2] = nodes[ list[6] ].getY();
						
						s[0] = nodes[ list[4] ].getNodeSlowness();
                        s[1] = nodes[ list[6] ].getNodeSlowness();
                        s[2] = nodes[ list[5] ].getNodeSlowness();
                        s[3] = nodes[ list[7] ].getNodeSlowness();
                        
						break;
						
					case 61:
						x[0] = nodes[n].getX();
                        y[0] = nodes[n].getZ();
						
						x[1] = nodes[ list[2] ].getX();
                        x[2] = nodes[ list[3] ].getX();
						
						y[1] = nodes[ list[2] ].getZ();
                        y[2] = nodes[ list[6] ].getZ();
						
						s[0] = nodes[ list[2] ].getNodeSlowness();
                        s[1] = nodes[ list[6] ].getNodeSlowness();
                        s[2] = nodes[ list[3] ].getNodeSlowness();
                        s[3] = nodes[ list[7] ].getNodeSlowness();
						
						break;
                }
                
                nodes[n].setNodeSlowness( Interpolator<T1>::bilinear(x, y, s) );
			}
		}
		return 0;
	}
	
    template<typename T1, typename T2>
    int Grid3Dri<T1,T2>::invDistInterpolation() const {
        
        std::vector<size_t>::iterator it;
        std::vector<size_t> list;
        
        std::vector<Node3Di<T1,T2>*> interpNodes;
        interpNodes.reserve(18);
        
        //Interpolation for the secondary nodes
		for (size_t n=0; n<nodes.size(); ++n ){
			T2 ind = nodes[n].getPrimary();
			if ( ind != 5 ) {
				
				//Lists the primary nodes around the secondary node
				list.resize(0);
                interpNodes.resize(0);
                for ( size_t n2=0; n2<nodes[n].getOwners().size(); ++n2) {
                    T2 cellno = nodes[n].getOwners()[n2];
                    for (size_t n3=0; n3 < neighbors[ cellno ].size(); n3++){
                        if ( nodes[neighbors[ cellno ][n3] ].getPrimary() == 5 ) {
                            list.push_back(neighbors[ cellno ][n3]);
                        }
                    }
                }
                sort( list.begin(), list.end() );
                it = unique( list.begin(), list.end() );
                list.resize( it - list.begin() );
                
                for ( size_t nn=0; nn<list.size(); ++nn )
                    interpNodes.push_back( &(nodes[list[nn] ]) );
                
                nodes[n].setNodeSlowness( Interpolator<T1>::inverseDistance( nodes[n], interpNodes ) );
                
            }
        }
        return 0;
	}
	
	// Linear interpolation functions for each axis
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::linX(const Node3Di<T1,T2>& node, const size_t C1, const size_t C2) const {
		return nodes[C1].getNodeSlowness()+((node.getX()- nodes[C1].getX())/(nodes[C2].getX()- nodes[C1].getX()) )*(nodes[C2].getNodeSlowness()- nodes[C1].getNodeSlowness());
	}
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::linX(const Node3Di<T1,T2>& node, const size_t C1, const size_t C2, const T1 S1, const T1 S2) const {
		return S1+((node.getX()- nodes[C1].getX())/(nodes[C2].getX()- nodes[C1].getX()) )*(S2-S1);
	}
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::linX(const sxyz<T1>& node, const size_t C1, const size_t C2) const {
		return nodes[C1].getNodeSlowness()+((node.x- nodes[C1].getX())/(nodes[C2].getX()- nodes[C1].getX()) )*(nodes[C2].getNodeSlowness()- nodes[C1].getNodeSlowness() );
	}
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::linX(const sxyz<T1>& node, const size_t C1, const size_t C2, const T1 S1, const T1 S2) const {
		return S1+((node.x- nodes[C1].getX())/(nodes[C2].getX()- nodes[C1].getX()) )*(S2-S1 );
	}
	
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::linY(const Node3Di<T1,T2>& node, const size_t C1, const size_t C2) const {
		return nodes[C1].getNodeSlowness()+((node.getY()- nodes[C1].getY())/(nodes[C2].getY()- nodes[C1].getY()) )*(nodes[C2].getNodeSlowness()- nodes[C1].getNodeSlowness() );
	}
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::linY(const Node3Di<T1,T2>& node, const size_t C1, const size_t C2, const T1 S1, const T1 S2) const {
		return S1+((node.getY()- nodes[C1].getY())/(nodes[C2].getY()- nodes[C1].getY()) )*(S2-S1);
	}
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::linY(const sxyz<T1>& node, const size_t C1, const size_t C2) const {
		return nodes[C1].getNodeSlowness()+((node.y- nodes[C1].getY())/(nodes[C2].getY()- nodes[C1].getY()) )*(nodes[C2].getNodeSlowness()- nodes[C1].getNodeSlowness() );
	}
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::linY(const sxyz<T1>& node, const size_t C1, const size_t C2, const T1 S1, const T1 S2) const {
		return S1+((node.y- nodes[C1].getY())/(nodes[C2].getY()- nodes[C1].getY()) )*(S2-S1 );
	}
	
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::linZ(const Node3Di<T1,T2>& node, const size_t C1, const size_t C2) const {
		return nodes[C1].getNodeSlowness()+((node.getZ()- nodes[C1].getZ())/(nodes[C2].getZ()- nodes[C1].getZ()) )*(nodes[C2].getNodeSlowness()- nodes[C1].getNodeSlowness() );
	}
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::linZ(const Node3Di<T1,T2>& node, const size_t C1, const size_t C2, const T1 S1, const T1 S2) const {
		return S1+((node.getZ()- nodes[C1].getZ())/(nodes[C2].getZ()- nodes[C1].getZ()) )*(S2-S1);
	}
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::linZ(const sxyz<T1>& node, const size_t C1, const size_t C2) const {
		return nodes[C1].getNodeSlowness()+((node.z- nodes[C1].getZ())/(nodes[C2].getZ()- nodes[C1].getZ()) )*(nodes[C2].getNodeSlowness()- nodes[C1].getNodeSlowness() );
	}
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::linZ(const sxyz<T1>& node, const size_t C1, const size_t C2,
							 const T1 S1, const T1 S2) const{
		return S1+((node.z- nodes[C1].getZ())/(nodes[C2].getZ()- nodes[C1].getZ()) )*( S2-S1 );
	}
	
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::linXY(const Node3Di<T1,T2>& node, const size_t N1, const size_t N2,
							  const size_t N3, const size_t N4) const {
		T1 iA1 = 1./fabs((node.getX()-nodes[N1].getX())*(node.getY()-nodes[N1].getY()));
		T1 iA2 = 1./fabs((node.getX()-nodes[N2].getX())*(node.getY()-nodes[N2].getY()));
		T1 iA3 = 1./fabs((node.getX()-nodes[N3].getX())*(node.getY()-nodes[N3].getY()));
		T1 iA4 = 1./fabs((node.getX()-nodes[N4].getX())*(node.getY()-nodes[N4].getY()));
		
		return (iA1*nodes[N1].getNodeSlowness() + iA2*nodes[N2].getNodeSlowness() +
				iA3*nodes[N3].getNodeSlowness() + iA4*nodes[N4].getNodeSlowness() ) /
		(iA1+iA2+iA3+iA4);
	}
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::linXY(const sxyz<T1>& node, const size_t N1, const size_t N2,
							  const size_t N3, const size_t N4) const {
		T1 iA1 = 1./fabs( ( node.x-nodes[N1].getX() )*( node.y-nodes[N1].getY() ) );
		T1 iA2 = 1./fabs( ( node.x-nodes[N2].getX() )*( node.y-nodes[N2].getY() ) );
		T1 iA3 = 1./fabs( ( node.x-nodes[N3].getX() )*( node.y-nodes[N3].getY() ) );
		T1 iA4 = 1./fabs( ( node.x-nodes[N4].getX() )*( node.y-nodes[N4].getY() ) );
		
		return (iA1*nodes[N1].getNodeSlowness() + iA2*nodes[N2].getNodeSlowness() +
				iA3*nodes[N3].getNodeSlowness() + iA4*nodes[N4].getNodeSlowness() ) /
		(iA1+iA2+iA3+iA4);
	}
	
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::linXZ(const Node3Di<T1,T2> &node, const size_t N1, const size_t N2,
							  const size_t N3, const size_t N4) const {
		T1 iA1 = 1./fabs((node.getX()-nodes[N1].getX())*(node.getZ()-nodes[N1].getZ()));
		T1 iA2 = 1./fabs((node.getX()-nodes[N2].getX())*(node.getZ()-nodes[N2].getZ()));
		T1 iA3 = 1./fabs((node.getX()-nodes[N3].getX())*(node.getZ()-nodes[N3].getZ()));
		T1 iA4 = 1./fabs((node.getX()-nodes[N4].getX())*(node.getZ()-nodes[N4].getZ()));
		
		return (iA1*nodes[N1].getNodeSlowness() + iA2*nodes[N2].getNodeSlowness() +
				iA3*nodes[N3].getNodeSlowness() + iA4*nodes[N4].getNodeSlowness() ) /
		(iA1+iA2+iA3+iA4);
	}
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::linXZ(const sxyz<T1>& node, const size_t N1, const size_t N2,
							  const size_t N3, const size_t N4) const {
		T1 iA1 = 1./fabs( ( node.x-nodes[N1].getX() )*( node.z-nodes[N1].getZ() ) );
		T1 iA2 = 1./fabs( ( node.x-nodes[N2].getX() )*( node.z-nodes[N2].getZ() ) );
		T1 iA3 = 1./fabs( ( node.x-nodes[N3].getX() )*( node.z-nodes[N3].getZ() ) );
		T1 iA4 = 1./fabs( ( node.x-nodes[N4].getX() )*( node.z-nodes[N4].getZ() ) );
		
		return (iA1*nodes[N1].getNodeSlowness() + iA2*nodes[N2].getNodeSlowness() +
				iA3*nodes[N3].getNodeSlowness() + iA4*nodes[N4].getNodeSlowness() ) /
		(iA1+iA2+iA3+iA4);
	}
	
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::linYZ(const Node3Di<T1,T2>& node, const size_t N1, const size_t N2,
							  const size_t N3, const size_t N4) const {
		T1 iA1 = 1./fabs((node.getY()-nodes[N1].getY())*(node.getZ()-nodes[N1].getZ()));
		T1 iA2 = 1./fabs((node.getY()-nodes[N2].getY())*(node.getZ()-nodes[N2].getZ()));
		T1 iA3 = 1./fabs((node.getY()-nodes[N3].getY())*(node.getZ()-nodes[N3].getZ()));
		T1 iA4 = 1./fabs((node.getY()-nodes[N4].getY())*(node.getZ()-nodes[N4].getZ()));
		
		return (iA1*nodes[N1].getNodeSlowness() + iA2*nodes[N2].getNodeSlowness() +
				iA3*nodes[N3].getNodeSlowness() + iA4*nodes[N4].getNodeSlowness() ) /
		(iA1+iA2+iA3+iA4);
	}
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::linYZ(const sxyz<T1>& node, const size_t N1, const size_t N2,
							  const size_t N3, const size_t N4) const {
		T1 iA1 = 1./fabs( ( node.y-nodes[N1].getY() )*( node.z-nodes[N1].getZ() ) );
		T1 iA2 = 1./fabs( ( node.y-nodes[N2].getY() )*( node.z-nodes[N2].getZ() ) );
		T1 iA3 = 1./fabs( ( node.y-nodes[N3].getY() )*( node.z-nodes[N3].getZ() ) );
		T1 iA4 = 1./fabs( ( node.y-nodes[N4].getY() )*( node.z-nodes[N4].getZ() ) );
		
		return (iA1*nodes[N1].getNodeSlowness() + iA2*nodes[N2].getNodeSlowness() +
				iA3*nodes[N3].getNodeSlowness() + iA4*nodes[N4].getNodeSlowness() ) /
		(iA1+iA2+iA3+iA4);
	}
	
	
	template<typename T1, typename T2>
	void Grid3Dri<T1,T2>::initQueue(const std::vector<sxyz<T1> >& Tx,
									const std::vector<T1>& t0,
									std::priority_queue<Node3Di<T1,T2>*,
									std::vector<Node3Di<T1,T2>*>,
									CompareNodePtr<T1> >& queue,
									std::vector<Node3Di<T1,T2> >& txNodes,
									std::vector<bool>& inQueue,
									std::vector<bool>& frozen,
									const size_t threadNo) const {
		
		//Find the starting nodes of the transmitters Tx and start the queue list
		for(size_t n=0; n<Tx.size(); ++n){
			bool found = false;
			for ( size_t nn=0; nn<nodes.size(); ++nn ) {
				if ( nodes[nn] == Tx[n] ) {
					found = true;
					nodes[nn].setTT( t0[n], threadNo );
					frozen[nn] = true;
					
					prepropagate(nodes[nn], queue, inQueue, frozen, threadNo); // See description in the function declaration
					
					//	queue.push( &(nodes[nn]) );   	//Don't use if prepropagate is used
					//	inQueue[nn] = true;				//Don't use if prepropagate is used
					
					break;
				}
			}
			if ( found==false ) {
				// If Tx[n] is not on a node, we create a new node and initialize the queue:
				txNodes.push_back( Node3Di<T1,T2>(t0[n], Tx[n].x, Tx[n].y, Tx[n].z, Grid3Dr<T1,T2>::nThreads, threadNo));
				txNodes.back().pushOwner( this->getCellNo(Tx[n]) );
				txNodes.back().setGridIndex( static_cast<T2>(nodes.size()+txNodes.size()-1) );
				frozen.push_back( true );
				T1 slo = computeSlowness( Tx[n] );
				txNodes.back().setNodeSlowness( slo );
				
				prepropagate(txNodes.back(), queue, inQueue, frozen, threadNo); // See description in the function declaration
				
				//	queue.push( &(txNodes.back()) );	//Don't use if prepropagate is used
				//	inQueue.push_back( true );			//Don't use if prepropagate is used
				
			}
		}
	}
	
	template<typename T1, typename T2>
	int Grid3Dri<T1,T2>::raytrace(const std::vector<sxyz<T1> >& Tx,
								  const std::vector<T1>& t0,
								  const std::vector<sxyz<T1> >& Rx,
								  std::vector<T1>& traveltimes,
								  const size_t threadNo) const {
		
		// Primary function
		
		// Checks if the points are in the grid
		if ( this->check_pts(Tx) == 1 ) return 1;
		if ( this->check_pts(Rx) == 1 ) return 1;
		
		for ( size_t n=0; n<nodes.size(); ++n ) {
			nodes[n].reinit( threadNo );
		}
		
		CompareNodePtr<T1> cmp(threadNo);
		std::priority_queue< Node3Di<T1,T2>*, std::vector<Node3Di<T1,T2>*>,
		CompareNodePtr<T1> > queue(cmp);
		// txNodes: Extra nodes if the sources points are not on an existing node
		std::vector<Node3Di<T1,T2> > txNodes;
		// inQueue lists the nodes waiting in the queue
		std::vector<bool> inQueue( nodes.size(), false );
		// Tx sources nodes are "frozen" and their traveltime can't be modified
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
	int Grid3Dri<T1,T2>::raytrace(const std::vector<sxyz<T1> >& Tx,
								  const std::vector<T1>& t0,
								  const std::vector<sxyz<T1> >& Rx,
								  std::vector<T1>& traveltimes,
								  std::vector<std::vector<sxyz<T1> > >& r_data,
								  const size_t threadNo) const {
		
		// Primary function
		
		// Checks if the points are in the grid
		if ( this->check_pts(Tx) == 1 ) return 1;
		if ( this->check_pts(Rx) == 1 ) return 1;
		
		for ( size_t n=0; n<nodes.size(); ++n ) {
			nodes[n].reinit( threadNo );
		}
		
		CompareNodePtr<T1> cmp(threadNo);
		std::priority_queue< Node3Di<T1,T2>*, std::vector<Node3Di<T1,T2>*>,
		CompareNodePtr<T1> > queue(cmp);
		// txNodes: Extra nodes if the sources points are not on an existing node
		std::vector<Node3Di<T1,T2> > txNodes;
		// inQueue lists the nodes waiting in the queue
		std::vector<bool> inQueue( nodes.size(), false );
		// Tx sources nodes are "frozen" and their traveltime can't be modified
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
			std::vector<Node3Di<T1,T2> > *node_p;
			node_p = &nodes;
			
			std::vector<sxyz<T1> > r_tmp;
			T2 iChild, iParent = nodeParentRx;
			sxyz<T1> child;
			//        siv<T1> cell;
			
			// store the son's coord
			child.x = Rx[n].x;
			child.y = Rx[n].y;
			child.z = Rx[n].z;
			//        cell.i = cellParentRx;
			while ( (*node_p)[iParent].getNodeParent(threadNo) != std::numeric_limits<T2>::max() ) {
				
				r_tmp.push_back( child );
				
				//			cell.v = (*node_p)[iParent].getDistance( child );
				//			bool found=false;
				//			for (size_t nc=0; nc<l_data[n].size(); ++nc) {
				//				if ( l_data[n][nc].i == cell.i ) {
				//					l_data[n][nc].v += cell.v;
				//					found = true;
				//				}
				//			}
				//			if ( found == false ) {
				//				l_data[n].push_back( cell );
				//			}
				
				// we now go up in time - parent becomes the child of grand'pa
				iChild = iParent;
				child.x = (*node_p)[iChild].getX();
				child.y = (*node_p)[iChild].getY();
				child.z = (*node_p)[iChild].getZ();
				//            cell.i = (*node_p)[iChild].getCellParent();
				
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
			
			//		cell.v = (*node_p)[iParent].getDistance( child );
			//		bool found=false;
			//		for (size_t nc=0; nc<l_data[n].size(); ++nc) {
			//			if ( l_data[n][nc].i == cell.i ) {
			//				l_data[n][nc].v += cell.v;
			//				found = true;
			//			}
			//		}
			//		if ( found == false ) {
			//			l_data[n].push_back( cell );
			//		}
			
			// finally, store Tx position
			child.x = (*node_p)[iParent].getX();
			child.y = (*node_p)[iParent].getY();
			child.z = (*node_p)[iParent].getZ();
			r_tmp.push_back( child );
			
			//  must be sorted to build matrix L
			//        sort(l_data[n].begin(), l_data[n].end(), CompareSiv_i<T1>());
			
			// the order should be from Tx to Rx, so we reorder...
			iParent = static_cast<T2>(r_tmp.size());
			r_data[n].resize( r_tmp.size() );
			for ( size_t nn=0; nn<r_data[n].size(); ++nn ) {
				r_data[n][nn].x = r_tmp[ iParent-1-nn ].x;
				r_data[n][nn].y = r_tmp[ iParent-1-nn ].y;
				r_data[n][nn].z = r_tmp[ iParent-1-nn ].z;
			}
			
		}
		return 0;
	}
	
	template<typename T1, typename T2>
	void Grid3Dri<T1,T2>::propagate( std::priority_queue<Node3Di<T1,T2>*,
									std::vector<Node3Di<T1,T2>*>,
									CompareNodePtr<T1> >& queue,
									std::vector<bool>& inQueue,
									std::vector<bool>& frozen,
									size_t threadNo) const {
		
		while ( !queue.empty() ) {
			const Node3Di<T1,T2>* source = queue.top();
			queue.pop();
			inQueue[ source->getGridIndex() ] = false;
			for ( size_t no=0; no<source->getOwners().size(); ++no ) {
				T2 cellNo = source->getOwners()[no];
				for ( size_t k=0; k< neighbors[cellNo].size(); ++k ) {
					T2 neibNo = neighbors[cellNo][k];
					if ( neibNo == source->getGridIndex() || frozen[neibNo] ) {
						continue;
					}
					
					T1 ttsource= source->getTT( threadNo );
					if (ttsource < nodes[neibNo].getTT(threadNo)){
						// Compute dt
						T1 dt = computeDt(*source, nodes[neibNo]);
						
						if ( ttsource +dt < nodes[neibNo].getTT( threadNo ) ) {
							nodes[neibNo].setTT( ttsource +dt, threadNo );
							nodes[neibNo].setnodeParent( source->getGridIndex(),
														threadNo );
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
	}
	
	template<typename T1, typename T2>
	void Grid3Dri<T1,T2>::prepropagate(const Node3Di<T1,T2>& node,
									   std::priority_queue<Node3Di<T1,T2>*,
									   std::vector<Node3Di<T1,T2>*>,
									   CompareNodePtr<T1> >& queue,
									   std::vector<bool>& inQueue,
									   std::vector<bool>& frozen,
									   const size_t threadNo) const {
		
		// This function can be used to "prepropagate" each Tx nodes one first time
		// during "initQueue", before running "propagate".
		// When a Tx source node seems to be lost in the queue and is not
		// propagated, corrupting the entire traveltime table,
		// this function force the propagation of every source points and can
		// solve the problem.
		
		for ( size_t no=0; no<node.getOwners().size(); ++no ) {
			T2 cellNo = node.getOwners()[no];
			for ( size_t k=0; k< neighbors[cellNo].size(); ++k ) {
				T2 neibNo = neighbors[cellNo][k];
				if ( neibNo == node.getGridIndex() || frozen[neibNo] ) {
					continue;
				}
				
				// compute dt
				T1 dt = computeDt(node, nodes[neibNo]);
				
				if ( node.getTT( threadNo )+dt < nodes[neibNo].getTT( threadNo ) ) {
					nodes[neibNo].setTT( node.getTT( threadNo )+dt, threadNo );
					nodes[neibNo].setnodeParent( node.getGridIndex(), threadNo );
					nodes[neibNo].setCellParent( cellNo, threadNo );
					
					if ( !inQueue[neibNo] ) {
						queue.push( &(nodes[neibNo]) );
						inQueue[neibNo] = true;
					}
				}
			}
		}
	}
	
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::getTraveltime(const sxyz<T1>& Rx,
									  const std::vector<Node3Di<T1,T2> >& nodes,
									  const size_t threadNo) const {
		
		// Calculate and return the traveltime for a Rx point.
		
		// If Rx is on a node:
		for ( size_t nn=0; nn<nodes.size(); ++nn ) {
			if ( nodes[nn] == Rx ) {
				return nodes[nn].getTT(threadNo);
			}
		}
		//If Rx is not on a node:
		T1 slo = computeSlowness( Rx );
		
		T2 cellNo = this->getCellNo( Rx );
		T2 neibNo = neighbors[cellNo][0];
		T1 dt = computeDt(nodes[neibNo], Rx, slo);
		
		T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
		for ( size_t k=1; k< neighbors[cellNo].size(); ++k ) {
			neibNo = neighbors[cellNo][k];
			dt = computeDt(nodes[neibNo], Rx, slo);
			if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
				traveltime =  nodes[neibNo].getTT(threadNo)+dt;
			}
		}
		return traveltime;
	}
	
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::getTraveltime(const sxyz<T1>& Rx,
									  const std::vector<Node3Di<T1,T2> >& nodes,
									  T2& nodeParentRx, T2& cellParentRx,
									  const size_t threadNo) const {
		
		// Calculate and return the traveltime for a Rx point.
		for ( size_t nn=0; nn<nodes.size(); ++nn ) {
			if ( nodes[nn] == Rx ) {
				nodeParentRx = nodes[nn].getNodeParent(threadNo);
				cellParentRx = nodes[nn].getCellParent(threadNo);
				return nodes[nn].getTT(threadNo);
			}
		}
		//If Rx is not on a node:
		T1 slo = computeSlowness( Rx );
		
		T2 cellNo = this->getCellNo( Rx );
		T2 neibNo = neighbors[cellNo][0];
		T1 dt = computeDt(nodes[neibNo], Rx, slo);
		
		T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
		nodeParentRx = neibNo;
		cellParentRx = cellNo;
		for ( size_t k=1; k< neighbors[cellNo].size(); ++k ) {
			neibNo = neighbors[cellNo][k];
			dt = computeDt(nodes[neibNo], Rx, slo);
			if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
				traveltime =  nodes[neibNo].getTT(threadNo)+dt;
				nodeParentRx = neibNo;
			}
		}
		return traveltime;
	}
	
	template<typename T1, typename T2>
	T1 Grid3Dri<T1,T2>::computeSlowness(const sxyz<T1>& Rx) const {

		
        //Calculate the slowness of any point that is not on a node
		
		T2 cellNo = this->getCellNo( Rx );
		
		//We calculate the Slowness at the point
		std::vector<T2> list;
		
		for (size_t n3=0; n3 < neighbors[ cellNo ].size(); n3++){
			if ( nodes[neighbors[ cellNo ][n3] ].getPrimary() == 5 ){
				list.push_back(neighbors[ cellNo ][n3]);
			}
		}
		
        if ( inverseDistance ) {
            
            std::vector<size_t>::iterator it;
            
            std::vector<Node3Di<T1,T2>*> interpNodes;
            
            for ( size_t nn=0; nn<list.size(); ++nn )
                interpNodes.push_back( &(nodes[list[nn] ]) );
            
            return Interpolator<T1>::inverseDistance( Rx, interpNodes );

        } else {

            // list elements are as following:
            //
            // list[0] = x_min, y_min, z_min
            // list[1] = x_max, y_min, z_min
            // list[2] = x_min, y_max, z_min
            // list[3] = x_max, y_max, z_min
            // list[4] = x_min, y_min, z_max
            // list[5] = x_max, y_min, z_max
            // list[6] = x_min, y_max, z_max
            // list[7] = x_max, y_max, z_max
            
            T1 x[3] = { Rx.x, nodes[list[0]].getX(), nodes[list[1]].getX() };
            T1 y[3] = { Rx.y, nodes[list[0]].getY(), nodes[list[2]].getY() };
            T1 z[3] = { Rx.z, nodes[list[0]].getZ(), nodes[list[4]].getZ() };
        
            T1 s[8] = { nodes[list[0]].getNodeSlowness(),
                nodes[list[4]].getNodeSlowness(),
                nodes[list[2]].getNodeSlowness(),
                nodes[list[6]].getNodeSlowness(),
                nodes[list[1]].getNodeSlowness(),
                nodes[list[5]].getNodeSlowness(),
                nodes[list[3]].getNodeSlowness(),
                nodes[list[7]].getNodeSlowness() };
        
            return Interpolator<T1>::trilinear(x, y, z, s);
        }
    }
	
	
	template<typename T1, typename T2>
	void Grid3Dri<T1,T2>::save(const char filename[]) const {
		std::ofstream fout( filename );
		
		fout << Grid3Dr<T1,T2>::dx << ' ' << Grid3Dr<T1,T2>::dy << ' ' << Grid3Dr<T1,T2>::dz << ' ' << Grid3Dr<T1,T2>::xmin << ' ' << Grid3Dr<T1,T2>::ymin << ' '
		<< Grid3Dr<T1,T2>::zmin << ' ' << Grid3Dr<T1,T2>::xmax << ' ' << Grid3Dr<T1,T2>::ymax << ' '<< Grid3Dr<T1,T2>::zmax << '\n';
		fout << Grid3Dr<T1,T2>::ncx << ' ' << Grid3Dr<T1,T2>::ncy << ' ' << Grid3Dr<T1,T2>::ncz << ' '
		<< Grid3Dr<T1,T2>::nsnx << ' ' << Grid3Dr<T1,T2>::nsny << ' ' << Grid3Dr<T1,T2>::nsnz << ' ' << '\n';
		
		fout << nodes.size() << '\n';
		for ( size_t n=0; n < nodes.size(); ++n ) {
			fout << nodes[n].getsize() ;
			for (size_t nt=0; nt< nodes[n].getsize(); nt++){
				fout << " " << nodes[n].getTT(nt) << " "
				<< nodes[n].getNodeParent(nt) << ' '<< nodes[n].getCellParent(nt);
			}
			fout << ' ' << nodes[n].getX() << ' ' << nodes[n].getY()
			<< ' ' << nodes[n].getZ() << ' '
			<< ' ' << nodes[n].getGridIndex();
			for (size_t no=0; no < nodes[n].getOwners().size(); ++no ) {
				fout << ' ' << nodes[n].getOwners()[no];
			}
			fout << '\n';
		}
		/*
		 fout << slowness.size() << '\n';
		 for ( size_t n=0; n < slowness.size(); ++n ) {
		 fout << slowness[n] << '\n';
		 }
		 fout << neighbors.size() << '\n';
		 for ( size_t n=0; n < neighbors.size(); ++n ) {
		 fout << neighbors[n].size();
		 for ( size_t nn=0; nn < neighbors[n].size(); ++nn ) {
		 fout << ' ' << neighbors[n][nn];
		 }
		 fout << '\n';
		 }
		 */
		fout.close();
	}
	
	template<typename T1, typename T2>
	void Grid3Dri<T1,T2>::dsave(const char filename[]) const {
		//Similar to 'save', with text information
		std::ofstream fout( filename );
		
		fout << "dx "<< Grid3Dr<T1,T2>::dx << "\t dy " << Grid3Dr<T1,T2>::dy << "\t dz " << Grid3Dr<T1,T2>::dz << "\t xmin " << Grid3Dr<T1,T2>::xmin
		<< "\t ymin " << Grid3Dr<T1,T2>::ymin << "\t zmin "	<< Grid3Dr<T1,T2>::zmin << "\t xmax " << Grid3Dr<T1,T2>::xmax
		<< "\t ymax " << Grid3Dr<T1,T2>::ymax << "\t zmax "<< Grid3Dr<T1,T2>::zmax << '\n';
		fout << "nCx " << Grid3Dr<T1,T2>::ncx << "\t nCy " << Grid3Dr<T1,T2>::ncy << "\t nCz " << Grid3Dr<T1,T2>::ncz << "\t nsnx "
		<< Grid3Dr<T1,T2>::nsnx << "\t nsny " << Grid3Dr<T1,T2>::nsny << "\t nsnz " << Grid3Dr<T1,T2>::nsnz << '\n';
		
		fout << "nb. nodes " << nodes.size() << '\n';
		for ( size_t n=0; n < nodes.size(); ++n ) {
			fout << "node " << nodes[n].getGridIndex() << "\t TT \t ";
			for (size_t nt=0; nt< nodes[n].getsize(); nt++){
				fout << nodes[n].getTT(nt) << "\t";
			}
			fout << " X " << nodes[n].getX() << "\t Y " << nodes[n].getY()
			<< "\t Z " << nodes[n].getZ() << "\t Slowness " << nodes[n].getNodeSlowness()<<
			"\t Ray Parent \t";
			for (size_t nt=0; nt< nodes[n].getsize(); nt++){
				fout << nodes[n].getNodeParent(nt) << '\t';
			}
			fout<< "Cell Parent \t";
			for (size_t nt=0; nt< nodes[n].getsize(); nt++){
				fout << nodes[n].getCellParent(nt) << '\t';
			}
			fout << "Owners: ";
			for (size_t no=0; no < nodes[n].getOwners().size(); ++no ) {
				fout << '\t' << nodes[n].getOwners()[no];
			}
			fout << '\n';
		}
		/*
		 fout << "neighbors size " << neighbors.size() << '\n';
		 for ( size_t n=0; n < neighbors.size(); ++n ) {
		 fout << "neighbors[" << n << "] size " << neighbors[n].size() << " :";
		 for ( size_t nn=0; nn < neighbors[n].size(); ++nn ) {
		 fout << '\t' << neighbors[n][nn];
		 }
		 fout << '\n';
		 }
		 */
		
		fout.close();
	}
	
	template<typename T1, typename T2>
	void Grid3Dri<T1,T2>::savefast(const char filename[]) const {
		std::ofstream fout( filename );
		
		for ( size_t n=0; n < nodes.size(); ++n ) {
			if ( floor((nodes[n].getX())/Grid3Dr<T1,T2>::dx)==(nodes[n].getX())/Grid3Dr<T1,T2>::dx &&
				floor((nodes[n].getZ())/Grid3Dr<T1,T2>::dz) == (nodes[n].getZ())/Grid3Dr<T1,T2>::dz &&
				floor((nodes[n].getY())/Grid3Dr<T1,T2>::dy) == (nodes[n].getY())/Grid3Dr<T1,T2>::dy )
			{
				//		fout <<  ((nodes[n].getX())/dx)+1 << '\t' << ((nodes[n].getY())/dy)+1
				//	    << '\t' << ((nodes[n].getZ())/dz)+1 ;
				for ( size_t nt=0; nt< nodes[n].getsize(); nt++ ) {
					fout.precision(9);
					fout //<< '\t'
					<< nodes[n].getTT(nt);
				}
				fout << '\n';
			}
		}
		
		fout.close();
	}
	
	template<typename T1, typename T2>
	void Grid3Dri<T1,T2>::savePrimary(const char filename[], const size_t nt,
									  const bool vtkFormat) const {
		
		if ( vtkFormat ) {
			
#ifdef VTK
			vtkSmartPointer<vtkDoubleArray> xCoords = vtkSmartPointer<vtkDoubleArray>::New();
			for (size_t ni=0; ni<=Grid3Dr<T1,T2>::ncx; ++ni)
				xCoords->InsertNextValue(Grid3Dr<T1,T2>::xmin + ni*Grid3Dr<T1,T2>::dx);
			vtkSmartPointer<vtkDoubleArray> yCoords = vtkSmartPointer<vtkDoubleArray>::New();
			for (size_t nj=0; nj<=Grid3Dr<T1,T2>::ncy; ++nj)
				yCoords->InsertNextValue(Grid3Dr<T1,T2>::ymin + nj*Grid3Dr<T1,T2>::dy);
			vtkSmartPointer<vtkDoubleArray> zCoords = vtkSmartPointer<vtkDoubleArray>::New();
			for (size_t nk=0; nk<=Grid3Dr<T1,T2>::ncz; ++nk)
				zCoords->InsertNextValue(Grid3Dr<T1,T2>::zmin + nk*Grid3Dr<T1,T2>::dz);
			
			vtkSmartPointer<vtkRectilinearGrid> rgrid = vtkSmartPointer<vtkRectilinearGrid>::New();
			rgrid->SetDimensions(Grid3Dr<T1,T2>::ncx, Grid3Dr<T1,T2>::ncy, Grid3Dr<T1,T2>::ncz);
			rgrid->SetXCoordinates(xCoords);
			rgrid->SetYCoordinates(yCoords);
			rgrid->SetZCoordinates(zCoords);
			
			vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
			data->SetName("Travel time");
			size_t n=0;
			for ( size_t nk=0; nk<=Grid3Dr<T1,T2>::ncz; ++nk ) {
				for ( size_t nj=0; nj<=Grid3Dr<T1,T2>::ncy; ++nj ) {
					for ( size_t ni=0; ni<=Grid3Dr<T1,T2>::ncx; ++ni ) {
						
						data->InsertNextValue( nodes[n++].getTT(nt) );
						
						// Secondary nodes on x edge
						if ( ni < Grid3Dr<T1,T2>::ncx ) {
							n += Grid3Dr<T1,T2>::nsnx;
						}
						
						// Secondary nodes on y edge
						if ( nj < Grid3Dr<T1,T2>::ncy ) {
							n += Grid3Dr<T1,T2>::nsny;
						}
						
						// Secondary nodes on z edge
						if ( nk < Grid3Dr<T1,T2>::ncz ) {
							n += Grid3Dr<T1,T2>::nsnz;
						}
						
						// Secondary nodes on the xy0 planes
						if ( ni < Grid3Dr<T1,T2>::ncx && nj < Grid3Dr<T1,T2>::ncy ) {
							n += Grid3Dr<T1,T2>::nsny*Grid3Dr<T1,T2>::nsnx;
						}
						
						// Secondary nodes on the x0z planes
						if ( ni < Grid3Dr<T1,T2>::ncx && nk < Grid3Dr<T1,T2>::ncz ) {
							n += Grid3Dr<T1,T2>::nsnz*Grid3Dr<T1,T2>::nsnx;
						}
						
						// Secondary nodes on the 0yz planes
						if ( nj < Grid3Dr<T1,T2>::ncy && nk < Grid3Dr<T1,T2>::ncz ) {
							n += Grid3Dr<T1,T2>::nsnz*Grid3Dr<T1,T2>::nsny;
						}
					}
				}
			}
			rgrid->GetPointData()->SetScalars( data );
			
			vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
			
			writer->SetFileName( filename );
			writer->SetInput( rgrid );
			writer->SetDataModeToBinary();
			writer->Update();
#else
			std::cerr << "VTK not included during compilation.\nNothing saved.\n";
#endif
		} else {
		
		std::ofstream fout( filename );
		fout.precision(9);
		
		size_t n=0;
		for ( size_t nk=0; nk<=Grid3Dr<T1,T2>::ncz; ++nk ) {
			
			for ( size_t nj=0; nj<=Grid3Dr<T1,T2>::ncy; ++nj ) {
				
				for (size_t ni=0; ni<=Grid3Dr<T1,T2>::ncx; ++ni ) {
					
					fout << nodes[n++].getTT(nt) << '\n';
					
					// Secondary nodes on x edge
					if ( ni < Grid3Dr<T1,T2>::ncx ) {
						n += Grid3Dr<T1,T2>::nsnx;
					}
					
					// Secondary nodes on y edge
					if ( nj < Grid3Dr<T1,T2>::ncy ) {
						n += Grid3Dr<T1,T2>::nsny;
					}
					
					// Secondary nodes on z edge
					if ( nk < Grid3Dr<T1,T2>::ncz ) {
						n += Grid3Dr<T1,T2>::nsnz;
					}
					
					// Secondary nodes on the xy0 planes
					if ( ni < Grid3Dr<T1,T2>::ncx && nj < Grid3Dr<T1,T2>::ncy ) {
						n += Grid3Dr<T1,T2>::nsny*Grid3Dr<T1,T2>::nsnx;
					}
					
					// Secondary nodes on the x0z planes
					if ( ni < Grid3Dr<T1,T2>::ncx && nk < Grid3Dr<T1,T2>::ncz ) {
						n += Grid3Dr<T1,T2>::nsnz*Grid3Dr<T1,T2>::nsnx;
					}
					
					// Secondary nodes on the 0yz planes
					if ( nj < Grid3Dr<T1,T2>::ncy && nk < Grid3Dr<T1,T2>::ncz ) {
						n += Grid3Dr<T1,T2>::nsnz*Grid3Dr<T1,T2>::nsny;
					}
				}
			}
		}
		fout.close();
	}
	}
	
#endif

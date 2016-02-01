//
//  Grid3Dri.h
//  ttcr.v2
//
//  Created by Giroux Bernard on 12-08-15.
//  Copyright (c) 2012 INRS-ETE. All rights reserved.
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

#ifndef __GRID3DRI_H__
#define __GRID3DRI_H__


#include <algorithm>
#include <cstring>
#include <iostream>
#include <fstream>
#include <queue>
#include <vector>
#include <ctime>

#include <boost/math/special_functions/sign.hpp>

#ifdef VTK
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkXMLRectilinearGridWriter.h"
#endif

#include "Grid3D.h"
#include "Interpolator.h"

template<typename T1, typename T2, typename NODE>
class Grid3Dri : public Grid3D<T1,T2> {
public:
    
    /* Constructor Format:
     Grid3Dri<T1,T2>::Grid3Dri(nb cells in x, nb cells in y, nb cells in z,
     x cells size, y cells size, z cells size,
     x origin, y origin, z origin,
     nb sec. cells in x, nb sec. cells in y, nb sec. cells in z,
     index of the thread)
     */
    Grid3Dri(const T2 nx, const T2 ny, const T2 nz,
             const T1 ddx, const T1 ddy, const T1 ddz,
             const T1 minx, const T1 miny, const T1 minz,
             const size_t nt=1, const bool invDist=false) :
    nThreads(nt),
    dx(ddx), dy(ddy), dz(ddz),
    xmin(minx), ymin(miny), zmin(minz),
    xmax(minx+nx*ddx), ymax(miny+ny*ddy), zmax(minz+nz*ddz),
    ncx(nx), ncy(ny), ncz(nz),
    nodes(std::vector<NODE>((nx+1)*(ny+1)*(nz+1), NODE(nt))),
    neighbors(std::vector<std::vector<T2>>(nx*ny*nz))
    {    }
    
    virtual ~Grid3Dri() {}
    
    T1 getSlowness(const size_t n) const {
        return nodes[n].getNodeSlowness();
    }
    
    void setSlowness(const T1 s) {
        for ( size_t n=0; n<nodes.size(); ++n ) {
            nodes[n].setNodeSlowness( s );
        }
    }
    
    virtual int setSlowness(const std::vector<T1>& s) {
        if ( nodes.size() != s.size() ) {
            std::cerr << "Error: slowness vectors of incompatible size.";
            return 1;
        }
        for ( size_t n=0; n<nodes.size(); ++n ) {
            nodes[n].setNodeSlowness( s[n] );
        }
        return 0;
    }
    
    size_t getNumberOfNodes() const { return nodes.size(); }
    
    virtual int raytrace(const std::vector<sxyz<T1>>& Tx,
                         const std::vector<T1>& t0,
                         const std::vector<sxyz<T1>>& Rx,
                         std::vector<T1>& traveltimes,
                         const size_t threadNo=0) const { return 0; }
    
    virtual int raytrace(const std::vector<sxyz<T1>>& Tx,
                         const std::vector<T1>& t0,
                         const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                         std::vector<std::vector<T1>*>& traveltimes,
                         const size_t=0) const { return 0; }
    
    virtual int raytrace(const std::vector<sxyz<T1>>& Tx,
                         const std::vector<T1>& t0,
                         const std::vector<sxyz<T1>>& Rx,
                         std::vector<T1>& traveltimes,
                         std::vector<std::vector<sxyz<T1>>>& r_data,
                         const size_t threadNo=0) const { return 0; }
    
    virtual int raytrace(const std::vector<sxyz<T1>>& Tx,
                         const std::vector<T1>& t0,
                         const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                         std::vector<std::vector<T1>*>& traveltimes,
                         std::vector<std::vector<std::vector<sxyz<T1>>>*>& r_data,
                         const size_t=0) const { return 0; }
    
    void saveSlownessXYZ(const char filename[]) const {
        //Saves the Slowness of the primary nodes
        std::ofstream fout( filename );
        for ( size_t n=0; n< nodes.size(); ++n ) {
            if (nodes[n].isPrimary() == 5 ){
                fout << nodes[n].getX() << "   "
                << nodes[n].getY() << "   "
                << nodes[n].getZ() << "   "
                << nodes[n].getNodeSlowness() << '\n';
            }
        }
        fout.close();
    }
    
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
    
    void saveTT(const std::string &, const int, const size_t nt=0,
                const bool vtkFormat=0) const;

    virtual const int get_niter() const { return 0; }
    virtual const int get_niterw() const { return 0; }

    const size_t getNthreads() const { return nThreads; }
    const T1 getXmin() const { return xmin; }
    const T1 getYmin() const { return ymin; }
    const T1 getZmin() const { return zmin; }
    const T1 getDx() const { return dx; }
    const T1 getDy() const { return dy; }
    const T1 getDz() const { return dz; }
    const T2 getNcx() const { return ncx; }
    const T2 getNcy() const { return ncy; }
    const T2 getNcz() const { return ncz; }

protected:
    size_t nThreads;	     // number of threads
    T1 dx;                   // cell size in x
    T1 dy;			         // cell size in y
    T1 dz;                   // cell size in z
    T1 xmin;                 // x origin of the grid
    T1 ymin;                 // y origin of the grid
    T1 zmin;                 // z origin of the grid
    T1 xmax;                 // x end of the grid
    T1 ymax;                 // y end of the grid
    T1 zmax;                 // z end of the grid
    T2 ncx;                  // number of cells in x
    T2 ncy;                  // number of cells in y
    T2 ncz;                  // number of cells in z

    bool inverseDistance;
    
    mutable std::vector<NODE> nodes;
    std::vector<std::vector<T2>> neighbors;  // nodes common to a cell
    
    void buildGridNeighbors();
    
    T2 getCellNo(const sxyz<T1>& pt) const {
        T1 x = xmax-pt.x < small ? xmax-.5*dx : pt.x;
        T1 y = ymax-pt.y < small ? ymax-.5*dy : pt.y;
        T1 z = zmax-pt.z < small ? zmax-.5*dz : pt.z;
        T2 nx = static_cast<T2>( small + (x-xmin)/dx );
        T2 ny = static_cast<T2>( small + (y-ymin)/dy );
        T2 nz = static_cast<T2>( small + (z-zmin)/dz );
        return ny*ncx + nz*(ncx*ncy) + nx;
    }
    
    T2 getCellNo(const NODE& node) const {
        T1 x = xmax-node.getX() < small ? xmax-.5*dx : node.getX();
        T1 y = ymax-node.getY() < small ? ymax-.5*dy : node.getY();
        T1 z = zmax-node.getZ() < small ? zmax-.5*dz : node.getZ();
        T2 nx = static_cast<T2>( small + (x-xmin)/dx );
        T2 ny = static_cast<T2>( small + (y-ymin)/dy );
        T2 nz = static_cast<T2>( small + (z-zmin)/dz );
        return ny*ncx + nz*(ncx*ncy) + nx;
    }
    
    void getIJK(const sxyz<T1>& pt, T2& i, T2& j, T2& k) const {
        i = static_cast<T2>( small + (pt.x-xmin)/dx );
        j = static_cast<T2>( small + (pt.y-ymin)/dy );
        k = static_cast<T2>( small + (pt.z-zmin)/dz );
    }
    
    void getIJK(const sxyz<T1>& pt, long long& i, long long& j, long long& k) const {
        i = static_cast<long long>( small + (pt.x-xmin)/dx );
        j = static_cast<long long>( small + (pt.y-ymin)/dy );
        k = static_cast<long long>( small + (pt.z-zmin)/dz );
    }
    
    int checkPts(const std::vector<sxyz<T1>>&) const;

    T1 computeDt(const NODE& source, const NODE& node) const {
        return (node.getNodeSlowness()+source.getNodeSlowness())/2. * source.getDistance( node );
    }
    
    T1 computeDt(const NODE& source, const sxyz<T1>& node, T1 slo) const {
        return (slo+source.getNodeSlowness())/2. * source.getDistance( node );
    }
    
    bool isNearInt( double value ) const {
        return ( remainder(value, 1.)  <= small );
    }
    
    T1 getTraveltime(const sxyz<T1> &pt, const size_t nt) const;
    
//    T1 getTraveltime(const sxyz<T1>& Rx,
//                     const std::vector<NODE>& nodes,
//                     const size_t threadNo) const;
    
    T1 getTraveltime(const sxyz<T1>& Rx,
                     T2&, T2&, const size_t threadNo) const;

//    T1 getTraveltime(const sxyz<T1>& Rx,
//                     const std::vector<NODE>& nodes,
//                     T2&, T2&, const size_t threadNo) const;
    
    void grad(sxyz<T1>& g, const size_t i, const size_t j, const size_t k,
              const size_t nt) const;

    void grad(sxyz<T1>& g, const sxyz<T1> &pt, const size_t nt) const;
    
    void getRaypath(const std::vector<sxyz<T1>>& Tx,
                    const sxyz<T1> &Rx,
                    std::vector<sxyz<T1>> &r_data,
                    const size_t threadNo=0) const;

    void getRaypath_old(const std::vector<sxyz<T1>>& Tx,
                        const sxyz<T1> &Rx,
                        std::vector<sxyz<T1>> &r_data,
                        const size_t threadNo=0) const;
    
    T1 computeSlowness(const sxyz<T1>& Rx ) const;

    void sweep(const std::vector<bool>& frozen,
               const size_t threadNo) const;
    void sweep_weno3(const std::vector<bool>& frozen,
                     const size_t threadNo) const;
    
    void update_node(const size_t, const size_t, const size_t, const size_t=0) const;
    void update_node_weno3(const size_t, const size_t, const size_t, const size_t=0) const;
    
    void initFSM(const std::vector<sxyz<T1>>& Tx,
                 const std::vector<T1>& t0,
                 std::vector<bool>& frozen,
                 const int npts,
                 const size_t threadNo) const;

private:
    Grid3Dri() {}
    Grid3Dri(const Grid3Dri<T1,T2,NODE>& g) {}
    Grid3Dri<T1,T2,NODE>& operator=(const Grid3Dri<T1,T2,NODE>& g) { return *this; }
    
};



template<typename T1, typename T2, typename NODE>
void Grid3Dri<T1,T2,NODE>::buildGridNeighbors() {
    
    //Index the neighbors nodes of each cell
    for ( T2 n=0; n<nodes.size(); ++n ) {
        for ( size_t n2=0; n2<nodes[n].getOwners().size(); ++n2) {
            T2 check = nodes[n].getOwners()[n2];
            neighbors[ check ].push_back(n);
        }
    }
}

template<typename T1, typename T2, typename NODE>
int Grid3Dri<T1,T2,NODE>::checkPts(const std::vector<sxyz<T1>>& pts) const {
    
    // Check if the points from a vector are in the grid
    for ( size_t n=0; n<pts.size(); ++n ) {
        if ( pts[n].x < xmin || pts[n].x > xmax ||
            pts[n].y < ymin || pts[n].y > ymax ||
            pts[n].z < zmin || pts[n].z > zmax ) {
            std::cerr << "Error: point no " << (n+1)
            << " outside the grid.\n";
            return 1;
        }
    }
    return 0;
}


template<typename T1, typename T2, typename NODE>
T1 Grid3Dri<T1,T2,NODE>::getTraveltime(const sxyz<T1> &pt, const size_t nt) const {
    
    static const size_t nnx = ncx+1;
    static const size_t nny = ncy+1;
    
    // trilinear interpolation if not on node
    
    T1 tt;
    T2 i, j, k;
    
    getIJK(pt, i, j, k);
    
    if ( fabs(pt.x - (xmin+i*dx))<small &&
        fabs(pt.y - (ymin+j*dy))<small &&
        fabs(pt.z - (zmin+k*dz))<small ) {
        // on node
        return nodes[(k*nny+j)*nnx+i].getTT(nt);
    } else if ( fabs(pt.x - (xmin+i*dx))<small &&
               fabs(pt.y - (ymin+j*dy))<small ) {
        // on edge
        T1 t1 = nodes[(    k*nny+j)*nnx+i].getTT(nt);
        T1 t2 = nodes[((k+1)*nny+j)*nnx+i].getTT(nt);
        
        T1 w1 = (zmin+(k+1)*dz - pt.z)/dz;
        T1 w2 = (pt.z - (zmin+k*dz))/dz;
        
        tt = t1*w1 + t2*w2;
        
    } else if ( fabs(pt.x - (xmin+i*dx))<small &&
               fabs(pt.z - (zmin+k*dz))<small ) {
        // on edge
        T1 t1 = nodes[(k*nny+j  )*nnx+i].getTT(nt);
        T1 t2 = nodes[(k*nny+j+1)*nnx+i].getTT(nt);
        
        T1 w1 = (ymin+(j+1)*dy - pt.y)/dy;
        T1 w2 = (pt.y - (ymin+j*dy))/dy;
        
        tt = t1*w1 + t2*w2;
        
    } else if ( fabs(pt.y - (ymin+j*dy))<small &&
               fabs(pt.z - (zmin+k*dz))<small ) {
        // on edge
        T1 t1 = nodes[(k*nny+j)*nnx+i  ].getTT(nt);
        T1 t2 = nodes[(k*nny+j)*nnx+i+1].getTT(nt);
        
        T1 w1 = (xmin+(i+1)*dx - pt.x)/dx;
        T1 w2 = (pt.x - (xmin+i*dx))/dx;
        
        tt = t1*w1 + t2*w2;
        
    } else if ( fabs(pt.x - (xmin+i*dx))<small ) {
        // on YZ face
        T1 t1 = nodes[(    k*nny+j  )*nnx+i].getTT(nt);
        T1 t2 = nodes[((k+1)*nny+j  )*nnx+i].getTT(nt);
        T1 t3 = nodes[(    k*nny+j+1)*nnx+i].getTT(nt);
        T1 t4 = nodes[((k+1)*nny+j+1)*nnx+i].getTT(nt);
        
        T1 w1 = (zmin+(k+1)*dz - pt.z)/dz;
        T1 w2 = (pt.z - (zmin+k*dz))/dz;
        
        t1 = t1*w1 + t2*w2;
        t2 = t3*w1 + t4*w2;
        
        w1 = (ymin+(j+1)*dy - pt.y)/dy;
        w2 = (pt.y - (ymin+j*dy))/dy;
        
        tt = t1*w1 + t2*w2;
        
    } else if ( fabs(pt.y - (ymin+j*dy))<small ) {
        // on XZ face
        T1 t1 = nodes[(    k*nny+j)*nnx+i  ].getTT(nt);
        T1 t2 = nodes[((k+1)*nny+j)*nnx+i  ].getTT(nt);
        T1 t3 = nodes[(    k*nny+j)*nnx+i+1].getTT(nt);
        T1 t4 = nodes[((k+1)*nny+j)*nnx+i+1].getTT(nt);
        
        T1 w1 = (zmin+(k+1)*dz - pt.z)/dz;
        T1 w2 = (pt.z - (zmin+k*dz))/dz;
        
        t1 = t1*w1 + t2*w2;
        t2 = t3*w1 + t4*w2;
        
        w1 = (xmin+(i+1)*dx - pt.x)/dx;
        w2 = (pt.x - (xmin+i*dx))/dx;
        
        tt = t1*w1 + t2*w2;
        
    } else if ( fabs(pt.z - (zmin+k*dz))<small ) {
        // on XY face
        T1 t1 = nodes[(k*nny+j  )*nnx+i  ].getTT(nt);
        T1 t2 = nodes[(k*nny+j+1)*nnx+i  ].getTT(nt);
        T1 t3 = nodes[(k*nny+j  )*nnx+i+1].getTT(nt);
        T1 t4 = nodes[(k*nny+j+1)*nnx+i+1].getTT(nt);
        
        T1 w1 = (ymin+(j+1)*dy - pt.y)/dy;
        T1 w2 = (pt.y - (ymin+j*dy))/dy;
        
        t1 = t1*w1 + t2*w2;
        t2 = t3*w1 + t4*w2;
        
        w1 = (xmin+(i+1)*dx - pt.x)/dx;
        w2 = (pt.x - (xmin+i*dx))/dx;
        
        tt = t1*w1 + t2*w2;
        
    } else {
        T1 t1 = nodes[(    k*nny+j  )*nnx+i  ].getTT(nt);
        T1 t2 = nodes[((k+1)*nny+j  )*nnx+i  ].getTT(nt);
        T1 t3 = nodes[(    k*nny+j+1)*nnx+i  ].getTT(nt);
        T1 t4 = nodes[((k+1)*nny+j+1)*nnx+i  ].getTT(nt);
        T1 t5 = nodes[(    k*nny+j  )*nnx+i+1].getTT(nt);
        T1 t6 = nodes[((k+1)*nny+j  )*nnx+i+1].getTT(nt);
        T1 t7 = nodes[(    k*nny+j+1)*nnx+i+1].getTT(nt);
        T1 t8 = nodes[((k+1)*nny+j+1)*nnx+i+1].getTT(nt);
        
        T1 w1 = (zmin+(k+1)*dz - pt.z)/dz;
        T1 w2 = (pt.z - (zmin+k*dz))/dz;
        
        t1 = t1*w1 + t2*w2;
        t2 = t3*w1 + t4*w2;
        t3 = t5*w1 + t6*w2;
        t4 = t7*w1 + t8*w2;
        
        w1 = (ymin+(j+1)*dy - pt.y)/dy;
        w2 = (pt.y - (ymin+j*dy))/dy;
        
        t1 = t1*w1 + t2*w2;
        t2 = t3*w1 + t4*w2;
        
        w1 = (xmin+(i+1)*dx - pt.x)/dx;
        w2 = (pt.x - (xmin+i*dx))/dx;
        
        tt = t1*w1 + t2*w2;
        
    }
    
    return tt;
}

//template<typename T1, typename T2, typename NODE>
//T1 Grid3Dri<T1,T2,NODE>::getTraveltime(const sxyz<T1>& Rx,
//                                    const std::vector<NODE>& nodes,
//                                    const size_t threadNo) const {
//    
//    // Calculate and return the traveltime for a Rx point.
//    
//    // If Rx is on a node:
//    for ( size_t nn=0; nn<nodes.size(); ++nn ) {
//        if ( nodes[nn] == Rx ) {
//            return nodes[nn].getTT(threadNo);
//        }
//    }
//    //If Rx is not on a node:
//    T1 slo = computeSlowness( Rx );
//    
//    T2 cellNo = getCellNo( Rx );
//    T2 neibNo = neighbors[cellNo][0];
//    T1 dt = computeDt(nodes[neibNo], Rx, slo);
//    
//    T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
//    for ( size_t k=1; k< neighbors[cellNo].size(); ++k ) {
//        neibNo = neighbors[cellNo][k];
//        dt = computeDt(nodes[neibNo], Rx, slo);
//        if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
//            traveltime =  nodes[neibNo].getTT(threadNo)+dt;
//        }
//    }
//    return traveltime;
//}

template<typename T1, typename T2, typename NODE>
T1 Grid3Dri<T1,T2,NODE>::getTraveltime(const sxyz<T1>& Rx,
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
    
    T2 cellNo = getCellNo( Rx );
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

//template<typename T1, typename T2, typename NODE>
//T1 Grid3Dri<T1,T2,NODE>::getTraveltime(const sxyz<T1>& Rx,
//                                       const std::vector<NODE>& nodes,
//                                       T2& nodeParentRx, T2& cellParentRx,
//                                       const size_t threadNo) const {
//    
//    // Calculate and return the traveltime for a Rx point.
//    for ( size_t nn=0; nn<nodes.size(); ++nn ) {
//        if ( nodes[nn] == Rx ) {
//            nodeParentRx = nodes[nn].getNodeParent(threadNo);
//            cellParentRx = nodes[nn].getCellParent(threadNo);
//            return nodes[nn].getTT(threadNo);
//        }
//    }
//    //If Rx is not on a node:
//    T1 slo = computeSlowness( Rx );
//    
//    T2 cellNo = getCellNo( Rx );
//    T2 neibNo = neighbors[cellNo][0];
//    T1 dt = computeDt(nodes[neibNo], Rx, slo);
//    
//    T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
//    nodeParentRx = neibNo;
//    cellParentRx = cellNo;
//    for ( size_t k=1; k< neighbors[cellNo].size(); ++k ) {
//        neibNo = neighbors[cellNo][k];
//        dt = computeDt(nodes[neibNo], Rx, slo);
//        if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
//            traveltime =  nodes[neibNo].getTT(threadNo)+dt;
//            nodeParentRx = neibNo;
//        }
//    }
//    return traveltime;
//}


template<typename T1, typename T2, typename NODE>
void Grid3Dri<T1,T2,NODE>::grad(sxyz<T1>& g, const size_t i, const size_t j, const size_t k,
                                const size_t nt) const {
    
    // compute average gradient for voxel (i,j,k)

    static const size_t nnx = ncx+1;
    static const size_t nny = ncy+1;

    g.x = 0.25*(nodes[(    k*nny+j  )*nnx+i+1].getTT(nt) - nodes[(    k*nny+j  )*nnx+i  ].getTT(nt) +
                nodes[(    k*nny+j+1)*nnx+i+1].getTT(nt) - nodes[(    k*nny+j+1)*nnx+i  ].getTT(nt) +
                nodes[((k+1)*nny+j  )*nnx+i+1].getTT(nt) - nodes[((k+1)*nny+j  )*nnx+i  ].getTT(nt) +
                nodes[((k+1)*nny+j+1)*nnx+i+1].getTT(nt) - nodes[((k+1)*nny+j+1)*nnx+i  ].getTT(nt))/dx;
    g.y = 0.25*(nodes[(    k*nny+j+1)*nnx+i  ].getTT(nt) - nodes[(    k*nny+j  )*nnx+i  ].getTT(nt) +
                nodes[(    k*nny+j+1)*nnx+i+1].getTT(nt) - nodes[(    k*nny+j  )*nnx+i+1].getTT(nt) +
                nodes[((k+1)*nny+j+1)*nnx+i  ].getTT(nt) - nodes[((k+1)*nny+j  )*nnx+i  ].getTT(nt) +
                nodes[((k+1)*nny+j+1)*nnx+i+1].getTT(nt) - nodes[((k+1)*nny+j  )*nnx+i+1].getTT(nt))/dy;
    g.z = 0.25*(nodes[((k+1)*nny+j  )*nnx+i  ].getTT(nt) - nodes[(    k*nny+j  )*nnx+i  ].getTT(nt) +
                nodes[((k+1)*nny+j  )*nnx+i+1].getTT(nt) - nodes[(    k*nny+j  )*nnx+i+1].getTT(nt) +
                nodes[((k+1)*nny+j+1)*nnx+i  ].getTT(nt) - nodes[(    k*nny+j+1)*nnx+i  ].getTT(nt) +
                nodes[((k+1)*nny+j+1)*nnx+i+1].getTT(nt) - nodes[(    k*nny+j+1)*nnx+i+1].getTT(nt))/dz;
}


template<typename T1, typename T2, typename NODE>
void Grid3Dri<T1,T2,NODE>::grad(sxyz<T1>& g, const sxyz<T1> &pt,
                                const size_t nt) const {
    
    // compute travel time gradient at point pt
    
    T1 p1 = pt.x - dx/2.0;
    T1 p2 = p1 + dx;
    g.x = (getTraveltime({p2, pt.y, pt.z}, nt) - getTraveltime({p1, pt.y, pt.z}, nt)) / dx;

    p1 = pt.y - dy/2.0;
    p2 = p1 + dy;
    g.y = (getTraveltime({pt.x, p2, pt.z}, nt) - getTraveltime({pt.x, p1, pt.z}, nt)) / dy;

    p1 = pt.z - dz/2.0;
    p2 = p1 + dz;
    g.z = (getTraveltime({pt.x, pt.y, p2}, nt) - getTraveltime({pt.x, pt.y, p1}, nt)) / dz;

}




template<typename T1, typename T2, typename NODE>
void Grid3Dri<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                      const sxyz<T1> &Rx,
                                      std::vector<sxyz<T1>> &r_data,
                                      const size_t threadNo) const {

    r_data.push_back( Rx );
    
    for ( size_t ns=0; ns<Tx.size(); ++ns ) {
        if ( Rx == Tx[ns] ) {
            return;
        }
    }
    
    sxyz<T1> curr_pt( Rx );
    // distance between opposite nodes of a voxel
    static const T1 maxDist = sqrt( dx*dx + dy*dy + dz*dz );
    sxyz<T1> g;
    
    bool reachedTx = false;
    while ( reachedTx == false ) {
        
        grad(g, curr_pt, threadNo);
        g *= -1.0;
        
        long long i, j, k;
        getIJK(curr_pt, i, j, k);
        
        // planes we will intersect
        T1 xp = xmin + dx*(i + boost::math::sign(g.x)>0.0 ? 1.0 : 0.0);
        T1 yp = ymin + dy*(j + boost::math::sign(g.y)>0.0 ? 1.0 : 0.0);
        T1 zp = zmin + dz*(k + boost::math::sign(g.z)>0.0 ? 1.0 : 0.0);
        
        if ( fabs(xp-curr_pt.x)<small) {
            xp += dx*boost::math::sign(g.x);
        }
        if ( fabs(yp-curr_pt.y)<small) {
            yp += dy*boost::math::sign(g.y);
        }
        if ( fabs(zp-curr_pt.z)<small) {
            zp += dz*boost::math::sign(g.z);
        }
        
        // dist to planes
        T1 tx = (xp - curr_pt.x)/g.x;
        T1 ty = (yp - curr_pt.y)/g.y;
        T1 tz = (zp - curr_pt.z)/g.z;
        
        if ( tx<ty && tx<tz ) { // closer to xp
            curr_pt += tx*g;
            curr_pt.x = xp;     // make sure we don't accumulate rounding errors
        } else if ( ty<tz ) {
            curr_pt += ty*g;
            curr_pt.y = yp;
        } else {
            curr_pt += tz*g;
            curr_pt.z = zp;
        }
        
        if ( curr_pt.x < xmin || curr_pt.x > xmax ||
            curr_pt.y < ymin || curr_pt.y > ymax ||
            curr_pt.z < zmin || curr_pt.z > zmax ) {
            //  we are going oustide the grid!
            std::cerr << "Error while computing raypaths: going outside grid!\n"
            << "  Stopping calculations, raypaths will be incomplete.\n" << std::endl;
            return;

        }

        r_data.push_back( curr_pt );
        
        // are we close enough to one the Tx nodes ?
        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( curr_pt.getDistance( Tx[ns] ) < maxDist ) {
                r_data.push_back( Tx[ns] );
                reachedTx = true;
            }
        }
    }
}


template<typename T1, typename T2, typename NODE>
void Grid3Dri<T1,T2,NODE>::getRaypath_old(const std::vector<sxyz<T1>>& Tx,
                                          const sxyz<T1> &Rx,
                                          std::vector<sxyz<T1>> &r_data,
                                          const size_t threadNo) const {
    
    r_data.push_back( Rx );
    
    for ( size_t ns=0; ns<Tx.size(); ++ns ) {
        if ( Rx == Tx[ns] ) {
            return;
        }
    }
    
    long long int iIn, jIn, kIn, iOut=-1, jOut=-1, kOut=-1; // Out for cell we are exiting; In for cell we are entering
    sxyz<T1> curr_pt( Rx );
    sxyz<T1> gOut = {0.0, 0.0, 0.0};
    
    // distance between opposite nodes of a voxel
    static const T1 maxDist = sqrt( dx*dx + dy*dy + dz*dz );

    getIJK(curr_pt, iIn, jIn, kIn);
    
    bool reachedTx = false;
    while ( reachedTx == false ) {
        
        bool onNode=false;
        bool onEdgeX=false;
        bool onEdgeY=false;
        bool onEdgeZ=false;
        
        if ( fabs(remainder(curr_pt.x,dx))<small &&
            fabs(remainder(curr_pt.y,dy))<small &&
            fabs(remainder(curr_pt.z,dz))<small ) {
            onNode = true;
        } else if ( fabs(remainder(curr_pt.y,dy))<small &&
                   fabs(remainder(curr_pt.z,dz))<small ) {
            onEdgeX = true;
        } else if ( fabs(remainder(curr_pt.x,dx))<small &&
                   fabs(remainder(curr_pt.z,dz))<small ) {
            onEdgeY = true;
        } else if ( fabs(remainder(curr_pt.x,dx))<small &&
                   fabs(remainder(curr_pt.y,dy))<small ) {
            onEdgeZ = true;
        }
        
        if ( onNode ) {
            
            T2 i, j, k;
            getIJK(curr_pt, i, j, k);
            std::vector<sijk<T2>> voxels;

            // find voxels touching node
            if ( i<=ncx && j<=ncy && k<=ncz )
                voxels.push_back( {i,j,k} );
            if ( i<=ncx && j<=ncy && k>0 )
                voxels.push_back( {i,j,k-1} );
            if ( i<=ncx && j>0 && k<=ncz )
                voxels.push_back( {i,j-1,k} );
            if ( i<=ncx && j>0 && k>0 )
                voxels.push_back( {i,j-1,k-1} );
            if ( i>0 && j<=ncy && k<=ncz )
                voxels.push_back( {i-1,j,k} );
            if ( i>0 && j<=ncy && k>0 )
                voxels.push_back( {i-1,j,k-1} );
            if ( i>0 && j>0 && k<=ncz )
                voxels.push_back( {i-1,j-1,k} );
            if ( i>0 && j>0 && k>0 )
                voxels.push_back( {i-1,j-1,k-1} );
            
            gOut = static_cast<T1>(0.0);
            sxyz<T1> g;
            T2 nc;
            for ( nc=0; nc<voxels.size(); ++nc ) {
                grad( g, voxels[nc].i, voxels[nc].j, voxels[nc].k, threadNo );
                g *= -1.0;
                gOut += g;
            }
            gOut /= nc;  // gOut holds average grad
            
            if ((gOut.x<0.0 && i==0) || (gOut.x>0.0 && i==ncx+1) ||
                (gOut.y<0.0 && j==0) || (gOut.y>0.0 && j==ncy+1) ||
                (gOut.z<0.0 && k==0) || (gOut.z>0.0 && k==ncz+1)) {
                //  we are going oustide the grid!
                std::cerr << "Error while computing raypaths: going outside grid!\n"
                << "  Stopping calculations, raypaths will be incomplete.\n" << std::endl;
                return;
            }
            
            iOut = boost::math::sign(gOut.x)<0.0 ? i-1 : i;
            jOut = boost::math::sign(gOut.y)<0.0 ? j-1 : j;
            kOut = boost::math::sign(gOut.z)<0.0 ? k-1 : k;

            // planes we will intersect
            T1 xp = xmin + dx*(i + boost::math::sign(gOut.x));
            T1 yp = ymin + dy*(j + boost::math::sign(gOut.y));
            T1 zp = zmin + dz*(k + boost::math::sign(gOut.z));
            
            // dist to planes
            T1 tx = (xp - curr_pt.x)/gOut.x;
            T1 ty = (yp - curr_pt.y)/gOut.y;
            T1 tz = (zp - curr_pt.z)/gOut.z;
            
            if ( tx<ty && tx<tz ) { // closer to xp
                curr_pt += tx*gOut;
                curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                iIn = iOut + boost::math::sign(gOut.x);
                jIn = jOut;
                kIn = kOut;
            } else if ( ty<tz ) {
                curr_pt += ty*gOut;
                curr_pt.y = yp;
                iIn = iOut;
                jIn = jOut + boost::math::sign(gOut.y);
                kIn = kOut;
            } else {
                curr_pt += tz*gOut;
                curr_pt.z = zp;
                iIn = iOut;
                jIn = jOut;
                kIn = kOut + boost::math::sign(gOut.z);
            }

        } else if ( onEdgeX ) {
            
            T2 i, j, k;
            getIJK(curr_pt, i, j, k);
            std::vector<sijk<T2>> voxels;
            
            // find voxels touching edge
            if ( i<=ncx && j<=ncy && k<=ncz )
                voxels.push_back( {i,j,k} );
            if ( i<=ncx && j<=ncy && k>0 )
                voxels.push_back( {i,j,k-1} );
            if ( i<=ncx && j>0 && k<=ncz )
                voxels.push_back( {i,j-1,k} );
            if ( i<=ncx && j>0 && k>0 )
                voxels.push_back( {i,j-1,k-1} );
            
            gOut = static_cast<T1>(0.0);
            sxyz<T1> g;
            T2 nc;
            for ( nc=0; nc<voxels.size(); ++nc ) {
                grad( g, voxels[nc].i, voxels[nc].j, voxels[nc].k, threadNo );
                g *= -1.0;
                gOut += g;
            }
            gOut /= nc;  // gOut holds average grad

            if ((gOut.x<0.0 && i==0) || (gOut.x>0.0 && i==ncx+1) ||
                (gOut.y<0.0 && j==0) || (gOut.y>0.0 && j==ncy+1) ||
                (gOut.z<0.0 && k==0) || (gOut.z>0.0 && k==ncz+1)) {
                //  we are going oustide the grid!
                std::cerr << "Error while computing raypaths: going outside grid!\n"
                << "  Stopping calculations, raypaths will be incomplete.\n" << std::endl;
                return;
            }
            
            iOut = i;
            jOut = boost::math::sign(gOut.y)<0.0 ? j-1 : j;
            kOut = boost::math::sign(gOut.z)<0.0 ? k-1 : k;
            
            // planes we will intersect
            T1 xp = xmin + dx*(i + boost::math::sign(gOut.x)>0.0 ? 1.0 : 0.0);
            T1 yp = ymin + dy*(j + boost::math::sign(gOut.y));
            T1 zp = zmin + dz*(k + boost::math::sign(gOut.z));
            
            if ( fabs(xp-curr_pt.x)<small) {
                xp += dx*boost::math::sign(gOut.x);
                iOut += boost::math::sign(gOut.x);
            }
            
            // dist to planes
            T1 tx = (xp - curr_pt.x)/gOut.x;
            T1 ty = (yp - curr_pt.y)/gOut.y;
            T1 tz = (zp - curr_pt.z)/gOut.z;
            
            if ( tx<ty && tx<tz ) { // closer to xp
                curr_pt += tx*gOut;
                curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                iIn = iOut + boost::math::sign(gOut.x);
                jIn = jOut;
                kIn = kOut;
            } else if ( ty<tz ) {
                curr_pt += ty*gOut;
                curr_pt.y = yp;
                iIn = iOut;
                jIn = jOut + boost::math::sign(gOut.y);
                kIn = kOut;
            } else {
                curr_pt += tz*gOut;
                curr_pt.z = zp;
                iIn = iOut;
                jIn = jOut;
                kIn = kOut + boost::math::sign(gOut.z);
            }

        } else if ( onEdgeY ) {
            
            T2 i, j, k;
            getIJK(curr_pt, i, j, k);
            std::vector<sijk<T2>> voxels;
            
            // find voxels touching node
            if ( i<=ncx && j<=ncy && k<=ncz )
                voxels.push_back( {i,j,k} );
            if ( i<=ncx && j<=ncy && k>0 )
                voxels.push_back( {i,j,k-1} );
            if ( i>0 && j<=ncy && k<=ncz )
                voxels.push_back( {i-1,j,k} );
            if ( i>0 && j<=ncy && k>0 )
                voxels.push_back( {i-1,j,k-1} );
            
            gOut = static_cast<T1>(0.0);
            sxyz<T1> g;
            T2 nc;
            for ( nc=0; nc<voxels.size(); ++nc ) {
                grad( g, voxels[nc].i, voxels[nc].j, voxels[nc].k, threadNo );
                g *= -1.0;
                gOut += g;
            }
            gOut /= nc;  // gOut holds average grad
            
            if ((gOut.x<0.0 && i==0) || (gOut.x>0.0 && i==ncx+1) ||
                (gOut.y<0.0 && j==0) || (gOut.y>0.0 && j==ncy+1) ||
                (gOut.z<0.0 && k==0) || (gOut.z>0.0 && k==ncz+1)) {
                //  we are going oustide the grid!
                std::cerr << "Error while computing raypaths: going outside grid!\n"
                << "  Stopping calculations, raypaths will be incomplete.\n" << std::endl;
                return;
            }

            iOut = boost::math::sign(gOut.x)<0.0 ? i-1 : i;
            jOut = j;
            kOut = boost::math::sign(gOut.z)<0.0 ? k-1 : k;
            
            // planes we will intersect
            T1 xp = xmin + dx*(i + boost::math::sign(gOut.x));
            T1 yp = ymin + dy*(j + boost::math::sign(gOut.y)>0.0 ? 1.0 : 0.0);
            T1 zp = zmin + dz*(k + boost::math::sign(gOut.z));
            
            if ( fabs(yp-curr_pt.y)<small) {
                yp += dy*boost::math::sign(gOut.y);
                jOut += boost::math::sign(gOut.y);
            }
            
            // dist to planes
            T1 tx = (xp - curr_pt.x)/gOut.x;
            T1 ty = (yp - curr_pt.y)/gOut.y;
            T1 tz = (zp - curr_pt.z)/gOut.z;
            
            if ( tx<ty && tx<tz ) { // closer to xp
                curr_pt += tx*gOut;
                curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                iIn = iOut + boost::math::sign(gOut.x);
                jIn = jOut;
                kIn = kOut;
            } else if ( ty<tz ) {
                curr_pt += ty*gOut;
                curr_pt.y = yp;
                iIn = iOut;
                jIn = jOut + boost::math::sign(gOut.y);
                kIn = kOut;
            } else {
                curr_pt += tz*gOut;
                curr_pt.z = zp;
                iIn = iOut;
                jIn = jOut;
                kIn = kOut + boost::math::sign(gOut.z);
            }

        } else if ( onEdgeZ ) {
            
            T2 i, j, k;
            getIJK(curr_pt, i, j, k);
            std::vector<sijk<T2>> voxels;
            
            // find voxels touching node
            if ( i<=ncx && j<=ncy && k<=ncz )
                voxels.push_back( {i,j,k} );
            if ( i<=ncx && j>0 && k<=ncz )
                voxels.push_back( {i,j-1,k} );
            if ( i>0 && j<=ncy && k<=ncz )
                voxels.push_back( {i-1,j,k} );
            if ( i>0 && j>0 && k<=ncz )
                voxels.push_back( {i-1,j-1,k} );
            
            gOut = static_cast<T1>(0.0);
            sxyz<T1> g;
            T2 nc;
            for ( nc=0; nc<voxels.size(); ++nc ) {
                grad( g, voxels[nc].i, voxels[nc].j, voxels[nc].k, threadNo );
                g *= -1.0;
                gOut += g;
            }
            gOut /= nc;  // gOut holds average grad
            
            if ((gOut.x<0.0 && i==0) || (gOut.x>0.0 && i==ncx+1) ||
                (gOut.y<0.0 && j==0) || (gOut.y>0.0 && j==ncy+1) ||
                (gOut.z<0.0 && k==0) || (gOut.z>0.0 && k==ncz+1)) {
                //  we are going oustide the grid!
                std::cerr << "Error while computing raypaths: going outside grid!\n"
                << "  Stopping calculations, raypaths will be incomplete.\n" << std::endl;
                return;
            }
            
            iOut = boost::math::sign(gOut.x)<0.0 ? i-1 : i;
            jOut = boost::math::sign(gOut.y)<0.0 ? j-1 : j;
            kOut = k;
            
            // planes we will intersect
            T1 xp = xmin + dx*(i + boost::math::sign(gOut.x));
            T1 yp = ymin + dy*(j + boost::math::sign(gOut.y));
            T1 zp = zmin + dz*(k + boost::math::sign(gOut.z)>0.0 ? 1.0 : 0.0);
            
            if ( fabs(zp-curr_pt.z)<small) {
                zp += dz*boost::math::sign(gOut.z);
                kOut += boost::math::sign(gOut.z);
            }
            
            // dist to planes
            T1 tx = (xp - curr_pt.x)/gOut.x;
            T1 ty = (yp - curr_pt.y)/gOut.y;
            T1 tz = (zp - curr_pt.z)/gOut.z;
            
            if ( tx<ty && tx<tz ) { // closer to xp
                curr_pt += tx*gOut;
                curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                iIn = iOut + boost::math::sign(gOut.x);
                jIn = jOut;
                kIn = kOut;
            } else if ( ty<tz ) {
                curr_pt += ty*gOut;
                curr_pt.y = yp;
                iIn = iOut;
                jIn = jOut + boost::math::sign(gOut.y);
                kIn = kOut;
            } else {
                curr_pt += tz*gOut;
                curr_pt.z = zp;
                iIn = iOut;
                jIn = jOut;
                kIn = kOut + boost::math::sign(gOut.z);
            }

        } else {
            
            sxyz<T1> gIn;
            grad( gIn, iIn, jIn, kIn, threadNo );
            gIn *= -1.0;
            
            if ( iIn == iOut && jIn == jOut && kIn == kOut) {
                // ray is returning to cell it was exiting
                // we might be at grazing incidence
                // check if gIn is significantly different from gOut
                
                sxyz<T1> diff = normalize(gOut)-normalize(gIn);
                if ( norm(diff) > small ) {
                    std::cerr << "Error while computing raypaths: raypath not converging!\n"
                    << "  Stopping calculations, raypaths will be incomplete.\n" << std::endl;
                    return;
                }
            }
            
            gOut = gIn;
            iOut = iIn;
            jOut = jIn;
            kOut = kIn;
            
            if ((gOut.x<0.0 && iOut==0) || (gOut.x>0.0 && iOut==ncx+1) ||
                (gOut.y<0.0 && jOut==0) || (gOut.y>0.0 && jOut==ncy+1) ||
                (gOut.z<0.0 && kOut==0) || (gOut.z>0.0 && kOut==ncz+1)) {
                //  we are going oustide the grid!
                std::cerr << "Error while computing raypaths: going outside grid!\n"
                << "  Stopping calculations, raypaths will be incomplete.\n" << std::endl;
                return;
            }
            
            // planes we will intersect
            T1 xp = xmin + dx*(iIn + boost::math::sign(gOut.x)>0.0 ? 1.0 : 0.0);
            T1 yp = ymin + dy*(jIn + boost::math::sign(gOut.y)>0.0 ? 1.0 : 0.0);
            T1 zp = zmin + dz*(kIn + boost::math::sign(gOut.z)>0.0 ? 1.0 : 0.0);
            
            if ( fabs(xp-curr_pt.x)<small) {
                xp += dx*boost::math::sign(gOut.x);
                iOut += boost::math::sign(gOut.x);
            }
            if ( fabs(yp-curr_pt.y)<small) {
                yp += dy*boost::math::sign(gOut.y);
                jOut += boost::math::sign(gOut.y);
            }
            if ( fabs(zp-curr_pt.z)<small) {
                zp += dz*boost::math::sign(gOut.z);
                kOut += boost::math::sign(gOut.z);
            }
            
            // dist to planes
            T1 tx = (xp - curr_pt.x)/gOut.x;
            T1 ty = (yp - curr_pt.y)/gOut.y;
            T1 tz = (zp - curr_pt.z)/gOut.z;
            
            if ( tx<ty && tx<tz ) { // closer to xp
                curr_pt += tx*gOut;
                curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                iIn = iOut + boost::math::sign(gOut.x);
                jIn = jOut;
                kIn = kOut;
            } else if ( ty<tz ) {
                curr_pt += ty*gOut;
                curr_pt.y = yp;
                iIn = iOut;
                jIn = jOut + boost::math::sign(gOut.y);
                kIn = kOut;
            } else {
                curr_pt += tz*gOut;
                curr_pt.z = zp;
                iIn = iOut;
                jIn = jOut;
                kIn = kOut + boost::math::sign(gOut.z);
            }

        }
        
        
        r_data.push_back( curr_pt );

        // are we close enough to one the of Tx nodes ?
        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( curr_pt.getDistance( Tx[ns] ) < maxDist ) {
                r_data.push_back( Tx[ns] );
                reachedTx = true;
            }
        }
    }
}

template<typename T1, typename T2, typename NODE>
T1 Grid3Dri<T1,T2,NODE>::computeSlowness(const sxyz<T1>& Rx) const {
    
    
    // Calculate the slowness of any point that is not on a node
    
    T2 cellNo = getCellNo( Rx );
    
    //We calculate the Slowness at the point
    std::vector<T2> list;
    
    for (size_t n3=0; n3 < neighbors[ cellNo ].size(); n3++){
        if ( nodes[neighbors[ cellNo ][n3] ].getPrimary() == 5 ){
            list.push_back(neighbors[ cellNo ][n3]);
        }
    }
    
    if ( inverseDistance ) {
        
        std::vector<size_t>::iterator it;
        
        std::vector<NODE*> interpNodes;
        
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

template<typename T1, typename T2, typename NODE>
void Grid3Dri<T1,T2,NODE>::saveTT(const std::string & fname, const int all,
                                  const size_t nt,
                                  const bool vtkFormat) const {
    
    if (vtkFormat) {
#ifdef VTK
        
        std::string filename = fname+".vtr";
        int nn[3] = {static_cast<int>(ncx+1), static_cast<int>(ncy+1), static_cast<int>(ncz+1)};
        
        vtkSmartPointer<vtkDoubleArray> xCoords = vtkSmartPointer<vtkDoubleArray>::New();
        for (size_t n=0; n<nn[0]; ++n)
            xCoords->InsertNextValue( xmin + n*dx );
        vtkSmartPointer<vtkDoubleArray> yCoords = vtkSmartPointer<vtkDoubleArray>::New();
        for (size_t n=0; n<nn[1]; ++n)
            yCoords->InsertNextValue( ymin + n*dy );
        vtkSmartPointer<vtkDoubleArray> zCoords = vtkSmartPointer<vtkDoubleArray>::New();
        for (size_t n=0; n<nn[2]; ++n)
            zCoords->InsertNextValue( zmin + n*dz );
        
        vtkSmartPointer<vtkRectilinearGrid> rgrid = vtkSmartPointer<vtkRectilinearGrid>::New();
        rgrid->SetDimensions( nn );
        rgrid->SetXCoordinates(xCoords);
        rgrid->SetYCoordinates(yCoords);
        rgrid->SetZCoordinates(zCoords);

        vtkSmartPointer<vtkDoubleArray> newScalars =
        vtkSmartPointer<vtkDoubleArray>::New();
        
        newScalars->SetName("Travel time");
        newScalars->SetNumberOfComponents(1);
        newScalars->SetNumberOfTuples( rgrid->GetNumberOfPoints() );
        
        for ( size_t n=0; n<nodes.size(); ++n ) {
            if ( nodes[n].isPrimary() == true ) {
                vtkIdType id = rgrid->FindPoint(nodes[n].getX(), nodes[n].getY(), nodes[n].getZ());
                newScalars->SetTuple1(id, nodes[n].getTT(nt) );
            }
        }
        rgrid->GetPointData()->SetScalars(newScalars);
        
        vtkSmartPointer<vtkXMLRectilinearGridWriter> writer =
        vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
        
        writer->SetFileName( filename.c_str() );
        writer->SetInputData( rgrid );
        writer->SetDataModeToBinary();
        writer->Update();
#else
        std::cerr << "VTK not included during compilation.\nNothing saved.\n";
#endif
    } else {
        std::ofstream fout(fname.c_str());
        fout.precision(12);
        for ( T2 n=0; n<nodes.size(); ++n ) {
            if ( nodes[n].isPrimary() == true || all==1 ) {
                fout << nodes[n].getX() << '\t'
                << nodes[n].getY() << '\t'
                << nodes[n].getZ() << '\t'
                << nodes[n].getTT(nt) << '\n';
            }
        }
        fout.close();
    }
}



template<typename T1, typename T2, typename NODE>
void Grid3Dri<T1,T2,NODE>::sweep(const std::vector<bool>& frozen,
                                 const size_t threadNo) const {
    
    // sweep first direction
    for ( size_t k=0; k<=ncz; ++k ) {
        for ( size_t j=0; j<=ncy; ++j ) {
            for ( size_t i=0; i<=ncx; ++i ) {
                if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                    update_node(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep second direction
    for ( size_t k=0; k<=ncz; ++k ) {
        for ( size_t j=0; j<=ncy; ++j ) {
            for ( long int i=ncx; i>=0; --i ) {
                if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                    update_node(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep third direction
    for ( size_t k=0; k<=ncz; ++k ) {
        for ( long int j=ncy; j>=0; --j ) {
            for ( size_t i=0; i<=ncx; ++i ) {
                if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                    update_node(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep fourth direction
    for ( size_t k=0; k<=ncz; ++k ) {
        for ( long int j=ncy; j>=0; --j ) {
            for ( long int i=ncx; i>=0; --i ) {
                if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                    update_node(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep fifth direction
    for ( long int k=ncz; k>=0; --k ) {
        for ( size_t j=0; j<=ncy; ++j ) {
            for ( size_t i=0; i<=ncx; ++i ) {
                if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                    update_node(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep sixth direction
    for ( long int k=ncz; k>=0; --k ) {
        for ( size_t j=0; j<=ncy; ++j ) {
            for ( long int i=ncx; i>=0; --i ) {
                if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                    update_node(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep seventh direction
    for ( long int k=ncz; k>=0; --k ) {
        for ( long int j=ncy; j>=0; --j ) {
            for ( size_t i=0; i<=ncx; ++i ) {
                if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                    update_node(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep eighth direction
    for ( long int k=ncz; k>=0; --k ) {
        for ( long int j=ncy; j>=0; --j ) {
            for ( long int i=ncx; i>=0; --i ) {
                if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                    update_node(i, j, k, threadNo);
                }
            }
        }
    }
}

template<typename T1, typename T2, typename NODE>
void Grid3Dri<T1,T2,NODE>::update_node(const size_t i, const size_t j, const size_t k,
                                   const size_t threadNo) const {
    T1 a1, a2, a3, t;
    
    if (k==0)
        a1 = nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
    else if (k==ncz)
        a1 = nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
    else {
        a1 = nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
        t  = nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
        a1 = a1<t ? a1 : t;
    }
    
    if (j==0)
        a2 = nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo);
    else if (j==ncy)
        a2 = nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo);
    else {
        a2 = nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo);
        t  = nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo);
        a2 = a2<t ? a2 : t;
    }
    
    if (i==0)
        a3 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo);
    else if (i==ncx)
        a3 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo);
    else {
        a3 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo);
        t  = nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo);
        a3 = a3<t ? a3 : t;
    }
    
    if ( a1>a2 ) std::swap(a1, a2);
    if ( a1>a3 ) std::swap(a1, a3);
    if ( a2>a3 ) std::swap(a2, a3);
    
    T1 fh = nodes[(k*(ncy+1)+j)*(ncx+1)+i].getNodeSlowness() * dx;
    
    t = a1 + fh;
    if ( t > a2 ) {
        
        t = 0.5*(a1+a2+sqrt(2.*fh*fh - (a1-a2)*(a1-a2)));
        
        if ( t > a3 ) {
            
            t = 1./3. * ((a1 + a2 + a3) + sqrt(-2.*a1*a1 + 2.*a1*a2 - 2.*a2*a2 +
                                               2.*a1*a3 + 2.*a2*a3 -
                                               2.*a3*a3 + 3.*fh*fh));
            
        }
    }
    
    if ( t<nodes[(k*(ncy+1)+j)*(ncx+1)+i].getTT(threadNo) )
        nodes[(k*(ncy+1)+j)*(ncx+1)+i].setTT(t,threadNo);
    
}

template<typename T1, typename T2, typename NODE>
void Grid3Dri<T1,T2,NODE>::sweep_weno3(const std::vector<bool>& frozen,
                                       const size_t threadNo) const {
    
    // sweep first direction
    for ( size_t k=0; k<=ncz; ++k ) {
        for ( size_t j=0; j<=ncy; ++j ) {
            for ( size_t i=0; i<=ncx; ++i ) {
                if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                    update_node_weno3(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep second direction
    for ( size_t k=0; k<=ncz; ++k ) {
        for ( size_t j=0; j<=ncy; ++j ) {
            for ( long int i=ncx; i>=0; --i ) {
                if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                    update_node_weno3(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep third direction
    for ( size_t k=0; k<=ncz; ++k ) {
        for ( long int j=ncy; j>=0; --j ) {
            for ( size_t i=0; i<=ncx; ++i ) {
                if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                    update_node_weno3(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep fourth direction
    for ( size_t k=0; k<=ncz; ++k ) {
        for ( long int j=ncy; j>=0; --j ) {
            for ( long int i=ncx; i>=0; --i ) {
                if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                    update_node_weno3(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep fifth direction
    for ( long int k=ncz; k>=0; --k ) {
        for ( size_t j=0; j<=ncy; ++j ) {
            for ( size_t i=0; i<=ncx; ++i ) {
                if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                    update_node_weno3(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep sixth direction
    for ( long int k=ncz; k>=0; --k ) {
        for ( size_t j=0; j<=ncy; ++j ) {
            for ( long int i=ncx; i>=0; --i ) {
                if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                    update_node_weno3(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep seventh direction
    for ( long int k=ncz; k>=0; --k ) {
        for ( long int j=ncy; j>=0; --j ) {
            for ( size_t i=0; i<=ncx; ++i ) {
                if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                    update_node_weno3(i, j, k, threadNo);
                }
            }
        }
    }
    // sweep eighth direction
    for ( long int k=ncz; k>=0; --k ) {
        for ( long int j=ncy; j>=0; --j ) {
            for ( long int i=ncx; i>=0; --i ) {
                if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                    update_node_weno3(i, j, k, threadNo);
                }
            }
        }
    }
}

template<typename T1, typename T2, typename NODE>
void Grid3Dri<T1,T2,NODE>::update_node_weno3(const size_t i,
                                             const size_t j,
                                             const size_t k,
                                             const size_t threadNo) const {
    T1 a1, a2, a3, t;
    
    if (k==0) {
        a1 = nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);  // first order
    } else if (k==1) {
        T1 num = nodes[ ((k+2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
        2.*nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
        nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
        num *= num;
        T1 den = nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
        2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
        nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
        den *= den;
        T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
        T1 w = 1./(1.+2.*r*r);
        
        T1 ap = (1.-w)*(nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo)-
                        nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx) +
        w*(-nodes[ ((k+2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
           4.*nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
           3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx);
        
        a1 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) + dx*ap;
        
        t = nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo); // first order for left
        a1 = a1<t ? a1 : t;
        
    } else if (k==ncz) {
        a1 = nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
    } else if (k==ncz-1) {
        T1 num = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
        2.*nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
        nodes[ ((k-2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
        num *= num;
        T1 den = nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
        2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
        nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
        den *= den;
        T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
        T1 w = 1./(1.+2.*r*r);
        
        T1 am = (1.-w)*(nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo)-
                        nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx) +
        w*(3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
           4.*nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
           nodes[ ((k-2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx);
        
        a1 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) - dx*am;
        
        t = nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo); // first order for right
        a1 = a1<t ? a1 : t;
        
    } else {
        T1 num = nodes[ ((k+2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
        2.*nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
        nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
        num *= num;
        T1 den = nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
        2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
        nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
        den *= den;
        T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
        T1 w = 1./(1.+2.*r*r);
        
        T1 ap = (1.-w)*(nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo)-
                        nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx) +
        w*(-nodes[ ((k+2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
           4.*nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
           3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx);
        
        a1 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) + dx*ap;

        num = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
        2.*nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
        nodes[ ((k-2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
        num *= num;
        r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
        w = 1./(1.+2.*r*r);
        
        T1 am = (1.-w)*(nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo)-
                        nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx) +
        w*(3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
           4.*nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
           nodes[ ((k-2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx);
        
        t = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) - dx*am;
        a1 = a1<t ? a1 : t;
        
    }
    
    if (j==0) {
        a2 = nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo);
    } else if (j==1) {
        T1 num = nodes[ (k*(ncy+1)+j+2)*(ncx+1)+i ].getTT(threadNo) -
        2.*nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo) +
        nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
        num *= num;
        T1 den = nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo) -
        2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
        nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo);
        den *= den;
        T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
        T1 w = 1./(1.+2.*r*r);
        
        T1 ap = (1.-w)*(nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo)-
                        nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo))/(2.*dx) +
        w*(-nodes[ (k*(ncy+1)+j+2)*(ncx+1)+i ].getTT(threadNo) +
           4.*nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo) -
           3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx);
        
        a2 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) + dx*ap;
        
        t = nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo); // first order for left
        a2 = a2<t ? a2 : t;

    } else if (j==ncy) {
        a2 = nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo);
    } else if (j==ncy-1) {
        T1 num = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
        2.*nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo) +
        nodes[ (k*(ncy+1)+j-2)*(ncx+1)+i ].getTT(threadNo);
        num *= num;
        T1 den = nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo) -
        2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
        nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo);
        den *= den;
        T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
        T1 w = 1./(1.+2.*r*r);
        
        T1 am = (1.-w)*(nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo)-
                        nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo))/(2.*dx) +
        w*(3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
           4.*nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo) +
           nodes[ (k*(ncy+1)+j-2)*(ncx+1)+i ].getTT(threadNo))/(2.*dx);
        
        a2 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) - dx*am;
        
        t = nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo); // first order for right
        a2 = a2<t ? a2 : t;
        
    } else {
        T1 num = nodes[ (k*(ncy+1)+j+2)*(ncx+1)+i ].getTT(threadNo) -
        2.*nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo) +
        nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
        num *= num;
        T1 den = nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo) -
        2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
        nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo);
        den *= den;
        T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
        T1 w = 1./(1.+2.*r*r);
        
        T1 ap = (1.-w)*(nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo)-
                        nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo))/(2.*dx) +
        w*(-nodes[ (k*(ncy+1)+j+2)*(ncx+1)+i ].getTT(threadNo) +
           4.*nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo) -
           3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx);
        
        a2 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) + dx*ap;
        
        num = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
        2.*nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo) +
        nodes[ (k*(ncy+1)+j-2)*(ncx+1)+i ].getTT(threadNo);
        num *= num;
        r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
        w = 1./(1.+2.*r*r);
        
        T1 am = (1.-w)*(nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo)-
                        nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo))/(2.*dx) +
        w*(3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
           4.*nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo) +
           nodes[ (k*(ncy+1)+j-2)*(ncx+1)+i ].getTT(threadNo))/(2.*dx);
        
        t = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) - dx*am;
        a2 = a2<t ? a2 : t;
        
    }
    
    if (i==0) {
        a3 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo);
    } else if (i==1) {
        T1 num = nodes[ (k*(ncy+1)+j)*(ncx+1)+i+2 ].getTT(threadNo) -
        2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo) +
        nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
        num *= num;
        T1 den = nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo) -
        2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
        nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo);
        den *= den;
        T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
        T1 w = 1./(1.+2.*r*r);
        
        T1 ap = (1.-w)*(nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo)-
                        nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo))/(2.*dx) +
        w*(-nodes[ (k*(ncy+1)+j)*(ncx+1)+i+2 ].getTT(threadNo) +
           4.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo) -
           3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx);
        
        a3 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) + dx*ap;
        
        t = nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo); // first order for left
        a3 = a3<t ? a3 : t;
        
    } else if (i==ncx) {
        a3 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo);
    } else if (i==ncx-1) {
        T1 num = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
        2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo) +
        nodes[ (k*(ncy+1)+j)*(ncx+1)+i-2 ].getTT(threadNo);
        num *= num;
        T1 den = nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo) -
        2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
        nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo);
        den *= den;
        T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
        T1 w = 1./(1.+2.*r*r);
        
        T1 am = (1.-w)*(nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo)-
                        nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo))/(2.*dx) +
        w*(3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
           4.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo) +
           nodes[ (k*(ncy+1)+j)*(ncx+1)+i-2 ].getTT(threadNo))/(2.*dx);
        
        a3 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) - dx*am;
        
        t = nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo); // first order for right
        a3 = a3<t ? a3 : t;

    } else {
        T1 num = nodes[ (k*(ncy+1)+j)*(ncx+1)+i+2 ].getTT(threadNo) -
        2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo) +
        nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
        num *= num;
        T1 den = nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo) -
        2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
        nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo);
        den *= den;
        T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
        T1 w = 1./(1.+2.*r*r);
        
        T1 ap = (1.-w)*(nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo)-
                        nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo))/(2.*dx) +
        w*(-nodes[ (k*(ncy+1)+j)*(ncx+1)+i+2 ].getTT(threadNo) +
           4.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo) -
           3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx);
        
        a3 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) + dx*ap;
        
        num = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
        2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo) +
        nodes[ (k*(ncy+1)+j)*(ncx+1)+i-2 ].getTT(threadNo);
        num *= num;
        r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
        w = 1./(1.+2.*r*r);
        
        T1 am = (1.-w)*(nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo)-
                        nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo))/(2.*dx) +
        w*(3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
           4.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo) +
           nodes[ (k*(ncy+1)+j)*(ncx+1)+i-2 ].getTT(threadNo))/(2.*dx);
        
        t = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) - dx*am;
        
        a3 = a3<t ? a3 : t;
    }
    
    if ( a1>a2 ) std::swap(a1, a2);
    if ( a1>a3 ) std::swap(a1, a3);
    if ( a2>a3 ) std::swap(a2, a3);
    
    T1 fh = nodes[(k*(ncy+1)+j)*(ncx+1)+i].getNodeSlowness() * dx;
    
    t = a1 + fh;
    if ( t > a2 ) {
        
        t = 0.5*(a1+a2+sqrt(2.*fh*fh - (a1-a2)*(a1-a2)));
        
        if ( t > a3 ) {
            
            t = 1./3. * ((a1 + a2 + a3) + sqrt(-2.*a1*a1 + 2.*a1*a2 -
                                               2.*a2*a2 + 2.*a1*a3 + 2.*a2*a3 -
                                               2.*a3*a3 + 3.*fh*fh));
            
        }
    }
    
    if ( t<nodes[(k*(ncy+1)+j)*(ncx+1)+i].getTT(threadNo) )
        nodes[(k*(ncy+1)+j)*(ncx+1)+i].setTT(t,threadNo);
    
}

template<typename T1, typename T2, typename NODE>
void Grid3Dri<T1,T2,NODE>::initFSM(const std::vector<sxyz<T1>>& Tx,
                                   const std::vector<T1>& t0,
                                   std::vector<bool>& frozen,
                                   const int npts,
                                   const size_t threadNo) const {
    
    for (size_t n=0; n<Tx.size(); ++n) {
        bool found = false;
        for ( long long nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == Tx[n] ) {
                found = true;
                nodes[nn].setTT( t0[n], threadNo );
                frozen[nn] = true;
                
                long long k = nn/((ncy+1)*(ncx+1));
                long long j = (nn-k*(ncy+1)*(ncx+1))/(ncx+1);
                long long i = nn - (k*(ncy+1)+j)*(ncx+1);
                
                for ( long long kk=k-npts; kk<=k+npts; ++kk ) {
                    if ( kk>=0 && kk<=ncz ) {
                        for ( long long jj=j-npts; jj<=j+npts; ++jj ) {
                            if ( jj>=0 && jj<=ncy ) {
                                for ( long long ii=i-npts; ii<=i+npts; ++ii ) {
                                    if ( ii>=0 && ii<=ncx && !(ii==i && jj==j && kk==k) ) {
                                        
                                        size_t nnn = (kk*(ncy+1)+jj)*(ncx+1)+ii;
                                        T1 tt = nodes[nnn].getDistance(Tx[n]) * nodes[nnn].getNodeSlowness();
                                        nodes[nnn].setTT( tt, threadNo );
                                        frozen[nnn] = true;
                                        
                                    }
                                }
                            }
                        }
                    }
                }
                
                break;
            }
        }
        if ( found==false ) {
            
            // find cell where Tx resides
            long long cellNo = getCellNo(Tx[n]);
            
            long long k = cellNo/(ncy*ncx);
            long long j = (cellNo-k*ncy*ncx)/ncx;
            long long i = cellNo - (k*ncy+j)*ncx;
            
            for ( long long kk=k-(npts-1); kk<=k+npts; ++kk ) {
                if ( kk>=0 && kk<=ncz ) {
                    for ( long long jj=j-(npts-1); jj<=j+npts; ++jj ) {
                        if ( jj>=0 && jj<=ncy ) {
                            for ( long long ii=i-(npts-1); ii<=i+npts; ++ii ) {
                                if ( ii>=0 && ii<=ncx && !(ii==i && jj==j && kk==k) ) {
                                    
                                    size_t nnn = (kk*(ncy+1)+jj)*(ncx+1)+ii;
                                    T1 tt = nodes[nnn].getDistance(Tx[n]) * nodes[nnn].getNodeSlowness();
                                    nodes[nnn].setTT( tt, threadNo );
                                    frozen[nnn] = true;
                                    
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}



#endif

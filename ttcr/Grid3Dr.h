//
//  Grid3Dr.h
//  ttcr.v2
//
//  Created by Giroux Bernard on 12-08-15.
//  Copyright (c) 2012 INRS-ETE. All rights reserved.
//

#ifndef __GRID3DR_H__
#define __GRID3DR_H__

#include "Grid3D.h"
#include "Node3Di.h"
#include "Node3Dc.h"

//namespace Grid3D {
//    size_t get_number_of_nodes3D_r(const size_t nx, const size_t ny,
//								   const size_t nz, const size_t nnx,
//								   const size_t nny, const size_t nnz) {
//        return nx*nnx*((ny+1)*(nz+1)) + // // secondary nodes on the edges
//        ny*nny*((nx+1)*(nz+1)) +
//        nz*nnz*((nx+1)*(ny+1)) +
//        // secondary nodes on the faces
//        (nnx*nny)*(nx*ny*(nz+1))+
//        (nnx*nnz)*(nx*nz*(ny+1))+
//        (nny*nnz)*(ny*nz*(nx+1))+
//        // primary nodes
//        (nx+1) * (ny+1) * (nz+1);
//    }
//};

template<typename T1, typename T2>
class Grid3Dr : public Grid3D<T1,T2> {
public:
    Grid3Dr(const T2 nx, const T2 ny, const T2 nz,
           const T1 ddx, const T1 ddy, const T1 ddz,
           const T1 minx, const T1 miny, const T1 minz,
           const T2 nnx, const T2 nny, const T2 nnz,
           const size_t nt=1) :
    nThreads(nt),
    dx(ddx), dy(ddy), dz(ddz),
    xmin(minx), ymin(miny), zmin(minz),
    xmax(minx+nx*ddx), ymax(miny+ny*ddy), zmax(minz+nz*ddz),
    ncx(nx), ncy(ny), ncz(nz),
    nsnx(nnx), nsny(nny), nsnz(nnz)
    { }
    
    virtual ~Grid3Dr() {}
    
    T1 getDx() const { return dx; }
    T1 getDy() const { return dy; }
    T1 getDz() const { return dz; }
    T1 getXmin() const { return xmin; }
    T1 getXmax() const { return xmax; }
    T1 getYmin() const { return ymin; }
    T1 getYmax() const { return ymax; }
    T1 getZmin() const { return zmin; }
    T1 getZmax() const { return zmax; }
    T2 getNcellx() const { return ncx; }
    T2 getNcelly() const { return ncy; }
    T2 getNcellz() const { return ncz; }
    T2 getNsnx() const { return nsnx; }
    T2 getNsny() const { return nsny; }
    T2 getNsnz() const { return nsnz; }

    T2 getNumberOfCells() const { return ncx*ncy*ncz; }
    size_t getNumberOfNodes() const { return 0; }
	
	virtual void setSlowness(const T1 s) {}
    virtual int setSlowness(const std::vector<T1>&) { return 0; }
    
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
    T2 nsnx;                 // number of secondary nodes in x
    T2 nsny;                 // number of secondary nodes in y
    T2 nsnz;                 // number of secondary nodes in z
    
    T2 getCellNo(const sxyz<T1>& pt) const {
        T1 x = xmax-pt.x < small ? xmax-.5*dx : pt.x;
        T1 y = ymax-pt.y < small ? ymax-.5*dy : pt.y;
        T1 z = zmax-pt.z < small ? zmax-.5*dz : pt.z;
        T2 nx = static_cast<T2>( small + (x-xmin)/dx );
        T2 ny = static_cast<T2>( small + (y-ymin)/dy );
        T2 nz = static_cast<T2>( small + (z-zmin)/dz );
//        T2 nn = ;
//        if ( nn > ncx*ncy*ncz ) {
//            std::cout << pt.x << ' ' << pt.y << ' ' << pt.z << '\n';
//            std::cout << x << ' ' << y << ' ' << z << '\n';
//            std::cout << xmax << ' ' << ymax << ' ' << zmax << '\n';
//            std::cout << xmin << ' ' << ymin << ' ' << zmin << '\n';
//            std::cout << dx << ' ' << dy << ' ' << dz << '\n';
//            std::cout << nx << ' ' << ny << ' ' << nz << '\n';
//            std::cout << ncx << ' ' << ncy << ' ' << ncz << '\n';
//        }
        return ny*ncx + nz*(ncx*ncy) + nx;
    }
	
    T2 getCellNo(const Node3Dc<T1,T2>& node) const {
        T1 x = xmax-node.getX() < small ? xmax-.5*dx : node.getX();
        T1 y = ymax-node.getY() < small ? ymax-.5*dy : node.getY();
        T1 z = zmax-node.getZ() < small ? zmax-.5*dz : node.getZ();
        T2 nx = static_cast<T2>( small + (x-xmin)/dx );
        T2 ny = static_cast<T2>( small + (y-ymin)/dy );
        T2 nz = static_cast<T2>( small + (z-zmin)/dz );
        return ny*ncx + nz*(ncx*ncy) + nx;
    }
	
    T2 getCellNo(const Node3Di<T1,T2>& node) const {
        T1 x = xmax-node.getX() < small ? xmax-.5*dx : node.getX();
        T1 y = ymax-node.getY() < small ? ymax-.5*dy : node.getY();
        T1 z = zmax-node.getZ() < small ? zmax-.5*dz : node.getZ();
        T2 nx = static_cast<T2>( small + (x-xmin)/dx );
        T2 ny = static_cast<T2>( small + (y-ymin)/dy );
        T2 nz = static_cast<T2>( small + (z-zmin)/dz );
        return ny*ncx + nz*(ncx*ncy) + nx;
    }

    int check_pts(const std::vector<sxyz<T1> >&) const;

};

template<typename T1, typename T2>
int Grid3Dr<T1,T2>::check_pts(const std::vector<sxyz<T1> >& pts) const {
    
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


#endif

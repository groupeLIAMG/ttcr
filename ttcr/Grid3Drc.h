//
//  Grid3Drc.h
//
//  Created by Bernard Giroux on 08-04-24.
//
//  Modified by Benoit Larouche on 12-07-20
//  	: now support parallel raytracing from many source points
//  	  on the same 3D grid simultaneously, using OpenMP.
//  	  Secondary nodes are placed on every edge and face of the grid cells.
//
//  	  The velocity model is sampled for each cell and is constant inside the
//        cell.
//
//

//
// Copyright (C) 2012 Bernard Giroux, Benoît Larouche.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#ifndef __GRID3DRC_H__
#define __GRID3DRC_H__


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

#include "Grid3D.h"


template<typename T1, typename T2, typename NODE, typename CELL>
class Grid3Drc : public Grid3D<T1,T2> {
public:
    
    Grid3Drc(const T2 nx, const T2 ny, const T2 nz,
             const T1 ddx, const T1 ddy, const T1 ddz,
             const T1 minx, const T1 miny, const T1 minz,
             const size_t nt=1) :
    nThreads(nt),
    dx(ddx), dy(ddy), dz(ddz),
    xmin(minx), ymin(miny), zmin(minz),
    xmax(minx+nx*ddx), ymax(miny+ny*ddy), zmax(minz+nz*ddz),
    ncx(nx), ncy(ny), ncz(nz),
    nodes(std::vector<NODE>((nx+1)*(ny+1)*(nz+1), NODE(nt))),
    cells(CELL(nx*ny*nz)),
    neighbors(std::vector<std::vector<T2>>(nx*ny*nz))
    { }
    
    virtual ~Grid3Drc() {}
    
//    T1 getSlowness(const size_t n) const { return slowness[n]; }
    
    void setSlowness(const T1 s) {
        cells.setSlowness( s );
    }
    
    int setSlowness(const T1 *s, const size_t ns) {
        return cells.setSlowness( s, ns );
    }
    
    int setSlowness(const std::vector<T1>& s) {
        return cells.setSlowness( s );
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
        std::ofstream fout( filename );
        
        for ( T2 nk=0, n=0; nk<=ncz; ++nk ) {
            T1 z = zmin + nk*dz;
            for ( T2 nj=0; nj<=ncy; ++nj ) {
                T1 y = ymin + nj*dy;
                for (T2 ni=0; ni<=ncx; ++ni, ++n ){
                    T1 x = xmin + ni*dx;
                    fout << x << "   " << y << "   " << z << "   "
                    << cells.getSlowness(n) << '\n';
                }
            }
        }
        fout.close();
    }
    
    void saveTT(const std::string &, const int, const size_t nt=0,
                const bool vtkFormat=0) const;
    
//    size_t getSlownessSize() const {
//        return slowness.size()*sizeof(T1);
//    }
    
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
    
    mutable std::vector<NODE> nodes;
    
    CELL cells;   // column-wise (z axis) slowness vector of the cells, NOT used by Grid3Dcinterp
    std::vector<std::vector<T2>> neighbors;  // nodes common to a cell

    
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
    
    void buildGridNeighbors();
    
    T1 getTraveltime(const sxyz<T1>& Rx,
                     const std::vector<NODE>& nodes,
                     const size_t threadNo) const;
    
    T1 getTraveltime(const sxyz<T1>& Rx,
                     const std::vector<NODE>& nodes,
                     T2&, T2& , const size_t threadNo) const;
    
private:
    Grid3Drc() {}
    Grid3Drc(const Grid3Drc<T1,T2,NODE,CELL>& g) {}
    Grid3Drc<T1,T2,NODE,CELL>& operator=(const Grid3Drc<T1,T2,NODE,CELL>& g) { return *this; }
    
};




template<typename T1, typename T2, typename NODE, typename CELL>
void Grid3Drc<T1,T2,NODE,CELL>::buildGridNeighbors() {
    
    // Index the neighbors nodes of each cell
    for ( T2 n=0; n<nodes.size(); ++n ) {
        for ( size_t n2=0; n2<nodes[n].getOwners().size(); ++n2) {
            neighbors[ nodes[n].getOwners()[n2] ].push_back(n);
        }
    }
}

template<typename T1, typename T2, typename NODE, typename CELL>
int Grid3Drc<T1,T2,NODE,CELL>::checkPts(const std::vector<sxyz<T1>>& pts) const {
    
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


template<typename T1, typename T2, typename NODE, typename CELL>
T1 Grid3Drc<T1,T2,NODE,CELL>::getTraveltime(const sxyz<T1>& Rx,
                                       const std::vector<NODE>& nodes,
                                       const size_t threadNo) const {
    
    // Calculate and return the traveltime for a Rx point.
    for ( size_t nn=0; nn<nodes.size(); ++nn ) {
        if ( nodes[nn] == Rx ) {
            return nodes[nn].getTT(threadNo);
        }
    }
    size_t cellNo = getCellNo( Rx );
    size_t neibNo = neighbors[cellNo][0];
    T1 dt = cells.computeDt(nodes[neibNo], Rx, cellNo);
    
    T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
    for ( size_t k=1; k< neighbors[cellNo].size(); ++k ) {
        neibNo = neighbors[cellNo][k];
        dt = cells.computeDt(nodes[neibNo], Rx, cellNo);
        if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
            traveltime =  nodes[neibNo].getTT(threadNo)+dt;
        }
    }
    return traveltime;
}

template<typename T1, typename T2, typename NODE, typename CELL>
T1 Grid3Drc<T1,T2,NODE,CELL>::getTraveltime(const sxyz<T1>& Rx,
                                       const std::vector<NODE>& nodes,
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
    T2 cellNo = getCellNo( Rx );
    T2 neibNo = neighbors[cellNo][0];
    T1 dt = cells.computeDt(nodes[neibNo], Rx, cellNo);
    
    T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
    nodeParentRx = neibNo;
    cellParentRx = cellNo;
    for ( size_t k=1; k< neighbors[cellNo].size(); ++k ) {
        neibNo = neighbors[cellNo][k];
        dt = cells.computeDt(nodes[neibNo], Rx, cellNo);
        if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
            traveltime =  nodes[neibNo].getTT(threadNo)+dt;
            nodeParentRx = neibNo;
        }
    }
    return traveltime;
}




template<typename T1, typename T2, typename NODE, typename CELL>
void Grid3Drc<T1,T2,NODE,CELL>::saveTT(const std::string &fname, const int all,
                                  const size_t nt, const bool vtkFormat) const {
    
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


#endif

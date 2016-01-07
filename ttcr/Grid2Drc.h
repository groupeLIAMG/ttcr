//
//  Grid2Drc.h
//  ttcr
//
//  Created by Bernard Giroux on 15-12-23.
//  Copyright © 2015 Bernard Giroux. All rights reserved.
//

#ifndef Grid2Drc_h
#define Grid2Drc_h

#include <algorithm>
#include <cstring>
#include <iostream>
#include <fstream>
#include <queue>
#include <vector>

#ifdef VTK
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkXMLRectilinearGridWriter.h"
#endif

#include "Grid2D.h"

template<typename T1, typename T2, typename NODE>
class Grid2Drc : public Grid2D<T1,T2,sxz<T1>> {
public:
    Grid2Drc(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
               const T1 minx, const T1 minz, const size_t nt=1);
    
    virtual ~Grid2Drc() {
    }

    
    T1 getDx() const { return dx; }
    T1 getDz() const { return dz; }
    T1 getXmin() const { return xmin; }
    T1 getXmax() const { return xmax; }
    T1 getZmin() const { return zmin; }
    T1 getZmax() const { return zmax; }
    T2 getNcellx() const { return nCellx; }
    T2 getNcellz() const { return nCellz; }

    T2 getNumberOfCells() const { return slowness.size(); }
    virtual T1 getSlowness(const size_t n) const { return slowness[n]; }
    
    virtual void setSlowness(const T1 s) {
        for ( size_t n=0; n<slowness.size(); ++n ) {
            slowness[n] = s;
        }
    }
    
    virtual int setSlowness(const T1 *s, const size_t ns) {
        if ( slowness.size() != ns ) {
            std::cerr << "Error: slowness vectors of incompatible size.";
            return 1;
        }
        for ( size_t n=0; n<slowness.size(); ++n ) {
            slowness[n] = s[n];
        }
        return 0;
    }
    
    virtual int setSlowness(const std::vector<T1>& s) {
        if ( slowness.size() != s.size() ) {
            std::cerr << "Error: slowness vectors of incompatible size.";
            return 1;
        }
        for ( size_t n=0; n<slowness.size(); ++n ) {
            slowness[n] = s[n];
        }
        return 0;
    }
    
    virtual int raytrace(const std::vector<sxz<T1>>& Tx,
                         const std::vector<T1>& t0,
                         const std::vector<sxz<T1>>& Rx,
                         std::vector<T1>& traveltimes,
                         const size_t threadNo=0) const { return 0; }
    
    virtual int raytrace(const std::vector<sxz<T1>>& Tx,
                         const std::vector<T1>& t0,
                         const std::vector<const std::vector<sxz<T1>>*>& Rx,
                         std::vector<std::vector<T1>*>& traveltimes,
                         const size_t threadNo=0) const { return 0; }
    
    virtual int raytrace(const std::vector<sxz<T1>>& Tx,
                         const std::vector<T1>& t0,
                         const std::vector<sxz<T1>>& Rx,
                         std::vector<T1>& traveltimes,
                         std::vector<std::vector<sxz<T1>>>& r_data,
                         const size_t threadNo=0) const { return 0; }
    
    virtual int raytrace(const std::vector<sxz<T1>>& Tx,
                         const std::vector<T1>& t0,
                         const std::vector<const std::vector<sxz<T1>>*>& Rx,
                         std::vector<std::vector<T1>*>& traveltimes,
                         std::vector<std::vector<std::vector<sxz<T1>>>*>& r_data,
                         const size_t threadNo=0) const { return 0; }

    size_t getNumberOfNodes() const { return nodes.size(); }
    
    T2 getCellNo(const sxz<T1>& pt) const {
        T1 x = xmax-pt.x < small ? xmax-.5*dx : pt.x;
        T1 z = zmax-pt.z < small ? zmax-.5*dz : pt.z;
        T2 nx = static_cast<T2>( small + (x-xmin)/dx );
        T2 nz = static_cast<T2>( small + (z-zmin)/dz );
        return nx*nCellz + nz;
    }
    
    void saveSlownessXYZ(const char filename[]) const {
        std::ofstream fout( filename );
        
        for ( T2 j=0, n=0; j<nCellx; ++j ) {
            T1 x = xmin + (0.5+j)*dx;
            for ( T2 i=0; i<nCellz; ++i, ++n ) {
                T1 z = zmin + (0.5+i)*dz;
                fout << x << "   " << z << "   " << slowness[n] << '\n';
            }
        }
        
        fout.close();
    }

    void saveTT(const std::string &, const int, const size_t nt=0,
                const bool vtkFormat=0) const;
    
    void saveTTgrad(const std::string &, const size_t nt=0,
                    const bool vtkFormat=0) const {}


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

    mutable std::vector<NODE> nodes;
    
    std::vector<T1> slowness;   // column-wise (z axis) slowness vector of the cells
    std::vector<std::vector<T2>> neighbors;  // nodes common to a cell

    void buildGridNeighbors();

    T1 computeDt(const NODE& source, const sxz<T1>& node,
                         const T2 cellNo) const {
        return slowness[cellNo] * source.getDistance( node );
    }
    
    T1 computeDt(const NODE& source, const NODE& node,
                         const T2 cellNo) const {
        return slowness[cellNo] * source.getDistance( node );
    }
    
    int check_pts(const std::vector<sxz<T1>>&) const;

    bool inPolygon(const sxz<T1>& p, const sxz<T1> poly[], const size_t N) const;
    

};

template<typename T1, typename T2, typename NODE>
Grid2Drc<T1,T2,NODE>::Grid2Drc(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
                               const T1 minx, const T1 minz, const size_t nt) : nThreads(nt),
dx(ddx), dz(ddz), xmin(minx), zmin(minz), xmax(minx+nx*ddx), zmax(minz+nz*ddz),
nCellx(nx), nCellz(nz),
nodes(std::vector<NODE>( (nCellx+1) * (nCellz+1), NODE(nt) )),
slowness(std::vector<T1>(nCellx*nCellz)),
neighbors(std::vector<std::vector<T2>>(nCellx*nCellz))
{ }

template<typename T1, typename T2, typename NODE>
void Grid2Drc<T1,T2,NODE>::buildGridNeighbors() {
    for ( T2 n=0; n<nodes.size(); ++n ) {
        for ( size_t n2=0; n2<nodes[n].getOwners().size(); ++n2) {
            neighbors[ nodes[n].getOwners()[n2] ].push_back(n);
        }
    }
}

template<typename T1, typename T2, typename NODE>
void Grid2Drc<T1,T2,NODE>::saveTT(const std::string& fname, const int all,
                                  const size_t nt,
                                  const bool vtkFormat) const {
    
    if (vtkFormat) {
#ifdef VTK
        
        std::string filename = fname+".vtr";
        int nn[3] = {static_cast<int>(nCellx+1), 1, static_cast<int>(nCellz+1)};
        
        vtkSmartPointer<vtkDoubleArray> xCoords = vtkSmartPointer<vtkDoubleArray>::New();
        for (size_t n=0; n<nn[0]; ++n)
            xCoords->InsertNextValue( xmin + n*dx );
        vtkSmartPointer<vtkDoubleArray> yCoords = vtkSmartPointer<vtkDoubleArray>::New();
        yCoords->InsertNextValue( 0.0 );
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
            vtkIdType id = rgrid->FindPoint(nodes[n].getX(), 0.0, nodes[n].getZ());
            newScalars->SetTuple1(id, nodes[n].getTT(nt) );
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
            fout << nodes[n].getX() << '\t'
            << nodes[n].getZ() << '\t'
            << nodes[n].getTT(nt) << '\n';
        }
        fout.close();
    }
}


template<typename T1, typename T2, typename NODE>
int Grid2Drc<T1,T2,NODE>::check_pts(const std::vector<sxz<T1>>& pts) const {
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


template<typename T1, typename T2, typename NODE>
bool Grid2Drc<T1,T2,NODE>::inPolygon(const sxz<T1>& p, const sxz<T1> poly[], const size_t N) const {
    bool c = false;
    for (size_t i = 0, j = N-1; i < N; j = i++) {
        if ((((poly[i].z <= p.z) && (p.z < poly[j].z)) ||
             ((poly[j].z <= p.z) && (p.z < poly[i].z))) &&
            (p.x < (poly[j].x - poly[i].x) * (p.z - poly[i].z) / (poly[j].z - poly[i].z) + poly[i].x))
            c = !c;
    }
    return c;
}





#endif /* Grid2Drc_h */

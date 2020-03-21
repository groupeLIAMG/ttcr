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
#include <exception>
#include <iostream>
#include <fstream>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifdef VTK
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkXMLRectilinearGridWriter.h"
#endif

#include "Grid2D.h"

namespace ttcr {
    
    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    class Grid2Drc : public Grid2D<T1,T2,S> {
    public:
        Grid2Drc(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
                 const T1 minx, const T1 minz, const size_t nt=1);
        
        virtual ~Grid2Drc() {
        }
        
        virtual void setSlowness(const std::vector<T1>& s) {
            try {
                cells.setSlowness( s );
            } catch (std::exception& e) {
                throw;
            }
        }
        void setXi(const std::vector<T1>& x) {
            try {
                cells.setXi( x );
            } catch (std::exception& e) {
                throw;
            }
        }
        void setTiltAngle(const std::vector<T1>& t) {
            try {
                cells.setTiltAngle( t );
            } catch (std::exception& e) {
                throw;
            }
        }
        void setVp0(const std::vector<T1>& s) {
            try {
                cells.setVp0(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        void setVs0(const std::vector<T1>& s) {
            try {
                cells.setVs0(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        void setDelta(const std::vector<T1>& s) {
            try {
                cells.setDelta(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        void setEpsilon(const std::vector<T1>& s) {
            try {
                cells.setEpsilon(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        void setGamma(const std::vector<T1>& s) {
            try {
                cells.setGamma(s);
            } catch (std::exception& e) {
                throw;
            }
        }

        size_t getNumberOfNodes(const bool primary=false) const {
            if ( primary ) {
                return (ncx+1) * (ncz+1);
            } else {
                return nodes.size();
            }
        }
        size_t getNumberOfCells() const { return ncx*ncz; }
        void getTT(std::vector<T1>& tt, const size_t threadNo=0) const final {
            size_t nPrimary = (ncx+1) * (ncz+1);
            tt.resize(nPrimary);
            for ( size_t n=0; n<nPrimary; ++n ) {
                tt[n] = nodes[n].getTT(threadNo);
            }
        }

        T2 getCellNo(const S& pt) const {
            T1 x = xmax-pt.x < small ? xmax-.5*dx : pt.x;
            T1 z = zmax-pt.z < small ? zmax-.5*dz : pt.z;
            T2 nx = static_cast<T2>( small + (x-xmin)/dx );
            T2 nz = static_cast<T2>( small + (z-zmin)/dz );
            return nx*ncz + nz;
        }
        
        void saveSlownessXYZ(const char filename[]) const {
            std::ofstream fout( filename );
            
            for ( T2 j=0, n=0; j<ncx; ++j ) {
                T1 x = xmin + (0.5+j)*dx;
                for ( T2 i=0; i<ncz; ++i, ++n ) {
                    T1 z = zmin + (0.5+i)*dz;
                    fout << x << "   " << z << "   " << cells.getSlowness(n) << '\n';
                }
            }
            
            fout.close();
        }
        
        void saveTT(const std::string &, const int, const size_t nt=0,
                    const int format=1) const;
        
        void saveTTgrad(const std::string &, const size_t nt=0,
                        const bool vtkFormat=0) const {}
        
        const T1 getXmin() const { return xmin; }
        const T1 getZmin() const { return zmin; }
        const T1 getDx() const { return dx; }
        const T1 getDz() const { return dz; }
        const T2 getNcx() const { return ncx; }
        const T2 getNcz() const { return ncz; }
        
        void dump_secondary(std::ofstream& os) const {
            size_t nPrimary = (ncx+1) * (ncz+1);
            for ( size_t n=nPrimary; n<nodes.size(); ++n ) {
                os << nodes[n].getX() << ' ' << nodes[n].getZ() << '\n';
            }
        }

    protected:
        T1 dx;           // cell size in x
        T1 dz;           // cell size in z
        T1 xmin;         // x origin of the grid
        T1 zmin;         // z origin of the grid
        T1 xmax;         // x end of the grid
        T1 zmax;         // z end of the grid
        T2 ncx;          // number of cells in x
        T2 ncz;          // number of cells in z
        
        mutable std::vector<NODE> nodes;
        
        CELL cells;   // column-wise (z axis) slowness vector of the cells
               
        void checkPts(const std::vector<S>&) const;
        
        bool inPolygon(const S& p, const S poly[], const size_t N) const;
        
        
    };
    
    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    Grid2Drc<T1,T2,S,NODE,CELL>::Grid2Drc(const T2 nx, const T2 nz,
                                          const T1 ddx, const T1 ddz,
                                          const T1 minx, const T1 minz,
                                          const size_t nt) :
    Grid2D<T1,T2,S>(nx*nz, nt),
    dx(ddx), dz(ddz), xmin(minx), zmin(minz),
    xmax(minx+nx*ddx), zmax(minz+nz*ddz),
    ncx(nx), ncz(nz),
    nodes(std::vector<NODE>( (ncx+1) * (ncz+1), NODE(nt) )),
    cells(ncx*ncz)
    { }
    
    
    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Drc<T1,T2,S,NODE,CELL>::saveTT(const std::string& fname,
                                             const int all,
                                             const size_t nt,
                                             const int format) const {
        
        if ( format == 1 ) {
            std::string filename = fname+".dat";
            std::ofstream fout(filename.c_str());
            fout.precision(12);
            for ( T2 n=0; n<nodes.size(); ++n ) {
                fout << nodes[n].getX() << '\t'
                << nodes[n].getZ() << '\t'
                << nodes[n].getTT(nt) << '\n';
            }
            fout.close();
        } else  if ( format == 2 ) {
#ifdef VTK
            
            std::string filename = fname+".vtr";
            int nn[3] = {static_cast<int>(ncx+1), 1, static_cast<int>(ncz+1)};
            
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
        } else if ( format == 3 ) {
            std::string filename = fname+".bin";
            std::ofstream fout(filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
            for ( T2 n=0; n<nodes.size(); ++n ) {
                T1 tmp[] = { nodes[n].getX(), nodes[n].getZ(), nodes[n].getTT(nt) };
                fout.write( (char*)tmp, 3*sizeof(T1) );
            }
            fout.close();
        } else {
            throw std::runtime_error("Unsupported format for saving traveltimes");
        }
    }
    
    
    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Drc<T1,T2,S,NODE,CELL>::checkPts(const std::vector<S>& pts) const {
        for (size_t n=0; n<pts.size(); ++n) {
            if ( pts[n].x < xmin || pts[n].x > xmax ||
                pts[n].z < zmin || pts[n].z > zmax ) {
                std::ostringstream msg;
                msg << "Error: Point (" << pts[n].x << ", "<< pts[n] .z << ") outside grid.";
                throw std::runtime_error(msg.str());
            }
        }
    }
    
    
    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    bool Grid2Drc<T1,T2,S,NODE,CELL>::inPolygon(const S& p, const S poly[],
                                                const size_t N) const {
        bool c = false;
        for (size_t i = 0, j = N-1; i < N; j = i++) {
            if ((((poly[i].z <= p.z) && (p.z < poly[j].z)) ||
                 ((poly[j].z <= p.z) && (p.z < poly[i].z))) &&
                (p.x < (poly[j].x - poly[i].x) * (p.z - poly[i].z) / (poly[j].z - poly[i].z) + poly[i].x))
                c = !c;
        }
        return c;
    }
    
}



#endif /* Grid2Drc_h */

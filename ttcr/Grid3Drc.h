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

#include "Grid3Dr.h"
//#include "Node3Dcsp.h"


template<typename T1, typename T2, typename NODE>
class Grid3Drc : public Grid3Dr<T1,T2> {
public:
    
    Grid3Drc(const T2 nx, const T2 ny, const T2 nz,
             const T1 ddx, const T1 ddy, const T1 ddz,
             const T1 minx, const T1 miny, const T1 minz,
             const T2 nnx, const T2 nny, const T2 nnz,
             const size_t nt=1) :
    Grid3Dr<T1,T2>(nx, ny, nz, ddx, ddy, ddz, minx, miny, minz, nnx, nny, nnz, nt),
    nodes(std::vector<NODE>(// secondary nodes on the edges
                            nx*nnx*((ny+1)*(nz+1)) +
                            ny*nny*((nx+1)*(nz+1)) +
                            nz*nnz*((nx+1)*(ny+1)) +
                            // secondary nodes on the faces
                            (nnx*nny)*(nx*ny*(nz+1))+
                            (nnx*nnz)*(nx*nz*(ny+1))+
                            (nny*nnz)*(ny*nz*(nx+1))+
                            // primary nodes
                            (nx+1) * (ny+1) * (nz+1),
                            NODE(nt) )),
    slowness(std::vector<T1>(nx*ny*nz)),
    neighbors(std::vector<std::vector<T2>>(nx*ny*nz))
    { }
    
    virtual ~Grid3Drc() {}
    
    T1 getSlowness(const size_t n) const { return slowness[n]; }
    
    void setSlowness(const T1 s) {
        for ( size_t n=0; n<slowness.size(); ++n ) {
            slowness[n] = s;
        }
    }
    
    int setSlowness(const T1 *s, const size_t ns) {
        if ( slowness.size() != ns ) {
            std::cerr << "Error: slowness vectors of incompatible size.";
            return 1;
        }
        for ( size_t n=0; n<slowness.size(); ++n ) {
            slowness[n] = s[n];
        }
        return 0;
    }
    
    int setSlowness(const std::vector<T1>& s) {
        if ( slowness.size() != s.size() ) {
            std::cerr << "Error: slowness vectors of incompatible size.";
            return 1;
        }
        for ( size_t n=0; n<slowness.size(); ++n ) {
            slowness[n] = s[n];
        }
        return 0;
    }
    
    size_t getNumberOfNodes() const { return nodes.size(); }
    
    virtual int raytrace(const std::vector<sxyz<T1>>& Tx,
                         const std::vector<T1>& t0,
                         const std::vector<sxyz<T1>>& Rx,
                         std::vector<T1>& traveltimes,
                         const size_t threadNo=0) const { return 0; }
    
    virtual int raytrace(const std::vector<sxyz<T1>>&,
                         const std::vector<T1>&,
                         const std::vector<const std::vector<sxyz<T1>>*>&,
                         std::vector<std::vector<T1>*>&,
                         const size_t=0) const { return 0; }
    
    virtual int raytrace(const std::vector<sxyz<T1>>& Tx,
                         const std::vector<T1>& t0,
                         const std::vector<sxyz<T1>>& Rx,
                         std::vector<T1>& traveltimes,
                         std::vector<std::vector<sxyz<T1>>>& r_data,
                         const size_t threadNo=0) const { return 0; }
    
    virtual int raytrace(const std::vector<sxyz<T1>>&,
                         const std::vector<T1>&,
                         const std::vector<const std::vector<sxyz<T1>>*>&,
                         std::vector<std::vector<T1>*>&,
                         std::vector<std::vector<std::vector<sxyz<T1>>>*>&,
                         const size_t=0) const { return 0; }
    
    virtual int raytrace2(const std::vector<sxyz<T1>>& Tx,
                          const std::vector<T1>& t0,
                          const std::vector<sxyz<T1>>& Rx,
                          std::vector<T1>& traveltimes,
                          const size_t threadNo=0) const { return 0; }
    
    
    void saveSlownessXYZ(const char filename[]) const {
        std::ofstream fout( filename );
        
        for ( T2 nk=0, n=0; nk<=Grid3Dr<T1,T2>::ncz; ++nk ) {
            T1 z = Grid3Dr<T1,T2>::zmin + nk*Grid3Dr<T1,T2>::dz;
            for ( T2 nj=0; nj<=Grid3Dr<T1,T2>::ncy; ++nj ) {
                T1 y = Grid3Dr<T1,T2>::ymin + nj*Grid3Dr<T1,T2>::dy;
                for (T2 ni=0; ni<=Grid3Dr<T1,T2>::ncx; ++ni, ++n ){
                    T1 x = Grid3Dr<T1,T2>::xmin + ni*Grid3Dr<T1,T2>::dx;
                    fout << x << "   " << y << "   " << z << "   "
                    << slowness[n] << '\n';
                }
            }
        }
        fout.close();
    }
    
    void save(const char filename[]) const;
    void dsave(const char filename[]) const;
    void savefast(const char filename[]) const;
    void savePrimary(const char filename[], const size_t nt=0,
                     const bool vtkFormat=0) const;
    
    void saveTT(const std::string &, const int, const size_t nt=0,
                const bool vtkFormat=0) const;
    
    size_t getSlownessSize() const {
        return slowness.size()*sizeof(T1);
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
    
protected:
    
    mutable std::vector<NODE> nodes;
    
    std::vector<T1> slowness;   // column-wise (z axis) slowness vector of the cells, NOT used by Grid3Dcinterp
    std::vector<std::vector<T2>> neighbors;  // nodes common to a cell
    
    void buildGridNeighbors();
    
    
    
    
    T1 computeDt(const NODE& source, const sxyz<T1>& node,
                 const size_t cellNo) const {
        return slowness[cellNo] * source.getDistance( node );
    }
    
    T1 computeDt(const NODE& source, const NODE& node,
                 const size_t cellNo) const {
        return slowness[cellNo] * source.getDistance( node );
    }
    
    T1 getTraveltime(const sxyz<T1>& Rx,
                     const std::vector<NODE>& nodes,
                     const size_t threadNo) const;
    
    T1 getTraveltime(const sxyz<T1>& Rx,
                     const std::vector<NODE>& nodes,
                     T2&, T2& , const size_t threadNo) const;
    
    Grid3Drc() {}
    Grid3Drc(const Grid3Drc<T1,T2,NODE>& g) {}
    Grid3Drc<T1,T2,NODE>& operator=(const Grid3Drc<T1,T2,NODE>& g) { return *this; }
    
};




template<typename T1, typename T2, typename NODE>
void Grid3Drc<T1,T2,NODE>::buildGridNeighbors() {
    
    // Index the neighbors nodes of each cell
    for ( T2 n=0; n<nodes.size(); ++n ) {
        for ( size_t n2=0; n2<nodes[n].getOwners().size(); ++n2) {
            neighbors[ nodes[n].getOwners()[n2] ].push_back(n);
        }
    }
}



template<typename T1, typename T2, typename NODE>
T1 Grid3Drc<T1,T2,NODE>::getTraveltime(const sxyz<T1>& Rx,
                                       const std::vector<NODE>& nodes,
                                       const size_t threadNo) const {
    
    // Calculate and return the traveltime for a Rx point.
    for ( size_t nn=0; nn<nodes.size(); ++nn ) {
        if ( nodes[nn] == Rx ) {
            return nodes[nn].getTT(threadNo);
        }
    }
    size_t cellNo = this->getCellNo( Rx );
    size_t neibNo = neighbors[cellNo][0];
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

template<typename T1, typename T2, typename NODE>
T1 Grid3Drc<T1,T2,NODE>::getTraveltime(const sxyz<T1>& Rx,
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
    T2 cellNo = this->getCellNo( Rx );
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

template<typename T1, typename T2, typename NODE>
void Grid3Drc<T1,T2,NODE>::save(const char filename[]) const {
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

template<typename T1, typename T2, typename NODE>
void Grid3Drc<T1,T2,NODE>::dsave(const char filename[]) const {
    //Similar to 'save', with text information
    std::ofstream fout( filename );
    
    fout << "dx "<< Grid3Dr<T1,T2>::dx << "\t dy " << Grid3Dr<T1,T2>::dy << "\t dz " << Grid3Dr<T1,T2>::dz
    << "\t xmin " << Grid3Dr<T1,T2>::xmin << "\t ymin " << Grid3Dr<T1,T2>::ymin << "\t zmin "	<< Grid3Dr<T1,T2>::zmin
    << "\t xmax " << Grid3Dr<T1,T2>::xmax<< "\t ymax " << Grid3Dr<T1,T2>::ymax << "\t zmax "<< Grid3Dr<T1,T2>::zmax
    << '\n';
    fout << "nCx " << Grid3Dr<T1,T2>::ncx << "\t nCy " << Grid3Dr<T1,T2>::ncy << "\t nCz " << Grid3Dr<T1,T2>::ncz
    << "\t nsnx " << Grid3Dr<T1,T2>::nsnx << "\t nsny " << Grid3Dr<T1,T2>::nsny << "\t nsnz " << Grid3Dr<T1,T2>::nsnz
    << '\n';
    
    fout << "nb. nodes " << nodes.size() << '\n';
    for ( size_t n=0; n < nodes.size(); ++n ) {
        fout << "node " << nodes[n].getGridIndex() << "\t TT \t ";
        for ( size_t nt=0; nt< nodes[n].getsize(); nt++ ) {
            fout << nodes[n].getTT(nt) << "\t";
        }
        fout << " X " << nodes[n].getX() << "\t Y " << nodes[n].getY()
        << "\t Z " << nodes[n].getZ() << "\t Ray Parent \t";
        for ( size_t nt=0; nt< nodes[n].getsize(); nt++ ) {
            fout << nodes[n].getNodeParent(nt) << '\t';
        }
        fout<< "Cell Parent \t";
        for ( size_t nt=0; nt< nodes[n].getsize(); nt++ ) {
            fout << nodes[n].getCellParent(nt) << '\t';
        }
        fout << "Owners: ";
        for ( size_t no=0; no < nodes[n].getOwners().size(); ++no ) {
            fout << '\t' << nodes[n].getOwners()[no];
        }
        fout << '\n';
    }
    /*
     fout << "slowness size " << slowness.size() << '\n';
     for ( size_t n=0; n < slowness.size(); ++n ) {
     fout << slowness[n] << '\n';
     }
     
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

template<typename T1, typename T2, typename NODE>
void Grid3Drc<T1,T2,NODE>::savefast(const char filename[]) const {
    
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


template<typename T1, typename T2, typename NODE>
void Grid3Drc<T1,T2,NODE>::savePrimary(const char filename[], const size_t nt,
                                       const bool vtkFormat) const {
    
    if ( vtkFormat ) {
        
#ifdef VTK
        
        std::string fname = std::string(filename)+".vtr";
        
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
        
        writer->SetFileName( fname.c_str() );
        //		writer->SetInputConnection( rgrid->GetProducerPort() );
        writer->SetInputData( rgrid );
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
                
                for ( size_t ni=0; ni<=Grid3Dr<T1,T2>::ncx; ++ni ) {
                    
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


template<typename T1, typename T2, typename NODE>
void Grid3Drc<T1,T2,NODE>::saveTT(const std::string &fname, const int all,
                                  const size_t nt, const bool vtkFormat) const {
    
    if (vtkFormat) {
#ifdef VTK
        if ( all == 1 )
            std::cout << "Warning, only primary nodes are save in VTK format\n";
        savePrimary(fname.c_str(), nt, vtkFormat);
#else
        std::cerr << "VTK not included during compilation.\nNothing saved.\n";
#endif
    } else {
        if ( all == 0 ) {
            savePrimary(fname.c_str(), nt, vtkFormat);
        } else {
            std::ofstream fout(fname.c_str());
            for ( T2 n=0; n<nodes.size(); ++n ) {
                fout << nodes[n].getX() << '\t'
                << nodes[n].getY() << '\t'
                << nodes[n].getZ() << '\t'
                << nodes[n].getTT(nt) << '\n';
            }
            fout.close();
        }
    }
}


#endif

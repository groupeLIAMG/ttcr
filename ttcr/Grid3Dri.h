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

#ifdef VTK
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkXMLRectilinearGridWriter.h"
#endif


template<typename T1, typename T2, typename NODE>
class Grid3Dri : public Grid3Dr<T1,T2> {
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
             const size_t nt=1) :
    Grid3Dr<T1,T2>(nx, ny, nz, ddx, ddy, ddz, minx, miny, minz, nt),
    nodes(std::vector<NODE>()),
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
    
    virtual int setSlowness(const std::vector<T1>& s) { return 0; }
    
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
    
    void save(const char filename[]) const;
    void dsave(const char filename[]) const;
    void savefast(const char filename[]) const;
    
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
    std::vector<std::vector<T2>> neighbors;  // nodes common to a cell
    
    void buildGridNeighbors();
    
    
    T1 computeDt(const NODE& source, const NODE& node) const {
        return (node.getNodeSlowness()+source.getNodeSlowness())/2. * source.getDistance( node );
    }
    
    T1 computeDt(const NODE& source, const sxyz<T1>& node, T1 slo) const {
        return (slo+source.getNodeSlowness())/2. * source.getDistance( node );
    }
    
    bool isNearInt( double value ) const {
        return ( remainder(value, 1.)  <= small );
    }
    
    Grid3Dri() {}
    Grid3Dri(const Grid3Dri<T1,T2,NODE>& g) {}
    Grid3Dri<T1,T2,NODE>& operator=(const Grid3Dri<T1,T2,NODE>& g) { return *this; }
    
};



template<typename T1, typename T2, typename NODE>
void Grid3Dri<T1,T2,NODE>::buildGridNeighbors() {
    
    //Index the neighbors nodes of each cell
    for ( T2 n=0; n<nodes.size(); ++n ) {
        for ( size_t n2=0; n2<nodes[n].getOwners().size(); ++n2) {
            neighbors[ nodes[n].getOwners()[n2] ].push_back(n);
        }
    }
}



template<typename T1, typename T2, typename NODE>
void Grid3Dri<T1,T2,NODE>::save(const char filename[]) const {
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
void Grid3Dri<T1,T2,NODE>::dsave(const char filename[]) const {
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

template<typename T1, typename T2, typename NODE>
void Grid3Dri<T1,T2,NODE>::savefast(const char filename[]) const {
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


#endif

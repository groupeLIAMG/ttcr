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
             const size_t nt=1, const bool invDist=false) :
    Grid3Dr<T1,T2>(nx, ny, nz, ddx, ddy, ddz, minx, miny, minz, nt),
    nodes(std::vector<NODE>(nx*ny*nz, NODE(nt))),
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
    bool inverseDistance;
    
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
    
    T1 getTraveltime(const sxyz<T1>& Rx,
                     const std::vector<NODE>& nodes,
                     const size_t threadNo) const;
    
    T1 getTraveltime(const sxyz<T1>& Rx,
                     const std::vector<NODE>& nodes,
                     T2&, T2&, const size_t threadNo) const;
    
    void grad(sxyz<T1>& g, const size_t i, const size_t j, const size_t k,
              const size_t nt) const;

    void grad(sxyz<T1>& g, const sxyz<T1> &pt, const size_t nt) const;
    
    T1 interpTT(const sxyz<T1> &pt, const size_t nt) const;

    void getRaypath(const std::vector<sxyz<T1>>& Tx,
                    const sxyz<T1> &Rx,
                    std::vector<sxyz<T1>> &r_data,
                    const size_t threadNo=0) const;

    void getRaypath_old(const std::vector<sxyz<T1>>& Tx,
                        const sxyz<T1> &Rx,
                        std::vector<sxyz<T1>> &r_data,
                        const size_t threadNo=0) const;
    
    T1 computeSlowness(const sxyz<T1>& Rx ) const;

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

template<typename T1, typename T2, typename NODE>
T1 Grid3Dri<T1,T2,NODE>::getTraveltime(const sxyz<T1>& Rx,
                                    const std::vector<NODE>& nodes,
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
    T2 neibNo = this->neighbors[cellNo][0];
    T1 dt = this->computeDt(nodes[neibNo], Rx, slo);
    
    T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
    for ( size_t k=1; k< this->neighbors[cellNo].size(); ++k ) {
        neibNo = this->neighbors[cellNo][k];
        dt = this->computeDt(nodes[neibNo], Rx, slo);
        if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
            traveltime =  nodes[neibNo].getTT(threadNo)+dt;
        }
    }
    return traveltime;
}

template<typename T1, typename T2, typename NODE>
T1 Grid3Dri<T1,T2,NODE>::getTraveltime(const sxyz<T1>& Rx,
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
    //If Rx is not on a node:
    T1 slo = computeSlowness( Rx );
    
    T2 cellNo = this->getCellNo( Rx );
    T2 neibNo = this->neighbors[cellNo][0];
    T1 dt = this->computeDt(nodes[neibNo], Rx, slo);
    
    T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
    nodeParentRx = neibNo;
    cellParentRx = cellNo;
    for ( size_t k=1; k< this->neighbors[cellNo].size(); ++k ) {
        neibNo = this->neighbors[cellNo][k];
        dt = this->computeDt(nodes[neibNo], Rx, slo);
        if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
            traveltime =  nodes[neibNo].getTT(threadNo)+dt;
            nodeParentRx = neibNo;
        }
    }
    return traveltime;
}


template<typename T1, typename T2, typename NODE>
void Grid3Dri<T1,T2,NODE>::grad(sxyz<T1>& g, const size_t i, const size_t j, const size_t k,
                                const size_t nt) const {
    
    // compute average gradient for voxel (i,j,k)

    static const size_t nnx = this->ncx+1;
    static const size_t nny = this->ncy+1;

    g.x = 0.25*(nodes[(    k*nny+j  )*nnx+i+1].getTT(nt) - nodes[(    k*nny+j  )*nnx+i  ].getTT(nt) +
                nodes[(    k*nny+j+1)*nnx+i+1].getTT(nt) - nodes[(    k*nny+j+1)*nnx+i  ].getTT(nt) +
                nodes[((k+1)*nny+j  )*nnx+i+1].getTT(nt) - nodes[((k+1)*nny+j  )*nnx+i  ].getTT(nt) +
                nodes[((k+1)*nny+j+1)*nnx+i+1].getTT(nt) - nodes[((k+1)*nny+j+1)*nnx+i  ].getTT(nt))/this->dx;
    g.y = 0.25*(nodes[(    k*nny+j+1)*nnx+i  ].getTT(nt) - nodes[(    k*nny+j  )*nnx+i  ].getTT(nt) +
                nodes[(    k*nny+j+1)*nnx+i+1].getTT(nt) - nodes[(    k*nny+j  )*nnx+i+1].getTT(nt) +
                nodes[((k+1)*nny+j+1)*nnx+i  ].getTT(nt) - nodes[((k+1)*nny+j  )*nnx+i  ].getTT(nt) +
                nodes[((k+1)*nny+j+1)*nnx+i+1].getTT(nt) - nodes[((k+1)*nny+j  )*nnx+i+1].getTT(nt))/this->dy;
    g.z = 0.25*(nodes[((k+1)*nny+j  )*nnx+i  ].getTT(nt) - nodes[(    k*nny+j  )*nnx+i  ].getTT(nt) +
                nodes[((k+1)*nny+j  )*nnx+i+1].getTT(nt) - nodes[(    k*nny+j  )*nnx+i+1].getTT(nt) +
                nodes[((k+1)*nny+j+1)*nnx+i  ].getTT(nt) - nodes[(    k*nny+j+1)*nnx+i  ].getTT(nt) +
                nodes[((k+1)*nny+j+1)*nnx+i+1].getTT(nt) - nodes[(    k*nny+j+1)*nnx+i+1].getTT(nt))/this->dz;
}


template<typename T1, typename T2, typename NODE>
void Grid3Dri<T1,T2,NODE>::grad(sxyz<T1>& g, const sxyz<T1> &pt,
                                const size_t nt) const {
    
    // compute travel time gradient at point pt
    
    T1 p1 = pt.x - this->dx/2.0;
    T1 p2 = p1 + this->dx;
    g.x = (interpTT({p2, pt.y, pt.z}, nt) - interpTT({p1, pt.y, pt.z}, nt)) / this->dx;

    p1 = pt.y - this->dy/2.0;
    p2 = p1 + this->dy;
    g.y = (interpTT({pt.x, p2, pt.z}, nt) - interpTT({pt.x, p1, pt.z}, nt)) / this->dy;

    p1 = pt.z - this->dz/2.0;
    p2 = p1 + this->dz;
    g.z = (interpTT({pt.x, pt.y, p2}, nt) - interpTT({pt.x, pt.y, p1}, nt)) / this->dz;

}


template<typename T1, typename T2, typename NODE>
T1 Grid3Dri<T1,T2,NODE>::interpTT(const sxyz<T1> &pt, const size_t nt) const {
    
    static const size_t nnx = this->ncx+1;
    static const size_t nny = this->ncy+1;

    // trilinear interpolation if not on node

    T1 tt;
    T2 i, j, k;
    
    this->getIJK(pt, i, j, k);

    if ( fabs(pt.x - (this->xmin+i*this->dx))<small &&
        fabs(pt.y - (this->ymin+j*this->dy))<small &&
        fabs(pt.z - (this->zmin+k*this->dz))<small ) {
        // on node
        return nodes[(k*nny+j)*nnx+i].getTT(nt);
    } else if ( fabs(pt.x - (this->xmin+i*this->dx))<small &&
               fabs(pt.y - (this->ymin+j*this->dy))<small ) {
        // on edge
        T1 t1 = nodes[(    k*nny+j)*nnx+i].getTT(nt);
        T1 t2 = nodes[((k+1)*nny+j)*nnx+i].getTT(nt);
        
		T1 w1 = (this->zmin+(k+1)*this->dz - pt.z)/this->dz;
		T1 w2 = (pt.z - (this->zmin+k*this->dz))/this->dz;
		
		tt = t1*w1 + t2*w2;
		
    } else if ( fabs(pt.x - (this->xmin+i*this->dx))<small &&
               fabs(pt.z - (this->zmin+k*this->dz))<small ) {
        // on edge
        T1 t1 = nodes[(k*nny+j  )*nnx+i].getTT(nt);
        T1 t2 = nodes[(k*nny+j+1)*nnx+i].getTT(nt);
        
        T1 w1 = (this->ymin+(j+1)*this->dy - pt.y)/this->dy;
		T1 w2 = (pt.y - (this->ymin+j*this->dy))/this->dy;
		
        tt = t1*w1 + t2*w2;

    } else if ( fabs(pt.y - (this->ymin+j*this->dy))<small &&
               fabs(pt.z - (this->zmin+k*this->dz))<small ) {
        // on edge
        T1 t1 = nodes[(k*nny+j)*nnx+i  ].getTT(nt);
        T1 t2 = nodes[(k*nny+j)*nnx+i+1].getTT(nt);
        
		T1 w1 = (this->xmin+(i+1)*this->dx - pt.x)/this->dx;
        T1 w2 = (pt.x - (this->xmin+i*this->dx))/this->dx;
		
        tt = t1*w1 + t2*w2;

    } else if ( fabs(pt.x - (this->xmin+i*this->dx))<small ) {
        // on YZ face
        T1 t1 = nodes[(    k*nny+j  )*nnx+i].getTT(nt);
        T1 t2 = nodes[((k+1)*nny+j  )*nnx+i].getTT(nt);
        T1 t3 = nodes[(    k*nny+j+1)*nnx+i].getTT(nt);
        T1 t4 = nodes[((k+1)*nny+j+1)*nnx+i].getTT(nt);

		T1 w1 = (this->zmin+(k+1)*this->dz - pt.z)/this->dz;
		T1 w2 = (pt.z - (this->zmin+k*this->dz))/this->dz;

		t1 = t1*w1 + t2*w2;
		t2 = t3*w1 + t4*w2;
		
		w1 = (this->ymin+(j+1)*this->dy - pt.y)/this->dy;
		w2 = (pt.y - (this->ymin+j*this->dy))/this->dy;
		
		tt = t1*w1 + t2*w2;
		
    } else if ( fabs(pt.y - (this->ymin+j*this->dy))<small ) {
        // on XZ face
        T1 t1 = nodes[(    k*nny+j)*nnx+i  ].getTT(nt);
        T1 t2 = nodes[((k+1)*nny+j)*nnx+i  ].getTT(nt);
        T1 t3 = nodes[(    k*nny+j)*nnx+i+1].getTT(nt);
        T1 t4 = nodes[((k+1)*nny+j)*nnx+i+1].getTT(nt);
        
		T1 w1 = (this->zmin+(k+1)*this->dz - pt.z)/this->dz;
		T1 w2 = (pt.z - (this->zmin+k*this->dz))/this->dz;
		
		t1 = t1*w1 + t2*w2;
		t2 = t3*w1 + t4*w2;

		w1 = (this->xmin+(i+1)*this->dx - pt.x)/this->dx;
		w2 = (pt.x - (this->xmin+i*this->dx))/this->dx;
		
		tt = t1*w1 + t2*w2;

	} else if ( fabs(pt.z - (this->zmin+k*this->dz))<small ) {
        // on XY face
        T1 t1 = nodes[(k*nny+j  )*nnx+i  ].getTT(nt);
        T1 t2 = nodes[(k*nny+j+1)*nnx+i  ].getTT(nt);
        T1 t3 = nodes[(k*nny+j  )*nnx+i+1].getTT(nt);
        T1 t4 = nodes[(k*nny+j+1)*nnx+i+1].getTT(nt);

		T1 w1 = (this->ymin+(j+1)*this->dy - pt.y)/this->dy;
		T1 w2 = (pt.y - (this->ymin+j*this->dy))/this->dy;

		t1 = t1*w1 + t2*w2;
		t2 = t3*w1 + t4*w2;
		
		w1 = (this->xmin+(i+1)*this->dx - pt.x)/this->dx;
		w2 = (pt.x - (this->xmin+i*this->dx))/this->dx;
		
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

		T1 w1 = (this->zmin+(k+1)*this->dz - pt.z)/this->dz;
		T1 w2 = (pt.z - (this->zmin+k*this->dz))/this->dz;
		
		t1 = t1*w1 + t2*w2;
		t2 = t3*w1 + t4*w2;
		t3 = t5*w1 + t6*w2;
		t4 = t7*w1 + t8*w2;
		
		w1 = (this->ymin+(j+1)*this->dy - pt.y)/this->dy;
		w2 = (pt.y - (this->ymin+j*this->dy))/this->dy;
		
		t1 = t1*w1 + t2*w2;
		t2 = t3*w1 + t4*w2;
		
		w1 = (this->xmin+(i+1)*this->dx - pt.x)/this->dx;
		w2 = (pt.x - (this->xmin+i*this->dx))/this->dx;
		
		tt = t1*w1 + t2*w2;
		
    }
    
    return tt;
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
    static const T1 maxDist = sqrt( this->dx*this->dx + this->dy*this->dy + this->dz*this->dz );
    sxyz<T1> g;
    
    bool reachedTx = false;
    while ( reachedTx == false ) {
        
        grad(g, curr_pt, threadNo);
        g *= -1.0;
        
        long long i, j, k;
        this->getIJK(curr_pt, i, j, k);
        
        // planes we will intersect
        T1 xp = this->xmin + this->dx*(i + boost::math::sign(g.x)>0.0 ? 1.0 : 0.0);
        T1 yp = this->ymin + this->dy*(j + boost::math::sign(g.y)>0.0 ? 1.0 : 0.0);
        T1 zp = this->zmin + this->dz*(k + boost::math::sign(g.z)>0.0 ? 1.0 : 0.0);
        
        if ( fabs(xp-curr_pt.x)<small) {
            xp += this->dx*boost::math::sign(g.x);
        }
        if ( fabs(yp-curr_pt.y)<small) {
            yp += this->dy*boost::math::sign(g.y);
        }
        if ( fabs(zp-curr_pt.z)<small) {
            zp += this->dz*boost::math::sign(g.z);
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
        
        if ( curr_pt.x < this->xmin || curr_pt.x > this->xmax ||
            curr_pt.y < this->ymin || curr_pt.y > this->ymax ||
            curr_pt.z < this->zmin || curr_pt.z > this->zmax ) {
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
    static const T1 maxDist = sqrt( this->dx*this->dx + this->dy*this->dy + this->dz*this->dz );

    this->getIJK(curr_pt, iIn, jIn, kIn);
    
    bool reachedTx = false;
    while ( reachedTx == false ) {
        
        bool onNode=false;
        bool onEdgeX=false;
        bool onEdgeY=false;
        bool onEdgeZ=false;
        
        if ( fabs(remainder(curr_pt.x,this->dx))<small &&
            fabs(remainder(curr_pt.y,this->dy))<small &&
            fabs(remainder(curr_pt.z,this->dz))<small ) {
            onNode = true;
        } else if ( fabs(remainder(curr_pt.y,this->dy))<small &&
                   fabs(remainder(curr_pt.z,this->dz))<small ) {
            onEdgeX = true;
        } else if ( fabs(remainder(curr_pt.x,this->dx))<small &&
                   fabs(remainder(curr_pt.z,this->dz))<small ) {
            onEdgeY = true;
        } else if ( fabs(remainder(curr_pt.x,this->dx))<small &&
                   fabs(remainder(curr_pt.y,this->dy))<small ) {
            onEdgeZ = true;
        }
        
        if ( onNode ) {
            
            T2 i, j, k;
            this->getIJK(curr_pt, i, j, k);
            std::vector<sijk<T2>> voxels;

            // find voxels touching node
            if ( i<=this->ncx && j<=this->ncy && k<=this->ncz )
                voxels.push_back( {i,j,k} );
            if ( i<=this->ncx && j<=this->ncy && k>0 )
                voxels.push_back( {i,j,k-1} );
            if ( i<=this->ncx && j>0 && k<=this->ncz )
                voxels.push_back( {i,j-1,k} );
            if ( i<=this->ncx && j>0 && k>0 )
                voxels.push_back( {i,j-1,k-1} );
            if ( i>0 && j<=this->ncy && k<=this->ncz )
                voxels.push_back( {i-1,j,k} );
            if ( i>0 && j<=this->ncy && k>0 )
                voxels.push_back( {i-1,j,k-1} );
            if ( i>0 && j>0 && k<=this->ncz )
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
            
            if ((gOut.x<0.0 && i==0) || (gOut.x>0.0 && i==this->ncx+1) ||
                (gOut.y<0.0 && j==0) || (gOut.y>0.0 && j==this->ncy+1) ||
                (gOut.z<0.0 && k==0) || (gOut.z>0.0 && k==this->ncz+1)) {
                //  we are going oustide the grid!
                std::cerr << "Error while computing raypaths: going outside grid!\n"
                << "  Stopping calculations, raypaths will be incomplete.\n" << std::endl;
                return;
            }
            
            iOut = boost::math::sign(gOut.x)<0.0 ? i-1 : i;
            jOut = boost::math::sign(gOut.y)<0.0 ? j-1 : j;
            kOut = boost::math::sign(gOut.z)<0.0 ? k-1 : k;

            // planes we will intersect
            T1 xp = this->xmin + this->dx*(i + boost::math::sign(gOut.x));
            T1 yp = this->ymin + this->dy*(j + boost::math::sign(gOut.y));
            T1 zp = this->zmin + this->dz*(k + boost::math::sign(gOut.z));
            
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
            this->getIJK(curr_pt, i, j, k);
            std::vector<sijk<T2>> voxels;
            
            // find voxels touching edge
            if ( i<=this->ncx && j<=this->ncy && k<=this->ncz )
                voxels.push_back( {i,j,k} );
            if ( i<=this->ncx && j<=this->ncy && k>0 )
                voxels.push_back( {i,j,k-1} );
            if ( i<=this->ncx && j>0 && k<=this->ncz )
                voxels.push_back( {i,j-1,k} );
            if ( i<=this->ncx && j>0 && k>0 )
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

            if ((gOut.x<0.0 && i==0) || (gOut.x>0.0 && i==this->ncx+1) ||
                (gOut.y<0.0 && j==0) || (gOut.y>0.0 && j==this->ncy+1) ||
                (gOut.z<0.0 && k==0) || (gOut.z>0.0 && k==this->ncz+1)) {
                //  we are going oustide the grid!
                std::cerr << "Error while computing raypaths: going outside grid!\n"
                << "  Stopping calculations, raypaths will be incomplete.\n" << std::endl;
                return;
            }
            
            iOut = i;
            jOut = boost::math::sign(gOut.y)<0.0 ? j-1 : j;
            kOut = boost::math::sign(gOut.z)<0.0 ? k-1 : k;
            
            // planes we will intersect
            T1 xp = this->xmin + this->dx*(i + boost::math::sign(gOut.x)>0.0 ? 1.0 : 0.0);
            T1 yp = this->ymin + this->dy*(j + boost::math::sign(gOut.y));
            T1 zp = this->zmin + this->dz*(k + boost::math::sign(gOut.z));
            
            if ( fabs(xp-curr_pt.x)<small) {
                xp += this->dx*boost::math::sign(gOut.x);
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
            this->getIJK(curr_pt, i, j, k);
            std::vector<sijk<T2>> voxels;
            
            // find voxels touching node
            if ( i<=this->ncx && j<=this->ncy && k<=this->ncz )
                voxels.push_back( {i,j,k} );
            if ( i<=this->ncx && j<=this->ncy && k>0 )
                voxels.push_back( {i,j,k-1} );
            if ( i>0 && j<=this->ncy && k<=this->ncz )
                voxels.push_back( {i-1,j,k} );
            if ( i>0 && j<=this->ncy && k>0 )
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
            
            if ((gOut.x<0.0 && i==0) || (gOut.x>0.0 && i==this->ncx+1) ||
                (gOut.y<0.0 && j==0) || (gOut.y>0.0 && j==this->ncy+1) ||
                (gOut.z<0.0 && k==0) || (gOut.z>0.0 && k==this->ncz+1)) {
                //  we are going oustide the grid!
                std::cerr << "Error while computing raypaths: going outside grid!\n"
                << "  Stopping calculations, raypaths will be incomplete.\n" << std::endl;
                return;
            }

            iOut = boost::math::sign(gOut.x)<0.0 ? i-1 : i;
            jOut = j;
            kOut = boost::math::sign(gOut.z)<0.0 ? k-1 : k;
            
            // planes we will intersect
            T1 xp = this->xmin + this->dx*(i + boost::math::sign(gOut.x));
            T1 yp = this->ymin + this->dy*(j + boost::math::sign(gOut.y)>0.0 ? 1.0 : 0.0);
            T1 zp = this->zmin + this->dz*(k + boost::math::sign(gOut.z));
            
            if ( fabs(yp-curr_pt.y)<small) {
                yp += this->dy*boost::math::sign(gOut.y);
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
            this->getIJK(curr_pt, i, j, k);
            std::vector<sijk<T2>> voxels;
            
            // find voxels touching node
            if ( i<=this->ncx && j<=this->ncy && k<=this->ncz )
                voxels.push_back( {i,j,k} );
            if ( i<=this->ncx && j>0 && k<=this->ncz )
                voxels.push_back( {i,j-1,k} );
            if ( i>0 && j<=this->ncy && k<=this->ncz )
                voxels.push_back( {i-1,j,k} );
            if ( i>0 && j>0 && k<=this->ncz )
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
            
            if ((gOut.x<0.0 && i==0) || (gOut.x>0.0 && i==this->ncx+1) ||
                (gOut.y<0.0 && j==0) || (gOut.y>0.0 && j==this->ncy+1) ||
                (gOut.z<0.0 && k==0) || (gOut.z>0.0 && k==this->ncz+1)) {
                //  we are going oustide the grid!
                std::cerr << "Error while computing raypaths: going outside grid!\n"
                << "  Stopping calculations, raypaths will be incomplete.\n" << std::endl;
                return;
            }
            
            iOut = boost::math::sign(gOut.x)<0.0 ? i-1 : i;
            jOut = boost::math::sign(gOut.y)<0.0 ? j-1 : j;
            kOut = k;
            
            // planes we will intersect
            T1 xp = this->xmin + this->dx*(i + boost::math::sign(gOut.x));
            T1 yp = this->ymin + this->dy*(j + boost::math::sign(gOut.y));
            T1 zp = this->zmin + this->dz*(k + boost::math::sign(gOut.z)>0.0 ? 1.0 : 0.0);
            
            if ( fabs(zp-curr_pt.z)<small) {
                zp += this->dz*boost::math::sign(gOut.z);
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
            
            if ((gOut.x<0.0 && iOut==0) || (gOut.x>0.0 && iOut==this->ncx+1) ||
                (gOut.y<0.0 && jOut==0) || (gOut.y>0.0 && jOut==this->ncy+1) ||
                (gOut.z<0.0 && kOut==0) || (gOut.z>0.0 && kOut==this->ncz+1)) {
                //  we are going oustide the grid!
                std::cerr << "Error while computing raypaths: going outside grid!\n"
                << "  Stopping calculations, raypaths will be incomplete.\n" << std::endl;
                return;
            }
            
            // planes we will intersect
            T1 xp = this->xmin + this->dx*(iIn + boost::math::sign(gOut.x)>0.0 ? 1.0 : 0.0);
            T1 yp = this->ymin + this->dy*(jIn + boost::math::sign(gOut.y)>0.0 ? 1.0 : 0.0);
            T1 zp = this->zmin + this->dz*(kIn + boost::math::sign(gOut.z)>0.0 ? 1.0 : 0.0);
            
            if ( fabs(xp-curr_pt.x)<small) {
                xp += this->dx*boost::math::sign(gOut.x);
                iOut += boost::math::sign(gOut.x);
            }
            if ( fabs(yp-curr_pt.y)<small) {
                yp += this->dy*boost::math::sign(gOut.y);
                jOut += boost::math::sign(gOut.y);
            }
            if ( fabs(zp-curr_pt.z)<small) {
                zp += this->dz*boost::math::sign(gOut.z);
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
    
    T2 cellNo = this->getCellNo( Rx );
    
    //We calculate the Slowness at the point
    std::vector<T2> list;
    
    for (size_t n3=0; n3 < this->neighbors[ cellNo ].size(); n3++){
        if ( this->nodes[this->neighbors[ cellNo ][n3] ].getPrimary() == 5 ){
            list.push_back(this->neighbors[ cellNo ][n3]);
        }
    }
    
    if ( inverseDistance ) {
        
        std::vector<size_t>::iterator it;
        
        std::vector<Node3Disp<T1, T2>*> interpNodes;
        
        for ( size_t nn=0; nn<list.size(); ++nn )
            interpNodes.push_back( &(this->nodes[list[nn] ]) );
        
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
        
        T1 x[3] = { Rx.x, this->nodes[list[0]].getX(), this->nodes[list[1]].getX() };
        T1 y[3] = { Rx.y, this->nodes[list[0]].getY(), this->nodes[list[2]].getY() };
        T1 z[3] = { Rx.z, this->nodes[list[0]].getZ(), this->nodes[list[4]].getZ() };
        
        T1 s[8] = { this->nodes[list[0]].getNodeSlowness(),
            this->nodes[list[4]].getNodeSlowness(),
            this->nodes[list[2]].getNodeSlowness(),
            this->nodes[list[6]].getNodeSlowness(),
            this->nodes[list[1]].getNodeSlowness(),
            this->nodes[list[5]].getNodeSlowness(),
            this->nodes[list[3]].getNodeSlowness(),
            this->nodes[list[7]].getNodeSlowness() };
        
        return Interpolator<T1>::trilinear(x, y, z, s);
    }
}

#endif

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

#ifndef ttcr_Grid3Drc_h
#define ttcr_Grid3Drc_h

#include <algorithm>
#include <cstring>
#include <exception>
#include <iostream>
#include <fstream>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <ctime>

#ifdef VTK
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkXMLRectilinearGridWriter.h"
#endif

#include <boost/math/special_functions/sign.hpp>

#include "Grid3D.h"

namespace ttcr {

    template<typename T1, typename T2, typename NODE, typename CELL>
    class Grid3Drc : public Grid3D<T1,T2> {
    public:

        Grid3Drc(const T2 nx, const T2 ny, const T2 nz,
                 const T1 ddx, const T1 ddy, const T1 ddz,
                 const T1 minx, const T1 miny, const T1 minz,
                 const bool ttrp, const size_t nt=1,
                 const bool _translateOrigin=false) :
        Grid3D<T1,T2>(ttrp, nx*ny*nz, nt, _translateOrigin),
        dx(ddx), dy(ddy), dz(ddz),
        xmin(minx), ymin(miny), zmin(minz),
        xmax(minx+nx*ddx), ymax(miny+ny*ddy), zmax(minz+nz*ddz),
        ncx(nx), ncy(ny), ncz(nz),
        nodes(std::vector<NODE>((nx+1)*(ny+1)*(nz+1), NODE(nt))),
        cells(CELL(nx*ny*nz))
        { }

        virtual ~Grid3Drc() {}

        void getSlowness(std::vector<T1>& slowness) const {
            if (slowness.size() != ncx*ncy*ncz) {
                slowness.resize(ncx*ncy*ncz);
            }
            for (size_t n=0; n<slowness.size(); ++n) {
                slowness[n] = cells.getSlowness(n);
            }
        }
        void setSlowness(const std::vector<T1>& s) {
            try {
                cells.setSlowness( s );
            } catch (std::exception& e) {
                throw;
            }
        }
        void setChi(const std::vector<T1>& x) {
            cells.setChi( x );
        }
        void setPsi(const std::vector<T1>& x) {
            cells.setPsi( x );
        }

        size_t getNumberOfNodes() const { return nodes.size(); }
        size_t getNumberOfCells() const { return ncx*ncy*ncz; }

        void getTT(std::vector<T1>& tt, const size_t threadNo=0) const final {
            size_t nPrimary = (ncx+1) * (ncy+1) * (ncz+1);
            tt.resize(nPrimary);
            for ( size_t n=0; n<nPrimary; ++n ) {
                tt[n] = nodes[n].getTT(threadNo);
            }
        }

        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const sxyz<T1> &Rx,
                        std::vector<sxyz<T1>> &r_data,
                        T1 &tt,
                        const size_t threadNo) const;

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
                    const int format=1) const;

        size_t getNeighborsSize() const {
            size_t n_elem = 0;
            for ( size_t n=0; n<this->neighbors.size(); ++n ) {
                n_elem += this->neighbors[n].size();
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

        const T1 getXmin() const { return xmin; }
        const T1 getYmin() const { return ymin; }
        const T1 getZmin() const { return zmin; }
        const T1 getDx() const { return dx; }
        const T1 getDy() const { return dy; }
        const T1 getDz() const { return dz; }
        const T2 getNcx() const { return ncx; }
        const T2 getNcy() const { return ncy; }
        const T2 getNcz() const { return ncz; }

        void dump_secondary(std::ofstream& os) const {
            size_t nPrimary = (ncx+1) * (ncy+1) * (ncz+1);
            for ( size_t n=nPrimary; n<nodes.size(); ++n ) {
                os << nodes[n].getX() << ' ' << nodes[n].getY() << ' ' << nodes[n].getZ() << '\n';
            }
        }

        T1 computeSlowness(const sxyz<T1>& pt) const {
            return cells.getSlowness(getCellNo(pt));
        }

    protected:
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

        void buildGridNodes(const T2 nsnx=0, const T2 nsny=0, const T2 nsnz=0);

        T2 getCellNo(const sxyz<T1>& pt) const {
            T1 x = xmax-pt.x < small2 ? xmax-.5*dx : pt.x;
            T1 y = ymax-pt.y < small2 ? ymax-.5*dy : pt.y;
            T1 z = zmax-pt.z < small2 ? zmax-.5*dz : pt.z;
            T2 nx = static_cast<T2>( small2 + (x-xmin)/dx );
            T2 ny = static_cast<T2>( small2 + (y-ymin)/dy );
            T2 nz = static_cast<T2>( small2 + (z-zmin)/dz );
            return ny*ncx + nz*(ncx*ncy) + nx;
        }


        T2 getCellNo(const NODE& node) const {
            T1 x = xmax-node.getX() < small2 ? xmax-.5*dx : node.getX();
            T1 y = ymax-node.getY() < small2 ? ymax-.5*dy : node.getY();
            T1 z = zmax-node.getZ() < small2 ? zmax-.5*dz : node.getZ();
            T2 nx = static_cast<T2>( small2 + (x-xmin)/dx );
            T2 ny = static_cast<T2>( small2 + (y-ymin)/dy );
            T2 nz = static_cast<T2>( small2 + (z-zmin)/dz );
            return ny*ncx + nz*(ncx*ncy) + nx;
        }

        void getCellIJK(const T2 cellNo, sijk<T2> &ind) const {
            ind.k = cellNo / (ncx*ncy);
            ind.j = (cellNo - ind.k*ncx*ncy) / ncx;
            ind.i = cellNo - ncx * ( ind.k*ncy + ind.j);
        }

        void getIJK(const sxyz<T1>& pt, T2& i, T2& j, T2& k) const {
            i = static_cast<T2>( small2 + (pt.x-xmin)/dx );
            j = static_cast<T2>( small2 + (pt.y-ymin)/dy );
            k = static_cast<T2>( small2 + (pt.z-zmin)/dz );
        }

        void getIJK(const sxyz<T1>& pt, long long& i, long long& j, long long& k) const {
            i = static_cast<long long>( small2 + (pt.x-xmin)/dx );
            j = static_cast<long long>( small2 + (pt.y-ymin)/dy );
            k = static_cast<long long>( small2 + (pt.z-zmin)/dz );
        }

        void checkPts(const std::vector<sxyz<T1>>&) const;

        T1 getTraveltime(const sxyz<T1>& Rx,
                         const std::vector<NODE>& nodes,
                         const size_t threadNo) const;

        T1 getTraveltime(const sxyz<T1>& Rx,
                         const std::vector<NODE>& nodes,
                         T2&, T2& , const size_t threadNo) const;

        T1 getTraveltimeFromRaypath(const std::vector<sxyz<T1>>& Tx,
                                    const std::vector<T1>& t0,
                                    const sxyz<T1>& Rx,
                                    const size_t threadNo) const;

        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const sxyz<T1> &Rx,
                        std::vector<siv<T1>> &l_data,
                        T1 &tt,
                        const size_t threadNo) const;

        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const sxyz<T1> &Rx,
                        std::vector<sxyz<T1>> &r_data,
                        std::vector<siv<T1>> &l_data,
                        T1 &tt,
                        const size_t threadNo) const;

    private:
        Grid3Drc() {}
        Grid3Drc(const Grid3Drc<T1,T2,NODE,CELL>& g) {}
        Grid3Drc<T1,T2,NODE,CELL>& operator=(const Grid3Drc<T1,T2,NODE,CELL>& g) { return *this; }

        T1 getTraveltime(const sxyz<T1>& pt,
                         const size_t threadNo) const final;

        void grad(sxyz<T1>& g, const sxyz<T1> &pt,
                  const size_t nt) const;

    };


    template<typename T1, typename T2, typename NODE, typename CELL>
    void Grid3Drc<T1,T2,NODE, CELL>::buildGridNodes(const T2 nsnx,
                                                    const T2 nsny,
                                                    const T2 nsnz) {

        if ( nsnx != 0 || nsny != 0 || nsnz != 0) {
            nodes.resize(// secondary nodes on the edges
                         ncx*nsnx*((ncy+1)*(ncz+1)) +
                         ncy*nsny*((ncx+1)*(ncz+1)) +
                         ncz*nsnz*((ncx+1)*(ncy+1)) +
                         // secondary nodes on the faces
                         (nsnx*nsny)*(ncx*ncy*(ncz+1))+
                         (nsnx*nsnz)*(ncx*ncz*(ncy+1))+
                         (nsny*nsnz)*(ncy*ncz*(ncx+1))+
                         // primary nodes
                         (ncx+1) * (ncy+1) * (ncz+1),
                         NODE(this->nThreads));
        }

        if ( this->translateOrigin ) {
            this->origin = {xmin, ymin, zmin};
            xmin = 0.0;
            ymin = 0.0;
            zmin = 0.0;
        } else {
            this->origin = {0.0, 0.0, 0.0};
        }

        // Create the grid, assign a number for each node and find the owners
        // Nodes and cells are first indexed in z, then y, and x.
        // Secondary nodes are placed on the faces and edges of every cells.
        // Ex: the node in "node[A]=(i,j,k)" is followed by the node in
        // "node[A+1]=(i+dx,j,k)"

        T2 cXmYmZm;     // cell in the (x-,y-,z-) direction from the node
        T2 cXpYmZm;     // cell in the (x+,y-,z-) direction from the node
        T2 cXmYpZm;
        T2 cXpYpZm;
        T2 cXmYmZp;
        T2 cXpYmZp;
        T2 cXmYpZp;
        T2 cXpYpZp;

        T2 n = 0;
        for ( T2 nk=0; nk<=ncz; ++nk ) {

            T1 z = zmin + nk*dz;

            for ( T2 nj=0; nj<=ncy; ++nj ) {

                T1 y = ymin + nj*dy;

                for ( T2 ni=0; ni<=ncx; ++ni, ++n ){

                    T1 x = xmin + ni*dx;

                    // Find the adjacent cells for each primary node

                    if (ni < ncx && nj < ncy && nk < ncz){
                        cXpYpZp = nj*ncx + nk*(ncx*ncy) + ni;
                    }
                    else {
                        cXpYpZp = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj < ncy && nk < ncz){
                        cXmYpZp = nj*ncx + nk*(ncx*ncy) + ni - 1;
                    }
                    else {
                        cXmYpZp = std::numeric_limits<T2>::max();
                    }

                    if (ni < ncx && nj > 0 && nk < ncz){
                        cXpYmZp = (nj-1)*ncx + nk*(ncx*ncy) + ni;
                    }
                    else {
                        cXpYmZp = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj > 0 && nk < ncz){
                        cXmYmZp = (nj-1)*ncx + nk*(ncx * ncy) + ni - 1;
                    }
                    else {
                        cXmYmZp = std::numeric_limits<T2>::max();
                    }

                    if (ni < ncx && nj < ncy && nk > 0){
                        cXpYpZm = nj*ncx + (nk-1)*(ncx*ncy) + ni;
                    }
                    else {
                        cXpYpZm = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj < ncy && nk > 0){
                        cXmYpZm = nj*ncx + (nk-1)*(ncx*ncy) + ni - 1;
                    }
                    else {
                        cXmYpZm = std::numeric_limits<T2>::max();
                    }

                    if (ni < ncx && nj > 0 && nk > 0){
                        cXpYmZm = (nj-1)*ncx + (nk-1)*(ncx*ncy) + ni;
                    }
                    else {
                        cXpYmZm = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj > 0 && nk > 0){
                        cXmYmZm = (nj-1)*ncx + (nk-1)*(ncx*ncy) + ni-1;
                    }
                    else {
                        cXmYmZm = std::numeric_limits<T2>::max();
                    }


                    // Index the primary nodes owners

                    if ( cXmYmZm != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cXmYmZm );
                    }
                    if ( cXpYmZm != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cXpYmZm );
                    }
                    if ( cXmYpZm != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cXmYpZm );
                    }
                    if ( cXpYpZm != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cXpYpZm );
                    }
                    if ( cXmYmZp != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cXmYmZp );
                    }
                    if ( cXpYmZp != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cXpYmZp );
                    }
                    if ( cXmYpZp != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cXmYpZp );
                    }
                    if ( cXpYpZp != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cXpYpZp );
                    }

                    nodes[n].setXYZindex( x, y, z, n );
                    nodes[n].setPrimary(true);
                }
            }
        }

        if ( nsnx != 0 || nsny != 0 || nsnz != 0) {
            T1 dxs = dx/(nsnx+1);     // distance between secondary nodes in x
            T1 dys = dy/(nsny+1);
            T1 dzs = dz/(nsnz+1);

            for ( T2 nk=0; nk<=ncz; ++nk ) {

                T1 z = zmin + nk*dz;

                for ( T2 nj=0; nj<=ncy; ++nj ) {

                    T1 y = ymin + nj*dy;

                    for ( T2 ni=0; ni<=ncx; ++ni ){

                        T1 x = xmin + ni*dx;

                        // Find the adjacent cells for each primary node

                        if (ni < ncx && nj < ncy && nk < ncz){
                            cXpYpZp = nj*ncx + nk*(ncx*ncy) + ni;
                        }
                        else {
                            cXpYpZp = std::numeric_limits<T2>::max();
                        }

                        if (ni > 0 && nj < ncy && nk < ncz){
                            cXmYpZp = nj*ncx + nk*(ncx*ncy) + ni - 1;
                        }
                        else {
                            cXmYpZp = std::numeric_limits<T2>::max();
                        }

                        if (ni < ncx && nj > 0 && nk < ncz){
                            cXpYmZp = (nj-1)*ncx + nk*(ncx*ncy) + ni;
                        }
                        else {
                            cXpYmZp = std::numeric_limits<T2>::max();
                        }

                        if (ni > 0 && nj > 0 && nk < ncz){
                            cXmYmZp = (nj-1)*ncx + nk*(ncx * ncy) + ni - 1;
                        }
                        else {
                            cXmYmZp = std::numeric_limits<T2>::max();
                        }

                        if (ni < ncx && nj < ncy && nk > 0){
                            cXpYpZm = nj*ncx + (nk-1)*(ncx*ncy) + ni;
                        }
                        else {
                            cXpYpZm = std::numeric_limits<T2>::max();
                        }

                        if (ni > 0 && nj < ncy && nk > 0){
                            cXmYpZm = nj*ncx + (nk-1)*(ncx*ncy) + ni - 1;
                        }
                        else {
                            cXmYpZm = std::numeric_limits<T2>::max();
                        }

                        if (ni < ncx && nj > 0 && nk > 0){
                            cXpYmZm = (nj-1)*ncx + (nk-1)*(ncx*ncy) + ni;
                        }
                        else {
                            cXpYmZm = std::numeric_limits<T2>::max();
                        }

                        if (ni > 0 && nj > 0 && nk > 0){
                            cXmYmZm = (nj-1)*ncx + (nk-1)*(ncx*ncy) + ni-1;
                        }
                        else {
                            cXmYmZm = std::numeric_limits<T2>::max();
                        }

                        // Secondary nodes on x edge
                        if ( ni < ncx ) {
                            for (T2 ns=0; ns< nsnx; ++ns, ++n ) {

                                T1 xsv = xmin + ni* dx + (ns+1)*dxs;

                                if ( cXpYmZm != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXpYmZm );
                                }
                                if ( cXpYpZm != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXpYpZm );
                                }
                                if ( cXpYmZp != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXpYmZp );
                                }
                                if ( cXpYpZp != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXpYpZp );
                                }
                                nodes[n].setXYZindex( xsv, y, z, n );
                            }
                        }

                        // Secondary nodes on y edge
                        if ( nj < ncy ) {
                            for (T2 ns=0; ns< nsny; ++ns, ++n ) {

                                T1 ysv = ymin + nj* dy + (ns+1)*dys;

                                if ( cXmYpZm != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXmYpZm );
                                }
                                if ( cXpYpZm != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXpYpZm );
                                }
                                if ( cXmYpZp != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXmYpZp );
                                }
                                if ( cXpYpZp != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXpYpZp );
                                }
                                nodes[n].setXYZindex( x, ysv, z, n );
                            }
                        }

                        // Secondary nodes on z edge
                        if ( nk < ncz ) {
                            for (T2 ns=0; ns< nsnz; ++ns, ++n ) {

                                T1 zsv = zmin + nk* dz + (ns+1)*dzs;

                                if ( cXmYmZp != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXmYmZp );
                                }
                                if ( cXpYmZp != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXpYmZp );
                                }
                                if ( cXmYpZp != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXmYpZp );
                                }
                                if ( cXpYpZp != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXpYpZp );
                                }
                                nodes[n].setXYZindex( x, y, zsv, n );
                            }
                        }

                        // Secondary nodes on the xy0 planes
                        if ( ni < ncx && nj < ncy ) {
                            for ( T2 sy=0; sy < nsny; ++sy ) {
                                for ( T2 sx=0; sx < nsnx; ++sx, n++ ) {

                                    T1 ysv = ymin+ nj* dy+ (sy+1)*dys;
                                    T1 xsv = xmin+ ni* dx+ (sx+1)*dxs;

                                    if ( cXpYpZm != std::numeric_limits<T2>::max() )
                                    {
                                        nodes[n].pushOwner( cXpYpZm );
                                    }
                                    if ( cXpYpZp != std::numeric_limits<T2>::max() )
                                    {
                                        nodes[n].pushOwner( cXpYpZp );
                                    }
                                    nodes[n].setXYZindex( xsv, ysv, z, n );
                                }
                            }
                        }

                        // Secondary nodes on the x0z planes
                        if ( ni < ncx && nk < ncz ) {
                            for ( T2 sz=0; sz < nsnz; ++sz ) {
                                for ( T2 sx=0; sx < nsnx; ++sx, n++ ) {

                                    T1 zsv = zmin+ nk* dz+ (sz+1)*dzs;
                                    T1 xsv = xmin+ ni* dx+ (sx+1)*dxs;

                                    if ( cXpYmZp != std::numeric_limits<T2>::max() )
                                    {
                                        nodes[n].pushOwner( cXpYmZp );
                                    }
                                    if ( cXpYpZp != std::numeric_limits<T2>::max() )
                                    {
                                        nodes[n].pushOwner( cXpYpZp );
                                    }
                                    nodes[n].setXYZindex( xsv, y, zsv, n );
                                }
                            }
                        }

                        // Secondary nodes on the 0yz planes
                        if ( nj < ncy && nk < ncz ) {
                            for ( T2 sz=0; sz < nsnz; ++sz ) {
                                for ( T2 sy=0; sy < nsny; ++sy, n++ ) {

                                    T1 zsv = zmin+ nk* dz+ (sz+1)*dzs;
                                    T1 ysv = ymin+ nj* dy+ (sy+1)*dys;

                                    if ( cXmYpZp != std::numeric_limits<T2>::max() )
                                    {
                                        nodes[n].pushOwner( cXmYpZp );
                                    }
                                    if ( cXpYpZp != std::numeric_limits<T2>::max() )
                                    {
                                        nodes[n].pushOwner( cXpYpZp );
                                    }
                                    nodes[n].setXYZindex( x, ysv, zsv, n );
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


    template<typename T1, typename T2, typename NODE, typename CELL>
    void Grid3Drc<T1,T2,NODE,CELL>::checkPts(const std::vector<sxyz<T1>>& pts) const {

        // Check if the points from a vector are in the grid
        for ( size_t n=0; n<pts.size(); ++n ) {
            if ( pts[n].x < xmin || pts[n].x > xmax ||
                pts[n].y < ymin || pts[n].y > ymax ||
                pts[n].z < zmin || pts[n].z > zmax ) {
                std::ostringstream msg;
                msg << "Error: Point (" << pts[n].x << ", " << pts[n].y << ", " << pts[n] .z << ") outside grid.";
                throw std::runtime_error(msg.str());
            }
        }
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
        size_t neibNo = this->neighbors[cellNo][0];
        T1 dt = cells.computeDt(nodes[neibNo], Rx, cellNo);

        T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
        for ( size_t k=1; k< this->neighbors[cellNo].size(); ++k ) {
            neibNo = this->neighbors[cellNo][k];
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
        T2 neibNo = this->neighbors[cellNo][0];
        T1 dt = cells.computeDt(nodes[neibNo], Rx, cellNo);

        T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
        nodeParentRx = neibNo;
        cellParentRx = cellNo;
        for ( size_t k=1; k< this->neighbors[cellNo].size(); ++k ) {
            neibNo = this->neighbors[cellNo][k];
            dt = cells.computeDt(nodes[neibNo], Rx, cellNo);
            if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
                traveltime =  nodes[neibNo].getTT(threadNo)+dt;
                nodeParentRx = neibNo;
            }
        }
        return traveltime;
    }


    template<typename T1, typename T2, typename NODE, typename CELL>
    T1 Grid3Drc<T1,T2,NODE,CELL>::getTraveltime(const sxyz<T1> &pt,
                                                const size_t nt) const {

        const size_t nnx = ncx+1;
        const size_t nny = ncy+1;

        // trilinear interpolation if not on node

        T1 tt;
        T2 i, j, k;

        getIJK(pt, i, j, k);

        if ( std::abs(pt.x - (xmin+i*dx))<small2 &&
            std::abs(pt.y - (ymin+j*dy))<small2 &&
            std::abs(pt.z - (zmin+k*dz))<small2 ) {
            // on node
            return nodes[(k*nny+j)*nnx+i].getTT(nt);
        } else if ( std::abs(pt.x - (xmin+i*dx))<small2 &&
                   std::abs(pt.y - (ymin+j*dy))<small2 ) {
            // on edge
            T1 t1 = nodes[(    k*nny+j)*nnx+i].getTT(nt);
            T1 t2 = nodes[((k+1)*nny+j)*nnx+i].getTT(nt);

            T1 w1 = (zmin+(k+1)*dz - pt.z)/dz;
            T1 w2 = (pt.z - (zmin+k*dz))/dz;

            tt = t1*w1 + t2*w2;

        } else if ( std::abs(pt.x - (xmin+i*dx))<small2 &&
                   std::abs(pt.z - (zmin+k*dz))<small2 ) {
            // on edge
            T1 t1 = nodes[(k*nny+j  )*nnx+i].getTT(nt);
            T1 t2 = nodes[(k*nny+j+1)*nnx+i].getTT(nt);

            T1 w1 = (ymin+(j+1)*dy - pt.y)/dy;
            T1 w2 = (pt.y - (ymin+j*dy))/dy;

            tt = t1*w1 + t2*w2;

        } else if ( std::abs(pt.y - (ymin+j*dy))<small2 &&
                   std::abs(pt.z - (zmin+k*dz))<small2 ) {
            // on edge
            T1 t1 = nodes[(k*nny+j)*nnx+i  ].getTT(nt);
            T1 t2 = nodes[(k*nny+j)*nnx+i+1].getTT(nt);

            T1 w1 = (xmin+(i+1)*dx - pt.x)/dx;
            T1 w2 = (pt.x - (xmin+i*dx))/dx;

            tt = t1*w1 + t2*w2;

        } else if ( std::abs(pt.x - (xmin+i*dx))<small2 ) {
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

        } else if ( std::abs(pt.y - (ymin+j*dy))<small2 ) {
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

        } else if ( std::abs(pt.z - (zmin+k*dz))<small2 ) {
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

    template<typename T1, typename T2, typename NODE, typename CELL>
    void Grid3Drc<T1,T2,NODE,CELL>::grad(sxyz<T1>& g, const sxyz<T1> &pt,
                                         const size_t nt) const {

        // compute travel time gradient at point pt

        T1 p1 = pt.x - dx/2.0;
        T1 p2 = p1 + dx;
        if ( p1 < xmin ) {  // check if on grid edge or out of grid
            p1 = xmin;  // shift pt to allow interpolating in getTraveltime
            p2 = p1 + dx;
        } else if ( p2 > xmax ) {
            p2 = xmax;
            p1 = p2 - dx;
        }
        g.x = (getTraveltime({p2, pt.y, pt.z}, nt) - getTraveltime({p1, pt.y, pt.z}, nt)) / dx;

        p1 = pt.y - dy/2.0;
        p2 = p1 + dy;
        if ( p1 < ymin ) {
            p1 = ymin;
            p2 = p1 + dy;
        } else if ( p2 > ymax ) {
            p2 = ymax;
            p1 = p2 - dy;
        }
        g.y = (getTraveltime({pt.x, p2, pt.z}, nt) - getTraveltime({pt.x, p1, pt.z}, nt)) / dy;

        p1 = pt.z - dz/2.0;
        p2 = p1 + dz;
        if ( p1 < zmin ) {
            p1 = zmin;
            p2 = p1 + dz;
        } else if ( p2 > zmax ) {
            p2 = zmax;
            p1 = p2 - dz;
        }
        g.z = (getTraveltime({pt.x, pt.y, p2}, nt) - getTraveltime({pt.x, pt.y, p1}, nt)) / dz;

    }

    template<typename T1, typename T2, typename NODE, typename CELL>
    T1 Grid3Drc<T1,T2,NODE,CELL>::getTraveltimeFromRaypath(const std::vector<sxyz<T1>>& Tx,
                                                           const std::vector<T1>& t0,
                                                           const sxyz<T1> &Rx,
                                                           const size_t threadNo) const {
        T1 tt = 0.0;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return tt;
            }
        }

        sxyz<T1> prev_pt( Rx );
        sxyz<T1> curr_pt( Rx );
        // distance between opposite nodes of a voxel
        const T1 maxDist = sqrt( dx*dx + dy*dy + dz*dz );
        sxyz<T1> g;

        bool reachedTx = false;
        while ( reachedTx == false ) {

            grad(g, curr_pt, threadNo);
            g *= -1.0;

            long long i, j, k;
            getIJK(curr_pt, i, j, k);

            // planes we will intersect
            T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
            T1 yp = ymin + dy*(j + (boost::math::sign(g.y)>0.0 ? 1.0 : 0.0));
            T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

            if ( std::abs(xp-curr_pt.x)<small2) {
                xp += dx*boost::math::sign(g.x);
            }
            if ( std::abs(yp-curr_pt.y)<small2) {
                yp += dy*boost::math::sign(g.y);
            }
            if ( std::abs(zp-curr_pt.z)<small2) {
                zp += dz*boost::math::sign(g.z);
            }

            // dist to planes
            T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
            T1 ty = g.y!=0.0 ? (yp - curr_pt.y)/g.y : std::numeric_limits<T1>::max();
            T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

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
                std::ostringstream msg;
                msg << "Error while computing raypaths: going outside grid \n\
                Rx: " << Rx << "\n\
                Tx: " << Tx[0] << "\n";
                for ( size_t ns=1; ns<Tx.size(); ++ns ) {
                    msg << "\
                    " << Tx[ns] << "\n";
                }
                throw std::runtime_error(msg.str());
            }
            sxyz<T1> mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
            T2 cellNo = getCellNo(mid_pt);
            tt += cells.computeDt(prev_pt, curr_pt, cellNo);
            prev_pt = curr_pt;

            // are we close enough to one the Tx nodes ?
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                T1 dist = curr_pt.getDistance( Tx[ns] );
                if ( dist < maxDist ) {

                    g = Tx[ns] - curr_pt;
                    // check if we intersect a plane between curr_pt & Tx

                    getIJK(curr_pt, i, j, k);

                    xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
                    yp = ymin + dy*(j + (boost::math::sign(g.y)>0.0 ? 1.0 : 0.0));
                    zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

                    if ( std::abs(xp-curr_pt.x)<small2) {
                        xp += dx*boost::math::sign(g.x);
                    }
                    if ( std::abs(yp-curr_pt.y)<small2) {
                        yp += dy*boost::math::sign(g.y);
                    }
                    if ( std::abs(zp-curr_pt.z)<small2) {
                        zp += dz*boost::math::sign(g.z);
                    }

                    // dist to planes
                    tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
                    ty = g.y!=0.0 ? (yp - curr_pt.y)/g.y : std::numeric_limits<T1>::max();
                    tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

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

                    if ( curr_pt.getDistance(prev_pt) > dist ||  // we do not intersect a plane
                        curr_pt == Tx[ns] ) {  // we have arrived
                        cellNo = getCellNo(Tx[ns]);
                        tt += cells.computeDt(Tx[ns], prev_pt, cellNo);
                    } else {
                        // to intersection
                        mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
                        cellNo = getCellNo(mid_pt);
                        tt += cells.computeDt(prev_pt, curr_pt, cellNo);
                        // to Tx
                        cellNo = getCellNo(Tx[ns]);
                        tt += cells.computeDt(Tx[ns], curr_pt, cellNo);
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
        return tt;
    }

    template<typename T1, typename T2, typename NODE, typename CELL>
    void Grid3Drc<T1,T2,NODE,CELL>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                               const std::vector<T1>& t0,
                                               const sxyz<T1> &Rx,
                                               std::vector<sxyz<T1>> &r_data,
                                               T1 &tt,
                                               const size_t threadNo) const {
        tt = 0.0;
        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        sxyz<T1> curr_pt( Rx );
        // distance between opposite nodes of a voxel
        const T1 maxDist = sqrt( dx*dx + dy*dy + dz*dz );
        sxyz<T1> g;
#ifdef DEBUG_RP
        std::cout << "Starting at " << curr_pt << '\n';
#endif

        bool reachedTx = false;
        while ( reachedTx == false ) {

            grad(g, curr_pt, threadNo);
            g *= -1.0;

            long long i, j, k;
            getIJK(curr_pt, i, j, k);

            // planes we will intersect
            T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
            T1 yp = ymin + dy*(j + (boost::math::sign(g.y)>0.0 ? 1.0 : 0.0));
            T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

            if ( std::abs(xp-curr_pt.x)<small2) {
                xp += dx*boost::math::sign(g.x);
            }
            if ( std::abs(yp-curr_pt.y)<small2) {
                yp += dy*boost::math::sign(g.y);
            }
            if ( std::abs(zp-curr_pt.z)<small2) {
                zp += dz*boost::math::sign(g.z);
            }

            // dist to planes
            T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
            T1 ty = g.y!=0.0 ? (yp - curr_pt.y)/g.y : std::numeric_limits<T1>::max();
            T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

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
#ifdef DEBUG_RP
            std::cout << "Grad: " << g << "\t going to: " << curr_pt << '\n';
#endif

            if ( curr_pt.x < xmin || curr_pt.x > xmax ||
                curr_pt.y < ymin || curr_pt.y > ymax ||
                curr_pt.z < zmin || curr_pt.z > zmax ) {
                //  we are going oustide the grid!
                std::ostringstream msg;
                msg << "Error while computing raypaths: going outside grid \n\
                Rx: " << Rx << "\n\
                Tx: " << Tx[0] << "\n";
                for ( size_t ns=1; ns<Tx.size(); ++ns ) {
                    msg << "\
                    " << Tx[ns] << "\n";
                }
                throw std::runtime_error(msg.str());
            }
            sxyz<T1> mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
            T2 cellNo = getCellNo(mid_pt);
            tt += cells.computeDt(r_data.back(), curr_pt, cellNo);
            r_data.push_back( curr_pt );

            // are we close enough to one the Tx nodes ?
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                T1 dist = curr_pt.getDistance( Tx[ns] );
                if ( dist < maxDist ) {

                    g = Tx[ns] - curr_pt;
                    // check if we intersect a plane between curr_pt & Tx

                    getIJK(curr_pt, i, j, k);

                    xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
                    yp = ymin + dy*(j + (boost::math::sign(g.y)>0.0 ? 1.0 : 0.0));
                    zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

                    if ( std::abs(xp-curr_pt.x)<small2) {
                        xp += dx*boost::math::sign(g.x);
                    }
                    if ( std::abs(yp-curr_pt.y)<small2) {
                        yp += dy*boost::math::sign(g.y);
                    }
                    if ( std::abs(zp-curr_pt.z)<small2) {
                        zp += dz*boost::math::sign(g.z);
                    }

                    // dist to planes
                    tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
                    ty = g.y!=0.0 ? (yp - curr_pt.y)/g.y : std::numeric_limits<T1>::max();
                    tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

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

                    if ( curr_pt.getDistance(r_data.back()) > dist ||  // we do not intersect a plane
                        curr_pt == Tx[ns] ) {  // we have arrived
                        cellNo = getCellNo(Tx[ns]);
                        tt += cells.computeDt(Tx[ns], r_data.back(), cellNo);
                        r_data.push_back( Tx[ns] );
                    } else {
                        // to intersection
                        mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
                        cellNo = getCellNo(mid_pt);
                        tt += cells.computeDt(r_data.back(), curr_pt, cellNo);
                        r_data.push_back( curr_pt );
                        // to Tx
                        cellNo = getCellNo(Tx[ns]);
                        tt += cells.computeDt(Tx[ns], curr_pt, cellNo);
                        r_data.push_back( Tx[ns] );
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
    }

    template<typename T1, typename T2, typename NODE, typename CELL>
    void Grid3Drc<T1,T2,NODE,CELL>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                               const std::vector<T1>& t0,
                                               const sxyz<T1> &Rx,
                                               std::vector<siv<T1>> &l_data,
                                               T1 &tt,
                                               const size_t threadNo) const {

        tt = 0.0;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        sxyz<T1> curr_pt( Rx ), prev_pt( Rx );
        // distance between opposite nodes of a voxel
        const T1 maxDist = sqrt( dx*dx + dy*dy + dz*dz );
        sxyz<T1> g;
#ifdef DEBUG_RP
        std::cout << "Starting at " << curr_pt << '\n';
#endif

        siv<T1> cell;
        bool reachedTx = false;
        while ( reachedTx == false ) {

            grad(g, curr_pt, threadNo);
            g *= -1.0;

            long long i, j, k;
            getIJK(curr_pt, i, j, k);

            // planes we will intersect
            T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
            T1 yp = ymin + dy*(j + (boost::math::sign(g.y)>0.0 ? 1.0 : 0.0));
            T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

            if ( std::abs(xp-curr_pt.x)<small2) {
                xp += dx*boost::math::sign(g.x);
            }
            if ( std::abs(yp-curr_pt.y)<small2) {
                yp += dy*boost::math::sign(g.y);
            }
            if ( std::abs(zp-curr_pt.z)<small2) {
                zp += dz*boost::math::sign(g.z);
            }

            // dist to planes
            T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
            T1 ty = g.y!=0.0 ? (yp - curr_pt.y)/g.y : std::numeric_limits<T1>::max();
            T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

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
#ifdef DEBUG_RP
            std::cout << "Grad: " << g << "\t going to: " << curr_pt << '\n';
#endif

            if ( curr_pt.x < xmin || curr_pt.x > xmax ||
                curr_pt.y < ymin || curr_pt.y > ymax ||
                curr_pt.z < zmin || curr_pt.z > zmax ) {
                //  we are going oustide the grid!
                std::ostringstream msg;
                msg << "Error while computing raypaths: going outside grid \n\
                Rx: " << Rx << "\n\
                Tx: " << Tx[0] << "\n";
                for ( size_t ns=1; ns<Tx.size(); ++ns ) {
                    msg << "\
                    " << Tx[ns] << "\n";
                }
                throw std::runtime_error(msg.str());
            }
            sxyz<T1> mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
            cell.i = getCellNo(mid_pt);
            cell.v = curr_pt.getDistance(prev_pt);
            l_data.push_back(cell);
            tt += cells.computeDt(prev_pt, curr_pt, cell.i);
            prev_pt = curr_pt;

            // are we close enough to one the Tx nodes ?
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                T1 dist = curr_pt.getDistance( Tx[ns] );
                if ( dist < maxDist ) {

                    g = Tx[ns] - curr_pt;
                    // check if we intersect a plane between curr_pt & Tx

                    getIJK(curr_pt, i, j, k);

                    xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
                    yp = ymin + dy*(j + (boost::math::sign(g.y)>0.0 ? 1.0 : 0.0));
                    zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

                    if ( std::abs(xp-curr_pt.x)<small2) {
                        xp += dx*boost::math::sign(g.x);
                    }
                    if ( std::abs(yp-curr_pt.y)<small2) {
                        yp += dy*boost::math::sign(g.y);
                    }
                    if ( std::abs(zp-curr_pt.z)<small2) {
                        zp += dz*boost::math::sign(g.z);
                    }

                    // dist to planes
                    tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
                    ty = g.y!=0.0 ? (yp - curr_pt.y)/g.y : std::numeric_limits<T1>::max();
                    tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

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

                    if ( curr_pt.getDistance(prev_pt) > dist ||  // we do not intersect a plane
                        curr_pt == Tx[ns] ) {  // we have arrived
                        cell.i = getCellNo(Tx[ns]);
                        cell.v = Tx[ns].getDistance(prev_pt);
                        l_data.push_back(cell);
                        tt += cells.computeDt(Tx[ns], prev_pt, cell.i);
                    } else {
                        // to intersection
                        mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
                        cell.i = getCellNo(mid_pt);
                        cell.v = curr_pt.getDistance(prev_pt);
                        l_data.push_back(cell);
                        tt += cells.computeDt(prev_pt, curr_pt, cell.i);

                        // to Tx
                        cell.i = getCellNo(Tx[ns]);
                        cell.v = Tx[ns].getDistance(curr_pt);
                        l_data.push_back(cell);
                        tt += cells.computeDt(Tx[ns], curr_pt, cell.i);
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
    }

    template<typename T1, typename T2, typename NODE, typename CELL>
    void Grid3Drc<T1,T2,NODE,CELL>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                               const std::vector<T1>& t0,
                                               const sxyz<T1> &Rx,
                                               std::vector<sxyz<T1>> &r_data,
                                               std::vector<siv<T1>> &l_data,
                                               T1 &tt,
                                               const size_t threadNo) const {

        tt = 0.0;
        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        sxyz<T1> curr_pt( Rx );
        // distance between opposite nodes of a voxel
        const T1 maxDist = sqrt( dx*dx + dy*dy + dz*dz );
        sxyz<T1> g;
#ifdef DEBUG_RP
        std::cout << "Starting at " << curr_pt << '\n';
#endif

        siv<T1> cell;
        bool reachedTx = false;
        while ( reachedTx == false ) {

            grad(g, curr_pt, threadNo);
            g *= -1.0;

            long long i, j, k;
            getIJK(curr_pt, i, j, k);

            // planes we will intersect
            T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
            T1 yp = ymin + dy*(j + (boost::math::sign(g.y)>0.0 ? 1.0 : 0.0));
            T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

            if ( std::abs(xp-curr_pt.x)<small2) {
                xp += dx*boost::math::sign(g.x);
            }
            if ( std::abs(yp-curr_pt.y)<small2) {
                yp += dy*boost::math::sign(g.y);
            }
            if ( std::abs(zp-curr_pt.z)<small2) {
                zp += dz*boost::math::sign(g.z);
            }

            // dist to planes
            T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
            T1 ty = g.y!=0.0 ? (yp - curr_pt.y)/g.y : std::numeric_limits<T1>::max();
            T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

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
#ifdef DEBUG_RP
            std::cout << "Grad: " << g << "\t going to: " << curr_pt << '\n';
#endif

            if ( curr_pt.x < xmin || curr_pt.x > xmax ||
                curr_pt.y < ymin || curr_pt.y > ymax ||
                curr_pt.z < zmin || curr_pt.z > zmax ) {
                //  we are going oustide the grid!
                std::ostringstream msg;
                msg << "Error while computing raypaths: going outside grid \n\
                Rx: " << Rx << "\n\
                Tx: " << Tx[0] << "\n";
                for ( size_t ns=1; ns<Tx.size(); ++ns ) {
                    msg << "\
                    " << Tx[ns] << "\n";
                }
                throw std::runtime_error(msg.str());
            }
            sxyz<T1> mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
            cell.i = getCellNo(mid_pt);
            cell.v = curr_pt.getDistance(r_data.back());
            l_data.push_back(cell);
            tt += cells.computeDt(r_data.back(), curr_pt, cell.i);
            r_data.push_back( curr_pt );

            // are we close enough to one the Tx nodes ?
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                T1 dist = curr_pt.getDistance( Tx[ns] );
                if ( dist < maxDist ) {

                    g = Tx[ns] - curr_pt;
                    // check if we intersect a plane between curr_pt & Tx

                    getIJK(curr_pt, i, j, k);

                    xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
                    yp = ymin + dy*(j + (boost::math::sign(g.y)>0.0 ? 1.0 : 0.0));
                    zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

                    if ( std::abs(xp-curr_pt.x)<small2) {
                        xp += dx*boost::math::sign(g.x);
                    }
                    if ( std::abs(yp-curr_pt.y)<small2) {
                        yp += dy*boost::math::sign(g.y);
                    }
                    if ( std::abs(zp-curr_pt.z)<small2) {
                        zp += dz*boost::math::sign(g.z);
                    }

                    // dist to planes
                    tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
                    ty = g.y!=0.0 ? (yp - curr_pt.y)/g.y : std::numeric_limits<T1>::max();
                    tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

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

                    if ( curr_pt.getDistance(r_data.back()) > dist ||  // we do not intersect a plane
                        curr_pt == Tx[ns] ) {  // we have arrived
                        cell.i = getCellNo(Tx[ns]);
                        cell.v = Tx[ns].getDistance(r_data.back());
                        l_data.push_back(cell);
                        tt += cells.computeDt(Tx[ns], r_data.back(), cell.i);
                        r_data.push_back( Tx[ns] );
                    } else {
                        // to intersection
                        mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
                        cell.i = getCellNo(mid_pt);
                        cell.v = curr_pt.getDistance(r_data.back());
                        l_data.push_back(cell);
                        tt += cells.computeDt(r_data.back(), curr_pt, cell.i);
                        r_data.push_back( curr_pt );
                        // to Tx
                        cell.i = getCellNo(Tx[ns]);
                        cell.v = Tx[ns].getDistance(r_data.back());
                        l_data.push_back(cell);
                        tt += cells.computeDt(Tx[ns], curr_pt, cell.i);
                        r_data.push_back( Tx[ns] );
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
    }



    template<typename T1, typename T2, typename NODE, typename CELL>
    void Grid3Drc<T1,T2,NODE,CELL>::saveTT(const std::string &fname,
                                           const int all,
                                           const size_t nt,
                                           const int format) const {

        if ( format == 1 ) {
            std::string filename = fname+".dat";
            std::ofstream fout(filename.c_str());
            fout.precision(12);
            for ( T2 n=0; n<nodes.size(); ++n ) {
                if ( nodes[n].isPrimary() || all==1 ) {
                    fout << nodes[n].getX() << '\t'
                    << nodes[n].getY() << '\t'
                    << nodes[n].getZ() << '\t'
                    << nodes[n].getTT(nt) << '\n';
                }
            }
            fout.close();
        } else if ( format == 2 ) {
#ifdef VTK

            std::string filename = fname+".vtr";
            int nn[3] = {static_cast<int>(ncx+1), static_cast<int>(ncy+1), static_cast<int>(ncz+1)};

            vtkSmartPointer<vtkDoubleArray> xCoords = vtkSmartPointer<vtkDoubleArray>::New();
            for (size_t n=0; n<nn[0]; ++n) {
                xCoords->InsertNextValue( xmin + n*dx );
            }
            vtkSmartPointer<vtkDoubleArray> yCoords = vtkSmartPointer<vtkDoubleArray>::New();
            for (size_t n=0; n<nn[1]; ++n) {
                yCoords->InsertNextValue( ymin + n*dy );
            }
            vtkSmartPointer<vtkDoubleArray> zCoords = vtkSmartPointer<vtkDoubleArray>::New();
            for (size_t n=0; n<nn[2]; ++n) {
                zCoords->InsertNextValue( zmin + n*dz );
            }

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
                if ( nodes[n].isPrimary() ) {
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
        } else if ( format == 3 ) {
            std::string filename = fname+".bin";
            std::ofstream fout(filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
            for ( T2 n=0; n<nodes.size(); ++n ) {
                if ( nodes[n].isPrimary() || all==1 ) {
                    T1 tmp[] = { nodes[n].getX(), nodes[n].getY(), nodes[n].getZ(), nodes[n].getTT(nt) };
                    fout.write( (char*)tmp, 4*sizeof(T1) );
                }
            }
            fout.close();
        } else {
            throw std::runtime_error("Unsupported format for saving traveltimes");
        }
    }

}

#endif

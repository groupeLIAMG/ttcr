//
//  Grid3Drn.h
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

#ifndef ttcr_Grid3Drn_h
#define ttcr_Grid3Drn_h

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <queue>
#include <sstream>
#include <stdexcept>
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

namespace ttcr {

    template<typename T1, typename T2, typename NODE>
    class Grid3Drn : public Grid3D<T1,T2> {
    public:

        /* Constructor Format:
         Grid3Drn<T1,T2>::Grid3Drn(nb cells in x, nb cells in y, nb cells in z,
         x cells size, y cells size, z cells size,
         x origin, y origin, z origin,
         index of the thread)
         */
        Grid3Drn(const T2 nx, const T2 ny, const T2 nz,
                 const T1 ddx, const T1 ddy, const T1 ddz,
                 const T1 minx, const T1 miny, const T1 minz,
                 const bool ttrp, const bool procVel, const size_t nt=1,
                 const bool _translateOrigin=false) :
        Grid3D<T1,T2>(ttrp, nx*ny*nz, nt, _translateOrigin),
        dx(ddx), dy(ddy), dz(ddz),
        xmin(minx), ymin(miny), zmin(minz),
        xmax(minx+nx*ddx), ymax(miny+ny*ddy), zmax(minz+nz*ddz),
        ncx(nx), ncy(ny), ncz(nz), processVel(procVel),
        nodes(std::vector<NODE>((nx+1)*(ny+1)*(nz+1), NODE(nt)))
        { }

        virtual ~Grid3Drn() {}

        virtual void setSlowness(const std::vector<T1>& s) {
            if ( nodes.size() != s.size() ) {
                throw std::length_error("Error: slowness vectors of incompatible size.");
            }
            for ( size_t n=0; n<nodes.size(); ++n ) {
                nodes[n].setNodeSlowness( s[n] );
            }
        }
        void getSlowness(std::vector<T1>& slowness) const {
            if (slowness.size() != (ncx+1) * (ncy+1) * (ncz+1)) {
                slowness.resize((ncx+1) * (ncy+1) * (ncz+1));
            }
            for ( size_t n=0; n<slowness.size(); ++n ) {
                slowness[n] = nodes[n].getNodeSlowness();
            }
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

        void saveSlownessXYZ(const char filename[]) const {
            //Saves the Slowness of the primary nodes
            std::ofstream fout( filename );
            for ( size_t n=0; n< nodes.size(); ++n ) {
                if (nodes[n].isPrimary() ){
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

        void saveTT(const std::string &, const int, const size_t nt=0,
                    const int format=1) const;
        void loadTT(const std::string &, const int, const size_t nt=0,
                    const int format=1) const;

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

        T1 computeSlowness(const sxyz<T1>&) const;

#ifdef VTK
        void saveModelVTR(const std::string &,
                          const bool saveSlowness=true) const;
#endif

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

        bool processVel;

        mutable std::vector<NODE> nodes;

        void buildGridNodes(const T2 nsnx=0, const T2 nsny=0, const T2 nsnz=0);

        void interpSecondary();

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

        void gradO2(sxyz<T1>& g, const sxyz<T1> &pt, const size_t nt) const;
        void grad(sxyz<T1>& g, const sxyz<T1> &pt, const size_t nt) const;

        T1 getTraveltimeFromRaypath(const std::vector<sxyz<T1>>& Tx,
                                    const std::vector<T1>& t0,
                                    const sxyz<T1> &Rx,
                                    const size_t threadNo) const;

        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const sxyz<T1> &Rx,
                        std::vector<sxyz<T1>> &r_data,
                        const size_t threadNo) const;

        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const sxyz<T1> &Rx,
                        std::vector<sxyz<T1>> &r_data,
                        T1 &tt,
                        const size_t threadNo) const;

        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const sxyz<T1>& Rx,
                        std::vector<sijv<T1>>& m_data,
                        T1 &tt,
                        const size_t RxNo,
                        const size_t threadNo) const;

        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const sxyz<T1>& Rx,
                        std::vector<sxyz<T1>>& r_data,
                        std::vector<sijv<T1>>& m_data,
                        T1 &tt,
                        const size_t RxNo,
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

        void getRaypath_old(const std::vector<sxyz<T1>>& Tx,
                            const sxyz<T1> &Rx,
                            std::vector<sxyz<T1>> &r_data,
                            const size_t threadNo=0) const;

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
        Grid3Drn() {}
        Grid3Drn(const Grid3Drn<T1,T2,NODE>& g) {}
        Grid3Drn<T1,T2,NODE>& operator=(const Grid3Drn<T1,T2,NODE>& g) { return *this; }

    };


    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::buildGridNodes(const T2 nsnx, const T2 nsny, const T2 nsnz) {

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

        // Create the grid, assign a number for each node, determine the type of the node and find the owners
        // Nodes and cells are first indexed in z, then y, and x.
        // Secondary nodes are placed on the faces and edges of every cells.
        // Ex: the node in "node[A]=(i,j,k)" is followed by the node in "node[A+1]=(i+dx,j,k)"

        T2 cell_XmYmZm;     // cell in the (x-,y-,z-) direction from the node
        T2 cell_XpYmZm;     // cell in the (x+,y-,z-) direction from the node
        T2 cell_XmYpZm;
        T2 cell_XpYpZm;
        T2 cell_XmYmZp;
        T2 cell_XpYmZp;
        T2 cell_XmYpZp;
        T2 cell_XpYpZp;

        T2 n=0;
        for ( T2 nk=0; nk<=ncz; ++nk ) {

            T1 z = zmin + nk*dz;

            for ( T2 nj=0; nj<=ncy; ++nj ) {

                T1 y = ymin + nj*dy;

                for ( T2 ni=0; ni<=ncx; ++ni, ++n ){

                    T1 x = xmin + ni*dx;

                    // Find the adjacent cells for each primary node

                    if (ni < ncx && nj < ncy && nk < ncz){
                        cell_XpYpZp = nj*ncx + nk*(ncx*ncy) + ni;
                    }
                    else {
                        cell_XpYpZp = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj < ncy && nk < ncz){
                        cell_XmYpZp = nj*ncx + nk*(ncx*ncy) + ni - 1;
                    }
                    else {
                        cell_XmYpZp = std::numeric_limits<T2>::max();
                    }

                    if (ni < ncx && nj > 0 && nk < ncz){
                        cell_XpYmZp = (nj-1)*ncx + nk*(ncx*ncy) + ni;
                    }
                    else {
                        cell_XpYmZp = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj > 0 && nk < ncz){
                        cell_XmYmZp = (nj-1)*ncx + nk*(ncx * ncy) + ni - 1;
                    }
                    else {
                        cell_XmYmZp = std::numeric_limits<T2>::max();
                    }

                    if (ni < ncx && nj < ncy && nk > 0){
                        cell_XpYpZm = nj*ncx + (nk-1)*(ncx*ncy) + ni;
                    }
                    else {
                        cell_XpYpZm = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj < ncy && nk > 0){
                        cell_XmYpZm = nj*ncx + (nk-1)*(ncx*ncy) + ni - 1;
                    }
                    else {
                        cell_XmYpZm = std::numeric_limits<T2>::max();
                    }

                    if (ni < ncx && nj > 0 && nk > 0){
                        cell_XpYmZm = (nj-1)*ncx + (nk-1)*(ncx*ncy) + ni;
                    }
                    else {
                        cell_XpYmZm = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj > 0 && nk > 0){
                        cell_XmYmZm = (nj-1)*ncx + (nk-1)*(ncx*ncy) + ni - 1;
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
                            cell_XpYpZp = nj*ncx + nk*(ncx*ncy) + ni;
                        }
                        else {
                            cell_XpYpZp = std::numeric_limits<T2>::max();
                        }

                        if (ni > 0 && nj < ncy && nk < ncz){
                            cell_XmYpZp = nj*ncx + nk*(ncx*ncy) + ni - 1;
                        }
                        else {
                            cell_XmYpZp = std::numeric_limits<T2>::max();
                        }

                        if (ni < ncx && nj > 0 && nk < ncz){
                            cell_XpYmZp = (nj-1)*ncx + nk*(ncx*ncy) + ni;
                        }
                        else {
                            cell_XpYmZp = std::numeric_limits<T2>::max();
                        }

                        if (ni > 0 && nj > 0 && nk < ncz){
                            cell_XmYmZp = (nj-1)*ncx + nk*(ncx * ncy) + ni - 1;
                        }
                        else {
                            cell_XmYmZp = std::numeric_limits<T2>::max();
                        }

                        if (ni < ncx && nj < ncy && nk > 0){
                            cell_XpYpZm = nj*ncx + (nk-1)*(ncx*ncy) + ni;
                        }
                        else {
                            cell_XpYpZm = std::numeric_limits<T2>::max();
                        }

                        if (ni > 0 && nj < ncy && nk > 0){
                            cell_XmYpZm = nj*ncx + (nk-1)*(ncx*ncy) + ni - 1;
                        }
                        else {
                            cell_XmYpZm = std::numeric_limits<T2>::max();
                        }

                        if (ni < ncx && nj > 0 && nk > 0){
                            cell_XpYmZm = (nj-1)*ncx + (nk-1)*(ncx*ncy) + ni;
                        }
                        else {
                            cell_XpYmZm = std::numeric_limits<T2>::max();
                        }

                        if (ni > 0 && nj > 0 && nk > 0){
                            cell_XmYmZm = (nj-1)*ncx + (nk-1)*(ncx*ncy) + ni - 1;
                        }
                        else {
                            cell_XmYmZm = std::numeric_limits<T2>::max();
                        }

                        // Secondary nodes on x edge
                        if ( ni < ncx ) {
                            for (T2 ns=0; ns< nsnx; ++ns, ++n ) {

                                T1 xsv = xmin + ni*dx + (ns+1)*dxs;

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
                                    nodes[n].setPrimary(false);
                                }
                                else if (nj==0 && nk>0){
                                    nodes[n].setPrimary(false);
                                }
                                else if (nj>0 && nk==0){
                                    nodes[n].setPrimary(false);
                                }
                                else if (nj==0 && nk==0){
                                    nodes[n].setPrimary(false);
                                }
                            }
                        }

                        // Secondary nodes on y edge
                        if ( nj < ncy ) {
                            for (T2 ns=0; ns< nsny; ++ns, ++n ) {

                                T1 ysv = ymin + nj* dy + (ns+1)*dys;

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
                                    nodes[n].setPrimary(false);
                                }
                                else if (ni>0 && nk==0){
                                    nodes[n].setPrimary(false);
                                }
                                else if (ni==0 && nk>0){
                                    nodes[n].setPrimary(false);
                                }
                                else if (ni==0 && nk==0){
                                    nodes[n].setPrimary(false);
                                }
                            }
                        }

                        // Secondary nodes on z edge
                        if ( nk < ncz ) {
                            for (T2 ns=0; ns< nsnz; ++ns, ++n ) {

                                T1 zsv = zmin + nk* dz + (ns+1)*dzs;

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
                                    nodes[n].setPrimary(false);
                                }
                                else if (ni>0 && nj==0){
                                    nodes[n].setPrimary(false);
                                }
                                else if (ni==0 && nj>0){
                                    nodes[n].setPrimary(false);
                                }
                                else if (ni==0 && nj==0){
                                    nodes[n].setPrimary(false);
                                }
                            }
                        }

                        // Secondary nodes on the xy0 planes
                        if ( ni < ncx && nj < ncy ) {
                            for (T2 sy=0; sy < nsny; ++sy){
                                for (T2 sx=0; sx < nsnx; ++sx, n++) {

                                    T1 ysv = ymin + nj*dy + (sy+1)*dys;
                                    T1 xsv = xmin + ni*dx + (sx+1)*dxs;

                                    if ( cell_XpYpZm != std::numeric_limits<T2>::max() ) {
                                        nodes[n].pushOwner( cell_XpYpZm );
                                    }
                                    if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
                                        nodes[n].pushOwner( cell_XpYpZp );
                                    }
                                    nodes[n].setXYZindex( xsv, ysv, z, n );

                                    if (nk>0){
                                        nodes[n].setPrimary(false);
                                    }
                                    else if (nk==0){
                                        nodes[n].setPrimary(false);
                                    }
                                }
                            }
                        }

                        // Secondary nodes on the x0z planes
                        if ( ni < ncx && nk < ncz ) {
                            for(T2 sz=0; sz < nsnz; ++sz){
                                for(T2 sx=0; sx < nsnx; ++sx, n++){

                                    T1 zsv = zmin + nk*dz + (sz+1)*dzs;
                                    T1 xsv = xmin + ni*dx + (sx+1)*dxs;

                                    if ( cell_XpYmZp != std::numeric_limits<T2>::max() ) {
                                        nodes[n].pushOwner( cell_XpYmZp );
                                    }
                                    if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
                                        nodes[n].pushOwner( cell_XpYpZp );
                                    }
                                    nodes[n].setXYZindex( xsv, y, zsv, n );

                                    if (nj>0){
                                        nodes[n].setPrimary(false);
                                    }
                                    else if (nj==0){
                                        nodes[n].setPrimary(false);
                                    }
                                }
                            }
                        }

                        // Secondary nodes on the 0yz planes
                        if ( nj < ncy && nk < ncz ) {
                            for(T2 sz=0; sz < nsnz; ++sz){
                                for(T2 sy=0; sy < nsny; ++sy, n++){

                                    T1 zsv = zmin + nk*dz + (sz+1)*dzs;
                                    T1 ysv = ymin + nj*dy + (sy+1)*dys;

                                    if ( cell_XmYpZp != std::numeric_limits<T2>::max() ) {
                                        nodes[n].pushOwner( cell_XmYpZp );
                                    }
                                    if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
                                        nodes[n].pushOwner( cell_XpYpZp );
                                    }
                                    nodes[n].setXYZindex( x, ysv, zsv, n );

                                    if (ni>0){
                                        nodes[n].setPrimary(false);
                                    }
                                    else if (ni==0){
                                        nodes[n].setPrimary(false);
                                    }
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

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::interpSecondary() {
        T2 nPrimary = (ncx+1)*(ncy+1)*(ncz+1);
        for ( size_t n=nPrimary; n<nodes.size(); ++n ) {
            nodes[n].setNodeSlowness( computeSlowness({nodes[n].getX(),
                nodes[n].getY(), nodes[n].getZ()}) );
        }
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::checkPts(const std::vector<sxyz<T1>>& pts) const {

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


    template<typename T1, typename T2, typename NODE>
    T1 Grid3Drn<T1,T2,NODE>::getTraveltime(const sxyz<T1> &pt, const size_t nt) const {

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


    template<typename T1, typename T2, typename NODE>
    T1 Grid3Drn<T1,T2,NODE>::getTraveltime(const sxyz<T1>& Rx,
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
        T2 neibNo = this->neighbors[cellNo][0];
        T1 dt = computeDt(nodes[neibNo], Rx, slo);

        T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
        nodeParentRx = neibNo;
        cellParentRx = cellNo;
        for ( size_t k=1; k< this->neighbors[cellNo].size(); ++k ) {
            neibNo = this->neighbors[cellNo][k];
            dt = computeDt(nodes[neibNo], Rx, slo);
            if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
                traveltime =  nodes[neibNo].getTT(threadNo)+dt;
                nodeParentRx = neibNo;
            }
        }
        return traveltime;
    }

    
    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::grad(sxyz<T1>& g, const size_t i, const size_t j, const size_t k,
                                    const size_t nt) const {

        // compute average gradient for voxel (i,j,k)

        const size_t nnx = ncx+1;
        const size_t nny = ncy+1;

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
    void Grid3Drn<T1,T2,NODE>::gradO2(sxyz<T1>& g, const sxyz<T1> &pt,
                                      const size_t nt) const {

        // compute travel time gradient (2nd order centered operator) at point pt

        T1 p1 = pt.x - dx/2.0;
        T1 p2 = p1 + dx;
        if ( p1 <= xmin ) {  // check if on grid edge or out of grid
            p1 = pt.x + dx/2.0;  // shift pt to allow interpolating in getTraveltime
            p2 = p1 + dx;
        } else if ( p2 >= xmax ) {
            p2 = pt.x - dx/2.0;
            p1 = p2 - dx;
        }
        g.x = (getTraveltime({p2, pt.y, pt.z}, nt) - getTraveltime({p1, pt.y, pt.z}, nt)) / dx;

        p1 = pt.y - dy/2.0;
        p2 = p1 + dy;
        if ( p1 <= ymin ) {
            p1 = pt.y + dy/2.0;
            p2 = p1 + dy;
        } else if ( p2 >= ymax ) {
            p2 = pt.y - dy/2.0;
            p1 = p2 - dy;
        }
        g.y = (getTraveltime({pt.x, p2, pt.z}, nt) - getTraveltime({pt.x, p1, pt.z}, nt)) / dy;

        p1 = pt.z - dz/2.0;
        p2 = p1 + dz;
        if ( p1 <= zmin ) {
            p1 = pt.z + dz/2.0;
            p2 = p1 + dz;
        } else if ( p2 >= zmax ) {
            p2 = pt.z - dz/2.0;
            p1 = p2 - dz;
        }
        g.z = (getTraveltime({pt.x, pt.y, p2}, nt) - getTraveltime({pt.x, pt.y, p1}, nt)) / dz;

    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::grad(sxyz<T1>& g, const sxyz<T1> &pt,
                                    const size_t nt) const {

        // compute travel time gradient (4th order centered operator) at point pt

        static const T1 k1 = 1./24.;
        static const T1 k2 = 9./8.;

        T1 p1 = pt.x - dx;
        T1 p2 = p1 + 0.5*dx;
        T1 p3 = p1 + 1.5*dx;
        T1 p4 = p1 + 2.0*dx;
        if ( p1 <= xmin ) {  // check if on grid edge or out of grid
            p1 = xmin;  // shift pt to allow interpolating in getTraveltime
            p2 = p1 + 0.5*dx;
            p3 = p1 + 1.5*dx;
            p4 = p1 + 2.0*dx;
        } else if ( p4 >= xmax ) {
            p4 = xmax;
            p3 = p4 - 0.5*dx;
            p2 = p4 - 1.5*dx;
            p1 = p4 - 2.0*dx;
        }
        g.x = (k1 * getTraveltime({p1, pt.y, pt.z}, nt) -
               k2 * getTraveltime({p2, pt.y, pt.z}, nt) +
               k2 * getTraveltime({p3, pt.y, pt.z}, nt) -
               k1 * getTraveltime({p4, pt.y, pt.z}, nt)) / dx;

        p1 = pt.y - dy/2.0;
        p2 = p1 + 0.5*dy;
        p3 = p1 + 1.5*dy;
        p4 = p1 + 2.0*dy;
        if ( p1 <= ymin ) {
            p1 = ymin;
            p2 = p1 + 0.5*dy;
            p3 = p1 + 1.5*dy;
            p4 = p1 + 2.0*dy;
        } else if ( p4 >= ymax ) {
            p4 = ymax;
            p3 = p4 - 0.5*dy;
            p2 = p4 - 1.5*dy;
            p1 = p4 - 2.0*dy;
        }
        g.y = (k1 * getTraveltime({pt.x, p1, pt.z}, nt) -
               k2 * getTraveltime({pt.x, p2, pt.z}, nt) +
               k2 * getTraveltime({pt.x, p3, pt.z}, nt) -
               k1 * getTraveltime({pt.x, p4, pt.z}, nt)) / dy;

        p1 = pt.z - dz/2.0;
        p2 = p1 + 0.5*dz;
        p3 = p1 + 1.5*dz;
        p4 = p1 + 2.0*dz;
        if ( p1 <= zmin ) {
            p1 = zmin;
            p2 = p1 + 0.5*dz;
            p3 = p1 + 1.5*dz;
            p4 = p1 + 2.0*dz;
        } else if ( p4 >= zmax ) {
            p4 = zmax;
            p3 = p4 - 0.5*dz;
            p2 = p4 - 1.5*dz;
            p1 = p4 - 2.0*dz;
        }
        g.z = (k1 * getTraveltime({pt.x, pt.y, p1}, nt) -
               k2 * getTraveltime({pt.x, pt.y, p2}, nt) +
               k2 * getTraveltime({pt.x, pt.y, p3}, nt) -
               k1 * getTraveltime({pt.x, pt.y, p4}, nt)) / dz;
    }

    template<typename T1, typename T2, typename NODE>
    T1 Grid3Drn<T1,T2,NODE>::getTraveltimeFromRaypath(const std::vector<sxyz<T1>>& Tx,
                                                      const std::vector<T1>& t0,
                                                      const sxyz<T1> &Rx,
                                                      const size_t threadNo) const {
        T1 tt = 0.0;
        T1 s1, s2;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return tt;
            }
        }

        sxyz<T1> prev_pt( Rx );
        sxyz<T1> curr_pt( Rx );
        s1 = computeSlowness( curr_pt );
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

            s2 = computeSlowness( curr_pt );
            tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
            s1 = s2;
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
                        s2 = computeSlowness( Tx[ns] );
                        tt += t0[ns] + 0.5*(s1 + s2) * prev_pt.getDistance( Tx[ns] );
                    } else {
                        // to intersection
                        s2 = computeSlowness( curr_pt );
                        tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
                        s1 = s2;
                        // to Tx
                        s2 = computeSlowness( Tx[ns] );
                        tt += t0[ns] + 0.5*(s1 + s2) * curr_pt.getDistance( Tx[ns] );
                    }

                    reachedTx = true;
                }
            }
        }

        return tt;
    }


    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
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
        const T1 maxDist = sqrt( dx*dx + dy*dy + dz*dz );
        sxyz<T1> g;
#ifdef DEBUG_RP
        std::cout << "\n\nStarting raypath computation\n  curr_pt = " << curr_pt << '\n';
#endif
        bool reachedTx = false;
        while ( reachedTx == false ) {

            grad(g, curr_pt, threadNo);
            g *= -1.0;
#ifdef DEBUG_RP
            std::cout << "  g = " << g << '\n';
#endif
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
            std::cout << "  curr_pt = " << curr_pt << '\n';
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
    void Grid3Drn<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                          const std::vector<T1>& t0,
                                          const sxyz<T1> &Rx,
                                          std::vector<sxyz<T1>> &r_data,
                                          T1 &tt,
                                          const size_t threadNo) const {
        tt = 0.0;
        T1 s1, s2;

        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        sxyz<T1> curr_pt( Rx );
        s1 = computeSlowness( curr_pt );
        // distance between opposite nodes of a voxel
        const T1 maxDist = sqrt( dx*dx + dy*dy + dz*dz );
        sxyz<T1> g;
#ifdef DEBUG_RP
        std::cout << '\n' << curr_pt << '\n';
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
#ifdef DEBUG_RP
            std::cout << curr_pt << '\n';
#endif
            s2 = computeSlowness( curr_pt );
            tt += 0.5*(s1 + s2) * r_data.back().getDistance( curr_pt );
            s1 = s2;
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

                    if ( curr_pt.getDistance(r_data.back()) > dist  ||  // we do not intersect a plane
                        curr_pt == Tx[ns] ) {  // we have arrived

                        s2 = computeSlowness( Tx[ns] );
                        tt += t0[ns] + 0.5*(s1 + s2) * r_data.back().getDistance( Tx[ns] );
                        r_data.push_back( Tx[ns] );
#ifdef DEBUG_RP
                        std::cout << Tx[ns] << "\n\n";
#endif
                    } else {
                        // to intersection
                        s2 = computeSlowness( curr_pt );
                        tt += 0.5*(s1 + s2) * r_data.back().getDistance( curr_pt );
                        r_data.push_back( curr_pt );
#ifdef DEBUG_RP
                        std::cout << curr_pt << '\n';
#endif
                        s1 = s2;
                        // to Tx
                        s2 = computeSlowness( Tx[ns] );
                        tt += t0[ns] + 0.5*(s1 + s2) * curr_pt.getDistance( Tx[ns] );
                        r_data.push_back( Tx[ns] );
#ifdef DEBUG_RP
                        std::cout << Tx[ns] << "\n\n";
#endif
                    }

                    reachedTx = true;
                }
            }
        }
    }


    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                          const std::vector<T1>& t0,
                                          const sxyz<T1> &Rx,
                                          std::vector<sijv<T1>>& m_data,
                                          T1 &tt,
                                          const size_t RxNo,
                                          const size_t threadNo) const {
        tt = 0.0;
        T1 s1, s2;

        const size_t nnx = ncx+1;
        const size_t nny = ncy+1;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                return;
            }
        }

        sxyz<T1> curr_pt( Rx ), prev_pt( Rx ), mid_pt;
        s1 = computeSlowness( curr_pt );
        sijv<T1> m;
        m.i = RxNo;

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

            s2 = computeSlowness( curr_pt );
            tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
            s1 = s2;
            prev_pt = curr_pt;

            // compute terms of matrix M
            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
            T1 s = computeSlowness(mid_pt);
            s *= s;
            T1 ds = curr_pt.getDistance( prev_pt );

            size_t ix = (mid_pt.x-xmin)/dx;
            size_t iy = (mid_pt.y-ymin)/dy;
            size_t iz = (mid_pt.z-zmin)/dz;
            for ( size_t ii=0; ii<2; ++ii ) {
                for ( size_t jj=0; jj<2; ++jj ) {
                    for ( size_t kk=0; kk<2; ++kk ) {
                        size_t iv = ix+ii;
                        size_t jv = iy+jj;
                        size_t kv = iz+kk;
                        T1 dvdv = (1. - std::abs(mid_pt.x - iv*dx)/dx) *
                        (1. - std::abs(mid_pt.y - jv*dy)/dy) *
                        (1. - std::abs(mid_pt.z - kv*dz)/dz);

                        m.j = (kv*nny+jv)*nnx+iv;
                        m.v = -s * ds * dvdv;

                        bool found = false;
                        for ( size_t nm=0; nm<m_data.size(); ++nm ) {
                            if ( m_data[nm].j == m.j ) {
                                m_data[nm].v += m.v;
                                found = true;
                                break;
                            }
                        }
                        if ( found == false ) {
                            m_data.push_back(m);
                        }
                    }
                }
            }


            // are we close enough to one of the Tx nodes ?
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

                    if ( curr_pt.getDistance(prev_pt) > dist  ||  // we do not intersect a plane
                        curr_pt == Tx[ns] ) {  // we have arrived

                        s2 = computeSlowness( Tx[ns] );
                        tt += t0[ns] + 0.5*(s1 + s2) * prev_pt.getDistance( Tx[ns] );

                        // compute terms of matrix M
                        mid_pt = static_cast<T1>(0.5)*(Tx[ns] + prev_pt);
                        s = computeSlowness(mid_pt);
                        s *= s;
                        ds = Tx[ns].getDistance( prev_pt );

                        ix = (mid_pt.x-xmin)/dx;
                        iy = (mid_pt.y-ymin)/dy;
                        iz = (mid_pt.z-zmin)/dz;
                        for ( size_t ii=0; ii<2; ++ii ) {
                            for ( size_t jj=0; jj<2; ++jj ) {
                                for ( size_t kk=0; kk<2; ++kk ) {
                                    size_t iv = ix+ii;
                                    size_t jv = iy+jj;
                                    size_t kv = iz+kk;
                                    T1 dvdv = (1. - std::abs(mid_pt.x - iv*dx)/dx) *
                                    (1. - std::abs(mid_pt.y - jv*dy)/dy) *
                                    (1. - std::abs(mid_pt.z - kv*dz)/dz);

                                    m.j = (kv*nny+jv)*nnx+iv;
                                    m.v = -s * ds * dvdv;

                                    bool found = false;
                                    for ( size_t nm=0; nm<m_data.size(); ++nm ) {
                                        if ( m_data[nm].j == m.j ) {
                                            m_data[nm].v += m.v;
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found == false ) {
                                        m_data.push_back(m);
                                    }
                                }
                            }
                        }
                    } else {
                        // to intersection
                        s2 = computeSlowness( curr_pt );
                        tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
                        s1 = s2;

                        // compute terms of matrix M
                        mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                        s = computeSlowness(mid_pt);
                        s *= s;
                        ds = curr_pt.getDistance( prev_pt );

                        ix = (mid_pt.x-xmin)/dx;
                        iy = (mid_pt.y-ymin)/dy;
                        iz = (mid_pt.z-zmin)/dz;
                        for ( size_t ii=0; ii<2; ++ii ) {
                            for ( size_t jj=0; jj<2; ++jj ) {
                                for ( size_t kk=0; kk<2; ++kk ) {
                                    size_t iv = ix+ii;
                                    size_t jv = iy+jj;
                                    size_t kv = iz+kk;
                                    T1 dvdv = (1. - std::abs(mid_pt.x - iv*dx)/dx) *
                                    (1. - std::abs(mid_pt.y - jv*dy)/dy) *
                                    (1. - std::abs(mid_pt.z - kv*dz)/dz);

                                    m.j = (kv*nny+jv)*nnx+iv;
                                    m.v = -s * ds * dvdv;

                                    bool found = false;
                                    for ( size_t nm=0; nm<m_data.size(); ++nm ) {
                                        if ( m_data[nm].j == m.j ) {
                                            m_data[nm].v += m.v;
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found == false ) {
                                        m_data.push_back(m);
                                    }

                                }
                            }
                        }

                        // to Tx
                        s2 = computeSlowness( Tx[ns] );
                        tt += t0[ns] + 0.5*(s1 + s2) * curr_pt.getDistance( Tx[ns] );

                        // compute terms of matrix M
                        mid_pt = static_cast<T1>(0.5)*(Tx[ns] + curr_pt);
                        s = computeSlowness(mid_pt);
                        s *= s;
                        ds = Tx[ns].getDistance( curr_pt );

                        ix = (mid_pt.x-xmin)/dx;
                        iy = (mid_pt.y-ymin)/dy;
                        iz = (mid_pt.z-zmin)/dz;
                        for ( size_t ii=0; ii<2; ++ii ) {
                            for ( size_t jj=0; jj<2; ++jj ) {
                                for ( size_t kk=0; kk<2; ++kk ) {
                                    size_t iv = ix+ii;
                                    size_t jv = iy+jj;
                                    size_t kv = iz+kk;
                                    T1 dvdv = (1. - std::abs(mid_pt.x - iv*dx)/dx) *
                                    (1. - std::abs(mid_pt.y - jv*dy)/dy) *
                                    (1. - std::abs(mid_pt.z - kv*dz)/dz);

                                    m.j = (kv*nny+jv)*nnx+iv;
                                    m.v = -s * ds * dvdv;

                                    bool found = false;
                                    for ( size_t nm=0; nm<m_data.size(); ++nm ) {
                                        if ( m_data[nm].j == m.j ) {
                                            m_data[nm].v += m.v;
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found == false ) {
                                        m_data.push_back(m);
                                    }
                                }
                            }
                        }
                    }

                    reachedTx = true;
                }
            }
        }
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                          const std::vector<T1>& t0,
                                          const sxyz<T1> &Rx,
                                          std::vector<siv<T1>> &l_data,
                                          T1 &tt,
                                          const size_t threadNo) const {
        tt = 0.0;
        T1 s1, s2;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        sxyz<T1> curr_pt( Rx ), prev_pt( Rx );
        s1 = computeSlowness( curr_pt );
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
            s2 = computeSlowness( curr_pt );
            tt += 0.5*(s1 + s2) * cell.v;
            s1 = s2;
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
                        s2 = computeSlowness( Tx[ns] );
                        tt += 0.5*(s1 + s2) * cell.v;

                    } else {
                        // to intersection
                        mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
                        cell.i = getCellNo(mid_pt);
                        cell.v = curr_pt.getDistance(prev_pt);
                        l_data.push_back(cell);
                        s2 = computeSlowness( curr_pt );
                        tt += 0.5*(s1 + s2) * cell.v;

                        s1 = s2;
                        // to Tx
                        cell.i = getCellNo(Tx[ns]);
                        cell.v = Tx[ns].getDistance(curr_pt);
                        l_data.push_back(cell);
                        s2 = computeSlowness( Tx[ns] );
                        tt += 0.5*(s1 + s2) * cell.v;
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                          const std::vector<T1>& t0,
                                          const sxyz<T1> &Rx,
                                          std::vector<sxyz<T1>> &r_data,
                                          std::vector<siv<T1>> &l_data,
                                          T1 &tt,
                                          const size_t threadNo) const {
        tt = 0.0;
        T1 s1, s2;
        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        sxyz<T1> curr_pt( Rx );
        s1 = computeSlowness( curr_pt );
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
            s2 = computeSlowness( curr_pt );
            tt += 0.5*(s1 + s2) * cell.v;
            s1 = s2;
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
                        s2 = computeSlowness( Tx[ns] );
                        tt += 0.5*(s1 + s2) * cell.v;
                        r_data.push_back( Tx[ns] );
                    } else {
                        // to intersection
                        mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
                        cell.i = getCellNo(mid_pt);
                        cell.v = curr_pt.getDistance(r_data.back());
                        l_data.push_back(cell);
                        s2 = computeSlowness( curr_pt );
                        tt += 0.5*(s1 + s2) * cell.v;
                        r_data.push_back( curr_pt );

                        s1 = s2;
                        // to Tx
                        cell.i = getCellNo(Tx[ns]);
                        cell.v = Tx[ns].getDistance(r_data.back());
                        l_data.push_back(cell);
                        s2 = computeSlowness( Tx[ns] );
                        tt += 0.5*(s1 + s2) * cell.v;
                        r_data.push_back( Tx[ns] );
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                          const std::vector<T1>& t0,
                                          const sxyz<T1> &Rx,
                                          std::vector<sxyz<T1>>& r_data,
                                          std::vector<sijv<T1>>& m_data,
                                          T1 &tt,
                                          const size_t RxNo,
                                          const size_t threadNo) const {
        tt = 0.0;
        T1 s1, s2;

        const size_t nnx = ncx+1;
        const size_t nny = ncy+1;

        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                return;
            }
        }

        sxyz<T1> curr_pt( Rx ), prev_pt, mid_pt;
        s1 = computeSlowness( curr_pt );
        sijv<T1> m;
        m.i = RxNo;

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

            prev_pt = r_data.back();
            s2 = computeSlowness( curr_pt );
            tt += 0.5*(s1 + s2) * r_data.back().getDistance( curr_pt );
            s1 = s2;
            r_data.push_back( curr_pt );

            // compute terms of matrix M
            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
            T1 s = computeSlowness(mid_pt);
            s *= s;
            T1 ds = curr_pt.getDistance( prev_pt );

            size_t ix = (mid_pt.x-xmin)/dx;
            size_t iy = (mid_pt.y-ymin)/dy;
            size_t iz = (mid_pt.z-zmin)/dz;
            for ( size_t ii=0; ii<2; ++ii ) {
                for ( size_t jj=0; jj<2; ++jj ) {
                    for ( size_t kk=0; kk<2; ++kk ) {
                        size_t iv = ix+ii;
                        size_t jv = iy+jj;
                        size_t kv = iz+kk;
                        T1 dvdv = (1. - std::abs(mid_pt.x - iv*dx)/dx) *
                        (1. - std::abs(mid_pt.y - jv*dy)/dy) *
                        (1. - std::abs(mid_pt.z - kv*dz)/dz);

                        m.j = (kv*nny+jv)*nnx+iv;
                        m.v = -s * ds * dvdv;

                        bool found = false;
                        for ( size_t nm=0; nm<m_data.size(); ++nm ) {
                            if ( m_data[nm].j == m.j ) {
                                m_data[nm].v += m.v;
                                found = true;
                                break;
                            }
                        }
                        if ( found == false ) {
                            m_data.push_back(m);
                        }
                    }
                }
            }


            // are we close enough to one of the Tx nodes ?
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

                    if ( curr_pt.getDistance(r_data.back()) > dist  ||  // we do not intersect a plane
                        curr_pt == Tx[ns] ) {  // we have arrived

                        prev_pt = r_data.back();
                        s2 = computeSlowness( Tx[ns] );
                        tt += t0[ns] + 0.5*(s1 + s2) * r_data.back().getDistance( Tx[ns] );
                        r_data.push_back( Tx[ns] );

                        // compute terms of matrix M
                        mid_pt = static_cast<T1>(0.5)*(Tx[ns] + prev_pt);
                        s = computeSlowness(mid_pt);
                        s *= s;
                        ds = Tx[ns].getDistance( prev_pt );

                        ix = (mid_pt.x-xmin)/dx;
                        iy = (mid_pt.y-ymin)/dy;
                        iz = (mid_pt.z-zmin)/dz;
                        for ( size_t ii=0; ii<2; ++ii ) {
                            for ( size_t jj=0; jj<2; ++jj ) {
                                for ( size_t kk=0; kk<2; ++kk ) {
                                    size_t iv = ix+ii;
                                    size_t jv = iy+jj;
                                    size_t kv = iz+kk;
                                    T1 dvdv = (1. - std::abs(mid_pt.x - iv*dx)/dx) *
                                    (1. - std::abs(mid_pt.y - jv*dy)/dy) *
                                    (1. - std::abs(mid_pt.z - kv*dz)/dz);

                                    m.j = (kv*nny+jv)*nnx+iv;
                                    m.v = -s * ds * dvdv;

                                    bool found = false;
                                    for ( size_t nm=0; nm<m_data.size(); ++nm ) {
                                        if ( m_data[nm].j == m.j ) {
                                            m_data[nm].v += m.v;
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found == false ) {
                                        m_data.push_back(m);
                                    }
                                }
                            }
                        }
                    } else {
                        // to intersection
                        s2 = computeSlowness( curr_pt );
                        tt += 0.5*(s1 + s2) * r_data.back().getDistance( curr_pt );
                        r_data.push_back( curr_pt );
                        s1 = s2;

                        prev_pt = r_data.back();
                        // compute terms of matrix M
                        mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                        s = computeSlowness(mid_pt);
                        s *= s;
                        ds = curr_pt.getDistance( prev_pt );

                        ix = (mid_pt.x-xmin)/dx;
                        iy = (mid_pt.y-ymin)/dy;
                        iz = (mid_pt.z-zmin)/dz;
                        for ( size_t ii=0; ii<2; ++ii ) {
                            for ( size_t jj=0; jj<2; ++jj ) {
                                for ( size_t kk=0; kk<2; ++kk ) {
                                    size_t iv = ix+ii;
                                    size_t jv = iy+jj;
                                    size_t kv = iz+kk;
                                    T1 dvdv = (1. - std::abs(mid_pt.x - iv*dx)/dx) *
                                    (1. - std::abs(mid_pt.y - jv*dy)/dy) *
                                    (1. - std::abs(mid_pt.z - kv*dz)/dz);

                                    m.j = (kv*nny+jv)*nnx+iv;
                                    m.v = -s * ds * dvdv;

                                    bool found = false;
                                    for ( size_t nm=0; nm<m_data.size(); ++nm ) {
                                        if ( m_data[nm].j == m.j ) {
                                            m_data[nm].v += m.v;
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found == false ) {
                                        m_data.push_back(m);
                                    }

                                }
                            }
                        }

                        // to Tx
                        prev_pt = r_data.back();
                        s2 = computeSlowness( Tx[ns] );
                        tt += t0[ns] + 0.5*(s1 + s2) * curr_pt.getDistance( Tx[ns] );
                        r_data.push_back( Tx[ns] );

                        // compute terms of matrix M
                        mid_pt = static_cast<T1>(0.5)*(Tx[ns] + prev_pt);
                        s = computeSlowness(mid_pt);
                        s *= s;
                        ds = Tx[ns].getDistance( prev_pt );

                        ix = (mid_pt.x-xmin)/dx;
                        iy = (mid_pt.y-ymin)/dy;
                        iz = (mid_pt.z-zmin)/dz;
                        for ( size_t ii=0; ii<2; ++ii ) {
                            for ( size_t jj=0; jj<2; ++jj ) {
                                for ( size_t kk=0; kk<2; ++kk ) {
                                    size_t iv = ix+ii;
                                    size_t jv = iy+jj;
                                    size_t kv = iz+kk;
                                    T1 dvdv = (1. - std::abs(mid_pt.x - iv*dx)/dx) *
                                    (1. - std::abs(mid_pt.y - jv*dy)/dy) *
                                    (1. - std::abs(mid_pt.z - kv*dz)/dz);

                                    m.j = (kv*nny+jv)*nnx+iv;
                                    m.v = -s * ds * dvdv;

                                    bool found = false;
                                    for ( size_t nm=0; nm<m_data.size(); ++nm ) {
                                        if ( m_data[nm].j == m.j ) {
                                            m_data[nm].v += m.v;
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found == false ) {
                                        m_data.push_back(m);
                                    }
                                }
                            }
                        }
                    }

                    reachedTx = true;
                }
            }
        }
    }


    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::getRaypath_old(const std::vector<sxyz<T1>>& Tx,
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
        const T1 maxDist = sqrt( dx*dx + dy*dy + dz*dz );

        getIJK(curr_pt, iIn, jIn, kIn);

        bool reachedTx = false;
        while ( reachedTx == false ) {

            bool onNode=false;
            bool onEdgeX=false;
            bool onEdgeY=false;
            bool onEdgeZ=false;

            if ( std::abs(remainder(curr_pt.x,dx))<small &&
                std::abs(remainder(curr_pt.y,dy))<small &&
                std::abs(remainder(curr_pt.z,dz))<small ) {
                onNode = true;
            } else if ( std::abs(remainder(curr_pt.y,dy))<small &&
                       std::abs(remainder(curr_pt.z,dz))<small ) {
                onEdgeX = true;
            } else if ( std::abs(remainder(curr_pt.x,dx))<small &&
                       std::abs(remainder(curr_pt.z,dz))<small ) {
                onEdgeY = true;
            } else if ( std::abs(remainder(curr_pt.x,dx))<small &&
                       std::abs(remainder(curr_pt.y,dy))<small ) {
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

                iOut = i;
                jOut = boost::math::sign(gOut.y)<0.0 ? j-1 : j;
                kOut = boost::math::sign(gOut.z)<0.0 ? k-1 : k;

                // planes we will intersect
                T1 xp = xmin + dx*(i + boost::math::sign(gOut.x)>0.0 ? 1.0 : 0.0);
                T1 yp = ymin + dy*(j + boost::math::sign(gOut.y));
                T1 zp = zmin + dz*(k + boost::math::sign(gOut.z));

                if ( std::abs(xp-curr_pt.x)<small) {
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

                iOut = boost::math::sign(gOut.x)<0.0 ? i-1 : i;
                jOut = j;
                kOut = boost::math::sign(gOut.z)<0.0 ? k-1 : k;

                // planes we will intersect
                T1 xp = xmin + dx*(i + boost::math::sign(gOut.x));
                T1 yp = ymin + dy*(j + boost::math::sign(gOut.y)>0.0 ? 1.0 : 0.0);
                T1 zp = zmin + dz*(k + boost::math::sign(gOut.z));

                if ( std::abs(yp-curr_pt.y)<small) {
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

                iOut = boost::math::sign(gOut.x)<0.0 ? i-1 : i;
                jOut = boost::math::sign(gOut.y)<0.0 ? j-1 : j;
                kOut = k;

                // planes we will intersect
                T1 xp = xmin + dx*(i + boost::math::sign(gOut.x));
                T1 yp = ymin + dy*(j + boost::math::sign(gOut.y));
                T1 zp = zmin + dz*(k + boost::math::sign(gOut.z)>0.0 ? 1.0 : 0.0);

                if ( std::abs(zp-curr_pt.z)<small) {
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
                        throw std::runtime_error("Error while computing raypaths: raypath not converging!");
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

                // planes we will intersect
                T1 xp = xmin + dx*(iIn + boost::math::sign(gOut.x)>0.0 ? 1.0 : 0.0);
                T1 yp = ymin + dy*(jIn + boost::math::sign(gOut.y)>0.0 ? 1.0 : 0.0);
                T1 zp = zmin + dz*(kIn + boost::math::sign(gOut.z)>0.0 ? 1.0 : 0.0);

                if ( std::abs(xp-curr_pt.x)<small) {
                    xp += dx*boost::math::sign(gOut.x);
                    iOut += boost::math::sign(gOut.x);
                }
                if ( std::abs(yp-curr_pt.y)<small) {
                    yp += dy*boost::math::sign(gOut.y);
                    jOut += boost::math::sign(gOut.y);
                }
                if ( std::abs(zp-curr_pt.z)<small) {
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
    T1 Grid3Drn<T1,T2,NODE>::computeSlowness(const sxyz<T1>& pt) const {

        const size_t nnx = ncx+1;
        const size_t nny = ncy+1;
        const size_t nnz = ncz+1;

        // are we on an node, an edge or a face?
        ptrdiff_t onX = -1;
        ptrdiff_t onY = -1;
        ptrdiff_t onZ = -1;
        for ( ptrdiff_t n=0; n<nnx; ++n ) {
            if ( std::abs(pt.x - (xmin+n*dx)) < small2 ) {
                onX = n;
                break;
            }
        }
        for ( ptrdiff_t n=0; n<nny; ++n ) {
            if ( std::abs(pt.y - (ymin+n*dy)) < small2 ) {
                onY = n;
                break;
            }
        }
        for ( ptrdiff_t n=0; n<nnz; ++n ) {
            if ( std::abs(pt.z - (zmin+n*dz)) < small2 ) {
                onZ = n;
                break;
            }
        }

        if ( onX!=-1 && onY!=-1 && onZ!=-1 ) {
            return nodes[(onZ*nny+onY)*nnx+onX].getNodeSlowness();
        } else if ( onX!=-1 && onY!=-1 ) {
            T2 k = static_cast<T2>( small + (pt.z-zmin)/dz );
            T1 s[2];
            T1 x[3];
            if ( processVel ) {
                s[0] = 1.0 / nodes[(k*nny+onY)*nnx+onX].getNodeSlowness();
                s[1] = 1.0 / nodes[((k+1)*nny+onY)*nnx+onX].getNodeSlowness();
            } else {
                s[0] = nodes[(k*nny+onY)*nnx+onX].getNodeSlowness();
                s[1] = nodes[((k+1)*nny+onY)*nnx+onX].getNodeSlowness();
            }
            x[0] = pt.z;
            x[1] = zmin + k*dz;
            x[2] = zmin + (k+1)*dz;

            if ( processVel )
                return 1.0 / Interpolator<T1>::linear(x, s);
            else
                return Interpolator<T1>::linear(x, s);

        } else if ( onX!=-1 && onZ!=-1 ) {
            T2 j = static_cast<T2>( small + (pt.y-ymin)/dy );
            T1 s[2];
            T1 x[3];
            if ( processVel ) {
                s[0] = 1.0 / nodes[(onZ*nny+j)*nnx+onX].getNodeSlowness();
                s[1] = 1.0 / nodes[(onZ*nny+j+1)*nnx+onX].getNodeSlowness();
            } else {
                s[0] = nodes[(onZ*nny+j)*nnx+onX].getNodeSlowness();
                s[1] = nodes[(onZ*nny+j+1)*nnx+onX].getNodeSlowness();
            }
            x[0] = pt.y;
            x[1] = ymin + j*dy;
            x[2] = ymin + (j+1)*dy;

            if ( processVel )
                return 1.0 / Interpolator<T1>::linear(x, s);
            else
                return Interpolator<T1>::linear(x, s);

        } else if ( onY!=-1 && onZ!=-1 ) {
            T2 i = static_cast<T2>( small + (pt.x-xmin)/dx );
            T1 s[2];
            T1 x[3];
            if ( processVel ) {
                s[0] = 1.0 / nodes[(onZ*nny+onY)*nnx+i].getNodeSlowness();
                s[1] = 1.0 / nodes[(onZ*nny+onY)*nnx+i+1].getNodeSlowness();
            } else {
                s[0] = nodes[(onZ*nny+onY)*nnx+i].getNodeSlowness();
                s[1] = nodes[(onZ*nny+onY)*nnx+i+1].getNodeSlowness();
            }
            x[0] = pt.x;
            x[1] = xmin + i*dx;
            x[2] = xmin + (i+1)*dx;

            if ( processVel )
                return 1.0 / Interpolator<T1>::linear(x, s);
            else
                return Interpolator<T1>::linear(x, s);
        } else if ( onX!=-1 ) {
            T2 j = static_cast<T2>( small + (pt.y-ymin)/dy );
            T2 k = static_cast<T2>( small + (pt.z-zmin)/dz );
            T1 s[4];
            T1 x[3];
            T1 y[3];
            if ( processVel ) {
                s[0] = 1.0 / nodes[((k  )*nny+j  )*nnx+onX].getNodeSlowness();
                s[1] = 1.0 / nodes[((k+1)*nny+j  )*nnx+onX].getNodeSlowness();
                s[2] = 1.0 / nodes[((k  )*nny+j+1)*nnx+onX].getNodeSlowness();
                s[3] = 1.0 / nodes[((k+1)*nny+j+1)*nnx+onX].getNodeSlowness();
            } else {
                s[0] = nodes[((k  )*nny+j  )*nnx+onX].getNodeSlowness();
                s[1] = nodes[((k+1)*nny+j  )*nnx+onX].getNodeSlowness();
                s[2] = nodes[((k  )*nny+j+1)*nnx+onX].getNodeSlowness();
                s[3] = nodes[((k+1)*nny+j+1)*nnx+onX].getNodeSlowness();
            }
            x[0] = pt.y;
            y[0] = pt.z;
            x[1] = ymin + j*dy;
            y[1] = zmin + k*dz;
            x[2] = ymin + (j+1)*dy;
            y[2] = zmin + (k+1)*dz;

            if ( processVel )
                return 1.0 / Interpolator<T1>::bilinear(x, y, s);
            else
                return Interpolator<T1>::bilinear(x, y, s);

        } else if ( onY!=-1 ) {
            T2 i = static_cast<T2>( small + (pt.x-xmin)/dx );
            T2 k = static_cast<T2>( small + (pt.z-zmin)/dz );
            T1 s[4];
            T1 x[3];
            T1 y[3];
            if ( processVel ) {
                s[0] = 1.0 / nodes[((k  )*nny+onY)*nnx+i  ].getNodeSlowness();
                s[1] = 1.0 / nodes[((k+1)*nny+onY)*nnx+i  ].getNodeSlowness();
                s[2] = 1.0 / nodes[((k  )*nny+onY)*nnx+i+1].getNodeSlowness();
                s[3] = 1.0 / nodes[((k+1)*nny+onY)*nnx+i+1].getNodeSlowness();
            } else {
                s[0] = nodes[((k  )*nny+onY)*nnx+i  ].getNodeSlowness();
                s[1] = nodes[((k+1)*nny+onY)*nnx+i  ].getNodeSlowness();
                s[2] = nodes[((k  )*nny+onY)*nnx+i+1].getNodeSlowness();
                s[3] = nodes[((k+1)*nny+onY)*nnx+i+1].getNodeSlowness();
            }
            x[0] = pt.x;
            y[0] = pt.z;
            x[1] = xmin + i*dx;
            y[1] = zmin + k*dz;
            x[2] = xmin + (i+1)*dx;
            y[2] = zmin + (k+1)*dz;

            if ( processVel )
                return 1.0 / Interpolator<T1>::bilinear(x, y, s);
            else
                return Interpolator<T1>::bilinear(x, y, s);

        } else if ( onZ!=-1 ) {
            T2 i = static_cast<T2>( small + (pt.x-xmin)/dx );
            T2 j = static_cast<T2>( small + (pt.y-ymin)/dy );
            T1 s[4];
            T1 x[3];
            T1 y[3];
            if ( processVel ) {
                s[0] = 1.0 / nodes[(onZ*nny+j  )*nnx+i  ].getNodeSlowness();
                s[1] = 1.0 / nodes[(onZ*nny+j+1)*nnx+i  ].getNodeSlowness();
                s[2] = 1.0 / nodes[(onZ*nny+j  )*nnx+i+1].getNodeSlowness();
                s[3] = 1.0 / nodes[(onZ*nny+j+1)*nnx+i+1].getNodeSlowness();
            } else {
                s[0] = nodes[(onZ*nny+j  )*nnx+i  ].getNodeSlowness();
                s[1] = nodes[(onZ*nny+j+1)*nnx+i  ].getNodeSlowness();
                s[2] = nodes[(onZ*nny+j  )*nnx+i+1].getNodeSlowness();
                s[3] = nodes[(onZ*nny+j+1)*nnx+i+1].getNodeSlowness();
            }
            x[0] = pt.x;
            y[0] = pt.y;
            x[1] = xmin + i*dx;
            y[1] = ymin + j*dy;
            x[2] = xmin + (i+1)*dx;
            y[2] = ymin + (j+1)*dy;

            if ( processVel )
                return 1.0 / Interpolator<T1>::bilinear(x, y, s);
            else
                return Interpolator<T1>::bilinear(x, y, s);

        } else {
            T2 i = static_cast<T2>( small + (pt.x-xmin)/dx );
            T2 j = static_cast<T2>( small + (pt.y-ymin)/dy );
            T2 k = static_cast<T2>( small + (pt.z-zmin)/dz );
            T1 s[8];
            T1 x[3];
            T1 y[3];
            T1 z[3];

            if ( processVel ) {
                s[0] = 1.0 / nodes[((k  )*nny+j  )*nnx+i  ].getNodeSlowness();
                s[1] = 1.0 / nodes[((k+1)*nny+j  )*nnx+i  ].getNodeSlowness();
                s[2] = 1.0 / nodes[((k  )*nny+j+1)*nnx+i  ].getNodeSlowness();
                s[3] = 1.0 / nodes[((k+1)*nny+j+1)*nnx+i  ].getNodeSlowness();
                s[4] = 1.0 / nodes[((k  )*nny+j  )*nnx+i+1].getNodeSlowness();
                s[5] = 1.0 / nodes[((k+1)*nny+j  )*nnx+i+1].getNodeSlowness();
                s[6] = 1.0 / nodes[((k  )*nny+j+1)*nnx+i+1].getNodeSlowness();
                s[7] = 1.0 / nodes[((k+1)*nny+j+1)*nnx+i+1].getNodeSlowness();
            } else {
                s[0] = nodes[((k  )*nny+j  )*nnx+i  ].getNodeSlowness();
                s[1] = nodes[((k+1)*nny+j  )*nnx+i  ].getNodeSlowness();
                s[2] = nodes[((k  )*nny+j+1)*nnx+i  ].getNodeSlowness();
                s[3] = nodes[((k+1)*nny+j+1)*nnx+i  ].getNodeSlowness();
                s[4] = nodes[((k  )*nny+j  )*nnx+i+1].getNodeSlowness();
                s[5] = nodes[((k+1)*nny+j  )*nnx+i+1].getNodeSlowness();
                s[6] = nodes[((k  )*nny+j+1)*nnx+i+1].getNodeSlowness();
                s[7] = nodes[((k+1)*nny+j+1)*nnx+i+1].getNodeSlowness();
            }
            x[0] = pt.x;
            y[0] = pt.y;
            z[0] = pt.z;
            x[1] = xmin + i*dx;
            y[1] = ymin + j*dy;
            z[1] = zmin + k*dz;
            x[2] = xmin + (i+1)*dx;
            y[2] = ymin + (j+1)*dy;
            z[2] = zmin + (k+1)*dz;

            if ( processVel )
                return 1.0 / Interpolator<T1>::trilinear(x, y, z, s);
            else
                return Interpolator<T1>::trilinear(x, y, z, s);

        }

    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::saveTT(const std::string & fname, const int all,
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

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::loadTT(const std::string & fname, const int all,
                                      const size_t nt,
                                      const int format) const {

        if ( format == 1 ) {
            std::string filename = fname+".dat";
            std::ifstream fin(filename.c_str());
            T1 x, y, z, tt;
            for ( T2 n=0; n<nodes.size(); ++n ) {
                if ( nodes[n].isPrimary() || all==1 ) {
                    fin >> x >> y >> z >> tt;
                    nodes[n].setTT(tt, nt);
                }
            }
            fin.close();
        } else if ( format == 2 ) {
#ifdef VTK
            std::string filename = fname+".vtr";

            vtkSmartPointer<vtkXMLRectilinearGridReader> reader =
            vtkSmartPointer<vtkXMLRectilinearGridReader>::New();

            reader->SetFileName( filename.c_str() );
            reader->Update();

            for ( size_t n=0; n<nodes.size(); ++n ) {
                if ( nodes[n].isPrimary() ) {
                    vtkIdType id = reader->GetOutput()->FindPoint(nodes[n].getX(), nodes[n].getY(), nodes[n].getZ());
                    nodes[n].setTT(reader->GetOutput()->GetPointData()->GetArray("Travel time")->GetTuple1(id), nt);
                }
            }
#else
            std::cerr << "VTK not included during compilation.\n";
#endif
        } else if ( format == 3 ) {
            std::string filename = fname+".bin";
            std::ifstream fin(filename.c_str(),std::ios::binary);
            for ( T2 n=0; n<nodes.size(); ++n ) {
                if ( nodes[n].isPrimary() || all==1 ) {
                    T1 tmp[4];
                    fin.read( (char*)tmp, 4*sizeof(T1) );
                    nodes[n].setTT(tmp[3], nt);
                }
            }
            fin.close();
        } else {
            throw std::runtime_error("Unsupported format for traveltimes");
        }
    }


    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::sweep(const std::vector<bool>& frozen,
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
    void Grid3Drn<T1,T2,NODE>::update_node(const size_t i, const size_t j, const size_t k,
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
    void Grid3Drn<T1,T2,NODE>::sweep_weno3(const std::vector<bool>& frozen,
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
    void Grid3Drn<T1,T2,NODE>::update_node_weno3(const size_t i,
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
    void Grid3Drn<T1,T2,NODE>::initFSM(const std::vector<sxyz<T1>>& Tx,
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
                                            T1 tt = t0[n] + nodes[nnn].getDistance(Tx[n]) * nodes[nnn].getNodeSlowness();
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
                                        T1 tt = t0[n] + nodes[nnn].getDistance(Tx[n]) * nodes[nnn].getNodeSlowness();
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

#ifdef VTK
    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::saveModelVTR(const std::string &fname,
                                            const bool saveSlowness) const {

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

        newScalars->SetNumberOfComponents(1);
        newScalars->SetNumberOfTuples( rgrid->GetNumberOfPoints() );

        if ( saveSlowness ) {
            newScalars->SetName("Slowness");
            for ( size_t n=0; n<nodes.size(); ++n ) {
                if ( nodes[n].isPrimary() ) {
                    vtkIdType id = rgrid->FindPoint(nodes[n].getX(), nodes[n].getY(), nodes[n].getZ());
                    newScalars->SetTuple1(id, nodes[n].getNodeSlowness() );
                }
            }
        } else {
            newScalars->SetName("Velocity");
            for ( size_t n=0; n<nodes.size(); ++n ) {
                if ( nodes[n].isPrimary() ) {
                    vtkIdType id = rgrid->FindPoint(nodes[n].getX(), nodes[n].getY(), nodes[n].getZ());
                    newScalars->SetTuple1(id, 1./nodes[n].getNodeSlowness() );
                }
            }
        }
        rgrid->GetPointData()->SetScalars(newScalars);

        vtkSmartPointer<vtkXMLRectilinearGridWriter> writer =
        vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();

        writer->SetFileName( fname.c_str() );
        writer->SetInputData( rgrid );
        writer->SetDataModeToBinary();
        writer->Update();
    }
#endif


}

#endif

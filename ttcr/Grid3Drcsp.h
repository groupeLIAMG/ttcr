//
//  Grid3drcsp.h
//  ttcr
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

//  Copyright (c) 2015 Bernard Giroux, Benoît Larouche.
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

#ifndef ttcr_Grid3Drcsp_h
#define ttcr_Grid3Drcsp_h

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <vector>

#include "Grid3Drc.h"
#include "Node3Dcsp.h"

namespace ttcr {
    
    template<typename T1, typename T2, typename CELL>
    class Grid3Drcsp : public Grid3Drc<T1,T2,Node3Dcsp<T1,T2>,CELL> {
    public:
        
        /* Constructor Format:
         Grid3Drc<T1,T2>::Grid3Drc(nb cells in x, nb cells in y, nb cells in z,
         x cells size, y cells size, z cells size,
         x origin, y origin, z origin,
         nb sec. cells in x, nb sec. cells in y, nb sec. cells in z,
         index of the thread)
         */
        Grid3Drcsp(const T2 nx, const T2 ny, const T2 nz,
                   const T1 ddx, const T1 ddy, const T1 ddz,
                   const T1 minx, const T1 miny, const T1 minz,
                   const T2 nnx, const T2 nny, const T2 nnz,
                   const bool ttrp, const size_t nt) :
        Grid3Drc<T1,T2,Node3Dcsp<T1,T2>,CELL>(nx, ny, nz, ddx, ddy, ddz, minx, miny, minz, ttrp, nt),
        nsnx(nnx), nsny(nny), nsnz(nnz)
        {
            buildGridNodes();
            this->template buildGridNeighbors<Node3Dcsp<T1,T2>>(this->nodes);
        }
        
        ~Grid3Drcsp() {
            
        }
        
        
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                     const std::vector<T1>& t0,
                     const std::vector<sxyz<T1>>& Rx,
                     std::vector<T1>& traveltimes,
                     const size_t threadNo=0) const;
        
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                     const std::vector<T1>& t0,
                     const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                     std::vector<std::vector<T1>*>& traveltimes,
                     const size_t threadNo=0) const;
        
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                     const std::vector<T1>& t0,
                     const std::vector<sxyz<T1>>& Rx,
                     std::vector<T1>& traveltimes,
                     std::vector<std::vector<sxyz<T1>>>& r_data,
                     const size_t threadNo=0) const;
        
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                     const std::vector<T1>& t0,
                     const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                     std::vector<std::vector<T1>*>& traveltimes,
                     std::vector<std::vector<std::vector<sxyz<T1>>>*>& r_data,
                     const size_t threadNo=0) const;
        
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                     const std::vector<T1>& t0,
                     const std::vector<sxyz<T1>>& Rx,
                     std::vector<T1>& traveltimes,
                     std::vector<std::vector<sxyz<T1>>>& r_data,
                     std::vector<std::vector<siv<T1>>>& l_data,
                     const size_t threadNo=0) const;

        void raytrace(const std::vector<sxyz<T1>>& Tx,
                     const std::vector<T1>& t0,
                     const std::vector<sxyz<T1>>& Rx,
                     std::vector<T1>& traveltimes,
                     std::vector<std::vector<siv<T1>>>& l_data,
                     const size_t threadNo=0) const;

        void raytrace2(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<sxyz<T1>>& Rx,
                      std::vector<T1>& traveltimes,
                      const size_t threadNo=0) const;
        
        const T2 getNsnx() const { return nsnx; }
        const T2 getNsny() const { return nsny; }
        const T2 getNsnz() const { return nsnz; }
        
    private:
        T2 nsnx;                 // number of secondary nodes in x
        T2 nsny;                 // number of secondary nodes in y
        T2 nsnz;                 // number of secondary nodes in z
        
        void buildGridNodes();
        
        void initQueue(const std::vector<sxyz<T1>>& Tx,
                       const std::vector<T1>& t0,
                       std::priority_queue<Node3Dcsp<T1,T2>*,
                       std::vector<Node3Dcsp<T1,T2>*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<Node3Dcsp<T1,T2>>& txNodes,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       const size_t threadNo) const;
        
        void propagate(std::priority_queue<Node3Dcsp<T1,T2>*,
                       std::vector<Node3Dcsp<T1,T2>*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       size_t threadNo) const;
        
        void propagate_lw(std::priority_queue<Node3Dcsp<T1,T2>*,
                          std::vector<Node3Dcsp<T1,T2>*>,
                          CompareNodePtr<T1>>& queue,
                          std::vector<bool>& inQueue,
                          std::vector<bool>& frozen,
                          size_t threadNo) const;
        
        void prepropagate(const Node3Dcsp<T1,T2>& node,
                          std::priority_queue<Node3Dcsp<T1,T2>*,
                          std::vector<Node3Dcsp<T1,T2>*>,
                          CompareNodePtr<T1>>& queue,
                          std::vector<bool>& inQueue,
                          std::vector<bool>& frozen,
                          size_t threadNo) const;
        
        void initQueue2(const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        std::priority_queue<Node3Dcsp<T1,T2>*,
                        std::deque<Node3Dcsp<T1,T2>*>,
                        CompareNodePtr<T1>>& queue,
                        std::vector<Node3Dcsp<T1,T2>>& txNodes,
                        std::vector<bool>& inQueue,
                        std::vector<bool>& frozen,
                        const size_t threadNo) const;
        
        void propagate2(std::priority_queue<Node3Dcsp<T1,T2>*,
                        std::deque<Node3Dcsp<T1,T2>*>,
                        CompareNodePtr<T1>>& queue,
                        std::vector<bool>& inQueue,
                        std::vector<bool>& frozen,
                        size_t threadNo) const;
        
        void prepropagate2(const Node3Dcsp<T1,T2>& node,
                           std::priority_queue<Node3Dcsp<T1,T2>*,
                           std::deque<Node3Dcsp<T1,T2>*>,
                           CompareNodePtr<T1>>& queue,
                           std::vector<bool>& inQueue,
                           std::vector<bool>& frozen,
                           size_t threadNo) const;
    };
    
    
    template<typename T1, typename T2, typename CELL>
    void Grid3Drcsp<T1,T2,CELL>::buildGridNodes() {
        
        
        this->nodes.resize(// secondary nodes on the edges
                           this->ncx*nsnx*((this->ncy+1)*(this->ncz+1)) +
                           this->ncy*nsny*((this->ncx+1)*(this->ncz+1)) +
                           this->ncz*nsnz*((this->ncx+1)*(this->ncy+1)) +
                           // secondary nodes on the faces
                           (nsnx*nsny)*(this->ncx*this->ncy*(this->ncz+1))+
                           (nsnx*nsnz)*(this->ncx*this->ncz*(this->ncy+1))+
                           (nsny*nsnz)*(this->ncy*this->ncz*(this->ncx+1))+
                           // primary nodes
                           (this->ncx+1) * (this->ncy+1) * (this->ncz+1),
                           Node3Dcsp<T1,T2>(this->nThreads));
        
        // Create the grid, assign a number for each node and find the owners
        // Nodes and cells are first indexed in z, then y, and x.
        // Secondary nodes are placed on the faces and edges of every cells.
        // Ex: the node in "node[A]=(i,j,k)" is followed by the node in
        // "node[A+1]=(i+dx,j,k)"
        
        T1 dxs = this->dx/(this->nsnx+1); 	// distance between secondary nodes in x
        T1 dys = this->dy/(this->nsny+1);
        T1 dzs = this->dz/(this->nsnz+1);
        
        T2 cXmYmZm; 	// cell in the (x-,y-,z-) direction from the node
        T2 cXpYmZm; 	// cell in the (x+,y-,z-) direction from the node
        T2 cXmYpZm;
        T2 cXpYpZm;
        T2 cXmYmZp;
        T2 cXpYmZp;
        T2 cXmYpZp;
        T2 cXpYpZp;
        
        T2 n = 0;
        for ( T2 nk=0; nk<=this->ncz; ++nk ) {
            
            T1 z = this->zmin + nk*this->dz;
            
            for ( T2 nj=0; nj<=this->ncy; ++nj ) {
                
                T1 y = this->ymin + nj*this->dy;
                
                for ( T2 ni=0; ni<=this->ncx; ++ni, ++n ){
                    
                    T1 x = this->xmin + ni*this->dx;
                    
                    // Find the adjacent cells for each primary node
                    
                    if (ni < this->ncx && nj < this->ncy && nk < this->ncz){
                        cXpYpZp = nj*this->ncx + nk*(this->ncx*this->ncy) + ni;
                    }
                    else {
                        cXpYpZp = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni > 0 && nj < this->ncy && nk < this->ncz){
                        cXmYpZp = nj*this->ncx + nk*(this->ncx*this->ncy) + ni - 1;
                    }
                    else {
                        cXmYpZp = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni < this->ncx && nj > 0 && nk < this->ncz){
                        cXpYmZp = (nj-1)*this->ncx + nk*(this->ncx*this->ncy) + ni;
                    }
                    else {
                        cXpYmZp = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni > 0 && nj > 0 && nk < this->ncz){
                        cXmYmZp = (nj-1)*this->ncx + nk*(this->ncx * this->ncy) + ni - 1;
                    }
                    else {
                        cXmYmZp = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni < this->ncx && nj < this->ncy && nk > 0){
                        cXpYpZm = nj*this->ncx + (nk-1)*(this->ncx*this->ncy) + ni;
                    }
                    else {
                        cXpYpZm = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni > 0 && nj < this->ncy && nk > 0){
                        cXmYpZm = nj*this->ncx + (nk-1)*(this->ncx*this->ncy) + ni - 1;
                    }
                    else {
                        cXmYpZm = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni < this->ncx && nj > 0 && nk > 0){
                        cXpYmZm = (nj-1)*this->ncx + (nk-1)*(this->ncx*this->ncy) + ni;
                    }
                    else {
                        cXpYmZm = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni > 0 && nj > 0 && nk > 0){
                        cXmYmZm = (nj-1)*this->ncx + (nk-1)*(this->ncx*this->ncy) + ni-1;
                    }
                    else {
                        cXmYmZm = std::numeric_limits<T2>::max();
                    }
                    
                    
                    // Index the primary nodes owners
                    
                    if ( cXmYmZm != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cXmYmZm );
                    }
                    if ( cXpYmZm != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cXpYmZm );
                    }
                    if ( cXmYpZm != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cXmYpZm );
                    }
                    if ( cXpYpZm != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cXpYpZm );
                    }
                    if ( cXmYmZp != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cXmYmZp );
                    }
                    if ( cXpYmZp != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cXpYmZp );
                    }
                    if ( cXmYpZp != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cXmYpZp );
                    }
                    if ( cXpYpZp != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cXpYpZp );
                    }
                    
                    this->nodes[n].setXYZindex( x, y, z, n );
                    this->nodes[n].setPrimary(true);
                }
            }
        }

        for ( T2 nk=0; nk<=this->ncz; ++nk ) {

            T1 z = this->zmin + nk*this->dz;

            for ( T2 nj=0; nj<=this->ncy; ++nj ) {

                T1 y = this->ymin + nj*this->dy;

                for ( T2 ni=0; ni<=this->ncx; ++ni ){

                    T1 x = this->xmin + ni*this->dx;

                    // Find the adjacent cells for each primary node

                    if (ni < this->ncx && nj < this->ncy && nk < this->ncz){
                        cXpYpZp = nj*this->ncx + nk*(this->ncx*this->ncy) + ni;
                    }
                    else {
                        cXpYpZp = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj < this->ncy && nk < this->ncz){
                        cXmYpZp = nj*this->ncx + nk*(this->ncx*this->ncy) + ni - 1;
                    }
                    else {
                        cXmYpZp = std::numeric_limits<T2>::max();
                    }

                    if (ni < this->ncx && nj > 0 && nk < this->ncz){
                        cXpYmZp = (nj-1)*this->ncx + nk*(this->ncx*this->ncy) + ni;
                    }
                    else {
                        cXpYmZp = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj > 0 && nk < this->ncz){
                        cXmYmZp = (nj-1)*this->ncx + nk*(this->ncx * this->ncy) + ni - 1;
                    }
                    else {
                        cXmYmZp = std::numeric_limits<T2>::max();
                    }

                    if (ni < this->ncx && nj < this->ncy && nk > 0){
                        cXpYpZm = nj*this->ncx + (nk-1)*(this->ncx*this->ncy) + ni;
                    }
                    else {
                        cXpYpZm = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj < this->ncy && nk > 0){
                        cXmYpZm = nj*this->ncx + (nk-1)*(this->ncx*this->ncy) + ni - 1;
                    }
                    else {
                        cXmYpZm = std::numeric_limits<T2>::max();
                    }

                    if (ni < this->ncx && nj > 0 && nk > 0){
                        cXpYmZm = (nj-1)*this->ncx + (nk-1)*(this->ncx*this->ncy) + ni;
                    }
                    else {
                        cXpYmZm = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj > 0 && nk > 0){
                        cXmYmZm = (nj-1)*this->ncx + (nk-1)*(this->ncx*this->ncy) + ni-1;
                    }
                    else {
                        cXmYmZm = std::numeric_limits<T2>::max();
                    }
                    
                    // Secondary nodes on x edge
                    if ( ni < this->ncx ) {
                        for (T2 ns=0; ns< this->nsnx; ++ns, ++n ) {
                            
                            T1 xsv = this->xmin + ni* this->dx + (ns+1)*dxs;
                            
                            if ( cXpYmZm != std::numeric_limits<T2>::max() )
                            {
                                this->nodes[n].pushOwner( cXpYmZm );
                            }
                            if ( cXpYpZm != std::numeric_limits<T2>::max() )
                            {
                                this->nodes[n].pushOwner( cXpYpZm );
                            }
                            if ( cXpYmZp != std::numeric_limits<T2>::max() )
                            {
                                this->nodes[n].pushOwner( cXpYmZp );
                            }
                            if ( cXpYpZp != std::numeric_limits<T2>::max() )
                            {
                                this->nodes[n].pushOwner( cXpYpZp );
                            }
                            this->nodes[n].setXYZindex( xsv, y, z, n );
                        }
                    }
                    
                    // Secondary nodes on y edge
                    if ( nj < this->ncy ) {
                        for (T2 ns=0; ns< this->nsny; ++ns, ++n ) {
                            
                            T1 ysv = this->ymin + nj* this->dy + (ns+1)*dys;
                            
                            if ( cXmYpZm != std::numeric_limits<T2>::max() )
                            {
                                this->nodes[n].pushOwner( cXmYpZm );
                            }
                            if ( cXpYpZm != std::numeric_limits<T2>::max() )
                            {
                                this->nodes[n].pushOwner( cXpYpZm );
                            }
                            if ( cXmYpZp != std::numeric_limits<T2>::max() )
                            {
                                this->nodes[n].pushOwner( cXmYpZp );
                            }
                            if ( cXpYpZp != std::numeric_limits<T2>::max() )
                            {
                                this->nodes[n].pushOwner( cXpYpZp );
                            }
                            this->nodes[n].setXYZindex( x, ysv, z, n );
                        }
                    }
                    
                    // Secondary nodes on z edge
                    if ( nk < this->ncz ) {
                        for (T2 ns=0; ns< this->nsnz; ++ns, ++n ) {
                            
                            T1 zsv = this->zmin + nk* this->dz + (ns+1)*dzs;
                            
                            if ( cXmYmZp != std::numeric_limits<T2>::max() )
                            {
                                this->nodes[n].pushOwner( cXmYmZp );
                            }
                            if ( cXpYmZp != std::numeric_limits<T2>::max() )
                            {
                                this->nodes[n].pushOwner( cXpYmZp );
                            }
                            if ( cXmYpZp != std::numeric_limits<T2>::max() )
                            {
                                this->nodes[n].pushOwner( cXmYpZp );
                            }
                            if ( cXpYpZp != std::numeric_limits<T2>::max() )
                            {
                                this->nodes[n].pushOwner( cXpYpZp );
                            }
                            this->nodes[n].setXYZindex( x, y, zsv, n );
                        }
                    }
                    
                    // Secondary nodes on the xy0 planes
                    if ( ni < this->ncx && nj < this->ncy ) {
                        for ( T2 sy=0; sy < this->nsny; ++sy ) {
                            for ( T2 sx=0; sx < this->nsnx; ++sx, n++ ) {
                                
                                T1 ysv = this->ymin+ nj* this->dy+ (sy+1)*dys;
                                T1 xsv = this->xmin+ ni* this->dx+ (sx+1)*dxs;
                                
                                if ( cXpYpZm != std::numeric_limits<T2>::max() )
                                {
                                    this->nodes[n].pushOwner( cXpYpZm );
                                }
                                if ( cXpYpZp != std::numeric_limits<T2>::max() )
                                {
                                    this->nodes[n].pushOwner( cXpYpZp );
                                }
                                this->nodes[n].setXYZindex( xsv, ysv, z, n );
                            }
                        }
                    }
                    
                    // Secondary nodes on the x0z planes
                    if ( ni < this->ncx && nk < this->ncz ) {
                        for ( T2 sz=0; sz < this->nsnz; ++sz ) {
                            for ( T2 sx=0; sx < this->nsnx; ++sx, n++ ) {
                                
                                T1 zsv = this->zmin+ nk* this->dz+ (sz+1)*dzs;
                                T1 xsv = this->xmin+ ni* this->dx+ (sx+1)*dxs;
                                
                                if ( cXpYmZp != std::numeric_limits<T2>::max() )
                                {
                                    this->nodes[n].pushOwner( cXpYmZp );
                                }
                                if ( cXpYpZp != std::numeric_limits<T2>::max() )
                                {
                                    this->nodes[n].pushOwner( cXpYpZp );
                                }
                                this->nodes[n].setXYZindex( xsv, y, zsv, n );
                            }
                        }
                    }
                    
                    // Secondary nodes on the 0yz planes
                    if ( nj < this->ncy && nk < this->ncz ) {
                        for ( T2 sz=0; sz < this->nsnz; ++sz ) {
                            for ( T2 sy=0; sy < this->nsny; ++sy, n++ ) {
                                
                                T1 zsv = this->zmin+ nk* this->dz+ (sz+1)*dzs;
                                T1 ysv = this->ymin+ nj* this->dy+ (sy+1)*dys;
                                
                                if ( cXmYpZp != std::numeric_limits<T2>::max() )
                                {
                                    this->nodes[n].pushOwner( cXmYpZp );
                                }
                                if ( cXpYpZp != std::numeric_limits<T2>::max() )
                                {
                                    this->nodes[n].pushOwner( cXpYpZp );
                                }
                                this->nodes[n].setXYZindex( x, ysv, zsv, n );
                            }
                        }
                    }
                }
            }
        }
        // sanity check
        if ( n != this->nodes.size() ) {
            std::cerr << "Error building grid, wrong number of nodes\n";
            abort();
        }
    }
    
    
    
    template<typename T1, typename T2, typename CELL>
    void Grid3Drcsp<T1,T2,CELL>::initQueue(const std::vector<sxyz<T1>>& Tx,
                                           const std::vector<T1>& t0,
                                           std::priority_queue<Node3Dcsp<T1,T2>*,
                                           std::vector<Node3Dcsp<T1,T2>*>,
                                           CompareNodePtr<T1>>& queue,
                                           std::vector<Node3Dcsp<T1,T2>>& txNodes,
                                           std::vector<bool>& inQueue,
                                           std::vector<bool>& frozen,
                                           const size_t threadNo) const {
        
        //Find the starting nodes of the transmitters Tx and start the queue list
        for ( size_t n=0; n<Tx.size(); ++n ) {
            bool found = false;
            for ( size_t nn=0; nn<this->nodes.size(); ++nn ) {
                if ( this->nodes[nn] == Tx[n] ) {
                    found = true;
                    this->nodes[nn].setTT( t0[n], threadNo );
                    frozen[nn] = true;
                    
                    prepropagate(this->nodes[nn], queue, inQueue, frozen, threadNo); // See description in the function declaration
                    
                    //	queue.push( &(this->nodes[nn]) );   	//Don't use if prepropagate is used
                    //	inQueue[nn] = true;				//Don't use if prepropagate is used
                    
                    break;
                }
            }
            if ( found==false ) {
                // If Tx[n] is not on a node, we create a new node and initialize the queue:
                txNodes.push_back( Node3Dcsp<T1,T2>(Tx[n].x, Tx[n].y, Tx[n].z,
                                                    static_cast<T2>(this->nodes.size()+txNodes.size()-1),
                                                    this->nThreads));
                txNodes.back().setTT(t0[n], threadNo);
                txNodes.back().pushOwner( this->getCellNo(Tx[n]) );
                frozen.push_back( true );
                
                prepropagate(txNodes.back(), queue, inQueue, frozen, threadNo); // See description in the function declaration
                
                //	queue.push( &(txNodes.back()) );	//Don't use if prepropagate is used
                //	inQueue.push_back( true );			//Don't use if prepropagate is used
                
            }
        }
    }
    
    template<typename T1, typename T2, typename CELL>
    void Grid3Drcsp<T1,T2,CELL>::propagate( std::priority_queue<Node3Dcsp<T1,T2>*,
                                           std::vector<Node3Dcsp<T1,T2>*>,
                                           CompareNodePtr<T1>>& queue,
                                           std::vector<bool>& inQueue,
                                           std::vector<bool>& frozen,
                                           size_t threadNo) const {
        
        while ( !queue.empty() ) {
            const Node3Dcsp<T1,T2>* source = queue.top();
            queue.pop();
            inQueue[ source->getGridIndex() ] = false;
            frozen[ source->getGridIndex() ] = true;
            
            for ( size_t no=0; no<source->getOwners().size(); ++no ) {
                T2 cellNo = source->getOwners()[no];
                for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                    size_t neibNo = this->neighbors[cellNo][k];
                    if ( neibNo == source->getGridIndex() || frozen[neibNo] ) {
                        continue;
                    }
                    
                    T1 ttsource= source->getTT( threadNo );
                    if (ttsource < this->nodes[neibNo].getTT(threadNo)){
                        // Compute dt
                        T1 dt = this->cells.computeDt(*source, this->nodes[neibNo], cellNo);
                        
                        if ( ttsource +dt < this->nodes[neibNo].getTT( threadNo ) ) {
                            this->nodes[neibNo].setTT( ttsource +dt, threadNo );
                            this->nodes[neibNo].setnodeParent( source->getGridIndex(),
                                                              threadNo );
                            this->nodes[neibNo].setCellParent( cellNo, threadNo );
                            
                            if ( !inQueue[neibNo] ) {
                                queue.push( &(this->nodes[neibNo]) );
                                inQueue[neibNo] = true;
                            }
                        }
                    }
                }
            }
        }
    }
    
    template<typename T1, typename T2, typename CELL>
    void Grid3Drcsp<T1,T2,CELL>::propagate_lw(std::priority_queue<Node3Dcsp<T1,T2>*,
                                              std::vector<Node3Dcsp<T1,T2>*>,
                                              CompareNodePtr<T1>>& queue,
                                              std::vector<bool>& inQueue,
                                              std::vector<bool>& frozen,
                                              size_t threadNo) const {
        // lightweight method where cell/node parent are not stored
        while ( !queue.empty() ) {
            const Node3Dcsp<T1,T2>* source = queue.top();
            queue.pop();
            inQueue[ source->getGridIndex() ] = false;
            frozen[ source->getGridIndex() ] = true;
            
            for ( size_t no=0; no<source->getOwners().size(); ++no ) {
                T2 cellNo = source->getOwners()[no];
                for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                    size_t neibNo = this->neighbors[cellNo][k];
                    if ( neibNo == source->getGridIndex() || frozen[neibNo] ) {
                        continue;
                    }
                    
                    T1 ttsource= source->getTT( threadNo );
//                    if (ttsource < this->nodes[neibNo].getTT(threadNo)){
                        // Compute dt
                        T1 dt = this->cells.computeDt(*source, this->nodes[neibNo], cellNo);
                        
                        if ( ttsource+dt < this->nodes[neibNo].getTT( threadNo ) ) {
                            this->nodes[neibNo].setTT( ttsource+dt, threadNo );
                            if ( !inQueue[neibNo] ) {
                                queue.push( &(this->nodes[neibNo]) );
                                inQueue[neibNo] = true;
                            }
                        }
//                    }
                }
            }
        }
    }
    
    template<typename T1, typename T2, typename CELL>
    void Grid3Drcsp<T1,T2,CELL>::prepropagate(const Node3Dcsp<T1,T2>& node,
                                              std::priority_queue<Node3Dcsp<T1,T2>*,
                                              std::vector<Node3Dcsp<T1,T2>*>,
                                              CompareNodePtr<T1>>& queue,
                                              std::vector<bool>& inQueue,
                                              std::vector<bool>& frozen,
                                              size_t threadNo) const {
        
        // This function can be used to "prepropagate" each Tx nodes one first time
        // during "initQueue", before running "propagate".
        // When a Tx source node seems to be lost in the queue and is not
        // propagated, corrupting the entire traveltime table,
        // this function force the propagation of every source points and can
        // solve the problem.
        
        for ( size_t no=0; no<node.getOwners().size(); ++no ) {
            T2 cellNo = node.getOwners()[no];
            for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                size_t neibNo = this->neighbors[cellNo][k];
                if ( neibNo == node.getGridIndex() || frozen[neibNo] ) {
                    continue;
                }
                
                // compute dt
                T1 dt = this->cells.computeDt(node, this->nodes[neibNo], cellNo);
                
                if ( node.getTT( threadNo )+dt < this->nodes[neibNo].getTT( threadNo ) ) {
                    this->nodes[neibNo].setTT( node.getTT( threadNo )+dt, threadNo );
                    this->nodes[neibNo].setnodeParent( node.getGridIndex(), threadNo );
                    this->nodes[neibNo].setCellParent( cellNo, threadNo );
                    
                    if ( !inQueue[neibNo] ) {
                        queue.push( &(this->nodes[neibNo]) );
                        inQueue[neibNo] = true;
                    }
                }
            }
        }
    }
    
    template<typename T1, typename T2, typename CELL>
    void Grid3Drcsp<T1,T2,CELL>::initQueue2(const std::vector<sxyz<T1>>& Tx,
                                            const std::vector<T1>& t0,
                                            std::priority_queue<Node3Dcsp<T1,T2>*,
                                            std::deque<Node3Dcsp<T1,T2>*>,
                                            CompareNodePtr<T1>>& queue,
                                            std::vector<Node3Dcsp<T1,T2>>& txNodes,
                                            std::vector<bool>& inQueue,
                                            std::vector<bool>& frozen,
                                            const size_t threadNo) const {
        
        //Find the starting nodes of the transmitters Tx and start the queue list
        for ( size_t n=0; n<Tx.size(); ++n ) {
            bool found = false;
            for ( size_t nn=0; nn<this->nodes.size(); ++nn ) {
                if ( this->nodes[nn] == Tx[n] ) {
                    found = true;
                    this->nodes[nn].setTT( t0[n], threadNo );
                    frozen[nn] = true;
                    
                    prepropagate2(this->nodes[nn], queue, inQueue, frozen, threadNo); // See description in the function declaration
                    
                    //	queue.push( &(this->nodes[nn]) );   	//Don't use if prepropagate is used
                    //	inQueue[nn] = true;				//Don't use if prepropagate is used
                    
                    break;
                }
            }
            if ( found==false ) {
                // If Tx[n] is not on a node, we create a new node and initialize the queue:
                txNodes.push_back( Node3Dcsp<T1,T2>(t0[n], Tx[n].x, Tx[n].y, Tx[n].z, this->nThreads, threadNo));
                txNodes.back().pushOwner( this->getCellNo(Tx[n]) );
                txNodes.back().setGridIndex( static_cast<T2>(this->nodes.size()+txNodes.size()-1) );
                frozen.push_back( true );
                
                prepropagate2(txNodes.back(), queue, inQueue, frozen, threadNo); // See description in the function declaration
                
                //	queue.push( &(txNodes.back()) );	//Don't use if prepropagate is used
                //	inQueue.push_back( true );			//Don't use if prepropagate is used
                
            }
        }
    }
    
    template<typename T1, typename T2, typename CELL>
    void Grid3Drcsp<T1,T2,CELL>::propagate2( std::priority_queue<Node3Dcsp<T1,T2>*,
                                            std::deque<Node3Dcsp<T1,T2>*>,
                                            CompareNodePtr<T1>>& queue,
                                            std::vector<bool>& inQueue,
                                            std::vector<bool>& frozen,
                                            size_t threadNo) const {
        
        while ( !queue.empty() ) {
            const Node3Dcsp<T1,T2>* source = queue.top();
            queue.pop();
            inQueue[ source->getGridIndex() ] = false;
            for ( size_t no=0; no<source->getOwners().size(); ++no ) {
                T2 cellNo = source->getOwners()[no];
                for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                    T2 neibNo = this->neighbors[cellNo][k];
                    if ( neibNo == source->getGridIndex() || frozen[neibNo] ) {
                        continue;
                    }
                    
                    T1 ttsource= source->getTT( threadNo );
                    if (ttsource < this->nodes[neibNo].getTT(threadNo)){
                        // Compute dt
                        T1 dt = this->cells.computeDt(source, &(this->nodes[neibNo]), cellNo);
                        
                        if ( ttsource +dt < this->nodes[neibNo].getTT( threadNo ) ) {
                            this->nodes[neibNo].setTT( ttsource +dt, threadNo );
                            this->nodes[neibNo].setnodeParent( source->getGridIndex(),
                                                              threadNo );
                            this->nodes[neibNo].setCellParent( cellNo, threadNo );
                            
                            if ( !inQueue[neibNo] ) {
                                queue.push( &(this->nodes[neibNo]) );
                                inQueue[neibNo] = true;
                            }
                        }
                    }
                    
                }
            }
        }
    }
    
    template<typename T1, typename T2, typename CELL>
    void Grid3Drcsp<T1,T2,CELL>::prepropagate2(const Node3Dcsp<T1,T2>& node,
                                               std::priority_queue<Node3Dcsp<T1,T2>*,
                                               std::deque<Node3Dcsp<T1,T2>*>,
                                               CompareNodePtr<T1>>& queue,
                                               std::vector<bool>& inQueue,
                                               std::vector<bool>& frozen,
                                               size_t threadNo) const {
        
        // This function can be used to "prepropagate" each Tx nodes one first time
        // during "initQueue", before running "propagate".
        // When a Tx source node seems to be lost in the queue and is not
        // propagated, corrupting the entire traveltime table,
        // this function force the propagation of every source points and can
        // solve the problem.
        
        for ( size_t no=0; no<node.getOwners().size(); ++no ) {
            T2 cellNo = node.getOwners()[no];
            for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                size_t neibNo = this->neighbors[cellNo][k];
                if ( neibNo == node.getGridIndex() || frozen[neibNo] ) {
                    continue;
                }
                
                // compute dt
                T1 dt = this->cells.computeDt(&node, &(this->nodes[neibNo]), cellNo);
                
                if ( node.getTT( threadNo )+dt < this->nodes[neibNo].getTT( threadNo ) ) {
                    this->nodes[neibNo].setTT( node.getTT( threadNo )+dt, threadNo );
                    this->nodes[neibNo].setnodeParent( node.getGridIndex(), threadNo );
                    this->nodes[neibNo].setCellParent( cellNo, threadNo );
                    
                    if ( !inQueue[neibNo] ) {
                        queue.push( &(this->nodes[neibNo]) );
                        inQueue[neibNo] = true;
                    }
                }
            }
        }
    }
    
    template<typename T1, typename T2, typename CELL>
    void Grid3Drcsp<T1,T2,CELL>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                          const std::vector<T1>& t0,
                                          const std::vector<sxyz<T1>>& Rx,
                                          std::vector<T1>& traveltimes,
                                          const size_t threadNo) const {
        
        // Primary function
        
        this->checkPts(Tx);
        this->checkPts(Rx);
        
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }
        
        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node3Dcsp<T1,T2>*, std::vector<Node3Dcsp<T1,T2>*>,
        CompareNodePtr<T1>> queue(cmp);
        // txNodes: Extra nodes if the sources points are not on an existing node
        std::vector<Node3Dcsp<T1,T2>> txNodes;
        // inQueue lists the nodes waiting in the queue
        std::vector<bool> inQueue( this->nodes.size(), false );
        // Tx sources nodes are "frozen" and their traveltime can't be modified
        std::vector<bool> frozen( this->nodes.size(), false );
        
        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
        
        propagate(queue, inQueue, frozen, threadNo);
        
        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        
//        if ( this->tt_from_rp ) {
//            for (size_t n=0; n<Rx.size(); ++n) {
//                traveltimes[n] = this->getTraveltimeFromRaypath(Tx, t0, Rx[n], threadNo);
//            }
//        } else {
            for (size_t n=0; n<Rx.size(); ++n) {
                bool onTx = false;
                for (size_t nt=0; nt<Tx.size(); ++nt) {
                    if ( Tx[nt] == Rx[n] ) {
                        onTx = true;
                        traveltimes[n] = t0[nt];
                        break;
                    }
                }
                if (onTx) { continue; }
                traveltimes[n] = this->getTraveltime(Rx[n], this->nodes, threadNo);
            }
//        }
    }
    
    template<typename T1, typename T2, typename CELL>
    void Grid3Drcsp<T1,T2,CELL>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                          const std::vector<T1>& t0,
                                          const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                                          std::vector<std::vector<T1>*>& traveltimes,
                                          const size_t threadNo) const {
        
        this->checkPts(Tx);
        for ( size_t n=0; n<Rx.size(); ++n )
            this->checkPts(*Rx[n]);
        
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }
        
        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node3Dcsp<T1,T2>*, std::vector<Node3Dcsp<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );
        
        std::vector<Node3Dcsp<T1,T2>> txNodes;
        std::vector<bool> inQueue( this->nodes.size(), false );
        std::vector<bool> frozen( this->nodes.size(), false );
        
        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
        
        propagate(queue, inQueue, frozen, threadNo);
        
        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        
//        if ( this->tt_from_rp ) {
//            for (size_t nr=0; nr<Rx.size(); ++nr) {
//                traveltimes[nr]->resize( Rx[nr]->size() );
//                for (size_t n=0; n<Rx[nr]->size(); ++n)
//                    (*traveltimes[nr])[n] = this->getTraveltimeFromRaypath(Tx, t0, (*Rx[nr])[n], threadNo);
//            }
//        } else {
            for (size_t nr=0; nr<Rx.size(); ++nr) {
                traveltimes[nr]->resize( Rx[nr]->size() );
                for (size_t n=0; n<Rx[nr]->size(); ++n) {
                    bool onTx = false;
                    for (size_t nt=0; nt<Tx.size(); ++nt) {
                        if ( Tx[nt] == (*Rx[nr])[n] ) {
                            onTx = true;
                            (*traveltimes[nr])[n] = t0[nt];
                            break;
                        }
                    }
                    if (onTx) { continue; }
                    (*traveltimes[nr])[n] = this->getTraveltime((*Rx[nr])[n], this->nodes, threadNo);
                }
            }
//        }
    }
    
    template<typename T1, typename T2, typename CELL>
    void Grid3Drcsp<T1,T2,CELL>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                          const std::vector<T1>& t0,
                                          const std::vector<sxyz<T1>>& Rx,
                                          std::vector<T1>& traveltimes,
                                          std::vector<std::vector<sxyz<T1>>>& r_data,
                                          const size_t threadNo) const {
        
        // Primary function
        
        this->checkPts(Tx);
        this->checkPts(Rx);
        
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }
        
        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node3Dcsp<T1,T2>*, std::vector<Node3Dcsp<T1,T2>*>,
        CompareNodePtr<T1>> queue(cmp);
        // txNodes: Extra nodes if the sources points are not on an existing node
        std::vector<Node3Dcsp<T1,T2>> txNodes;
        // inQueue lists the nodes waiting in the queue
        std::vector<bool> inQueue( this->nodes.size(), false );
        // Tx sources nodes are "frozen" and their traveltime can't be modified
        std::vector<bool> frozen( this->nodes.size(), false );
        
        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
        
        propagate(queue, inQueue, frozen, threadNo);
        
        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        if ( r_data.size() != Rx.size() ) {
            r_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<r_data.size(); ++ni ) {
            r_data[ni].resize( 0 );
        }
        T2 nodeParentRx;
        T2 cellParentRx;
        
        for (size_t n=0; n<Rx.size(); ++n) {
            bool onTx = false;
            for (size_t nt=0; nt<Tx.size(); ++nt) {
                if ( Tx[nt] == Rx[n] ) {
                    onTx = true;
                    traveltimes[n] = t0[nt];
                    r_data[n].push_back(Tx[nt]);
                    break;
                }
            }
            if (onTx) { continue; }
            
            traveltimes[n] = this->getTraveltime(Rx[n], this->nodes, nodeParentRx, cellParentRx,
                                                 threadNo);
            
            // Rx are in nodes (not txNodes)
            std::vector<Node3Dcsp<T1,T2>> *node_p;
            node_p = &(this->nodes);
            
            std::vector<sxyz<T1>> r_tmp;
            T2 iChild, iParent = nodeParentRx;
            sxyz<T1> child;
            
            // store the son's coord
            child.x = Rx[n].x;
            child.y = Rx[n].y;
            child.z = Rx[n].z;
            while ( (*node_p)[iParent].getNodeParent(threadNo) != std::numeric_limits<T2>::max() ) {
                
                r_tmp.push_back( child );
                
                // we now go up in time - parent becomes the child of grand'pa
                iChild = iParent;
                child.x = (*node_p)[iChild].getX();
                child.y = (*node_p)[iChild].getY();
                child.z = (*node_p)[iChild].getZ();
                
                // grand'pa is now papa
                iParent = (*node_p)[iChild].getNodeParent(threadNo);
                if ( iParent >= this->nodes.size() ) {
                    node_p = &txNodes;
                    iParent -= this->nodes.size();
                }
                else {
                    node_p = &(this->nodes);
                }
            }
            
            // parent is now at Tx
            r_tmp.push_back( child );
            
            // finally, store Tx position
            child.x = (*node_p)[iParent].getX();
            child.y = (*node_p)[iParent].getY();
            child.z = (*node_p)[iParent].getZ();
            r_tmp.push_back( child );
            
            // the order should be from Tx to Rx, so we reorder...
            iParent = static_cast<T2>(r_tmp.size());
            r_data[n].resize( r_tmp.size() );
            for ( size_t nn=0; nn<r_data[n].size(); ++nn ) {
                r_data[n][nn].x = r_tmp[ iParent-1-nn ].x;
                r_data[n][nn].y = r_tmp[ iParent-1-nn ].y;
                r_data[n][nn].z = r_tmp[ iParent-1-nn ].z;
            }
        }
    }
    
    template<typename T1, typename T2, typename CELL>
    void Grid3Drcsp<T1,T2,CELL>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                          const std::vector<T1>& t0,
                                          const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                                          std::vector<std::vector<T1>*>& traveltimes,
                                          std::vector<std::vector<std::vector<sxyz<T1>>>*>& r_data,
                                          const size_t threadNo) const {
        
        this->checkPts(Tx);
        for ( size_t n=0; n<Rx.size(); ++n )
            this->checkPts(*Rx[n]);
        
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }
        
        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node3Dcsp<T1,T2>*, std::vector<Node3Dcsp<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );
        
        std::vector<Node3Dcsp<T1,T2>> txNodes;
        std::vector<bool> inQueue( this->nodes.size(), false );
        std::vector<bool> frozen( this->nodes.size(), false );
        
        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
        
        propagate(queue, inQueue, frozen, threadNo);
        
        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        if ( r_data.size() != Rx.size() ) {
            r_data.resize( Rx.size() );
        }
        
        for (size_t nr=0; nr<Rx.size(); ++nr) {
            
            traveltimes[nr]->resize( Rx[nr]->size() );
            r_data[nr]->resize( Rx[nr]->size() );
            for ( size_t ni=0; ni<r_data[nr]->size(); ++ni ) {
                (*r_data[nr])[ni].resize( 0 );
            }
            
            T2 nodeParentRx;
            T2 cellParentRx;
            
            for (size_t n=0; n<Rx[nr]->size(); ++n) {
                
                (*traveltimes[nr])[n] = this->getTraveltime((*Rx[nr])[n], this->nodes,
                                                            nodeParentRx, cellParentRx,
                                                            threadNo);
                
                bool flag=false;
                for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                    if ( (*Rx[nr])[n] == Tx[ns] ) {
                        
                        (*r_data[nr])[n].resize( 1 );
                        (*r_data[nr])[n][0] = (*Rx[nr])[n];
                        
                        flag = true;
                        break;
                    }
                }
                if ( flag ) continue;
                
                // Rx are in nodes (not txNodes)
                std::vector<Node3Dcsp<T1,T2>> *node_p;
                node_p = &(this->nodes);
                
                std::vector<sxyz<T1>> r_tmp;
                T2 iChild, iParent = nodeParentRx;
                sxyz<T1> child;
                
                // store the son's coord
                child = (*Rx[nr])[n];
                while ( (*node_p)[iParent].getNodeParent(threadNo) !=
                       std::numeric_limits<T2>::max() ) {
                    
                    r_tmp.push_back( child );
                    
                    // we now go up in time - parent becomes the child of grand'pa
                    iChild = iParent;
                    child = (*node_p)[iChild];
                    
                    // grand'pa is now papa
                    iParent = (*node_p)[iChild].getNodeParent(threadNo);
                    if ( iParent >= this->nodes.size() ) {
                        node_p = &txNodes;
                        iParent -= this->nodes.size();
                    }
                    else {
                        node_p = &(this->nodes);
                    }
                }
                
                // parent is now at Tx
                r_tmp.push_back( child );
                
                // finally, store Tx position
                child = (*node_p)[iParent];
                r_tmp.push_back( child );
                
                // the order should be from Tx to Rx, so we reorder...
                iParent = static_cast<T2>(r_tmp.size());
                (*r_data[nr])[n].resize( r_tmp.size() );
                for ( size_t nn=0; nn<(*r_data[nr])[n].size(); ++nn ) {
                    (*r_data[nr])[n][nn] = r_tmp[ iParent-1-nn ];
                }
            }
        }
    }
    
    template<typename T1, typename T2, typename CELL>
    void Grid3Drcsp<T1,T2,CELL>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                          const std::vector<T1>& t0,
                                          const std::vector<sxyz<T1>>& Rx,
                                          std::vector<T1>& traveltimes,
                                          std::vector<std::vector<sxyz<T1>>>& r_data,
                                          std::vector<std::vector<siv<T1>>>& l_data,
                                          const size_t threadNo) const {
        
        // Primary function
        
        this->checkPts(Tx);
        this->checkPts(Rx);
        
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }
        
        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node3Dcsp<T1,T2>*, std::vector<Node3Dcsp<T1,T2>*>,
        CompareNodePtr<T1>> queue(cmp);
        // txNodes: Extra nodes if the sources points are not on an existing node
        std::vector<Node3Dcsp<T1,T2>> txNodes;
        // inQueue lists the nodes waiting in the queue
        std::vector<bool> inQueue( this->nodes.size(), false );
        // Tx sources nodes are "frozen" and their traveltime can't be modified
        std::vector<bool> frozen( this->nodes.size(), false );
        
        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
        
        propagate(queue, inQueue, frozen, threadNo);
        
        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        if ( r_data.size() != Rx.size() ) {
            r_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<r_data.size(); ++ni ) {
            r_data[ni].resize( 0 );
        }
        if ( l_data.size() != Rx.size() ) {
            l_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<l_data.size(); ++ni ) {
            l_data[ni].resize( 0 );
        }
        
        T2 nodeParentRx;
        T2 cellParentRx;
        
        for (size_t n=0; n<Rx.size(); ++n) {
            traveltimes[n] = this->getTraveltime(Rx[n], this->nodes, nodeParentRx, cellParentRx,
                                                 threadNo);
            
            // Rx are in nodes (not txNodes)
            std::vector<Node3Dcsp<T1,T2>> *node_p;
            node_p = &(this->nodes);
            
            std::vector<sxyz<T1>> r_tmp;
            T2 iChild, iParent = nodeParentRx;
            sxyz<T1> child;
            siv<T1> cell;
            
            // store the son's coord
            child.x = Rx[n].x;
            child.y = Rx[n].y;
            child.z = Rx[n].z;
            cell.i = cellParentRx;
            while ( (*node_p)[iParent].getNodeParent(threadNo) != std::numeric_limits<T2>::max() ) {
                
                r_tmp.push_back( child );
                
                cell.v = (*node_p)[iParent].getDistance( child );
                bool found=false;
                for (size_t nc=0; nc<l_data[n].size(); ++nc) {
                    if ( l_data[n][nc].i == cell.i ) {
                        l_data[n][nc].v += cell.v;  // must add in case we pass through secondary nodes along edge
                        found = true;
                        break;
                    }
                }
                if ( found == false ) {
                    l_data[n].push_back( cell );
                }
                
                // we now go up in time - parent becomes the child of grand'pa
                iChild = iParent;
                child.x = (*node_p)[iChild].getX();
                child.y = (*node_p)[iChild].getY();
                child.z = (*node_p)[iChild].getZ();
                cell.i = (*node_p)[iChild].getCellParent(threadNo);
                
                // grand'pa is now papa
                iParent = (*node_p)[iChild].getNodeParent(threadNo);
                if ( iParent >= this->nodes.size() ) {
                    node_p = &txNodes;
                    iParent -= this->nodes.size();
                }
                else {
                    node_p = &(this->nodes);
                }
            }
            
            // parent is now at Tx
            r_tmp.push_back( child );
            
            cell.v = (*node_p)[iParent].getDistance( child );
            bool found=false;
            for (size_t nc=0; nc<l_data[n].size(); ++nc) {
                if ( l_data[n][nc].i == cell.i ) {
                    l_data[n][nc].v += cell.v;  // must add in case we pass through secondary nodes along edge
                    found = true;
                    break;
                }
            }
            if ( found == false ) {
                l_data[n].push_back( cell );
            }
            
            // finally, store Tx position
            child.x = (*node_p)[iParent].getX();
            child.y = (*node_p)[iParent].getY();
            child.z = (*node_p)[iParent].getZ();
            r_tmp.push_back( child );
            
            //  must be sorted to build matrix L
            sort(l_data[n].begin(), l_data[n].end(), CompareSiv_i<T1>());
            
            // the order should be from Tx to Rx, so we reorder...
            iParent = static_cast<T2>(r_tmp.size());
            r_data[n].resize( r_tmp.size() );
            for ( size_t nn=0; nn<r_data[n].size(); ++nn ) {
                r_data[n][nn].x = r_tmp[ iParent-1-nn ].x;
                r_data[n][nn].y = r_tmp[ iParent-1-nn ].y;
                r_data[n][nn].z = r_tmp[ iParent-1-nn ].z;
            }
        }
    }
    
    template<typename T1, typename T2, typename CELL>
    void Grid3Drcsp<T1,T2,CELL>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                          const std::vector<T1>& t0,
                                          const std::vector<sxyz<T1>>& Rx,
                                          std::vector<T1>& traveltimes,
                                          std::vector<std::vector<siv<T1>>>& l_data,
                                          const size_t threadNo) const {

        // Primary function

        this->checkPts(Tx);
        this->checkPts(Rx);

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node3Dcsp<T1,T2>*, std::vector<Node3Dcsp<T1,T2>*>,
        CompareNodePtr<T1>> queue(cmp);
        // txNodes: Extra nodes if the sources points are not on an existing node
        std::vector<Node3Dcsp<T1,T2>> txNodes;
        // inQueue lists the nodes waiting in the queue
        std::vector<bool> inQueue( this->nodes.size(), false );
        // Tx sources nodes are "frozen" and their traveltime can't be modified
        std::vector<bool> frozen( this->nodes.size(), false );

        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);

        propagate(queue, inQueue, frozen, threadNo);

        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        if ( l_data.size() != Rx.size() ) {
            l_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<l_data.size(); ++ni ) {
            l_data[ni].resize( 0 );
        }

        T2 nodeParentRx;
        T2 cellParentRx;

        for (size_t n=0; n<Rx.size(); ++n) {
            traveltimes[n] = this->getTraveltime(Rx[n], this->nodes, nodeParentRx, cellParentRx,
                                                 threadNo);

            // Rx are in nodes (not txNodes)
            std::vector<Node3Dcsp<T1,T2>> *node_p;
            node_p = &(this->nodes);

            std::vector<sxyz<T1>> r_tmp;
            T2 iChild, iParent = nodeParentRx;
            sxyz<T1> child;
            siv<T1> cell;

            // store the son's coord
            child.x = Rx[n].x;
            child.y = Rx[n].y;
            child.z = Rx[n].z;
            cell.i = cellParentRx;
            while ( (*node_p)[iParent].getNodeParent(threadNo) != std::numeric_limits<T2>::max() ) {

                r_tmp.push_back( child );

                cell.v = (*node_p)[iParent].getDistance( child );
                bool found=false;
                for (size_t nc=0; nc<l_data[n].size(); ++nc) {
                    if ( l_data[n][nc].i == cell.i ) {
                        l_data[n][nc].v += cell.v;  // must add in case we pass through secondary nodes along edge
                        found = true;
                        break;
                    }
                }
                if ( found == false ) {
                    l_data[n].push_back( cell );
                }

                // we now go up in time - parent becomes the child of grand'pa
                iChild = iParent;
                child.x = (*node_p)[iChild].getX();
                child.y = (*node_p)[iChild].getY();
                child.z = (*node_p)[iChild].getZ();
                cell.i = (*node_p)[iChild].getCellParent(threadNo);

                // grand'pa is now papa
                iParent = (*node_p)[iChild].getNodeParent(threadNo);
                if ( iParent >= this->nodes.size() ) {
                    node_p = &txNodes;
                    iParent -= this->nodes.size();
                }
                else {
                    node_p = &(this->nodes);
                }
            }

            // parent is now at Tx
            r_tmp.push_back( child );

            cell.v = (*node_p)[iParent].getDistance( child );
            bool found=false;
            for (size_t nc=0; nc<l_data[n].size(); ++nc) {
                if ( l_data[n][nc].i == cell.i ) {
                    l_data[n][nc].v += cell.v;  // must add in case we pass through secondary nodes along edge
                    found = true;
                    break;
                }
            }
            if ( found == false ) {
                l_data[n].push_back( cell );
            }

            // finally, store Tx position
            child.x = (*node_p)[iParent].getX();
            child.y = (*node_p)[iParent].getY();
            child.z = (*node_p)[iParent].getZ();
            r_tmp.push_back( child );

            //  must be sorted to build matrix L
            sort(l_data[n].begin(), l_data[n].end(), CompareSiv_i<T1>());
        }
    }

    template<typename T1, typename T2, typename CELL>
    void Grid3Drcsp<T1,T2,CELL>::raytrace2(const std::vector<sxyz<T1>>& Tx,
                                           const std::vector<T1>& t0,
                                           const std::vector<sxyz<T1>>& Rx,
                                           std::vector<T1>& traveltimes,
                                           const size_t threadNo) const {
        
        // Primary function
        
        this->checkPts(Tx);
        this->checkPts(Rx);
        
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }
        
        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node3Dcsp<T1,T2>*, std::deque<Node3Dcsp<T1,T2>*>,
        CompareNodePtr<T1>> queue(cmp);
        // txNodes: Extra nodes if the sources points are not on an existing node
        std::vector<Node3Dcsp<T1,T2>> txNodes;
        // inQueue lists the nodes waiting in the queue
        std::vector<bool> inQueue( this->nodes.size(), false );
        // Tx sources nodes are "frozen" and their traveltime can't be modified
        std::vector<bool> frozen( this->nodes.size(), false );
        
        initQueue2(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
        
        propagate2(queue, inQueue, frozen, threadNo);
        
        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        
        for (size_t n=0; n<Rx.size(); ++n) {
            traveltimes[n] = this->getTraveltime(Rx[n], this->nodes, threadNo);
        }
    }    
}

#endif

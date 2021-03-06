//
//  Grid3Drnsp.h
//  ttcr
//
//  Created by Bernard Giroux on 2015-07-03.
//  Copyright (c) 2015 Bernard Giroux. All rights reserved.
//

#ifndef ttcr_Grid3Drnsp_h
#define ttcr_Grid3Drnsp_h

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <ostream>
#include <map>
#include <queue>
#include <vector>

#include "Grid3Drn.h"
#include "Node3Dnsp.h"
#include "utils.h"

#include "Interpolator.h"

namespace ttcr {
    
    template<typename T1, typename T2>
    class Grid3Drnsp : public Grid3Drn<T1,T2,Node3Dnsp<T1,T2>> {
    public:
        Grid3Drnsp(const T2 nx, const T2 ny, const T2 nz,
                   const T1 ddx, const T1 ddy, const T1 ddz,
                   const T1 minx, const T1 miny, const T1 minz,
                   const T2 nnx, const T2 nny, const T2 nnz, const bool ttrp,
                   const bool intVel, const size_t nt=1) :
        Grid3Drn<T1,T2,Node3Dnsp<T1,T2>>(nx, ny, nz, ddx, ddy, ddz, minx, miny, minz, ttrp, intVel, nt),
        nsnx(nnx), nsny(nny), nsnz(nnz)
        {
            buildGridNodes();
            this->template buildGridNeighbors<Node3Dnsp<T1,T2>>(this->nodes);
        }
        
        ~Grid3Drnsp() {
            
        }
        
        void setSlowness(const std::vector<T1>& s);
        
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
        
        void savePrimary(const char filename[], const size_t nt=0,
                         const bool vtkFormat=0) const;
        
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
                       std::priority_queue<Node3Dnsp<T1,T2>*,
                       std::vector<Node3Dnsp<T1,T2>*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<Node3Dnsp<T1,T2>>& txNodes,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       const size_t threadNo) const;
        
        void propagate(std::priority_queue<Node3Dnsp<T1,T2>*, std::vector<Node3Dnsp<T1,T2>*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       size_t threadNo) const;
        
        void prepropagate(const Node3Dnsp<T1,T2>& node,
                          std::priority_queue<Node3Dnsp<T1,T2>*, std::vector<Node3Dnsp<T1,T2>*>,
                          CompareNodePtr<T1>>& queue,
                          std::vector<bool>& inQueue,
                          std::vector<bool>& frozen,
                          size_t threadNo) const;
        
    };
    
    
    template<typename T1, typename T2>
    void Grid3Drnsp<T1,T2>::buildGridNodes() {
        
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
                           Node3Dnsp<T1,T2>(this->nThreads));
        
        // Create the grid, assign a number for each node, determine the type of the node and find the owners
        // Nodes and cells are first indexed in z, then y, and x.
        // Secondary nodes are placed on the faces and edges of every cells.
        // Ex: the node in "node[A]=(i,j,k)" is followed by the node in "node[A+1]=(i+dx,j,k)"
        
        T1 dxs = this->dx/(this->nsnx+1); 	// distance between secondary nodes in x
        T1 dys = this->dy/(this->nsny+1);
        T1 dzs = this->dz/(this->nsnz+1);
        
        T2 cell_XmYmZm; 	// cell in the (x-,y-,z-) direction from the node
        T2 cell_XpYmZm; 	// cell in the (x+,y-,z-) direction from the node
        T2 cell_XmYpZm;
        T2 cell_XpYpZm;
        T2 cell_XmYmZp;
        T2 cell_XpYmZp;
        T2 cell_XmYpZp;
        T2 cell_XpYpZp;
        
        T2 n=0;
        for ( T2 nk=0; nk<=this->ncz; ++nk ) {
            
            T1 z = this->zmin + nk*this->dz;
            
            for ( T2 nj=0; nj<=this->ncy; ++nj ) {
                
                T1 y = this->ymin + nj*this->dy;
                
                for ( T2 ni=0; ni<=this->ncx; ++ni, ++n ){
                    
                    T1 x = this->xmin + ni*this->dx;
                    
                    // Find the adjacent cells for each primary node
                    
                    if (ni < this->ncx && nj < this->ncy && nk < this->ncz){
                        cell_XpYpZp = nj*this->ncx + nk*(this->ncx*this->ncy) + ni;
                    }
                    else {
                        cell_XpYpZp = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni > 0 && nj < this->ncy && nk < this->ncz){
                        cell_XmYpZp = nj*this->ncx + nk*(this->ncx*this->ncy) + ni - 1;
                    }
                    else {
                        cell_XmYpZp = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni < this->ncx && nj > 0 && nk < this->ncz){
                        cell_XpYmZp = (nj-1)*this->ncx + nk*(this->ncx*this->ncy) + ni;
                    }
                    else {
                        cell_XpYmZp = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni > 0 && nj > 0 && nk < this->ncz){
                        cell_XmYmZp = (nj-1)*this->ncx + nk*(this->ncx * this->ncy) + ni - 1;
                    }
                    else {
                        cell_XmYmZp = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni < this->ncx && nj < this->ncy && nk > 0){
                        cell_XpYpZm = nj*this->ncx + (nk-1)*(this->ncx*this->ncy) + ni;
                    }
                    else {
                        cell_XpYpZm = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni > 0 && nj < this->ncy && nk > 0){
                        cell_XmYpZm = nj*this->ncx + (nk-1)*(this->ncx*this->ncy) + ni - 1;
                    }
                    else {
                        cell_XmYpZm = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni < this->ncx && nj > 0 && nk > 0){
                        cell_XpYmZm = (nj-1)*this->ncx + (nk-1)*(this->ncx*this->ncy) + ni;
                    }
                    else {
                        cell_XpYmZm = std::numeric_limits<T2>::max();
                    }
                    
                    if (ni > 0 && nj > 0 && nk > 0){
                        cell_XmYmZm = (nj-1)*this->ncx + (nk-1)*(this->ncx*this->ncy) + ni - 1;
                    }
                    else {
                        cell_XmYmZm = std::numeric_limits<T2>::max();
                    }
                    
                    
                    // Index the primary nodes owners
                    
                    if (cell_XmYmZm != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cell_XmYmZm );
                    }
                    if (cell_XpYmZm != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cell_XpYmZm );
                    }
                    if (cell_XmYpZm != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cell_XmYpZm );
                    }
                    if (cell_XpYpZm != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cell_XpYpZm );
                    }
                    if (cell_XmYmZp != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cell_XmYmZp );
                    }
                    if (cell_XpYmZp != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cell_XpYmZp );
                    }
                    if (cell_XmYpZp != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cell_XmYpZp );
                    }
                    if (cell_XpYpZp != std::numeric_limits<T2>::max() ) {
                        this->nodes[n].pushOwner( cell_XpYpZp );
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
                        cell_XpYpZp = nj*this->ncx + nk*(this->ncx*this->ncy) + ni;
                    }
                    else {
                        cell_XpYpZp = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj < this->ncy && nk < this->ncz){
                        cell_XmYpZp = nj*this->ncx + nk*(this->ncx*this->ncy) + ni - 1;
                    }
                    else {
                        cell_XmYpZp = std::numeric_limits<T2>::max();
                    }

                    if (ni < this->ncx && nj > 0 && nk < this->ncz){
                        cell_XpYmZp = (nj-1)*this->ncx + nk*(this->ncx*this->ncy) + ni;
                    }
                    else {
                        cell_XpYmZp = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj > 0 && nk < this->ncz){
                        cell_XmYmZp = (nj-1)*this->ncx + nk*(this->ncx * this->ncy) + ni - 1;
                    }
                    else {
                        cell_XmYmZp = std::numeric_limits<T2>::max();
                    }

                    if (ni < this->ncx && nj < this->ncy && nk > 0){
                        cell_XpYpZm = nj*this->ncx + (nk-1)*(this->ncx*this->ncy) + ni;
                    }
                    else {
                        cell_XpYpZm = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj < this->ncy && nk > 0){
                        cell_XmYpZm = nj*this->ncx + (nk-1)*(this->ncx*this->ncy) + ni - 1;
                    }
                    else {
                        cell_XmYpZm = std::numeric_limits<T2>::max();
                    }

                    if (ni < this->ncx && nj > 0 && nk > 0){
                        cell_XpYmZm = (nj-1)*this->ncx + (nk-1)*(this->ncx*this->ncy) + ni;
                    }
                    else {
                        cell_XpYmZm = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj > 0 && nk > 0){
                        cell_XmYmZm = (nj-1)*this->ncx + (nk-1)*(this->ncx*this->ncy) + ni - 1;
                    }
                    else {
                        cell_XmYmZm = std::numeric_limits<T2>::max();
                    }

                    // Secondary nodes on x edge
                    if ( ni < this->ncx ) {
                        for (T2 ns=0; ns< this->nsnx; ++ns, ++n ) {
                            
                            T1 xsv = this->xmin + ni*this->dx + (ns+1)*dxs;
                            
                            if ( cell_XpYmZm != std::numeric_limits<T2>::max() ) {
                                this->nodes[n].pushOwner( cell_XpYmZm );
                            }
                            if ( cell_XpYpZm != std::numeric_limits<T2>::max() ) {
                                this->nodes[n].pushOwner( cell_XpYpZm );
                            }
                            if ( cell_XpYmZp != std::numeric_limits<T2>::max() ) {
                                this->nodes[n].pushOwner( cell_XpYmZp );
                            }
                            if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
                                this->nodes[n].pushOwner( cell_XpYpZp );
                            }
                            this->nodes[n].setXYZindex( xsv, y, z, n );
                            
                            if (nj >0 && nk>0){
                                this->nodes[n].setPrimary(false);
                            }
                            else if (nj==0 && nk>0){
                                this->nodes[n].setPrimary(false);
                            }
                            else if (nj>0 && nk==0){
                                this->nodes[n].setPrimary(false);
                            }
                            else if (nj==0 && nk==0){
                                this->nodes[n].setPrimary(false);
                            }
                        }
                    }
                    
                    // Secondary nodes on y edge
                    if ( nj < this->ncy ) {
                        for (T2 ns=0; ns< this->nsny; ++ns, ++n ) {
                            
                            T1 ysv = this->ymin + nj* this->dy + (ns+1)*dys;
                            
                            if ( cell_XmYpZm != std::numeric_limits<T2>::max() ) {
                                this->nodes[n].pushOwner( cell_XmYpZm );
                            }
                            if ( cell_XpYpZm != std::numeric_limits<T2>::max() ) {
                                this->nodes[n].pushOwner( cell_XpYpZm );
                            }
                            if ( cell_XmYpZp != std::numeric_limits<T2>::max() ) {
                                this->nodes[n].pushOwner( cell_XmYpZp );
                            }
                            if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
                                this->nodes[n].pushOwner( cell_XpYpZp );
                            }
                            this->nodes[n].setXYZindex( x, ysv, z, n );
                            
                            if (ni >0 && nk>0){
                                this->nodes[n].setPrimary(false);
                            }
                            else if (ni>0 && nk==0){
                                this->nodes[n].setPrimary(false);
                            }
                            else if (ni==0 && nk>0){
                                this->nodes[n].setPrimary(false);
                            }
                            else if (ni==0 && nk==0){
                                this->nodes[n].setPrimary(false);
                            }
                        }
                    }
                    
                    // Secondary nodes on z edge
                    if ( nk < this->ncz ) {
                        for (T2 ns=0; ns< this->nsnz; ++ns, ++n ) {
                            
                            T1 zsv = this->zmin + nk* this->dz + (ns+1)*dzs;
                            
                            if ( cell_XmYmZp != std::numeric_limits<T2>::max() ) {
                                this->nodes[n].pushOwner( cell_XmYmZp );
                            }
                            if ( cell_XpYmZp != std::numeric_limits<T2>::max() ) {
                                this->nodes[n].pushOwner( cell_XpYmZp );
                            }
                            if ( cell_XmYpZp != std::numeric_limits<T2>::max() ) {
                                this->nodes[n].pushOwner( cell_XmYpZp );
                            }
                            if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
                                this->nodes[n].pushOwner( cell_XpYpZp );
                            }
                            this->nodes[n].setXYZindex( x, y, zsv, n );
                            
                            if (ni >0 && nj>0){
                                this->nodes[n].setPrimary(false);
                            }
                            else if (ni>0 && nj==0){
                                this->nodes[n].setPrimary(false);
                            }
                            else if (ni==0 && nj>0){
                                this->nodes[n].setPrimary(false);
                            }
                            else if (ni==0 && nj==0){
                                this->nodes[n].setPrimary(false);
                            }
                        }
                    }
                    
                    // Secondary nodes on the xy0 planes
                    if ( ni < this->ncx && nj < this->ncy ) {
                        for(T2 sy=0; sy < this->nsny; ++sy){
                            for(T2 sx=0; sx < this->nsnx; ++sx, n++){
                                
                                T1 ysv= this->ymin+ nj* this->dy+ (sy+1)*dys;
                                T1 xsv= this->xmin+ ni* this->dx+ (sx+1)*dxs;
                                
                                if ( cell_XpYpZm != std::numeric_limits<T2>::max() ) {
                                    this->nodes[n].pushOwner( cell_XpYpZm );
                                }
                                if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
                                    this->nodes[n].pushOwner( cell_XpYpZp );
                                }
                                this->nodes[n].setXYZindex( xsv, ysv, z, n );
                                
                                if (nk>0){
                                    this->nodes[n].setPrimary(false);
                                }
                                else if (nk==0){
                                    this->nodes[n].setPrimary(false);
                                }
                            }
                        }
                    }
                    
                    // Secondary nodes on the x0z planes
                    if ( ni < this->ncx && nk < this->ncz ) {
                        for(T2 sz=0; sz < this->nsnz; ++sz){
                            for(T2 sx=0; sx < this->nsnx; ++sx, n++){
                                
                                T1 zsv= this->zmin+ nk* this->dz+ (sz+1)*dzs;
                                T1 xsv= this->xmin+ ni* this->dx+ (sx+1)*dxs;
                                
                                if ( cell_XpYmZp != std::numeric_limits<T2>::max() ) {
                                    this->nodes[n].pushOwner( cell_XpYmZp );
                                }
                                if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
                                    this->nodes[n].pushOwner( cell_XpYpZp );
                                }
                                this->nodes[n].setXYZindex( xsv, y, zsv, n );
                                
                                if (nj>0){
                                    this->nodes[n].setPrimary(false);
                                }
                                else if (nj==0){
                                    this->nodes[n].setPrimary(false);
                                }
                            }
                        }
                    }
                    
                    // Secondary nodes on the 0yz planes
                    if ( nj < this->ncy && nk < this->ncz ) {
                        for(T2 sz=0; sz < this->nsnz; ++sz){
                            for(T2 sy=0; sy < this->nsny; ++sy, n++){
                                
                                T1 zsv= this->zmin + nk*this->dz + (sz+1)*dzs;
                                T1 ysv= this->ymin + nj*this->dy + (sy+1)*dys;
                                
                                if ( cell_XmYpZp != std::numeric_limits<T2>::max() ) {
                                    this->nodes[n].pushOwner( cell_XmYpZp );
                                }
                                if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
                                    this->nodes[n].pushOwner( cell_XpYpZp );
                                }
                                this->nodes[n].setXYZindex( x, ysv, zsv, n );
                                
                                if (ni>0){
                                    this->nodes[n].setPrimary(false);
                                }
                                else if (ni==0){
                                    this->nodes[n].setPrimary(false);
                                }
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
    
    template<typename T1, typename T2>
    void Grid3Drnsp<T1,T2>::setSlowness(const std::vector<T1>& s) {
        
        if ( ((this->ncx+1)*(this->ncy+1)*(this->ncz+1)) != s.size() ) {
            throw std::length_error("Error: slowness vector of incompatible size.");
        }
        //Set the slowness for primary nodes
        size_t i=0;
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            if (this->nodes[n].isPrimary()){
                this->nodes[n].setNodeSlowness( s[i] );
                i++;
            }
        }
        this->interpSecondary();

    }

    
    template<typename T1, typename T2>
    void Grid3Drnsp<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
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
        std::priority_queue< Node3Dnsp<T1,T2>*, std::vector<Node3Dnsp<T1,T2>*>,
        CompareNodePtr<T1>> queue(cmp);
        // txNodes: Extra nodes if the sources points are not on an existing node
        std::vector<Node3Dnsp<T1,T2>> txNodes;
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
//            if ( verbose ) {
//                std::cout << "done.\n  Updating traveltimes from raypaths ... ";
//            }
//            for (size_t n=0; n<Rx.size(); ++n) {
//                if ( verbose > 1 ) {
//                    std::cout << "\n    Rx no " << n << std::flush;
//                }
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
                traveltimes[n] = this->getTraveltime(Rx[n], threadNo);
            }
//        }
    }
    
    template<typename T1, typename T2>
    void Grid3Drnsp<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
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
        std::priority_queue< Node3Dnsp<T1,T2>*, std::vector<Node3Dnsp<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );
        
        std::vector<Node3Dnsp<T1,T2>> txNodes;
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
                    
                    (*traveltimes[nr])[n] = this->getTraveltime((*Rx[nr])[n], threadNo);
                }
            }
//        }
    }
    
    template<typename T1, typename T2>
    void Grid3Drnsp<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
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
        std::priority_queue< Node3Dnsp<T1,T2>*, std::vector<Node3Dnsp<T1,T2>*>,
        CompareNodePtr<T1>> queue(cmp);
        // txNodes: Extra nodes if the sources points are not on an existing node
        std::vector<Node3Dnsp<T1,T2>> txNodes;
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
            traveltimes[n] = this->getTraveltime(Rx[n], nodeParentRx, cellParentRx,
                                                 threadNo);
            
            // Rx are in nodes (not txNodes)
            std::vector<Node3Dnsp<T1,T2>> *node_p;
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
    
    template<typename T1, typename T2>
    void Grid3Drnsp<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
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
        std::priority_queue< Node3Dnsp<T1,T2>*, std::vector<Node3Dnsp<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );
        
        std::vector<Node3Dnsp<T1,T2>> txNodes;
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
                bool onTx = false;
                for (size_t nt=0; nt<Tx.size(); ++nt) {
                    if ( Tx[nt] == (*Rx[nr])[n] ) {
                        onTx = true;
                        (*traveltimes[nr])[n] = t0[nt];
                        (*r_data[nr])[n].push_back(Tx[nt]);
                        break;
                    }
                }
                if (onTx) { continue; }
                
                (*traveltimes[nr])[n] = this->getTraveltime((*Rx[nr])[n],
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
                std::vector<Node3Dnsp<T1,T2>> *node_p;
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
    
    template<typename T1, typename T2>
    void Grid3Drnsp<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
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
        std::priority_queue< Node3Dnsp<T1,T2>*, std::vector<Node3Dnsp<T1,T2>*>,
        CompareNodePtr<T1>> queue(cmp);
        // txNodes: Extra nodes if the sources points are not on an existing node
        std::vector<Node3Dnsp<T1,T2>> txNodes;
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
            traveltimes[n] = this->getTraveltime(Rx[n], nodeParentRx, cellParentRx,
                                                 threadNo);
            
            // Rx are in nodes (not txNodes)
            std::vector<Node3Dnsp<T1,T2>> *node_p;
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
    
    template<typename T1, typename T2>
    void Grid3Drnsp<T1,T2>::initQueue(const std::vector<sxyz<T1>>& Tx,
                                      const std::vector<T1>& t0,
                                      std::priority_queue<Node3Dnsp<T1,T2>*,
                                      std::vector<Node3Dnsp<T1,T2>*>,
                                      CompareNodePtr<T1>>& queue,
                                      std::vector<Node3Dnsp<T1,T2>>& txNodes,
                                      std::vector<bool>& inQueue,
                                      std::vector<bool>& frozen,
                                      const size_t threadNo) const {
        
        //Find the starting nodes of the transmitters Tx and start the queue list
        for (size_t n=0; n<Tx.size(); ++n){
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
                txNodes.push_back( Node3Dnsp<T1,T2>(t0[n], Tx[n].x, Tx[n].y, Tx[n].z, this->nThreads, threadNo));
                txNodes.back().pushOwner( this->getCellNo(Tx[n]) );
                txNodes.back().setGridIndex( static_cast<T2>(this->nodes.size()+txNodes.size()-1) );
                frozen.push_back( true );
                T1 slo = this->computeSlowness( Tx[n] );
                txNodes.back().setNodeSlowness( slo );
                
                prepropagate(txNodes.back(), queue, inQueue, frozen, threadNo); // See description in the function declaration
                
                //	queue.push( &(txNodes.back()) );	//Don't use if prepropagate is used
                //	inQueue.push_back( true );			//Don't use if prepropagate is used
                
            }
        }
    }
    
    
    template<typename T1, typename T2>
    void Grid3Drnsp<T1,T2>::propagate(std::priority_queue<Node3Dnsp<T1,T2>*,
                                      std::vector<Node3Dnsp<T1,T2>*>,
                                      CompareNodePtr<T1>>& queue,
                                      std::vector<bool>& inQueue,
                                      std::vector<bool>& frozen,
                                      size_t threadNo) const {
        
        while ( !queue.empty() ) {
            const Node3Dnsp<T1,T2>* source = queue.top();
            queue.pop();
            inQueue[ source->getGridIndex() ] = false;
            frozen[ source->getGridIndex() ] = true;
            
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
                        T1 dt = this->computeDt(*source, this->nodes[neibNo]);
                        
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
    
    template<typename T1, typename T2>
    void Grid3Drnsp<T1,T2>::prepropagate(const Node3Dnsp<T1,T2>& node,
                                         std::priority_queue<Node3Dnsp<T1,T2>*,
                                         std::vector<Node3Dnsp<T1,T2>*>,
                                         CompareNodePtr<T1>>& queue,
                                         std::vector<bool>& inQueue,
                                         std::vector<bool>& frozen,
                                         const size_t threadNo) const {
        
        // This function can be used to "prepropagate" each Tx nodes one first time
        // during "initQueue", before running "propagate".
        // When a Tx source node seems to be lost in the queue and is not
        // propagated, corrupting the entire traveltime table,
        // this function force the propagation of every source points and can
        // solve the problem.
        
        for ( size_t no=0; no<node.getOwners().size(); ++no ) {
            T2 cellNo = node.getOwners()[no];
            for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                T2 neibNo = this->neighbors[cellNo][k];
                if ( neibNo == node.getGridIndex() || frozen[neibNo] ) {
                    continue;
                }
                
                // compute dt
                T1 dt = this->computeDt(node, this->nodes[neibNo]);
                
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
    
    
    template<typename T1, typename T2>
    void Grid3Drnsp<T1,T2>::savePrimary(const char filename[], const size_t nt,
                                        const bool vtkFormat) const {
        
        if ( vtkFormat ) {
            
#ifdef VTK
            vtkSmartPointer<vtkDoubleArray> xCoords = vtkSmartPointer<vtkDoubleArray>::New();
            for (size_t ni=0; ni<=this->ncx; ++ni)
                xCoords->InsertNextValue(this->xmin + ni*this->dx);
            vtkSmartPointer<vtkDoubleArray> yCoords = vtkSmartPointer<vtkDoubleArray>::New();
            for (size_t nj=0; nj<=this->ncy; ++nj)
                yCoords->InsertNextValue(this->ymin + nj*this->dy);
            vtkSmartPointer<vtkDoubleArray> zCoords = vtkSmartPointer<vtkDoubleArray>::New();
            for (size_t nk=0; nk<=this->ncz; ++nk)
                zCoords->InsertNextValue(this->zmin + nk*this->dz);
            
            vtkSmartPointer<vtkRectilinearGrid> rgrid = vtkSmartPointer<vtkRectilinearGrid>::New();
            rgrid->SetDimensions(this->ncx, this->ncy, this->ncz);
            rgrid->SetXCoordinates(xCoords);
            rgrid->SetYCoordinates(yCoords);
            rgrid->SetZCoordinates(zCoords);
            
            vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
            data->SetName("Travel time");
            size_t n=0;
            for ( size_t nk=0; nk<=this->ncz; ++nk ) {
                for ( size_t nj=0; nj<=this->ncy; ++nj ) {
                    for ( size_t ni=0; ni<=this->ncx; ++ni ) {
                        
                        data->InsertNextValue( this->nodes[n++].getTT(nt) );
                        
                        // Secondary nodes on x edge
                        if ( ni < this->ncx ) {
                            n += this->nsnx;
                        }
                        
                        // Secondary nodes on y edge
                        if ( nj < this->ncy ) {
                            n += this->nsny;
                        }
                        
                        // Secondary nodes on z edge
                        if ( nk < this->ncz ) {
                            n += this->nsnz;
                        }
                        
                        // Secondary nodes on the xy0 planes
                        if ( ni < this->ncx && nj < this->ncy ) {
                            n += this->nsny*this->nsnx;
                        }
                        
                        // Secondary nodes on the x0z planes
                        if ( ni < this->ncx && nk < this->ncz ) {
                            n += this->nsnz*this->nsnx;
                        }
                        
                        // Secondary nodes on the 0yz planes
                        if ( nj < this->ncy && nk < this->ncz ) {
                            n += this->nsnz*this->nsny;
                        }
                    }
                }
            }
            rgrid->GetPointData()->SetScalars( data );
            
            vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
            
            writer->SetFileName( filename );
            //			writer->SetInputConnection( rgrid->GetProducerPort() );
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
            for ( size_t nk=0; nk<=this->ncz; ++nk ) {
                
                for ( size_t nj=0; nj<=this->ncy; ++nj ) {
                    
                    for (size_t ni=0; ni<=this->ncx; ++ni ) {
                        
                        fout << this->nodes[n++].getTT(nt) << '\n';
                        
                        // Secondary nodes on x edge
                        if ( ni < this->ncx ) {
                            n += this->nsnx;
                        }
                        
                        // Secondary nodes on y edge
                        if ( nj < this->ncy ) {
                            n += this->nsny;
                        }
                        
                        // Secondary nodes on z edge
                        if ( nk < this->ncz ) {
                            n += this->nsnz;
                        }
                        
                        // Secondary nodes on the xy0 planes
                        if ( ni < this->ncx && nj < this->ncy ) {
                            n += this->nsny*this->nsnx;
                        }
                        
                        // Secondary nodes on the x0z planes
                        if ( ni < this->ncx && nk < this->ncz ) {
                            n += this->nsnz*this->nsnx;
                        }
                        
                        // Secondary nodes on the 0yz planes
                        if ( nj < this->ncy && nk < this->ncz ) {
                            n += this->nsnz*this->nsny;
                        }
                    }
                }
            }
            fout.close();
        }
    }
    
}

#endif

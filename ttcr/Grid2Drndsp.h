//
//  Grid2Drndsp.h
//  ttcr
//
//  Created by Bernard Giroux on 2021-02-27.
//  Copyright © 2021 Bernard Giroux. All rights reserved.
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

#ifndef ttcr_Grid2Drndsp_h
#define ttcr_Grid2Drndsp_h

#include "Grid2Drn.h"
#include "Node2Dn.h"
#include "Node2Dnd.h"

namespace ttcr {

    template<typename T1, typename T2, typename S>
    class Grid2Drndsp : public Grid2Drn<T1,T2,S,Node2Dn<T1,T2>> {
    public:
        Grid2Drndsp(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
                    const T1 minx, const T1 minz, const T2 ns, const T2 nd,
                    const T1 drad, const bool ttrp, const bool useEdgeLength=true,
                    const size_t nt=1) :
        Grid2Drn<T1,T2,S,Node2Dn<T1,T2>>(nx,nz,ddx,ddz,minx,minz,ttrp,nt),
        nSecondary(ns), nTertiary(nd), nPermanent(0),
        dynRadius(drad),
        tempNodes(std::vector<std::vector<Node2Dnd<T1,T2>>>(nt)),
        tempNeighbors(std::vector<std::vector<std::vector<T2>>>(nt))
        {
            buildGridNodes();
            this->template buildGridNeighbors<Node2Dn<T1,T2>>(this->nodes);
            nPermanent = static_cast<T2>(this->nodes.size());
            for ( size_t n=0; n<nt; ++n ) {
                tempNeighbors[n].resize(this->ncx * this->ncz);
            }
            if (useEdgeLength) {
                T1 dx = 0.5 * (ddx+ddz);
                dynRadius *= dx;
            }
        }

        ~Grid2Drndsp() {
        }

        void setSlowness(const std::vector<T1>& s) {
            T2 nPrimary = (this->ncx+1) * (this->ncz+1);
            if ( nPrimary != s.size() ) {
                throw std::length_error("Error: slowness vectors of incompatible size.");
            }
            for ( size_t n=0; n<nPrimary; ++n ) {
                this->nodes[n].setNodeSlowness( s[n] );
            }
            interpSlownessSecondary();
        }

    private:
        T2 nSecondary;                 // number of secondary nodes
        T2 nTertiary;                  // number of tertiary nodes
        T2 nPermanent;                 // total nb of primary & secondary
        T1 dynRadius;

        // we will store temporary nodes in a separate container.  This is to
        // allow threaded computations with different Tx (location of temp
        // nodes vary from one Tx to the other)
        mutable std::vector<std::vector<Node2Dnd<T1,T2>>> tempNodes;
        mutable std::vector<std::vector<std::vector<T2>>> tempNeighbors;

        void buildGridNodes();

        void interpSlownessSecondary();

        void addTemporaryNodes(const std::vector<S>&, const size_t) const;

        void initQueue(const std::vector<S>& Tx,
                       const std::vector<T1>& t0,
                       std::priority_queue<Node2Dn<T1,T2>*,
                       std::vector<Node2Dn<T1,T2>*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<Node2Dnd<T1,T2>>& txNodes,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       const size_t threadNo) const;

        void propagate(std::priority_queue<Node2Dn<T1,T2>*,
                       std::vector<Node2Dn<T1,T2>*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       const size_t threadNo) const;

        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<S>& Rx,
                      const size_t threadNo=0) const;

        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<const std::vector<S>*>& Rx,
                      const size_t threadNo=0) const;
    };

    template<typename T1, typename T2, typename S>
    void Grid2Drndsp<T1,T2,S>::buildGridNodes() {

        this->nodes.resize( // noeuds secondaires
                           this->ncx*nSecondary*(this->ncz+1) +
                           this->ncz*nSecondary*(this->ncx+1) +
                           // noeuds primaires
                           (this->ncx+1) * (this->ncz+1),
                           Node2Dn<T1,T2>(this->nThreads));

        T1 dxs = this->dx/(nSecondary+1);
        T1 dzs = this->dz/(nSecondary+1);

        T2 cell_upLeft = std::numeric_limits<T2>::max();
        T2 cell_upRight = std::numeric_limits<T2>::max();
        T2 cell_downLeft = 0;
        T2 cell_downRight = 0;

        T2 n = 0;
        // start with primary nodes
        for ( T2 nc=0; nc<=this->ncx; ++nc ) {

            T1 x = this->xmin + nc*this->dx;

            for ( T2 nr=0; nr<=this->ncz; ++nr ) {

                T1 z = this->zmin + nr*this->dz;

                if ( nr < this->ncz && nc < this->ncx ) {
                    cell_downRight = nc*this->ncz + nr;
                }
                else {
                    cell_downRight = std::numeric_limits<T2>::max();
                }

                if ( nr > 0 && nc < this->ncx ) {
                    cell_upRight = nc*this->ncz + nr - 1;
                }
                else {
                    cell_upRight = std::numeric_limits<T2>::max();
                }

                if ( nr < this->ncz && nc > 0 ) {
                    cell_downLeft = (nc-1)*this->ncz + nr;
                }
                else {
                    cell_downLeft = std::numeric_limits<T2>::max();
                }

                if ( nr > 0 && nc > 0 ) {
                    cell_upLeft = (nc-1)*this->ncz + nr - 1;
                }
                else {
                    cell_upLeft = std::numeric_limits<T2>::max();
                }

                if ( cell_upLeft != std::numeric_limits<T2>::max() ) {
                    this->nodes[n].pushOwner( cell_upLeft );
                }
                if ( cell_downLeft != std::numeric_limits<T2>::max() ) {
                    this->nodes[n].pushOwner( cell_downLeft );
                }
                if ( cell_upRight != std::numeric_limits<T2>::max() ) {
                    this->nodes[n].pushOwner( cell_upRight );
                }
                if ( cell_downRight != std::numeric_limits<T2>::max() ) {
                    this->nodes[n].pushOwner( cell_downRight );
                }

                this->nodes[n].setX( x );
                this->nodes[n].setZ( z );
                this->nodes[n].setGridIndex( n );
                this->nodes[n].setPrimary(true);

                ++n;
            }
        }

        // continue with secondary nodes
        for ( T2 nc=0; nc<=this->ncx; ++nc ) {

            T1 x = this->xmin + nc*this->dx;

            for ( T2 nr=0; nr<=this->ncz; ++nr ) {

                T1 z = this->zmin + nr*this->dz;

                if ( nr < this->ncz && nc < this->ncx ) {
                    cell_downRight = nc*this->ncz + nr;
                }
                else {
                    cell_downRight = std::numeric_limits<T2>::max();
                }

                if ( nr > 0 && nc < this->ncx ) {
                    cell_upRight = nc*this->ncz + nr - 1;
                }
                else {
                    cell_upRight = std::numeric_limits<T2>::max();
                }

                if ( nr < this->ncz && nc > 0 ) {
                    cell_downLeft = (nc-1)*this->ncz + nr;
                }
                else {
                    cell_downLeft = std::numeric_limits<T2>::max();
                }

                if ( nr > 0 && nc > 0 ) {
                    cell_upLeft = (nc-1)*this->ncz + nr - 1;
                }
                else {
                    cell_upLeft = std::numeric_limits<T2>::max();
                }

                // secondary nodes on the vertical
                if ( nr < this->ncz ) {
                    for (T2 ns=0; ns<nSecondary; ++ns, ++n ) {

                        T1 zsv = this->zmin + nr*this->dz + (ns+1)*dzs;

                        if ( cell_downLeft != std::numeric_limits<T2>::max() ) {
                            this->nodes[n].pushOwner( cell_downLeft );
                        }
                        if ( cell_downRight != std::numeric_limits<T2>::max() ) {
                            this->nodes[n].pushOwner( cell_downRight );
                        }

                        this->nodes[n].setX( x );
                        this->nodes[n].setZ( zsv );
                        this->nodes[n].setGridIndex( n );
                        this->nodes[n].setPrimary(false);
                    }
                }

                // secondary nodes on the horizontal
                if ( nc < this->ncx ) {
                    for ( T2 ns=0; ns<nSecondary; ++ns, ++n ) {

                        T1 xsh = this->xmin + nc*this->dx + (ns+1)*dxs;

                        if ( cell_upRight != std::numeric_limits<T2>::max() ) {
                            this->nodes[n].pushOwner( cell_upRight );
                        }
                        if ( cell_downRight != std::numeric_limits<T2>::max() ) {
                            this->nodes[n].pushOwner( cell_downRight );
                        }

                        this->nodes[n].setX( xsh );
                        this->nodes[n].setZ( z );
                        this->nodes[n].setGridIndex( n );
                        this->nodes[n].setPrimary(false);
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

    template<typename T1, typename T2, typename S>
    void Grid2Drndsp<T1,T2,S>::interpSlownessSecondary() {

        T1 dxs = this->dx/(nSecondary+1);
        T1 dzs = this->dz/(nSecondary+1);

        T2 nnz = this->ncz+1;
        T2 n = (this->ncx+1) * nnz;
        for ( T2 nc=0; nc<=this->ncx; ++nc ) {
            for ( T2 nr=0; nr<=this->ncz; ++nr ) {

                T2 np1 = nc * nnz + nr;
                T2 np2v = np1 + 1;
                T2 np2h = np1 + nnz;

                // secondary nodes on the vertical
                if ( nr < this->ncz ) {
                    T1 slope = (this->nodes[np2v].getNodeSlowness()-this->nodes[np1].getNodeSlowness())/this->dz;
                    for (T2 ns=0; ns<nSecondary; ++ns, ++n ) {

                        T1 s = this->nodes[np1].getNodeSlowness() + slope * (ns+1)*dzs;

                        this->nodes[n].setNodeSlowness( s );
                    }
                }
                // secondary nodes on the horizontal
                if ( nc < this->ncx ) {
                    T1 slope = (this->nodes[np2h].getNodeSlowness()-this->nodes[np1].getNodeSlowness())/this->dx;
                    for ( T2 ns=0; ns<nSecondary; ++ns, ++n ) {

                        T1 s = this->nodes[np1].getNodeSlowness() + slope * (ns+1)*dxs;

                        this->nodes[n].setNodeSlowness( s );
                    }
                }
            }
        }

        //#ifdef VTK
        //    vtkSmartPointer<vtkPolyData> polydata1 = vtkSmartPointer<vtkPolyData>::New();
        //    vtkSmartPointer<vtkPoints> pts1 = vtkSmartPointer<vtkPoints>::New();
        //    vtkSmartPointer<vtkDoubleArray> slow1 =
        //    vtkSmartPointer<vtkDoubleArray>::New();
        //    slow1->SetName("Slowness");
        //    T2 nPrimary = (this->ncx+1) * nnz;
        //    for ( size_t n=0; n<nPrimary; ++n ) {
        //        pts1->InsertNextPoint(this->nodes[n].getX(),
        //                             0.0,
        //                             this->nodes[n].getZ());
        //        slow1->InsertValue(n, this->nodes[n].getNodeSlowness() );
        //    }
        //    polydata1->SetPoints(pts1);
        //    polydata1->GetPointData()->SetScalars(slow1);
        //
        //    vtkSmartPointer<vtkXMLPolyDataWriter> writer1 = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        //    writer1->SetFileName( "tempNodesPrimary.vtp" );
        //    writer1->SetInputData( polydata1 );
        //    writer1->SetDataModeToBinary();
        //    writer1->Update();
        //
        //    vtkSmartPointer<vtkPolyData> polydata2 = vtkSmartPointer<vtkPolyData>::New();
        //    vtkSmartPointer<vtkPoints> pts2 = vtkSmartPointer<vtkPoints>::New();
        //    vtkSmartPointer<vtkDoubleArray> slow2 =
        //    vtkSmartPointer<vtkDoubleArray>::New();
        //    slow2->SetName("Slowness");
        //    for ( size_t n=nPrimary; n<this->nodes.size(); ++n ) {
        //        pts2->InsertNextPoint(this->nodes[n].getX(),
        //                              0.0,
        //                              this->nodes[n].getZ());
        //        slow2->InsertValue(n-nPrimary, this->nodes[n].getNodeSlowness() );
        //    }
        //    polydata2->SetPoints(pts2);
        //    polydata2->GetPointData()->SetScalars(slow2);
        //
        //    vtkSmartPointer<vtkXMLPolyDataWriter> writer2 = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        //    writer2->SetFileName( "tempNodesSecondary.vtp" );
        //    writer2->SetInputData( polydata2 );
        //    writer2->SetDataModeToBinary();
        //    writer2->Update();
        //#endif

    }

    template<typename T1, typename T2, typename S>
    void Grid2Drndsp<T1,T2,S>::addTemporaryNodes(const std::vector<S>& Tx,
                                                 const size_t threadNo) const {

        // clear previously assigned nodes
        tempNodes[threadNo].clear();
        for ( size_t nt=0; nt<tempNeighbors[threadNo].size(); ++nt ) {
            tempNeighbors[threadNo][nt].clear();
        }

        // find cells surrounding Tx
        std::set<T2> txCells;
        for (size_t n=0; n<Tx.size(); ++n) {
            long long i, k;
            this->getIJ(Tx[n], i, k);

            T2 nsx = dynRadius / this->dx;
            T2 nsz = dynRadius / this->dz;

            T1 xstart = this->xmin + (i-nsx-1)*this->dx;
            xstart = xstart < this->xmin ? this->xmin : xstart;
            T1 xstop  = this->xmin + (i+nsx+2)*this->dx;
            xstop = xstop > this->xmax ? this->xmax : xstop;

            T1 zstart = this->zmin + (k-nsz-1)*this->dz;
            zstart = zstart < this->zmin ? this->zmin : zstart;
            T1 zstop  = this->zmin + (k+nsz+2)*this->dz;
            zstop = zstop > this->zmax ? this->zmax : zstop;

            S p;
            for ( p.x=xstart+this->dx/2.; p.x<xstop; p.x+=this->dx ) {
                for ( p.z=zstart+this->dz/2.; p.z<zstop; p.z+=this->dz ) {
                    if ( Tx[n].getDistance(p) < dynRadius ) {
                        txCells.insert( this->getCellNo(p) );
                    }
                }
            }
        }
        if ( verbose )
            std::cout << "\n  *** thread no " << threadNo << ": found " << txCells.size() << " cells within radius ***" << std::endl;

        std::set<T2> adjacentCells(txCells.begin(), txCells.end());

        T2 nTemp = nTertiary * (nSecondary+1);

        T1 dxDyn = this->dx / (nTemp + nSecondary + 1);
        T1 dzDyn = this->dz / (nTemp + nSecondary + 1);

        std::map<std::array<T2,2>,std::vector<T2>> lineMap;
        std::array<T2,2> lineKey;
        typename std::map<std::array<T2,2>,std::vector<T2>>::iterator lineIt;

        T2 nnz = this->ncz+1;

        T1 slope, islown;

        // edge nodes
        T2 nTmpNodes = 0;
        Node2Dnd<T1,T2> tmpNode;

        sij<T2> ind;
        for ( auto cell=txCells.begin(); cell!=txCells.end(); cell++ ) {
            this->getCellIJ(*cell, ind);

            if ( ind.i > 0 )
                adjacentCells.insert( (ind.i-1) * this->ncz + ind.j );
            if ( ind.i < this->ncx-1 )
                adjacentCells.insert( (ind.i+1) * this->ncz + ind.j );
            if ( ind.j > 0 )
                adjacentCells.insert( ind.i * this->ncz + ind.j-1 );
            if ( ind.j < this->ncz-1 )
                adjacentCells.insert( ind.i * this->ncz + ind.j+1 );

            if ( ind.i > 0 && ind.j > 0 )
                adjacentCells.insert( (ind.i-1) * this->ncz + ind.j-1 );
            if ( ind.i < this->ncx-1 && ind.j > 0 )
                adjacentCells.insert( (ind.i+1) * this->ncz + ind.j-1 );
            if ( ind.i > 0 && ind.j < this->ncz-1 )
                adjacentCells.insert( (ind.i-1) * this->ncz + ind.j+1 );
            if ( ind.i < this->ncx-1 && ind.j < this->ncz-1 )
                adjacentCells.insert( (ind.i+1) * this->ncz + ind.j+1 );

            T1 x0 = this->xmin + ind.i*this->dx;
            T1 z0 = this->zmin + ind.j*this->dz;

            //
            // along X
            //
            for ( T2 k=0; k<2; ++k ) {

                lineKey = { ind.i * nnz + ind.j+k, (ind.i+1) * nnz + ind.j+k };
                std::sort(lineKey.begin(), lineKey.end());

                lineIt = lineMap.find( lineKey );
                if ( lineIt == lineMap.end() ) {
                    // not found, insert new pair
                    lineMap[ lineKey ] = std::vector<T2>(nTemp);

                    slope = (this->nodes[lineKey[1]].getNodeSlowness() - this->nodes[lineKey[0]].getNodeSlowness())/this->dx;

                    size_t nd = 0;
                    for ( size_t n2=0; n2<nSecondary+1; ++n2 ) {
                        for ( size_t n3=0; n3<nTertiary; ++n3 ) {
                            tmpNode.setXZindex(x0 + (1+n2*(nTertiary+1)+n3)*dxDyn,
                                               z0 + k*this->dz,
                                               nPermanent+nTmpNodes );
                            islown = this->nodes[lineKey[0]].getNodeSlowness() + slope * tmpNode.getDistance(this->nodes[lineKey[0]]);
                            tmpNode.setNodeSlowness( islown );

                            lineMap[lineKey][nd++] = nTmpNodes++;
                            tempNodes[threadNo].push_back( tmpNode );
                            tempNodes[threadNo].back().pushOwner( *cell );
                        }
                    }
                } else {
                    for ( size_t n=0; n<lineIt->second.size(); ++n ) {
                        // setting owners
                        tempNodes[threadNo][ lineIt->second[n] ].pushOwner( *cell );
                    }
                }
            }

            //
            // along Z
            //
            for ( T2 i=0; i<2; ++i ) {

                lineKey = { (ind.i+i) * nnz + ind.j, (ind.i+i) * nnz + ind.j + 1 };
                std::sort(lineKey.begin(), lineKey.end());

                slope = (this->nodes[lineKey[1]].getNodeSlowness() - this->nodes[lineKey[0]].getNodeSlowness())/this->dz;

                lineIt = lineMap.find( lineKey );
                if ( lineIt == lineMap.end() ) {
                    // not found, insert new pair
                    lineMap[ lineKey ] = std::vector<T2>(nTemp);

                    size_t nd = 0;
                    for ( size_t n2=0; n2<nSecondary+1; ++n2 ) {
                        for ( size_t n3=0; n3<nTertiary; ++n3 ) {
                            tmpNode.setXZindex(x0 + i*this->dx,
                                               z0 + (1+n2*(nTertiary+1)+n3)*dzDyn,
                                               nPermanent+nTmpNodes );
                            islown = this->nodes[lineKey[0]].getNodeSlowness() + slope * tmpNode.getDistance(this->nodes[lineKey[0]]);
                            tmpNode.setNodeSlowness( islown );

                            lineMap[lineKey][nd++] = nTmpNodes++;
                            tempNodes[threadNo].push_back( tmpNode );
                            tempNodes[threadNo].back().pushOwner( *cell );
                        }
                    }
                } else {
                    for ( size_t n=0; n<lineIt->second.size(); ++n ) {
                        // setting owners
                        tempNodes[threadNo][ lineIt->second[n] ].pushOwner( *cell );
                    }
                }
            }
        }

        for ( auto cell=txCells.begin(); cell!=txCells.end(); ++cell ) {
            adjacentCells.erase(*cell);
        }
        for ( auto adj=adjacentCells.begin(); adj!=adjacentCells.end(); ++adj ) {
            this->getCellIJ(*adj, ind);

            //
            // along X
            //
            for ( T2 k=0; k<2; ++k ) {

                lineKey = { ind.i * nnz + ind.j+k, (ind.i+1) * nnz + ind.j+k };
                std::sort(lineKey.begin(), lineKey.end());

                lineIt = lineMap.find( lineKey );
                if ( lineIt != lineMap.end() ) {
                    for ( size_t n=0; n<lineIt->second.size(); ++n ) {
                        // setting owners
                        tempNodes[threadNo][ lineIt->second[n] ].pushOwner( *adj );
                    }
                }
            }

            //
            // along Z
            //
            for ( T2 i=0; i<2; ++i ) {

                lineKey = { (ind.i+i) * nnz + ind.j, (ind.i+i) * nnz + ind.j + 1 };
                std::sort(lineKey.begin(), lineKey.end());

                lineIt = lineMap.find( lineKey );
                if ( lineIt != lineMap.end() ) {
                    for ( size_t n=0; n<lineIt->second.size(); ++n ) {
                        // setting owners
                        tempNodes[threadNo][ lineIt->second[n] ].pushOwner( *adj );
                    }
                }
            }
        }
        for ( T2 n=0; n<tempNodes[threadNo].size(); ++n ) {
            for ( size_t n2=0; n2<tempNodes[threadNo][n].getOwners().size(); ++n2) {
                tempNeighbors[threadNo][ tempNodes[threadNo][n].getOwners()[n2] ].push_back(n);
            }
        }
        if ( verbose )
            std::cout << "  *** thread no " << threadNo << ": " << tempNodes[threadNo].size() << " dynamic nodes were added ***" << std::endl;
    }

    template<typename T1, typename T2, typename S>
    void Grid2Drndsp<T1,T2,S>::initQueue(const std::vector<S>& Tx,
                                         const std::vector<T1>& t0,
                                         std::priority_queue<Node2Dn<T1,T2>*,
                                         std::vector<Node2Dn<T1,T2>*>,
                                         CompareNodePtr<T1>>& queue,
                                         std::vector<Node2Dnd<T1,T2>>& txNodes,
                                         std::vector<bool>& inQueue,
                                         std::vector<bool>& frozen,
                                         const size_t threadNo) const {
        for (size_t n=0; n<Tx.size(); ++n) {
            bool found = false;
            for ( size_t nn=0; nn<this->nodes.size(); ++nn ) {
                if ( this->nodes[nn] == Tx[n] ) {
                    found = true;
                    this->nodes[nn].setTT( t0[n], threadNo );
                    queue.push( &(this->nodes[nn]) );
                    inQueue[nn] = true;
                    frozen[nn] = true;
                    break;
                }
            }
            if ( found==false ) {
                for ( size_t nn=0; nn<tempNodes[threadNo].size(); ++nn ) {
                    if ( tempNodes[threadNo][nn] == Tx[n] ) {
                        found = true;
                        tempNodes[threadNo][nn].setTT( t0[n], 0 );
                        queue.push( &(tempNodes[threadNo][nn]) );
                        inQueue[nPermanent+nn] = true;
                        frozen[nPermanent+nn] = true;
                        break;
                    }
                }
            }
            if ( found==false ) {
                // If Tx[n] is not on a node, we create a new node and initialize the queue:
                txNodes.push_back( Node2Dnd<T1,T2>(t0[n], Tx[n].x, Tx[n].z, 1, 0));
                txNodes.back().pushOwner( this->getCellNo(Tx[n]) );
                txNodes.back().setGridIndex( static_cast<T2>(this->nodes.size()+
                                                             tempNodes.size()+
                                                             txNodes.size()-1) );
                frozen.push_back( true );

                T1 slo = this->getSlowness( Tx[n] );
                txNodes.back().setNodeSlowness( slo );

                queue.push( &(txNodes.back()) );
                inQueue.push_back( true );

            }
        }
    }

    template<typename T1, typename T2, typename S>
    void Grid2Drndsp<T1,T2,S>::propagate(std::priority_queue<Node2Dn<T1,T2>*,
                                         std::vector<Node2Dn<T1,T2>*>,
                                         CompareNodePtr<T1>>& queue,
                                         std::vector<bool>& inQueue,
                                         std::vector<bool>& frozen,
                                         const size_t threadNo) const {
        while ( !queue.empty() ) {
            const Node2Dn<T1,T2>* src = queue.top();
            queue.pop();
            inQueue[ src->getGridIndex() ] = false;
            frozen[ src->getGridIndex() ] = true;

            for ( size_t no=0; no<src->getOwners().size(); ++no ) {

                T2 cellNo = src->getOwners()[no];

                for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                    T2 neibNo = this->neighbors[cellNo][k];
                    if ( neibNo == src->getGridIndex() || frozen[neibNo] ) {
                        continue;
                    }

                    // compute dt
                    T1 dt = this->computeDt(*src, this->nodes[neibNo]);

                    if (src->getTT(threadNo)+dt < this->nodes[neibNo].getTT(threadNo)) {
                        this->nodes[neibNo].setTT( src->getTT(threadNo)+dt, threadNo );

                        if ( !inQueue[neibNo] ) {
                            queue.push( &(this->nodes[neibNo]) );
                            inQueue[neibNo] = true;
                        }
                    }
                }

                for ( size_t k=0; k < tempNeighbors[threadNo][cellNo].size(); ++k ) {
                    T2 neibNo = tempNeighbors[threadNo][cellNo][k];
                    if ( neibNo == src->getGridIndex()-nPermanent || frozen[nPermanent+neibNo] ) {
                        continue;
                    }

                    // compute dt
                    T1 dt = this->computeDt(*src, tempNodes[threadNo][neibNo]);

                    if (src->getTT(threadNo)+dt < tempNodes[threadNo][neibNo].getTT(0)) {
                        tempNodes[threadNo][neibNo].setTT( src->getTT(threadNo)+dt, 0 );

                        if ( !inQueue[nPermanent+neibNo] ) {
                            queue.push( &(tempNodes[threadNo][neibNo]) );
                            inQueue[nPermanent+neibNo] = true;
                        }
                    }
                }
            }
        }
    }

    template<typename T1, typename T2, typename S>
    void Grid2Drndsp<T1,T2,S>::raytrace(const std::vector<S>& Tx,
                                        const std::vector<T1>& t0,
                                        const std::vector<S>& Rx,
                                        const size_t threadNo) const {
        this->checkPts(Tx);
        this->checkPts(Rx);

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node2Dn<T1,T2>*, std::vector<Node2Dn<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );

        addTemporaryNodes(Tx, threadNo);

        std::vector<Node2Dnd<T1,T2>> txNodes;
        std::vector<bool> inQueue( this->nodes.size()+tempNodes[threadNo].size(), false );
        std::vector<bool> frozen( this->nodes.size()+tempNodes[threadNo].size(), false );

        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);

        propagate(queue, inQueue, frozen, threadNo);
    }

    template<typename T1, typename T2, typename S>
    void Grid2Drndsp<T1,T2,S>::raytrace(const std::vector<S>& Tx,
                                        const std::vector<T1>& t0,
                                        const std::vector<const std::vector<S>*>& Rx,
                                        const size_t threadNo) const {
        this->checkPts(Tx);
        for ( size_t n=0; n<Rx.size(); ++n ) {
            this->checkPts(*Rx[n]);
        }

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node2Dn<T1,T2>*, std::vector<Node2Dn<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );

        addTemporaryNodes(Tx, threadNo);

        std::vector<Node2Dnd<T1,T2>> txNodes;
        std::vector<bool> inQueue( this->nodes.size()+tempNodes[threadNo].size(), false );
        std::vector<bool> frozen( this->nodes.size()+tempNodes[threadNo].size(), false );

        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);

        propagate(queue, inQueue, frozen, threadNo);
    }

}
#endif /* Grid2Drndsp_h */

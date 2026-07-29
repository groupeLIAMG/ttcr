//
//  Grid3Drcdsp.h
//  ttcr
//
//  Created by Bernard Giroux on 2018-11-23.
//  Copyright © 2018 Bernard Giroux. All rights reserved.
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

/**
 * @file Grid3Drcdsp.h
 * @brief Dynamic shortest-path solver on a 3-D rectilinear grid with cell-based
 *        slowness.
 *
 * Declares ttcr::Grid3Drcdsp, which concentrates its node refinement near the
 * source rather than spreading it uniformly as ttcr::Grid3Drcsp does. See
 * @ref g3drcdsp_dyn.
 *
 * @sa Grid3Drc.h, Grid3Drcsp.h, Grid2Drcdsp.h, Node3Dcd.h
 */

#ifndef ttcr_Grid3Drcdsp_h
#define ttcr_Grid3Drcdsp_h

#ifdef VTK
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkXMLPolyDataWriter.h"
#endif

#include "Grid3Drc.h"
#include "Node3Dc.h"
#include "Node3Dcd.h"

namespace ttcr {

    /**
     * @brief Dynamic shortest-path solver on a rectilinear cell-slowness grid.
     *
     * @tparam T1   floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2   integer type of node and cell indices.
     * @tparam CELL cell policy from Cell.h.
     *
     * @section g3drcdsp_dyn Why "dynamic"
     * The plain shortest-path method (ttcr::Grid3Drcsp) buys accuracy with
     * secondary nodes spread uniformly over the grid, even though the
     * traveltime error concentrates near the source where wavefront curvature is
     * highest. This class keeps a modest number of permanent secondary nodes
     * everywhere and inserts extra **tertiary** nodes only within
     * @ref dynRadius of each source, for that source's solve alone.
     *
     * Since the tertiary nodes depend on where the source is, they cannot live
     * in the shared node vector — a second thread solving a different source
     * would need different ones. They are held in @ref tempNodes, indexed by
     * thread, which is why the class has a per-thread node container.
     *
     * @note Unlike ttcr::Grid2Drcdsp — which is hard-wired to an isotropic
     *       ttcr::Cell — this 3-D version stays templated on @c CELL and so
     *       supports the anisotropic policies.
     *
     * @sa Grid3Drc.h, Grid3Drcsp.h, Grid2Drcdsp.h, Node3Dcd.h
     */
    template<typename T1, typename T2, typename CELL>
    class Grid3Drcdsp : public Grid3Drc<T1,T2,Node3Dc<T1,T2>,CELL> {
    public:
        /**
         * @brief Build the grid and its permanent nodes.
         *
         * @param nx   number of cells along x.
         * @param ny   number of cells along y.
         * @param nz   number of cells along z.
         * @param ddx  cell size along x.
         * @param ddy  cell size along y.
         * @param ddz  cell size along z.
         * @param minx x coordinate of the grid origin.
         * @param miny y coordinate of the grid origin.
         * @param minz z coordinate of the grid origin.
         * @param ns   number of permanent secondary nodes per cell edge, applied
         *             to all three axes.
         * @param ttrp recompute receiver traveltimes along the raypath.
         * @param nd   number of tertiary nodes added per edge near a source
         *             (ttcr::input_parameters::nTertiary).
         * @param drad radius around a source within which tertiary nodes are
         *             inserted (ttcr::input_parameters::radius_tertiary_nodes).
         * @param useEdgeLength if true, @p drad is a **multiple of the mean cell
         *             edge** @f$(ddx+ddy+ddz)/3@f$ rather than an absolute
         *             distance; the constructor rescales it accordingly.
         * @param nt   number of threads; sizes the per-thread temporary
         *             containers.
         * @param _translateOrigin shift the grid origin to (0,0,0).
         *
         * @post @ref nPermanent records the node count before any tertiary nodes
         *       are added, which is how the solver later distinguishes permanent
         *       nodes from temporary ones.
         * @note The parameter order differs from ttcr::Grid2Drcdsp, where
         *       @p ttrp comes after @p drad rather than between @p ns and @p nd.
         */
        Grid3Drcdsp(const T2 nx, const T2 ny, const T2 nz,
                    const T1 ddx, const T1 ddy, const T1 ddz,
                    const T1 minx, const T1 miny, const T1 minz,
                    const T2 ns, const bool ttrp, const T2 nd, const T1 drad,
                    const bool useEdgeLength=true, const size_t nt=1,
                    const bool _translateOrigin=false) :
        Grid3Drc<T1,T2,Node3Dc<T1,T2>,CELL>(nx, ny, nz, ddx, ddy, ddz, minx, miny, minz, ttrp, nt, _translateOrigin),
        nSecondary(ns), nTertiary(nd), nPermanent(0),
        dynRadius(drad),
        tempNodes(std::vector<std::vector<Node3Dcd<T1,T2>>>(nt)),
        tempNeighbors(std::vector<std::vector<std::vector<T2>>>(nt))
        {
            this->buildGridNodes(nSecondary, nSecondary, nSecondary);
            this->template buildGridNeighbors<Node3Dc<T1,T2>>(this->nodes);
            nPermanent = static_cast<T2>(this->nodes.size());
            for ( size_t n=0; n<nt; ++n ) {
                tempNeighbors[n].resize(static_cast<size_t>(this->ncx) * this->ncy * this->ncz);
            }
            if (useEdgeLength) {
                T1 dx = 0.3333333 * (ddx+ddy+ddz);
                dynRadius *= dx;
            }
        }

        /// Destructor.
        ~Grid3Drcdsp() {
        }

    private:
        T2 nSecondary;                 ///< number of permanent secondary nodes per edge
        T2 nTertiary;                  ///< number of tertiary nodes added per edge near a source
        T2 nPermanent;                 ///< total nb of primary & permanent secondary; nodes beyond this index are temporary
        T1 dynRadius;                  ///< radius around a source within which tertiary nodes are inserted, already scaled by the mean edge length if requested

        /// Temporary (tertiary) nodes, one set per thread. Kept out of the
        /// shared node vector because their positions depend on the source.
        /// @sa @ref g3drcdsp_dyn
        // we will store temporary nodes in a separate container.  This is to
        // allow threaded computations with different Tx (location of temp
        // nodes vary from one Tx to the other)
        mutable std::vector<std::vector<Node3Dcd<T1,T2>>> tempNodes;
        /// Per-thread, per-cell lists of the temporary nodes falling in each cell.
        mutable std::vector<std::vector<std::vector<T2>>> tempNeighbors;

        /**
         * @brief Insert tertiary nodes around the sources, for one thread.
         * @param Tx       source positions.
         * @param threadNo thread whose @ref tempNodes set to fill.
         * @post That thread's temporary containers are cleared and repopulated,
         *       so the call is safe to repeat for a new source.
         */
        void addTemporaryNodes(const std::vector<sxyz<T1>>& Tx, const size_t threadNo) const;

        /**
         * @brief Seed the priority queue with the source nodes.
         * @param[in]  Tx       source positions.
         * @param[in]  t0       origin time of each source.
         * @param[out] queue    min-queue ordered by traveltime.
         * @param[out] txNodes  nodes created at the source positions themselves,
         *                      when a source does not coincide with a grid node.
         * @param[out] inQueue  per-node flag, whether it is already queued.
         * @param[out] frozen   per-node flag, whether its traveltime is final.
         * @param[in]  threadNo thread to work on.
         */
        void initQueue(const std::vector<sxyz<T1>>& Tx,
                       const std::vector<T1>& t0,
                       std::priority_queue<Node3Dc<T1,T2>*,
                       std::vector<Node3Dc<T1,T2>*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<Node3Dcd<T1,T2>>& txNodes,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       const size_t threadNo) const;

        /**
         * @brief Relax the graph until the queue empties.
         * @param[in,out] queue    min-queue of nodes to expand.
         * @param[in,out] inQueue  per-node queued flag.
         * @param[in,out] frozen   per-node finalised flag.
         * @param[in]     threadNo thread to work on.
         * @note Dijkstra expansion: pop the least-traveltime node, freeze it,
         *       relax its neighbours.
         */
        void propagate(std::priority_queue<Node3Dc<T1,T2>*,
                       std::vector<Node3Dc<T1,T2>*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       const size_t threadNo) const;

        /**
         * @brief Propagate the traveltime field and evaluate it at the receivers.
         * @param Tx       source positions.
         * @param t0       origin time of each source.
         * @param Rx       receiver positions.
         * @param threadNo thread to compute on.
         * @note Adds this thread's tertiary nodes first, then runs
         *       @ref initQueue and @ref propagate.
         */
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<sxyz<T1>>& Rx,
                      const size_t threadNo=0) const;

        /**
         * @brief Propagate once and evaluate at several receiver sets.
         * @param Tx       source positions.
         * @param t0       origin time of each source.
         * @param Rx       the receiver sets.
         * @param threadNo thread to compute on.
         */
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      const size_t threadNo=0) const;
    };


    template<typename T1, typename T2, typename CELL>
    void Grid3Drcdsp<T1,T2,CELL>::addTemporaryNodes(const std::vector<sxyz<T1>>& Tx,
                                                    const size_t threadNo) const {

        // clear previously assigned nodes
        tempNodes[threadNo].clear();
        for ( size_t nt=0; nt<tempNeighbors[threadNo].size(); ++nt ) {
            tempNeighbors[threadNo][nt].clear();
        }

        // find cells surrounding Tx
        std::set<T2> txCells;
        for (size_t n=0; n<Tx.size(); ++n) {
            ptrdiff_t i, j, k;
            this->getIJK(Tx[n], i, j, k);

            T2 nsx = dynRadius / this->dx;
            T2 nsy = dynRadius / this->dy;
            T2 nsz = dynRadius / this->dz;

            T1 xstart = this->xmin + (i-nsx-1)*this->dx;
            xstart = xstart < this->xmin ? this->xmin : xstart;
            T1 xstop  = this->xmin + (i+nsx+2)*this->dx;
            xstop = xstop > this->xmax ? this->xmax : xstop;

            T1 ystart = this->ymin + (j-nsy-1)*this->dy;
            ystart = ystart < this->ymin ? this->ymin : ystart;
            T1 ystop  = this->ymin + (j+nsy+2)*this->dy;
            ystop = ystop > this->ymax ? this->ymax : ystop;

            T1 zstart = this->zmin + (k-nsz-1)*this->dz;
            zstart = zstart < this->zmin ? this->zmin : zstart;
            T1 zstop  = this->zmin + (k+nsz+2)*this->dz;
            zstop = zstop > this->zmax ? this->zmax : zstop;

            sxyz<T1> p;
            for ( p.x=xstart+this->dx/2.; p.x<xstop; p.x+=this->dx ) {
                for ( p.y=ystart+this->dy/2.; p.y<ystop; p.y+=this->dy ) {
                    for ( p.z=zstart+this->dz/2.; p.z<zstop; p.z+=this->dz ) {
                        if ( Tx[n].getDistance(p) < dynRadius ) {
                            txCells.insert( this->getCellNo(p) );
                        }
                    }
                }
            }
        }
        if ( verbose )
            std::cout << "\n  *** thread no " << threadNo << ": found " << txCells.size() << " cells within radius ***" << std::endl;

        std::set<T2> adjacentCells(txCells.begin(), txCells.end());

        T2 nTemp = nTertiary * (nSecondary+1);

        T1 dxDyn = this->dx / (nTemp + nSecondary + 1);
        T1 dyDyn = this->dy / (nTemp + nSecondary + 1);
        T1 dzDyn = this->dz / (nTemp + nSecondary + 1);

        std::map<std::array<T2,2>,std::vector<T2>> lineMap;
        std::array<T2,2> lineKey;
        typename std::map<std::array<T2,2>,std::vector<T2>>::iterator lineIt;

        T2 nnx = this->ncx+1;
        T2 nny = this->ncy+1;

        // edge nodes
        T2 nTmpNodes = 0;
        Node3Dcd<T1,T2> tmpNode;

        sijk<T2> ind;
        for ( auto cell=txCells.begin(); cell!=txCells.end(); cell++ ) {
            this->getCellIJK(*cell, ind);

            if ( ind.i > 0 )
                adjacentCells.insert( (ind.k * this->ncy + ind.j) * this->ncx + ind.i-1);
            if ( ind.i < this->ncx-1 )
                adjacentCells.insert( (ind.k * this->ncy + ind.j) * this->ncx + ind.i+1);
            if ( ind.j > 0 )
                adjacentCells.insert( (ind.k * this->ncy + ind.j-1) * this->ncx + ind.i);
            if ( ind.j < this->ncy-1 )
                adjacentCells.insert( (ind.k * this->ncy + ind.j+1) * this->ncx + ind.i);
            if ( ind.k > 0 )
                adjacentCells.insert( ((ind.k-1) * this->ncy + ind.j) * this->ncx + ind.i);
            if ( ind.k < this->ncz-1 )
                adjacentCells.insert( ((ind.k+1) * this->ncy + ind.j) * this->ncx + ind.i);

            if ( ind.i > 0 && ind.j > 0 )
                adjacentCells.insert( (ind.k * this->ncy + ind.j-1) * this->ncx + ind.i-1);
            if ( ind.i < this->ncx-1 && ind.j > 0 )
                adjacentCells.insert( (ind.k * this->ncy + ind.j-1) * this->ncx + ind.i+1);
            if ( ind.i > 0 && ind.j < this->ncy-1 )
                adjacentCells.insert( (ind.k * this->ncy + ind.j+1) * this->ncx + ind.i-1);
            if ( ind.i < this->ncx-1 && ind.j < this->ncy-1 )
                adjacentCells.insert( (ind.k * this->ncy + ind.j+1) * this->ncx + ind.i+1);

            if ( ind.i > 0 && ind.k > 0 )
                adjacentCells.insert( ((ind.k-1) * this->ncy + ind.j) * this->ncx + ind.i-1);
            if ( ind.i < this->ncx-1 && ind.k > 0 )
                adjacentCells.insert( ((ind.k-1) * this->ncy + ind.j) * this->ncx + ind.i+1);
            if ( ind.i > 0 && ind.k < this->ncz-1 )
                adjacentCells.insert( ((ind.k+1) * this->ncy + ind.j) * this->ncx + ind.i-1);
            if ( ind.i < this->ncx-1 && ind.k < this->ncz-1 )
                adjacentCells.insert( ((ind.k+1) * this->ncy + ind.j) * this->ncx + ind.i+1);

            if ( ind.j > 0 && ind.k > 0 )
                adjacentCells.insert( ((ind.k-1) * this->ncy + ind.j-1) * this->ncx + ind.i);
            if ( ind.j < this->ncy-1 && ind.k > 0 )
                adjacentCells.insert( ((ind.k-1) * this->ncy + ind.j+1) * this->ncx + ind.i);
            if ( ind.j > 0 && ind.k < this->ncz-1 )
                adjacentCells.insert( ((ind.k+1) * this->ncy + ind.j-1) * this->ncx + ind.i);
            if ( ind.j < this->ncy-1 && ind.k < this->ncz-1 )
                adjacentCells.insert( ((ind.k+1) * this->ncy + ind.j+1) * this->ncx + ind.i);

            T1 x0 = this->xmin + ind.i*this->dx;
            T1 y0 = this->ymin + ind.j*this->dy;
            T1 z0 = this->zmin + ind.k*this->dz;

            //
            // along X
            //
            for ( T2 j=0; j<2; ++j ) {
                for ( T2 k=0; k<2; ++k ) {

                    lineKey = { ((ind.k+k)*nny + ind.j+j) * nnx + ind.i,
                        ((ind.k+k)*nny + ind.j+j) * nnx + ind.i+1 };
                    std::sort(lineKey.begin(), lineKey.end());

                    lineIt = lineMap.find( lineKey );
                    if ( lineIt == lineMap.end() ) {
                        // not found, insert new pair
                        lineMap[ lineKey ] = std::vector<T2>(nTemp);

                        size_t nd = 0;
                        for ( size_t n2=0; n2<nSecondary+1; ++n2 ) {
                            for ( size_t n3=0; n3<nTertiary; ++n3 ) {
                                tmpNode.setXYZindex(x0 + (1+n2*(nTertiary+1)+n3)*dxDyn,
                                                    y0 + j*this->dy,
                                                    z0 + k*this->dz,
                                                    nPermanent+nTmpNodes );

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

            //
            // along Y
            //
            for ( T2 i=0; i<2; ++i ) {
                for ( T2 k=0; k<2; ++k ) {

                    lineKey = { ((ind.k+k)*nny + ind.j)   * nnx + ind.i+i,
                        ((ind.k+k)*nny + ind.j+1) * nnx + ind.i+i };
                    std::sort(lineKey.begin(), lineKey.end());

                    lineIt = lineMap.find( lineKey );
                    if ( lineIt == lineMap.end() ) {
                        // not found, insert new pair
                        lineMap[ lineKey ] = std::vector<T2>(nTemp);

                        size_t nd = 0;
                        for ( size_t n2=0; n2<nSecondary+1; ++n2 ) {
                            for ( size_t n3=0; n3<nTertiary; ++n3 ) {
                                tmpNode.setXYZindex(x0 + i*this->dx,
                                                    y0 + (1+n2*(nTertiary+1)+n3)*dyDyn,
                                                    z0 + k*this->dz,
                                                    nPermanent+nTmpNodes );

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

            //
            // along Z
            //
            for ( T2 i=0; i<2; ++i ) {
                for ( T2 j=0; j<2; ++j ) {

                    lineKey = { ((ind.k)  *nny + ind.j+j) * nnx + ind.i+i,
                        ((ind.k+1)*nny + ind.j+j) * nnx + ind.i+i };
                    std::sort(lineKey.begin(), lineKey.end());

                    lineIt = lineMap.find( lineKey );
                    if ( lineIt == lineMap.end() ) {
                        // not found, insert new pair
                        lineMap[ lineKey ] = std::vector<T2>(nTemp);

                        size_t nd = 0;
                        for ( size_t n2=0; n2<nSecondary+1; ++n2 ) {
                            for ( size_t n3=0; n3<nTertiary; ++n3 ) {
                                tmpNode.setXYZindex(x0 + i*this->dx,
                                                    y0 + j*this->dy,
                                                    z0 + (1+n2*(nTertiary+1)+n3)*dzDyn,
                                                    nPermanent+nTmpNodes );

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
        }

        std::map<std::array<T2,4>,std::vector<T2>> faceMap;
        std::array<T2,4> faceKey;
        typename std::map<std::array<T2,4>,std::vector<T2>>::iterator faceIt;
        nTemp = (nTertiary * (nSecondary+1) + nSecondary) * (nTertiary * (nSecondary+1) + nSecondary) - nSecondary*nSecondary;

        for ( auto cell=txCells.begin(); cell!=txCells.end(); cell++ ) {
            this->getCellIJK(*cell, ind);

            //
            // XY faces
            //
            for ( T2 k=0; k<2; ++k ) {
                faceKey = { ((ind.k+k)*nny + ind.j) * nnx + ind.i,
                    ((ind.k+k)*nny + ind.j) * nnx + ind.i+1,
                    ((ind.k+k)*nny + ind.j+1) * nnx + ind.i,
                    ((ind.k+k)*nny + ind.j+1) * nnx + ind.i+1 };
                std::sort(faceKey.begin(), faceKey.end());

                faceIt = faceMap.find( faceKey );
                if ( faceIt == faceMap.end() ) {
                    // not found, insert new pair
                    faceMap[ faceKey ] = std::vector<T2>(nTemp);
                } else {
                    for ( size_t n=0; n<faceIt->second.size(); ++n ) {
                        // setting owners
                        tempNodes[threadNo][ faceIt->second[n] ].pushOwner( *cell );
                    }
                    continue;
                }

                T1 x0 = this->xmin + ind.i*this->dx;
                T1 y0 = this->ymin + ind.j*this->dy;
                T1 z0 = this->zmin + ind.k+k*this->dz;

                size_t ifn = 0;
                size_t n0 = 0;
                while (n0 < nTertiary*(nSecondary+1)+nSecondary) {
                    for ( size_t n1=0; n1<nTertiary; ++n1 ) {
                        y0 += dyDyn;
                        for ( size_t n2=0; n2<nTertiary*(nSecondary+1)+nSecondary; ++n2 ) {
                            tmpNode.setXYZindex(x0 + (n2+1)*dxDyn,
                                                y0,
                                                z0,
                                                nPermanent+nTmpNodes );

                            faceMap[faceKey][ifn++] = nTmpNodes++;
                            tempNodes[threadNo].push_back( tmpNode );
                            tempNodes[threadNo].back().pushOwner( *cell );
                        }
                        n0++;
                    }
                    if (n0 == nTertiary*(nSecondary+1)+nSecondary) {
                        break;
                    }
                    y0 += dyDyn;
                    for ( size_t n2=0; n2<nSecondary+1; ++n2 ) {
                        for ( size_t n3=0; n3<nTertiary; ++n3 ) {
                            tmpNode.setXYZindex(x0 + (1+n2*(nTertiary+1)+n3)*dxDyn,
                                                y0,
                                                z0,
                                                nPermanent+nTmpNodes );

                            faceMap[faceKey][ifn++] = nTmpNodes++;
                            tempNodes[threadNo].push_back( tmpNode );
                            tempNodes[threadNo].back().pushOwner( *cell );
                        }
                    }
                    n0++;
                }
            }

            //
            // XZ faces
            //
            for ( T2 j=0; j<2; ++j ) {
                faceKey = { ((ind.k)*nny + ind.j+j) * nnx + ind.i,
                    ((ind.k)*nny + ind.j+j) * nnx + ind.i+1,
                    ((ind.k+1)*nny + ind.j+j) * nnx + ind.i,
                    ((ind.k+1)*nny + ind.j+j) * nnx + ind.i+1 };
                std::sort(faceKey.begin(), faceKey.end());

                faceIt = faceMap.find( faceKey );
                if ( faceIt == faceMap.end() ) {
                    // not found, insert new pair
                    faceMap[ faceKey ] = std::vector<T2>(nTemp);
                } else {
                    for ( size_t n=0; n<faceIt->second.size(); ++n ) {
                        // setting owners
                        tempNodes[threadNo][ faceIt->second[n] ].pushOwner( *cell );
                    }
                    continue;
                }

                T1 x0 = this->xmin + ind.i*this->dx;
                T1 y0 = this->ymin + ind.j+j*this->dy;
                T1 z0 = this->zmin + ind.k*this->dz;

                size_t ifn = 0;
                size_t n0 = 0;
                while (n0 < nTertiary*(nSecondary+1)+nSecondary) {
                    for ( size_t n1=0; n1<nTertiary; ++n1 ) {
                        z0 += dzDyn;
                        for ( size_t n2=0; n2<nTertiary*(nSecondary+1)+nSecondary; ++n2 ) {
                            tmpNode.setXYZindex(x0 + (n2+1)*dxDyn,
                                                y0,
                                                z0,
                                                nPermanent+nTmpNodes );

                            faceMap[faceKey][ifn++] = nTmpNodes++;
                            tempNodes[threadNo].push_back( tmpNode );
                            tempNodes[threadNo].back().pushOwner( *cell );
                        }
                        n0++;
                    }
                    if (n0 == nTertiary*(nSecondary+1)+nSecondary) {
                        break;
                    }
                    z0 += dzDyn;
                    for ( size_t n2=0; n2<nSecondary+1; ++n2 ) {
                        for ( size_t n3=0; n3<nTertiary; ++n3 ) {
                            tmpNode.setXYZindex(x0 + (1+n2*(nTertiary+1)+n3)*dxDyn,
                                                y0,
                                                z0,
                                                nPermanent+nTmpNodes );

                            faceMap[faceKey][ifn++] = nTmpNodes++;
                            tempNodes[threadNo].push_back( tmpNode );
                            tempNodes[threadNo].back().pushOwner( *cell );
                        }
                    }
                    n0++;
                }
            }

            //
            // YZ faces
            //
            for ( T2 i=0; i<2; ++i ) {
                faceKey = { ((ind.k)*nny + ind.j) * nnx + ind.i+i,
                    ((ind.k)*nny + ind.j+1) * nnx + ind.i+i,
                    ((ind.k+1)*nny + ind.j) * nnx + ind.i+i,
                    ((ind.k+1)*nny + ind.j+1) * nnx + ind.i+i };
                std::sort(faceKey.begin(), faceKey.end());

                faceIt = faceMap.find( faceKey );
                if ( faceIt == faceMap.end() ) {
                    // not found, insert new pair
                    faceMap[ faceKey ] = std::vector<T2>(nTemp);
                } else {
                    for ( size_t n=0; n<faceIt->second.size(); ++n ) {
                        // setting owners
                        tempNodes[threadNo][ faceIt->second[n] ].pushOwner( *cell );
                    }
                    continue;
                }

                T1 x0 = this->xmin + ind.i+i*this->dx;
                T1 y0 = this->ymin + ind.j*this->dy;
                T1 z0 = this->zmin + ind.k*this->dz;

                size_t ifn = 0;
                size_t n0 = 0;
                while (n0 < nTertiary*(nSecondary+1)+nSecondary) {
                    for ( size_t n1=0; n1<nTertiary; ++n1 ) {
                        z0 += dzDyn;
                        for ( size_t n2=0; n2<nTertiary*(nSecondary+1)+nSecondary; ++n2 ) {
                            tmpNode.setXYZindex(x0,
                                                y0 + (n2+1)*dyDyn,
                                                z0,
                                                nPermanent+nTmpNodes );

                            faceMap[faceKey][ifn++] = nTmpNodes++;
                            tempNodes[threadNo].push_back( tmpNode );
                            tempNodes[threadNo].back().pushOwner( *cell );
                        }
                        n0++;
                    }
                    if (n0 == nTertiary*(nSecondary+1)+nSecondary) {
                        break;
                    }
                    z0 += dzDyn;
                    for ( size_t n2=0; n2<nSecondary+1; ++n2 ) {
                        for ( size_t n3=0; n3<nTertiary; ++n3 ) {
                            tmpNode.setXYZindex(x0,
                                                y0 + (1+n2*(nTertiary+1)+n3)*dyDyn,
                                                z0,
                                                nPermanent+nTmpNodes );

                            faceMap[faceKey][ifn++] = nTmpNodes++;
                            tempNodes[threadNo].push_back( tmpNode );
                            tempNodes[threadNo].back().pushOwner( *cell );
                        }
                    }
                    n0++;
                }
            }
        }

        for ( auto cell=txCells.begin(); cell!=txCells.end(); ++cell ) {
            adjacentCells.erase(*cell);
        }
        for ( auto adj=adjacentCells.begin(); adj!=adjacentCells.end(); ++adj ) {
            this->getCellIJK(*adj, ind);

            //
            // along X
            //
            for ( T2 j=0; j<2; ++j ) {
                for ( T2 k=0; k<2; ++k ) {

                    lineKey = { ((ind.k+k)*nny + ind.j+j) * nnx + ind.i,
                        ((ind.k+k)*nny + ind.j+j) * nnx + ind.i+1 };
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

            //
            // along Y
            //
            for ( T2 i=0; i<2; ++i ) {
                for ( T2 k=0; k<2; ++k ) {

                    lineKey = { ((ind.k+k)*nny + ind.j)   * nnx + ind.i+i,
                        ((ind.k+k)*nny + ind.j+1) * nnx + ind.i+i };
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

            //
            // along Z
            //
            for ( T2 i=0; i<2; ++i ) {
                for ( T2 j=0; j<2; ++j ) {

                    lineKey = { ((ind.k)  *nny + ind.j+j) * nnx + ind.i+i,
                        ((ind.k+1)*nny + ind.j+j) * nnx + ind.i+i };
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

            //
            // XY faces
            //
            for ( T2 k=0; k<2; ++k ) {
                faceKey = { ((ind.k+k)*nny + ind.j) * nnx + ind.i,
                    ((ind.k+k)*nny + ind.j) * nnx + ind.i+1,
                    ((ind.k+k)*nny + ind.j+1) * nnx + ind.i,
                    ((ind.k+k)*nny + ind.j+1) * nnx + ind.i+1 };
                std::sort(faceKey.begin(), faceKey.end());

                faceIt = faceMap.find( faceKey );
                if ( faceIt != faceMap.end() ) {
                    for ( size_t n=0; n<faceIt->second.size(); ++n ) {
                        // setting owners
                        tempNodes[threadNo][ faceIt->second[n] ].pushOwner( *adj );
                    }
                }
            }

            //
            // XZ faces
            //
            for ( T2 j=0; j<2; ++j ) {
                faceKey = { ((ind.k)*nny + ind.j+j) * nnx + ind.i,
                    ((ind.k)*nny + ind.j+j) * nnx + ind.i+1,
                    ((ind.k+1)*nny + ind.j+j) * nnx + ind.i,
                    ((ind.k+1)*nny + ind.j+j) * nnx + ind.i+1 };
                std::sort(faceKey.begin(), faceKey.end());

                faceIt = faceMap.find( faceKey );
                if ( faceIt != faceMap.end() ) {
                    for ( size_t n=0; n<faceIt->second.size(); ++n ) {
                        // setting owners
                        tempNodes[threadNo][ faceIt->second[n] ].pushOwner( *adj );
                    }
                }
            }

            //
            // YZ faces
            //
            for ( T2 i=0; i<2; ++i ) {
                faceKey = { ((ind.k)*nny + ind.j) * nnx + ind.i+i,
                    ((ind.k)*nny + ind.j+1) * nnx + ind.i+i,
                    ((ind.k+1)*nny + ind.j) * nnx + ind.i+i,
                    ((ind.k+1)*nny + ind.j+1) * nnx + ind.i+i };
                std::sort(faceKey.begin(), faceKey.end());

                faceIt = faceMap.find( faceKey );
                if ( faceIt != faceMap.end() ) {
                    for ( size_t n=0; n<faceIt->second.size(); ++n ) {
                        // setting owners
                        tempNodes[threadNo][ faceIt->second[n] ].pushOwner( *adj );
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

    template<typename T1, typename T2, typename CELL>
    void Grid3Drcdsp<T1,T2,CELL>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                           const std::vector<T1>& t0,
                                           const std::vector<sxyz<T1>>& Rx,
                                           const size_t threadNo) const {
        this->checkPts(Tx, true);
        this->checkPts(Rx, true);

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node3Dc<T1,T2>*, std::vector<Node3Dc<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );

        addTemporaryNodes(Tx, threadNo);

        std::vector<Node3Dcd<T1,T2>> txNodes;
        std::vector<bool> inQueue( this->nodes.size()+tempNodes[threadNo].size(), false );
        std::vector<bool> frozen( this->nodes.size()+tempNodes[threadNo].size(), false );

        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);

        propagate(queue, inQueue, frozen, threadNo);
    }

    template<typename T1, typename T2, typename CELL>
    void Grid3Drcdsp<T1,T2,CELL>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                           const std::vector<T1>& t0,
                                           const std::vector<std::vector<sxyz<T1>>>& Rx,
                                           const size_t threadNo) const {
        this->checkPts(Tx, true);
        for ( size_t n=0; n<Rx.size(); ++n ) {
            this->checkPts(Rx[n], true);
        }

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node3Dc<T1,T2>*, std::vector<Node3Dc<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );

        addTemporaryNodes(Tx, threadNo);

        std::vector<Node3Dcd<T1,T2>> txNodes;
        std::vector<bool> inQueue( this->nodes.size()+tempNodes[threadNo].size(), false );
        std::vector<bool> frozen( this->nodes.size()+tempNodes[threadNo].size(), false );

        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);

        propagate(queue, inQueue, frozen, threadNo);
    }


    template<typename T1, typename T2, typename CELL>
    void Grid3Drcdsp<T1,T2,CELL>::initQueue(const std::vector<sxyz<T1>>& Tx,
                                            const std::vector<T1>& t0,
                                            std::priority_queue<Node3Dc<T1,T2>*,
                                            std::vector<Node3Dc<T1,T2>*>,
                                            CompareNodePtr<T1>>& queue,
                                            std::vector<Node3Dcd<T1,T2>>& txNodes,
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
                txNodes.push_back( Node3Dcd<T1,T2>(t0[n], Tx[n].x, Tx[n].y, Tx[n].z, 1, 0));
                txNodes.back().pushOwner( this->getCellNo(Tx[n]) );
                txNodes.back().setGridIndex( static_cast<T2>(this->nodes.size()+
                                                             tempNodes.size()+
                                                             txNodes.size()-1) );
                frozen.push_back( true );

                queue.push( &(txNodes.back()) );
                inQueue.push_back( true );

            }
        }
    }

    template<typename T1, typename T2, typename CELL>
    void Grid3Drcdsp<T1,T2,CELL>::propagate(std::priority_queue<Node3Dc<T1,T2>*,
                                            std::vector<Node3Dc<T1,T2>*>,
                                            CompareNodePtr<T1>>& queue,
                                            std::vector<bool>& inQueue,
                                            std::vector<bool>& frozen,
                                            const size_t threadNo) const {

        while ( !queue.empty() ) {
            const Node3Dc<T1,T2>* src = queue.top();
            queue.pop();
            inQueue[ src->getGridIndex() ] = false;
            frozen[ src->getGridIndex() ] = true;

            T1 srcTT;
            if ( src->getGridIndex() >= nPermanent )
                srcTT = src->getTT(0);
            else
                srcTT = src->getTT(threadNo);

            for ( size_t no=0; no<src->getOwners().size(); ++no ) {

                T2 cellNo = src->getOwners()[no];

                for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                    T2 neibNo = this->neighbors[cellNo][k];
                    if ( neibNo == src->getGridIndex() || frozen[neibNo] ) {
                        continue;
                    }

                    // compute dt
                    T1 dt = this->cells.computeDt(*src, this->nodes[neibNo], cellNo);

                    if (srcTT+dt < this->nodes[neibNo].getTT(threadNo)) {
                        this->nodes[neibNo].setTT( srcTT+dt, threadNo );

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
                    T1 dt = this->cells.computeDt(*src, tempNodes[threadNo][neibNo], cellNo);
                    if (srcTT+dt < tempNodes[threadNo][neibNo].getTT(0)) {
                        tempNodes[threadNo][neibNo].setTT(srcTT+dt,0);

                        if ( !inQueue[nPermanent+neibNo] ) {
                            queue.push( &(tempNodes[threadNo][neibNo]) );
                            inQueue[nPermanent+neibNo] = true;
                        }
                    }
                }
            }
        }
    }

}
#endif /* Grid3Drcdsp_h */

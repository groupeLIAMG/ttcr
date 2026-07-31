//
//  Grid2Ducdsp.h
//  ttcr
//
//  Created by Bernard Giroux on 2021-02-23.
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

/**
 * @file Grid2Ducdsp.h
 * @brief Dynamic shortest-path solver on a 2-D triangular mesh with cell-based
 *        slowness.
 *
 * Declares ttcr::Grid2Ducdsp, which concentrates node refinement near the
 * source rather than spreading it uniformly as ttcr::Grid2Ducsp does.
 *
 * @sa Grid2Duc.h, Grid2Drcdsp.h, Grid2Ducsp.h, Node2Dcd.h
 */

#ifndef ttcr_Grid2Ducdsp_h
#define ttcr_Grid2Ducdsp_h

#include "Grid2Duc.h"
#include "Node2Dc.h"
#include "Node2Dcd.h"

namespace ttcr {

    template<typename T1, typename T2, typename S>
    /**
     * @brief Dynamic shortest-path solver on a triangular mesh, slowness per cell.
     *
     * @tparam T1 floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2 integer type of node and cell indices.
     * @tparam S  point type, @ref sxz or @ref sxyz.
     *
     * The unstructured counterpart of ttcr::Grid2Drcdsp. Rather than spreading
     * secondary nodes uniformly as ttcr::Grid2Ducsp does, it keeps a modest
     * permanent set everywhere and inserts extra **tertiary** nodes only within
     * @ref dyn_radius of each source, where wavefront curvature — and hence the
     * traveltime error — is greatest.
     *
     * Tertiary nodes depend on the source position, so they live per-thread in
     * @ref tempNodes rather than in the shared node vector, letting different
     * workers solve different sources concurrently.
     *
     * @note On a mesh there is no cell size to scale the radius against, so
     *       @ref getAverageEdgeLength supplies the length scale — see the
     *       @c useEdgeLength constructor argument.
     * @note Not templated on a @c CELL policy, so isotropic only.
     *
     * @sa Grid2Duc.h, Grid2Drcdsp.h, Grid2Ducsp.h, Node2Dcd.h
     */
    class Grid2Ducdsp : public Grid2Duc<T1,T2,S,Node2Dc<T1,T2>,Cell<T1,Node2Dc<T1,T2>,sxz<T1>>> {
    public:
        /**
         * @brief Build the mesh and its permanent nodes.
         *
         * @param no   node coordinates.
         * @param tri  triangles, each naming three node indices.
         * @param ns   number of permanent secondary nodes per triangle edge.
         * @param nd   number of tertiary nodes added per edge near a source
         *             (ttcr::input_parameters::nTertiary).
         * @param drad radius around a source within which tertiary nodes are
         *             inserted (ttcr::input_parameters::radius_tertiary_nodes).
         * @param ttrp recompute receiver traveltimes along the raypath.
         * @param useEdgeLength if true, @p drad is a **multiple of the mesh's
         *             average edge length** rather than an absolute distance;
         *             the constructor rescales it accordingly.
         * @param nt   number of threads; sizes the per-thread temporary
         *             containers.
         *
         * @post @ref nPermanent records the node count before any tertiary nodes
         *       are added, which is how the solver later tells permanent nodes
         *       from temporary ones. Slowness is not set.
         */
        Grid2Ducdsp(const std::vector<S>& no,
                    const std::vector<triangleElem<T2>>& tri,
                    const T2 ns, const int nd, const T1 drad, const bool ttrp,
                    const bool useEdgeLength=true, const size_t nt=1) :
        Grid2Duc<T1,T2,S,Node2Dc<T1,T2>,Cell<T1, Node2Dc<T1, T2>, sxz<T1>>>(no, tri, ttrp, nt),
        nSecondary(ns), nTertiary(nd), nPermanent(0),
        dyn_radius(drad),
        tempNodes(std::vector<std::vector<Node2Dcd<T1,T2>>>(nt)),
        tempNeighbors(std::vector<std::vector<std::vector<T2>>>(nt))
        {
            this->buildGridNodes(no, ns, nt);
            this->template buildGridNeighbors<Node2Dc<T1,T2>>(this->nodes);
            nPermanent = static_cast<T2>(this->nodes.size());
            for ( size_t n=0; n<nt; ++n ) {
                tempNeighbors[n].resize(tri.size());
            }
            if (useEdgeLength) dyn_radius *= this->getAverageEdgeLength();
        }
        /// Destructor.
        ~Grid2Ducdsp() {
        }

    private:
        T2 nSecondary;   ///< number of permanent secondary nodes per triangle edge
        T2 nTertiary;    ///< number of tertiary nodes added per edge near a source
        T2 nPermanent;   ///< total nb of primary & permanent secondary; nodes beyond this index are temporary
        T1 dyn_radius;   ///< radius around a source within which tertiary nodes are inserted, already scaled by the average edge length if requested

        /// Temporary (tertiary) nodes, one set per thread. Kept out of the
        /// shared node vector because their positions depend on the source.
        // we will store temporary nodes in a separate container.  This is to
        // allow threaded computations with different Tx (location of temp
        // nodes vary from one Tx to the other)
        mutable std::vector<std::vector<Node2Dcd<T1,T2>>> tempNodes;
        /// Per-thread, per-triangle lists of the temporary nodes in each triangle.
        mutable std::vector<std::vector<std::vector<T2>>> tempNeighbors;

        /**
         * @brief Insert tertiary nodes around the sources, for one thread.
         * @param Tx       source positions.
         * @param threadNo thread whose @ref tempNodes set to fill.
         * @post That thread's temporary containers are cleared and repopulated,
         *       so the call is safe to repeat for a new source.
         */
        void addTemporaryNodes(const std::vector<S>& Tx, const size_t threadNo) const;

        /**
         * @brief Seed the priority queue with the source nodes.
         * @param[in]  Tx       source positions.
         * @param[in]  t0       origin time of each source.
         * @param[out] queue    min-queue ordered by traveltime.
         * @param[out] txNodes  nodes created at source positions that do not
         *                      coincide with a mesh node.
         * @param[out] inQueue  per-node flag, whether it is already queued.
         * @param[out] frozen   per-node flag, whether its traveltime is final.
         * @param[in]  threadNo thread to work on.
         */
        void initQueue(const std::vector<S>& Tx,
                       const std::vector<T1>& t0,
                       std::priority_queue<Node2Dc<T1,T2>*,
                       std::vector<Node2Dc<T1,T2>*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<Node2Dcd<T1,T2>>& txNodes,
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
        void propagate(std::priority_queue<Node2Dc<T1,T2>*,
                       std::vector<Node2Dc<T1,T2>*>,
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
        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<S>& Rx,
                      const size_t threadNo=0) const;

        /**
         * @brief Propagate once and evaluate at several receiver sets.
         * @param Tx       source positions.
         * @param t0       origin time of each source.
         * @param Rx       pointers to the receiver sets.
         * @param threadNo thread to compute on.
         */
        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<const std::vector<S>*>& Rx,
                      const size_t threadNo=0) const;
    };

    template<typename T1, typename T2, typename S>
    void Grid2Ducdsp<T1,T2,S>::addTemporaryNodes(const std::vector<S>& Tx,
                                                 const size_t threadNo) const {

        // clear previously assigned nodes
        tempNodes[threadNo].clear();
        for ( size_t nt=0; nt<this->triangles.size(); ++nt ) {
            tempNeighbors[threadNo][nt].clear();
        }

        // find cells surrounding Tx
        std::set<T2> txCells;
        for ( size_t nt=0; nt<this->triangles.size(); ++nt ) {
            // centroid of tet
            S cent = S(this->nodes[this->triangles[nt].i[0]] +
                       this->nodes[this->triangles[nt].i[1]]);
            cent += this->nodes[this->triangles[nt].i[2]];
            cent *= 0.33333333;
            for (size_t n=0; n<Tx.size(); ++n) {
                if ( cent.getDistance(Tx[n]) <= this->dyn_radius ) {
                    txCells.insert(static_cast<T2>(nt));
                }
            }
        }
        if ( verbose )
            std::cout << "\n  *** thread no " << threadNo << ": found " << txCells.size() << " cells within radius ***" << std::endl;

        std::set<T2> adjacentCells(txCells.begin(), txCells.end());

        std::map<std::array<T2,2>,std::vector<T2>> lineMap;
        std::array<T2,2> lineKey;
        typename std::map<std::array<T2,2>,std::vector<T2>>::iterator lineIt;

        // edge nodes
        T2 nTmpNodes = 0;
        Node2Dcd<T1,T2> tmpNode;
        size_t nDynTot = static_cast<size_t>(nTertiary) * (nSecondary+1);  // total number of dynamic nodes on edges

        for ( auto cell=txCells.begin(); cell!=txCells.end(); cell++ ) {

            //  adjacent cells to the tetrahedron where new nodes will be added
            for ( size_t i=0; i<4; ++i ) {
                T2 vertex = this->neighbors[*cell][i];
                for ( size_t c=0; c<this->nodes[vertex].getOwners().size(); ++c ) {
                    adjacentCells.insert(this->nodes[vertex].getOwners()[c]);
                }
            }

            for ( size_t nl=0; nl<3; ++nl ) {

                lineKey = {this->triangles[*cell].i[ nl ],
                    this->triangles[*cell].i[ (nl+1)%3 ]};
                std::sort(lineKey.begin(), lineKey.end());

                lineIt = lineMap.find( lineKey );
                if ( lineIt == lineMap.end() ) {
                    // not found, insert new pair
                    lineMap[ lineKey ] = std::vector<T2>(nDynTot);
                } else {
                    for ( size_t n=0; n<lineIt->second.size(); ++n ) {
                        // setting owners
                        tempNodes[threadNo][ lineIt->second[n] ].pushOwner( *cell );
                    }
                    continue;
                }

                S d = (this->nodes[lineKey[1]]-this->nodes[lineKey[0]])/static_cast<T1>(nDynTot+nSecondary+1);

                size_t nd = 0;
                for ( size_t n2=0; n2<nSecondary+1; ++n2 ) {
                    for ( size_t n3=0; n3<nTertiary; ++n3 ) {
                        tmpNode.setXZindex(this->nodes[lineKey[0]].getX()+(1+n2*(nTertiary+1)+n3)*d.x,
                                           this->nodes[lineKey[0]].getZ()+(1+n2*(nTertiary+1)+n3)*d.z,
                                           nPermanent+nTmpNodes );

                        lineMap[lineKey][nd++] = nTmpNodes++;
                        tempNodes[threadNo].push_back( tmpNode );
                        tempNodes[threadNo].back().pushOwner( *cell );
                    }
                }
            }

        }

        for ( auto cell=txCells.begin(); cell!=txCells.end(); ++cell ) {
            adjacentCells.erase(*cell);
        }
        for ( auto adj=adjacentCells.begin(); adj!=adjacentCells.end(); ++adj ) {

            for ( size_t nl=0; nl<3; ++nl ) {

                lineKey = {this->triangles[*adj].i[ nl ],
                    this->triangles[*adj].i[ (nl+1)%3 ]};
                std::sort(lineKey.begin(), lineKey.end());

                lineIt = lineMap.find( lineKey );
                if ( lineIt != lineMap.end() ) {
                    // setting owners
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
    void Grid2Ducdsp<T1,T2,S>::initQueue(const std::vector<S>& Tx,
                                         const std::vector<T1>& t0,
                                         std::priority_queue<Node2Dc<T1,T2>*,
                                         std::vector<Node2Dc<T1,T2>*>,
                                         CompareNodePtr<T1>>& queue,
                                         std::vector<Node2Dcd<T1,T2>>& txNodes,
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
                txNodes.push_back( Node2Dcd<T1,T2>(t0[n], Tx[n].x, Tx[n].z, 1, 0) );
                // we belong to cell index no
                txNodes.back().pushOwner( this->getCellNo(Tx[n]) );
                txNodes.back().setGridIndex( static_cast<T2>(this->nodes.size()+
                                                             txNodes.size()-1) );

                queue.push( &(txNodes.back()) );
                inQueue.push_back( true );
                frozen.push_back( true );
            }
        }
    }

    template<typename T1, typename T2, typename S>
    void Grid2Ducdsp<T1,T2,S>::propagate(std::priority_queue<Node2Dc<T1,T2>*,
                                         std::vector<Node2Dc<T1,T2>*>,
                                         CompareNodePtr<T1>>& queue,
                                         std::vector<bool>& inQueue,
                                         std::vector<bool>& frozen,
                                         const size_t threadNo) const {

        while ( !queue.empty() ) {
            const Node2Dc<T1,T2>* src = queue.top();
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

    template<typename T1, typename T2, typename S>
    void Grid2Ducdsp<T1,T2,S>::raytrace(const std::vector<S>& Tx,
                                        const std::vector<T1>& t0,
                                        const std::vector<S>& Rx,
                                        const size_t threadNo) const {
        this->checkPts(Tx);
        this->checkPts(Rx);

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );

        addTemporaryNodes(Tx, threadNo);

        std::vector<Node2Dcd<T1,T2>> txNodes;
        std::vector<bool> inQueue( this->nodes.size()+tempNodes[threadNo].size(), false );
        std::vector<bool> frozen( this->nodes.size()+tempNodes[threadNo].size(), false );

        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);

        propagate(queue, inQueue, frozen, threadNo);
    }

    template<typename T1, typename T2, typename S>
    void Grid2Ducdsp<T1,T2,S>::raytrace(const std::vector<S>& Tx,
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
        std::priority_queue< Node2Dc<T1,T2>*, std::vector<Node2Dc<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );

        addTemporaryNodes(Tx, threadNo);

        std::vector<Node2Dcd<T1,T2>> txNodes;
        std::vector<bool> inQueue( this->nodes.size()+tempNodes[threadNo].size(), false );
        std::vector<bool> frozen( this->nodes.size()+tempNodes[threadNo].size(), false );

        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);

        propagate(queue, inQueue, frozen, threadNo);
    }

}
#endif /* Grid2Ducdsp_h */

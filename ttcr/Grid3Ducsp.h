//
//  Grid3Ducsp.h
//  ttcr
//
//  Created by Bernard Giroux on 2012-09-22.
//  Copyright (c) 2012 Bernard Giroux. All rights reserved.
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
 * @file Grid3Ducsp.h
 * @brief Shortest-path solver on a 3-D tetrahedral mesh with cell-based slowness.
 *
 * Declares ttcr::Grid3Ducsp, the 3-D counterpart of ttcr::Grid2Ducsp.
 *
 * @sa Grid3Duc.h, Grid2Ducsp.h, Node3Dcsp.h, Grid3Ducdsp.h
 */

#ifndef ttcr_Grid3Ducsp_h
#define ttcr_Grid3Ducsp_h

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <vector>

#include "Grid3Duc.h"
#include "Node3Dcsp.h"
#include "utils.h"

namespace ttcr {

    template<typename T1, typename T2, typename CELL>
    /**
     * @brief Shortest-path eikonal solver on a tetrahedral mesh, slowness per cell.
     *
     * @tparam T1 floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2 integer type of node and cell indices.
     *
     * The 3-D counterpart of ttcr::Grid2Ducsp: the mesh is treated as a graph
     * over primary and secondary nodes and relaxed with a Dijkstra-style
     * propagation. Secondary nodes are placed along the tetrahedron edges to
     * widen the set of ray directions available.
     *
     * Being graph-based, it does not use the local eikonal update of
     * ttcr::Grid3Duc::localUpdate3D at all — that belongs to the sweeping and
     * marching solvers.
     *
     * @sa Grid3Duc.h, Grid2Ducsp.h, Node3Dcsp.h, Grid3Ducdsp.h
     */
    class Grid3Ducsp : public Grid3Duc<T1,T2,Node3Dcsp<T1,T2>,CELL> {
    public:
        /**
         * @brief Build the mesh and its secondary nodes.
         *
         * @param no   node coordinates.
         * @param tet  tetrahedra, each naming four node indices.
         * @param ns   number of secondary nodes per tetrahedron edge
         *             (ttcr::input_parameters::nn).
         * @param rptt recompute receiver traveltimes along the raypath.
         * @param md   minimum step retained when integrating a raypath.
         * @param nt   number of threads.
         * @param _translateOrigin shift the mesh origin to (0,0,0).
         *
         * @post Primary and secondary nodes are built and their neighbour lists
         *       populated. Slowness is **not** set.
         * @note The base's raypath method is hard-wired to **1**
         *       (ttcr::LS_SO, second-order least squares) rather than taken from
         *       ttcr::input_parameters::raypath_method, so that setting has no
         *       effect on this solver — unlike ttcr::Grid3Ducfs and
         *       ttcr::Grid3Ducdsp, which pass it through.
         */
        Grid3Ducsp(const std::vector<sxyz<T1>>& no,
                   const std::vector<tetrahedronElem<T2>>& tet,
                   const int ns, const bool rptt, const T1 md,
                   const size_t nt=1, const bool _translateOrigin=false) :
        Grid3Duc<T1,T2,Node3Dcsp<T1,T2>,CELL>(no, tet, 1, rptt, md, nt, _translateOrigin)
        {
            this->buildGridNodes(no, ns, nt);
            this->template buildGridNeighbors<Node3Dcsp<T1,T2>>(this->nodes);
        }

        ~Grid3Ducsp() {
        }

        void raytrace(const std::vector<sxyz<T1>>&,
                      const std::vector<T1>&,
                      const std::vector<sxyz<T1>>&,
                      std::vector<T1>&,
                      const size_t=0) const;

        void raytrace(const std::vector<sxyz<T1>>&,
                      const std::vector<T1>&,
                      const std::vector<std::vector<sxyz<T1>>>&,
                      std::vector<std::vector<T1>*>&,
                      const size_t=0) const;

        void raytrace(const std::vector<sxyz<T1>>&,
                      const std::vector<T1>& ,
                      const std::vector<sxyz<T1>>&,
                      std::vector<T1>&,
                      std::vector<std::vector<sxyz<T1>>>&,
                      const size_t=0) const;

        void raytrace(const std::vector<sxyz<T1>>&,
                      const std::vector<T1>&,
                      const std::vector<std::vector<sxyz<T1>>>&,
                      std::vector<std::vector<T1>*>&,
                      std::vector<std::vector<std::vector<sxyz<T1>>>*>&,
                      const size_t=0) const;

        /// @name Traveltimes and per-cell sensitivities
        ///
        /// The cells report the sensitivity of the traveltime to their medium
        /// parameters, one value per parameter, so the container widens with
        /// the number of parameters the medium has: @ref siv for an isotropic
        /// cell, @ref siv2 or @ref siv4 for the anisotropic ones.  All of
        /// these walk the same shortest-path tree, through
        /// raytraceSensitivity().
        /// @{
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<sxyz<T1>>& Rx,
                      std::vector<T1>& traveltimes,
                      std::vector<std::vector<sxyz<T1>>>& r_data,
                      std::vector<std::vector<siv<T1>>>& l_data,
                      const size_t threadNo=0) const {
            raytraceSensitivity(Tx, t0, Rx, traveltimes, &r_data, l_data, threadNo);
        }

        void raytrace(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<sxyz<T1>>& Rx,
                      std::vector<T1>& traveltimes,
                      std::vector<std::vector<sxyz<T1>>>& r_data,
                      std::vector<std::vector<siv2<T1>>>& l_data,
                      const size_t threadNo=0) const {
            raytraceSensitivity(Tx, t0, Rx, traveltimes, &r_data, l_data, threadNo);
        }

        void raytrace(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<sxyz<T1>>& Rx,
                      std::vector<T1>& traveltimes,
                      std::vector<std::vector<sxyz<T1>>>& r_data,
                      std::vector<std::vector<siv4<T1>>>& l_data,
                      const size_t threadNo=0) const {
            raytraceSensitivity(Tx, t0, Rx, traveltimes, &r_data, l_data, threadNo);
        }

        /// The same, without the raypaths.
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<sxyz<T1>>& Rx,
                      std::vector<T1>& traveltimes,
                      std::vector<std::vector<siv<T1>>>& l_data,
                      const size_t threadNo=0) const {
            raytraceSensitivity(Tx, t0, Rx, traveltimes, nullptr, l_data, threadNo);
        }

        void raytrace(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<sxyz<T1>>& Rx,
                      std::vector<T1>& traveltimes,
                      std::vector<std::vector<siv2<T1>>>& l_data,
                      const size_t threadNo=0) const {
            raytraceSensitivity(Tx, t0, Rx, traveltimes, nullptr, l_data, threadNo);
        }

        void raytrace(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<sxyz<T1>>& Rx,
                      std::vector<T1>& traveltimes,
                      std::vector<std::vector<siv4<T1>>>& l_data,
                      const size_t threadNo=0) const {
            raytraceSensitivity(Tx, t0, Rx, traveltimes, nullptr, l_data, threadNo);
        }
        /// @}

    private:

        /**
         * @brief Traveltimes and per-cell sensitivities, whatever the container
         *
         * Walks the shortest-path tree from each receiver back to the source,
         * asking the cells for the sensitivity of the traveltime of every
         * segment and accumulating it.  The six public overloads differ only
         * in the container the cells report into, and in whether the raypath
         * is wanted, so they share this one body.  Two cases have no
         * counterpart on a rectilinear grid and are kept: a receiver lying on
         * a source, and a receiver whose parent is itself a source node.
         *
         * @param Tx          source coordinates
         * @param t0          source excitation times
         * @param Rx          receiver coordinates
         * @param traveltimes traveltimes at the receivers
         * @param r_data      raypaths, or nullptr if they are not wanted
         * @param l_data      per-cell sensitivities
         * @param threadNo    thread on which to run
         */
        template<typename SIV>
        void raytraceSensitivity(const std::vector<sxyz<T1>>& Tx,
                                 const std::vector<T1>& t0,
                                 const std::vector<sxyz<T1>>& Rx,
                                 std::vector<T1>& traveltimes,
                                 std::vector<std::vector<sxyz<T1>>>* r_data,
                                 std::vector<std::vector<SIV>>& l_data,
                                 const size_t threadNo) const;

        /// @brief Add the contribution of a segment to the cell it crosses
        ///
        /// A ray may cross the same cell more than once, passing through the
        /// secondary nodes along an edge, so the contributions are summed.
        template<typename SIV>
        static void accumulate(std::vector<SIV>& row, const SIV& cell) {
            for ( size_t nc=0; nc<row.size(); ++nc ) {
                if ( row[nc].i == cell.i ) {
                    row[nc] += cell;
                    return;
                }
            }
            row.push_back( cell );
        }

        void initQueue(const std::vector<sxyz<T1>>& Tx,
                       const std::vector<T1>& t0,
                       std::priority_queue<Node3Dcsp<T1,T2>*,
                       std::vector<Node3Dcsp<T1,T2>*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<Node3Dcsp<T1,T2>>& txNodes,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       const size_t threadNo) const;

        void prepropagate(const Node3Dcsp<T1,T2>& node,
                          std::priority_queue<Node3Dcsp<T1,T2>*, std::vector<Node3Dcsp<T1,T2>*>,
                          CompareNodePtr<T1>>& queue,
                          std::vector<bool>& inQueue,
                          std::vector<bool>& frozen,
                          size_t threadNo) const;

        void propagate(std::priority_queue<Node3Dcsp<T1,T2>*,
                       std::vector<Node3Dcsp<T1,T2>*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       const size_t threadNo) const;

    };

    template<typename T1, typename T2, typename CELL>
    void Grid3Ducsp<T1,T2,CELL>::raytrace(const std::vector<sxyz<T1>>& _Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<sxyz<T1>>& _Rx,
                                     std::vector<T1>& traveltimes,
                                     const size_t threadNo) const {

        std::vector<sxyz<T1>> Tx = _Tx;
        std::vector<sxyz<T1>> Rx = _Rx;
        if ( this->translateOrigin ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                Tx[n] -= this->origin;
            }
            for ( size_t n=0; n<Rx.size(); ++n ) {
                Rx[n] -= this->origin;
            }
        }

        this->checkPts(Tx, true);
        this->checkPts(Rx, true);

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

        for (size_t n=0; n<Rx.size(); ++n) {
            traveltimes[n] = this->getTraveltime(Rx[n], this->nodes, threadNo);
        }
    }

    template<typename T1, typename T2, typename CELL>
    void Grid3Ducsp<T1,T2,CELL>::raytrace(const std::vector<sxyz<T1>>& _Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<std::vector<sxyz<T1>>>& _Rx,
                                     std::vector<std::vector<T1>*>& traveltimes,
                                     const size_t threadNo) const {

        std::vector<sxyz<T1>> Tx = _Tx;
        std::vector<std::vector<sxyz<T1>>> Rx = _Rx;
        if ( this->translateOrigin ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                Tx[n] -= this->origin;
            }
            for ( size_t n=0; n<Rx.size(); ++n ) {
                for ( size_t nn=0; nn<Rx[n].size(); ++nn ) {
                    Rx[n][nn] =  _Rx[n][nn] - this->origin;
                }
            }
        }

        this->checkPts(Tx, true);
        for ( size_t n=0; n<Rx.size(); ++n ) {
            this->checkPts(Rx[n], true);
        }

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

        for (size_t nr=0; nr<Rx.size(); ++nr) {
            traveltimes[nr]->resize( Rx[nr].size() );
            for (size_t n=0; n<Rx[nr].size(); ++n) {
                (*traveltimes[nr])[n] = this->getTraveltime(Rx[nr][n], this->nodes, threadNo);
            }
        }
    }

    template<typename T1, typename T2, typename CELL>
    void Grid3Ducsp<T1,T2,CELL>::raytrace(const std::vector<sxyz<T1>>& _Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<sxyz<T1>>& _Rx,
                                     std::vector<T1>& traveltimes,
                                     std::vector<std::vector<sxyz<T1>>>& r_data,
                                     const size_t threadNo) const {

        std::vector<sxyz<T1>> Tx = _Tx;
        std::vector<sxyz<T1>> Rx = _Rx;
        if ( this->translateOrigin ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                Tx[n] -= this->origin;
            }
            for ( size_t n=0; n<Rx.size(); ++n ) {
                Rx[n] -= this->origin;
            }
        }

        this->checkPts(Tx, true);
        this->checkPts(Rx, true);

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
        for ( size_t ni=0; ni<r_data.size(); ++ni ) {
            r_data[ni].resize( 0 );
        }
        T2 nodeParentRx;
        T2 cellParentRx;

        for (size_t n=0; n<Rx.size(); ++n) {

            traveltimes[n] = this->getTraveltime(Rx[n], this->nodes, nodeParentRx, cellParentRx,
                                                 threadNo);

            bool flag=false;
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                if ( Rx[n] == Tx[ns] ) {

                    r_data[n].resize( 1 );
                    r_data[n][0] = Rx[n];

                    flag = true;
                    break;
                }
            }
            if ( flag ) continue;
            for ( size_t ns=0; ns<txNodes.size(); ++ns ) {
                if ( nodeParentRx == txNodes[ns].getGridIndex() ) {
                    // insert Tx at begining
                    r_data[n].push_back(sxyz<T1>(txNodes[ns]));
                    r_data[n].push_back(Rx[n]);
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
            child = Rx[n];
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
            r_data[n].resize( r_tmp.size() );
            for ( size_t nn=0; nn<r_data[n].size(); ++nn ) {
                r_data[n][nn] = r_tmp[ iParent-1-nn ];
            }
        }
        if ( this->translateOrigin ) {
            for (size_t n=0; n<r_data.size(); ++n) {
                for (size_t nn=0; nn<r_data[n].size(); ++nn) {
                    r_data[n][nn] += this->origin;
                }
            }
        }
    }

    template<typename T1, typename T2, typename CELL>
    void Grid3Ducsp<T1,T2,CELL>::raytrace(const std::vector<sxyz<T1>>& _Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<std::vector<sxyz<T1>>>& _Rx,
                                     std::vector<std::vector<T1>*>& traveltimes,
                                     std::vector<std::vector<std::vector<sxyz<T1>>>*>& r_data,
                                     const size_t threadNo) const {

        std::vector<sxyz<T1>> Tx = _Tx;
        std::vector<std::vector<sxyz<T1>>> Rx = _Rx;
        if ( this->translateOrigin ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                Tx[n] -= this->origin;
            }
            for ( size_t n=0; n<Rx.size(); ++n ) {
                for ( size_t nn=0; nn<Rx[n].size(); ++nn ) {
                    Rx[n][nn] =  _Rx[n][nn] - this->origin;
                }
            }
        }

        this->checkPts(Tx, true);
        for ( size_t n=0; n<Rx.size(); ++n ) {
            this->checkPts(Rx[n], true);
        }

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

            traveltimes[nr]->resize( Rx[nr].size() );
            r_data[nr]->resize( Rx[nr].size() );
            for ( size_t ni=0; ni<r_data[nr]->size(); ++ni ) {
                (*r_data[nr])[ni].resize( 0 );
            }

            T2 nodeParentRx;
            T2 cellParentRx;

            for (size_t n=0; n<Rx[nr].size(); ++n) {

                (*traveltimes[nr])[n] = this->getTraveltime(Rx[nr][n], this->nodes,
                                                            nodeParentRx, cellParentRx,
                                                            threadNo);

                bool flag=false;
                for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                    if ( Rx[nr][n] == Tx[ns] ) {

                        (*r_data[nr])[n].resize( 1 );
                        (*r_data[nr])[n][0] = Rx[nr][n];

                        flag = true;
                        break;
                    }
                }
                if ( flag ) continue;
                for ( size_t ns=0; ns<txNodes.size(); ++ns ) {
                    if ( nodeParentRx == txNodes[ns].getGridIndex() ) {
                        // insert Tx at begining
                        (*r_data[nr])[n].push_back(sxyz<T1>(txNodes[ns]));
                        (*r_data[nr])[n].push_back(Rx[nr][n]);
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
                child = Rx[nr][n];
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
        if ( this->translateOrigin ) {
            for (size_t n=0; n<r_data.size(); ++n) {
                for (size_t nn=0; nn<r_data[n]->size(); ++nn) {
                    for (size_t nnn=0; nnn<(*r_data[n])[nn].size(); ++nnn) {
                        (*r_data[n])[nn][nnn] += this->origin;
                    }
                }
            }
        }
    }

    template<typename T1, typename T2, typename CELL>
    void Grid3Ducsp<T1,T2,CELL>::initQueue(const std::vector<sxyz<T1>>& Tx,
                                      const std::vector<T1>& t0,
                                      std::priority_queue<Node3Dcsp<T1,T2>*,
                                      std::vector<Node3Dcsp<T1,T2>*>,
                                      CompareNodePtr<T1>>& queue,
                                      std::vector<Node3Dcsp<T1,T2>>& txNodes,
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
                // If Tx[n] is not on a node, we create a new node and initialize the queue:
                txNodes.push_back( Node3Dcsp<T1,T2>(t0[n], Tx[n].x, Tx[n].y, Tx[n].z,
                                                    this->nThreads, threadNo));
                txNodes.back().pushOwner( this->getCellNo(Tx[n]) );
                txNodes.back().setGridIndex( static_cast<T2>(this->nodes.size()+
                                                             txNodes.size()-1) );
                frozen.push_back( true );

                prepropagate(txNodes.back(), queue, inQueue, frozen, threadNo); // See description in the function declaration

                //	queue.push( &(txNodes.back()) );	//Don't use if prepropagate is used
                //	inQueue.push_back( true );			//Don't use if prepropagate is used

            }
        }
    }


    template<typename T1, typename T2, typename CELL>
    void Grid3Ducsp<T1,T2,CELL>::prepropagate(const Node3Dcsp<T1,T2>& node,
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
                T1 dt = this->computeDt(node, this->nodes[neibNo], cellNo);

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
    void Grid3Ducsp<T1,T2,CELL>::propagate(std::priority_queue<Node3Dcsp<T1,T2>*,
                                      std::vector<Node3Dcsp<T1,T2>*>,
                                      CompareNodePtr<T1>>& queue,
                                      std::vector<bool>& inQueue,
                                      std::vector<bool>& frozen,
                                      const size_t threadNo) const {

        while ( !queue.empty() ) {
            const Node3Dcsp<T1,T2>* src = queue.top();
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
                    T1 dt = this->computeDt(*src, this->nodes[neibNo], cellNo);

                    if (src->getTT(threadNo)+dt < this->nodes[neibNo].getTT(threadNo)) {
                        this->nodes[neibNo].setTT( src->getTT(threadNo)+dt, threadNo );
                        this->nodes[neibNo].setnodeParent(src->getGridIndex(),threadNo);
                        this->nodes[neibNo].setCellParent(cellNo, threadNo );

                        if ( !inQueue[neibNo] ) {
                            queue.push( &(this->nodes[neibNo]) );
                            inQueue[neibNo] = true;
                        }
                    }
                }
            }
        }
    }

    template<typename T1, typename T2, typename CELL>
    template<typename SIV>
    void Grid3Ducsp<T1,T2,CELL>::raytraceSensitivity(const std::vector<sxyz<T1>>& _Tx,
                                                     const std::vector<T1>& t0,
                                                     const std::vector<sxyz<T1>>& _Rx,
                                                     std::vector<T1>& traveltimes,
                                                     std::vector<std::vector<sxyz<T1>>>* r_data,
                                                     std::vector<std::vector<SIV>>& l_data,
                                                     const size_t threadNo) const {

        if constexpr ( !cell_reports_into<CELL, SIV, Node3Dcsp<T1,T2>, sxyz<T1>>::value ) {
            throw std::logic_error("Error: these cells do not report their "
                                   "sensitivity in the container requested.");
        } else {

        std::vector<sxyz<T1>> Tx = _Tx;
        std::vector<sxyz<T1>> Rx = _Rx;
        if ( this->translateOrigin ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                Tx[n] -= this->origin;
            }
            for ( size_t n=0; n<Rx.size(); ++n ) {
                Rx[n] -= this->origin;
            }
        }

        this->checkPts(Tx, true);
        this->checkPts(Rx, true);

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
        if ( l_data.size() != Rx.size() ) {
            l_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<l_data.size(); ++ni ) {
            l_data[ni].resize( 0 );
        }
        if ( r_data != nullptr ) {
            if ( r_data->size() != Rx.size() ) {
                r_data->resize( Rx.size() );
            }
            for ( size_t ni=0; ni<r_data->size(); ++ni ) {
                (*r_data)[ni].resize( 0 );
            }
        }

        T2 nodeParentRx;
        T2 cellParentRx;

        for (size_t n=0; n<Rx.size(); ++n) {

            traveltimes[n] = this->getTraveltime(Rx[n], this->nodes, nodeParentRx,
                                                 cellParentRx, threadNo);

            // the receiver sits on a source: the ray has no length
            bool flag=false;
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                if ( Rx[n] == Tx[ns] ) {
                    if ( r_data != nullptr ) {
                        (*r_data)[n].resize( 1 );
                        (*r_data)[n][0] = Rx[n];
                    }
                    // no need to update l_data: ray length is zero
                    flag = true;
                }
            }
            if ( flag ) continue;

            SIV cell;

            // the receiver hangs directly off a source node
            for ( size_t ns=0; ns<txNodes.size(); ++ns ) {
                if ( nodeParentRx == txNodes[ns].getGridIndex() ) {
                    if ( r_data != nullptr ) {
                        // insert Tx at begining
                        (*r_data)[n].push_back(sxyz<T1>(txNodes[ns]));
                        (*r_data)[n].push_back(Rx[n]);
                    }
                    cell.i = cellParentRx;
                    this->cells.computeDistance( txNodes[ns], Rx[n], cell );
                    l_data[n].push_back( cell );
                    flag = true;
                    break;
                }
            }
            if ( flag ) continue;

            // Rx are in nodes (not txNodes)
            std::vector<Node3Dcsp<T1,T2>> *node_p;
            node_p = &this->nodes;

            std::vector<sxyz<T1>> r_tmp;
            T2 iChild, iParent = nodeParentRx;
            sxyz<T1> child;

            // store the son's coord
            child = Rx[n];
            cell.i = cellParentRx;
            while ( (*node_p)[iParent].getNodeParent(threadNo) !=
                   std::numeric_limits<T2>::max() ) {

                if ( r_data != nullptr ) r_tmp.push_back( child );

                this->cells.computeDistance( (*node_p)[iParent], child, cell );
                accumulate( l_data[n], cell );

                // we now go up in time - parent becomes the child of grand'pa
                iChild = iParent;
                child = (*node_p)[iChild];
                cell.i = (*node_p)[iChild].getCellParent(threadNo);

                // grand'pa is now papa
                iParent = (*node_p)[iChild].getNodeParent(threadNo);
                if ( iParent >= this->nodes.size() ) {
                    node_p = &txNodes;
                    iParent -= this->nodes.size();
                }
                else {
                    node_p = &this->nodes;
                }
            }

            // parent is now at Tx
            if ( r_data != nullptr ) r_tmp.push_back( child );

            this->cells.computeDistance( (*node_p)[iParent], child, cell );
            accumulate( l_data[n], cell );

            // finally, store Tx position
            child = (*node_p)[iParent];
            if ( r_data != nullptr ) r_tmp.push_back( child );

            //  must be sorted to build matrix L
            std::sort(l_data[n].begin(), l_data[n].end(),
                      [](const SIV& a, const SIV& b) { return a.i < b.i; });

            if ( r_data != nullptr ) {
                // the order should be from Tx to Rx, so we reorder...
                iParent = static_cast<T2>(r_tmp.size());
                (*r_data)[n].resize( r_tmp.size() );
                for ( size_t nn=0; nn<(*r_data)[n].size(); ++nn ) {
                    (*r_data)[n][nn] = r_tmp[ iParent-1-nn ];
                }
            }
        }

        if ( this->translateOrigin && r_data != nullptr ) {
            for ( size_t n=0; n<r_data->size(); ++n ) {
                for ( size_t nn=0; nn<(*r_data)[n].size(); ++nn ) {
                    (*r_data)[n][nn] += this->origin;
                }
            }
        }
        }
    }

}

#endif

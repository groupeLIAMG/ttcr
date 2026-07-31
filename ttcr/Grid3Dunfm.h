//
//  Grid3Dunfm.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-04-21.
//  Copyright (c) 2014 Bernard Giroux. All rights reserved.
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
 * @file Grid3Dunfm.h
 * @brief Fast-marching solver on a 3-D tetrahedral mesh with node slowness.
 *
 * Declares ttcr::Grid3Dunfm, which solves the eikonal equation by advancing a
 * narrow band outwards from the sources in order of increasing traveltime
 * (Sethian, 1996). It is the node-slowness counterpart of ttcr::Grid3Ducfm and
 * the 3-D counterpart of ttcr::Grid2Dunfm.
 *
 * @section g3dunfm_vs_fs Marching versus sweeping
 * Both solve the same equation on the same mesh; they differ in how they order
 * the work.
 *
 * - **Marching (here)** — a priority queue holds the band of nodes adjacent to
 *   the solved region. The one of least traveltime is frozen and its neighbours
 *   updated. One pass suffices, since a node is only frozen once every path
 *   that could reach it has been considered.
 * - **Sweeping (ttcr::Grid3Dunfs)** — fixed node orderings are swept
 *   repeatedly until the field stops changing. No queue, but several passes.
 *
 * Marching is the more predictable of the two — a single pass, no convergence
 * tolerance to choose — at the cost of maintaining the queue.
 *
 * @sa Grid3Dun.h, Grid3Ducfm.h, Grid2Dunfm.h, Grid3Dunfs.h
 */

#ifndef ttcr_Grid3Dunfm_h
#define ttcr_Grid3Dunfm_h

#include <cmath>
#include <fstream>
#include <queue>
#include <vector>

#include "Grid3Dun.h"
#include "Node3Dn.h"

namespace ttcr {

    /**
     * @brief Fast-marching traveltimes on a tetrahedral mesh with node slowness.
     *
     * @tparam T1 floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2 unsigned integer type used for node and cell indices.
     *
     * Uses ttcr::Node3Dn and inserts no secondary nodes; the mesh connectivity
     * is all the algorithm needs.
     *
     * @sa Grid3Dun.h, Grid3Dunfs.h, Grid3Ducfm.h
     */
    template<typename T1, typename T2>
    class Grid3Dunfm : public Grid3Dun<T1,T2,Node3Dn<T1,T2>> {
    public:
        /**
         * @brief Build the mesh.
         * @param no   primary node coordinates.
         * @param tet  tetrahedra, as quadruples of indices into @p no.
         * @param rp   raypath gradient estimator; see @ref g3dun_raypath.
         * @param iv   interpolate velocity rather than slowness.
         * @param rptt compute traveltimes by integrating along the raypath.
         * @param md   minimum distance for snapping a raypath point onto a node.
         * @param nt   number of threads.
         * @param _translateOrigin shift the mesh to the origin.
         *
         */
        Grid3Dunfm(const std::vector<sxyz<T1>>& no,
                   const std::vector<tetrahedronElem<T2>>& tet,
                   const int rp, const bool iv, const bool rptt, const T1 md,
                   const size_t nt=1, const bool _translateOrigin=false) :
        Grid3Dun<T1,T2,Node3Dn<T1,T2>>(no, tet, rp, iv, rptt, md, nt, _translateOrigin)
        {
            this->buildGridNodes(no, nt);
            this->template buildGridNeighbors<Node3Dn<T1,T2>>(this->nodes);
        }

        ~Grid3Dunfm() {
        }

    private:

        /**
         * @brief Seed the narrow band from the sources.
         * @param Tx source coordinates.
         * @param t0 source excitation times.
         * @param[out] narrow_band the priority queue, ordered by traveltime.
         * @param[out] inQueue     whether each node is currently in the band.
         * @param[out] frozen      whether each node's traveltime is final.
         * @param threadNo         thread number.
         *
         * Sets the source nodes' traveltimes, freezes them, and puts their
         * neighbours into the band. How far the seeding reaches is governed by
         * @c source_radius; see ttcr::Grid3Dun::setSourceRadius.
         */
        void initBand(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      std::priority_queue<Node3Dn<T1,T2>*,
                      std::vector<Node3Dn<T1,T2>*>,
                      CompareNodePtr<T1>>& narrow_band,
                      std::vector<bool>& inQueue,
                      std::vector<bool>& frozen,
                      const size_t threadNo) const;

        /**
         * @brief Advance the narrow band until it is empty.
         * @param[in,out] narrow_band the priority queue.
         * @param[in,out] inQueue     band membership flags.
         * @param[in,out] frozen      finalised flags.
         * @param threadNo            thread number.
         *
         * Takes the band node of least traveltime, freezes it, updates its
         * unfrozen neighbours through ttcr::Grid3Dun::localUpdate3D, and adds
         * any newly reachable node to the band.
         */
        void propagate(std::priority_queue<Node3Dn<T1,T2>*,
                       std::vector<Node3Dn<T1,T2>*>,
                       CompareNodePtr<T1>>& narrow_band,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       const size_t threadNo) const;

        /**
         * @brief Solve the traveltime field by marching.
         * @param Tx       source coordinates, already origin-translated.
         * @param t0       source excitation times.
         * @param Rx       receiver coordinates; checked to lie in the mesh.
         * @param threadNo thread number.
         *
         * Overrides the base class's solver hook, which the public @c raytrace
         * overloads call; see ttcr::Grid3D.
         */
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<sxyz<T1>>& Rx,
                      const size_t threadNo=0) const;

        /**
         * @brief Solver hook, grouped-receiver form.
         * @param Tx       source coordinates.
         * @param t0       source excitation times.
         * @param Rx       receiver coordinates, one group per source.
         * @param threadNo thread number.
         */
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      const size_t threadNo=0) const;

    };

    template<typename T1, typename T2>
    void Grid3Dunfm<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<sxyz<T1>>& Rx,
                                     const size_t threadNo) const {

        this->checkPts(Tx, true);
        this->checkPts(Rx, true);

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node3Dn<T1,T2>*, std::vector<Node3Dn<T1,T2>*>,
        CompareNodePtr<T1>> narrow_band( cmp );

        std::vector<bool> inQueue( this->nodes.size(), false );
        std::vector<bool> frozen( this->nodes.size(), false );

        initBand(Tx, t0, narrow_band, inQueue, frozen, threadNo);

        propagate(narrow_band, inQueue, frozen, threadNo);
    }

    template<typename T1, typename T2>
    void Grid3Dunfm<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
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
        std::priority_queue< Node3Dn<T1,T2>*, std::vector<Node3Dn<T1,T2>*>,
        CompareNodePtr<T1>> narrow_band( cmp );

        std::vector<bool> inBand( this->nodes.size(), false );
        std::vector<bool> frozen( this->nodes.size(), false );

        initBand(Tx, t0, narrow_band, inBand, frozen, threadNo);

        propagate(narrow_band, inBand, frozen, threadNo);
    }

    template<typename T1, typename T2>
    void Grid3Dunfm<T1,T2>::initBand(const std::vector<sxyz<T1>>& Tx,
                                     const std::vector<T1>& t0,
                                     std::priority_queue<Node3Dn<T1,T2>*,
                                     std::vector<Node3Dn<T1,T2>*>,
                                     CompareNodePtr<T1>>& narrow_band,
                                     std::vector<bool>& inBand,
                                     std::vector<bool>& frozen,
                                     const size_t threadNo) const {

        for (size_t n=0; n<Tx.size(); ++n) {
            bool found = false;
            for ( size_t nn=0; nn<this->nodes.size(); ++nn ) {
                if ( this->nodes[nn] == Tx[n] ) {
                    found = true;
                    this->nodes[nn].setTT( t0[n], threadNo );
                    narrow_band.push( &(this->nodes[nn]) );
                    inBand[nn] = true;
                    frozen[nn] = true;

                    if ( Tx.size()==1 ) {
                        if ( Grid3Dun<T1,T2,Node3Dn<T1,T2>>::source_radius == 0.0 ) {
                            // populate around Tx
                            for ( size_t no=0; no<this->nodes[nn].getOwners().size(); ++no ) {

                                T2 cellNo = this->nodes[nn].getOwners()[no];
                                for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                                    T2 neibNo = this->neighbors[cellNo][k];
                                    if ( neibNo == nn ) continue;
                                    T1 dt = this->computeDt(this->nodes[nn], this->nodes[neibNo]);

                                    if ( t0[n]+dt < this->nodes[neibNo].getTT(threadNo) ) {
                                        this->nodes[neibNo].setTT( t0[n]+dt, threadNo );

                                        if ( !inBand[neibNo] ) {
                                            narrow_band.push( &(this->nodes[neibNo]) );
                                            inBand[neibNo] = true;
                                            frozen[neibNo] = true;
                                        }
                                    }
                                }
                            }
                        } else {

                            // find nodes within source radius
                            size_t nodes_added = 0;
                            for ( size_t no=0; no<this->nodes.size(); ++no ) {

                                if ( no == nn ) continue;

                                T1 d = this->nodes[nn].getDistance( this->nodes[no] );
                                if ( d <= Grid3Dun<T1,T2,Node3Dn<T1,T2>>::source_radius ) {

                                    T1 dt = this->computeDt(this->nodes[nn], this->nodes[no] );

                                    if ( t0[n]+dt < this->nodes[no].getTT(threadNo) ) {
                                        this->nodes[no].setTT( t0[n]+dt, threadNo );

                                        if ( !inBand[no] ) {
                                            narrow_band.push( &(this->nodes[no]) );
                                            inBand[no] = true;
                                            frozen[no] = true;
                                            nodes_added++;
                                        }
                                    }
                                }
                            }
                            if ( nodes_added == 0 ) {
                                std::cerr << "Error: no nodes found within source radius, aborting" << std::endl;
                                abort();
                            } else {
                                std::cout << "(found " << nodes_added << " nodes around Tx point)\n";
                            }
                        }

                    }

                    break;
                }
            }
            if ( found==false ) {

                T2 sTx = this->computeSlowness( Tx[n], true );

                T2 cellNo = this->getCellNo(Tx[n]);
                if ( Grid3Dun<T1,T2,Node3Dn<T1,T2>>::source_radius == 0.0 ) {
                    for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                        T2 neibNo = this->neighbors[cellNo][k];

                        // compute dt
                        T1 dt = this->computeDt(this->nodes[neibNo], Tx[n], sTx);

                        this->nodes[neibNo].setTT( t0[n]+dt, threadNo );
                        narrow_band.push( &(this->nodes[neibNo]) );
                        inBand[neibNo] = true;
                        frozen[neibNo] = true;

                    }
                } else if ( Tx.size()==1 ) { // look into source radius only for point sources

                    // find nodes within source radius
                    size_t nodes_added = 0;
                    for ( size_t no=0; no<this->nodes.size(); ++no ) {

                        T1 d = this->nodes[no].getDistance( Tx[n] );
                        if ( d <= Grid3Dun<T1,T2,Node3Dn<T1,T2>>::source_radius ) {

                            T1 dt = this->computeDt(this->nodes[no], Tx[n], sTx);

                            if ( t0[n]+dt < this->nodes[no].getTT(threadNo) ) {
                                this->nodes[no].setTT( t0[n]+dt, threadNo );

                                if ( !inBand[no] ) {
                                    narrow_band.push( &(this->nodes[no]) );
                                    inBand[no] = true;
                                    frozen[no] = true;
                                    nodes_added++;
                                }
                            }
                        }
                    }
                    if ( nodes_added == 0 ) {
                        std::cerr << "Error: no nodes found within source radius, aborting" << std::endl;
                        abort();
                    } else {
                        std::cout << "(found " << nodes_added << " nodes around Tx point)\n";
                    }
                }
            }
        }
    }

    template<typename T1, typename T2>
    void Grid3Dunfm<T1,T2>::propagate(std::priority_queue<Node3Dn<T1,T2>*,
                                      std::vector<Node3Dn<T1,T2>*>,
                                      CompareNodePtr<T1>>& narrow_band,
                                      std::vector<bool>& inNarrowBand,
                                      std::vector<bool>& frozen,
                                      const size_t threadNo) const {

        while ( !narrow_band.empty() ) {

            const Node3Dn<T1,T2>* source = narrow_band.top();
            narrow_band.pop();
            inNarrowBand[ source->getGridIndex() ] = false;
            frozen[ source->getGridIndex() ] = true;   // marked as known

            for ( size_t no=0; no<source->getOwners().size(); ++no ) {

                T2 cellNo = source->getOwners()[no];

                for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                    T2 neibNo = this->neighbors[cellNo][k];
                    if ( neibNo == source->getGridIndex() || frozen[neibNo] ) {
                        continue;
                    }

                    this->localUpdate3D( &(this->nodes[neibNo]), threadNo );

                    if ( !inNarrowBand[neibNo] ) {
                        narrow_band.push( &(this->nodes[neibNo]) );
                        inNarrowBand[neibNo] = true;
                    }
                }
            }
        }
    }

}

#endif

//
//  Grid3Dunsp.h
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

#ifndef ttcr_Grid3Dunsp_h
#define ttcr_Grid3Dunsp_h

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <stdexcept>
#include <vector>

#include "Grid3Dun.h"
#include "Interpolator.h"
#include "Node3Dnsp.h"
#include "utils.h"

namespace ttcr {

    template<typename T1, typename T2>
    class Grid3Dunsp : public Grid3Dun<T1,T2,Node3Dnsp<T1,T2>> {
    public:
        Grid3Dunsp(const std::vector<sxyz<T1>>& no,
                   const std::vector<tetrahedronElem<T2>>& tet,
                   const int ns, const bool iv, const bool rptt, const T1 md,
                   const size_t nt=1, const bool _translateOrigin=false) :
        Grid3Dun<T1,T2,Node3Dnsp<T1,T2>>(no, tet, 1, iv, rptt, md, nt, _translateOrigin),
        nSecondary(ns)
        {
            this->buildGridNodes(no, ns, nt);
            this->template buildGridNeighbors<Node3Dnsp<T1,T2>>(this->nodes);
        }

        ~Grid3Dunsp() {
        }

        void setSlowness(const std::vector<T1>& s) {
            if ( this->nPrimary != s.size() ) {
                throw std::length_error("Error: slowness vectors of incompatible size.");
            }
            for ( size_t n=0; n<this->nPrimary; ++n ) {
                this->nodes[n].setNodeSlowness( s[n] );
            }
            if ( nSecondary>0 ) {
                if ( this->processVel )
                    this->interpVelocitySecondary(nSecondary);
                else
                    this->interpSlownessSecondary(nSecondary);
            }
        }


        void setSlowness(const T1 *s, const size_t ns) {
            if ( this->nPrimary != ns ) {
                throw std::length_error("Error: slowness vectors of incompatible size.");
            }
            for ( size_t n=0; n<this->nPrimary; ++n ) {
                this->nodes[n].setNodeSlowness( s[n] );
            }
            if ( nSecondary>0 ) {
                if ( this->processVel )
                    this->interpVelocitySecondary(nSecondary);
                else
                    this->interpSlownessSecondary(nSecondary);
            }
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

        void raytrace(const std::vector<sxyz<T1>>&,
                      const std::vector<T1>& ,
                      const std::vector<sxyz<T1>>&,
                      std::vector<T1>&,
                      std::vector<std::vector<sxyz<T1>>>&,
                      std::vector<std::vector<siv<T1>>>&,
                      const size_t=0) const;


    private:
        T2 nSecondary;

        void initQueue(const std::vector<sxyz<T1>>& Tx,
                       const std::vector<T1>& t0,
                       std::priority_queue<Node3Dnsp<T1,T2>*,
                       std::vector<Node3Dnsp<T1,T2>*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<Node3Dnsp<T1,T2>>& txNodes,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       const size_t threadNo) const;

        void prepropagate(const Node3Dnsp<T1,T2>& node,
                          std::priority_queue<Node3Dnsp<T1,T2>*, std::vector<Node3Dnsp<T1,T2>*>,
                          CompareNodePtr<T1>>& queue,
                          std::vector<bool>& inQueue,
                          std::vector<bool>& frozen,
                          size_t threadNo) const;

        void propagate(std::priority_queue<Node3Dnsp<T1,T2>*,
                       std::vector<Node3Dnsp<T1,T2>*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       const size_t threadNo) const;

        T1 getTraveltime(const sxyz<T1>& Rx,
                         const std::vector<Node3Dnsp<T1,T2>>& nodes,
                         const size_t threadNo) const;

        T1 getTraveltime(const sxyz<T1>& Rx,
                         const std::vector<Node3Dnsp<T1,T2>>& nodes,
                         T2& nodeParentRx,
                         T2& cellParentRx,
                         const size_t threadNo) const;
    };

    template<typename T1, typename T2>
    void Grid3Dunsp<T1,T2>::raytrace(const std::vector<sxyz<T1>>& _Tx,
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

        this->checkPts(Tx);
        this->checkPts(Rx);

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

        for (size_t n=0; n<Rx.size(); ++n) {
            traveltimes[n] = this->getTraveltime(Rx[n], this->nodes, threadNo);
        }
    }

    template<typename T1, typename T2>
    void Grid3Dunsp<T1,T2>::raytrace(const std::vector<sxyz<T1>>& _Tx,
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

        this->checkPts(Tx);
        for ( size_t n=0; n<Rx.size(); ++n ) {
            this->checkPts(Rx[n]);
        }

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

        for (size_t nr=0; nr<Rx.size(); ++nr) {
            traveltimes[nr]->resize( Rx[nr].size() );
            for (size_t n=0; n<Rx[nr].size(); ++n) {
                (*traveltimes[nr])[n] = this->getTraveltime(Rx[nr][n], this->nodes, threadNo);
            }
        }
    }


    template<typename T1, typename T2>
    void Grid3Dunsp<T1,T2>::raytrace(const std::vector<sxyz<T1>>& _Tx,
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

        this->checkPts(Tx);
        this->checkPts(Rx);

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
            std::vector<Node3Dnsp<T1,T2>> *node_p;
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

    template<typename T1, typename T2>
    void Grid3Dunsp<T1,T2>::raytrace(const std::vector<sxyz<T1>>& _Tx,
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

        this->checkPts(Tx);
        for ( size_t n=0; n<Rx.size(); ++n ) {
            this->checkPts(Rx[n]);
        }

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
                std::vector<Node3Dnsp<T1,T2>> *node_p;
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

    template<typename T1, typename T2>
    void Grid3Dunsp<T1,T2>::raytrace(const std::vector<sxyz<T1>>& _Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<sxyz<T1>>& _Rx,
                                     std::vector<T1>& traveltimes,
                                     std::vector<std::vector<sxyz<T1>>>& r_data,
                                     std::vector<std::vector<siv<T1>>>& l_data,
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

        this->checkPts(Tx);
        this->checkPts(Rx);

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
        if ( l_data.size() != Rx.size() ) {
            l_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<l_data.size(); ++ni ) {
            l_data[ni].resize( 0 );
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

            traveltimes[n] = getTraveltime(Rx[n], this->nodes, nodeParentRx, cellParentRx,
                                           threadNo);

            bool flag=false;
            for ( size_t ns=0; ns<Tx.size(); ++n ) {
                if ( Rx[n] == Tx[ns] ) {

                    r_data[n].resize( 1 );
                    r_data[n][0] = Rx[n];

                    // no need to update l_data: ray length is zero

                    flag = true;
                }
            }
            if ( flag ) continue;

            siv<T1> cell;
            for ( size_t ns=0; ns<txNodes.size(); ++ns ) {
                if ( nodeParentRx == txNodes[ns].getGridIndex() ) {
                    // insert Tx at begining
                    r_data[n].push_back(sxyz<T1>(txNodes[ns]));
                    r_data[n].push_back(Rx[n]);
                    cell.i = cellParentRx;
                    cell.v = Rx[n].getDistance(txNodes[ns]);
                    l_data[n].push_back( cell );
                    flag = true;
                    break;
                }
            }
            if ( flag ) continue;

            // Rx are in nodes (not txNodes)
            std::vector<Node3Dnsp<T1,T2>> *node_p;
            node_p = &this->nodes;

            std::vector<sxyz<T1>> r_tmp;
            T2 iChild, iParent = nodeParentRx;
            sxyz<T1> child;

            // store the son's coord
            child = Rx[n];
            cell.i = cellParentRx;
            while ( (*node_p)[iParent].getNodeParent(threadNo) !=
                   std::numeric_limits<T2>::max() ) {

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
            child = (*node_p)[iParent];
            r_tmp.push_back( child );

            //  must be sorted to build matrix L
            std::sort(l_data[n].begin(), l_data[n].end(), CompareSiv_i<T1>());

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

    template<typename T1, typename T2>
    void Grid3Dunsp<T1,T2>::initQueue(const std::vector<sxyz<T1>>& Tx,
                                      const std::vector<T1>& t0,
                                      std::priority_queue<Node3Dnsp<T1,T2>*,
                                      std::vector<Node3Dnsp<T1,T2>*>,
                                      CompareNodePtr<T1>>& queue,
                                      std::vector<Node3Dnsp<T1,T2>>& txNodes,
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
                txNodes.push_back( Node3Dnsp<T1,T2>(t0[n], Tx[n].x, Tx[n].y, Tx[n].z,
                                                    this->nThreads, threadNo));
                T2 cn = this->getCellNo(Tx[n]);
                txNodes.back().pushOwner( cn );
                txNodes.back().setGridIndex( static_cast<T2>(this->nodes.size()+
                                                             txNodes.size()-1) );
                frozen.push_back( true );

                T1 s;
                if ( this->processVel )
                    s = Interpolator<T1>::trilinearTriangleVel(txNodes.back(),
                                                               this->nodes[this->neighbors[cn][0]],
                                                               this->nodes[this->neighbors[cn][1]],
                                                               this->nodes[this->neighbors[cn][2]],
                                                               this->nodes[this->neighbors[cn][3]]);
                else
                    s = Interpolator<T1>::trilinearTriangle(txNodes.back(),
                                                            this->nodes[this->neighbors[cn][0]],
                                                            this->nodes[this->neighbors[cn][1]],
                                                            this->nodes[this->neighbors[cn][2]],
                                                            this->nodes[this->neighbors[cn][3]]);
                txNodes.back().setNodeSlowness(s);


                // prepropagate(txNodes.back(), queue, inQueue, frozen, threadNo); // See description in the function declaration

                queue.push( &(txNodes.back()) );	//Don't use if prepropagate is used
                inQueue.push_back( true );			//Don't use if prepropagate is used

            }
        }
    }


    template<typename T1, typename T2>
    void Grid3Dunsp<T1,T2>::prepropagate(const Node3Dnsp<T1,T2>& node,
                                         std::priority_queue<Node3Dnsp<T1,T2>*,
                                         std::vector<Node3Dnsp<T1,T2>*>,
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
    void Grid3Dunsp<T1,T2>::propagate(std::priority_queue<Node3Dnsp<T1,T2>*,
                                      std::vector<Node3Dnsp<T1,T2>*>,
                                      CompareNodePtr<T1>>& queue,
                                      std::vector<bool>& inQueue,
                                      std::vector<bool>& frozen,
                                      const size_t threadNo) const {

        while ( !queue.empty() ) {
            const Node3Dnsp<T1,T2>* src = queue.top();
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

    template<typename T1, typename T2>
    T1 Grid3Dunsp<T1,T2>::getTraveltime(const sxyz<T1>& Rx,
                                        const std::vector<Node3Dnsp<T1,T2>>& nodes,
                                        const size_t threadNo) const {

        for ( size_t nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == Rx ) {
                return nodes[nn].getTT(threadNo);
            }
        }
        //If Rx is not on a node:
        T1 slo = this->computeSlowness( Rx );

        T2 cellNo = this->getCellNo( Rx );
        T2 neibNo = this->neighbors[cellNo][0];
        T1 dt = this->computeDt(nodes[neibNo], Rx, slo);

        T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
        for ( size_t k=1; k< this->neighbors[cellNo].size(); ++k ) {
            neibNo = this->neighbors[cellNo][k];
            dt = this->computeDt(nodes[neibNo], Rx, slo);
            if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
                traveltime = nodes[neibNo].getTT(threadNo)+dt;
            }
        }
        return traveltime;
    }

    template<typename T1, typename T2>
    T1 Grid3Dunsp<T1,T2>::getTraveltime(const sxyz<T1>& Rx,
                                        const std::vector<Node3Dnsp<T1,T2>>& nodes,
                                        T2& nodeParentRx, T2& cellParentRx,
                                        const size_t threadNo) const {

        for ( size_t nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == Rx ) {
                nodeParentRx = nodes[nn].getNodeParent(threadNo);
                cellParentRx = nodes[nn].getCellParent(threadNo);
                return nodes[nn].getTT(threadNo);
            }
        }
        //If Rx is not on a node:
        T1 slo = this->computeSlowness( Rx );

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
                traveltime = nodes[neibNo].getTT(threadNo)+dt;
                nodeParentRx = neibNo;
            }
        }
        return traveltime;
    }

}

#endif

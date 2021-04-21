//
//  Grid2Ducsp.h
//  ttcr
//
//  Created by Bernard Giroux on 2012-09-20.
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

#ifndef ttcr_Grid2Ducsp_h
#define ttcr_Grid2Ducsp_h

#include <array>
#include <fstream>
#include <queue>

#include "Grid2Duc.h"

namespace ttcr {

    template<typename T1, typename T2, typename NODE, typename S>
    class Grid2Ducsp : public Grid2Duc<T1,T2,NODE,S> {
    public:
        Grid2Ducsp(const std::vector<S>& no,
                   const std::vector<triangleElem<T2>>& tri,
                   const T2 ns, const bool ttrp, const size_t nt=1) :
        Grid2Duc<T1,T2,NODE,S>(no, tri, ttrp, nt)
        {
            this->buildGridNodes(no, ns, nt);
            this->template buildGridNeighbors<NODE>(this->nodes);
        }

        ~Grid2Ducsp() {
        }

        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<S>& Rx,
                      std::vector<T1>& traveltimes,
                      const size_t threadNo=0) const;

        void raytrace(const std::vector<S>&,
                      const std::vector<T1>&,
                      const std::vector<const std::vector<S>*>&,
                      std::vector<std::vector<T1>*>&,
                      const size_t=0) const;

        void raytrace(const std::vector<S>&,
                      const std::vector<T1>& ,
                      const std::vector<S>&,
                      std::vector<T1>&,
                      std::vector<std::vector<S>>&,
                      const size_t=0) const;

        void raytrace(const std::vector<S>&,
                      const std::vector<T1>&,
                      const std::vector<const std::vector<S>*>&,
                      std::vector<std::vector<T1>*>&,
                      std::vector<std::vector<std::vector<S>>*>&,
                      const size_t=0) const;

        void raytrace(const std::vector<S>&,
                      const std::vector<T1>& ,
                      const std::vector<S>&,
                      std::vector<T1>&,
                      std::vector<std::vector<S>>&,
                      std::vector<std::vector<siv<T1>>>&,
                      const size_t=0) const;

    private:

        void initQueue(const std::vector<S>& Tx,
                       const std::vector<T1>& t0,
                       std::priority_queue<NODE*,
                       std::vector<NODE*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<NODE>& txNodes,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       const size_t threadNo) const;

        void propagate(std::priority_queue<NODE*,
                       std::vector<NODE*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       const size_t threadNo) const;

        void raytrace(const std::vector<S>&,
                      const std::vector<T1>&,
                      const std::vector<S>&,
                      const size_t=0) const;
    };

    template<typename T1, typename T2, typename NODE, typename S>
    void Grid2Ducsp<T1,T2,NODE,S>::raytrace(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            const std::vector<S>& Rx,
                                            const size_t threadNo) const {

        this->checkPts(Tx);
        this->checkPts(Rx);

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< NODE*, std::vector<NODE*>,
        CompareNodePtr<T1>> queue( cmp );

        std::vector<NODE> txNodes;
        std::vector<bool> inQueue( this->nodes.size(), false );
        std::vector<bool> frozen( this->nodes.size(), false );

        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);

        propagate(queue, inQueue, frozen, threadNo);
    }

    template<typename T1, typename T2, typename NODE, typename S>
    void Grid2Ducsp<T1,T2,NODE,S>::raytrace(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            const std::vector<S>& Rx,
                                            std::vector<T1>& traveltimes,
                                            const size_t threadNo) const {

        this->checkPts(Tx);
        this->checkPts(Rx);

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< NODE*, std::vector<NODE*>,
        CompareNodePtr<T1>> queue( cmp );

        std::vector<NODE> txNodes;
        std::vector<bool> inQueue( this->nodes.size(), false );
        std::vector<bool> frozen( this->nodes.size(), false );

        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);

        propagate(queue, inQueue, frozen, threadNo);

        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }

        for (size_t n=0; n<Rx.size(); ++n) {
            traveltimes[n] = this->getTraveltime(Rx[n], threadNo);
        }
    }

    template<typename T1, typename T2, typename NODE, typename S>
    void Grid2Ducsp<T1,T2,NODE,S>::raytrace(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            const std::vector<const std::vector<S>*>& Rx,
                                            std::vector<std::vector<T1>*>& traveltimes,
                                            const size_t threadNo) const {

        this->checkPts(Tx);
        for ( size_t n=0; n<Rx.size(); ++n ) {
            this->checkPts(*Rx[n]);
        }

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< NODE*, std::vector<NODE*>,
        CompareNodePtr<T1>> queue( cmp );

        std::vector<NODE> txNodes;
        std::vector<bool> inQueue( this->nodes.size(), false );
        std::vector<bool> frozen( this->nodes.size(), false );

        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);

        propagate(queue, inQueue, frozen, threadNo);

        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }

        for (size_t nr=0; nr<Rx.size(); ++nr) {
            traveltimes[nr]->resize( Rx[nr]->size() );
            for (size_t n=0; n<Rx[nr]->size(); ++n) {
                (*traveltimes[nr])[n] = this->getTraveltime((*Rx[nr])[n], threadNo);
            }
        }
    }

    template<typename T1, typename T2, typename NODE, typename S>
    void Grid2Ducsp<T1,T2,NODE,S>::raytrace(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            const std::vector<S>& Rx,
                                            std::vector<T1>& traveltimes,
                                            std::vector<std::vector<S>>& r_data,
                                            const size_t threadNo) const {

        this->checkPts(Tx);
        this->checkPts(Rx);

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< NODE*, std::vector<NODE*>,
        CompareNodePtr<T1>> queue( cmp );

        std::vector<NODE> txNodes;
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

            traveltimes[n] = this->getTraveltime(Rx[n], nodeParentRx,
                                                 cellParentRx, threadNo);

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
                    r_data[n].push_back(S(txNodes[ns]));
                    r_data[n].push_back(Rx[n]);
                    flag = true;
                    break;
                }
            }
            if ( flag ) continue;

            // Rx are in nodes (not txNodes)
            std::vector<NODE> *node_p;
            node_p = &(this->nodes);

            std::vector<S> r_tmp;
            T2 iChild, iParent = nodeParentRx;
            S child;

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
    }

    template<typename T1, typename T2, typename NODE, typename S>
    void Grid2Ducsp<T1,T2,NODE,S>::raytrace(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            const std::vector<const std::vector<S>*>& Rx,
                                            std::vector<std::vector<T1>*>& traveltimes,
                                            std::vector<std::vector<std::vector<S>>*>& r_data,
                                            const size_t threadNo) const {

        this->checkPts(Tx);
        for ( size_t n=0; n<Rx.size(); ++n ) {
            this->checkPts(*Rx[n]);
        }

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< NODE*, std::vector<NODE*>,
        CompareNodePtr<T1>> queue( cmp );

        std::vector<NODE> txNodes;
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
                for ( size_t ns=0; ns<txNodes.size(); ++ns ) {
                    if ( nodeParentRx == txNodes[ns].getGridIndex() ) {
                        // insert Tx at begining
                        (*r_data[nr])[n].push_back(S(txNodes[ns]));
                        (*r_data[nr])[n].push_back((*Rx[nr])[n]);
                        flag = true;
                        break;
                    }
                }
                if ( flag ) continue;

                // Rx are in nodes (not txNodes)
                std::vector<NODE> *node_p;
                node_p = &(this->nodes);

                std::vector<S> r_tmp;
                T2 iChild, iParent = nodeParentRx;
                S child;

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

    template<typename T1, typename T2, typename NODE, typename S>
    void Grid2Ducsp<T1,T2,NODE,S>::raytrace(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            const std::vector<S>& Rx,
                                            std::vector<T1>& traveltimes,
                                            std::vector<std::vector<S>>& r_data,
                                            std::vector<std::vector<siv<T1>>>& l_data,
                                            const size_t threadNo) const {

        this->checkPts(Tx);
        this->checkPts(Rx);

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< NODE*, std::vector<NODE*>,
        CompareNodePtr<T1>> queue( cmp );

        std::vector<NODE> txNodes;
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
            traveltimes[n] = this->getTraveltime(Rx[n], nodeParentRx, cellParentRx,
                                                 threadNo);

            bool flag=false;
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
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
                    r_data[n].push_back(S(txNodes[ns]));
                    r_data[n].push_back(Rx[n]);
                    cell.i = cellParentRx;
                    cell.v = Rx[n].getDistance(txNodes[ns]);
                    l_data[n].push_back( cell );
                    flag = true;
                    flag = true;
                    break;
                }
            }
            if ( flag ) continue;

            // Rx are in nodes (not txNodes)
            std::vector<NODE> *node_p;
            node_p = &(this->nodes);

            std::vector<S> r_tmp;
            T2 iChild, iParent = nodeParentRx;
            S child;

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
            child = (*node_p)[iParent];
            r_tmp.push_back( child );

            //  must be sorted to build matrix L
            sort(l_data[n].begin(), l_data[n].end(), CompareSiv_i<T1>());

            // the order should be from Tx to Rx, so we reorder...
            iParent = static_cast<T2>(r_tmp.size());
            r_data[n].resize( r_tmp.size() );
            for ( size_t nn=0; nn<r_data[n].size(); ++nn ) {
                r_data[n][nn] = r_tmp[ iParent-1-nn ];
            }
        }
    }


    template<typename T1, typename T2, typename NODE, typename S>
    void Grid2Ducsp<T1,T2,NODE,S>::initQueue(const std::vector<S>& Tx,
                                             const std::vector<T1>& t0,
                                             std::priority_queue<NODE*,
                                             std::vector<NODE*>,
                                             CompareNodePtr<T1>>& queue,
                                             std::vector<NODE>& txNodes,
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
                txNodes.push_back( NODE(t0[n], Tx[n], this->nThreads, threadNo) );
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


    template<typename T1, typename T2, typename NODE, typename S>
    void Grid2Ducsp<T1,T2,NODE,S>::propagate(std::priority_queue<NODE*,
                                             std::vector<NODE*>,
                                             CompareNodePtr<T1>>& queue,
                                             std::vector<bool>& inQueue,
                                             std::vector<bool>& frozen,
                                             const size_t threadNo) const {

#ifdef DEBUG_OF
        std::string fname;
        size_t n=1;
#endif

        while ( !queue.empty() ) {

            const NODE* source = queue.top();
            queue.pop();
            inQueue[ source->getGridIndex() ] = false;
            frozen[ source->getGridIndex() ] = true;

#ifdef DEBUG_OF
            char fname[80];
            sprintf(fname, "src_node%04zd", n);
            ofstream of(fname);
            of << source->getGridIndex() << '\n';
            of.close();

            sprintf(fname, "in_queue%04zd", n);
            of.open(fname);
            for ( size_t ni=0; ni<inQueue.size(); ++ni ) {
                if ( inQueue[ni]==true ) {
                    of << ni << '\n';
                }
            }
            of.close();

            sprintf(fname, "timed%04zd", n);
            of.open(fname);
            for ( size_t ni=0; ni<this->nodes.size(); ++ni ) {
                if ( this->nodes[ni].getTT(threadNo)<std::numeric_limits<T1>::max() && inQueue[ni]==false && source->getGridIndex()!=ni ) {
                    of << ni << '\n';
                }
            }
            of.close();

            sprintf(fname, "neighbors%04zd", n);
            of.open(fname);
#endif
            for ( size_t no=0; no<source->getOwners().size(); ++no ) {

                T2 cellNo = source->getOwners()[no];

                for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                    T2 neibNo = this->neighbors[cellNo][k];
                    if ( neibNo == source->getGridIndex() || frozen[neibNo] ) {
                        continue;
                    }
#ifdef DEBUG_OF
                    of << neibNo << '\n';
#endif
                    // compute dt
                    T1 dt = this->computeDt(*source, this->nodes[neibNo], cellNo);

                    if (source->getTT(threadNo)+dt < this->nodes[neibNo].getTT(threadNo)) {
                        this->nodes[neibNo].setTT( source->getTT(threadNo)+dt, threadNo );
                        this->nodes[neibNo].setnodeParent(source->getGridIndex(),threadNo);
                        this->nodes[neibNo].setCellParent(cellNo, threadNo );

                        if ( !inQueue[neibNo] ) {
                            queue.push( &(this->nodes[neibNo]) );
                            inQueue[neibNo] = true;
                        }
                    }
                }
            }
#ifdef DEBUG_OF
            of.close();
            n++;
#endif
        }
    }

}

#endif

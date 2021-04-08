//
//  Grid2Dunfm.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-04-15.
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

#ifndef ttcr_Grid2Dunfm_h
#define ttcr_Grid2Dunfm_h

#include <fstream>
#include <queue>

#include "Grid2Dun.h"

namespace ttcr {

    template<typename T1, typename T2, typename NODE, typename S>
    class Grid2Dunfm : public Grid2Dun<T1,T2,NODE,S> {
    public:
        Grid2Dunfm(const std::vector<S>& no,
                   const std::vector<triangleElem<T2>>& tri,
                   const bool ttrp, const size_t nt=1,
                   const bool procObtuse=true) :
        Grid2Dun<T1, T2,NODE,S>(no, tri, ttrp, nt)
        {
            this->buildGridNodes(no, nt);
            this->template buildGridNeighbors<NODE>(this->nodes);
            if ( procObtuse ) this->processObtuse();
        }

        ~Grid2Dunfm() {
        }

    private:

        void initBand(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      std::priority_queue<NODE*,
                      std::vector<NODE*>,
                      CompareNodePtr<T1>>&,
                      std::vector<NODE>&,
                      std::vector<bool>&,
                      std::vector<bool>&,
                      const size_t) const;

        void propagate(std::priority_queue<NODE*,
                       std::vector<NODE*>,
                       CompareNodePtr<T1>>&,
                       std::vector<bool>&,
                       std::vector<bool>&,
                       const size_t) const;

        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<S>& Rx,
                      const size_t threadNo) const;

        void raytrace(const std::vector<S>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<const std::vector<S>*>& Rx,
                      const size_t threadNo) const;
    };

    template<typename T1, typename T2, typename NODE, typename S>
    void Grid2Dunfm<T1,T2,NODE,S>::raytrace(const std::vector<S>& Tx,
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
        CompareNodePtr<T1>> narrow_band( cmp );

        std::vector<NODE> txNodes;
        std::vector<bool> inQueue( this->nodes.size(), false );
        std::vector<bool> frozen( this->nodes.size(), false );

        initBand(Tx, t0, narrow_band, txNodes, inQueue, frozen, threadNo);

        propagate(narrow_band, inQueue, frozen, threadNo);
    }

    template<typename T1, typename T2, typename NODE, typename S>
    void Grid2Dunfm<T1,T2,NODE,S>::raytrace(const std::vector<S>& Tx,
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
        std::priority_queue< NODE*, std::vector<NODE*>,
        CompareNodePtr<T1>> narrow_band( cmp );

        std::vector<NODE> txNodes;
        std::vector<bool> inBand( this->nodes.size(), false );
        std::vector<bool> frozen( this->nodes.size(), false );

        initBand(Tx, t0, narrow_band, txNodes, inBand, frozen, threadNo);

        propagate(narrow_band, inBand, frozen, threadNo);
    }

    template<typename T1, typename T2, typename NODE, typename S>
    void Grid2Dunfm<T1,T2,NODE,S>::initBand(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            std::priority_queue<NODE*,
                                            std::vector<NODE*>,
                                            CompareNodePtr<T1>>& narrow_band,
                                            std::vector<NODE>& txNodes,
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
                    }

                    break;
                }
            }
            if ( found==false ) {

                T2 cellNo = this->getCellNo(Tx[n]);
                for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                    T2 neibNo = this->neighbors[cellNo][k];

                    // compute dt
                    T1 dt = this->nodes[neibNo].getDistance(Tx[n])*this->nodes[neibNo].getNodeSlowness();

                    this->nodes[neibNo].setTT( t0[n]+dt, threadNo );
                    narrow_band.push( &(this->nodes[neibNo]) );
                    inBand[neibNo] = true;
                    frozen[neibNo] = true;

                }
            }
        }
    }

    template<typename T1, typename T2, typename NODE, typename S>
    void Grid2Dunfm<T1,T2,NODE,S>::propagate(std::priority_queue<NODE*,
                                             std::vector<NODE*>,
                                             CompareNodePtr<T1>>& narrow_band,
                                             std::vector<bool>& inNarrowBand,
                                             std::vector<bool>& frozen,
                                             const size_t threadNo) const {

        //    size_t n=1;
        while ( !narrow_band.empty() ) {

            const NODE* source = narrow_band.top();
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

                    this->localSolver( &(this->nodes[neibNo]), threadNo );

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

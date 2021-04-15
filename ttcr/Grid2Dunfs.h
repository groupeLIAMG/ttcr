//
//  Grid2Dunfs.h
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

#ifndef ttcr_Grid2Dunfs_h
#define ttcr_Grid2Dunfs_h

#include <cmath>
#include <fstream>
#include <queue>

#include "Grid2Dun.h"
#include "Metric.h"

namespace ttcr {

    template<typename T1, typename T2, typename NODE, typename S>
    class Grid2Dunfs : public Grid2Dun<T1,T2,NODE,S> {
    public:
        Grid2Dunfs(const std::vector<S>& no,
                   const std::vector<triangleElem<T2>>& tri,
                   const T1 eps, const int maxit, const bool ttrp,
                   const size_t nt=1,
                   const bool procObtuse=true) :
        Grid2Dun<T1,T2,NODE,S>(no, tri, ttrp, nt),
        epsilon(eps), nitermax(maxit), niter_final(0), sorted()
        {
            this->buildGridNodes(no, nt);
            this->template buildGridNeighbors<NODE>(this->nodes);
            if ( procObtuse ) this->processObtuse();
        }

        Grid2Dunfs(const std::vector<S>& no,
                   const std::vector<triangleElem<T2>>& tri,
                   const T1 eps, const int maxit,
                   const std::vector<S>& refPts, const int order,
                   const bool ttrp, const size_t nt=1,
                   const bool procObtuse=true) :
        Grid2Dun<T1,T2,NODE,S>(no, tri, ttrp, nt),
        epsilon(eps), nitermax(maxit), niter_final(0), sorted()
        {
            this->buildGridNodes(no, nt);
            this->template buildGridNeighbors<NODE>(this->nodes);
            if ( procObtuse ) this->processObtuse();
            initOrdering(refPts, order);
        }

        ~Grid2Dunfs() {
        }

        void initOrdering(const std::vector<S>& refPts, const int order);

        const int get_niter() const { return niter_final; }

    private:
        T1 epsilon;
        int nitermax;
        mutable int niter_final;
        std::vector<std::vector<NODE*>> sorted;

        void initTx(const std::vector<S>& Tx, const std::vector<T1>& t0,
                    std::vector<bool>& frozen, const size_t threadNo) const;

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
    void Grid2Dunfs<T1,T2,NODE,S>::initOrdering(const std::vector<S>& refPts,
                                                const int order) {
        sorted.resize( refPts.size() );

        Metric<T1> *m;
        if ( order == 1 )
            m = new Metric1<T1>();
        else
            m = new Metric2<T1>();

        std::priority_queue<siv<T1>,std::vector<siv<T1>>,CompareSiv_vr<T1>> queue;

        for ( size_t np=0; np<refPts.size(); ++np ) {

            for ( size_t n=0; n<this->nodes.size(); ++n ) {
                queue.push( {n, m->l(this->nodes[n], refPts[np])} );
            }

            while ( !queue.empty() ) {
                siv<T1> s = queue.top();
                queue.pop();
                sorted[np].push_back( &(this->nodes[s.i]) );
            }
        }

        delete m;
    }


    template<typename T1, typename T2, typename NODE, typename S>
    void Grid2Dunfs<T1,T2,NODE,S>::raytrace(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            const std::vector<S>& Rx,
                                            const size_t threadNo) const {

        this->checkPts(Tx);
        this->checkPts(Rx);

        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        std::vector<bool> frozen( this->nodes.size(), false );
        initTx(Tx, t0, frozen, threadNo);

        std::vector<T1> times( this->nodes.size() );
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            times[n] = this->nodes[n].getTT( threadNo );
        }

        int niter=0;
        T1 change = std::numeric_limits<T1>::max();
        while ( change >= epsilon && niter<nitermax ) {

            for ( size_t i=0; i<sorted.size(); ++i ) {

                // ascending
                for ( auto vertexC=sorted[i].begin(); vertexC!=sorted[i].end(); ++vertexC ) {
                    if ( !frozen[(*vertexC)->getGridIndex()] )
                        this->localSolver(*vertexC, threadNo);
                }

                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = std::abs( times[n] - this->nodes[n].getTT(threadNo) );

                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                if ( change < epsilon ) {
                    niter++;
                    break;
                }

                // descending
                for ( auto vertexC=sorted[i].rbegin(); vertexC!=sorted[i].rend(); ++vertexC ) {
                    if ( !frozen[(*vertexC)->getGridIndex()] )
                        this->localSolver(*vertexC, threadNo);
                }

                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = std::abs( times[n] - this->nodes[n].getTT(threadNo) );

                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                if ( change < epsilon ) {
                    niter++;
                    break;
                }
            }
            niter++;
        }
        niter_final = niter;
    }


    template<typename T1, typename T2, typename NODE, typename S>
    void Grid2Dunfs<T1,T2,NODE,S>::raytrace(const std::vector<S>& Tx,
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

        std::vector<bool> frozen( this->nodes.size(), false );
        initTx(Tx, t0, frozen, threadNo);

        std::vector<T1> times( this->nodes.size() );
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            times[n] = this->nodes[n].getTT( threadNo );
        }

        int niter=0;
        T1 change = std::numeric_limits<T1>::max();
        while ( change >= epsilon && niter<nitermax ) {

            for ( size_t i=0; i<sorted.size(); ++i ) {

                // ascending
                for ( auto vertexC=sorted[i].begin(); vertexC!=sorted[i].end(); ++vertexC ) {
                    if ( !frozen[(*vertexC)->getGridIndex()] )
                        this->localSolver(*vertexC, threadNo);
                }

                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = std::abs( times[n] - this->nodes[n].getTT(threadNo) );

                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                if ( change < epsilon ) {
                    niter++;
                    break;
                }

                // descending
                for ( auto vertexC=sorted[i].rbegin(); vertexC!=sorted[i].rend(); ++vertexC ) {
                    if ( !frozen[(*vertexC)->getGridIndex()] )
                        this->localSolver(*vertexC, threadNo);
                }

                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = std::abs( times[n] - this->nodes[n].getTT(threadNo) );

                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                if ( change < epsilon ) {
                    niter++;
                    break;
                }

            }
            niter++;
        }
        niter_final = niter;
    }

    template<typename T1, typename T2, typename NODE, typename S>
    void Grid2Dunfs<T1,T2,NODE,S>::initTx(const std::vector<S>& Tx,
                                          const std::vector<T1>& t0,
                                          std::vector<bool>& frozen,
                                          const size_t threadNo) const {

        for (size_t n=0; n<Tx.size(); ++n) {
            bool found = false;
            for ( size_t nn=0; nn<this->nodes.size(); ++nn ) {
                if ( this->nodes[nn] == Tx[n] ) {
                    found = true;
                    this->nodes[nn].setTT( t0[n], threadNo );
                    frozen[nn] = true;

                    // populate around Tx
                    for ( size_t no=0; no<this->nodes[nn].getOwners().size(); ++no ) {

                        T2 cellNo = this->nodes[nn].getOwners()[no];
                        for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                            T2 neibNo = this->neighbors[cellNo][k];
                            if ( neibNo == nn ) continue;
                            T1 dt = this->computeDt(this->nodes[nn], this->nodes[neibNo]);

                            if ( t0[n]+dt < this->nodes[neibNo].getTT(threadNo) ) {
                                this->nodes[neibNo].setTT( t0[n]+dt, threadNo );
                                // this->nodes[neibNo].setnodeParent(this->nodes[nn].getGridIndex(),threadNo);
                                // this->nodes[neibNo].setCellParent(cellNo, threadNo );
                                // frozen[neibNo] = true;
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
                    // first, get slowness at point
                    T1 slowness_pt = this->computeSlowness(Tx[n], cellNo);
                    T1 dt = this->computeDt(this->nodes[neibNo], Tx[n], slowness_pt);

                    this->nodes[neibNo].setTT( t0[n]+dt, threadNo );
                    frozen[neibNo] = true;

                }
            }
        }
    }
}

#endif

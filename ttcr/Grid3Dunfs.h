//
//  Grid3Dunfs.h
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

#ifndef ttcr_Grid3Dunfs_h
#define ttcr_Grid3Dunfs_h

#include <cmath>
#include <fstream>
#include <queue>
#include <vector>

#include "Grid3Dun.h"
#include "Node3Dn.h"
#include "Metric.h"

namespace ttcr {

    template<typename T1, typename T2>
    class Grid3Dunfs : public Grid3Dun<T1,T2,Node3Dn<T1,T2>> {
    public:
        Grid3Dunfs(const std::vector<sxyz<T1>>& no,
                   const std::vector<tetrahedronElem<T2>>& tet,
                   const T1 eps,
                   const int maxit,
                   const int rp,
                   const bool iv,
                   const bool rptt,
                   const T1 md,
                   const size_t nt=1,
                   const bool _translateOrigin=false) :
        Grid3Dun<T1,T2,Node3Dn<T1,T2>>(no, tet, rp, iv, rptt, md, nt, _translateOrigin),
        epsilon(eps), nitermax(maxit), S(), niter_final(0)
        {
            this->buildGridNodes(no, nt);
            this->template buildGridNeighbors<Node3Dn<T1,T2>>(this->nodes);
        }
        Grid3Dunfs(const std::vector<sxyz<T1>>& no,
                   const std::vector<tetrahedronElem<T2>>& tet,
                   const T1 eps,
                   const int maxit,
                   const std::vector<sxyz<T1>>& refPts,
                   const int order,
                   const int rp,
                   const bool iv,
                   const bool rptt,
                   const T1 md,
                   const size_t nt=1,
                   const bool _translateOrigin=false) :
        Grid3Dun<T1,T2,Node3Dn<T1,T2>>(no, tet, rp, iv, rptt, md, nt, _translateOrigin),
        epsilon(eps), nitermax(maxit), S(), niter_final(0)
        {
            this->buildGridNodes(no, nt);
            this->buildGridNeighbors(this->nodes);
            initOrdering(refPts, order);
        }

        ~Grid3Dunfs() {
        }

        void initOrdering(const std::vector<sxyz<T1>>& refPts, const int order);

        const int get_niter() const { return niter_final; }

    private:
        T1 epsilon;
        int nitermax;
        std::vector<std::vector<Node3Dn<T1,T2>*>> S;
        mutable int niter_final;

        void initTx(const std::vector<sxyz<T1>>& Tx, const std::vector<T1>& t0,
                    std::vector<bool>& frozen, const size_t threadNo) const;

        void initBand(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      std::priority_queue<Node3Dn<T1,T2>*,
                      std::vector<Node3Dn<T1,T2>*>,
                      CompareNodePtr<T1>>&,
                      std::vector<Node3Dn<T1,T2>>&,
                      std::vector<bool>&,
                      std::vector<bool>&,
                      const size_t) const;

        void propagate(std::priority_queue<Node3Dn<T1,T2>*,
                       std::vector<Node3Dn<T1,T2>*>,
                       CompareNodePtr<T1>>&,
                       std::vector<bool>&,
                       std::vector<bool>&,
                       const size_t) const;

        void raytrace(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<sxyz<T1>>& Rx,
                      const size_t threadNo=0) const;

        void raytrace(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      const size_t threadNo=0) const;

    };

    template<typename T1, typename T2>
    void Grid3Dunfs<T1,T2>::initOrdering(const std::vector<sxyz<T1>>& refPts,
                                         const int order) {
        S.resize( refPts.size() );

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
                S[np].push_back( &(this->nodes[s.i]) );
            }
        }

        delete m;
    }

    template<typename T1, typename T2>
    void Grid3Dunfs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<sxyz<T1>>& Rx,
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

        int niter = 0;
        T1 change = std::numeric_limits<T1>::max();
        while ( change >= epsilon && niter<nitermax ) {

            for ( size_t i=0; i<S.size(); ++i ) {

                // ascending
                for ( auto vertexC=S[i].begin(); vertexC!=S[i].end(); ++vertexC ) {
                    if ( !frozen[(*vertexC)->getGridIndex()] )
                        //                    this->local3Dsolver(*vertexC, threadNo);
                        this->localUpdate3D(*vertexC, threadNo);
                }

                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = std::abs( times[n] - this->nodes[n].getTT(threadNo) );

                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                if ( change < epsilon ) {
                    break;
                }

                // descending
                for ( auto vertexC=S[i].rbegin(); vertexC!=S[i].rend(); ++vertexC ) {
                    if ( !frozen[(*vertexC)->getGridIndex()] )
                        //                    this->local3Dsolver(*vertexC, threadNo);
                        this->localUpdate3D(*vertexC, threadNo);
                }

                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = std::abs( times[n] - this->nodes[n].getTT(threadNo) );

                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                if ( change < epsilon ) {
                    break;
                }
            }
            niter++;
        }
        niter_final = niter;
    }

    template<typename T1, typename T2>
    void Grid3Dunfs<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                     const std::vector<T1>& t0,
                                     const std::vector<std::vector<sxyz<T1>>>& Rx,
                                     const size_t threadNo) const {

        this->checkPts(Tx);
        for ( size_t n=0; n<Rx.size(); ++n ) {
            this->checkPts(Rx[n]);
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

        int niter = 0;
        T1 change = std::numeric_limits<T1>::max();
        while ( change >= epsilon && niter<nitermax ) {

            for ( size_t i=0; i<S.size(); ++i ) {

                // ascending
                for ( auto vertexC=S[i].begin(); vertexC!=S[i].end(); ++vertexC ) {
                    if ( !frozen[(*vertexC)->getGridIndex()] )
                        //                    this->local3Dsolver(*vertexC, threadNo);
                        this->localUpdate3D(*vertexC, threadNo);
                }

                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = std::abs( times[n] - this->nodes[n].getTT(threadNo) );

                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                if ( change < epsilon ) {
                    break;
                }

                // descending
                for ( auto vertexC=S[i].rbegin(); vertexC!=S[i].rend(); ++vertexC ) {
                    if ( !frozen[(*vertexC)->getGridIndex()] )
                        //                    this->local3Dsolver(*vertexC, threadNo);
                        this->localUpdate3D(*vertexC, threadNo);
                }

                change = 0.0;
                for ( size_t n=0; n<this->nodes.size(); ++n ) {
                    T1 dt = std::abs( times[n] - this->nodes[n].getTT(threadNo) );

                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                if ( change < epsilon ) {
                    break;
                }

            }
            niter++;
        }
        niter_final = niter;
    }


    template<typename T1, typename T2>
    void Grid3Dunfs<T1,T2>::initTx(const std::vector<sxyz<T1>>& Tx,
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
                                    //frozen[neibNo] = true;
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
                                    if ( this->nodes[no].getTT(threadNo) == std::numeric_limits<T1>::max() ) nodes_added++;
                                    this->nodes[no].setTT( t0[n]+dt, threadNo );
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

                    break;
                }
            }
            if ( found==false ) {

                T2 sTx = this->computeSlowness( Tx[n] );

                T2 cellNo = this->getCellNo(Tx[n]);
                if ( Grid3Dun<T1,T2,Node3Dn<T1,T2>>::source_radius == 0.0 ) {
                    for ( size_t k=0; k< this->neighbors[cellNo].size(); ++k ) {
                        T2 neibNo = this->neighbors[cellNo][k];
                        // compute dt
                        T1 dt = this->computeDt(this->nodes[neibNo], Tx[n], sTx);

                        this->nodes[neibNo].setTT( t0[n]+dt, threadNo );
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
                                if ( this->nodes[no].getTT(threadNo) == std::numeric_limits<T1>::max() ) nodes_added++;
                                this->nodes[no].setTT( t0[n]+dt, threadNo );
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

}

#endif

//
//  Grid3Ducdsp.h
//  ttcr
//
//  Created by Bernard Giroux on 18-10-11.
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

#ifndef ttcr_Grid3Ducdsp_h
#define ttcr_Grid3Ducdsp_h

#include "Grid3Duc.h"
#include "Node3Dc.h"
#include "Node3Dcd.h"

namespace ttcr {

    template<typename T1, typename T2>
    class Grid3Ducdsp : public Grid3Duc<T1,T2,Node3Dc<T1,T2>> {
    public:
        Grid3Ducdsp(const std::vector<sxyz<T1>>& no,
                    const std::vector<tetrahedronElem<T2>>& tet,
                    const int ns, const int nd, const T1 rad,
                    const int rp, const bool rptt, const T1 min_dist,
                    const T1 drad, const bool useEdgeLength=true,
                    const size_t nt=1, const bool _translateOrigin=false) :
        Grid3Duc<T1,T2,Node3Dc<T1,T2>>(no, tet, rp, rptt, min_dist, nt, _translateOrigin),
        nSecondary(ns), nTertiary(nd), nPermanent(0),
        dyn_radius(drad),
        tempNodes(std::vector<std::vector<Node3Dcd<T1,T2>>>(nt)),
        tempNeighbors(std::vector<std::vector<std::vector<T2>>>(nt))
        {
            this->buildGridNodes(no, ns, nt);
            this->template buildGridNeighbors<Node3Dc<T1,T2>>(this->nodes);
            this->source_radius = rad;
            nPermanent = static_cast<T2>(this->nodes.size());
            for ( size_t n=0; n<nt; ++n ) {
                tempNeighbors[n].resize(tet.size());
            }
            if (useEdgeLength) dyn_radius *= this->getAverageEdgeLength();
        }

        ~Grid3Ducdsp() {
        }

    private:
        T2 nSecondary;
        T2 nTertiary;
        T2 nPermanent;
        T1 dyn_radius;

        // we will store temporary nodes in a separate container.  This is to
        // allow threaded computations with different Tx (location of temp
        // nodes vary from one Tx to the other)
        mutable std::vector<std::vector<Node3Dcd<T1,T2>>> tempNodes;
        mutable std::vector<std::vector<std::vector<T2>>> tempNeighbors;

        void addTemporaryNodes(const std::vector<sxyz<T1>>&, const size_t) const;

        void initQueue(const std::vector<sxyz<T1>>& Tx,
                       const std::vector<T1>& t0,
                       std::priority_queue<Node3Dc<T1,T2>*,
                       std::vector<Node3Dc<T1,T2>*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<Node3Dcd<T1,T2>>& txNodes,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       const size_t threadNo) const;

        void propagate(std::priority_queue<Node3Dc<T1,T2>*,
                       std::vector<Node3Dc<T1,T2>*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       const size_t threadNo) const;

        void raytrace(const std::vector<sxyz<T1>>&,
                      const std::vector<T1>&,
                      const std::vector<sxyz<T1>>&,
                      const size_t=0) const;

        void raytrace(const std::vector<sxyz<T1>>&,
                      const std::vector<T1>&,
                      const std::vector<std::vector<sxyz<T1>>>&,
                      const size_t=0) const;

    };

    template<typename T1, typename T2>
    void Grid3Ducdsp<T1,T2>::addTemporaryNodes(const std::vector<sxyz<T1>>& Tx,
                                               const size_t threadNo) const {

        // clear previously assigned nodes
        tempNodes[threadNo].clear();
        for ( size_t nt=0; nt<this->tetrahedra.size(); ++nt ) {
            tempNeighbors[threadNo][nt].clear();
        }

        // find cells surrounding Tx
        std::set<T2> txCells;
        for ( size_t nt=0; nt<this->tetrahedra.size(); ++nt ) {
            // centroid of tet
            sxyz<T1> cent = sxyz<T1>(this->nodes[this->tetrahedra[nt].i[0]] +
                                     this->nodes[this->tetrahedra[nt].i[1]]);
            cent += this->nodes[this->tetrahedra[nt].i[2]] +
            this->nodes[this->tetrahedra[nt].i[3]];
            cent *= 0.25;
            for (size_t n=0; n<Tx.size(); ++n) {
                if ( cent.getDistance(Tx[n]) <= this->dyn_radius ) {
                    txCells.insert(static_cast<T2>(nt));
                }
            }
        }
        if ( verbose )
            std::cout << "\n  *** thread no " << threadNo << ": found " << txCells.size() << " cells within radius ***" << std::endl;

        std::set<T2> adjacentCells(txCells.begin(), txCells.end());

        T2 iNodes[4][3] = {
            {0,1,2},  // (relative) indices of nodes of 1st triangle
            {1,2,3},  // (relative) indices of nodes of 2nd triangle
            {0,2,3},  // (relative) indices of nodes of 3rd triangle
            {0,1,3}   // (relative) indices of nodes of 4th triangle
        };
        std::map<std::array<T2,2>,std::vector<T2>> lineMap;
        std::array<T2,2> lineKey;
        typename std::map<std::array<T2,2>,std::vector<T2>>::iterator lineIt;

        std::map<std::array<T2,3>,std::vector<T2>> faceMap;
        std::array<T2,3> faceKey;
        typename std::map<std::array<T2,3>,std::vector<T2>>::iterator faceIt;

        // edge nodes
        T2 nTmpNodes = 0;
        Node3Dcd<T1,T2> tmpNode;
        size_t nDynTot = (nSecondary+1) * nTertiary;  // total number of dynamic nodes on edges

        for ( auto cell=txCells.begin(); cell!=txCells.end(); cell++ ) {

            //  adjacent cells to the tetrahedron where new nodes will be added
            for ( size_t i=0; i<4; ++i ) {
                T2 vertex = this->neighbors[*cell][i];
                for ( size_t c=0; c<this->nodes[vertex].getOwners().size(); ++c ) {
                    adjacentCells.insert(this->nodes[vertex].getOwners()[c]);
                }
            }

            // for each triangle
            for ( T2 ntri=0; ntri<4; ++ntri ) {

                // start from ntri to avoid redundancy
                for ( size_t nl=ntri; nl<3; ++nl ) {

                    lineKey = {this->tetrahedra[*cell].i[ iNodes[ntri][nl] ],
                        this->tetrahedra[*cell].i[ iNodes[ntri][(nl+1)%3] ]};
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

                    sxyz<T1> d = (this->nodes[lineKey[1]]-this->nodes[lineKey[0]])/static_cast<T1>(nDynTot+nSecondary+1);

                    size_t nd = 0;
                    for ( size_t n2=0; n2<nSecondary+1; ++n2 ) {
                        for ( size_t n3=0; n3<nTertiary; ++n3 ) {
                            tmpNode.setXYZindex(this->nodes[lineKey[0]].getX()+(1+n2*(nTertiary+1)+n3)*d.x,
                                                this->nodes[lineKey[0]].getY()+(1+n2*(nTertiary+1)+n3)*d.y,
                                                this->nodes[lineKey[0]].getZ()+(1+n2*(nTertiary+1)+n3)*d.z,
                                                nPermanent+nTmpNodes );

                            lineMap[lineKey][nd++] = nTmpNodes++;
                            tempNodes[threadNo].push_back( tmpNode );
                            tempNodes[threadNo].back().pushOwner( *cell );
                        }
                    }
                }
            }
        }
        // on faces
        size_t ncut = nDynTot + nSecondary - 1;
        size_t nSecNodes = 0;
        for ( int n=1; n<=(nSecondary-1); ++n ) {
            nSecNodes += n;
        }
        size_t nFaceNodes = 0;
        for ( int n=1; n<=ncut; ++n ) {
            nFaceNodes += n;
        }
        nFaceNodes -= nSecNodes;

        for ( auto cell=txCells.begin(); cell!=txCells.end(); cell++ ) {

            // for each triangle
            for ( T2 ntri=0; ntri<4; ++ntri ) {

                faceKey = {this->tetrahedra[*cell].i[ iNodes[ntri][0] ],
                    this->tetrahedra[*cell].i[ iNodes[ntri][1] ],
                    this->tetrahedra[*cell].i[ iNodes[ntri][2] ]};
                std::sort(faceKey.begin(), faceKey.end());

                faceIt = faceMap.find( faceKey );
                if ( faceIt == faceMap.end() ) {
                    // not found, insert new pair
                    faceMap[ faceKey ] = std::vector<T2>(nFaceNodes);
                } else {
                    for ( size_t n=0; n<faceIt->second.size(); ++n ) {
                        // setting owners
                        tempNodes[threadNo][ faceIt->second[n] ].pushOwner( *cell );
                    }
                    continue;
                }

                sxyz<T1> d1 = (this->nodes[faceKey[1]]-this->nodes[faceKey[0]])/static_cast<T1>(nDynTot+nSecondary+1);
                sxyz<T1> d2 = (this->nodes[faceKey[1]]-this->nodes[faceKey[2]])/static_cast<T1>(nDynTot+nSecondary+1);

                size_t ifn = 0;
                size_t n = 0;
                for ( int n2=nSecondary; n2>-1; --n2 ) {
                    for ( size_t n3=0; n3<nTertiary; ++n3 ) {
                        sxyz<T1> pt1 = this->nodes[faceKey[0]]+static_cast<T1>(1+n)*d1;
                        sxyz<T1> pt2 = this->nodes[faceKey[2]]+static_cast<T1>(1+n)*d2;

                        size_t nseg = ncut + 1 - n;
                        sxyz<T1> d = (pt2-pt1)/static_cast<T1>(nseg);
                        for ( size_t n4=0; n4<nseg-1; ++n4 ) {
                            tmpNode.setXYZindex(pt1.x+(1+n4)*d.x,
                                                pt1.y+(1+n4)*d.y,
                                                pt1.z+(1+n4)*d.z,
                                                nPermanent+nTmpNodes );

                            faceMap[faceKey][ifn++] = nTmpNodes++;
                            tempNodes[threadNo].push_back( tmpNode );
                            tempNodes[threadNo].back().pushOwner( *cell );
                        }
                        n++;
                    }
                    sxyz<T1> pt1 = this->nodes[faceKey[0]]+static_cast<T1>(1+n)*d1;
                    sxyz<T1> pt2 = this->nodes[faceKey[2]]+static_cast<T1>(1+n)*d2;

                    size_t nseg = ncut + 1 - n;
                    if ( nseg == 0 ) break;

                    sxyz<T1> d = (pt2-pt1)/static_cast<T1>(nseg);
                    size_t n5 = 0;
                    for ( size_t n4=0; n4<n2; ++n4 ) {
                        for ( size_t n3=0; n3<nTertiary; ++n3 ) {
                            tmpNode.setXYZindex(pt1.x+(1+n5)*d.x,
                                                pt1.y+(1+n5)*d.y,
                                                pt1.z+(1+n5)*d.z,
                                                nPermanent+nTmpNodes );
                            n5++;

                            faceMap[faceKey][ifn++] = nTmpNodes++;
                            tempNodes[threadNo].push_back( tmpNode );
                            tempNodes[threadNo].back().pushOwner( *cell );
                        }
                        n5++;
                    }
                    n++;
                }
            }
        }

        for ( auto cell=txCells.begin(); cell!=txCells.end(); ++cell ) {
            adjacentCells.erase(*cell);
        }
        for ( auto adj=adjacentCells.begin(); adj!=adjacentCells.end(); ++adj ) {
            for ( T2 ntri=0; ntri<4; ++ntri ) {
                for ( size_t nl=ntri; nl<3; ++nl ) {

                    lineKey = {this->tetrahedra[*adj].i[ iNodes[ntri][nl] ],
                        this->tetrahedra[*adj].i[ iNodes[ntri][(nl+1)%3] ]};
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

                faceKey = {this->tetrahedra[*adj].i[ iNodes[ntri][0] ],
                    this->tetrahedra[*adj].i[ iNodes[ntri][1] ],
                    this->tetrahedra[*adj].i[ iNodes[ntri][2] ]};
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

    template<typename T1, typename T2>
    void Grid3Ducdsp<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                      const std::vector<T1>& t0,
                                      const std::vector<sxyz<T1>>& Rx,
                                      const size_t threadNo) const {
        this->checkPts(Tx);
        this->checkPts(Rx);

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

    template<typename T1, typename T2>
    void Grid3Ducdsp<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
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

    template<typename T1, typename T2>
    void Grid3Ducdsp<T1,T2>::initQueue(const std::vector<sxyz<T1>>& Tx,
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

                //                prepropagate(txNodes.back(), queue, inQueue, frozen, threadNo); // See description in the function declaration

                queue.push( &(txNodes.back()) );    //Don't use if prepropagate is used
                inQueue.push_back( true );            //Don't use if prepropagate is used

            }
        }
    }

    template<typename T1, typename T2>
    void Grid3Ducdsp<T1,T2>::propagate(std::priority_queue<Node3Dc<T1,T2>*,
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
                    T1 dt = this->computeDt(*src, this->nodes[neibNo], cellNo);

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
                    T1 dt = this->computeDt(*src, tempNodes[threadNo][neibNo], cellNo);
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
#endif /* Grid3Ducdsp_h */

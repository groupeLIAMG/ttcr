//
//  Grid3Duisp.h
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

#ifndef ttcr_Grid3Duisp_h
#define ttcr_Grid3Duisp_h

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <vector>

#include "Grid3Dui.h"
#include "Interpolator.h"
#include "Node3Disp.h"
#include "utils.h"

namespace ttcr {
    
    template<typename T1, typename T2>
    class Grid3Duisp : public Grid3Dui<T1,T2,Node3Disp<T1,T2>> {
    public:
        Grid3Duisp(const std::vector<sxyz<T1>>& no,
                   const std::vector<tetrahedronElem<T2>>& tet,
                   const int ns, const size_t nt=1, const int verbose=0) :
        Grid3Dui<T1,T2,Node3Disp<T1,T2>>(no, tet, nt), nsecondary(ns)
        {
            buildGridNodes(no, ns, nt, verbose);
            this->buildGridNeighbors();
            //cout<<this->nodes.size();
        }
        
        ~Grid3Duisp() {
        }
        
        int setSlowness(const std::vector<T1>& s) {
            if ( this->nPrimary != s.size() ) {
                std::cerr << "Error: slowness vectors of incompatible size.";
                return 1;
            }
            for ( size_t n=0; n<this->nPrimary; ++n ) {
                this->nodes[n].setNodeSlowness( s[n] );
            }
            if ( nsecondary>0 ) interpSlownessSecondary();
            return 0;
        }
        size_t getNumberOfNodes() const { return this->nPrimary; }
        int setSlowness(const T1 *s, const size_t ns) {
            if ( this->nPrimary != ns ) {
                std::cerr << "Error: slowness vectors of incompatible size.";
                return 1;
            }
            for ( size_t n=0; n<this->nPrimary; ++n ) {
                this->nodes[n].setNodeSlowness( s[n] );
            }
            if ( nsecondary>0 ) interpSlownessSecondary();
            return 0;
        }
        T1 AddNewSecondaryNodes(const sxyz<T1> & source,const size_t &  Number, const T1 & Radius1,
                                  const T1 & Radius2=0);
        T1 AddNodesInside(const sxyz<T1> & source,const size_t &  Number, const T1 & Radius1)const;
        
        void DelateNodes();
        int raytrace(const std::vector<sxyz<T1>>&,
                     const std::vector<T1>&,
                     const std::vector<sxyz<T1>>&,
                     std::vector<T1>&,
                     const size_t=0) const;
        
        int raytrace(const std::vector<sxyz<T1>>&,
                     const std::vector<T1>&,
                     const std::vector<const std::vector<sxyz<T1>>*>&,
                     std::vector<std::vector<T1>*>&,
                     const size_t=0) const;
        
        int raytrace(const std::vector<sxyz<T1>>&,
                     const std::vector<T1>& ,
                     const std::vector<sxyz<T1>>&,
                     std::vector<T1>&,
                     std::vector<std::vector<sxyz<T1>>>&,
                     const size_t=0) const;
        
        int raytrace(const std::vector<sxyz<T1>>&,
                     const std::vector<T1>&,
                     const std::vector<const std::vector<sxyz<T1>>*>&,
                     std::vector<std::vector<T1>*>&,
                     std::vector<std::vector<std::vector<sxyz<T1>>>*>&,
                     const size_t=0) const;
        
        int raytrace(const std::vector<sxyz<T1>>&,
                     const std::vector<T1>& ,
                     const std::vector<sxyz<T1>>&,
                     std::vector<T1>&,
                     std::vector<std::vector<sxyz<T1>>>&,
                     std::vector<std::vector<siv<T1>>>&,
                     const size_t=0) const;
        int raytrace(const std::vector<sxyz<T1>>&,
                     const std::vector<T1>& ,
                     const std::vector<sxyz<T1>>&,
                     std::vector<T1>&,
                     std::vector<std::vector<sxyz<T1>>>&,
                     std::vector<std::vector<sijv<T1>>>&,
                     const size_t=0) const;
        int raytrace(const std::vector<sxyz<T1>>&,
                     const std::vector<T1>& ,
                     const std::vector<sxyz<T1>>&,
                     std::vector<T1>&,
                     std::vector<std::vector<sxyz<T1>>>&,
                     std::vector<std::vector<sijv<T1>>>&, T1 &,
                     const size_t=0) const;
        int raytrace(const std::vector<sxyz<T1>>& Tx,
                 const std::vector<T1>& t0,
                 const std::vector<sxyz<T1>>& Rx,
                 std::vector<T1>& traveltimes,
                 std::vector<std::vector<sxyz<T1>>>& r_data,
                     T1 & v0,
                 const size_t=0)const;
        int raytrace(const std::vector<sxyz<T1>>& Tx,
                     const std::vector<T1>& t0,
                     const std::vector<sxyz<T1>>& Rx,
                     std::vector<T1>& traveltimes,
                     std::vector<std::vector<sxyz<T1>>>& r_data,
                     T1 & v0,
                     std::vector<std::vector<sijv<T1>>>& m_data,
                     const size_t threadNo=0) const;
        int computeD(const std::vector<sxyz<T1>>& Pts, std::vector<std::vector<sijv<T1>>>& d_data)const;
    private:
        T2 nsecondary;
        
        void buildGridNodes(const std::vector<sxyz<T1>>&,
                            const int, const size_t, const int);
        
        void interpSlownessSecondary();
        
        void initQueue(const std::vector<sxyz<T1>>& Tx,
                       const std::vector<T1>& t0,
                       std::priority_queue<Node3Disp<T1,T2>*,
                       std::vector<Node3Disp<T1,T2>*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<Node3Disp<T1,T2>>& txNodes,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       const size_t threadNo) const;
        
        void prepropagate(const Node3Disp<T1,T2>& node,
                          std::priority_queue<Node3Disp<T1,T2>*, std::vector<Node3Disp<T1,T2>*>,
                          CompareNodePtr<T1>>& queue,
                          std::vector<bool>& inQueue,
                          std::vector<bool>& frozen,
                          size_t threadNo) const;
        
        void propagate(std::priority_queue<Node3Disp<T1,T2>*,
                       std::vector<Node3Disp<T1,T2>*>,
                       CompareNodePtr<T1>>& queue,
                       std::vector<bool>& inQueue,
                       std::vector<bool>& frozen,
                       const size_t threadNo) const;
        
        T1 getTraveltime(const sxyz<T1>& Rx,
                         const std::vector<Node3Disp<T1,T2>>& nodes,
                         const size_t threadNo) const;
        
        T1 getTraveltime(const sxyz<T1>& Rx,
                         const std::vector<Node3Disp<T1,T2>>& nodes,
                         T2& nodeParentRx,
                         T2& cellParentRx,
                         const size_t threadNo) const;
    };
    
    
    template<typename T1, typename T2>
    void Grid3Duisp<T1,T2>::buildGridNodes(const std::vector<sxyz<T1>>& no,
                                           const int nsecondary,
                                           const size_t nt, const int verbose) {
        
        // primary nodes
        for ( T2 n=0; n<no.size(); ++n ) {
            this->nodes[n].setXYZindex( no[n].x, no[n].y, no[n].z, n );
            this->nodes[n].setPrimary(5);
        }
        T2 nNodes = static_cast<T2>(this->nodes.size());
        
        size_t nFaceNodes = 0;
        for ( int n=1; n<=(nsecondary-1); ++n ) nFaceNodes += n;//nsecondary*(nsecondary-1)/2
        
        std::map<std::array<T2,2>,std::vector<T2>> lineMap;
        std::array<T2,2> lineKey;
        typename std::map<std::array<T2,2>,std::vector<T2>>::iterator lineIt;
        
        size_t estLineNo = (this->tetrahedra.size()+this->tetrahedra.size()/10) * 6/2;
        size_t estFaceNo = (this->tetrahedra.size()+this->tetrahedra.size()/10) * 4/2;
        this->nodes.reserve( nNodes + estLineNo*nsecondary + estFaceNo*nFaceNodes );
        
        
        T2 iNodes[4][3] = {
            {0,1,2},  // (relative) indices of nodes of 1st triangle
            {1,2,3},  // (relative) indices of nodes of 2nd triangle
            {0,2,3},  // (relative) indices of nodes of 3rd triangle
            {0,1,3}   // (relative) indices of nodes of 4th triangle
        };
        
        //
        //              1
        //            ,/|`\
        //          ,/  |  `\
        //        ,0    '.   `4
        //      ,/       1     `\
        //    ,/         |       `\
        //   0-----5-----'.--------3
        //    `\.        |      ,/
        //       `\.     |     3
        //          `2.  '. ,/
        //             `\|/
        //               2
        //
        //
        //  triangle 0:  0-1  1-2  2-0     (first occurence of segment underlined)
        //               ---  ---  ---
        //  triangle 1:  1-2  2-3  3-1
        //                    ---  ---
        //  triangle 2:  0-2  2-3  3-0
        //                         ---
        //  triangle 3:  0-1  1-3  3-0
        
        if ( verbose>1 && nsecondary > 0 ) {
            std::cout << '\n';
        }
        
        // edge nodes
        Node3Disp<T1,T2> tmpNode(nt);
        for ( T2 ntet=0; ntet<this->tetrahedra.size(); ++ntet ) {
            
            if ( verbose>1 && nsecondary > 0 ) {
                std::cout << "\r  Building edge nodes: " << (100*ntet)/this->tetrahedra.size() << "%";
                std::cout.flush();
            }
            
            // for each triangle
            for ( T2 ntri=0; ntri<4; ++ntri ) {
                
                // push owner for primary nodes
                this->nodes[ this->tetrahedra[ntet].i[ntri] ].pushOwner( ntet );
                
                if ( nsecondary > 0 ) {
                    // start from ntri to avoid redundancy
                    for ( size_t nl=ntri; nl<3; ++nl ) {
                        
                        lineKey = {this->tetrahedra[ntet].i[ iNodes[ntri][nl] ],
                            this->tetrahedra[ntet].i[ iNodes[ntri][(nl+1)%3] ]};
                        std::sort(lineKey.begin(), lineKey.end());
                        
                        lineIt = lineMap.find( lineKey );
                        if ( lineIt == lineMap.end() ) {
                            // not found, insert new pair
                            lineMap[ lineKey ] = std::vector<T2>(nsecondary);
                        } else {
                            for ( size_t n=0; n<lineIt->second.size(); ++n ) {
                                // setting owners
                                this->nodes[ lineIt->second[n] ].pushOwner( ntet );
                            }
                            continue;
                        }
                        
                        sxyz<T1> d = (no[lineKey[1]]-no[lineKey[0]])/static_cast<T1>(nsecondary+1);
                        
                        for ( size_t n2=0; n2<nsecondary; ++n2 ) {
                            tmpNode.setXYZindex(no[lineKey[0]].x+(1+n2)*d.x,
                                                no[lineKey[0]].y+(1+n2)*d.y,
                                                no[lineKey[0]].z+(1+n2)*d.z,
                                                nNodes );
                            tmpNode.setPrimary(10);
                            tmpNode.SetPrincipals(lineKey[0]);
                            tmpNode.SetPrincipals(lineKey[1]);
                            lineMap[lineKey][n2] = nNodes++;
                            this->nodes.push_back( tmpNode );
                            this->nodes.back().pushOwner( ntet );
                        }
                    }
                }
            }
        }
        
        if ( verbose>1 && nsecondary > 0 ) {
            std::cout << "\r  Building edge nodes: 100%\n";
        }
        
        if ( nsecondary > 1 ) {
            
            std::map<std::array<T2,3>,std::vector<T2>> faceMap;
            std::array<T2,3> faceKey;
            typename std::map<std::array<T2,3>,std::vector<T2>>::iterator faceIt;
            
            int ncut = nsecondary - 1;
            
            for ( T2 ntet=0; ntet<this->tetrahedra.size(); ++ntet ) {
                
                if ( verbose>1 ) {
                    std::cout << "\r  Building face nodes: " << (100*ntet)/this->tetrahedra.size() << "%";
                    std::cout.flush();
                }
                
                // for each triangle
                for ( T2 ntri=0; ntri<4; ++ntri ) {
                    
                    faceKey = {this->tetrahedra[ntet].i[ iNodes[ntri][0] ],
                        this->tetrahedra[ntet].i[ iNodes[ntri][1] ],
                        this->tetrahedra[ntet].i[ iNodes[ntri][2] ]};
                    std::sort(faceKey.begin(), faceKey.end());
                    
                    
                    faceIt = faceMap.find( faceKey );
                    if ( faceIt == faceMap.end() ) {
                        // not found, insert new pair
                        faceMap[ faceKey ] = std::vector<T2>(nFaceNodes);
                    } else {
                        for ( size_t n=0; n<faceIt->second.size(); ++n ) {
                            // setting owners
                            this->nodes[ faceIt->second[n] ].pushOwner( ntet );
                        }
                        continue;
                    }
                    sxyz<T1> d1 = (no[faceKey[1]]-no[faceKey[0]])/static_cast<T1>(nsecondary+1);
                    sxyz<T1> d2 = (no[faceKey[1]]-no[faceKey[2]])/static_cast<T1>(nsecondary+1);
                    
                    size_t ifn = 0;
                    for ( size_t n=0; n<ncut; ++n ) {
                        
                        sxyz<T1> pt1 = no[faceKey[0]]+static_cast<T1>(1+n)*d1;
                        sxyz<T1> pt2 = no[faceKey[2]]+static_cast<T1>(1+n)*d2;
                        
                        size_t nseg = ncut+1-n;
                        
                        sxyz<T1> d = (pt2-pt1)/static_cast<T1>(nseg);
                        
                        for ( size_t n2=0; n2<nseg-1; ++n2 ) {
                            tmpNode.setXYZindex(pt1.x+(1+n2)*d.x,
                                                pt1.y+(1+n2)*d.y,
                                                pt1.z+(1+n2)*d.z,
                                                nNodes );
                            tmpNode.setPrimary(10);
                            tmpNode.SetPrincipals(faceKey[0]);
                            tmpNode.SetPrincipals(faceKey[1]);
                            tmpNode.SetPrincipals(faceKey[2]);
                            faceMap[faceKey][ifn++] = nNodes++;
                            this->nodes.push_back( tmpNode );
                            this->nodes.back().pushOwner( ntet );
                        }
                    }
                }
            }
        }
        if ( verbose>1 && nsecondary > 0 ) {
            std::cout << "\r  Building face nodes: 100%\n";
        }
        
        this->nodes.shrink_to_fit();
    }
    
    
    template<typename T1, typename T2>
    void Grid3Duisp<T1,T2>::interpSlownessSecondary() {
        
        T2 nNodes = this->nPrimary;
        
        std::map<std::array<T2,2>,std::vector<T2>> lineMap;
        std::array<T2,2> lineKey;
        typename std::map<std::array<T2,2>,std::vector<T2>>::iterator lineIt;
        
        size_t nFaceNodes = 0;
        for ( int n=1; n<=(nsecondary-1); ++n ) nFaceNodes += n;
        
        for ( T2 ntet=0; ntet<this->tetrahedra.size(); ++ntet ) {
            
            // for each triangle
            for ( T2 ntri=0; ntri<4; ++ntri ) {
                
                // start from ntri to avoid redundancy
                for ( size_t nl=ntri; nl<3; ++nl ) {
                    
                    lineKey = {this->tetrahedra[ntet].i[ iNodes[ntri][nl] ],
                        this->tetrahedra[ntet].i[ iNodes[ntri][(nl+1)%3] ]};
                    std::sort(lineKey.begin(), lineKey.end());
                    
                    lineIt = lineMap.find( lineKey );
                    if ( lineIt == lineMap.end() ) {
                        // not found, insert new pair
                        lineMap[ lineKey ] = std::vector<T2>(nsecondary);
                    } else {
                        continue;
                    }
                    T1 slope=(1.0/this->nodes[lineKey[1]].getNodeSlowness() - 1.0/this->nodes[lineKey[0]].getNodeSlowness())/this->nodes[lineKey[1]].getDistance(this->nodes[lineKey[0]]);
                    for ( size_t n2=0; n2<nsecondary; ++n2 ) {
                        T1 s=1.0/(1.0/this->nodes[lineKey[0]].getNodeSlowness() + slope * this->nodes[nNodes].getDistance(this->nodes[lineKey[0]]));
                        this->nodes[nNodes].setNodeSlowness( s );
                        lineMap[lineKey][n2] = nNodes++;
                    }
                }
            }
        }
        
        
        
        if ( nsecondary > 1 ) {
            
            std::map<std::array<T2,3>,std::vector<T2>> faceMap;
            std::array<T2,3> faceKey;
            typename std::map<std::array<T2,3>,std::vector<T2>>::iterator faceIt;
            
            int ncut = nsecondary - 1;
            
            for ( T2 ntet=0; ntet<this->tetrahedra.size(); ++ntet ) {
                
                // for each triangle
                for ( T2 ntri=0; ntri<4; ++ntri ) {
                    
                    faceKey = {this->tetrahedra[ntet].i[ iNodes[ntri][0] ],
                        this->tetrahedra[ntet].i[ iNodes[ntri][1] ],
                        this->tetrahedra[ntet].i[ iNodes[ntri][2] ]};
                    std::sort(faceKey.begin(), faceKey.end());
                    
                    
                    faceIt = faceMap.find( faceKey );
                    if ( faceIt == faceMap.end() ) {
                        // not found, insert new pair
                        faceMap[ faceKey ] = std::vector<T2>(nFaceNodes);
                    } else {
                        continue;
                    }
                    
                    std::vector<Node3Disp<T1,T2>*> inodes;
                    inodes.push_back( &(this->nodes[faceKey[0]]) );
                    inodes.push_back( &(this->nodes[faceKey[1]]) );
                    inodes.push_back( &(this->nodes[faceKey[2]]) );
                    
                    size_t ifn = 0;
                    for ( size_t n=0; n<ncut; ++n ) {
                        size_t nseg = ncut+1-n;
                        for ( size_t n2=0; n2<nseg-1; ++n2 ) {
                            T1 s = Interpolator<T1>::bilinearTriangleVel(this->nodes[nNodes], this->nodes[faceKey[0]], this->nodes[faceKey[1]], this->nodes[faceKey[2]]);
                            this->nodes[nNodes].setNodeSlowness( s );
                            faceMap[faceKey][ifn++] = nNodes++;
                            
                        }
                    }
                }
            }
        }
    }
    template<typename T1, typename T2>
    T1 Grid3Duisp<T1,T2>::AddNewSecondaryNodes(const sxyz<T1> & source,const size_t & SecondNodesN, const T1 & Radius1,const T1 & Radius2){
        std::vector<T2> Tetrahedrons;
        for(T2 tet=0; tet<this->neighbors.size();++tet){
            sxyz<T1> Centroid=this->nodes[this->neighbors[tet][0]]+this->nodes[this->neighbors[tet][1]];
            Centroid+=this->nodes[this->neighbors[tet][2]]+this->nodes[this->neighbors[tet][3]];
            Centroid/=4;
            if (source.getDistance(Centroid)<=Radius1 && source.getDistance(Centroid)>=Radius2)
                Tetrahedrons.push_back(tet);
        }
        if (Tetrahedrons.size()==0)
            return 0.0;
        std::set<T2> AdjacentCells(Tetrahedrons.begin(),Tetrahedrons.end());
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
        
        // edges vector
        T1 EdgesCoord [(nsecondary+1)*SecondNodesN];
        T2 n=0;
        for (size_t i=0;i<((nsecondary+1)*(SecondNodesN+1)-1);++i){
            if ((i+1)%(SecondNodesN+1)!=0){
                EdgesCoord[n++]=(i+1.0)/((nsecondary+1)*(SecondNodesN+1));
            }
        }
        T2 nFaceNodes=0.5*(pow(((nsecondary+1)*SecondNodesN+nsecondary),2)-(nsecondary+1)*SecondNodesN-pow(nsecondary,2));
        T1 FacesCoord [nFaceNodes][2];
        n=0;
        for (size_t i=0;i<((nsecondary+1)*(SecondNodesN+1)-2);++i){
            for (size_t j=0;j<((nsecondary+1)*(SecondNodesN+1)-i-2);++j){
                if (((i+1)%(SecondNodesN+1)!=0)||((j+1)%(SecondNodesN+1)!=0)){
                    FacesCoord[n][0]=(i+1.0)/((nsecondary+1)*(SecondNodesN+1));
                    FacesCoord[n++][1]=(j+1.0)/((nsecondary+1)*(SecondNodesN+1));
                }
            }
        }
        T2 nNodes = static_cast<T2>(this->nodes.size());
        for (auto tetN=Tetrahedrons.begin();tetN!=Tetrahedrons.end();tetN++){
            //  adjacent cells to the tetrahedron where new nodes will be added
            for(size_t i=0;i<4;++i){
                T2 Vertex=this->neighbors[*tetN][i];
                for(size_t c=0;c<this->nodes[Vertex].getOwners().size();++c){
                    AdjacentCells.insert(this->nodes[Vertex].getOwners()[c]);
                }
            }

            // for each triangle
            for ( T2 ntri=0; ntri<4; ++ntri ) {
                for ( size_t nl=ntri; nl<3; ++nl ) {
                    lineKey = {this->tetrahedra[*tetN].i[ iNodes[ntri][nl] ],
                        this->tetrahedra[*tetN].i[ iNodes[ntri][(nl+1)%3] ]};
                    std::sort(lineKey.begin(), lineKey.end());
                    
                    lineIt = lineMap.find( lineKey );
                    if ( lineIt == lineMap.end() ) {
                        // not found, insert new pair
                        lineMap[lineKey] = std::vector<T2>((nsecondary+1)*SecondNodesN);
                    } else {
                        for ( size_t n=0; n<lineIt->second.size(); ++n ) {
                            // setting owners
                            this->nodes[ lineIt->second[n] ].pushOwner(*tetN);
                            this->neighbors[*tetN].push_back(lineIt->second[n]);
                        }
                        continue;
                    }
                    sxyz<T1> Vect =sxyz<T1>(this->nodes[lineKey[1]])-sxyz<T1>(this->nodes[lineKey[0]]);
                    for ( size_t n2=0; n2<((nsecondary+1)*SecondNodesN); ++n2 ) {
                        Node3Disp<T1,T2> tmpNode(this->nThreads);
                        tmpNode.setXYZindex(this->nodes[lineKey[0]].getX()+EdgesCoord[n2]*Vect.x,
                                            this->nodes[lineKey[0]].getY()+EdgesCoord[n2]*Vect.y,
                                            this->nodes[lineKey[0]].getZ()+EdgesCoord[n2]*Vect.z,
                                            nNodes );
                        tmpNode.setPrimary(10);
                        tmpNode.SetPrincipals(lineKey[0]);
                        tmpNode.SetPrincipals(lineKey[1]);
                        this->neighbors[*tetN].push_back(nNodes);
                        lineMap[lineKey][n2] = nNodes++;
                        this->nodes.push_back( tmpNode );
                        this->nodes.back().pushOwner( *tetN );
                        // interpolate slowness
                        T1 d1,d0;
                        d0=this->nodes[lineKey[0]].getDistance(this->nodes.back());
                        d1=this->nodes[lineKey[1]].getDistance(this->nodes.back());
                        T1 s=1.0/((d0*1.0/this->nodes[lineKey[1]].getNodeSlowness()+d1*1.0/this->nodes[lineKey[0]].getNodeSlowness())/(d1+d0));
                        this->nodes.back().setNodeSlowness(s);
                    }
                    
                }
                faceKey = {this->tetrahedra[*tetN].i[ iNodes[ntri][0] ],
                    this->tetrahedra[*tetN].i[ iNodes[ntri][1] ],
                    this->tetrahedra[*tetN].i[ iNodes[ntri][2] ]};
                std::sort(faceKey.begin(), faceKey.end());
                
                
                faceIt = faceMap.find( faceKey );
                if ( faceIt == faceMap.end() ) {
                    // not found, insert new pair
                    faceMap[ faceKey ] = std::vector<T2>(nFaceNodes);
                } else {
                    for ( size_t n=0; n<faceIt->second.size(); ++n ) {
                        // setting owners
                        this->nodes[ faceIt->second[n] ].pushOwner(*tetN);
                        this->neighbors[*tetN].push_back(faceIt->second[n]);
                    }
                    continue;
                }
                sxyz<T1> Vect1 =sxyz<T1>(this->nodes[faceKey[2]])-sxyz<T1>(this->nodes[faceKey[0]]);
                sxyz<T1> Vect2 =sxyz<T1>(this->nodes[faceKey[1]])-sxyz<T1>(this->nodes[faceKey[0]]);
                for ( size_t n=0; n<nFaceNodes; ++n ) {
                    //std::cout<<FacesCoord[n][0]<<" ;"<< FacesCoord[n][1]<<std::endl;
                    Node3Disp<T1,T2> tmpNode(this->nThreads);
                    tmpNode.setXYZindex(this->nodes[faceKey[0]].getX()+FacesCoord[n][0]*Vect1.x+FacesCoord[n][1]*Vect2.x,
                                        this->nodes[faceKey[0]].getY()+FacesCoord[n][0]*Vect1.y+FacesCoord[n][1]*Vect2.y,
                                        this->nodes[faceKey[0]].getZ()+FacesCoord[n][0]*Vect1.z+FacesCoord[n][1]*Vect2.z,
                                        nNodes );
                    tmpNode.setPrimary(10);
                    tmpNode.SetPrincipals(faceKey[0]);
                    tmpNode.SetPrincipals(faceKey[1]);
                    tmpNode.SetPrincipals(faceKey[2]);
                    this->neighbors[*tetN].push_back(nNodes);
                    this->nodes.push_back(tmpNode);
                    this->nodes.back().pushOwner(*tetN);
                    std::vector<Node3Disp<T1,T2>*> inodes;
                    inodes.push_back( &(this->nodes[faceKey[0]]) );
                    inodes.push_back( &(this->nodes[faceKey[1]]) );
                    inodes.push_back( &(this->nodes[faceKey[2]]) );
                    T1 s = Interpolator<T1>::bilinearTriangleVel(this->nodes[nNodes],(this->nodes[faceKey[0]]) , (this->nodes[faceKey[1]]) , (this->nodes[faceKey[2]]));
                    this->nodes[nNodes].setNodeSlowness( s );
                    faceMap[faceKey][n] = nNodes++;
                }
            }
            
        }
        for (auto tet=Tetrahedrons.begin();tet!=Tetrahedrons.end();++tet){
            AdjacentCells.erase(*tet);
        }
        for(auto Adj=AdjacentCells.begin(); Adj!=AdjacentCells.end();++Adj){
            for ( T2 ntri=0; ntri<4; ++ntri ) {
                for ( size_t nl=ntri; nl<3; ++nl ) {
                    lineKey = {this->tetrahedra[*Adj].i[ iNodes[ntri][nl] ],
                        this->tetrahedra[*Adj].i[ iNodes[ntri][(nl+1)%3] ]};
                    std::sort(lineKey.begin(), lineKey.end());
                    
                    lineIt = lineMap.find( lineKey );
                    if ( lineIt != lineMap.end() ) {
                        // setting owners
                        for ( size_t n=0; n<lineIt->second.size(); ++n ) {
                            // setting owners
                            this->nodes[ lineIt->second[n] ].pushOwner(*Adj);
                            this->neighbors[*Adj].push_back(lineIt->second[n]);
                        }
                    }
                }
                faceKey = {this->tetrahedra[*Adj].i[ iNodes[ntri][0] ],
                    this->tetrahedra[*Adj].i[ iNodes[ntri][1] ],
                    this->tetrahedra[*Adj].i[ iNodes[ntri][2] ]};
                std::sort(faceKey.begin(), faceKey.end());
                
                
                faceIt = faceMap.find( faceKey );
                if ( faceIt != faceMap.end() ) {
                    for ( size_t n=0; n<faceIt->second.size(); ++n ) {
                        // setting owners
                        this->nodes[ faceIt->second[n] ].pushOwner(*Adj);
                        this->neighbors[*Adj].push_back(faceIt->second[n]);
                    }
                   
                }
            }
        }
        //cout<<this->nodes.size();
        return 0.0;
    }
    
    template<typename T1, typename T2>
    void Grid3Duisp<T1,T2>::DelateNodes(){
        T2 PermanentNodes=2*nsecondary*(nsecondary-1)+6*nsecondary+4;
        T2 firstAddedNodes=static_cast<T2>(this->nodes.size());
        for (auto tet=this->neighbors.begin();tet!=this->neighbors.end();++tet){
            if ((*tet).size()<=PermanentNodes)
                continue;
            firstAddedNodes=(*tet)[PermanentNodes]<firstAddedNodes?(*tet)[PermanentNodes]:firstAddedNodes;
            (*tet).erase((*tet).begin()+PermanentNodes,(*tet).end());
        }
        this->nodes.erase(this->nodes.begin()+firstAddedNodes,this->nodes.end());
    }
    
    
    template<typename T1,typename T2>
    T1 Grid3Duisp<T1,T2>::AddNodesInside(const sxyz<T1> & source,const size_t &  Number, const T1 & Radius) const{
        std::vector<T1> RelativeCoordinates(Number);
        for (size_t p=0;p<Number;p++)
            RelativeCoordinates[p]=(p+1.0)/(Number+1.0);
        RelativeCoordinates.erase(std::remove(RelativeCoordinates.begin(),RelativeCoordinates.end(),(2.0/3.0)),RelativeCoordinates.end());
        T2 nNodes = static_cast<T2>(this->nodes.size());
        for(T2 tet=0; tet<this->neighbors.size();++tet){
            sxyz<T1> Centroid=this->nodes[this->neighbors[tet][0]]+this->nodes[this->neighbors[tet][1]];
            Centroid+=this->nodes[this->neighbors[tet][2]]+this->nodes[this->neighbors[tet][3]];
            Centroid/=4;
            if (source.getDistance(Centroid)<=Radius){
                Node3Disp<T1,T2> tmpNode (this->nThreads);
                tmpNode.setXYZindex(Centroid.x,Centroid.y,Centroid.z,nNodes);
                tmpNode.setPrimary(10);
                tmpNode.SetPrincipals(this->neighbors[tet][0]);
                tmpNode.SetPrincipals(this->neighbors[tet][1]);
                tmpNode.SetPrincipals(this->neighbors[tet][2]);
                tmpNode.SetPrincipals(this->neighbors[tet][3]);
                this->neighbors[tet].push_back(nNodes);
                this->nodes.push_back(tmpNode);
                this->nodes.back().pushOwner(tet);
                std::vector<Node3Disp<T1,T2>*> inodes;
                inodes.push_back( &(this->nodes[this->neighbors[tet][0]]) );
                inodes.push_back( &(this->nodes[this->neighbors[tet][1]]) );
                inodes.push_back( &(this->nodes[this->neighbors[tet][2]]) );
                inodes.push_back( &(this->nodes[this->neighbors[tet][3]]) );
                T1 s = Interpolator<T1>::inverseDistance(this->nodes[nNodes], inodes);
                this->nodes[nNodes++].setNodeSlowness( s );
                
                for (size_t i=0;i<4;i++){
                    sxyz<T1> Median={Centroid.x-this->nodes[this->neighbors[tet][i]].getX(),
                        Centroid.y-this->nodes[this->neighbors[tet][i]].getY(),Centroid.z-this->nodes[this->neighbors[tet][i]].getZ()};
                    Median*=1.5;
                    for(auto nf=RelativeCoordinates.begin();nf!=RelativeCoordinates.end();nf++){
                        Node3Disp<T1,T2> tmpNode(this->nThreads);
                        tmpNode.setXYZindex(this->nodes[this->neighbors[tet][i]].getX()+(*nf)*Median.x,this->nodes[this->neighbors[tet][i]].getY()+(*nf)*Median.y,this->nodes[this->neighbors[tet][i]].getZ()+(*nf)*Median.z,nNodes);
                        tmpNode.setPrimary(10);
                        tmpNode.SetPrincipals(this->neighbors[tet][0]);
                        tmpNode.SetPrincipals(this->neighbors[tet][1]);
                        tmpNode.SetPrincipals(this->neighbors[tet][2]);
                        tmpNode.SetPrincipals(this->neighbors[tet][3]);
                        this->neighbors[tet].push_back(nNodes);
                        this->nodes.push_back(tmpNode);
                        this->nodes.back().pushOwner(tet);
                        T1 s = Interpolator<T1>::inverseDistance(this->nodes.back(), inodes);
                        this->nodes[nNodes++].setNodeSlowness( s );
                        
                    }
                    
                }
//                for (auto nei=this->neighbors[tet].begin();nei!=this->neighbors[tet].end();++nei){
//                    std::cout<<std::setprecision(10)<<this->nodes[*nei].getX()<<","<<this->nodes[*nei].getY()<<","<<this->nodes[*nei].getZ()<<";"<<std::endl;
//                }
                
            }
            
        }
        return 0.0;
    }
    
    template<typename T1, typename T2>
    int Grid3Duisp<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                    const std::vector<T1>& t0,
                                    const std::vector<sxyz<T1>>& Rx,
                                    std::vector<T1>& traveltimes,
                                    const size_t threadNo) const {
        
        if ( this->checkPts(Tx) == 1 ) return 1;
        if ( this->checkPts(Rx) == 1 ) return 1;
        
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }
        
        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node3Disp<T1,T2>*, std::vector<Node3Disp<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );
        
        std::vector<Node3Disp<T1,T2>> txNodes;
        std::vector<bool> inQueue( this->nodes.size(), false );
        std::vector<bool> frozen( this->nodes.size(), false );
        
        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
        
        propagate(queue, inQueue, frozen, threadNo);
        
        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        
//        for (size_t n=0; n<Rx.size(); ++n) {
//            traveltimes[n] = this->getTraveltime(Rx[n], this->nodes, threadNo);
//        }
        T2 StandartSPM=1.0;
        T1 Slowness_Srcs[Tx.size()];
        T2 Cells_Tx [Tx.size()];
        if (StandartSPM==1.0){
            for(size_t i=0;i<Tx.size();++i){
                Cells_Tx[i] =this->getCellNo(Tx[i]);
                Slowness_Srcs[i]=this->computeSlowness(Tx[i]);
            }
        }
        
        for (size_t n=0; n<Rx.size(); ++n) {
            T1 t=0.0;
            if (StandartSPM==1.0){
                bool find=false;
                for(size_t i=0;i<Tx.size();++i){
                    if (Cells_Tx[i]==this->getCellNo(Rx[n])){
                        t=Tx[i].getDistance(Rx[n])*0.5*(Slowness_Srcs[i]+this->computeSlowness(Rx[n]));
                        find=true;
                        break;
                    }
                }
                if (!find){
                    t= this->getTraveltime(Rx[n], this->nodes, threadNo);
                }
            }else{
                std::vector<sxyz<T1>> r_data;
                this->getRaypathABM(Tx, Rx[n],r_data,threadNo);
                T1 s2=this->computeSlowness(r_data[0]);
                for(size_t no=0;no<r_data.size()-1;++no){
                    T1 s1=this->computeSlowness(r_data[no+1]);
                    t+=r_data[no].getDistance(r_data[no+1])*0.5*(s1+s2);
                    s2=s1;
                }
                
            }
            traveltimes[n] =t;
        }

        return 0;
    }
    
    template<typename T1, typename T2>
    int Grid3Duisp<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                    const std::vector<T1>& t0,
                                    const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                                    std::vector<std::vector<T1>*>& traveltimes,
                                    const size_t threadNo) const {
        
        if ( this->checkPts(Tx) == 1 ) return 1;
        for ( size_t n=0; n<Rx.size(); ++n )
            if ( this->checkPts(*Rx[n]) == 1 ) return 1;
        
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }
        
        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node3Disp<T1,T2>*, std::vector<Node3Disp<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );
        
        std::vector<Node3Disp<T1,T2>> txNodes;
        std::vector<bool> inQueue( this->nodes.size(), false );
        std::vector<bool> frozen( this->nodes.size(), false );
        
        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
        
        propagate(queue, inQueue, frozen, threadNo);
        
        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        
        for (size_t nr=0; nr<Rx.size(); ++nr) {
            traveltimes[nr]->resize( Rx[nr]->size() );
            for (size_t n=0; n<Rx[nr]->size(); ++n)
                (*traveltimes[nr])[n] = this->getTraveltime((*Rx[nr])[n], this->nodes, threadNo);
        }
        
        return 0;
    }
    
    
    template<typename T1, typename T2>
    int Grid3Duisp<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                    const std::vector<T1>& t0,
                                    const std::vector<sxyz<T1>>& Rx,
                                    std::vector<T1>& traveltimes,
                                    std::vector<std::vector<sxyz<T1>>>& r_data,
                                    const size_t threadNo) const {
        
        if ( this->checkPts(Tx) == 1 ) return 1;
        if ( this->checkPts(Rx) == 1 ) return 1;
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }
        
        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node3Disp<T1,T2>*, std::vector<Node3Disp<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );
        
        std::vector<Node3Disp<T1,T2>> txNodes;
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
        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        /////////////////
       
        std::ofstream myfile;
        myfile.open ("StandartSPM.dat");
        myfile.close();
        //myfile.open ("CreatedSPM.dat");
        //myfile.close();
        //////////////////
        T2 StandartSPM=0.0;
        T1 Slowness_Srcs[Tx.size()];
        T2 Cells_Tx [Tx.size()];
        if (StandartSPM==1.0){
            for(size_t i=0;i<Tx.size();++i){
                Cells_Tx[i] =this->getCellNo(Tx[i]);
                Slowness_Srcs[i]=this->computeSlowness(Tx[i]);
            }
        }

        for (size_t n=0; n<Rx.size(); ++n) {
        
            //this->getRaypathBLIT(Tx, Rx[n], r_data[n],threadNo);
            //this->getRaypath(Tx, Rx[n], r_data[n],threadNo);
            //this->getRaypathABM(Tx, Rx[n], r_data[n],threadNo);
           //this->getRaypath_ho(Tx, Rx[n], r_data[n],threadNo);
            ///////////////////////////////////////////////////////////////////////////////////////////////
            T1 t=0.0;
            if (StandartSPM==1.0){
                bool find=false;
                for(size_t i=0;i<Tx.size();++i){
                    if (Cells_Tx[i]==this->getCellNo(Rx[n])){
                        t=Tx[i].getDistance(Rx[n])*0.5*(Slowness_Srcs[i]+this->computeSlowness(Rx[n]));
                        find=true;
                        break;
                    }
                }
                if (!find){
                    //cout<<" in cell source"<< endl;
                    t= this->getTraveltime(Rx[n], this->nodes, threadNo);
                }
                
            }else{
                this->getRaypathABM(Tx, Rx[n], r_data[n],threadNo);
                T1 s2=this->computeSlowness(r_data[n][0]);
        
                //T1 s2=this->Highorder_computeSlowness(r_data[n][0]);
                for(size_t no=0;no<r_data[n].size()-1;++no){
                    //T1 s1=this->Highorder_computeSlowness(r_data[n][no+1]);
                    T1 s1=this->computeSlowness(r_data[n][no+1]);
                    t+=r_data[n][no].getDistance(r_data[n][no+1])*0.5*(s1+s2);
                    s2=s1;
                }

            }
//            for(T2 sr=0;sr<Tx.size();++sr){
//                if(r_data[n][0]==Tx[sr])
//                        t=t0[sr];
//                }
            traveltimes[n] =t;
            ///////////////////////////////////////////////////////////////////////////////////////////////
            std::ofstream outfile;


            outfile.open("CreatedSPM.dat", std::ios_base::app);
            outfile << t<<"\n";
            outfile.close();

            ///////////////////////////////////////////////////////////////////////////////////////////////
        }
        return 0;
    }
    template<typename T1, typename T2>
    int Grid3Duisp<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                    const std::vector<T1>& t0,
                                    const std::vector<sxyz<T1>>& Rx,
                                    std::vector<T1>& traveltimes,
                                    std::vector<std::vector<sxyz<T1>>>& r_data, T1 &v0,
                                    const size_t threadNo) const {

        if ( this->checkPts(Tx) == 1 ) return 1;
        if ( this->checkPts(Rx) == 1 ) return 1;
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }

        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node3Disp<T1,T2>*, std::vector<Node3Disp<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );

        std::vector<Node3Disp<T1,T2>> txNodes;
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
        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        v0 = 0.0;
        for ( size_t n=0; n<Tx.size(); ++n ) {
            v0 += this->computeSlowness( Tx[n] );
        }
        v0 = Tx.size() / v0;
        /////////////////
        T2 StandartSPM=0.0;
        T1 Slowness_Srcs[Tx.size()];
        T2 Cells_Tx [Tx.size()];
        if (StandartSPM==1.0){
            for(size_t i=0;i<Tx.size();++i){
                Cells_Tx[i] =this->getCellNo(Tx[i]);
                Slowness_Srcs[i]=this->computeSlowness(Tx[i]);
            }
        }
        
        for (size_t n=0; n<Rx.size(); ++n) {

            T1 t=0.0;
            if (StandartSPM==1.0){
                bool find=false;
                for(size_t i=0;i<Tx.size();++i){
                    if (Cells_Tx[i]==this->getCellNo(Rx[n])){
                        t=Tx[i].getDistance(Rx[n])*0.5*(Slowness_Srcs[i]+this->computeSlowness(Rx[n]));
                        find=true;
                        break;
                    }
                }
                if (!find){
            
                    t= this->getTraveltime(Rx[n], this->nodes, threadNo);
                }
                
            }else{
                this->getRaypathABM(Tx, Rx[n], r_data[n],threadNo);
                T1 s2=this->computeSlowness(r_data[n][0]);

                for(size_t no=0;no<r_data[n].size()-1;++no){
                    T1 s1=this->computeSlowness(r_data[n][no+1]);
                    t+=r_data[n][no].getDistance(r_data[n][no+1])*0.5*(s1+s2);
                    s2=s1;
                }
                
            }
            traveltimes[n] =t;
        }
        return 0;
    }
    template<typename T1, typename T2>
    int Grid3Duisp<T1,T2>::computeD(const std::vector<sxyz<T1>> &Pts,
                                    std::vector<std::vector<sijv<T1>>> &d_data) const{
        
        if (d_data.size()!=Pts.size()) d_data.resize(Pts.size());
        for (size_t i=0;i<Pts.size();++i){
            d_data[i].resize(0);
        }
        for(size_t np=0;np<Pts.size();++np){
            bool found=false;
            for(size_t nn=0;nn<this->nodes.size();++nn){
                if (this->nodes[nn].getDistance(Pts[np])<small){
                    found=true;
                    d_data[np].push_back({np,nn,1.0});
                }
            }
            if (!found){
                T2 CellNO=this->getCellNo(Pts[np]);
                if (CellNO==std::numeric_limits<T2>::max()) return 1;
                std::array<T1,4> Weights;
                T1 sum (0.0);
                for(size_t n=0;n<4;++n){
                    Weights[n]=1.0/this->nodes[this->neighbors[CellNO][n]].getDistance(Pts[np]);
                    sum+=Weights[n];
                }
                for(size_t n=0;n<4;++n){
                    Weights[n]/=sum;
                    d_data[np].push_back({np,this->neighbors[CellNO][n],Weights[n]});
                }
            }
        }
        
        return 0;
    }
    
    template<typename T1, typename T2>
    int Grid3Duisp<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                    const std::vector<T1>& t0,
                                    const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                                    std::vector<std::vector<T1>*>& traveltimes,
                                    std::vector<std::vector<std::vector<sxyz<T1>>>*>& r_data,
                                    const size_t threadNo) const {
        
        if ( this->checkPts(Tx) == 1 ) return 1;
        for ( size_t n=0; n<Rx.size(); ++n )
            if ( this->checkPts(*Rx[n]) == 1 ) return 1;
        
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }
        
        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node3Disp<T1,T2>*, std::vector<Node3Disp<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );
        
        std::vector<Node3Disp<T1,T2>> txNodes;
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
                
                (*traveltimes[nr])[n] = this->getTraveltime((*Rx[nr])[n], this->nodes,
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
                
                // Rx are in nodes (not txNodes)
                std::vector<Node3Disp<T1,T2>> *node_p;
                node_p = &(this->nodes);
                
                std::vector<sxyz<T1>> r_tmp;
                T2 iChild, iParent = nodeParentRx;
                sxyz<T1> child;
                
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
        return 0;
    }
    
    template<typename T1, typename T2>
    int Grid3Duisp<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                    const std::vector<T1>& t0,
                                    const std::vector<sxyz<T1>>& Rx,
                                    std::vector<T1>& traveltimes,
                                    std::vector<std::vector<sxyz<T1>>>& r_data,
                                    std::vector<std::vector<siv<T1>>>& l_data,
                                    const size_t threadNo) const {
        
        //	std::cout << "   running in thread no " << threadNo << std::endl;
        if ( checkPts(Tx) == 1 ) return 1;
        if ( checkPts(Rx) == 1 ) return 1;
        
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }
        
        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node3Disp<T1,T2>*, std::vector<Node3Disp<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );
        
        std::vector<Node3Disp<T1,T2>> txNodes;
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
            
            // Rx are in nodes (not txNodes)
            std::vector<Node3Disp<T1,T2>> *node_p;
            node_p = &this->nodes;
            
            std::vector<sxyz<T1>> r_tmp;
            T2 iChild, iParent = nodeParentRx;
            sxyz<T1> child;
            siv<T1> cell;
            
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
        return 0;
    }
    
    template<typename T1, typename T2>
    int Grid3Duisp<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                    const std::vector<T1>& t0,
                                    const std::vector<sxyz<T1>>& Rx,
                                    std::vector<T1>& traveltimes,
                                    std::vector<std::vector<sxyz<T1>>>& r_data,
                                    std::vector<std::vector<sijv<T1>>>& m_data,
                                    const size_t threadNo) const {
        
        //    std::cout << "   running in thread no " << threadNo << std::endl;
        if (this->checkPts(Tx) == 1 ) return 1;
        if (this->checkPts(Rx) == 1 ) return 1;
        
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }
        
        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node3Disp<T1,T2>*, std::vector<Node3Disp<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );
        
        std::vector<Node3Disp<T1,T2>> txNodes;
        std::vector<bool> inQueue( this->nodes.size(), false );
        std::vector<bool> frozen( this->nodes.size(), false );
        
        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
        
        propagate(queue, inQueue, frozen, threadNo);
        
        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        if ( m_data.size() != Rx.size() ) {
            m_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<m_data.size(); ++ni ) {
            m_data[ni].resize( 0 );
        }
        if ( r_data.size() != Rx.size() ) {
            r_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<r_data.size(); ++ni ) {
            r_data[ni].resize( 0 );
        }

        
        for (size_t n=0; n<Rx.size(); ++n) {
            traveltimes[n] = getTraveltime(Rx[n], this->nodes,threadNo);
            this->getRaypath_ho(Tx, Rx[n], r_data[n], m_data[n], n, threadNo);
        }
        
        return 0;
    }
    template<typename T1, typename T2>
    int Grid3Duisp<T1,T2>::raytrace(const std::vector<sxyz<T1>>& Tx,
                                    const std::vector<T1>& t0,
                                    const std::vector<sxyz<T1>>& Rx,
                                    std::vector<T1>& traveltimes,
                                    std::vector<std::vector<sxyz<T1>>>& r_data,
                                     T1 & v0,
                                    std::vector<std::vector<sijv<T1>>>& m_data,
                                    const size_t threadNo) const {
        
        //    std::cout << "   running in thread no " << threadNo << std::endl;
        if (this->checkPts(Tx) == 1 ) return 1;
        if (this->checkPts(Rx) == 1 ) return 1;
        
        for ( size_t n=0; n<this->nodes.size(); ++n ) {
            this->nodes[n].reinit( threadNo );
        }
        
        CompareNodePtr<T1> cmp(threadNo);
        std::priority_queue< Node3Disp<T1,T2>*, std::vector<Node3Disp<T1,T2>*>,
        CompareNodePtr<T1>> queue( cmp );
        
        std::vector<Node3Disp<T1,T2>> txNodes;
        std::vector<bool> inQueue( this->nodes.size(), false );
        std::vector<bool> frozen( this->nodes.size(), false );
        
        initQueue(Tx, t0, queue, txNodes, inQueue, frozen, threadNo);
        
        propagate(queue, inQueue, frozen, threadNo);
        
        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        if ( m_data.size() != Rx.size() ) {
            m_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<m_data.size(); ++ni ) {
            m_data[ni].resize( 0 );
        }
        if ( r_data.size() != Rx.size() ) {
            r_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<r_data.size(); ++ni ) {
            r_data[ni].resize( 0 );
        }
        v0 = 0.0;
        for ( size_t n=0; n<Tx.size(); ++n ) {
            v0 += this->computeSlowness( Tx[n] );
        }
        v0 = Tx.size() / v0;

        
        T2 StandartSPM=0.0;
        T1 Slowness_Srcs[Tx.size()];
        T2 Cells_Tx [Tx.size()];
        if (StandartSPM==1.0){
            for(size_t i=0;i<Tx.size();++i){
                Cells_Tx[i] =this->getCellNo(Tx[i]);
                Slowness_Srcs[i]=this->computeSlowness(Tx[i]);
            }
        }
        
        for (size_t n=0; n<Rx.size(); ++n) {
            
            T1 t=0.0;
            if (StandartSPM==1.0){
                bool find=false;
                for(size_t i=0;i<Tx.size();++i){
                    if (Cells_Tx[i]==this->getCellNo(Rx[n])){
                        t=Tx[i].getDistance(Rx[n])*0.5*(Slowness_Srcs[i]+this->computeSlowness(Rx[n]));
                        find=true;
                        break;
                    }
                }
                if (!find){
                    
                    t= this->getTraveltime(Rx[n], this->nodes, threadNo);
                }
                
            }else{
                this->getRaypath_ho(Tx, Rx[n], r_data[n],m_data[n],n,threadNo);
                T1 s2=this->computeSlowness(r_data[n][0]);
                for(size_t no=0;no<r_data[n].size()-1;++no){
                    T1 s1=this->computeSlowness(r_data[n][no+1]);
                    t+=r_data[n][no].getDistance(r_data[n][no+1])*0.5*(s1+s2);
                    s2=s1;
                }
                
            }
            traveltimes[n] =t;
        }
        
        return 0;
    }
    template<typename T1, typename T2>
    void Grid3Duisp<T1,T2>::initQueue(const std::vector<sxyz<T1>>& Tx,
                                      const std::vector<T1>& t0,
                                      std::priority_queue<Node3Disp<T1,T2>*,
                                      std::vector<Node3Disp<T1,T2>*>,
                                      CompareNodePtr<T1>>& queue,
                                      std::vector<Node3Disp<T1,T2>>& txNodes,
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
                txNodes.push_back( Node3Disp<T1,T2>(t0[n], Tx[n].x, Tx[n].y, Tx[n].z,
                                                    this->nThreads, threadNo));
                T2 CellNeighbor=this->getCellNo(Tx[n]);
                txNodes.back().pushOwner( CellNeighbor );
                txNodes.back().setGridIndex( static_cast<T2>(this->nodes.size()+
                                                             txNodes.size()-1) );
                T1 s=Interpolator<T1>::TrilinearTriangleVel(txNodes.back(), this->nodes[this->neighbors[CellNeighbor][0]], this->nodes[this->neighbors[CellNeighbor][1]], this->nodes[this->neighbors[CellNeighbor][2]], this->nodes[this->neighbors[CellNeighbor][3]]);
                txNodes.back().setNodeSlowness(s);
                frozen.push_back( true );
                
                //prepropagate(txNodes.back(), queue, inQueue, frozen, threadNo); // See description in the function declaration
                
                	queue.push( &(txNodes.back()) );	//Don't use if prepropagate is used
                	inQueue.push_back( true );			//Don't use if prepropagate is used
                
            }
        }
    }
    
    
    template<typename T1, typename T2>
    void Grid3Duisp<T1,T2>::prepropagate(const Node3Disp<T1,T2>& node,
                                         std::priority_queue<Node3Disp<T1,T2>*,
                                         std::vector<Node3Disp<T1,T2>*>,
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
    void Grid3Duisp<T1,T2>::propagate(std::priority_queue<Node3Disp<T1,T2>*,
                                      std::vector<Node3Disp<T1,T2>*>,
                                      CompareNodePtr<T1>>& queue,
                                      std::vector<bool>& inQueue,
                                      std::vector<bool>& frozen,
                                      const size_t threadNo) const {
        
        while ( !queue.empty() ) {
            const Node3Disp<T1,T2>* src = queue.top();
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
    T1 Grid3Duisp<T1,T2>::getTraveltime(const sxyz<T1>& Rx,
                                        const std::vector<Node3Disp<T1,T2>>& nodes,
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
    T1 Grid3Duisp<T1,T2>::getTraveltime(const sxyz<T1>& Rx,
                                        const std::vector<Node3Disp<T1,T2>>& nodes,
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

//
//  Grid3Dun.h
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

#ifndef ttcr_Grid3Dun_h
#define ttcr_Grid3Dun_h

#include <cassert>

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifdef VTK
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkProbeFilter.h"
#include "vtkRectilinearGrid.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSmartPointer.h"
#include "vtkTetra.h"
#include "vtkXMLRectilinearGridWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"
#endif


#include "Grad.h"
#include "Grid3D.h"
#include "Interpolator.h"
#include "utils.h"

namespace ttcr {
    
    template<typename T1, typename T2, typename NODE>
    class Grid3Dun : public Grid3D<T1,T2> {
    public:
        Grid3Dun(const std::vector<sxyz<T1>>& no,
                 const std::vector<tetrahedronElem<T2>>& tet,
                 const int rp,
                 const bool iv,
                 const bool rptt,
                 const T1 md,
                 const size_t nt=1) :
        nThreads(nt), rp_method(rp), interpVel(iv), tt_from_rp(rptt),
        nPrimary(static_cast<T2>(no.size())),
        source_radius(0.0), min_dist(md),
        nodes(std::vector<NODE>(no.size(), NODE(nt))),
        neighbors(std::vector<std::vector<T2>>(tet.size())),
        tetrahedra(tet)
        {}
        
        virtual ~Grid3Dun() {}
        
        void setSlowness(const T1 s) {
            for ( size_t n=0; n<nodes.size(); ++n ) {
                nodes[n].setNodeSlowness( s[n] );
            }
        }
        
        void setSlowness(const T1 *s, const size_t ns) {
            if ( nodes.size() != ns ) {
                throw std::length_error("Error: slowness vectors of incompatible size.");
            }
            for ( size_t n=0; n<nodes.size(); ++n ) {
                nodes[n].setNodeSlowness( s[n] );
            }
        }
        
        void setSlowness(const std::vector<T1>& s) {
            if ( nodes.size() != s.size() ) {
                throw std::length_error("Error: slowness vectors of incompatible size.");
            }
            for ( size_t n=0; n<nodes.size(); ++n ) {
                nodes[n].setNodeSlowness( s[n] );
            }
        }
        
        void setSourceRadius(const double r) { source_radius = r; }
        
        void setTT(const T1 tt, const size_t nn, const size_t nt=0) {
            nodes[nn].setTT(tt, nt);
        }
        
        size_t getNumberOfNodes() const { return nodes.size(); }
        
        const T1 getXmin() const {
            T1 xmin = nodes[0].getX();
            for ( auto it=nodes.begin(); it!=nodes.end(); ++it )
                xmin = xmin<it->getX() ? xmin : it->getX();
            return xmin;
        }
        const T1 getXmax() const {
            T1 xmax = nodes[0].getX();
            for ( auto it=nodes.begin(); it!=nodes.end(); ++it )
                xmax = xmax>it->getX() ? xmax : it->getX();
            return xmax;
        }
        const T1 getYmin() const {
            T1 ymin = nodes[0].getY();
            for ( auto it=nodes.begin(); it!=nodes.end(); ++it )
                ymin = ymin<it->getY() ? ymin : it->getY();
            return ymin;
        }
        const T1 getYmax() const {
            T1 ymax = nodes[0].getY();
            for ( auto it=nodes.begin(); it!=nodes.end(); ++it )
                ymax = ymax>it->getY() ? ymax : it->getY();
            return ymax;
        }
        const T1 getZmin() const {
            T1 zmin = nodes[0].getZ();
            for ( auto it=nodes.begin(); it!=nodes.end(); ++it )
                zmin = zmin<it->getZ() ? zmin : it->getZ();
            return zmin;
        }
        const T1 getZmax() const {
            T1 zmax = nodes[0].getZ();
            for ( auto it=nodes.begin(); it!=nodes.end(); ++it )
                zmax = zmax>it->getZ() ? zmax : it->getZ();
            return zmax;
        }
        
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              const size_t threadNo=0) const {}
        
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                              std::vector<std::vector<T1>*>& traveltimes,
                              const size_t=0) const {}
        
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<sxyz<T1>>>& r_data,
                              const size_t threadNo=0) const {}
        
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                              std::vector<std::vector<T1>*>& traveltimes,
                              std::vector<std::vector<std::vector<sxyz<T1>>>*>& r_data,
                              const size_t=0) const {}
        
        void saveTT(const std::string &, const int, const size_t nt=0,
                    const bool vtkFormat=0) const;
        
#ifdef VTK
        void saveModelVTU(const std::string &, const bool saveSlowness=true,
                          const bool savePhysicalEntity=false) const;
#endif
        
        const size_t getNthreads() const { return nThreads; }
        
        void dump_secondary(std::ofstream& os) const {
            for ( size_t n=nPrimary; n<nodes.size(); ++n ) {
                os << nodes[n].getX() << ' ' << nodes[n].getY() << ' ' << nodes[n].getZ() << '\n';
            }
        }

    protected:
        const size_t nThreads;
        int rp_method;
        bool interpVel;
        bool tt_from_rp;
        T2 nPrimary;
        T1 source_radius;
        T1 min_dist;
        mutable std::vector<NODE> nodes;
        std::vector<std::vector<T2>> neighbors;  // nodes common to a cell
        std::vector<tetrahedronElem<T2>> tetrahedra;
        
        T1 computeDt(const NODE& source, const NODE& node) const {
            return (node.getNodeSlowness()+source.getNodeSlowness())/2 * source.getDistance( node );
        }
        
        T1 computeDt(const NODE& source, const sxyz<T1>& node, T1 slo) const {
            return (slo+source.getNodeSlowness())/2 * source.getDistance( node );
        }
        
        T1 getTraveltime(const sxyz<T1>& Rx,
                         const std::vector<NODE>& nodes,
                         const size_t threadNo) const;
        
        void checkPts(const std::vector<sxyz<T1>>&) const;
        
        bool insideTetrahedron(const sxyz<T1>&, const T2) const;
        bool insideTetrahedron_old(const sxyz<T1>&, const T2) const;
        
        T2 getCellNo(const sxyz<T1>& pt) const {
            for ( T2 n=0; n<tetrahedra.size(); ++n ) {
                if ( insideTetrahedron(pt, n) ) {
                    return n;
                }
            }
            return -1;
        }
        
        void buildGridNodes(const std::vector<sxyz<T1>>&, const size_t);
        void buildGridNodes(const std::vector<sxyz<T1>>&,
                            const int, const size_t, const int);
        
        void buildGridNeighbors() {
            // Index the neighbors nodes of each cell
            for ( T2 n=0; n<nodes.size(); ++n ) {
                for ( size_t n2=0; n2<nodes[n].getOwners().size(); ++n2) {
                    neighbors[ nodes[n].getOwners()[n2] ].push_back(n);
                }
            }
        }
        
        void localUpdate3D(NODE *vertexC, const size_t threadNo) const;
        
        T1 localUpdate2D(const NODE *vertexA,
                         const NODE *vertexB,
                         const NODE *vertexC,
                         const T2 tetraNo,
                         const size_t threadNo) const;
        
        void local3Dsolver(NODE *vertexC, const size_t threadNo) const;
        
        T1 local2Dsolver(const NODE *vertexA,
                         const NODE *vertexB,
                         const NODE *vertexC,
                         const T2 tetraNo,
                         const size_t threadNo) const;
        
        int solveEq23(const T1 a[], const T1 b[], T1 n[][3]) const;
        
        bool testInTriangle(const NODE *vertexA,
                            const NODE *vertexB,
                            const NODE *vertexC,
                            const sxyz<T1> &E) const;
        
        void barycentric(const NODE *a, const NODE *b, const NODE *c,
                         const sxyz<T1> &p, T1 &u, T1 &v, T1 &w) const;
        
        T1 getTraveltimeFromRaypath(const std::vector<sxyz<T1>>& Tx,
                                    const std::vector<T1>& t0,
                                    const sxyz<T1>& Rx,
                                    const size_t threadNo) const;
        
        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const sxyz<T1>& Rx,
                        std::vector<sxyz<T1>>& r_data,
                        const size_t threadNo) const;
        
        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const sxyz<T1>& Rx,
                        std::vector<sxyz<T1>>& r_data,
                        T1 &tt,
                        const size_t threadNo) const;
        
        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const sxyz<T1>& Rx,
                        std::vector<sxyz<T1>>& r_data,
                        std::vector<sijv<T1>>& m_data,
                        const size_t RxNo,
                        const size_t threadNo) const;
        
        void getRaypath_blti(const std::vector<sxyz<T1>>& Tx,
                             const std::vector<T1>& t0,
                             const sxyz<T1>& Rx,
                             std::vector<sxyz<T1>>& r_data,
                             T1 &tt,
                             const size_t threadNo) const;
        
        bool blti_solver_around_source(const sxyz<T1>& Source,
                                      const sxyz<T1>& curr_pt,
                                      const std::array<T2,3>& face,
                                      std::array<T1,3>& barycenters) const;
        bool blti_raytrace(const sxyz<T1>& curr_pt,
                           const std::array<T2, 3>& faces,
                           sxyz<T1>& next_pt,
                           const size_t threadNo,
                           const T1& s) const;
        bool blti2D_raytrace(const sxyz<T1>& curr_pt,
                             const T2& node1,
                             const T2& node2,
                             sxyz<T1>& next_pt,
                             const size_t threadNo,
                             const T1& s) const;

        
        bool check_pt_location(sxyz<T1>& curr_pt,
                               const std::array<T2,3>& ind,
                               bool& onNode,
                               T2& nodeNo,
                               bool& onEdge,
                               std::array<T2,2>& edgeNodes,
                               bool& onFace,
                               std::array<T2,3>& faceNodes) const;
        
        bool check_pt_location(sxyz<T1>& curr_pt,
                               const std::vector<T2>& ind1,
                               const std::array<T2,3>& ind2,
                               bool& onNode,
                               T2& nodeNo,
                               bool& onEdge,
                               std::array<T2,2>& edgeNodes,
                               bool& onFace,
                               std::array<T2,3>& faceNodes) const;
        
        void update_m_data(std::vector<sijv<T1>>& m_data,
                           sijv<T1>& m,
                           const std::set<T2>& allNodes,
                           const sxyz<T1>& mid_pt,
                           const T1 s,
                           const T1 ds) const ;
        
        bool intersectVecTriangle(const T2 iO, const sxyz<T1> &vec,
                                  const T2 iA, T2 iB, T2 iC,
                                  sxyz<T1> &pt_i) const;
        bool intersectVecTriangle(const sxyz<T1> &O, const sxyz<T1> &vec,
                                  const T2 iA, T2 iB, T2 iC,
                                  sxyz<T1> &pt_i) const;
        
        bool areCollinear(const sxyz<T1> &pt, const T2 i0, const T2 i1) const;
        bool areCoplanar(const sxyz<T1> &pt, const T2 i0, const T2 i1, const T2 i2) const;
        bool areCoplanar(const sxyz<T1> &pt1, const sxyz<T1> &pt2,
                         const sxyz<T1> &pt3, const sxyz<T1> &pt4) const;
        
        T2 findAdjacentCell1(const std::array<T2,3> &faceNodes, const T2 nodeNo) const;
        T2 findAdjacentCell2(const std::array<T2,3> &faceNodes, const T2 cellNo) const;
        T2 findAdjacentCell2(const std::array<T2,3> &faceNodes,
                             const T2 & cellNo,
                             const sxyz<T1>& curr_pt ) const;
        
        void getNeighborNodes(const T2, std::set<NODE*>&) const;
        void getNeighborNodesAB(const std::vector<NODE*>&,
                                std::vector<std::vector<std::array<NODE*,3>>>&) const;

        void plotCell(const T2 cellNo, const sxyz<T1> &pt, const sxyz<T1> &g) const;
        
        T1 computeSlowness( const sxyz<T1>& Rx ) const;
        T1 computeSlowness(const sxyz<T1>& curr_pt,
                           const bool onNode,
                           const T2 nodeNo,
                           const bool onEdge,
                           const std::array<T2,2>& edgeNodes,
                           const std::array<T2,3>& faceNodes) const;

        void printRaypathData(const sxyz<T1>& curr_pt,
                              const sxyz<T1>& g,
                              const bool onNode,
                              const bool onEdge,
                              const bool onFace,
                              const T2 cellNo,
                              const T2 nodeNo,
                              const std::array<T2,2> &edgeNodes,
                              const std::array<T2,3> &faceNodes) const;
        
        std::array<T2,4> getPrimary(const T2 cellNo) const {
            size_t i = 0;
            std::array<T2,4> tmp;
            for (size_t n=0; n<neighbors[cellNo].size(); ++n) {
                if ( nodes[neighbors[cellNo][n]].isPrimary() )
                    tmp[i++] = neighbors[cellNo][n];
                if ( i == 4 ) break;
            }
            return tmp;
        }
        
        void checkCloseToTx(const sxyz<T1>& curr_pt,
                       sxyz<T1>& g,
                       const T2 cellNo,
                       const std::vector<sxyz<T1>>& Tx,
                       const std::vector<T2>& txCell) const {
            
            for (size_t nt=0; nt<Tx.size(); ++nt) {
                std::array<T2,4> itmp = getPrimary(txCell[nt]);
                for (size_t n=0; n<4; ++n){
                    for (auto nc=nodes[itmp[n]].getOwners().begin(); nc!=nodes[itmp[n]].getOwners().end(); ++nc){
                        if ( *nc==cellNo ) {
                            g = Tx[nt]-curr_pt;
                            return;
                        }
                    }
                }
            }
        }
        
    };
    
    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::buildGridNodes(const std::vector<sxyz<T1>>& no,
                                              const size_t nt) {
        
        // primary nodes
        for ( T2 n=0; n<no.size(); ++n ) {
            nodes[n].setXYZindex( no[n].x, no[n].y, no[n].z, n );
            nodes[n].setPrimary(5);
        }
        
        //
        //              1
        //            ,/|`\
        //          ,/  |  `\
        //        ,0    '.   `4
        //      ,/       1     `\
        //    ,/         |       `\
        //   0-----5-----'.--------3
        //    `\.         |      ,/
        //       `\.      |     3
        //          `2.   '. ,/
        //             `\. |/
        //                `2
        //
        //
        //  triangle 0:  0-1  1-2  2-0     (first occurence of segment underlined)
        //               ---  ---  ---
        //  triangle 1:  1-2  2-3  3-1
        //                    ---  ---
        //  triangle 2:  0-2  2-3  3-0
        //                         ---
        //  triangle 3:  0-1  1-3  3-0
        
        
        for ( T2 ntet=0; ntet<tetrahedra.size(); ++ntet ) {
            
            // for each triangle
            for ( T2 ntri=0; ntri<4; ++ntri ) {
                
                // push owner for primary nodes
                nodes[ tetrahedra[ntet].i[ntri] ].pushOwner( ntet );
            }
        }
    }
    
    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::buildGridNodes(const std::vector<sxyz<T1>>& no,
                                              const int nsecondary,
                                              const size_t nt, const int verbose) {
        
        // primary nodes
        for ( T2 n=0; n<no.size(); ++n ) {
            nodes[n].setXYZindex( no[n].x, no[n].y, no[n].z, n );
            nodes[n].setPrimary(5);
        }
        T2 nNodes = static_cast<T2>(nodes.size());
        
        size_t nFaceNodes = 0;
        for ( int n=1; n<=(nsecondary-1); ++n ) nFaceNodes += n;
        
        std::map<std::array<T2,2>,std::vector<T2>> lineMap;
        std::array<T2,2> lineKey;
        typename std::map<std::array<T2,2>,std::vector<T2>>::iterator lineIt;
        
        size_t estLineNo = (tetrahedra.size()+tetrahedra.size()/10) * 6/2;
        size_t estFaceNo = (tetrahedra.size()+tetrahedra.size()/10) * 4/2;
        nodes.reserve( nNodes + estLineNo*nsecondary + estFaceNo*nFaceNodes );
        
        
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
        //    `\.         |      ,/
        //       `\.      |     3
        //          `2.   '. ,/
        //             `\. |/
        //                `2
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
        NODE tmpNode(nt);
        for ( T2 ntet=0; ntet<tetrahedra.size(); ++ntet ) {
            
            if ( verbose>1 && nsecondary > 0 ) {
                std::cout << "\r  Building edge nodes: " << (100*ntet)/tetrahedra.size() << "%";
                std::cout.flush();
            }
            
            // for each triangle
            for ( T2 ntri=0; ntri<4; ++ntri ) {
                
                // push owner for primary nodes
                nodes[ tetrahedra[ntet].i[ntri] ].pushOwner( ntet );
                
                if ( nsecondary > 0 ) {
                    // start from ntri to avoid redundancy
                    for ( size_t nl=ntri; nl<3; ++nl ) {
                        
                        lineKey = {tetrahedra[ntet].i[ iNodes[ntri][nl] ],
                            tetrahedra[ntet].i[ iNodes[ntri][(nl+1)%3] ]};
                        std::sort(lineKey.begin(), lineKey.end());
                        
                        lineIt = lineMap.find( lineKey );
                        if ( lineIt == lineMap.end() ) {
                            // not found, insert new pair
                            lineMap[ lineKey ] = std::vector<T2>(nsecondary);
                        } else {
                            for ( size_t n=0; n<lineIt->second.size(); ++n ) {
                                // setting owners
                                nodes[ lineIt->second[n] ].pushOwner( ntet );
                            }
                            continue;
                        }
                        
                        sxyz<T1> d = (no[lineKey[1]]-no[lineKey[0]])/static_cast<T1>(nsecondary+1);
                        
                        for ( size_t n2=0; n2<nsecondary; ++n2 ) {
                            tmpNode.setXYZindex(no[lineKey[0]].x+(1+n2)*d.x,
                                                no[lineKey[0]].y+(1+n2)*d.y,
                                                no[lineKey[0]].z+(1+n2)*d.z,
                                                nNodes );
                            lineMap[lineKey][n2] = nNodes++;
                            nodes.push_back( tmpNode );
                            nodes.back().pushOwner( ntet );
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
            
            for ( T2 ntet=0; ntet<tetrahedra.size(); ++ntet ) {
                
                if ( verbose>1 ) {
                    std::cout << "\r  Building face nodes: " << (100*ntet)/tetrahedra.size() << "%";
                    std::cout.flush();
                }
                
                // for each triangle
                for ( T2 ntri=0; ntri<4; ++ntri ) {
                    
                    faceKey = {tetrahedra[ntet].i[ iNodes[ntri][0] ],
                        tetrahedra[ntet].i[ iNodes[ntri][1] ],
                        tetrahedra[ntet].i[ iNodes[ntri][2] ]};
                    std::sort(faceKey.begin(), faceKey.end());
                    
                    faceIt = faceMap.find( faceKey );
                    if ( faceIt == faceMap.end() ) {
                        // not found, insert new pair
                        faceMap[ faceKey ] = std::vector<T2>(nFaceNodes);
                    } else {
                        for ( size_t n=0; n<faceIt->second.size(); ++n ) {
                            // setting owners
                            nodes[ faceIt->second[n] ].pushOwner( ntet );
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
                            faceMap[faceKey][ifn++] = nNodes++;
                            nodes.push_back( tmpNode );
                            nodes.back().pushOwner( ntet );
                        }
                    }
                }
            }
        }
        if ( verbose>1 && nsecondary > 0 ) {
            std::cout << "\r  Building face nodes: 100%\n";
        }
        
        nodes.shrink_to_fit();
    }
    
    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dun<T1,T2,NODE>::getTraveltime(const sxyz<T1>& Rx,
                                           const std::vector<NODE>& nodes,
                                           const size_t threadNo) const {
        
        for ( size_t nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == Rx ) {
                return nodes[nn].getTT(threadNo);
            }
        }
        //If Rx is not on a node:
        T1 slo = computeSlowness( Rx );
        
        T2 cellNo = getCellNo( Rx );
        
        T2 neibNo = neighbors[cellNo][0];
        T1 dt = computeDt(nodes[neibNo], Rx, slo);
        
        T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
        for ( size_t k=1; k< neighbors[cellNo].size(); ++k ) {
            neibNo = neighbors[cellNo][k];
            dt = computeDt(nodes[neibNo], Rx, slo);
            if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
                traveltime =  nodes[neibNo].getTT(threadNo)+dt;
            }
        }
        return traveltime;
    }
    
    
    
    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::checkPts(const std::vector<sxyz<T1>>& pts) const {
        
        for (size_t n=0; n<pts.size(); ++n) {
            bool found = false;
            // check first if point is on a node
            for ( T2 nt=0; nt<nodes.size(); ++nt ) {
                if ( nodes[nt] == pts[n]) {
                    found = true;
                    break;
                }
            }
            if ( found == false ) {
                // check if inside tetrahedra
                for ( T2 nt=0; nt<tetrahedra.size(); ++nt ) {
                    if ( insideTetrahedron(pts[n], nt) ) {
                        found = true;
                        break;
                    }
                }
            }
            if ( found == false ) {
                std::ostringstream msg;
                msg << "Error: Point (" << pts[n].x << ", " << pts[n].y << ", " << pts[n] .z << ") outside grid.";
                throw std::runtime_error(msg.str());
            }
        }
    }
    
    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::insideTetrahedron_old(const sxyz<T1>& p, const T2 nt) const {
        
        sxyz<T1> v1 = { nodes[ tetrahedra[nt].i[0] ].getX(),
            nodes[ tetrahedra[nt].i[0] ].getY(),
            nodes[ tetrahedra[nt].i[0] ].getZ() };
        
        sxyz<T1> v2 = { nodes[ tetrahedra[nt].i[1] ].getX(),
            nodes[ tetrahedra[nt].i[1] ].getY(),
            nodes[ tetrahedra[nt].i[1] ].getZ() };
        
        sxyz<T1> v3 = { nodes[ tetrahedra[nt].i[2] ].getX(),
            nodes[ tetrahedra[nt].i[2] ].getY(),
            nodes[ tetrahedra[nt].i[2] ].getZ() };
        
        sxyz<T1> v4 = { nodes[ tetrahedra[nt].i[3] ].getX(),
            nodes[ tetrahedra[nt].i[3] ].getY(),
            nodes[ tetrahedra[nt].i[3] ].getZ() };
        
        return sameSide(v1, v2, v3, v4, p) && sameSide(v2, v3, v4, v1, p) &&
            sameSide(v3, v4, v1, v2, p) && sameSide(v4, v1, v2, v3, p);
    }
    
    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::insideTetrahedron(const sxyz<T1>& v, const T2 nt) const {
        
        
        // from http://steve.hollasch.net/cgindex/geometry/ptintet.html
        
        sxyz<T1> v1 = { nodes[ tetrahedra[nt].i[0] ].getX(),
            nodes[ tetrahedra[nt].i[0] ].getY(),
            nodes[ tetrahedra[nt].i[0] ].getZ() };
        
        sxyz<T1> v2 = { nodes[ tetrahedra[nt].i[1] ].getX(),
            nodes[ tetrahedra[nt].i[1] ].getY(),
            nodes[ tetrahedra[nt].i[1] ].getZ() };
        
        sxyz<T1> v3 = { nodes[ tetrahedra[nt].i[2] ].getX(),
            nodes[ tetrahedra[nt].i[2] ].getY(),
            nodes[ tetrahedra[nt].i[2] ].getZ() };
        
        sxyz<T1> v4 = { nodes[ tetrahedra[nt].i[3] ].getX(),
            nodes[ tetrahedra[nt].i[3] ].getY(),
            nodes[ tetrahedra[nt].i[3] ].getZ() };
        
        T1 D0 = det4(v1, v2, v3, v4);
        T1 D1 = det4( v, v2, v3, v4);
        T1 D2 = det4(v1,  v, v3, v4);
        T1 D3 = det4(v1, v2,  v, v4);
        T1 D4 = det4(v1, v2, v3,  v);
        
        int t1, t2, t3, t4;
        
        if ( fabs(D1)<small2 ) {
            // points are coplanar, check if pt is inside triangle
            if ( testInTriangle(&(nodes[ tetrahedra[nt].i[1] ]),
                                &(nodes[ tetrahedra[nt].i[2] ]),
                                &(nodes[ tetrahedra[nt].i[3] ]), v))
                return 1;
            else
                return 0;
            
        } else {
            t1 = (signum(D0)==signum(D1));
        }
        
        if ( fabs(D2)<small2 ) {
            if ( testInTriangle(&(nodes[ tetrahedra[nt].i[0] ]),
                                &(nodes[ tetrahedra[nt].i[2] ]),
                                &(nodes[ tetrahedra[nt].i[3] ]), v))
                return 1;
            else
                return 0;
        } else {
            t2 = (signum(D0)==signum(D2));
        }
        
        if ( fabs(D3)<small2 ) {
            if ( testInTriangle(&(nodes[ tetrahedra[nt].i[0] ]),
                                &(nodes[ tetrahedra[nt].i[1] ]),
                                &(nodes[ tetrahedra[nt].i[3] ]), v))
                return 1;
            else
                return 0;
        } else {
            t3 = (signum(D0)==signum(D3));
        }
        
        if ( fabs(D4)<small2 ) {
            if ( testInTriangle(&(nodes[ tetrahedra[nt].i[0] ]),
                                &(nodes[ tetrahedra[nt].i[1] ]),
                                &(nodes[ tetrahedra[nt].i[2] ]), v))
                return 1;
            else
                return 0;
        } else {
            t4 = (signum(D0)==signum(D4));
        }
        
        return t1 && t2 && t3 && t4;
    }
    
    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::saveTT(const std::string &fname, const int all,
                                      const size_t nt, const bool vtkFormat) const {
        
        if (vtkFormat) {
#ifdef VTK
            std::string filename = fname+".vtu";
            
            vtkSmartPointer<vtkUnstructuredGrid> ugrid =
            vtkSmartPointer<vtkUnstructuredGrid>::New();
            
            vtkSmartPointer<vtkPoints> newPts =
            vtkSmartPointer<vtkPoints>::New();
            vtkSmartPointer<vtkDoubleArray> newScalars =
            vtkSmartPointer<vtkDoubleArray>::New();
            
            newScalars->SetName("Travel time");
            
            double xyz[3];
            T2 nMax = nPrimary;  // only primary are saved
            for (size_t n=0; n<nMax; ++n) {
                xyz[0] = nodes[n].getX();
                xyz[1] = nodes[n].getY();
                xyz[2] = nodes[n].getZ();
                newPts->InsertPoint(n, xyz);
                newScalars->InsertValue(n, nodes[n].getTT(nt) );
            }
            
            ugrid->SetPoints(newPts);
            ugrid->GetPointData()->SetScalars(newScalars);
            
            vtkSmartPointer<vtkTetra> tet =
            vtkSmartPointer<vtkTetra>::New();
            for (size_t n=0; n<tetrahedra.size(); ++n) {
                tet->GetPointIds()->SetId(0, tetrahedra[n].i[0] );
                tet->GetPointIds()->SetId(1, tetrahedra[n].i[1] );
                tet->GetPointIds()->SetId(2, tetrahedra[n].i[2] );
                tet->GetPointIds()->SetId(3, tetrahedra[n].i[3] );
                
                ugrid->InsertNextCell( tet->GetCellType(), tet->GetPointIds() );
            }
            vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
            vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
            
            writer->SetFileName( filename.c_str() );
            //		writer->SetInputConnection( ugrid->GetProducerPort() );
            writer->SetInputData( ugrid );
            writer->SetDataModeToBinary();
            writer->Update();
#else
            std::cerr << "VTK not included during compilation.\nNothing saved.\n";
#endif
        } else {
            std::string filename = fname+".dat";
            std::ofstream fout(filename.c_str());
            T2 nMax = nPrimary;
            if ( all == 1 ) {
                nMax = static_cast<T2>(nodes.size());
            }
            for ( T2 n=0; n<nMax; ++n ) {
                fout << nodes[n].getX() << '\t'
                << nodes[n].getY() << '\t'
                << nodes[n].getZ() << '\t'
                << nodes[n].getTT(nt) << '\n';
            }
            fout.close();
        }
    }
    
#ifdef VTK
    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::saveModelVTU(const std::string &fname,
                                            const bool saveSlowness,
                                            const bool savePhysicalEntity) const {
        
        vtkSmartPointer<vtkUnstructuredGrid> ugrid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
        
        vtkSmartPointer<vtkPoints> newPts =
        vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkDoubleArray> newScalars =
        vtkSmartPointer<vtkDoubleArray>::New();
        
        double xyz[3];
        T2 nMax = nPrimary;  // only primary are saved
        if ( saveSlowness ) {
            newScalars->SetName("Slowness");
            
            for (size_t n=0; n<nMax; ++n) {
                xyz[0] = nodes[n].getX();
                xyz[1] = nodes[n].getY();
                xyz[2] = nodes[n].getZ();
                newPts->InsertPoint(n, xyz);
                newScalars->InsertValue(n, nodes[n].getNodeSlowness() );
            }
        } else {
            newScalars->SetName("Velocity");
            
            for (size_t n=0; n<nMax; ++n) {
                xyz[0] = nodes[n].getX();
                xyz[1] = nodes[n].getY();
                xyz[2] = nodes[n].getZ();
                newPts->InsertPoint(n, xyz);
                newScalars->InsertValue(n, static_cast<T1>(1.0)/nodes[n].getNodeSlowness() );
            }
        }
        
        ugrid->SetPoints(newPts);
        ugrid->GetPointData()->SetScalars(newScalars);
        
        vtkSmartPointer<vtkTetra> tet =
        vtkSmartPointer<vtkTetra>::New();
        
        for (size_t n=0; n<tetrahedra.size(); ++n) {
            tet->GetPointIds()->SetId(0, tetrahedra[n].i[0] );
            tet->GetPointIds()->SetId(1, tetrahedra[n].i[1] );
            tet->GetPointIds()->SetId(2, tetrahedra[n].i[2] );
            tet->GetPointIds()->SetId(3, tetrahedra[n].i[3] );
            
            ugrid->InsertNextCell( tet->GetCellType(), tet->GetPointIds() );
        }
        
        vtkSmartPointer<vtkIntArray> data_pe = vtkSmartPointer<vtkIntArray>::New();
        if ( savePhysicalEntity ) {
            data_pe->SetName("Physical entity");
            for (size_t n=0; n<tetrahedra.size(); ++n) {
                data_pe->InsertNextValue(tetrahedra[n].physical_entity );
            }
            ugrid->GetCellData()->AddArray(data_pe);
        }
        
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        
        writer->SetFileName( fname.c_str() );
        writer->SetInputData( ugrid );
        writer->SetDataModeToBinary();
        writer->Update();
        
    }
    
#endif
    
    
    
    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::localUpdate3D(NODE *vertexD,
                                             const size_t threadNo) const {
        
        // method of Lelievre et al. 2011
        
        T2 iA, iB, iC, iD;
        NODE *vertexA, *vertexB, *vertexC;
        
        for ( size_t no=0; no<vertexD->getOwners().size(); ++no ) {
            
            T2 tetNo = vertexD->getOwners()[no];
            
            for ( iD=0; iD<4; ++iD ) {
                if ( vertexD->getGridIndex() == tetrahedra[tetNo].i[iD] ) break;
            }
            
            iA = (iD+1)%4;
            iB = (iD+2)%4;
            iC = (iD+3)%4;
            vertexA = &(nodes[tetrahedra[tetNo].i[iA]]);
            vertexB = &(nodes[tetrahedra[tetNo].i[iB]]);
            vertexC = &(nodes[tetrahedra[tetNo].i[iC]]);
            
            if ( vertexA->getTT(threadNo) > vertexB->getTT(threadNo) ) {
                std::swap(iA, iB);
                std::swap(vertexA, vertexB);
            }
            if ( vertexA->getTT(threadNo) > vertexC->getTT(threadNo) ) {
                std::swap(iA, iC);
                std::swap(vertexA, vertexC);
            }
            if ( vertexB->getTT(threadNo) > vertexC->getTT(threadNo) ) {
                std::swap(iB, iC);
                std::swap(vertexB, vertexC);
            }

            if ( vertexA->getTT(threadNo) == std::numeric_limits<T1>::max() ) {
                continue;
            }
            
            T1 tABC = std::numeric_limits<T1>::max();
            
            if ( vertexB->getTT(threadNo) != std::numeric_limits<T1>::max() &&
                vertexC->getTT(threadNo) != std::numeric_limits<T1>::max() ) {
                
                T1 u = vertexB->getTT(threadNo) - vertexA->getTT(threadNo);
                T1 v = vertexC->getTT(threadNo) - vertexA->getTT(threadNo);
                
                sxyz<T1> v_b = { vertexC->getX() - vertexA->getX(),
                    vertexC->getY() - vertexA->getY(),
                    vertexC->getZ() - vertexA->getZ() };
                sxyz<T1> v_c = { vertexB->getX() - vertexA->getX(),
                    vertexB->getY() - vertexA->getY(),
                    vertexB->getZ() - vertexA->getZ() };
                
                // vector normal to plane
                sxyz<T1> v_n = cross(v_b, v_c);
                
                T1 b = norm( v_b );
                T1 c = norm( v_c );
                T1 d2 = dot(v_b, v_c);
                
                T1 alpha = acos( d2 / (b*c) );
                
                T1 phi = c*b*sin(alpha);
                
                // check for negative value
                T1 w_tilde = vertexD->getNodeSlowness()*vertexD->getNodeSlowness()*phi*phi -
                                  u*u*b*b - v*v*c*c + 2.*u*v*d2;
                if ( w_tilde > 0.0 ) {
                
                    w_tilde = sqrt( w_tilde );
                    
                    // Point (ξ_0 , ζ_0 ) is the normalized projection of node D onto face ABC
                    // project D on plane
                    
                    T1 d_tmp = -vertexA->getX()*v_n.x - vertexA->getY()*v_n.y - vertexA->getZ()*v_n.z;
                    
                    T1 k = -(d_tmp + v_n.x*vertexD->getX() + v_n.y*vertexD->getY() + v_n.z*vertexD->getZ())/
                    norm2(v_n);
                    
                    sxyz<T1> pt;   // -> Point (ξ_0 , ζ_0 )
                    pt.x = vertexD->getX() + k*v_n.x;
                    pt.y = vertexD->getY() + k*v_n.y;
                    pt.z = vertexD->getZ() + k*v_n.z;
                    
                    T1 rho0 = vertexD->getDistance( pt );
                    
                    sxyz<T1> v_pt = {pt.x-vertexA->getX(), pt.y-vertexA->getY(), pt.z-vertexA->getZ()};
                    
                    T1 xi0;
                    T1 zeta0;
                    projNorm(v_b/b, v_c/c, v_pt, xi0, zeta0);
                    if ( xi0 < 0.0 || zeta0 < 0.0 ) {
                        // this should not happen unless we have incorrect triangle
                        continue;
                    }
                    
                    T1 beta = u*b*b - v*d2;
                    T1 gamma = v*c*c - u*d2;
                    
                    T1 xi_tilde = -fabs(beta)*rho0/(phi*w_tilde);
                    T1 zeta_tilde = -fabs(gamma)*rho0/(phi*w_tilde);
                    
                    T1 xi = xi_tilde + xi0;
                    T1 zeta = zeta_tilde + zeta0;
                    
                    if ( 0.<xi && xi<1. && 0.<zeta && zeta<1. && 0.<(xi+zeta) && (xi+zeta)<1. ) {
                        tABC = vertexA->getTT(threadNo) + u*xi0 + v*zeta0 + w_tilde*rho0/phi;
                    }
                }
            }
            
            T1 t = vertexA->getTT(threadNo) + vertexD->getNodeSlowness() * vertexD->getDistance( *vertexA );
            if ( t < tABC ) tABC = t;
            t = vertexB->getTT(threadNo) + vertexD->getNodeSlowness() * vertexD->getDistance( *vertexB );
            if ( t < tABC ) tABC = t;
            t = vertexC->getTT(threadNo) + vertexD->getNodeSlowness() * vertexD->getDistance( *vertexC );
            if ( t < tABC ) tABC = t;
            
            t = localUpdate2D(vertexA, vertexB, vertexD, tetNo, threadNo);
            if ( t < tABC ) tABC = t;
            t = localUpdate2D(vertexA, vertexC, vertexD, tetNo, threadNo);
            if ( t < tABC ) tABC = t;
            t = localUpdate2D(vertexB, vertexC, vertexD, tetNo, threadNo);
            if ( t < tABC ) tABC = t;
            
            if ( tABC<vertexD->getTT(threadNo) )
                vertexD->setTT(tABC, threadNo);
            
        }
    }
    
    
    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dun<T1,T2,NODE>::localUpdate2D(const NODE *vertexA,
                                           const NODE *vertexB,
                                           const NODE *vertexC,
                                           const T2 tetNo,
                                           const size_t threadNo) const {
        
        if ( vertexB->getTT(threadNo)==std::numeric_limits<T1>::max() &&
            vertexA->getTT(threadNo)==std::numeric_limits<T1>::max() ) {
            return std::numeric_limits<T1>::max();
        }
        T1 t;
        
        T1 u = vertexB->getTT(threadNo) - vertexA->getTT(threadNo);
        
        sxyz<T1> v_b = { vertexC->getX() - vertexA->getX(),
            vertexC->getY() - vertexA->getY(),
            vertexC->getZ() - vertexA->getZ() };
        sxyz<T1> v_c = { vertexB->getX() - vertexA->getX(),
            vertexB->getY() - vertexA->getY(),
            vertexB->getZ() - vertexA->getZ() };
        
        T1 c = norm( v_c );
        
        T1 w2 = vertexC->getNodeSlowness()*vertexC->getNodeSlowness()*c*c - u*u;
        if ( w2 < 0.0 ) return std::numeric_limits<T1>::max();
        
        T1 w = sqrt( w2 );
        
        T1 k = dot(v_b,v_c)/dot(v_c,v_c);
        sxyz<T1> pt;
        pt.x = vertexA->getX() + k*v_c.x;
        pt.y = vertexA->getY() + k*v_c.y;
        pt.z = vertexA->getZ() + k*v_c.z;
        
        T1 rho0 = vertexC->getDistance( pt );
//        T1 xi0 = vertexA->getDistance( pt )/c;
        T1 xi0 = k;
        
        T1 xi = xi0 - u*rho0/(w*c);
        
        if ( 0.<xi && xi<1. ) {
            t = vertexA->getTT(threadNo) + u*xi0 + w*rho0/c;
        } else {
            t = std::numeric_limits<T1>::max();
        }
        
        return t;
    }
    
    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::local3Dsolver(NODE *vertexD,
                                             const size_t threadNo) const {
        
        // Méthode de Qian et al. 2007
        
        T2 iA, iB, iC, iD;
        NODE *vertexA, *vertexB, *vertexC;
        T1 AB, AC;
        
        for ( size_t no=0; no<vertexD->getOwners().size(); ++no ) {
            
            T2 tetNo = vertexD->getOwners()[no];
            
            for ( iD=0; iD<4; ++iD ) {
                if ( vertexD->getGridIndex() == tetrahedra[tetNo].i[iD] ) break;
            }
            
            iA = (iD+1)%4;
            iB = (iD+2)%4;
            iC = (iD+3)%4;
            vertexA = &(nodes[tetrahedra[tetNo].i[iA]]);
            vertexB = &(nodes[tetrahedra[tetNo].i[iB]]);
            vertexC = &(nodes[tetrahedra[tetNo].i[iC]]);
            
            if ( vertexA->getTT(threadNo) > vertexB->getTT(threadNo) ) {
                std::swap(iA, iB);
                std::swap(vertexA, vertexB);
            }
            if ( vertexA->getTT(threadNo) > vertexC->getTT(threadNo) ) {
                std::swap(iA, iC);
                std::swap(vertexA, vertexC);
            }
            
            if ( vertexA->getTT(threadNo) == std::numeric_limits<T1>::max() ) {
                continue;
            }
            
            AB = vertexA->getDistance( *vertexB );
            AC = vertexA->getDistance( *vertexC );
            
            bool apply2Dsolvers = true;
            
            if (fabs(vertexB->getTT(threadNo)-vertexA->getTT(threadNo))<=AB*vertexD->getNodeSlowness() &&
                fabs(vertexC->getTT(threadNo)-vertexA->getTT(threadNo))<=AC*vertexD->getNodeSlowness()) {
                
                // Qian et al, 2007, eq 2.3
                
                T1 ab[4], ac[4], n[2][3];
                
                // vec(AB)
                ab[0] = vertexB->getX()-vertexA->getX();
                ab[1] = vertexB->getY()-vertexA->getY();
                ab[2] = vertexB->getZ()-vertexA->getZ();
                
                ab[3] = (vertexB->getTT(threadNo)-vertexA->getTT(threadNo)) / vertexD->getNodeSlowness();
                
                // vec(AC)
                ac[0] = vertexC->getX()-vertexA->getX();
                ac[1] = vertexC->getY()-vertexA->getY();
                ac[2] = vertexC->getZ()-vertexA->getZ();
                
                ac[3] = (vertexC->getTT(threadNo)-vertexA->getTT(threadNo)) / vertexD->getNodeSlowness();
                
                int rv = solveEq23(ab, ac, n);
                
                if ( rv == 1 ) {
                    
                    for ( size_t ns=0; ns<2; ++ns ) {
                        //
                        // find pt E
                        //
                        
                        // plane vec(AB) cross vec(AC) passing by A: ax + by + cz + d = 0
                        
                        T1 a = ab[1]*ac[2] - ac[1]*ab[2];
                        T1 b = ac[0]*ab[2] - ab[0]*ac[2];
                        T1 c = ab[0]*ac[1] - ac[0]*ab[1];
                        
                        T1 d = -vertexA->getX()*a - vertexA->getY()*b - vertexA->getZ()*c;
                        
                        //				T1 k = -(d + a*vertexA->getX() + b*vertexA->getY() + c*vertexA->getZ())/
                        //				(a*n[0] + b*n[1] + c*n[2]);
                        T1 k = -(d + a*vertexD->getX() + b*vertexD->getY() + c*vertexD->getZ())/  // TODO check here if vertex D
                        (a*n[ns][0] + b*n[ns][1] + c*n[ns][2]);
                        
                        sxyz<T1> E;
                        E.x = vertexD->getX() + k*n[ns][0];
                        E.y = vertexD->getY() + k*n[ns][1];
                        E.z = vertexD->getZ() + k*n[ns][2];
                        
                        if ( testInTriangle(vertexA, vertexB, vertexC, E) ) {
                            
                            // find point on wavefront plane
                            
                            a = n[ns][0];
                            b = n[ns][1];
                            c = n[ns][2];
                            d = -vertexA->getX()*a - vertexA->getY()*b - vertexA->getZ()*c;
                            
                            k = -(d + a*vertexD->getX() + b*vertexD->getY() + c*vertexD->getZ())/
                            (a*n[ns][0] + b*n[ns][1] + c*n[ns][2]);
                            
                            sxyz<T1> pt;
                            pt.x = vertexD->getX() + k*n[ns][0];
                            pt.y = vertexD->getY() + k*n[ns][1];
                            pt.z = vertexD->getZ() + k*n[ns][2];
                            
                            sxyz<T1> AD;
                            AD.x = vertexD->getX() - vertexA->getX();
                            AD.y = vertexD->getY() - vertexA->getY();
                            AD.z = vertexD->getZ() - vertexA->getZ();
                            
                            T1 d2 = vertexD->getDistance( E );
                            T1 d3 = vertexD->getDistance( pt );
                            T1 d4 = fabs( AD.x*n[ns][0] + AD.y*n[ns][1] + AD.z*n[ns][2] );
                            
                            if ( fabs(d3-d4)>small ) {
                                std::cout << " d3 ne d4: " << d3 << '\t' << d4 << '\t' << d2 << '\n';
                            }
                            
                            T1 t = vertexA->getTT(threadNo) +
                            d3*vertexD->getNodeSlowness();
                            
                            if ( t<vertexD->getTT(threadNo) )
                                vertexD->setTT(t, threadNo);
                            
                            apply2Dsolvers = false;
                            break;
                        }
                    }
                }
            }
            
            if ( apply2Dsolvers ) {
                T1 tABD = local2Dsolver(vertexA, vertexB, vertexD, tetNo, threadNo);
                T1 tACD = local2Dsolver(vertexA, vertexC, vertexD, tetNo, threadNo);
                T1 tBCD = local2Dsolver(vertexB, vertexC, vertexD, tetNo, threadNo);
                
                T1 t = tABD < tACD ? tABD : tACD;
                t = t < tBCD ? t : tBCD;
                
                if ( t<vertexD->getTT(threadNo) )
                    vertexD->setTT(t, threadNo);
            }
            
        }
    }
    
    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dun<T1,T2,NODE>::local2Dsolver(const NODE *vertexA,
                                           const NODE *vertexB,
                                           const NODE *vertexC,
                                           const T2 tetraNo,
                                           const size_t threadNo) const {
        static const double pi2 = pi / 2.;
        
        if ( vertexB->getTT(threadNo)==std::numeric_limits<T1>::max() &&
            vertexA->getTT(threadNo)==std::numeric_limits<T1>::max() ) {
            return std::numeric_limits<T1>::max();
        }
        
        T1 t;
        
        T1 a = vertexB->getDistance( *vertexC );
        T1 b = vertexA->getDistance( *vertexC );
        T1 c = vertexA->getDistance( *vertexB );
        if ( fabs(vertexB->getTT(threadNo)-vertexA->getTT(threadNo))<= c*vertexC->getNodeSlowness() ) {
            
            T1 theta = asin( fabs(vertexB->getTT(threadNo)-vertexA->getTT(threadNo))/
                            (c*vertexC->getNodeSlowness()) );
            
            T1 gamma = acos((a*a + b*b - c*c)/(2.*a*b));
            
            if ( gamma > pi2 ) {
                std::cout << "*** Obtuse angle: " << gamma*57.2957795 << " ***\n";
            } else {
                std::cout << "Accute angle: " << gamma*57.2957795 << " \n";
            }
            
            T1 beta  = acos((b*b + c*c - a*a)/(2.*b*c));
            T1 alpha = acos((a*a + c*c - b*b)/(2.*a*c));
            
            if ( ((0.>alpha-pi2?0.:alpha-pi2)<=theta && theta<=(pi2-beta) ) ||
                ((alpha-pi2)<=theta && theta<=(0.<pi2-beta?0.:pi2-beta)) ) {
                T1 h = a*sin(alpha-theta);
                T1 H = b*sin(beta+theta);
                t = 0.5*(h*vertexC->getNodeSlowness() + vertexB->getTT(threadNo)) +
                0.5 *(H*vertexC->getNodeSlowness() + vertexA->getTT(threadNo));
                
            } else {
                t = vertexA->getTT(threadNo) + b*vertexC->getNodeSlowness();
                t = t < vertexB->getTT(threadNo) + a*vertexC->getNodeSlowness() ? t : vertexB->getTT(threadNo) + a*vertexC->getNodeSlowness();
            }
        } else {
            t = vertexA->getTT(threadNo) + b*vertexC->getNodeSlowness();
            t = t < vertexB->getTT(threadNo) + a*vertexC->getNodeSlowness() ? t : vertexB->getTT(threadNo) + a*vertexC->getNodeSlowness();
        }
        t = t<vertexC->getTT(threadNo) ? t : vertexC->getTT(threadNo);
        
        return t;
    }
    
    
    template<typename T1, typename T2, typename NODE>
    int Grid3Dun<T1,T2,NODE>::solveEq23(const T1 a[], const T1 b[], T1 n[][3]) const {
        // returns 0 if no solution
        //         1 if solutions exist
        
        T1 a02 = a[0]*a[0];
        T1 a12 = a[1]*a[1];
        T1 a22 = a[2]*a[2];
        T1 a32 = a[3]*a[3];
        T1 b02 = b[0]*b[0];
        T1 b12 = b[1]*b[1];
        T1 b22 = b[2]*b[2];
        T1 b32 = b[3]*b[3];
        T1 a23 = a[2]*a[2]*a[2];
        
        T1 s1 = (a[2]*b[1] - a[1]*b[2])*(a[2]*b[1] - a[1]*b[2])*
        (a02*(b12 + b22) - a32*(b02 + b12 + b22) + 2*a[0]*a[3]*b[0]*b[3] - a02*b32 -
         2*a[1]*b[1]*(a[0]*b[0] + a[2]*b[2] - a[3]*b[3]) + 2*a[2]*b[2]*
         (-(a[0]*b[0]) + a[3]*b[3]) + a22*(b02 + b12 - b32) + a12*(b02 + b22 - b32));
        
        if ( s1 < 0.0 ) {
            return 0;
        } else {
            
            T1 d1 = (a22*(b02 + b12) - 2*a[0]*a[2]*b[0]*b[2] -
                     2*a[1]*b[1]*(a[0]*b[0] + a[2]*b[2])  +
                     a12*(b02 + b22) + a02*(b12 + b22));
            T1 d2 = ((a[2]*b[1] - a[1]*b[2]) *(a22*(b02 + b12) -
                                               2*a[0]*a[2]*b[0]*b[2] -
                                               2*a[1]*b[1]*(a[0]*b[0] + a[2]*b[2]) +
                                               a12*(b02 + b22) + a02*(b12 + b22)));
            
            if ( d1==0.0 || d2==0.0 ) return 0;
            
            s1 = sqrt(s1);
            
            n[0][0] = (a[0]*a[3]*(b12 + b22) + a12*b[0]*b[3] - a[0]*a[2]*b[2]*b[3] -
                       a[1]*b[1]*(a[3]*b[0] + a[0]*b[3]) +
                       a[2]*b[0]*(-(a[3]*b[2]) + a[2]*b[3])  - s1) / d1;
            
            n[0][1] = (a[2]*a[3]*b[1]*(-(a[0]*b[0]*b[1])  + a[1]*(b02 + 2*b22)) +
                       a23*b12*b[3] + a[2]*(a[0]*b[1]*(-(a[1]*b[0])  + a[0]*b[1]) +
                                            a12*b22)*b[3] -
                       a22*b[1]*b[2]*(a[3]*b[1] + 2*a[1]*b[3]) -
                       a[1]*b[2]*(a[1]*a[3]*(b02 + b22) - a[0]*a[1]*b[0]*b[3] +
                                  a[0]*b[1]*(-(a[3]*b[0])  + a[0]*b[3]) ) +
                       a[2]*b[0]*s1 - a[0]*b[2]*s1) / d2;
            
            n[0][2] = (a[1]*b22*(a[3]*(a[0]*b[0] + a[1]*b[1])  - (a02 + a12)*b[3]) +
                       a22*b[1]*(a[3]*(b02 + b12) - (a[0]*b[0] + a[1]*b[1]) *b[3]) +
                       a[2]*b[2]*(-(a[1]*a[3]*(b02 + 2*b12)) + a[0]*a[1]*b[0]*b[3] +
                                  2*a12*b[1]*b[3] + a[0]*b[1]*(-(a[3]*b[0]) + a[0]*b[3]) ) -
                       a[1]*b[0]*s1 + a[0]*b[1]*s1)/ d2;
            
            n[1][0] = (a[0]*a[3]*(b12 + b22) + a12*b[0]*b[3] - a[0]*a[2]*b[2]*b[3] -
                       a[1]*b[1]*(a[3]*b[0] + a[0]*b[3]) +
                       a[2]*b[0]*(-(a[3]*b[2]) + a[2]*b[3])  + s1) / d1;
            
            n[1][1] = (a[2]*a[3]*b[1]*(-(a[0]*b[0]*b[1])  + a[1]*(b02 + 2*b22)) +
                       a23*b12*b[3] + a[2]*(a[0]*b[1]*(-(a[1]*b[0])  + a[0]*b[1]) +
                                            a12*b22)*b[3] - a22*b[1]*b[2]*(a[3]*b[1] +
                                                                           2*a[1]*b[3]) -
                       a[1]*b[2]*(a[1]*a[3]*(b02 + b22) - a[0]*a[1]*b[0]*b[3] +
                                  a[0]*b[1]*(-(a[3]*b[0])  + a[0]*b[3]) ) -
                       a[2]*b[0]*s1 + a[0]*b[2]*s1) / d2;
            
            n[1][2] = (a[1]*b22*(a[3]*(a[0]*b[0] + a[1]*b[1])  - (a02 + a12)*b[3]) +
                       a22*b[1]*(a[3]*(b02 + b12) - (a[0]*b[0] + a[1]*b[1]) *b[3]) +
                       a[2]*b[2]*(-(a[1]*a[3]*(b02 + 2*b12)) + a[0]*a[1]*b[0]*b[3] +
                                  2*a12*b[1]*b[3] + a[0]*b[1]*(-(a[3]*b[0])  + a[0]*b[3]) ) +
                       a[1]*b[0]*s1 - a[0]*b[1]*s1) / d2;
        }
        return 1;
    }

    
    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dun<T1,T2,NODE>::getTraveltimeFromRaypath(const std::vector<sxyz<T1>>& Tx,
                                                      const std::vector<T1>& t0,
                                                      const sxyz<T1> &Rx,
                                                      const size_t threadNo) const {
        T1 tt = 0.0;
        
        T1 minDist = small;
        T1 s1, s2;
        
        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                return t0[ns];
            }
        }
        
        std::vector<bool> txOnNode( Tx.size(), false );
        std::vector<T2> txNode( Tx.size() );
        std::vector<T2> txCell( Tx.size() );
        std::vector<std::vector<T2>> txNeighborCells( Tx.size() );
        for ( size_t nt=0; nt<Tx.size(); ++nt ) {
            for ( T2 nn=0; nn<nodes.size(); ++nn ) {
                if ( nodes[nn] == Tx[nt] ) {
                    txOnNode[nt] = true;
                    txNode[nt] = nn;
                    break;
                }
            }
        }
#ifdef DEBUG_RP
        std::cout << "\n*** RP debug data - Source\n";
#endif
        for ( size_t nt=0; nt<Tx.size(); ++nt ) {
            if ( !txOnNode[nt] ) {
                txCell[nt] = getCellNo( Tx[nt] );
                
                std::array<T2,4> itmp = getPrimary(txCell[nt]);
                // find adjacent cells
                const T2 ind[6][2] = {
                    {itmp[0], itmp[1]},
                    {itmp[0], itmp[2]},
                    {itmp[0], itmp[3]},
                    {itmp[1], itmp[2]},
                    {itmp[1], itmp[3]},
                    {itmp[2], itmp[3]} };
                
                for ( size_t nedge=0; nedge<6; ++nedge ) {
                    for ( auto nc0=nodes[ind[nedge][0]].getOwners().begin(); nc0!=nodes[ind[nedge][0]].getOwners().end(); ++nc0 ) {
                        if ( std::find(nodes[ind[nedge][1]].getOwners().begin(), nodes[ind[nedge][1]].getOwners().end(), *nc0)!=nodes[ind[nedge][1]].getOwners().end() )
                            txNeighborCells[nt].push_back( *nc0 );
                    }
                }
            }
#ifdef DEBUG_RP
            std::cout << "   src no: " << nt << '\n';
            if ( txOnNode[nt] ) {
                std::cout << "     onNode\n"
                << "\t    txNode: " << txNode[nt] << '\n'
                << "\t    coords: " << nodes[txNode[nt]] << '\n';
            } else {
                std::array<T2,4> itmp = getPrimary(txCell[nt]);
                std::cout << "     inCell\n"
                << "\t    cellNo: " << txCell[nt] << '\n'
                << "\t  vertices: " << nodes[itmp[0]] << '\n'
                << "\t          : " << nodes[itmp[1]] << '\n'
                << "\t          : " << nodes[itmp[2]] << '\n'
                << "\t          : " << nodes[itmp[3]] << '\n';
            }
            std::cout << '\n';
#endif
        }
        
        T2 cellNo, nodeNo;
        sxyz<T1> curr_pt( Rx ), prev_pt( Rx );
        bool atRx = true;
        
        bool onNode = false;
        bool onEdge = false;
        bool onFace = false;
        std::array<T2,2> edgeNodes{ {0, 0} };
        std::array<T2,3> faceNodes{ {0, 0, 0} };
        Grad3D<T1,NODE>* grad3d = nullptr;
        if ( rp_method == 0 ) {
            grad3d = new Grad3D_ls_fo<T1,NODE>();
        } else if ( rp_method == 1 ) {
            grad3d = new Grad3D_ls_so<T1,NODE>();
        } else if ( rp_method == 2 ) {
            grad3d = new Grad3D_ab<T1,NODE>();
        }
        bool reachedTx = false;
        
        for ( T2 nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == curr_pt ) {
                nodeNo = nn;
                onNode = true;
                s1 = nodes[nodeNo].getNodeSlowness();
                break;
            }
        }
        if ( !onNode ) {
            cellNo = getCellNo( curr_pt );
            
            std::array<T2,4> itmp = getPrimary(cellNo);
            // find adjacent cells
            const T2 ind[6][2] = {
                {itmp[0], itmp[1]},
                {itmp[0], itmp[2]},
                {itmp[0], itmp[3]},
                {itmp[1], itmp[2]},
                {itmp[1], itmp[3]},
                {itmp[2], itmp[3]} };
            
            for ( size_t n=0; n<6; ++n ) {
                if ( areCollinear(curr_pt, ind[n][0], ind[n][1]) ) {
                    onEdge = true;
                    edgeNodes[0] = ind[n][0];
                    edgeNodes[1] = ind[n][1];
                    if ( interpVel )
                        s1 = Interpolator<T1>::linearVel(curr_pt,
                                                         nodes[edgeNodes[0]],
                                                         nodes[edgeNodes[1]]);
                    else
                        s1 = Interpolator<T1>::linear(curr_pt,
                                                      nodes[edgeNodes[0]],
                                                      nodes[edgeNodes[1]]);
                    break;
                }
            }
            if ( !onEdge ) {
                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                
                for ( size_t n=0; n<4; ++n ) {
                    if ( areCoplanar(curr_pt, ind[n][0], ind[n][1], ind[n][2]) ) {
                        onFace = true;
                        faceNodes[0] = ind[n][0];
                        faceNodes[1] = ind[n][1];
                        faceNodes[2] = ind[n][2];
                        if ( interpVel )
                            s1 = Interpolator<T1>::bilinearTriangleVel(curr_pt,
                                                                       nodes[faceNodes[0]],
                                                                       nodes[faceNodes[1]],
                                                                       nodes[faceNodes[2]]);
                        else
                            s1 = Interpolator<T1>::bilinearTriangle(curr_pt,
                                                                    nodes[faceNodes[0]],
                                                                    nodes[faceNodes[1]],
                                                                    nodes[faceNodes[2]]);
                        //  faceNodes shoud not be assigned, face was not intersected
                        faceNodes = {0, 0, 0};
                        break;
                    }
                }
                if ( !onFace ) {
                    if ( interpVel )
                        s1 = Interpolator<T1>::trilinearTriangleVel(curr_pt,
                                                                    nodes[itmp[0]],
                                                                    nodes[itmp[1]],
                                                                    nodes[itmp[2]],
                                                                    nodes[itmp[3]]);
                    else
                        s1 = Interpolator<T1>::trilinearTriangle(curr_pt,
                                                                 nodes[itmp[0]],
                                                                 nodes[itmp[1]],
                                                                 nodes[itmp[2]],
                                                                 nodes[itmp[3]]);

                }
            }
        }
        
        for ( auto nt=0; nt<txCell.size(); ++nt ) {
            if ( cellNo == txCell[nt] ) {
                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( interpVel )
                    s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                nodes[itmp[0]],
                                                                nodes[itmp[1]],
                                                                nodes[itmp[2]],
                                                                nodes[itmp[3]]);
                else
                    s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                             nodes[itmp[0]],
                                                             nodes[itmp[1]],
                                                             nodes[itmp[2]],
                                                             nodes[itmp[3]]);
                
                tt += t0[nt] + 0.5*(s1 + s2) * curr_pt.getDistance( Tx[nt] );
                reachedTx = true;
                break;
            }
        }
        
        sxyz<T1> g;
        while ( reachedTx == false ) {
            
            if ( onNode ) {
#ifdef DEBUG_RP
                printRaypathData(curr_pt, g, onNode, onEdge, onFace, cellNo,
                                 nodeNo, edgeNodes, faceNodes);
#endif
                if ( rp_method < 2 ) {
                    // find cells common to edge
                    std::set<NODE*> nnodes;
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                        getNeighborNodes(*nc, nnodes);
                    }
                    // compute gradient with nodes from all common cells
                    g = grad3d->compute(curr_pt, nodes[nodeNo].getTT(threadNo), nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(1);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[nodeNo]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);
                
                // find cell for which gradient intersect opposing face
                bool foundIntersection = false;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                    
                    std::array<T2,3> nb;
                    size_t n=0;
                    for (auto nn=neighbors[*nc].begin(); nn!=neighbors[*nc].end(); ++nn ) {
                        if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                            nb[n++] = *nn;
                        }
                    }
                    std::sort(nb.begin(), nb.end());
                    
                    foundIntersection = intersectVecTriangle( nodeNo, g, nb[0], nb[1], nb[2], curr_pt);
                    if ( !foundIntersection ) {
                        continue;
                    }
                    
                    bool break_flag = check_pt_location(curr_pt, nb, onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);
                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
                    s1 = s2;
                    prev_pt = curr_pt;
                    
                    if ( break_flag ) break;
                    
                    // find next cell
                    cellNo = findAdjacentCell1(faceNodes, nodeNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        reachedTx = true;
                    }
                    break;
                }
                
                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath on node failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    tt = 0.0;
                    reachedTx = true;
                }
                
            } else if ( onEdge ) {
#ifdef DEBUG_RP
                printRaypathData(curr_pt, g, onNode, onEdge, onFace, cellNo,
                                 nodeNo, edgeNodes, faceNodes);
#endif
                // find cells common to edge
                std::vector<T2> cells;
                for ( auto nc0=nodes[edgeNodes[0]].getOwners().begin(); nc0!=nodes[edgeNodes[0]].getOwners().end(); ++nc0 ) {
                    if ( std::find(nodes[edgeNodes[1]].getOwners().begin(), nodes[edgeNodes[1]].getOwners().end(), *nc0)!=nodes[edgeNodes[1]].getOwners().end() ) {
                        cells.push_back( *nc0 );
                    }
                }
                if ( rp_method < 2 ) {
                    std::set<NODE*> nnodes;
                    for (size_t n=0; n<cells.size(); ++n ) {
                        getNeighborNodes(cells[n], nnodes);
                    }
                    T1 d01 = nodes[edgeNodes[0]].getDistance(nodes[edgeNodes[1]]);
                    T1 w0 = curr_pt.getDistance(nodes[edgeNodes[1]]) / d01;
                    T1 w1 = curr_pt.getDistance(nodes[edgeNodes[0]]) / d01;
                    T1 curr_t = nodes[edgeNodes[0]].getTT(threadNo)*w0 + nodes[edgeNodes[1]].getTT(threadNo)*w1;
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(2);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[edgeNodes[0]]);
                    ref_pt[1] = &(nodes[edgeNodes[1]]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);
                
                bool foundIntersection = false;
                for (size_t n=0; n<cells.size(); ++n ) {
                    
                    cellNo = cells[n];
                    
                    // there are 2 faces that might be intersected
                    std::array<T2,2> edgeNodes2;
                    size_t n2=0;
                    for ( auto nn=neighbors[cellNo].begin(); nn!= neighbors[cellNo].end(); ++nn ) {
                        if ( *nn!=edgeNodes[0] && *nn!=edgeNodes[1] && nodes[*nn].isPrimary() ) {
                            edgeNodes2[n2++] = *nn;
                        }
                    }
                    
                    sxyz<T1> pt_i;
                    T2 itmpNode;
                    foundIntersection = intersectVecTriangle(curr_pt, g,
                                                             edgeNodes[0],
                                                             edgeNodes2[0],
                                                             edgeNodes2[1], pt_i);
                    itmpNode = edgeNodes[0];
                    if ( !foundIntersection ) {
                        foundIntersection = intersectVecTriangle(curr_pt, g,
                                                                 edgeNodes[1],
                                                                 edgeNodes2[0],
                                                                 edgeNodes2[1], pt_i);
                        itmpNode = edgeNodes[1];
                    }
                    if ( !foundIntersection ) {
                        continue;
                    }
                    
                    curr_pt = pt_i;
                    
                    bool break_flag = check_pt_location(curr_pt, neighbors[cellNo],
                                                        {itmpNode, edgeNodes2[0], edgeNodes2[1]},
                                                        onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
                    s1 = s2;
                    prev_pt = curr_pt;
                    
                    if ( break_flag ) break;
                    
                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        reachedTx = true;
                    }
                    break;
                }
                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath on edge failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    tt = 0.0;
                    reachedTx = true;
                }
                
            } else if ( onFace ) { // on Face
#ifdef DEBUG_RP
                printRaypathData(curr_pt, g, onNode, onEdge, onFace, cellNo,
                                 nodeNo, edgeNodes, faceNodes);
#endif
                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( rp_method < 2 ) {
                    std::set<NODE*> nnodes;
                    getNeighborNodes(cellNo, nnodes);
                    T1 curr_t;
                    if ( atRx ) {
                        curr_t = Interpolator<T1>::trilinearTime(curr_pt,
                                                                 nodes[itmp[0]],
                                                                 nodes[itmp[1]],
                                                                 nodes[itmp[2]],
                                                                 nodes[itmp[3]],
                                                                 threadNo);
                        atRx = false;
                    } else {
                        curr_t = Interpolator<T1>::bilinearTime(curr_pt,
                                                                nodes[faceNodes[0]],
                                                                nodes[faceNodes[1]],
                                                                nodes[faceNodes[2]],
                                                                threadNo);
                    }
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(3);
                    if ( atRx ) {
                        ref_pt[0] = &(nodes[itmp[0]]);
                        ref_pt[1] = &(nodes[itmp[1]]);
                        ref_pt[2] = &(nodes[itmp[2]]);
                        ref_pt.push_back( &(nodes[itmp[3]]) );
                        atRx = false;
                    } else {
                        ref_pt[0] = &(nodes[faceNodes[0]]);
                        ref_pt[1] = &(nodes[faceNodes[1]]);
                        ref_pt[2] = &(nodes[faceNodes[2]]);
                    }
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                // there are 3 faces that might be intersected
                
                bool foundIntersection = false;
                for ( size_t n=0; n<4; ++n ) {
                    if ( ind[n] == faceNodes ) continue;
                    
                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);
                    
                    if ( !foundIntersection )
                        continue;
                    
                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);
                    tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
                    s1 = s2;
                    prev_pt = curr_pt;
                    
                    if ( break_flag ) break;
                    
                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        reachedTx = true;
                    }
                    break;
                }
                
                if ( foundIntersection == false ) {
                    
                    // we must be on an face with gradient pointing slightly outward tetrahedron
                    // return in other cell but keep gradient
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    
                    std::array<T2,4> itmp = getPrimary(cellNo);
                    ind[0] = { { itmp[0], itmp[1], itmp[2] } };
                    ind[1] = { { itmp[0], itmp[1], itmp[3] } };
                    ind[2] = { { itmp[0], itmp[2], itmp[3] } };
                    ind[3] = { { itmp[1], itmp[2], itmp[3] } };
                    
                    for ( size_t n=0; n<4; ++n )
                        std::sort( ind[n].begin(), ind[n].end() );
                    
                    for ( size_t n=0; n<4; ++n ) {
                        if ( ind[n] == faceNodes ) continue;
                        
                        sxyz<T1> pt_i;
                        foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                                 ind[n][1], ind[n][2],
                                                                 pt_i);
                        
                        if ( !foundIntersection ) {
                            continue;
                        }
                        curr_pt = pt_i;

                        bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                            nodeNo, onEdge, edgeNodes,
                                                            onFace, faceNodes);

                        s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                             faceNodes);

                        tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
                        s1 = s2;
                        prev_pt = curr_pt;
                        
                        if ( break_flag ) break;
                        
                        // find next cell
                        cellNo = findAdjacentCell2(faceNodes, cellNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                            << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                            tt = 0.0;
                            reachedTx = true;
                        }
                        break;
                    }
                }
                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath on face failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    tt = 0.0;
                    reachedTx = true;
                }
            } else { // at Rx, somewhere in a tetrahedron
#ifdef DEBUG_RP
                printRaypathData(curr_pt, g, onNode, onEdge, onFace, cellNo,
                                 nodeNo, edgeNodes, faceNodes);
#endif
                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( rp_method < 2 ) {
                    std::set<NODE*> nnodes;
                    getNeighborNodes(cellNo, nnodes);
                    T1 curr_t = Interpolator<T1>::trilinearTime(curr_pt,
                                                                nodes[itmp[0]],
                                                                nodes[itmp[1]],
                                                                nodes[itmp[2]],
                                                                nodes[itmp[3]],
                                                                threadNo);
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(4);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[itmp[0]]);
                    ref_pt[1] = &(nodes[itmp[1]]);
                    ref_pt[2] = &(nodes[itmp[2]]);
                    ref_pt[3] = &(nodes[itmp[3]]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);
                
                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                
                bool foundIntersection = false;
                for ( size_t n=0; n<4; ++n ) {
                    
                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);
                    
                    if ( !foundIntersection )
                        continue;
                    
                    curr_pt = pt_i;
                    
                    bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);
                    
                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);
                    tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
                    s1 = s2;
                    prev_pt = curr_pt;
                    
                    if ( break_flag ) break;
                    
                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        reachedTx = true;
                    }
                    break;
                }
                
                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath in cell failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    tt = 0.0;
                    reachedTx = true;
                }
            }
            
            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( curr_pt.getDistance( Tx[nt] ) < minDist ) {
                        tt += t0[nt];
                        reachedTx = true;
                        break;
                    }
                }
            } else {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        for ( auto nc=nodes[txNode[nt]].getOwners().begin();
                             nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                            if ( cellNo == *nc ) {
                                std::array<T2,4> itmp = getPrimary(cellNo);

                                if ( interpVel )
                                    s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                nodes[itmp[0]],
                                                                                nodes[itmp[1]],
                                                                                nodes[itmp[2]],
                                                                                nodes[itmp[3]]);
                                else
                                    s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                             nodes[itmp[0]],
                                                                             nodes[itmp[1]],
                                                                             nodes[itmp[2]],
                                                                             nodes[itmp[3]]);

                                tt += t0[nt] + 0.5*(s1 + s2) * prev_pt.getDistance( Tx[nt] );
                                reachedTx = true;
                                break;
                            }
                        }
                    } else {
                        if ( cellNo == txCell[nt] ) {
                            std::array<T2,4> itmp = getPrimary(cellNo);
                            
                            if ( interpVel )
                                s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                            nodes[itmp[0]],
                                                                            nodes[itmp[1]],
                                                                            nodes[itmp[2]],
                                                                            nodes[itmp[3]]);
                            else
                                s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                         nodes[itmp[0]],
                                                                         nodes[itmp[1]],
                                                                         nodes[itmp[2]],
                                                                         nodes[itmp[3]]);

                            tt += t0[nt] + 0.5*(s1 + s2) * prev_pt.getDistance( Tx[nt] );
                            reachedTx = true;
                        } else {
                            for ( size_t nn=0; nn<txNeighborCells[nt].size(); ++nn ) {
                                if ( cellNo == txNeighborCells[nt][nn] ) {
                                    std::array<T2,4> itmp = getPrimary(txCell[nt]);
                                    
                                    bool found = false;
                                    for ( T2 i=0; i<4; ++i ) {
                                        T2 iA = itmp[(i+1)%4];
                                        T2 iB = itmp[(i+2)%4];
                                        T2 iC = itmp[(i+3)%4];
                                        if (nodes[iA].getDistance(curr_pt)+nodes[iB].getDistance(curr_pt)-nodes[iB].getDistance(nodes[iA])<minDist*minDist||
                                            nodes[iA].getDistance(curr_pt)+nodes[iC].getDistance(curr_pt)-nodes[iC].getDistance(nodes[iA])<minDist*minDist||
                                            nodes[iC].getDistance(curr_pt)+nodes[iB].getDistance(curr_pt)-nodes[iB].getDistance(nodes[iC])<minDist*minDist){
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found ) {
                                        std::array<T2,4> itmp = getPrimary(cellNo);
                                        
                                        if ( interpVel )
                                            s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                        nodes[itmp[0]],
                                                                                        nodes[itmp[1]],
                                                                                        nodes[itmp[2]],
                                                                                        nodes[itmp[3]]);
                                        else
                                            s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                                     nodes[itmp[0]],
                                                                                     nodes[itmp[1]],
                                                                                     nodes[itmp[2]],
                                                                                     nodes[itmp[3]]);
                                        
                                        tt += t0[nt] + 0.5*(s1 + s2) * prev_pt.getDistance( Tx[nt] );
                                        reachedTx = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if ( reachedTx ) break;
                }
            }
        }
        delete grad3d;
        return tt;
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::printRaypathData(const sxyz<T1>& curr_pt,
                                                const sxyz<T1>& g,
                                                const bool onNode,
                                                const bool onEdge,
                                                const bool onFace,
                                                const T2 cellNo,
                                                const T2 nodeNo,
                                                const std::array<T2,2> &edgeNodes,
                                                const std::array<T2,3> &faceNodes) const {
        std::array<T2,4> itmp = getPrimary(cellNo);
        std::cout << "\n*** RP debug data\n   curr_pt: " << curr_pt << '\n'
        << "         g: " << g << '\n'
        << "    cellNo: " << cellNo << '\n'
        << "  vertices: " << nodes[itmp[0]] << '\n'
        << "          : " << nodes[itmp[1]] << '\n'
        << "          : " << nodes[itmp[2]] << '\n'
        << "          : " << nodes[itmp[3]] << '\n';
        if ( onNode ) {
            std::cout << "\tonNode\n"
            << "\t    nodeNo: " << nodeNo << '\n'
            << "\t    coords: " << nodes[nodeNo] << '\n';
        }
        if ( onEdge ) {
            std::cout << "\tonEdge\n"
            << "\t    edgeNo: " << edgeNodes[0] << ' ' << edgeNodes[1] << '\n'
            << "\t    coords: " << nodes[edgeNodes[0]] << '\n'
            << "\t    coords: " << nodes[edgeNodes[1]] << '\n';
        }
        if ( onFace ) {
            std::cout << "\tonFace\n"
            << "\t    faceNo: " << faceNodes[0] << ' ' << faceNodes[1] << ' ' << faceNodes[2] << '\n'
            << "\t    coords: " << nodes[faceNodes[0]] << '\n'
            << "\t    coords: " << nodes[faceNodes[1]] << '\n'
            << "\t    coords: " << nodes[faceNodes[2]] << '\n';
        }
        std::cout << std::endl;
    }
    
    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                          const sxyz<T1> &Rx,
                                          std::vector<sxyz<T1>> &r_data,
                                          const size_t threadNo) const {

        T1 minDist = small;
        std::vector<sxyz<T1>> r_tmp;
        r_tmp.emplace_back( Rx );
        
        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                return;
            }
        }
        
        std::vector<bool> txOnNode( Tx.size(), false );
        std::vector<T2> txNode( Tx.size() );
        std::vector<T2> txCell( Tx.size() );
        std::vector<std::vector<T2>> txNeighborCells( Tx.size() );
        for ( size_t nt=0; nt<Tx.size(); ++nt ) {
            for ( T2 nn=0; nn<nodes.size(); ++nn ) {
                if ( nodes[nn] == Tx[nt] ) {
                    txOnNode[nt] = true;
                    txNode[nt] = nn;
                    break;
                }
            }
        }
        for ( size_t nt=0; nt<Tx.size(); ++nt ) {
            if ( !txOnNode[nt] ) {
                txCell[nt] = getCellNo( Tx[nt] );
                
                std::array<T2,4> itmp = getPrimary(txCell[nt]);
                // find adjacent cells
                const T2 ind[6][2] = {
                    {itmp[0], itmp[1]},
                    {itmp[0], itmp[2]},
                    {itmp[0], itmp[3]},
                    {itmp[1], itmp[2]},
                    {itmp[1], itmp[3]},
                    {itmp[2], itmp[3]} };
                
                for ( size_t nedge=0; nedge<6; ++nedge ) {
                    for ( auto nc0=nodes[ind[nedge][0]].getOwners().begin(); nc0!=nodes[ind[nedge][0]].getOwners().end(); ++nc0 ) {
                        if ( std::find(nodes[ind[nedge][1]].getOwners().begin(), nodes[ind[nedge][1]].getOwners().end(), *nc0)!=nodes[ind[nedge][1]].getOwners().end() )
                            txNeighborCells[nt].push_back( *nc0 );
                    }
                }
            }
        }
        
        T2 cellNo, nodeNo;
        sxyz<T1> curr_pt( Rx );
        
        bool onNode = false;
        bool onEdge = false;
        bool onFace = false;
        std::array<T2,2> edgeNodes{ {0, 0} };
        std::array<T2,3> faceNodes{ {0, 0, 0} };
        Grad3D<T1,NODE>* grad3d = nullptr;
        if ( rp_method == 0 ) {
            grad3d = new Grad3D_ls_fo<T1,NODE>();
        } else if ( rp_method == 1 ) {
            grad3d = new Grad3D_ls_so<T1,NODE>();
        } else if ( rp_method == 2 ) {
            grad3d = new Grad3D_ab<T1,NODE>();
        }
        bool reachedTx = false;
        
        for ( T2 nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == curr_pt ) {
                nodeNo = nn;
                onNode = true;
                break;
            }
        }
        if ( !onNode ) {
            cellNo = getCellNo( curr_pt );
            
            std::array<T2,4> itmp = getPrimary(cellNo);
            T2 ind[6][2] = {
                {itmp[0], itmp[1]},
                {itmp[0], itmp[2]},
                {itmp[0], itmp[3]},
                {itmp[1], itmp[2]},
                {itmp[1], itmp[3]},
                {itmp[2], itmp[3]} };
            
            for ( size_t n=0; n<6; ++n ) {
                if ( areCollinear(curr_pt, ind[n][0], ind[n][1]) ) {
                    onEdge = true;
                    edgeNodes[0] = ind[n][0];
                    edgeNodes[1] = ind[n][1];
                    break;
                }
            }
            if ( !onEdge ) {
                std::array<T2,4> itmp = getPrimary(cellNo);
                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                
                for ( size_t n=0; n<4; ++n ) {
                    if ( areCoplanar(curr_pt, ind[n][0], ind[n][1], ind[n][2]) ) {
                        onFace = true;
                        //  faceNodes shoud not be assigned, face was not intersected
                        break;
                    }
                }
            }
        }
        
        for ( auto nt=0; nt<txCell.size(); ++nt ) {
            if ( cellNo == txCell[nt] ) {
                r_tmp.push_back( Tx[nt] );
                reachedTx = true;
                break;
            }
        }
        
        sxyz<T1> g;
        while ( reachedTx == false ) {
            
            if ( onNode ) {
                
                if ( rp_method < 2 ) {
                    // find cells common to edge
                    std::set<NODE*> nnodes;
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                        getNeighborNodes(*nc, nnodes);
                    }
                    // compute gradient with nodes from all common cells
                    g = grad3d->compute(curr_pt, nodes[nodeNo].getTT(threadNo), nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(1);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[nodeNo]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                // find cell for which gradient intersect opposing face
                bool foundIntersection = false;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                    
                    std::array<T2,3> nb;
                    size_t n=0;
                    for (auto nn=neighbors[*nc].begin(); nn!=neighbors[*nc].end(); ++nn ) {
                        if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                            nb[n++] = *nn;
                        }
                    }
                    std::sort(nb.begin(), nb.end());
                    
                    foundIntersection = intersectVecTriangle( nodeNo, g, nb[0], nb[1], nb[2], curr_pt);
                    if ( !foundIntersection ) {
                        continue;
                    }
                    
                    bool break_flag = check_pt_location(curr_pt, nb, onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    r_tmp.push_back( curr_pt );
                    
                    if ( break_flag ) break;
                    
                    // find next cell
                    cellNo = findAdjacentCell1(faceNodes, nodeNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }
                
                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath on node failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    reachedTx = true;
                }
                
            } else if ( onEdge ) {
                
                // find cells common to edge
                std::vector<T2> cells;
                for ( auto nc0=nodes[edgeNodes[0]].getOwners().begin(); nc0!=nodes[edgeNodes[0]].getOwners().end(); ++nc0 ) {
                    if ( std::find(nodes[edgeNodes[1]].getOwners().begin(), nodes[edgeNodes[1]].getOwners().end(), *nc0)!=nodes[edgeNodes[1]].getOwners().end() ) {
                        cells.push_back( *nc0 );
                    }
                }
                if ( rp_method < 2 ) {
                    std::set<NODE*> nnodes;
                    for (size_t n=0; n<cells.size(); ++n ) {
                        getNeighborNodes(cells[n], nnodes);
                    }
                    T1 d01 = nodes[edgeNodes[0]].getDistance(nodes[edgeNodes[1]]);
                    T1 w0 = curr_pt.getDistance(nodes[edgeNodes[1]]) / d01;
                    T1 w1 = curr_pt.getDistance(nodes[edgeNodes[0]]) / d01;
                    T1 curr_t = nodes[edgeNodes[0]].getTT(threadNo)*w0 + nodes[edgeNodes[1]].getTT(threadNo)*w1;
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(2);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[edgeNodes[0]]);
                    ref_pt[1] = &(nodes[edgeNodes[1]]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                bool foundIntersection = false;
                for (size_t n=0; n<cells.size(); ++n ) {
                    
                    cellNo = cells[n];
                    
                    // there are 2 faces that might be intersected
                    std::array<T2,2> edgeNodes2;
                    size_t n2=0;
                    for ( auto nn=neighbors[cellNo].begin(); nn!= neighbors[cellNo].end(); ++nn ) {
                        if ( *nn!=edgeNodes[0] && *nn!=edgeNodes[1] && nodes[*nn].isPrimary() ) {
                            edgeNodes2[n2++] = *nn;
                        }
                    }
                    
                    sxyz<T1> pt_i;
                    T2 itmpNode;
                    foundIntersection = intersectVecTriangle(curr_pt, g,
                                                             edgeNodes[0],
                                                             edgeNodes2[0],
                                                             edgeNodes2[1], pt_i);
                    itmpNode = edgeNodes[0];
                    if ( !foundIntersection ) {
                        foundIntersection = intersectVecTriangle(curr_pt, g,
                                                                 edgeNodes[1],
                                                                 edgeNodes2[0],
                                                                 edgeNodes2[1], pt_i);
                        itmpNode = edgeNodes[1];
                    }
                    if ( !foundIntersection ) {
                        continue;
                    }
                    
                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, neighbors[cellNo],
                                                        {itmpNode, edgeNodes2[0], edgeNodes2[1]},
                                                        onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    r_tmp.push_back( curr_pt );
                    
                    if ( break_flag ) break;
                    
                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }
                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath on edge failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    reachedTx = true;
                }
                
            } else if ( onFace ) { // on Face
                
                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( rp_method < 2 ) {
                    std::set<NODE*> nnodes;
                    getNeighborNodes(cellNo, nnodes);
                    T1 curr_t;
                    if ( r_tmp.size <= 1 ) {
                        curr_t = Interpolator<T1>::trilinearTime(curr_pt,
                                                                 nodes[itmp[0]],
                                                                 nodes[itmp[1]],
                                                                 nodes[itmp[2]],
                                                                 nodes[itmp[3]],
                                                                 threadNo);
                    } else {
                        curr_t = Interpolator<T1>::bilinearTime(curr_pt,
                                                                nodes[faceNodes[0]],
                                                                nodes[faceNodes[1]],
                                                                nodes[faceNodes[2]],
                                                                threadNo);
                    }
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(3);
                    if ( r_tmp.size <= 1 ) {
                        ref_pt[0] = &(nodes[itmp[0]]);
                        ref_pt[1] = &(nodes[itmp[1]]);
                        ref_pt[2] = &(nodes[itmp[2]]);
                        ref_pt.push_back( &(nodes[itmp[3]]) );
                    } else {
                        ref_pt[0] = &(nodes[faceNodes[0]]);
                        ref_pt[1] = &(nodes[faceNodes[1]]);
                        ref_pt[2] = &(nodes[faceNodes[2]]);
                    }
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                // there are 3 faces that might be intersected
                
                bool foundIntersection = false;
                for ( size_t n=0; n<4; ++n ) {
                    if ( ind[n] == faceNodes ) continue;
                    
                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);
                    
                    if ( !foundIntersection )
                        continue;
                    
                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    r_tmp.push_back( curr_pt );
                    
                    if ( break_flag ) break;
                    
                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }
                
                if ( foundIntersection == false ) {
                    
                    // we must be on an face with gradient pointing slightly outward tetrahedron
                    // return in other cell but keep gradient
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    
                    std::array<T2,4> itmp = getPrimary(cellNo);
                    ind[0] = { { itmp[0], itmp[1], itmp[2] } };
                    ind[1] = { { itmp[0], itmp[1], itmp[3] } };
                    ind[2] = { { itmp[0], itmp[2], itmp[3] } };
                    ind[3] = { { itmp[1], itmp[2], itmp[3] } };
                    
                    for ( size_t n=0; n<4; ++n )
                        std::sort( ind[n].begin(), ind[n].end() );
                    
                    for ( size_t n=0; n<4; ++n ) {
                        if ( ind[n] == faceNodes ) continue;
                        
                        sxyz<T1> pt_i;
                        foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                                 ind[n][1], ind[n][2],
                                                                 pt_i);
                        
                        if ( !foundIntersection ) {
                            continue;
                        }
                        curr_pt = pt_i;

                        bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                            nodeNo, onEdge, edgeNodes,
                                                            onFace, faceNodes);

                        r_tmp.push_back( curr_pt );
                        
                        if ( break_flag ) break;
                        
                        // find next cell
                        cellNo = findAdjacentCell2(faceNodes, cellNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                            << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                            r_tmp.resize(1);
                            r_tmp[0] = Rx;
                            reachedTx = true;
                        }
                        break;
                    }
                }
                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath on face failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    reachedTx = true;
                }
            } else { // at Rx, somewhere in a tetrahedron
                
                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( rp_method < 2 ) {
                    std::set<NODE*> nnodes;
                    getNeighborNodes(cellNo, nnodes);
                    T1 curr_t = Interpolator<T1>::trilinearTime(curr_pt,
                                                                nodes[itmp[0]],
                                                                nodes[itmp[1]],
                                                                nodes[itmp[2]],
                                                                nodes[itmp[3]],
                                                                threadNo);
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(4);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[itmp[0]]);
                    ref_pt[1] = &(nodes[itmp[1]]);
                    ref_pt[2] = &(nodes[itmp[2]]);
                    ref_pt[3] = &(nodes[itmp[3]]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                // there are 3 faces that might be intersected
                
                bool foundIntersection = false;
                for ( size_t n=0; n<4; ++n ) {
                    
                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);
                    
                    if ( !foundIntersection )
                        continue;
                    
                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    r_tmp.push_back( curr_pt );
                    
                    if ( break_flag ) break;
                    
                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }
                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath within cell failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_data.resize(1);
                    r_data[0] = Rx;
                    reachedTx = true;
                }
            }
            
            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( curr_pt.getDistance( Tx[nt] ) < minDist ) {
                        reachedTx = true;
                        break;
                    }
                }
            } else {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        for ( auto nc=nodes[txNode[nt]].getOwners().begin();
                             nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                            if ( cellNo == *nc ) {
                                r_tmp.push_back( Tx[nt] );
                                reachedTx = true;
                                break;
                            }
                        }
                    } else {
                        if ( cellNo == txCell[nt] ) {
                            r_tmp.push_back( Tx[nt] );
                            reachedTx = true;
                        } else {
                            for ( size_t nn=0; nn<txNeighborCells[nt].size(); ++nn ) {
                                if ( cellNo == txNeighborCells[nt][nn] ) {
                                    std::array<T2,4> itmp = getPrimary(txCell[nt]);
                                    
                                    bool found = false;
                                    for ( T2 i=0; i<4; ++i ) {
                                        T2 iA = itmp[(i+1)%4];
                                        T2 iB = itmp[(i+2)%4];
                                        T2 iC = itmp[(i+3)%4];
                                        if (nodes[iA].getDistance(curr_pt)+nodes[iB].getDistance(curr_pt)-nodes[iB].getDistance(nodes[iA])<minDist*minDist||
                                            nodes[iA].getDistance(curr_pt)+nodes[iC].getDistance(curr_pt)-nodes[iC].getDistance(nodes[iA])<minDist*minDist||
                                            nodes[iC].getDistance(curr_pt)+nodes[iB].getDistance(curr_pt)-nodes[iB].getDistance(nodes[iC])<minDist*minDist){
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found ) {
                                        r_tmp.push_back( Tx[nt] );
                                        reachedTx = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if ( reachedTx ) break;
                }
            }
        }
        // for inversion, the order should be from Tx to Rx, so we reorder...
        size_t npts = r_tmp.size();
        r_data.resize( npts );
        for ( size_t nn=0; nn<npts; ++nn ) {
            r_data[nn] = r_tmp[ npts-1-nn ];
        }
        
        delete grad3d;
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                          const std::vector<T1>& t0,
                                          const sxyz<T1> &Rx,
                                          std::vector<sxyz<T1>> &r_data,
                                          T1 &tt,
                                          const size_t threadNo) const {
        
        T1 minDist = small;
        std::vector<sxyz<T1>> r_tmp;
        r_tmp.emplace_back( Rx );
        tt = 0.0;
        T1 s1, s2;
        
        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }
        
        std::vector<bool> txOnNode( Tx.size(), false );
        std::vector<T2> txNode( Tx.size() );
        std::vector<T2> txCell( Tx.size() );
        std::vector<std::vector<T2>> txNeighborCells( Tx.size() );
        for ( size_t nt=0; nt<Tx.size(); ++nt ) {
            for ( T2 nn=0; nn<nodes.size(); ++nn ) {
                if ( nodes[nn] == Tx[nt] ) {
                    txOnNode[nt] = true;
                    txNode[nt] = nn;
                    break;
                }
            }
        }
#ifdef DEBUG_RP
        std::cout << "\n*** RP debug data - Source\n";
#endif
        for ( size_t nt=0; nt<Tx.size(); ++nt ) {
            if ( !txOnNode[nt] ) {
                txCell[nt] = getCellNo( Tx[nt] );
                
                std::array<T2,4> itmp = getPrimary(txCell[nt]);
                // find adjacent cells
                T2 ind[6][2] = {
                    {itmp[0], itmp[1]},
                    {itmp[0], itmp[2]},
                    {itmp[0], itmp[3]},
                    {itmp[1], itmp[2]},
                    {itmp[1], itmp[3]},
                    {itmp[2], itmp[3]} };
                
                for ( size_t nedge=0; nedge<6; ++nedge ) {
                    for ( auto nc0=nodes[ind[nedge][0]].getOwners().begin(); nc0!=nodes[ind[nedge][0]].getOwners().end(); ++nc0 ) {
                        if ( std::find(nodes[ind[nedge][1]].getOwners().begin(), nodes[ind[nedge][1]].getOwners().end(), *nc0)!=nodes[ind[nedge][1]].getOwners().end() )
                            txNeighborCells[nt].push_back( *nc0 );
                    }
                }
            }
#ifdef DEBUG_RP
            std::cout << "   src no: " << nt << '\n';
            if ( txOnNode[nt] ) {
                std::cout << "     onNode\n"
                << "\t    txNode: " << txNode[nt] << '\n'
                << "\t    coords: " << nodes[txNode[nt]] << '\n';
            } else {
                std::array<T2,4> itmp = getPrimary(txCell[nt]);
                std::cout << "     inCell\n"
                << "\t    cellNo: " << txCell[nt] << '\n'
                << "\t  vertices: " << nodes[itmp[0]] << '\n'
                << "\t          : " << nodes[itmp[1]] << '\n'
                << "\t          : " << nodes[itmp[2]] << '\n'
                << "\t          : " << nodes[itmp[3]] << '\n';
            }
            std::cout << '\n';
#endif
        }
        
        T2 cellNo, nodeNo;
        sxyz<T1> curr_pt( Rx );
        
        bool onNode = false;
        bool onEdge = false;
        bool onFace = false;
        std::array<T2,2> edgeNodes{ {0, 0} };
        std::array<T2,3> faceNodes{ {0, 0, 0} };
        Grad3D<T1,NODE>* grad3d = nullptr;
        if ( rp_method == 0 ) {
            grad3d = new Grad3D_ls_fo<T1,NODE>();
        } else if ( rp_method == 1 ) {
            grad3d = new Grad3D_ls_so<T1,NODE>();
        } else if ( rp_method == 2 ) {
            grad3d = new Grad3D_ab<T1,NODE>();
        }
        bool reachedTx = false;
        
        for ( T2 nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == curr_pt ) {
                nodeNo = nn;
                onNode = true;
                s1 = nodes[nodeNo].getNodeSlowness();
                break;
            }
        }
        if ( !onNode ) {
            cellNo = getCellNo( curr_pt );
            
            std::array<T2,4> itmp = getPrimary(cellNo);
            T2 ind[6][2] = {
                {itmp[0], itmp[1]},
                {itmp[0], itmp[2]},
                {itmp[0], itmp[3]},
                {itmp[1], itmp[2]},
                {itmp[1], itmp[3]},
                {itmp[2], itmp[3]} };
            
            for ( size_t n=0; n<6; ++n ) {
                if ( areCollinear(curr_pt, ind[n][0], ind[n][1]) ) {
                    onEdge = true;
                    edgeNodes[0] = ind[n][0];
                    edgeNodes[1] = ind[n][1];
                    
                    if ( interpVel )
                        s1 = Interpolator<T1>::linearVel(r_tmp.back(),
                                                         nodes[edgeNodes[0]],
                                                         nodes[edgeNodes[1]]);
                    else
                        s1 = Interpolator<T1>::linear(r_tmp.back(),
                                                      nodes[edgeNodes[0]],
                                                      nodes[edgeNodes[1]]);

                    break;
                }
            }
            if ( !onEdge ) {
                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                
                for ( size_t n=0; n<4; ++n ) {
                    if ( areCoplanar(curr_pt, ind[n][0], ind[n][1], ind[n][2]) ) {
                        onFace = true;
                        faceNodes[0] = ind[n][0];
                        faceNodes[1] = ind[n][1];
                        faceNodes[2] = ind[n][2];

                        if ( interpVel )
                            s1 = Interpolator<T1>::bilinearTriangleVel(curr_pt,
                                                                       nodes[faceNodes[0]],
                                                                       nodes[faceNodes[1]],
                                                                       nodes[faceNodes[2]]);
                        else
                            s1 = Interpolator<T1>::bilinearTriangle(curr_pt,
                                                                    nodes[faceNodes[0]],
                                                                    nodes[faceNodes[1]],
                                                                    nodes[faceNodes[2]]);
                        //  faceNodes shoud not be assigned, face was not intersected
                        faceNodes = {0, 0, 0};
                        break;
                    }
                }
                if ( !onFace ) {
                    if ( interpVel )
                        s1 = Interpolator<T1>::trilinearTriangleVel(curr_pt,
                                                                    nodes[itmp[0]],
                                                                    nodes[itmp[1]],
                                                                    nodes[itmp[2]],
                                                                    nodes[itmp[3]]);
                    else
                        s1 = Interpolator<T1>::trilinearTriangle(curr_pt,
                                                                 nodes[itmp[0]],
                                                                 nodes[itmp[1]],
                                                                 nodes[itmp[2]],
                                                                 nodes[itmp[3]]);

                }
            }
        }
        
        for ( auto nt=0; nt<txCell.size(); ++nt ) {
            if ( cellNo == txCell[nt] ) {
                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( interpVel )
                    s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                nodes[itmp[0]],
                                                                nodes[itmp[1]],
                                                                nodes[itmp[2]],
                                                                nodes[itmp[3]]);
                else
                    s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                             nodes[itmp[0]],
                                                             nodes[itmp[1]],
                                                             nodes[itmp[2]],
                                                             nodes[itmp[3]]);
                
                tt += t0[nt] + 0.5*(s1 + s2) * r_tmp.back().getDistance( Tx[nt] );
                r_tmp.push_back( Tx[nt] );
                reachedTx = true;
                break;
            }
        }
        
        sxyz<T1> g;
        while ( reachedTx == false ) {
            
            if ( onNode ) {
#ifdef DEBUG_RP
                printRaypathData(curr_pt, g, onNode, onEdge, onFace, cellNo,
                                 nodeNo, edgeNodes, faceNodes);
#endif
                if ( rp_method < 2 ) {
                    // find cells common to edge
                    std::set<NODE*> nnodes;
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                        getNeighborNodes(*nc, nnodes);
                    }
                    // compute gradient with nodes from all common cells
                    g = grad3d->compute(curr_pt, nodes[nodeNo].getTT(threadNo), nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(1);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[nodeNo]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                // find cell for which gradient intersect opposing face
                bool foundIntersection = false;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                    
                    std::array<T2,3> nb;
                    size_t n=0;
                    for (auto nn=neighbors[*nc].begin(); nn!=neighbors[*nc].end(); ++nn ) {
                        if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                            nb[n++] = *nn;
                        }
                    }
                    std::sort(nb.begin(), nb.end());
                    
                    foundIntersection = intersectVecTriangle( nodeNo, g, nb[0], nb[1], nb[2], curr_pt);
                    if ( !foundIntersection ) {
                        continue;
                    }
                    
                    bool break_flag = check_pt_location(curr_pt, nb, onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    tt += 0.5*(s1 + s2) * r_tmp.back().getDistance( curr_pt );
//                    std::cout << curr_pt << '\t' << tt << '\t' << s1 << '\t' << s2 << '\n';
                    s1 = s2;
                    r_tmp.push_back( curr_pt );
                    
                    if ( break_flag ) break;
                    
                    // find next cell
                    cellNo = findAdjacentCell1(faceNodes, nodeNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }
                
                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath on node failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    reachedTx = true;
                }
                
            } else if ( onEdge ) {
#ifdef DEBUG_RP
                printRaypathData(curr_pt, g, onNode, onEdge, onFace, cellNo,
                                 nodeNo, edgeNodes, faceNodes);
#endif
                // find cells common to edge
                std::vector<T2> cells;
                for ( auto nc0=nodes[edgeNodes[0]].getOwners().begin(); nc0!=nodes[edgeNodes[0]].getOwners().end(); ++nc0 ) {
                    if ( std::find(nodes[edgeNodes[1]].getOwners().begin(), nodes[edgeNodes[1]].getOwners().end(), *nc0)!=nodes[edgeNodes[1]].getOwners().end() ) {
                        cells.push_back( *nc0 );
                    }
                }
                if ( rp_method < 2 ) {
                    std::set<NODE*> nnodes;
                    for (size_t n=0; n<cells.size(); ++n ) {
                        getNeighborNodes(cells[n], nnodes);
                    }
                    T1 d01 = nodes[edgeNodes[0]].getDistance(nodes[edgeNodes[1]]);
                    T1 w0 = curr_pt.getDistance(nodes[edgeNodes[1]]) / d01;
                    T1 w1 = curr_pt.getDistance(nodes[edgeNodes[0]]) / d01;
                    T1 curr_t = nodes[edgeNodes[0]].getTT(threadNo)*w0 + nodes[edgeNodes[1]].getTT(threadNo)*w1;
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(2);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[edgeNodes[0]]);
                    ref_pt[1] = &(nodes[edgeNodes[1]]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                bool foundIntersection = false;
                for (size_t n=0; n<cells.size(); ++n ) {
                    
                    cellNo = cells[n];
                    
                    // there are 2 faces that might be intersected
                    std::array<T2,2> edgeNodes2;
                    size_t n2=0;
                    for ( auto nn=neighbors[cellNo].begin(); nn!= neighbors[cellNo].end(); ++nn ) {
                        if ( *nn!=edgeNodes[0] && *nn!=edgeNodes[1] && nodes[*nn].isPrimary() ) {
                            edgeNodes2[n2++] = *nn;
                        }
                    }
                    
                    sxyz<T1> pt_i;
                    T2 itmpNode;
                    foundIntersection = intersectVecTriangle(curr_pt, g,
                                                             edgeNodes[0],
                                                             edgeNodes2[0],
                                                             edgeNodes2[1], pt_i);
                    itmpNode = edgeNodes[0];
                    if ( !foundIntersection ) {
                        foundIntersection = intersectVecTriangle(curr_pt, g,
                                                                 edgeNodes[1],
                                                                 edgeNodes2[0],
                                                                 edgeNodes2[1], pt_i);
                        itmpNode = edgeNodes[1];
                    }
                    if ( !foundIntersection ) {
                        continue;
                    }
                    
                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, neighbors[cellNo],
                                                        {itmpNode, edgeNodes2[0], edgeNodes2[1]},
                                                        onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    tt += 0.5*(s1 + s2) * r_tmp.back().getDistance( curr_pt );
//                    std::cout << curr_pt << '\t' << tt << '\t' << s1 << '\t' << s2 << '\n';

                    s1 = s2;
                    r_tmp.push_back( curr_pt );
                    
                    if ( break_flag ) break;
                    
                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }
                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath on edge failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    reachedTx = true;
                }
                
            } else if ( onFace ) { // on Face
#ifdef DEBUG_RP
                printRaypathData(curr_pt, g, onNode, onEdge, onFace, cellNo,
                                 nodeNo, edgeNodes, faceNodes);
#endif
                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( rp_method < 2 ) {
                    std::set<NODE*> nnodes;
                    getNeighborNodes(cellNo, nnodes);
                    T1 curr_t;
                    if ( r_tmp.size() <= 1 ) {
                        curr_t = Interpolator<T1>::trilinearTime(curr_pt,
                                                                 nodes[itmp[0]],
                                                                 nodes[itmp[1]],
                                                                 nodes[itmp[2]],
                                                                 nodes[itmp[3]],
                                                                 threadNo);
                    } else {
                        curr_t = Interpolator<T1>::bilinearTime(curr_pt,
                                                                nodes[faceNodes[0]],
                                                                nodes[faceNodes[1]],
                                                                nodes[faceNodes[2]],
                                                                threadNo);
                    }
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(3);
                    if ( r_tmp.size() <= 1 ) {
                        ref_pt[0] = &(nodes[itmp[0]]);
                        ref_pt[1] = &(nodes[itmp[1]]);
                        ref_pt[2] = &(nodes[itmp[2]]);
                        ref_pt.push_back( &(nodes[itmp[3]]) );
                    } else {
                        ref_pt[0] = &(nodes[faceNodes[0]]);
                        ref_pt[1] = &(nodes[faceNodes[1]]);
                        ref_pt[2] = &(nodes[faceNodes[2]]);
                    }
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                // there are 3 faces that might be intersected
                
                bool foundIntersection = false;
                for ( size_t n=0; n<4; ++n ) {
                    if ( ind[n] == faceNodes ) continue;
                    
                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);
                    
                    if ( !foundIntersection )
                        continue;
                    
                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    tt += 0.5*(s1 + s2) * r_tmp.back().getDistance( curr_pt );
//                    std::cout << curr_pt << '\t' << tt << '\t' << s1 << '\t' << s2 << '\n';
                    s1 = s2;
                    r_tmp.push_back( curr_pt );
                    
                    if ( break_flag ) break;
                    
                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }
                
                if ( foundIntersection == false ) {
                    
                    // we must be on an face with gradient pointing slightly outward tetrahedron
                    // return in other cell but keep gradient
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    
                    std::array<T2,4> itmp = getPrimary(cellNo);
                    ind[0] = { { itmp[0], itmp[1], itmp[2] } };
                    ind[1] = { { itmp[0], itmp[1], itmp[3] } };
                    ind[2] = { { itmp[0], itmp[2], itmp[3] } };
                    ind[3] = { { itmp[1], itmp[2], itmp[3] } };
                    
                    for ( size_t n=0; n<4; ++n )
                        std::sort( ind[n].begin(), ind[n].end() );
                    
                    for ( size_t n=0; n<4; ++n ) {
                        if ( ind[n] == faceNodes ) continue;
                        
                        sxyz<T1> pt_i;
                        foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                                 ind[n][1], ind[n][2],
                                                                 pt_i);
                        
                        if ( !foundIntersection ) {
                            continue;
                        }
                        curr_pt = pt_i;

                        bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                            nodeNo, onEdge, edgeNodes,
                                                            onFace, faceNodes);

                        s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                             faceNodes);

                        tt += 0.5*(s1 + s2) * r_tmp.back().getDistance( curr_pt );
//                        std::cout << curr_pt << '\t' << tt << '\t' << s1 << '\t' << s2 << '\n';

                        s1 = s2;
                        r_tmp.push_back( curr_pt );
                        
                        if ( break_flag ) break;
                        
                        // find next cell
                        cellNo = findAdjacentCell2(faceNodes, cellNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                            << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                            r_tmp.resize(1);
                            r_tmp[0] = Rx;
                            reachedTx = true;
                        }
                        break;
                    }
                }
                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath on face failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    reachedTx = true;
                }
            } else { // at Rx, somewhere in a tetrahedron
#ifdef DEBUG_RP
                printRaypathData(curr_pt, g, onNode, onEdge, onFace, cellNo,
                                 nodeNo, edgeNodes, faceNodes);
#endif
                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( rp_method < 2 ) {
                    std::set<NODE*> nnodes;
                    getNeighborNodes(cellNo, nnodes);
                    T1 curr_t = Interpolator<T1>::trilinearTime(curr_pt,
                                                                nodes[itmp[0]],
                                                                nodes[itmp[1]],
                                                                nodes[itmp[2]],
                                                                nodes[itmp[3]],
                                                                threadNo);
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(4);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[itmp[0]]);
                    ref_pt[1] = &(nodes[itmp[1]]);
                    ref_pt[2] = &(nodes[itmp[2]]);
                    ref_pt[3] = &(nodes[itmp[3]]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                // there are 3 faces that might be intersected
                
                bool foundIntersection = false;
                for ( size_t n=0; n<4; ++n ) {
                    
                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);
                    
                    if ( !foundIntersection )
                        continue;
                    
                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    tt += 0.5*(s1 + s2) * r_tmp.back().getDistance( curr_pt );
//                    std::cout << curr_pt << '\t' << tt << '\t' << s1 << '\t' << s2 << '\n';

                    s1 = s2;
                    r_tmp.push_back( curr_pt );
                    
                    if ( break_flag ) break;
                    
                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }
                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath within cell failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_data.resize(1);
                    r_data[0] = Rx;
                    reachedTx = true;
                }
            }
            
            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( curr_pt.getDistance( Tx[nt] ) < minDist ) {
                        tt += t0[nt];
                        reachedTx = true;
                        break;
                    }
                }
            } else {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        for ( auto nc=nodes[txNode[nt]].getOwners().begin();
                             nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                            if ( cellNo == *nc ) {
                                std::array<T2,4> itmp = getPrimary(cellNo);
                                
                                if ( interpVel )
                                    s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                nodes[itmp[0]],
                                                                                nodes[itmp[1]],
                                                                                nodes[itmp[2]],
                                                                                nodes[itmp[3]]);
                                else
                                    s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                             nodes[itmp[0]],
                                                                             nodes[itmp[1]],
                                                                             nodes[itmp[2]],
                                                                             nodes[itmp[3]]);

                                tt += t0[nt] + 0.5*(s1 + s2) * r_tmp.back().getDistance( Tx[nt] );
                                r_tmp.push_back( Tx[nt] );
                                reachedTx = true;
                                break;
                            }
                        }
                    } else {
                        if ( cellNo == txCell[nt] ) {
                            std::array<T2,4> itmp = getPrimary(cellNo);
                            if ( interpVel )
                                s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                            nodes[itmp[0]],
                                                                            nodes[itmp[1]],
                                                                            nodes[itmp[2]],
                                                                            nodes[itmp[3]]);
                            else
                                s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                         nodes[itmp[0]],
                                                                         nodes[itmp[1]],
                                                                         nodes[itmp[2]],
                                                                         nodes[itmp[3]]);

                            tt += t0[nt] + 0.5*(s1 + s2) * r_tmp.back().getDistance( Tx[nt] );
                            r_tmp.push_back( Tx[nt] );
                            reachedTx = true;
                            break;
                        } else {
                            for ( size_t nn=0; nn<txNeighborCells[nt].size(); ++nn ) {
                                if ( cellNo == txNeighborCells[nt][nn] ) {
                                    std::array<T2,4> itmp = getPrimary(txCell[nt]);
                                    
                                    bool found = false;
                                    for ( T2 i=0; i<4; ++i ) {
                                        T2 iA = itmp[(i+1)%4];
                                        T2 iB = itmp[(i+2)%4];
                                        T2 iC = itmp[(i+3)%4];
                                        if (nodes[iA].getDistance(curr_pt)+nodes[iB].getDistance(curr_pt)-nodes[iB].getDistance(nodes[iA])<minDist*minDist||
                                            nodes[iA].getDistance(curr_pt)+nodes[iC].getDistance(curr_pt)-nodes[iC].getDistance(nodes[iA])<minDist*minDist||
                                            nodes[iC].getDistance(curr_pt)+nodes[iB].getDistance(curr_pt)-nodes[iB].getDistance(nodes[iC])<minDist*minDist){
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found ) {
                                        itmp = getPrimary(cellNo);
                                        if ( interpVel )
                                            s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                        nodes[itmp[0]],
                                                                                        nodes[itmp[1]],
                                                                                        nodes[itmp[2]],
                                                                                        nodes[itmp[3]]);
                                        else
                                            s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                                     nodes[itmp[0]],
                                                                                     nodes[itmp[1]],
                                                                                     nodes[itmp[2]],
                                                                                     nodes[itmp[3]]);
                                        
                                        tt += t0[nt] + 0.5*(s1 + s2) * r_tmp.back().getDistance( Tx[nt] );
                                        r_tmp.push_back( Tx[nt] );
                                        reachedTx = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if ( reachedTx ) break;
                }
            }
        }
        // for inversion, the order should be from Tx to Rx, so we reorder...
        size_t npts = r_tmp.size();
        r_data.resize( npts );
        for ( size_t nn=0; nn<npts; ++nn ) {
            r_data[nn] = r_tmp[ npts-1-nn ];
        }
        
        delete grad3d;
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                          const sxyz<T1> &Rx,
                                          std::vector<sxyz<T1>> &r_data,
                                          std::vector<sijv<T1>>& m_data,
                                          const size_t RxNo,
                                          const size_t threadNo) const {
        
        T1 minDist = small;
        std::vector<sxyz<T1>> r_tmp;
        r_tmp.emplace_back( Rx );
        
        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                return;
            }
        }
        
        std::vector<bool> txOnNode( Tx.size(), false );
        std::vector<T2> txNode( Tx.size() );
        std::vector<T2> txCell( Tx.size() );
        std::vector<std::vector<T2>> txNeighborCells( Tx.size() );
        for ( size_t nt=0; nt<Tx.size(); ++nt ) {
            for ( T2 nn=0; nn<nodes.size(); ++nn ) {
                if ( nodes[nn] == Tx[nt] ) {
                    txOnNode[nt] = true;
                    txNode[nt] = nn;
                    break;
                }
            }
        }
        for ( size_t nt=0; nt<Tx.size(); ++nt ) {
            if ( !txOnNode[nt] ) {
                txCell[nt] = getCellNo( Tx[nt] );
                
                std::array<T2,4> itmp = getPrimary(txCell[nt]);
                // find adjacent cells
                T2 ind[6][2] = {
                    {itmp[0], itmp[1]},
                    {itmp[0], itmp[2]},
                    {itmp[0], itmp[3]},
                    {itmp[1], itmp[2]},
                    {itmp[1], itmp[3]},
                    {itmp[2], itmp[3]} };
                
                for ( size_t nedge=0; nedge<6; ++nedge ) {
                    for ( auto nc0=nodes[ind[nedge][0]].getOwners().begin(); nc0!=nodes[ind[nedge][0]].getOwners().end(); ++nc0 ) {
                        if ( std::find(nodes[ind[nedge][1]].getOwners().begin(), nodes[ind[nedge][1]].getOwners().end(), *nc0)!=nodes[ind[nedge][1]].getOwners().end() )
                            txNeighborCells[nt].push_back( *nc0 );
                    }
                }
            }
        }
        
        T2 cellNo, nodeNo, nodeNoPrev;
        sxyz<T1> curr_pt( Rx ), mid_pt, prev_pt( Rx );
        sijv<T1> m;
        m.i = RxNo;
        
        bool onNode = false;
        bool onEdge = false;
        bool onFace = false;
        bool onNodePrev = false;
        bool onEdgePrev = false;
        bool onFacePrev = false;
        std::array<T2,2> edgeNodes, edgeNodesPrev;
        std::array<T2,3> faceNodes, faceNodesPrev;
        Grad3D<T1,NODE>* grad3d = nullptr;
        if ( rp_method == 0 ) {
            grad3d = new Grad3D_ls_fo<T1,NODE>();
        } else if ( rp_method == 1 ) {
            grad3d = new Grad3D_ls_so<T1,NODE>();
        }
        bool reachedTx = false;
        
        for ( T2 nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == curr_pt ) {
                nodeNo = nn;
                onNode = true;
                break;
            }
        }
        if ( !onNode ) {
            cellNo = getCellNo( curr_pt );
            
            std::array<T2,4> itmp = getPrimary(cellNo);
            T2 ind[6][2] = {
                {itmp[0], itmp[1]},
                {itmp[0], itmp[2]},
                {itmp[0], itmp[3]},
                {itmp[1], itmp[2]},
                {itmp[1], itmp[3]},
                {itmp[2], itmp[3]} };
            
            for ( size_t n=0; n<6; ++n ) {
                if ( areCollinear(curr_pt, ind[n][0], ind[n][1]) ) {
                    onEdge = true;
                    edgeNodes[0] = ind[n][0];
                    edgeNodes[1] = ind[n][1];
                    break;
                }
            }
            if ( !onEdge ) {
                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                
                for ( size_t n=0; n<4; ++n ) {
                    if ( areCoplanar(curr_pt, ind[n][0], ind[n][1], ind[n][2]) ) {
                        onFace = true;
                        //  faceNodes shoud not be assigned, face was not intersected
                        break;
                    }
                }
            }
        }
        
        for ( auto nt=0; nt<txCell.size(); ++nt ) {
            if ( cellNo == txCell[nt] ) {
                r_tmp.push_back( Tx[nt] );
                reachedTx = true;
                break;
            }
        }
        
        T1 s, ds;
        while ( reachedTx == false ) {
            
            if ( onNode ) {
                
                // find cells common to edge
                std::set<NODE*> nnodes;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                    getNeighborNodes(*nc, nnodes);
                }
                
                // compute gradient with nodes from all common cells
                sxyz<T1> g = grad3d->compute(curr_pt, nodes[nodeNo].getTT(threadNo), nnodes, threadNo);
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                nodeNoPrev = nodeNo;
                onNodePrev = onNode;
                
                // find cell for which gradient intersect opposing face
                bool foundIntersection = false;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                    
                    std::array<T2,3> nb;
                    size_t n=0;
                    for (auto nn=neighbors[*nc].begin(); nn!=neighbors[*nc].end(); ++nn ) {
                        if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                            nb[n++] = *nn;
                        }
                    }
                    std::sort(nb.begin(), nb.end());
                    
                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle( nodeNo, g, nb[0], nb[1], nb[2], pt_i);
                    if ( !foundIntersection ) {
                        continue;
                    }
                    
                    prev_pt = curr_pt;
                    curr_pt = pt_i;
                    r_tmp.push_back( curr_pt );
                    
                    if ( r_tmp.size() > 1 ) {
                        // compute terms of matrix M
                        mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                        s = computeSlowness(mid_pt);
                        s *= s;
                        ds = curr_pt.getDistance( prev_pt );
                    }
                    
                    bool break_flag=false;
                    for ( n=0; n<3; ++n ) {
                        if ( nodes[ nb[n] ].getDistance( curr_pt ) < small ) {
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            nodeNo = nb[n];
                            onNode = true;
                            onEdge = false;
                            onFace = false;
                            
                            if ( r_tmp.size() > 1) {
                                std::set<T2> allNodes;
                                allNodes.insert(nodeNoPrev);
                                allNodes.insert(nodeNo);
                                
                                update_m_data(m_data, m, allNodes, mid_pt, s,  ds);
                            }
                            
                            break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;
                    
                    for ( size_t n1=0; n1<3; ++n1 ) {
                        size_t n2 = (n1+1)%3;
                        if ( areCollinear(curr_pt, nb[n1], nb[n2]) ) {
                            edgeNodesPrev = edgeNodes;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            edgeNodes[0] = nb[n1];
                            edgeNodes[1] = nb[n2];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            
                            if ( r_tmp.size() > 1) {
                                std::set<T2> allNodes;
                                allNodes.insert(nodeNoPrev);
                                allNodes.insert( edgeNodes[0] );
                                allNodes.insert( edgeNodes[1] );
                                
                                update_m_data(m_data, m, allNodes, mid_pt, s,  ds);
                            }
                            
                            break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;
                    
                    onEdgePrev = onEdge;
                    onFacePrev = onFace;
                    onNode = false;
                    onEdge = false;
                    onFace = true;
                    
                    faceNodesPrev = faceNodes;
                    faceNodes = nb;
                    
                    if ( r_tmp.size() > 1) {
                        std::set<T2> allNodes;
                        allNodes.insert(nodeNoPrev);
                        allNodes.insert( faceNodes[0] );
                        allNodes.insert( faceNodes[1] );
                        allNodes.insert( faceNodes[2] );
                        
                        update_m_data(m_data, m, allNodes, mid_pt, s,  ds);
                    }
                    
                    // find next cell
                    cellNo = findAdjacentCell1(faceNodes, nodeNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }
                
                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath on node failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    reachedTx = true;
                }
                
            } else if ( onEdge ) {
                
                // find cells common to edge
                std::vector<T2> cells;
                std::set<NODE*> nnodes;
                for ( auto nc0=nodes[edgeNodes[0]].getOwners().begin(); nc0!=nodes[edgeNodes[0]].getOwners().end(); ++nc0 ) {
                    if ( std::find(nodes[edgeNodes[1]].getOwners().begin(), nodes[edgeNodes[1]].getOwners().end(), *nc0)!=nodes[edgeNodes[1]].getOwners().end() ) {
                        cells.push_back( *nc0 );
                        getNeighborNodes(*nc0, nnodes);
                    }
                }
                T1 d01 = nodes[edgeNodes[0]].getDistance(nodes[edgeNodes[1]]);
                T1 w0 = curr_pt.getDistance(nodes[edgeNodes[1]]) / d01;
                T1 w1 = curr_pt.getDistance(nodes[edgeNodes[0]]) / d01;
                T1 curr_t = nodes[edgeNodes[0]].getTT(threadNo)*w0 + nodes[edgeNodes[1]].getTT(threadNo)*w1;
                sxyz<T1> g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                edgeNodesPrev = edgeNodes;
                onEdgePrev = onEdge;
                
                bool foundIntersection = false;
                for (size_t n=0; n<cells.size(); ++n ) {
                    
                    cellNo = cells[n];
                    
                    // there are 2 faces that might be intersected
                    std::array<T2,2> edgeNodes2;
                    size_t n2=0;
                    for ( auto nn=neighbors[cellNo].begin(); nn!= neighbors[cellNo].end(); ++nn ) {
                        if ( *nn!=edgeNodes[0] && *nn!=edgeNodes[1] && nodes[*nn].isPrimary() ) {
                            edgeNodes2[n2++] = *nn;
                        }
                    }
                    
                    sxyz<T1> pt_i;
                    T2 itmpNode;
                    foundIntersection = intersectVecTriangle(curr_pt, g,
                                                             edgeNodes[0],
                                                             edgeNodes2[0],
                                                             edgeNodes2[1], pt_i);
                    itmpNode = edgeNodes[0];
                    if ( !foundIntersection ) {
                        foundIntersection = intersectVecTriangle(curr_pt, g,
                                                                 edgeNodes[1],
                                                                 edgeNodes2[0],
                                                                 edgeNodes2[1], pt_i);
                        itmpNode = edgeNodes[1];
                    }
                    if ( !foundIntersection ) {
                        continue;
                    }
                    
                    prev_pt = curr_pt;
                    curr_pt = pt_i;
                    r_tmp.push_back( curr_pt );
                    
                    if (r_tmp.size() > 1 ) {
                        // compute terms of matrix M
                        mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                        s = computeSlowness(mid_pt);
                        s *= s;
                        ds = curr_pt.getDistance( prev_pt );
                    }
                    
                    bool break_flag = false;
                    for ( size_t n2=0; n2<4; ++n2 ) {
                        if ( nodes[ neighbors[cellNo][n2] ].getDistance( curr_pt ) < small ) {
                            onNodePrev = onNode;
                            onFacePrev = onFace;
                            
                            nodeNo = neighbors[cellNo][n2];
                            onNode = true;
                            onEdge = false;
                            onFace = false;
                            
                            if ( r_tmp.size() > 1) {
                                std::set<T2> allNodes;
                                allNodes.insert(nodeNo);
                                allNodes.insert( edgeNodesPrev[0] );
                                allNodes.insert( edgeNodesPrev[1] );
                                
                                update_m_data(m_data, m, allNodes, mid_pt, s,  ds);
                            }
                            
                            break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;
                    
                    if ( areCollinear(curr_pt, itmpNode, edgeNodes2[0]) ) {
                        onNodePrev = onNode;
                        onFacePrev = onFace;
                        
                        edgeNodes[0] = itmpNode;
                        edgeNodes[1] = edgeNodes2[0];
                        onNode = false;
                        onEdge = true;
                        onFace = false;
                        
                        if ( r_tmp.size() > 1) {
                            std::set<T2> allNodes;
                            allNodes.insert( edgeNodesPrev[0] );
                            allNodes.insert( edgeNodesPrev[1] );
                            allNodes.insert( edgeNodes[0] );
                            allNodes.insert( edgeNodes[1] );
                            
                            update_m_data(m_data, m, allNodes, mid_pt, s,  ds);
                        }
                        
                        break_flag = true;
                        break;
                    } else if ( areCollinear(curr_pt, itmpNode, edgeNodes2[1]) ) {
                        onNodePrev = onNode;
                        onFacePrev = onFace;
                        
                        edgeNodes[0] = itmpNode;
                        edgeNodes[1] = edgeNodes2[1];
                        onNode = false;
                        onEdge = true;
                        onFace = false;
                        
                        if ( r_tmp.size() > 1) {
                            std::set<T2> allNodes;
                            allNodes.insert( edgeNodesPrev[0] );
                            allNodes.insert( edgeNodesPrev[1] );
                            allNodes.insert( edgeNodes[0] );
                            allNodes.insert( edgeNodes[1] );
                            
                            update_m_data(m_data, m, allNodes, mid_pt, s,  ds);
                        }
                        
                        break_flag = true;
                        break;
                    } else if ( areCollinear(curr_pt, edgeNodes2[0], edgeNodes2[1]) ) {
                        onNodePrev = onNode;
                        onFacePrev = onFace;
                        
                        edgeNodes[0] = edgeNodes2[0];
                        edgeNodes[1] = edgeNodes2[1];
                        onNode = false;
                        onEdge = true;
                        onFace = false;
                        
                        if ( r_tmp.size() > 1) {
                            std::set<T2> allNodes;
                            allNodes.insert( edgeNodesPrev[0] );
                            allNodes.insert( edgeNodesPrev[1] );
                            allNodes.insert( edgeNodes[0] );
                            allNodes.insert( edgeNodes[1] );
                            
                            update_m_data(m_data, m, allNodes, mid_pt, s,  ds);
                        }
                        
                        break_flag = true;
                        break;
                    }
                    if ( break_flag ) break;
                    
                    onNodePrev = onNode;
                    onFacePrev = onFace;
                    onNode = false;
                    onEdge = false;
                    onFace = true;
                    
                    faceNodesPrev = faceNodes;
                    faceNodes[0] = itmpNode;
                    faceNodes[1] = edgeNodes2[0];
                    faceNodes[2] = edgeNodes2[1];
                    std::sort(faceNodes.begin(), faceNodes.end());
                    
                    if ( r_tmp.size() > 1) {
                        std::set<T2> allNodes;
                        allNodes.insert( edgeNodesPrev[0] );
                        allNodes.insert( edgeNodesPrev[1] );
                        allNodes.insert( faceNodes[0] );
                        allNodes.insert( faceNodes[1] );
                        allNodes.insert( faceNodes[2] );
                        
                        update_m_data(m_data, m, allNodes, mid_pt, s,  ds);
                    }
                    
                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }
                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath on edge failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    reachedTx = true;
                }
                
            } else { // on Face
                
                std::set<NODE*> nnodes;
                getNeighborNodes(cellNo, nnodes);
                
                std::array<T2,4> itmp = getPrimary(cellNo);
                T1 curr_t = Interpolator<T1>::trilinearTime(curr_pt,
                                                            nodes[itmp[0]],
                                                            nodes[itmp[1]],
                                                            nodes[itmp[2]],
                                                            nodes[itmp[3]],
                                                            threadNo);
                sxyz<T1> g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                // there are 3 faces that might be intersected
                
                onFacePrev = onFace;
                faceNodesPrev = faceNodes;
                
                bool foundIntersection = false;
                for ( size_t n=0; n<4; ++n ) {
                    if ( ind[n] == faceNodes ) continue;
                    
                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);
                    
                    if ( !foundIntersection )
                        continue;
                    
                    prev_pt = curr_pt;
                    curr_pt = pt_i;
                    r_tmp.push_back( curr_pt );
                    
                    if (r_tmp.size() > 1 ) {
                        // compute terms of matrix M
                        mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                        s = computeSlowness(mid_pt);
                        s *= s;
                        ds = curr_pt.getDistance( prev_pt );
                    }
                    
                    bool break_flag = false;
                    for ( size_t n2=0; n2<3; ++n2 ) {
                        if ( nodes[ ind[n][n2] ].getDistance( curr_pt ) < small ) {
                            nodeNoPrev = nodeNo;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            
                            nodeNo = ind[n][n2];
                            onNode = true;
                            onEdge = false;
                            onFace = false;
                            
                            if ( r_tmp.size() > 1) {
                                std::set<T2> allNodes;
                                allNodes.insert(nodeNo);
                                allNodes.insert( faceNodesPrev[0] );
                                allNodes.insert( faceNodesPrev[1] );
                                allNodes.insert( faceNodesPrev[2] );
                                
                                update_m_data(m_data, m, allNodes, mid_pt, s,  ds);
                            }
                            
                            break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;
                    
                    for ( size_t n1=0; n1<3; ++n1 ) {
                        size_t n2 = (n1+1)%3;
                        if ( areCollinear(curr_pt, ind[n][n1], ind[n][n2]) ) {
                            edgeNodesPrev = edgeNodes;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            
                            edgeNodes[0] = ind[n][n1];
                            edgeNodes[1] = ind[n][n2];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            
                            if ( r_tmp.size() > 1) {
                                std::set<T2> allNodes;
                                allNodes.insert( edgeNodes[0] );
                                allNodes.insert( edgeNodes[1] );
                                allNodes.insert( faceNodesPrev[0] );
                                allNodes.insert( faceNodesPrev[1] );
                                allNodes.insert( faceNodesPrev[2] );
                                
                                update_m_data(m_data, m, allNodes, mid_pt, s,  ds);
                            }
                            
                            break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;
                    
                    onNodePrev = onNode;
                    onEdgePrev = onEdge;
                    onNode = false;
                    onEdge = false;
                    onFace = true;
                    
                    faceNodes = ind[n];
                    
                    if ( r_tmp.size() > 1) {
                        std::set<T2> allNodes;
                        allNodes.insert( faceNodesPrev[0] );
                        allNodes.insert( faceNodesPrev[1] );
                        allNodes.insert( faceNodesPrev[2] );
                        allNodes.insert( faceNodes[0] );
                        allNodes.insert( faceNodes[1] );
                        allNodes.insert( faceNodes[2] );
                        
                        update_m_data(m_data, m, allNodes, mid_pt, s,  ds);
                    }
                    
                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }
                
                if ( foundIntersection == false ) {
                    
                    // we must be on an face with gradient pointing slightly outward tetrahedron
                    // return in other cell but keep gradient
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    
                    std::array<T2,4> itmp = getPrimary(cellNo);
                    ind[0] = { { itmp[0], itmp[1], itmp[2] } };
                    ind[1] = { { itmp[0], itmp[1], itmp[3] } };
                    ind[2] = { { itmp[0], itmp[2], itmp[3] } };
                    ind[3] = { { itmp[1], itmp[2], itmp[3] } };
                    
                    for ( size_t n=0; n<4; ++n )
                        std::sort( ind[n].begin(), ind[n].end() );
                    
                    for ( size_t n=0; n<4; ++n ) {
                        if ( ind[n] == faceNodes ) continue;
                        
                        sxyz<T1> pt_i;
                        foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                                 ind[n][1], ind[n][2],
                                                                 pt_i);
                        
                        if ( !foundIntersection ) {
                            continue;
                        }
                        prev_pt = curr_pt;
                        curr_pt = pt_i;
                        r_tmp.push_back( curr_pt );
                        
                        if (r_tmp.size() > 1 ) {
                            // compute terms of matrix M
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            s = computeSlowness(mid_pt);
                            s *= s;
                            ds = curr_pt.getDistance( prev_pt );
                        }
                        
                        bool break_flag = false;
                        for ( size_t n2=0; n2<3; ++n2 ) {
                            if ( nodes[ ind[n][n2] ].getDistance( curr_pt ) < small ) {
                                nodeNoPrev = nodeNo;
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                
                                nodeNo = ind[n][n2];
                                onNode = true;
                                onEdge = false;
                                onFace = false;
                                
                                if ( r_tmp.size() > 1) {
                                    std::set<T2> allNodes;
                                    allNodes.insert(nodeNo);
                                    allNodes.insert( faceNodesPrev[0] );
                                    allNodes.insert( faceNodesPrev[1] );
                                    allNodes.insert( faceNodesPrev[2] );
                                    
                                    update_m_data(m_data, m, allNodes, mid_pt, s,  ds);
                                }
                                
                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;
                        
                        for ( size_t n1=0; n1<3; ++n1 ) {
                            size_t n2 = (n1+1)%3;
                            if ( areCollinear(curr_pt, ind[n][n1], ind[n][n2]) ) {
                                edgeNodesPrev = edgeNodes;
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                
                                edgeNodes[0] = ind[n][n1];
                                edgeNodes[1] = ind[n][n2];
                                onNode = false;
                                onEdge = true;
                                onFace = false;
                                
                                if ( r_tmp.size() > 1) {
                                    std::set<T2> allNodes;
                                    allNodes.insert( edgeNodes[0] );
                                    allNodes.insert( edgeNodes[1] );
                                    allNodes.insert( faceNodesPrev[0] );
                                    allNodes.insert( faceNodesPrev[1] );
                                    allNodes.insert( faceNodesPrev[2] );
                                    
                                    update_m_data(m_data, m, allNodes, mid_pt, s,  ds);
                                }
                                
                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;
                        
                        onNodePrev = onNode;
                        onEdgePrev = onEdge;
                        onNode = false;
                        onEdge = false;
                        onFace = true;
                        
                        faceNodes = ind[n];
                        
                        if ( r_tmp.size() > 1) {
                            std::set<T2> allNodes;
                            allNodes.insert( faceNodesPrev[0] );
                            allNodes.insert( faceNodesPrev[1] );
                            allNodes.insert( faceNodesPrev[2] );
                            allNodes.insert( faceNodes[0] );
                            allNodes.insert( faceNodes[1] );
                            allNodes.insert( faceNodes[2] );
                            
                            update_m_data(m_data, m, allNodes, mid_pt, s,  ds);
                        }
                        
                        // find next cell
                        cellNo = findAdjacentCell2(faceNodes, cellNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                            << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                            r_tmp.resize(1);
                            r_tmp[0] = Rx;
                            reachedTx = true;
                        }
                        break;
                    }
                }
                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath on face failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    reachedTx = true;
                }
            }
            
            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( curr_pt.getDistance( Tx[nt] ) < minDist ) {
                        reachedTx = true;
                        break;
                    }
                }
            } else {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        for ( auto nc=nodes[txNode[nt]].getOwners().begin();
                             nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                            if ( cellNo == *nc ) {
                                r_tmp.push_back( Tx[nt] );
                                reachedTx = true;
                                break;
                            }
                        }
                    } else {
                        if ( cellNo == txCell[nt] ) {
                            r_tmp.push_back( Tx[nt] );
                            reachedTx = true;
                        } else {
                            for ( size_t nn=0; nn<txNeighborCells[nt].size(); ++nn ) {
                                if ( cellNo == txNeighborCells[nt][nn] ) {
                                    std::array<T2,4> itmp = getPrimary(txCell[nt]);
                                    
                                    bool found = false;
                                    for ( T2 i=0; i<4; ++i ) {
                                        T2 iA = itmp[(i+1)%4];
                                        T2 iB = itmp[(i+2)%4];
                                        T2 iC = itmp[(i+3)%4];
                                        if (nodes[iA].getDistance(curr_pt)+nodes[iB].getDistance(curr_pt)-nodes[iB].getDistance(nodes[iA])<minDist*minDist||
                                            nodes[iA].getDistance(curr_pt)+nodes[iC].getDistance(curr_pt)-nodes[iC].getDistance(nodes[iA])<minDist*minDist||
                                            nodes[iC].getDistance(curr_pt)+nodes[iB].getDistance(curr_pt)-nodes[iB].getDistance(nodes[iC])<minDist*minDist){
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found ) {
                                        r_tmp.push_back( Tx[nt] );
                                        reachedTx = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if ( reachedTx ) break;
                }
            }
        }
        // for inversion, the order should be from Tx to Rx, so we reorder...
        size_t npts = r_tmp.size();
        r_data.resize( npts );
        for ( size_t nn=0; nn<npts; ++nn ) {
            r_data[nn] = r_tmp[ npts-1-nn ];
        }
        
        delete grad3d;
    }

    
    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::getRaypath_blti(const std::vector<sxyz<T1>>& Tx,
                                               const std::vector<T1>& t0,
                                               const sxyz<T1> &Rx,
                                               std::vector<sxyz<T1>> &r_data,
                                               T1& tt,
                                               const size_t threadNo) const {
        
        T1 minDist = small;
        std::vector<sxyz<T1>> r_tmp;
        r_tmp.emplace_back( Rx );
        tt = 0.0;
        T1 s1, s2;
        
        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                return;
            }
        }
        
        std::vector<bool> txOnNode( Tx.size(), false );
        std::vector<T2> txNode( Tx.size() );
        std::vector<T2> txCell( Tx.size() );
        std::vector<std::vector<T2>> txNeighborCells( Tx.size() );
        for ( size_t nt=0; nt<Tx.size(); ++nt ) {
            for ( T2 nn=0; nn<nodes.size(); ++nn ) {
                if ( nodes[nn] == Tx[nt] ) {
                    txOnNode[nt] = true;
                    txNode[nt] = nn;
                    break;
                }
            }
        }
        for ( size_t nt=0; nt<Tx.size(); ++nt ) {
            if ( !txOnNode[nt] ) {
                txCell[nt] = getCellNo( Tx[nt] );
                
                std::array<T2,4> itmp = getPrimary(txCell[nt]);
                // find adjacent cells
                T2 ind[6][2] = {
                    {itmp[0], itmp[1]},
                    {itmp[0], itmp[2]},
                    {itmp[0], itmp[3]},
                    {itmp[1], itmp[2]},
                    {itmp[1], itmp[3]},
                    {itmp[2], itmp[3]} };
                
                for ( size_t nedge=0; nedge<6; ++nedge ) {
                    for ( auto nc0=nodes[ind[nedge][0]].getOwners().begin(); nc0!=nodes[ind[nedge][0]].getOwners().end(); ++nc0 ) {
                        if ( std::find(nodes[ind[nedge][1]].getOwners().begin(), nodes[ind[nedge][1]].getOwners().end(), *nc0)!=nodes[ind[nedge][1]].getOwners().end() )
                            txNeighborCells[nt].push_back( *nc0 );
                    }
                }
            }
        }
        
        T2 cellNo, nodeNo, nodeNoPrev;
        sxyz<T1> curr_pt( Rx ), prev_pt( Rx );
        bool onNode = false;
        bool onEdge = false;
        bool onFace = false;
        bool onNodePrev = false;
        bool onEdgePrev = false;
        bool onFacePrev = false;
        std::array<T2,2> edgeNodes, edgeNodesPrev;
        std::array<T2,3> faceNodes={{0,0,0}};
        std::array<T2,3> faceNodesPrev={{0,0,0}};
        bool reachedTx = false;
        
        for ( T2 nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == curr_pt ) {
                nodeNo = nn;
                onNode = true;
                s1 = nodes[nodeNo].getNodeSlowness();
                break;
            }
        }
        if ( !onNode ) {
            cellNo = getCellNo( curr_pt );
            
            std::array<T2,4> itmp = getPrimary(cellNo);
            T2 ind[6][2] = {
                {itmp[0], itmp[1]},
                {itmp[0], itmp[2]},
                {itmp[0], itmp[3]},
                {itmp[1], itmp[2]},
                {itmp[1], itmp[3]},
                {itmp[2], itmp[3]} };
            
            for ( size_t n=0; n<6; ++n ) {
                if ( areCollinear(curr_pt, ind[n][0], ind[n][1]) ) {
                    onEdge = true;
                    edgeNodes[0] = ind[n][0];
                    edgeNodes[1] = ind[n][1];
                    
                    if ( interpVel )
                        s1 = Interpolator<T1>::linearVel(r_tmp.back(),
                                                         nodes[edgeNodes[0]],
                                                         nodes[edgeNodes[1]]);
                    else
                        s1 = Interpolator<T1>::linear(r_tmp.back(),
                                                      nodes[edgeNodes[0]],
                                                      nodes[edgeNodes[1]]);

                    break;
                }
            }
            if ( !onEdge ) {
                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };

                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                
                for ( size_t n=0; n<4; ++n ) {
                    if ( areCoplanar(curr_pt, ind[n][0], ind[n][1], ind[n][2]) ) {
                        onFace = true;
                        faceNodes[0] = ind[n][0];
                        faceNodes[1] = ind[n][1];
                        faceNodes[2] = ind[n][2];
                        
                        if ( interpVel )
                            s1 = Interpolator<T1>::bilinearTriangleVel(curr_pt,
                                                                       nodes[faceNodes[0]],
                                                                       nodes[faceNodes[1]],
                                                                       nodes[faceNodes[2]]);
                        else
                            s1 = Interpolator<T1>::bilinearTriangle(curr_pt,
                                                                    nodes[faceNodes[0]],
                                                                    nodes[faceNodes[1]],
                                                                    nodes[faceNodes[2]]);
                        break;
                    }
                }
            }
        }
        T1 time=std::numeric_limits<T1>::max();
        if (!onEdge && ! onNode && ! onFace){
            onFace=true;
            std::array<T2,4> itmp = getPrimary(cellNo);
            time=Interpolator<T1>::trilinearTime(curr_pt,
                                                 nodes[itmp[0]],
                                                 nodes[itmp[1]],
                                                 nodes[itmp[2]],
                                                 nodes[itmp[3]], threadNo);
            
            if ( interpVel )
                s1 = Interpolator<T1>::trilinearTriangleVel(curr_pt,
                                                            nodes[itmp[0]],
                                                            nodes[itmp[1]],
                                                            nodes[itmp[2]],
                                                            nodes[itmp[3]]);
            else
                s1 = Interpolator<T1>::trilinearTriangle(curr_pt,
                                                         nodes[itmp[0]],
                                                         nodes[itmp[1]],
                                                         nodes[itmp[2]],
                                                         nodes[itmp[3]]);
        }
        for(auto nt=0;nt<txCell.size();++nt){
            if (getCellNo( Rx )==txCell[nt]){
                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( interpVel )
                    s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                nodes[itmp[0]],
                                                                nodes[itmp[1]],
                                                                nodes[itmp[2]],
                                                                nodes[itmp[3]]);
                else
                    s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                             nodes[itmp[0]],
                                                             nodes[itmp[1]],
                                                             nodes[itmp[2]],
                                                             nodes[itmp[3]]);
                
                tt += t0[nt] + 0.5*(s1 + s2) * r_tmp.back().getDistance( Tx[nt] );
                r_tmp.push_back(Tx[nt]);
                reachedTx=true;
                break;
            }
        }
        T2 N=0;
        while ( reachedTx == false && N<500) {
            ++N;
            //            if (N==29   )
            //               cout<<"stop";
            sxyz<T1> NodeSource;
            bool NearSource=false;
            for(size_t nt=0;nt<Tx.size();++nt){
                for(auto n=neighbors[txCell[nt]].begin();n!=neighbors[txCell[nt]].begin()+4;++n){
                    for(auto nc=nodes[*n].getOwners().begin();nc!=nodes[*n].getOwners().end();++nc){
                        if(*nc==cellNo){
                            NearSource=true;
                            NodeSource=Tx[nt];
                            break;
                        }
                    }
                    if (NearSource)
                        break;
                }
                if (curr_pt.getDistance(Tx[nt])<50.0*1.e-3){
                    NearSource=true;
                    NodeSource=Tx[nt];
                    break;
                }
                
            }
            
            if ( onNode ) {
                
                T1 t_i=std::numeric_limits<T1>::max();
                T1 Slow;
                sxyz<T1> pt_i;
                T2 newNode;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                    T2 cellNoi = *nc;
                    std::array<T2, 3>ind;
                    for (T2 i=0;i<4;++i){
                        if (neighbors[cellNoi][i]==nodeNo){
                            ind[0]=neighbors[cellNoi][(i+1)%4];
                            ind[1]=neighbors[cellNoi][(i+2)%4];
                            ind[2]=neighbors[cellNoi][(i+3)%4];
                            break;
                        }
                    }
                    std::sort(ind.begin(), ind.end());
                    if (NearSource){
                        std::array<T1,3> Barycenter;
                        if (blti_solver_around_source(NodeSource, curr_pt, ind, Barycenter)==true){
                            pt_i.x=Barycenter[0]*nodes[ind[0]].getX()+Barycenter[1]*nodes[ind[1]].getX()+Barycenter[2]*nodes[ind[2]].getX();
                            pt_i.y=Barycenter[0]*nodes[ind[0]].getY()+Barycenter[1]*nodes[ind[1]].getY()+Barycenter[2]*nodes[ind[2]].getY();
                            pt_i.z=Barycenter[0]*nodes[ind[0]].getZ()+Barycenter[1]*nodes[ind[1]].getZ()+Barycenter[2]*nodes[ind[2]].getZ();
                            if(NodeSource.getDistance(pt_i)>NodeSource.getDistance(curr_pt))
                                continue;
                            for(T2 ni=0;ni<3;++ni){
                                if(abs(1.0-Barycenter[ni])<minDist*minDist){
                                    onNodePrev = onNode;
                                    onEdgePrev = onEdge;
                                    onFacePrev = onFace;
                                    onNode = true;
                                    onEdge = false;
                                    onFace = false;
                                    nodeNoPrev=nodeNo;
                                    newNode=ind[ni];
                                    cellNo=cellNoi;
                                    break;
                                }
                                if(abs(Barycenter[ni])<minDist*minDist){
                                    onNodePrev = onNode;
                                    onEdgePrev = onEdge;
                                    onFacePrev = onFace;
                                    onNode = false;
                                    onEdge = true;
                                    onFace = false;
                                    nodeNoPrev=nodeNo;
                                    edgeNodes[0] = ind[(ni+1)%3];
                                    edgeNodes[1] = ind[(ni+2)%3];
                                    cellNo=cellNoi;
                                    break;
                                }
                            }
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            onNode = false;
                            onEdge = false;
                            onFace = true;
                            nodeNoPrev=nodeNo;
                            faceNodes=ind;
                            cellNo=cellNoi;
                            break;
                        }
                        continue;
                    }
                    T1 s=0.25*(nodes[neighbors[cellNoi][0]].getNodeSlowness()+
                               nodes[neighbors[cellNoi][1]].getNodeSlowness()+
                               nodes[neighbors[cellNoi][2]].getNodeSlowness()+
                               nodes[neighbors[cellNoi][3]].getNodeSlowness());
                    std::array<T2,2> indEdges[3] = {{{ ind[0],ind[1]}},{{ind[0],ind[2]}},{{ ind[1],ind[2]}}};
                    for (T2 i=0;i<4;++i){
                        if (neighbors[cellNoi][i]==nodeNo)
                            continue;
                        T1 t=nodes[neighbors[cellNoi][i]].getTT(threadNo);
                        if (t>time)
                            continue;
                        t+=s*nodes[neighbors[cellNoi][i]].getDistance(curr_pt);
                        if (t<t_i){
                            t_i=t;
                            Slow=s;
                            pt_i=sxyz<T1>(nodes[neighbors[cellNoi][i]]);
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            onNode = true;
                            onEdge = false;
                            onFace = false;
                            newNode=neighbors[cellNoi][i];
                            cellNo=cellNoi;
                        }
                    }
                    for(T2 n=0;n<3;++n){
                        
                        sxyz<T1> pt;
                        if(blti2D_raytrace(curr_pt, indEdges[n][0], indEdges[n][1], pt, threadNo, s) ==false) continue;
                        T1 dist0=pt.getDistance(nodes[indEdges[n][0]]);
                        T1 dist1=pt.getDistance(nodes[indEdges[n][1]]);
                        T1 t= (dist1*nodes[indEdges[n][0]].getTT(threadNo)+dist0*nodes[indEdges[n][1]].getTT(threadNo))/(dist0+dist1);
                        if (t>time)
                            continue;
                        t+=curr_pt.getDistance(pt)*s;
                        if (t<t_i ){
                            t_i=t;
                            Slow=s;
                            pt_i=pt;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            edgeNodes[0] = indEdges[n][0];
                            edgeNodes[1] = indEdges[n][1];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            nodeNoPrev=nodeNo;
                            cellNo=cellNoi;
                            
                        }
                        
                    }
                    sxyz<T1> pt;
                    if(blti_raytrace(curr_pt, ind, pt, threadNo, s)==false) continue;
                    T1 t=Interpolator<T1>::bilinearTime(pt, nodes[ind[0]], nodes[ind[1]], nodes[ind[2]], threadNo);
                    if (t>time)
                        continue;
                    t+=curr_pt.getDistance(pt)*s;
                    if (t<t_i){
                        t_i=t;
                        Slow=s;
                        pt_i=pt;
                        onNodePrev = onNode;
                        onEdgePrev = onEdge;
                        onFacePrev = onFace;
                        onNode = false;
                        onEdge = false;
                        onFace = true;
                        nodeNoPrev=nodeNo;
                        faceNodes = ind;
                        cellNo=cellNoi;
                    }
                    
                }
                
                prev_pt = curr_pt;
                curr_pt = pt_i;
                s2 = computeSlowness(curr_pt);
                
                tt += 0.5*(s1 + s2) * r_tmp.back().getDistance( curr_pt );
                s1 = s2;
                r_tmp.push_back( curr_pt );
                time=t_i-curr_pt.getDistance(prev_pt)*Slow;
                if (onNode){
                    nodeNoPrev=nodeNo;
                    nodeNo=newNode;
                }
                if (!NearSource){
                    if ( t_i==std::numeric_limits<T1>::max()) {
                        std::cout << "\n\nWarning: finding raypath on edge failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                }
                
            } else if ( onEdge ) {
                
                // find cells common to edge
                std::vector<T2> cells;
                T1 Slow;
                for ( auto nc0=nodes[edgeNodes[0]].getOwners().begin(); nc0!=nodes[edgeNodes[0]].getOwners().end(); ++nc0 ) {
                    if ( std::find(nodes[edgeNodes[1]].getOwners().begin(), nodes[edgeNodes[1]].getOwners().end(), *nc0)!=nodes[edgeNodes[1]].getOwners().end() ) {
                        cells.push_back( *nc0 );
                    }
                }
                std::array<T2,2> edgeNodestmp;
                T1 t_i=std::numeric_limits<T1>::max();
                sxyz<T1> pt_i;
                for (T2 nc=0; nc<cells.size(); ++nc ) {
                    T2 cellNoi = cells[nc];
                    // there are 2 faces that might be intersected
                    std::array<T2,2> edgeNodes2;
                    size_t n2=0;
                    for ( auto nn=neighbors[cellNoi].begin(); nn!= neighbors[cellNoi].begin()+4; ++nn ) {
                        if ( *nn!=edgeNodes[0] && *nn!=edgeNodes[1] ) {
                            edgeNodes2[n2++] = *nn;
                        }
                    }
                    std::array<T2,3> ind[2] = {
                        { { edgeNodes[0], edgeNodes2[0], edgeNodes2[1]}},
                        { { edgeNodes[1], edgeNodes2[0], edgeNodes2[1]}}};
                    std::sort( ind[0].begin(), ind[0].end());
                    std::sort( ind[1].begin(), ind[1].end());
                    std::array<T2,2> indEdges[5] = {
                        {{ edgeNodes[0],edgeNodes2[0]}},
                        {{ edgeNodes[0],edgeNodes2[1]}},
                        {{ edgeNodes[1],edgeNodes2[0]}},
                        {{ edgeNodes[1],edgeNodes2[1]}},
                        {{ edgeNodes2[0],edgeNodes2[1]}}};
                    T1 s=0.25*(nodes[edgeNodes[0]].getNodeSlowness()+
                               nodes[edgeNodes[1]].getNodeSlowness()+
                               nodes[edgeNodes2[0]].getNodeSlowness()+
                               nodes[edgeNodes2[1]].getNodeSlowness());
                    bool flag=false;
                    for ( size_t n=0; n<2; ++n ) {
                        if (NearSource){
                            std::array<T1,3> Barycenter;
                            if (blti_solver_around_source(NodeSource, curr_pt, ind[n], Barycenter)==true){
                                pt_i.x=Barycenter[0]*nodes[ind[n][0]].getX()+Barycenter[1]*nodes[ind[n][1]].getX()+Barycenter[2]*nodes[ind[n][2]].getX();
                                pt_i.y=Barycenter[0]*nodes[ind[n][0]].getY()+Barycenter[1]*nodes[ind[n][1]].getY()+Barycenter[2]*nodes[ind[n][2]].getY();
                                pt_i.z=Barycenter[0]*nodes[ind[n][0]].getZ()+Barycenter[1]*nodes[ind[n][1]].getZ()+Barycenter[2]*nodes[ind[n][2]].getZ();
                                if(NodeSource.getDistance(pt_i)>NodeSource.getDistance(curr_pt))
                                    continue;
                                for(T2 ni=0;ni<3;++ni){
                                    if(abs(1.0-Barycenter[ni])<minDist*minDist){
                                        onNodePrev = onNode;
                                        onEdgePrev = onEdge;
                                        onFacePrev = onFace;
                                        onNode = true;
                                        onEdge = false;
                                        onFace = false;
                                        edgeNodesPrev=edgeNodes;
                                        nodeNo=ind[n][ni];
                                        cellNo=cellNoi;
                                        flag=true;
                                        break;
                                    }
                                    if(abs(Barycenter[ni])<minDist*minDist){
                                        onNodePrev = onNode;
                                        onEdgePrev = onEdge;
                                        onFacePrev = onFace;
                                        onNode = false;
                                        onEdge = true;
                                        onFace = false;
                                        edgeNodes[0] = ind[n][(ni+1)%3];
                                        edgeNodes[1] = ind[n][(ni+2)%3];
                                        cellNo=cellNoi;
                                        flag=true;
                                        break;
                                    }
                                }
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;
                                onNode = false;
                                onEdge = false;
                                onFace = true;
                                edgeNodesPrev=edgeNodes;
                                faceNodes=ind[n];
                                cellNo=cellNoi;
                                flag=true;
                                break;
                            }
                            continue;
                        }
                        sxyz<T1> pt;
                        if(blti_raytrace(curr_pt, ind[n], pt, threadNo, s)==false) continue;
                        T1 t=Interpolator<T1>::bilinearTime(pt, nodes[ind[n][0]], nodes[ind[n][1]], nodes[ind[n][2]], threadNo);
                        if (t>time)
                            continue;
                        t+=curr_pt.getDistance(pt)*s;
                        if (t<t_i){
                            t_i=t;
                            Slow=s;
                            pt_i=pt;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            onNode = false;
                            onEdge = false;
                            onFace = true;
                            edgeNodesPrev = edgeNodes;
                            faceNodes = ind[n];
                            cellNo=cellNoi;
                        }
                    }
                    if (flag)
                        break;
                    if (NearSource)
                        continue;
                    for(T2 n=0;n<5;++n){
                        sxyz<T1> pt;
                        if(blti2D_raytrace(curr_pt, indEdges[n][0], indEdges[n][1], pt, threadNo, s) ==false) continue;
                        T1 dist0=pt.getDistance(nodes[indEdges[n][0]]);
                        T1 dist1=pt.getDistance(nodes[indEdges[n][1]]);
                        T1 t= (dist1*nodes[indEdges[n][0]].getTT(threadNo)+dist0*nodes[indEdges[n][1]].getTT(threadNo))/(dist0+dist1);
                        if (t>time)
                            continue;
                        t+=curr_pt.getDistance(pt)*s;
                        if (t<t_i ){
                            t_i=t;
                            Slow=s;
                            pt_i=pt;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            edgeNodestmp[0] = indEdges[n][0];
                            edgeNodestmp[1] = indEdges[n][1];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            cellNo=cellNoi;
                            
                        }
                        
                    }
                    for (T2 i=0;i<4;++i){
                        T1 t=nodes[neighbors[cellNoi][i]].getTT(threadNo);
                        if (t>time)
                            continue;
                        t+=s*nodes[neighbors[cellNoi][i]].getDistance(curr_pt);
                        if (t<t_i){
                            t_i=t;
                            Slow=s;
                            pt_i=sxyz<T1>(nodes[neighbors[cellNoi][i]]);
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            onNode = true;
                            onEdge = false;
                            onFace = false;
                            edgeNodesPrev = edgeNodes;
                            nodeNo=neighbors[cellNoi][i];
                            cellNo=cellNoi;
                        }
                    }
                    
                    // find next cell
                }
                if (!NearSource){
                    if ( t_i==std::numeric_limits<T1>::max()) {
                        std::cout << "\n\nWarning: finding raypath on edge failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                }
                if (onEdge){
                    edgeNodesPrev = edgeNodes;
                    edgeNodes[0] = edgeNodestmp[0];
                    edgeNodes[1] = edgeNodestmp[1];
                }
                prev_pt = curr_pt;
                curr_pt = pt_i;
                s2 = computeSlowness(curr_pt);
                
                tt += 0.5*(s1 + s2) * r_tmp.back().getDistance( curr_pt );
                s1 = s2;
                r_tmp.push_back( curr_pt );
                time=t_i-curr_pt.getDistance(prev_pt)*Slow;
            } else{ // on Face
                //////////////////////////
                // cout<<N<<endl;
                T2 cellNo1=findAdjacentCell2(faceNodes, cellNo,curr_pt);
                if (r_tmp.size()==1)
                    cellNo1=cellNo;
                sxyz<T1> pt_i;
                T1 Slow;
                ///////////////////////////
                std::array<T2,3> ind[8] = {
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][2] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][2], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][1], neighbors[cellNo][2], neighbors[cellNo][3] } },
                    {{ neighbors[cellNo1][0], neighbors[cellNo1][1], neighbors[cellNo1][2] }},
                    { { neighbors[cellNo1][0], neighbors[cellNo1][1], neighbors[cellNo1][3]}},
                    { { neighbors[cellNo1][0], neighbors[cellNo1][2], neighbors[cellNo1][3] } },
                    { { neighbors[cellNo1][1], neighbors[cellNo1][2], neighbors[cellNo1][3] } }
                };
                for ( size_t n=0; n<8; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                // there are 3 faces that might be intersected
                
                T1 t_i=std::numeric_limits<T1>::max();
                T2 face;
                T1 s1=0.25*(nodes[neighbors[cellNo][0]].getNodeSlowness()+
                            nodes[neighbors[cellNo][1]].getNodeSlowness()+
                            nodes[neighbors[cellNo][2]].getNodeSlowness()+
                            nodes[neighbors[cellNo][3]].getNodeSlowness());
                T1 s2=0.25*(nodes[neighbors[cellNo1][0]].getNodeSlowness()+
                            nodes[neighbors[cellNo1][1]].getNodeSlowness()+
                            nodes[neighbors[cellNo1][2]].getNodeSlowness()+
                            nodes[neighbors[cellNo1][3]].getNodeSlowness());
                
                for ( size_t n=0; n<8; ++n ) {
                    if ( ind[n] == faceNodes) continue;
                    if (NearSource){
                        // plotCell(cellNo, curr_pt, NodeSource-curr_pt);
                        //plotCell(cellNo1, curr_pt, sxyz<T1>(0.0,0.0,0.0));
                        std::array<T1,3> Barycenter;
                        if (blti_solver_around_source(NodeSource, curr_pt, ind[n], Barycenter)==true){
                            pt_i.x=Barycenter[0]*nodes[ind[n][0]].getX()+Barycenter[1]*nodes[ind[n][1]].getX()+Barycenter[2]*nodes[ind[n][2]].getX();
                            pt_i.y=Barycenter[0]*nodes[ind[n][0]].getY()+Barycenter[1]*nodes[ind[n][1]].getY()+Barycenter[2]*nodes[ind[n][2]].getY();
                            pt_i.z=Barycenter[0]*nodes[ind[n][0]].getZ()+Barycenter[1]*nodes[ind[n][1]].getZ()+Barycenter[2]*nodes[ind[n][2]].getZ();
                            
                            if(NodeSource.getDistance(pt_i)>NodeSource.getDistance(curr_pt))
                                continue;
                            for(T2 ni=0;ni<3;++ni){
                                if(abs(1.0-Barycenter[ni])<minDist*minDist){
                                    onNodePrev = onNode;
                                    onEdgePrev = onEdge;
                                    onFacePrev = onFace;
                                    onNode = true;
                                    onEdge = false;
                                    onFace = false;
                                    faceNodesPrev = faceNodes;
                                    nodeNo=ind[n][ni];
                                    break;
                                }
                                if(abs(Barycenter[ni])<minDist*minDist){
                                    onNodePrev = onNode;
                                    onEdgePrev = onEdge;
                                    onFacePrev = onFace;
                                    onNode = false;
                                    onEdge = true;
                                    onFace = false;
                                    faceNodesPrev = faceNodes;
                                    edgeNodes[0] = ind[n][(ni+1)%3];
                                    edgeNodes[1] = ind[n][(ni+2)%3];
                                    break;
                                }
                            }
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            onNode = false;
                            onEdge = false;
                            onFace = true;
                            faceNodesPrev = faceNodes;
                            face=static_cast<T1>(n);
                            break;
                        }
                        continue;
                    }
                    if ( ind[n] == faceNodesPrev) continue;
                    T1 s;
                    sxyz<T1> pt;
                    s=(n<4)?s1:s2;
                    if(blti_raytrace(curr_pt, ind[n], pt, threadNo, s)==false) continue;
                    T1 t=Interpolator<T1>::bilinearTime(pt, nodes[ind[n][0]], nodes[ind[n][1]], nodes[ind[n][2]], threadNo);
                    if (t>time)
                        continue;
                    t+=curr_pt.getDistance(pt)*s;
                    if (t<t_i){
                        t_i=t;
                        Slow=s;
                        face=static_cast<T1>(n);
                        pt_i=pt;
                        onNodePrev = onNode;
                        onEdgePrev = onEdge;
                        onFacePrev = onFace;
                        onNode = false;
                        onEdge = false;
                        onFace = true;
                    }
                }
                if (!NearSource){
                    std::array<T2,2> indEdges[9] = {
                        {{ neighbors[cellNo][0],neighbors[cellNo][1]}},
                        {{ neighbors[cellNo][0],neighbors[cellNo][2]}},
                        {{ neighbors[cellNo][0],neighbors[cellNo][3]}},
                        {{ neighbors[cellNo][1],neighbors[cellNo][2]}},
                        {{ neighbors[cellNo][1],neighbors[cellNo][3]}},
                        {{ neighbors[cellNo][2],neighbors[cellNo][3]}}
                    };
                    T2 forthNode=0;
                    for (T2 i=0; i<4;++i){
                        if (neighbors[cellNo1][i]!=neighbors[cellNo][0] &&
                            neighbors[cellNo1][i]!=neighbors[cellNo][1]&&
                            neighbors[cellNo1][i]!=neighbors[cellNo][2]&&
                            neighbors[cellNo1][i]!=neighbors[cellNo][3]){
                            forthNode=i;
                            break;
                        }
                    }
                    indEdges[6]={neighbors[cellNo1][forthNode],neighbors[cellNo1][(forthNode+1)%4]};
                    indEdges[7]={neighbors[cellNo1][forthNode],neighbors[cellNo1][(forthNode+2)%4]};
                    indEdges[8]={neighbors[cellNo1][forthNode],neighbors[cellNo1][(forthNode+3)%4]};
                    for(T2 n=0;n<9;++n){
                        T1 s;
                        sxyz<T1> pt;
                        s=(n<6)?s1:s2;
                        if(blti2D_raytrace(curr_pt, indEdges[n][0], indEdges[n][1], pt, threadNo, s) ==false) continue;
                        T1 dist0=pt.getDistance(nodes[indEdges[n][0]]);
                        T1 dist1=pt.getDistance(nodes[indEdges[n][1]]);
                        T1 t= (dist1*nodes[indEdges[n][0]].getTT(threadNo)+dist0*nodes[indEdges[n][1]].getTT(threadNo))/(dist0+dist1);
                        if (t>time)
                            continue;
                        t+=curr_pt.getDistance(pt)*s;
                        if (t<t_i ){
                            t_i=t;
                            Slow=s;
                            pt_i=pt;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            edgeNodes[0] = indEdges[n][0];
                            edgeNodes[1] = indEdges[n][1];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            faceNodesPrev = faceNodes;
                            
                        }
                        
                    }
                    for (T2 n=0;n<4;++n){
                        T1 t=nodes[neighbors[cellNo][n]].getTT(threadNo);
                        if(t<time){
                            sxyz<T1> pt;
                            t+=s1*nodes[neighbors[cellNo][n]].getDistance(curr_pt);
                            if(t<t_i){
                                t_i=t;
                                Slow=s1;
                                pt_i=sxyz<T1>(nodes[neighbors[cellNo][n]]);
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;
                                onNode = true;
                                onEdge = false;
                                onFace = false;
                                faceNodesPrev = faceNodes;
                                nodeNo=neighbors[cellNo][n];
                            }
                        }
                    }
                    if(nodes[neighbors[cellNo1][forthNode]].getTT(threadNo)<time){
                        sxyz<T1> pt;
                        if(nodes[neighbors[cellNo1][forthNode]].getTT(threadNo)+s2*nodes[neighbors[cellNo1][forthNode]].getDistance(curr_pt)<t_i){
                            t_i=nodes[neighbors[cellNo1][forthNode]].getTT(threadNo)+s2*nodes[neighbors[cellNo1][forthNode]].getDistance(curr_pt);
                            Slow=s2;
                            pt_i=sxyz<T1>(nodes[neighbors[cellNo1][forthNode]]);
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            onNode = true;
                            onEdge = false;
                            onFace = false;
                            faceNodesPrev = faceNodes;
                            nodeNo=neighbors[cellNo1][forthNode];
                        }
                    }
                    if (t_i==std::numeric_limits<T1>::max()){
                        std::cout << "\n\nWarning: finding raypath on face failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                        continue;
                    }
                }
                prev_pt = curr_pt;
                curr_pt = pt_i;
                s2 = computeSlowness(curr_pt);
                
                tt += 0.5*(s1 + s2) * r_tmp.back().getDistance( curr_pt );
                s1 = s2;
                r_tmp.push_back( curr_pt );
                time=t_i-curr_pt.getDistance(prev_pt)*Slow;
                if (onFace){
                    faceNodesPrev = faceNodes;
                    
                    faceNodes = ind[face];
                    if (face<4)
                        cellNo=findAdjacentCell2(faceNodes, cellNo,curr_pt);
                    else
                        cellNo=findAdjacentCell2(faceNodes, cellNo1,curr_pt);
                }
                
            }
            //            std::cout<<"distance= "<<1000.0*curr_pt.getDistance(prev_pt)<<endl;
            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( curr_pt.getDistance( Tx[nt] ) < minDist ) {
                        tt += t0[nt];
                        reachedTx = true;
                        break;
                    }
                    for(size_t n=0;n<4;++n){
                        sxyz<T1> NearTx={nodes[neighbors[txCell[nt]][n]]};
                        if (curr_pt.getDistance(NearTx)< minDist){
                            tt += t0[nt];
                            reachedTx = true;
                            r_tmp.push_back( Tx[nt] );
                            break;
                        }
                    }
                    
                }
            } else {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        for ( auto nc=nodes[txNode[nt]].getOwners().begin();
                             nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                            if ( cellNo == *nc ) {
                                std::array<T2,4> itmp = getPrimary(cellNo);
                                
                                if ( interpVel )
                                    s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                nodes[itmp[0]],
                                                                                nodes[itmp[1]],
                                                                                nodes[itmp[2]],
                                                                                nodes[itmp[3]]);
                                else
                                    s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                             nodes[itmp[0]],
                                                                             nodes[itmp[1]],
                                                                             nodes[itmp[2]],
                                                                             nodes[itmp[3]]);
                                
                                tt += t0[nt] + 0.5*(s1 + s2) * r_tmp.back().getDistance( Tx[nt] );
                                r_tmp.push_back( Tx[nt] );
                                reachedTx = true;
                                break;
                            }
                        }
                    } else {
                        if ( cellNo == txCell[nt] ) {
                            std::array<T2,4> itmp = getPrimary(cellNo);
                            if ( interpVel )
                                s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                            nodes[itmp[0]],
                                                                            nodes[itmp[1]],
                                                                            nodes[itmp[2]],
                                                                            nodes[itmp[3]]);
                            else
                                s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                         nodes[itmp[0]],
                                                                         nodes[itmp[1]],
                                                                         nodes[itmp[2]],
                                                                         nodes[itmp[3]]);
                            
                            tt += t0[nt] + 0.5*(s1 + s2) * r_tmp.back().getDistance( Tx[nt] );
                            r_tmp.push_back( Tx[nt] );
                            reachedTx = true;
                        } else {
                            for ( size_t nn=0; nn<txNeighborCells[nt].size(); ++nn ) {
                                if ( cellNo == txNeighborCells[nt][nn] ) {
                                    std::array<T2,4> itmp = getPrimary(cellNo);
                                    if ( interpVel )
                                        s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                    nodes[itmp[0]],
                                                                                    nodes[itmp[1]],
                                                                                    nodes[itmp[2]],
                                                                                    nodes[itmp[3]]);
                                    else
                                        s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                                 nodes[itmp[0]],
                                                                                 nodes[itmp[1]],
                                                                                 nodes[itmp[2]],
                                                                                 nodes[itmp[3]]);
                                    
                                    tt += t0[nt] + 0.5*(s1 + s2) * r_tmp.back().getDistance( Tx[nt] );
                                    r_tmp.push_back( Tx[nt] );
                                    reachedTx = true;
                                    break;
                                }
                            }
                        }
                    }
                    if ( reachedTx ) break;
                }
            }
        }
        // for inversion, the order should be from Tx to Rx, so we reorder...
        size_t npts = r_tmp.size();
        r_data.resize( npts );
        for ( size_t nn=0; nn<npts; ++nn ) {
            r_data[nn] = r_tmp[ npts-1-nn ];
        }
    }
    
    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::blti_raytrace(const sxyz<T1>& curr_pt,
                                             const std::array<T2,3>&face,
                                             sxyz<T1>& next_pt,
                                             const size_t threadNo,
                                             const T1& s)const{
        
        NODE *vertexA=&(nodes[face[0]]);
        NODE *vertexB=&(nodes[face[1]]);
        NODE *vertexC=&(nodes[face[2]]);
        
        if ( vertexA->getTT(threadNo) > vertexB->getTT(threadNo) )
            std::swap(vertexA, vertexB);
        if ( vertexA->getTT(threadNo) > vertexC->getTT(threadNo) )
            std::swap(vertexA, vertexC);
        
        T1 u = vertexB->getTT(threadNo) - vertexA->getTT(threadNo);
        T1 v = vertexC->getTT(threadNo) - vertexA->getTT(threadNo);
        sxyz<T1> v_b = { vertexC->getX() - vertexA->getX(),
            vertexC->getY() - vertexA->getY(),
            vertexC->getZ() - vertexA->getZ() };
        sxyz<T1> v_c = { vertexB->getX() - vertexA->getX(),
            vertexB->getY() - vertexA->getY(),
            vertexB->getZ() - vertexA->getZ() };
        
        sxyz<T1> v_n = cross(v_b, v_c);
        
        T1 b = norm( v_b );
        T1 c = norm( v_c );
        T1 d2 = dot(v_b, v_c);
        
        //       T1 alpha = acos(d2/(b*c) );
        
        T1 phi=norm(v_n);
        T1 d_tmp = -vertexA->getX()*v_n.x - vertexA->getY()*v_n.y - vertexA->getZ()*v_n.z;
        
        T1 k = -(d_tmp + v_n.x*curr_pt.x + v_n.y*curr_pt.y+ v_n.z*curr_pt.z)/norm2(v_n);
        
        sxyz<T1> pt;
        pt.x = curr_pt.x+ k*v_n.x;
        pt.y = curr_pt.y + k*v_n.y;
        pt.z = curr_pt.z + k*v_n.z;
        //        bool test=areCoplanar(pt, face[0], face[1], face[2]);
        //        sxyz<T1> vect=pt-curr_pt;
        //        T1 dd=dot(v_b,vect);
        T1 rho0 = curr_pt.getDistance( pt );
        
        // project point on AB
        sxyz<T1> v_pt = {pt.x-vertexA->getX(), pt.y-vertexA->getY(), pt.z-vertexA->getZ()};
        //// decomposition of Ap
        sxz<T1> AtA_Vect1={b*b,d2};
        sxz<T1> AtA_Vect2={d2,c*c};
        sxz<T1> Atb={dot(v_b,v_pt),dot(v_c,v_pt)};
        T1 DeT=det(AtA_Vect1,AtA_Vect2);
        T1 xi0=det(AtA_Vect1,Atb)/DeT;
        T1 zeta0=det(Atb,AtA_Vect2)/DeT;
        //        std::array<T1, 3> Weights;
        //      Interpolator<T1>::bilinearTriangleWeight(pt, *vertexA, *vertexB, *vertexC, Weights);
        
        T1 beta = u*b*b - v*d2;
        T1 gamma = v*c*c - u*d2;
        T1 w_tilde2 = (s*s*phi*phi-u*u*b*b-v*v*c*c+2.0*u*v*d2 );
        if (w_tilde2>0.0){
            T1 xi_tilde = -fabs(beta)*rho0/(phi*sqrt(w_tilde2));
            T1 zeta_tilde = -fabs(gamma)*rho0/(phi*sqrt(w_tilde2));
            T1 xi = xi_tilde + xi0;
            T1 zeta = zeta_tilde + zeta0;
            if ( 0.<xi && xi<1. && 0.<zeta && zeta<1. && 0.<(xi+zeta) && (xi+zeta)<1. ){
                next_pt.x=xi*vertexB->getX()+zeta*vertexC->getX()+(1-xi-zeta)*vertexA->getX();
                next_pt.y=xi*vertexB->getY()+zeta*vertexC->getY()+(1-xi-zeta)*vertexA->getY();
                next_pt.z=xi*vertexB->getZ()+zeta*vertexC->getZ()+(1-xi-zeta)*vertexA->getZ();
                return true;
            }else
                return false;
            
        }else
            return false;
    }
    
    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::blti2D_raytrace(const sxyz<T1>& curr_pt,
                                               const T2& node1,
                                               const T2& node2,
                                               sxyz<T1>& next_pt,
                                               const size_t threadNo,
                                               const T1& s) const{
        
        NODE *vertexA=&(nodes[node1]);
        NODE *vertexB=&(nodes[node2]);
        if ( vertexA->getTT(threadNo) > vertexB->getTT(threadNo) ) {
            std::swap(vertexA, vertexB);
        }
        T1 u = vertexB->getTT(threadNo) - vertexA->getTT(threadNo);
        
        sxyz<T1> v_b = { curr_pt.x- vertexA->getX(),
            curr_pt.y- vertexA->getY(),
            curr_pt.z- vertexA->getZ() };
        sxyz<T1> v_c = { vertexB->getX() - vertexA->getX(),
            vertexB->getY() - vertexA->getY(),
            vertexB->getZ() - vertexA->getZ() };
        
        T1 c = norm( v_c );
        
        T1 w2 = (s*s*c*c - u*u);
        if (w2<0.0)
            return false;
        T1 w = sqrt( w2 );
        T1 k = dot(v_b,v_c)/dot(v_c,v_c);
        sxyz<T1> pt;
        pt.x = vertexA->getX() + k*v_c.x;
        pt.y = vertexA->getY() + k*v_c.y;
        pt.z = vertexA->getZ() + k*v_c.z;

        T1 rho0 = curr_pt.getDistance( pt );

        T1 xi0 = k;
        T1 xi = xi0 -u*rho0/(w*c);
        
        if ( 0.0<xi && xi<1.0 ) {
            next_pt.x=(1-xi)*vertexA->getX()+xi*vertexB->getX();
            next_pt.y=(1-xi)*vertexA->getY()+xi*vertexB->getY();
            next_pt.z=(1-xi)*vertexA->getZ()+xi*vertexB->getZ();
            return true;
        }
        return false;
    }
    
    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::blti_solver_around_source(const sxyz<T1>& source,
                                                         const sxyz<T1>& curr_pt,
                                                         const std::array<T2,3>& face,
                                                         std::array<T1,3>& barycenters) const {
        T2 node1=face[0];
        T2 node2=face[1];
        T2 node3=face[2];
        sxyz<T1> PQ=source-curr_pt;
        sxyz<T1>PA={nodes[node1].getX()-curr_pt.x,nodes[node1].getY()-curr_pt.y,nodes[node1].getZ()-curr_pt.z};
        sxyz<T1>PB={nodes[node2].getX()-curr_pt.x,nodes[node2].getY()-curr_pt.y,nodes[node2].getZ()-curr_pt.z};
        sxyz<T1>PC={nodes[node3].getX()-curr_pt.x,nodes[node3].getY()-curr_pt.y,nodes[node3].getZ()-curr_pt.z};
        sxyz<T1> m=cross(PQ,PC);
        T1 u=dot(PB,m);
        T1 v=-dot(PA,m);
        if(signum(u)!=signum(v))
            return false;
        T1 w=tripleScalar(PQ, PB, PA);
        if(signum(u)!=signum(w))
            return false;
        T1 denom=1.0/(u+v+w);
        barycenters[0]=denom*u;
        barycenters[1]=denom*v;
        barycenters[2]=denom*w;
        return true;
    }
    

    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::check_pt_location(sxyz<T1> &curr_pt,
                                                 const std::array<T2,3> &ind,
                                                 bool &onNode,
                                                 T2 &nodeNo,
                                                 bool &onEdge,
                                                 std::array<T2,2> &edgeNodes,
                                                 bool &onFace,
                                                 std::array<T2,3> &faceNodes) const {
        // check if point is on vertex, edge or face
        
        bool break_flag = false;
        
        for ( size_t n=0; n<3; ++n ) {
            if ( nodes[ ind[n] ].getDistance( curr_pt ) < min_dist) {
                curr_pt = nodes[ ind[n] ];
                nodeNo = ind[n];
                onNode = true;
                onEdge = false;
                onFace = false;
                break_flag = true;
                return break_flag;
            }
        }
        
        for ( size_t n1=0; n1<3; ++n1 ) {
            size_t n2 = (n1+1)%3;
            if ( areCollinear(curr_pt, ind[n1], ind[n2]) ) {
                edgeNodes[0] = ind[n1];
                edgeNodes[1] = ind[n2];
                onNode = false;
                onEdge = true;
                onFace = false;
                break_flag = true;
                return break_flag;
            }
        }
        
        onNode = false;
        onEdge = false;
        onFace = true;
        
        faceNodes = ind;
        
        return break_flag;
    }
    

    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::check_pt_location(sxyz<T1> &curr_pt,
                                                 const std::vector<T2> &ind1,
                                                 const std::array<T2,3> &ind2,
                                                 bool &onNode,
                                                 T2 &nodeNo,
                                                 bool &onEdge,
                                                 std::array<T2,2> &edgeNodes,
                                                 bool &onFace,
                                                 std::array<T2,3> &faceNodes) const {
        // chech if point is on vertex, edge or face
        
        bool break_flag = false;
        
        for ( size_t n=0; n<4; ++n ) {
            if ( nodes[ ind1[n] ].getDistance( curr_pt ) < min_dist ) {
                curr_pt = nodes[ ind1[n] ];
                nodeNo = ind1[n];
                onNode = true;
                onEdge = false;
                onFace = false;
                break_flag = true;
                return break_flag;
            }
        }
        
        for ( size_t n1=0; n1<3; ++n1 ) {
            size_t n2 = (n1+1)%3;
            if ( areCollinear(curr_pt, ind2[n1], ind2[n2]) ) {
                edgeNodes[0] = ind2[n1];
                edgeNodes[1] = ind2[n2];
                onNode = false;
                onEdge = true;
                onFace = false;
                break_flag = true;
                return break_flag;
            }
        }
        
        onNode = false;
        onEdge = false;
        onFace = true;
        
        faceNodes = ind2;
        std::sort(faceNodes.begin(), faceNodes.end());
        
        return break_flag;
    }
    
    
    
    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::update_m_data(std::vector<sijv<T1>>& m_data,
                                             sijv<T1>& m,
                                             const std::set<T2>& allNodes,
                                             const sxyz<T1>& mid_pt,
                                             const T1 s,
                                             const T1 ds) const {
        std::vector<T1> w;
        T1 sum_w = 0.0;
        for ( auto it=allNodes.begin(); it!=allNodes.end(); ++it ) {
            w.push_back( 1./nodes[*it].getDistance( mid_pt ) );
            sum_w += w.back();
        }
        size_t nn=0;
        for ( auto it=allNodes.begin(); it!=allNodes.end(); ++it ) {
            m.j = *it;
            m.v = -s * ds * w[nn++]/sum_w;
            bool found = false;
            for ( size_t nm=0; nm<m_data.size(); ++nm ) {
                if ( m_data[nm].j == m.j ) {
                    m_data[nm].v += m.v;
                    found = true;
                    break;
                }
            }
            if ( found == false ) {
                m_data.push_back(m);
            }
        }
    }
    
    
    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::intersectVecTriangle(const T2 iO, const sxyz<T1> &vec,
                                                    const T2 iA, T2 iB, T2 iC,
                                                    sxyz<T1> &pt_i) const {
        
        sxyz<T1> OA = {nodes[iA].getX()-nodes[iO].getX(), nodes[iA].getY()-nodes[iO].getY(), nodes[iA].getZ()-nodes[iO].getZ()};
        // check if counterclockwise
        sxyz<T1> AB = {nodes[iB].getX()-nodes[iA].getX(),
            nodes[iB].getY()-nodes[iA].getY(),
            nodes[iB].getZ()-nodes[iA].getZ()};
        sxyz<T1> AC = {nodes[iC].getX()-nodes[iA].getX(),
            nodes[iC].getY()-nodes[iA].getY(),
            nodes[iC].getZ()-nodes[iA].getZ()};
        sxyz<T1> n = cross(AB, AC);
        if ( dot(OA, n) > 0. ) std::swap(iB, iC);
        
        sxyz<T1> OB = {nodes[iB].getX()-nodes[iO].getX(), nodes[iB].getY()-nodes[iO].getY(), nodes[iB].getZ()-nodes[iO].getZ()};
        sxyz<T1> OC = {nodes[iC].getX()-nodes[iO].getX(), nodes[iC].getY()-nodes[iO].getY(), nodes[iC].getZ()-nodes[iO].getZ()};
        
        T1 u, v, w;
        u = tripleScalar(vec, OC, OB);
        if ( u<0.0 ) return false;
        v = tripleScalar(vec, OA, OC);
        if ( v<0.0 ) return false;
        w = tripleScalar(vec, OB, OA);
        if ( w<0.0 ) return false;
        
        T1 den = 1./(u+v+w);
        u *= den;
        v *= den;
        w *= den;
        
        pt_i.x = u*nodes[iA].getX() + v*nodes[iB].getX() + w*nodes[iC].getX();
        pt_i.y = u*nodes[iA].getY() + v*nodes[iB].getY() + w*nodes[iC].getY();
        pt_i.z = u*nodes[iA].getZ() + v*nodes[iB].getZ() + w*nodes[iC].getZ();
        
        return true;
    }
    
    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::intersectVecTriangle(const sxyz<T1> &O, const sxyz<T1> &vec,
                                                    const T2 iA, T2 iB, T2 iC,
                                                    sxyz<T1> &pt_i) const {
        
        sxyz<T1> OA = {nodes[iA].getX()-O.x, nodes[iA].getY()-O.y, nodes[iA].getZ()-O.z};
        // check if counterclockwise
        sxyz<T1> AB = {nodes[iB].getX()-nodes[iA].getX(),
            nodes[iB].getY()-nodes[iA].getY(),
            nodes[iB].getZ()-nodes[iA].getZ()};
        sxyz<T1> AC = {nodes[iC].getX()-nodes[iA].getX(),
            nodes[iC].getY()-nodes[iA].getY(),
            nodes[iC].getZ()-nodes[iA].getZ()};
        sxyz<T1> n = cross(AB, AC);
        if ( dot(OA, n) > 0. ) std::swap(iB, iC);
        
        sxyz<T1> OB = {nodes[iB].getX()-O.x, nodes[iB].getY()-O.y, nodes[iB].getZ()-O.z};
        sxyz<T1> OC = {nodes[iC].getX()-O.x, nodes[iC].getY()-O.y, nodes[iC].getZ()-O.z};
        
        T1 u, v, w;
        u = tripleScalar(vec, OC, OB);
        if ( u<0.0 ) return false;
        v = tripleScalar(vec, OA, OC);
        if ( v<0.0 ) return false;
        w = tripleScalar(vec, OB, OA);
        if ( w<0.0 ) return false;
        
        T1 den = 1./(u+v+w);
        u *= den;
        v *= den;
        w *= den;
        
        pt_i.x = u*nodes[iA].getX() + v*nodes[iB].getX() + w*nodes[iC].getX();
        pt_i.y = u*nodes[iA].getY() + v*nodes[iB].getY() + w*nodes[iC].getY();
        pt_i.z = u*nodes[iA].getZ() + v*nodes[iB].getZ() + w*nodes[iC].getZ();
        
        return true;
    }
    
    
    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::areCollinear(const sxyz<T1> &pt, const T2 i0, const T2 i1) const {
        
        // http://mathworld.wolfram.com/Collinear.html
        //
        sxyz<T1> v1 = {pt.x-nodes[i0].getX(), pt.y-nodes[i0].getY(), pt.z-nodes[i0].getZ()};
        sxyz<T1> v2 = {pt.x-nodes[i1].getX(), pt.y-nodes[i1].getY(), pt.z-nodes[i1].getZ()};
        sxyz<T1> v3 = cross(v1, v2);
        return norm(v3)<small2;
        
    }
    
    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::areCoplanar(const sxyz<T1> &x1, const T2 i0, const T2 i1, const T2 i2) const {
        
        // http://mathworld.wolfram.com/Coplanar.html
        //
        sxyz<T1> x2 = {nodes[i0].getX(), nodes[i0].getY(), nodes[i0].getZ()};
        sxyz<T1> x3 = {nodes[i1].getX(), nodes[i1].getY(), nodes[i1].getZ()};
        sxyz<T1> x4 = {nodes[i2].getX(), nodes[i2].getY(), nodes[i2].getZ()};
        
        return fabs( dot( x3-x1, cross(x2-x1, x4-x3) ) )<small2;
    }
    
    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::areCoplanar(const sxyz<T1> &x1, const sxyz<T1> &x2,
                                           const sxyz<T1> &x3, const sxyz<T1> &x4) const {
        return (fabs( dot( x3-x1, cross(x2-x1, x4-x3) ) )<small2);
    }

    
    template<typename T1, typename T2, typename NODE>
    T2 Grid3Dun<T1,T2,NODE>::findAdjacentCell1(const std::array<T2,3> &faceNodes,
                                               const T2 nodeNo) const {
        
        std::vector<T2> cells;
        for ( auto nc0=nodes[faceNodes[0]].getOwners().begin(); nc0!=nodes[faceNodes[0]].getOwners().end(); ++nc0 ) {
            if ( std::find(nodes[faceNodes[1]].getOwners().begin(), nodes[faceNodes[1]].getOwners().end(), *nc0)!=nodes[faceNodes[1]].getOwners().end() &&
                std::find(nodes[faceNodes[2]].getOwners().begin(), nodes[faceNodes[2]].getOwners().end(), *nc0)!=nodes[faceNodes[2]].getOwners().end() )
                cells.push_back( *nc0 );
        }
        if ( cells.size() == 1 ) {
            return cells[0];
        }
        for ( auto nc0=nodes[nodeNo].getOwners().begin(); nc0!=nodes[nodeNo].getOwners().end(); ++nc0 ) {
            if ( *nc0 == cells[0] ) {
                return cells[1];
            } else if ( *nc0 == cells[1] ) {
                return cells[0];
            }
        }
        return std::numeric_limits<T2>::max();
    }
    
    template<typename T1, typename T2, typename NODE>
    T2 Grid3Dun<T1,T2,NODE>::findAdjacentCell2(const std::array<T2,3> &faceNodes,
                                               const T2 cellNo) const {
        
        std::vector<T2> cells;
        for ( auto nc0=nodes[faceNodes[0]].getOwners().begin(); nc0!=nodes[faceNodes[0]].getOwners().end(); ++nc0 ) {
            if ( std::find(nodes[faceNodes[1]].getOwners().begin(), nodes[faceNodes[1]].getOwners().end(), *nc0)!=nodes[faceNodes[1]].getOwners().end() &&
                std::find(nodes[faceNodes[2]].getOwners().begin(), nodes[faceNodes[2]].getOwners().end(), *nc0)!=nodes[faceNodes[2]].getOwners().end() )
                cells.push_back( *nc0 );
        }
        if ( cells.size() == 1 ) {
            return cells[0];
        }
        if ( cellNo == cells[0] ) {
            return cells[1];
        } else if ( cellNo == cells[1] ) {
            return cells[0];
        }
        return std::numeric_limits<T2>::max();
    }
    
    template<typename T1, typename T2, typename NODE>
    T2 Grid3Dun<T1,T2,NODE>::findAdjacentCell2(const std::array<T2,3> &faceNodes,
                                               const T2 & cellNo,
                                               const sxyz<T1>& curr_pt ) const {
        
        if (findAdjacentCell2(faceNodes, cellNo)!=cellNo &&
            findAdjacentCell2(faceNodes, cellNo)!=std::numeric_limits<T2>::max())
            return (findAdjacentCell2(faceNodes, cellNo));
        std::set<T2> AdjacentCells;
        for (size_t n1=0;n1<3;++n1){
            size_t n2((n1+1)%3),n3((n1+2)%3);
            for ( auto nc0=nodes[faceNodes[n1]].getOwners().begin(); nc0!=nodes[faceNodes[n1]].getOwners().end(); ++nc0 ) {
                for(size_t N=0;N<4;++N){
                    for ( auto nc=nodes[neighbors[(*nc0)][N]].getOwners().begin(); nc!=nodes[neighbors[(*nc0)][N]].getOwners().end(); ++nc ) {
                        for(size_t iD=0;iD<4;++iD){
                            size_t iA((iD+1)%4),iB((iD+2)%4),iC((iD+3)%4);
                            bool coarsetofine=(testInTriangle( &nodes[faceNodes[n1]],  &nodes[faceNodes[n2]],  &nodes[faceNodes[n3]], sxyz<T1>(nodes[neighbors[(*nc)][iA]])))&&(testInTriangle( &nodes[faceNodes[n1]],  &nodes[faceNodes[n2]],  &nodes[faceNodes[n3]], sxyz<T1>(nodes[neighbors[(*nc)][iB]])))&&
                            (testInTriangle( &nodes[faceNodes[n1]],  &nodes[faceNodes[n2]],  &nodes[faceNodes[n3]], sxyz<T1>(nodes[neighbors[(*nc)][iC]])))&&(*nc!=cellNo);
                            
                            bool fintocoarse=(testInTriangle( &nodes[neighbors[(*nc)][iA]],  &nodes[neighbors[(*nc)][iB]],  &nodes[neighbors[(*nc)][iC]], sxyz<T1>(nodes[faceNodes[n1]])))&&
                            (testInTriangle( &nodes[neighbors[(*nc)][iA]],  &nodes[neighbors[(*nc)][iB]],  &nodes[neighbors[(*nc)][iC]], sxyz<T1>(nodes[faceNodes[n2]])))&&
                            (testInTriangle( &nodes[neighbors[(*nc)][iA]],  &nodes[neighbors[(*nc)][iB]],  &nodes[neighbors[(*nc)][iC]], sxyz<T1>(nodes[faceNodes[n3]])))&&(*nc!=cellNo);
                            
                            if(fintocoarse || coarsetofine){
                                AdjacentCells.insert(*(nc));
                                break;
                            }
                        }
                    }
                }
            }
        }
        for(auto nc=AdjacentCells.begin();nc!=AdjacentCells.end();++nc){
            for(T2 iD=0;iD<4;++iD){
                T2 iA((iD+1)%4),iB((iD+2)%4),iC((iD+3)%4);
                if(testInTriangle(& nodes[neighbors[*nc][iA]], & nodes[neighbors[*nc][iB]]
                                  ,& nodes[neighbors[*nc][iC]], curr_pt)){
                    return (*nc);
                }
            }
        }
        return -1;
    }

    
    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::plotCell(const T2 cellNo, const sxyz<T1> &pt, const sxyz<T1> &g) const {
        
        
        if ( true ) {
            std::array<T2,4> itmp = getPrimary(cellNo);
            T2 i0 = itmp[0];
            T2 i1 = itmp[1];
            T2 i2 = itmp[2];
            T2 i3 = itmp[3];
            
            std::cout << "\nplot3(["<<nodes[ i0 ].getX()<<' ' << nodes[ i1 ].getX() <<"],["
            <<nodes[ i0 ].getY()<<' ' << nodes[ i1 ].getY() <<"],["
            <<nodes[ i0 ].getZ()<<' ' << nodes[ i1 ].getZ() <<"]); hold on;\n";
            std::cout << "plot3(["<<nodes[ i0 ].getX()<<' ' << nodes[ i2 ].getX() <<"],["
            <<nodes[ i0 ].getY()<<' ' << nodes[ i2 ].getY() <<"],["
            <<nodes[ i0 ].getZ()<<' ' << nodes[ i2 ].getZ() <<"])\n";
            std::cout << "plot3(["<<nodes[ i0 ].getX()<<' ' << nodes[ i3 ].getX() <<"],["
            <<nodes[ i0 ].getY()<<' ' << nodes[ i3 ].getY() <<"],["
            <<nodes[ i0 ].getZ()<<' ' << nodes[ i3 ].getZ() <<"])\n";
            std::cout << "plot3(["<<nodes[ i1 ].getX()<<' '<<nodes[ i2 ].getX()<<' '<<nodes[ i3 ].getX()<<' '<<nodes[ i1 ].getX()<<"],["
            <<nodes[ i1 ].getY()<<' '<<nodes[ i2 ].getY()<<' '<<nodes[ i3 ].getY()<<' '<<nodes[ i1 ].getY()<<"],["
            <<nodes[ i1 ].getZ()<<' '<<nodes[ i2 ].getZ()<<' '<<nodes[ i3 ].getZ()<<' '<<nodes[ i1 ].getZ()<<"])\n";
            std::cout << "plot3(["<<pt.x<< ' ' << pt.x+g.x<<"],["<<pt.y<< ' ' << pt.y+g.y<<"],["<<pt.z<< ' ' << pt.z+g.z<<"],'r')\naxis equal\n\n";
        }
    }
    
    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::getNeighborNodes(const T2 cellNo,
                                                std::set<NODE*> &nnodes) const {
        
        for ( size_t n=0; n<neighbors[cellNo].size(); ++n ) {
            if ( nodes[neighbors[cellNo][n]].isPrimary() ) {
                T2 nodeNo = neighbors[cellNo][n];
                nnodes.insert( &(nodes[nodeNo]) );
                if ( rp_method == 1 ) {  // second-order
                    for ( auto nc=nodes[nodeNo].getOwners().cbegin(); nc!=nodes[nodeNo].getOwners().cend(); ++nc ) {
                        for ( size_t nn=0; nn<neighbors[*nc].size(); ++nn ) {
                            if ( nodes[ neighbors[*nc][nn] ].isPrimary() ) {
                                nnodes.insert( &(nodes[ neighbors[*nc][nn] ]) );
                            }
                        }
                    }
                }
            }
        }
    }
    
    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::getNeighborNodesAB(const std::vector<NODE*>& ref_pt,
                                                  std::vector<std::vector<std::array<NODE*,3>>>& opp_pts) const {
        opp_pts.resize( ref_pt.size() );
        for ( size_t nr=0; nr<ref_pt.size(); ++nr ) {
            opp_pts[nr].resize( ref_pt[nr]->getOwners().size() );
            for ( size_t nc=0; nc<ref_pt[nr]->getOwners().size(); ++nc ) {
                T2 cellNo = ref_pt[nr]->getOwners()[nc];
                size_t ind=0;
                for ( size_t nn=0; nn<neighbors[cellNo].size(); ++nn ) {
                    if ( &(nodes[ neighbors[cellNo][nn] ]) != ref_pt[nr] && nodes[ neighbors[cellNo][nn] ].isPrimary() ) {
                        opp_pts[nr][nc][ind++] = &(nodes[ neighbors[cellNo][nn] ]);
                    }
                }
            }
        }
    }
    
    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::testInTriangle(const NODE *vertexA,
                                              const NODE *vertexB,
                                              const NODE *vertexC,
                                              const sxyz<T1> &E) const {
        if (vertexA->getDistance(E)+vertexB->getDistance(E)-vertexA->getDistance(*vertexB)<small2)
            return true;
        if (vertexA->getDistance(E)+vertexC->getDistance(E)-vertexA->getDistance(*vertexC)<small2)
            return true;
        if (vertexC->getDistance(E)+vertexB->getDistance(E)-vertexC->getDistance(*vertexB)<small2)
            return true;
        if (!areCoplanar(sxyz<T1>((*vertexA)), sxyz<T1>((*vertexB)), sxyz<T1>((*vertexC)), E))
            return false;
        T1 u, v, w;
        barycentric(vertexA, vertexB, vertexC, E, u, v, w);
        return v >= 0.0 && w >= 0.0 && (v + w) <= 1.0;
    }
    
    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::barycentric(const NODE *a,
                                           const NODE *b,
                                           const NODE *c,
                                           const sxyz<T1> &p,
                                           T1 &u, T1 &v, T1 &w) const {
        
        sxyz<T1> ab = {b->getX()-a->getX(), b->getY()-a->getY(), b->getZ()-a->getZ()};
        sxyz<T1> ac = {c->getX()-a->getX(), c->getY()-a->getY(), c->getZ()-a->getZ()};
        
        // Unnormalized triangle normal
        sxyz<T1> m = cross(ab, ac);
        
        // Nominators and one-over-denominator for u and v ratios
        T1 nu, nv, ood;
        
        // Absolute components for determining projection plane
        T1 x = fabs(m.x), y = fabs(m.y), z = fabs(m.z);
        
        // Compute areas in plane of largest projection
        if (x >= y && x >= z) {
            // x is largest, project to the yz plane
            nu = triangleArea2D(p.y, p.z, b->getY(), b->getZ(), c->getY(), c->getZ()); // Area of PBC in yz plane
            nv = triangleArea2D(p.y, p.z, c->getY(), c->getZ(), a->getY(), a->getZ()); // Area of PCA in yz plane
            ood = 1.0 / m.x; // 1/(2*area of ABC in yz plane)
        } else if (y >= x && y >= z) {
            // y is largest, project to the xz plane
            nu = triangleArea2D(p.x, p.z, b->getX(), b->getZ(), c->getX(), c->getZ());
            nv = triangleArea2D(p.x, p.z, c->getX(), c->getZ(), a->getX(), a->getZ());
            ood = 1.0 / -m.y;
        } else {
            // z is largest, project to the xy plane
            nu = triangleArea2D(p.x, p.y, b->getX(), b->getY(), c->getX(), c->getY());
            nv = triangleArea2D(p.x, p.y, c->getX(), c->getY(), a->getX(), a->getY());
            ood = 1.0 / m.z;
        }
        u = nu * ood;
        v = nv * ood;
        w = 1.0 - u - v;
    }
    
    
    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dun<T1,T2,NODE>::computeSlowness( const sxyz<T1>& Rx ) const {
        
        //Calculate the slowness of any point that is not on a node
        
        T2 cellNo = this->getCellNo( Rx );

        std::vector<NODE*> interpNodes;
        
        for (size_t n=0; n < neighbors[ cellNo ].size(); n++){
            if ( nodes[neighbors[ cellNo ][n] ].isPrimary() ){
                interpNodes.push_back( &(nodes[neighbors[ cellNo ][n] ]) );
            }
        }
        if ( interpVel )
            return Interpolator<T1>::trilinearTriangleVel( Rx, interpNodes );
        else
            return Interpolator<T1>::trilinearTriangle( Rx, interpNodes );
    }
    
    
    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dun<T1,T2,NODE>::computeSlowness(const sxyz<T1>& curr_pt,
                                             const bool onNode,
                                             const T2 nodeNo,
                                             const bool onEdge,
                                             const std::array<T2,2>& edgeNodes,
                                             const std::array<T2,3>& faceNodes) const {
        if ( onNode ) {
            return nodes[nodeNo].getNodeSlowness();
        } else {
            if ( interpVel ) {
                if ( onEdge ) {
                    return Interpolator<T1>::linearVel(curr_pt,
                                                       nodes[edgeNodes[0]],
                                                       nodes[edgeNodes[1]]);
                } else {
                    return Interpolator<T1>::bilinearTriangleVel(curr_pt,
                                                                 nodes[faceNodes[0]],
                                                                 nodes[faceNodes[1]],
                                                                 nodes[faceNodes[2]]);
                }
            } else {
                if ( onEdge ) {
                    return Interpolator<T1>::linear(curr_pt,
                                                    nodes[edgeNodes[0]],
                                                    nodes[edgeNodes[1]]);
                } else {
                    return Interpolator<T1>::bilinearTriangle(curr_pt,
                                                              nodes[faceNodes[0]],
                                                              nodes[faceNodes[1]],
                                                              nodes[faceNodes[2]]);
                }
            }
        }
    }
    
}

#endif

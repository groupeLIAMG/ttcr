 
//
//  Grid3Dui.h
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

#ifndef ttcr_Grid3Dui_h
#define ttcr_Grid3Dui_h

#include <cassert>

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>
#include <math.h>
#include <Eigen/Dense>

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
    class Grid3Dui : public Grid3D<T1,T2> {
    public:
        Grid3Dui(const std::vector<sxyz<T1>>& no,
                 const std::vector<tetrahedronElem<T2>>& tet,
                 const size_t nt=1) :
        nThreads(nt),
        nPrimary(static_cast<T2>(no.size())),
        source_radius(0.0),
        nodes(std::vector<NODE>(no.size(), NODE(nt))),
        neighbors(std::vector<std::vector<T2>>(tet.size())),
        tetrahedra(tet)
        {}

        virtual ~Grid3Dui() {}

        void setSlowness(const T1 s) {
            for ( size_t n=0; n<nodes.size(); ++n ) {
                nodes[n].setNodeSlowness( s[n] );
            }
        }

        int setSlowness(const T1 *s, const size_t ns) {
            if ( nodes.size() != ns ) {
                std::cerr << "Error: slowness vectors of incompatible size.";
                return 1;
            }
            for ( size_t n=0; n<nodes.size(); ++n ) {
                nodes[n].setNodeSlowness( s[n] );
            }
            return 0;
        }

        int setSlowness(const std::vector<T1>& s) {
            if ( nodes.size() != s.size() ) {
                std::cerr << "Error: slowness vectors of incompatible size.";
                return 1;
            }
            for ( size_t n=0; n<nodes.size(); ++n ) {
                nodes[n].setNodeSlowness( s[n] );
            }
            return 0;
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
        T2 getCellNo(const sxyz<T1>& pt) const {
            for ( T2 n=0; n<tetrahedra.size(); ++n ) {
                if ( insideTetrahedron(pt, n) ) {
                    return n;
                }
            }
            // if no cell is found
            T2 ClosestNode=0;
            T1 distance=std::numeric_limits<T1>::max();
            for(auto Node=nodes.begin();Node!=nodes.begin()+nPrimary;++Node){
                T1 dist=pt.getDistance(*Node);
                if(dist<distance){
                    distance=dist;
                    ClosestNode=Node->getGridIndex();
                }
            }
            T1 MinVolumeDiff=std::numeric_limits<T1>::max();
            T2 Cell;
                 
            for(auto tet=nodes[ClosestNode].getOwners().begin();tet!=nodes[ClosestNode].getOwners().end();++tet){
                sxyz<T1> v1 = { nodes[ neighbors[*tet][0] ]};
                sxyz<T1> v2 = { nodes[ neighbors[*tet][1] ]};
                sxyz<T1> v3 ={ nodes[ neighbors[*tet][2] ]};
                sxyz<T1> v4 = { nodes[ neighbors[*tet][3] ]};
                
                T1 D0 = 1.e6*det4(v1, v2, v3, v4);
                T1 D1 = 1.e6*det4( pt, v2, v3, v4);
                T1 D2 = 1.e6*det4(v1,  pt, v3, v4);
                T1 D3 = 1.e6*det4(v1, v2,  pt, v4);
                T1 D4 = 1.e6*det4(v1, v2, v3,  pt);
  
                T1 VolumeDiff=abs(abs(D0)-abs(D1)-abs(D2)-abs(D3)-abs(D4));
                 plotCell(*tet, pt, sxyz<T1>(0.0,0.0,0.00));
                if(VolumeDiff<MinVolumeDiff){
                    MinVolumeDiff=VolumeDiff;
                    Cell=*tet;
                }
            }
            plotCell(Cell, pt, sxyz<T1>(0.0,0.0,0.00));
            return Cell;
        }
         int checkPts(const std::vector<sxyz<T1>>&) const;

        virtual int raytrace(const std::vector<sxyz<T1>>& Tx,
                             const std::vector<T1>& t0,
                             const std::vector<sxyz<T1>>& Rx,
                             std::vector<T1>& traveltimes,
                             const size_t threadNo=0) const { return 0; }

        virtual int raytrace(const std::vector<sxyz<T1>>& Tx,
                             const std::vector<T1>& t0,
                             const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                             std::vector<std::vector<T1>*>& traveltimes,
                             const size_t=0) const { return 0; }

        virtual int raytrace(const std::vector<sxyz<T1>>& Tx,
                             const std::vector<T1>& t0,
                             const std::vector<sxyz<T1>>& Rx,
                             std::vector<T1>& traveltimes,
                             std::vector<std::vector<sxyz<T1>>>& r_data,
                             const size_t threadNo=0) const { return 0; }

        virtual int raytrace(const std::vector<sxyz<T1>>& Tx,
                             const std::vector<T1>& t0,
                             const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                             std::vector<std::vector<T1>*>& traveltimes,
                             std::vector<std::vector<std::vector<sxyz<T1>>>*>& r_data,
                             const size_t=0) const { return 0; }
        T2 DivideTetraHedron (const std::vector<sxyz<T1>>& Tetrahedron,
                              std::vector<std::vector<sxyz<T1>>>& NewTetrahedron)const;
        T2 addNodes(const sxyz<T1>& Pnt,std::vector<tetrahedronElem<T2>>& Delete_CeLL,const size_t & nt) const;
        void saveTT(const std::string &, const int, const size_t nt=0,
                    const bool vtkFormat=0) const;
        void delateCells(const std::vector<tetrahedronElem<T2>> & delatedCells) const;
        sxyz<T1> reflectedGradient(const sxyz<T1>& g1, const std::array<T2, 3>& faces,const T2 & cellNo) const;
#ifdef VTK
        void saveModelVTU(const std::string &, const bool saveSlowness=true,
                          const bool savePhysicalEntity=false) const;
#endif

        const size_t getNthreads() const { return nThreads; }

    protected:
        const size_t nThreads;
        T2 nPrimary;
        T1 source_radius;
        mutable std::vector<NODE> nodes;
        mutable std::vector<std::vector<T2>> neighbors;  // nodes common to a cell
        mutable std::vector<tetrahedronElem<T2>> tetrahedra;

        T1 computeDt(const NODE& source, const NODE& node) const;

        T1 computeDt(const NODE& source, const sxyz<T1>& node, T1 slo) const;
        
        bool BLTISolver_ArroundSource(const sxyz<T1>& Source,const sxyz<T1>& curr_pt,const std::array<T2,3>& face, std::array<T1,3>& barycenters) const;

        T1 getTraveltime(const sxyz<T1>& Rx,
                         const std::vector<NODE>& nodes,
                         const size_t threadNo) const;
        bool insideTetrahedron(const sxyz<T1>&, const T2) const;


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
        bool blti_raytrace(const sxyz<T1> & curr_pt,const std::array<T2, 3> &faces, sxyz<T1> & next_pt,const size_t threadNo, const T1 & s) const;
        bool blti2D_raytrace(const sxyz<T1> & curr_pt,const T2 & node1,const T2 & node2, sxyz<T1> & next_pt,const size_t threadNo, const T1 & s) const;
        T1 local2Dsolver(const NODE *vertexA,
                         const NODE *vertexB,
                         const NODE *vertexC,
                         const size_t threadNo) const;
        void ObtusAngle2Dsolver(const NODE *vertexA,
                              const NODE *vertexB,
                              const NODE *vertexC,
                              std::vector<std::vector<T2>> & Vertexes) const;
        int solveEq23(const T1 a[], const T1 b[], T1 n[][3]) const;

        void barycentric(const NODE *a, const NODE *b, const NODE *c,
                         const sxyz<T1> &p, T1 &u, T1 &v, T1 &w) const;

        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const sxyz<T1> &Rx,
                        std::vector<sxyz<T1>> &r_data,
                        const size_t threadNo) const;
        void getRaypathBLIT(const std::vector<sxyz<T1>>& Tx,
                        const sxyz<T1> &Rx,
                        std::vector<sxyz<T1>> &r_data,
                        const size_t threadNo) const;

        void getRaypath_ho(const std::vector<sxyz<T1>>& Tx,
                           const sxyz<T1> &Rx,
                           std::vector<sxyz<T1>> &r_data,
                           const size_t threadNo) const;

        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const sxyz<T1> &Rx,
                        std::vector<sxyz<T1>> &r_data,
                        std::vector<sijv<T1>>& m_data,
                        const size_t RxNo,
                        const size_t threadNo) const;

        void getRaypath_ho(const std::vector<sxyz<T1>>& Tx,
                           const sxyz<T1> &Rx,
                           std::vector<sxyz<T1>> &r_data,
                           std::vector<sijv<T1>>& m_data,
                           const size_t RxNo,
                           const size_t threadNo) const;
        void getRaypathABM(const std::vector<sxyz<T1>>& Tx,
                         const sxyz<T1> &Rx,std::vector<sxyz<T1>> &r_data,
                           const size_t threadNo)const;
        bool intersectVecTriangle(const T2 iO, const sxyz<T1> &vec,
                                  const T2 iA, T2 iB, T2 iC,
                                  sxyz<T1> &pt_i) const;
        bool intersectVecTriangle(const sxyz<T1> &O, const sxyz<T1> &vec,
                                  const T2 iA, T2 iB, T2 iC,
                                  sxyz<T1> &pt_i) const;

        bool areCollinear(const sxyz<T1> &pt, const T2 i0, const T2 i1) const;
        bool areCoplanar(const sxyz<T1> &pt, const T2 i0, const T2 i1, const T2 i2) const;
        bool areCoplanar(const sxyz<T1> &pt1, const sxyz<T1> &pt2,
                         const sxyz<T1> &pt3,const sxyz<T1> &pt4) const;
        T2 findAdjacentCell1(const std::array<T2,3> &faceNodes, const T2 nodeNo) const;
        T2 findAdjacentCell2(const std::array<T2,3> &faceNodes, const T2 cellNo) const;
        T2 findAdjacentCell2(const std::array<T2,3> &faceNodes, const T2 &cellNo,const sxyz<T1>& curr_pt) const;
        void findAdjacentCell2(const std::array<T2,3> &faceNodes, const T2& cellNo,std::set<T2>& AdjacentCells) const;
       void getNeighborNodes(const T2 cellNo, std::set<NODE*> &nnodes) const;
        void getNeighborNodes2(const T2 cellNo, std::set<NODE*> &nnodes) const;
        void plotCell(const T2 cellNo, const sxyz<T1> &pt, const sxyz<T1> &g) const;
        bool testInTriangle(const NODE *vertexA,
                            const NODE *vertexB,
                            const NODE *vertexC,
                            const sxyz<T1> &E) const;
        T1 AreaTriangle(const NODE *vertexA,
                            const NODE *vertexB,
                            const NODE *vertexC) const;
        T1 AreaTriangle(const NODE *vertexA,
                        const NODE *vertexB,
                        const sxyz<T1>&E) const;
        bool areInLimits(const T2 & Node) const;

        T1 computeSlowness( const sxyz<T1>& Rx ) const;
        T1 Highorder_computeSlowness( const sxyz<T1>& Rx ) const;
        T1 computeSlowness( const sxyz<T1>& Rx,const T2 & cell ) const;
        bool AreSameFace(const std::array<T2, 3>& Face1,const std::array<T2, 3>& Face2)const;
    };
  template<typename T1,typename T2,typename NODE>
    T2  Grid3Dui<T1,T2,NODE>::DivideTetraHedron (const std::vector<sxyz<T1>>& Tetrahedron,
                                                   std::vector<std::vector<sxyz<T1>>>& NewTetrahedron) const {

        
        
//            for(size_t iD=0;iD<4;++iD){
//                T2 iB = (iD+1)%4, iC=(iD+2)%4, iA=(iD+3)%4;
//                T1 ang1=dot(Tetrahedron[iB]-Tetrahedron[iA],Tetrahedron[iC]-Tetrahedron[iA])/(norm(Tetrahedron[iB]-Tetrahedron[iA])*norm(Tetrahedron[iC]-Tetrahedron[iA]));
//                ang1=acos(ang1)*180.0/(4*atan(1));
//                T1 ang2=dot(Tetrahedron[iA]-Tetrahedron[iB],Tetrahedron[iC]-Tetrahedron[iB])/(norm(Tetrahedron[iA]-Tetrahedron[iB])*norm(Tetrahedron[iC]-Tetrahedron[iB]));
//                ang2=acos(ang2)*180.0/(4*atan(1));
//                T1 ang3=dot(Tetrahedron[iB]-Tetrahedron[iC],Tetrahedron[iA]-Tetrahedron[iC])/(norm(Tetrahedron[iB]-Tetrahedron[iC])*norm(Tetrahedron[iA]-Tetrahedron[iC]));
//                ang3=acos(ang3)*180.0/(4*atan(1));
//                if(ang1>92.0||ang2>92.0||ang3>92.0){
//                    std::cout<<"we have an obtuse angle in the origin tetrahedron: angle1  : "<<ang1<<" ; angle 2 : ";
//                    std::cout<<ang2<<" ; angle 3 : "<< ang3<<std::endl;
//                }
//            }

        
        std::vector<sxyz<T1>> Tetra;
        std::vector<sxyz<T1>> MidEdges;
        sxyz<T1> midEdge;
        for(size_t iA=0;iA<4;++iA){
            Tetra.push_back(Tetrahedron[iA]);
            T2 iB = (iA+1)%4, iC=(iA+2)%4, iD=(iA+3)%4;
            midEdge=(Tetrahedron[iA]+Tetrahedron[iB]);midEdge/=2.0;Tetra.push_back(midEdge);
            if(std::find(MidEdges.begin(),MidEdges.end(),midEdge)==MidEdges.end()) MidEdges.push_back(midEdge);
            midEdge=(Tetrahedron[iA]+Tetrahedron[iC]);midEdge/=2.0;Tetra.push_back(midEdge);
            if(std::find(MidEdges.begin(),MidEdges.end(),midEdge)==MidEdges.end()) MidEdges.push_back(midEdge);
            midEdge=(Tetrahedron[iA]+Tetrahedron[iD]);midEdge/=2.0;Tetra.push_back(midEdge);
            if(std::find(MidEdges.begin(),MidEdges.end(),midEdge)==MidEdges.end()) MidEdges.push_back(midEdge);
            NewTetrahedron.push_back(Tetra);Tetra.clear();
        }
        auto it = MidEdges.cbegin();
        if(areCoplanar(*it, *(std::next(it,2)), *(std::next(it,5)), *(std::next(it,3)))){
            sxyz<T1> vect1(*it-*(std::next(it,2)));
            sxyz<T1> vect2(*it-*(std::next(it,3)));
            if(dot(vect1,vect2)<=0){
                // fifth tetrahedron
                Tetra.push_back(*it);Tetra.push_back(*(std::next(it,2)));
                Tetra.push_back(*(std::next(it,4)));Tetra.push_back(*(std::next(it,5)));
                NewTetrahedron.push_back(Tetra);Tetra.clear();
                // sixth tetrahedron
                Tetra.push_back(*it);Tetra.push_back(*(std::next(it,4)));
                Tetra.push_back(*(std::next(it,3)));Tetra.push_back(*(std::next(it,5)));
                NewTetrahedron.push_back(Tetra);Tetra.clear();
                // seventh tetrahedron
                Tetra.push_back(*it);Tetra.push_back(*(std::next(it,5)));
                Tetra.push_back(*(std::next(it,1)));Tetra.push_back(*(std::next(it,2)));
                NewTetrahedron.push_back(Tetra);Tetra.clear();
                // eigthth tetrahedron
                Tetra.push_back(*it);Tetra.push_back(*(std::next(it,3)));
                Tetra.push_back(*(std::next(it,5)));Tetra.push_back(*(std::next(it,1)));
                NewTetrahedron.push_back(Tetra);Tetra.clear();
            }
            vect1=(*(std::next(it,2))-*(std::next(it,5)));
            vect2=(*(std::next(it,2))-*it);
            if(dot(vect1,vect2)<=0){
                // fifth tetrahedron
                Tetra.push_back(*it);Tetra.push_back(*(std::next(it,2)));
                Tetra.push_back(*(std::next(it,3)));Tetra.push_back(*(std::next(it,4)));
                NewTetrahedron.push_back(Tetra);Tetra.clear();
                // sixth tetrahedron
                Tetra.push_back(*(std::next(it,2)));Tetra.push_back(*(std::next(it,3)));
                Tetra.push_back(*(std::next(it,5)));Tetra.push_back(*(std::next(it,4)));
                NewTetrahedron.push_back(Tetra);Tetra.clear();
                // seventh tetrahedron
                Tetra.push_back(*(std::next(it,1)));Tetra.push_back(*(std::next(it,2)));
                Tetra.push_back(*(std::next(it,3)));Tetra.push_back(*(std::next(it,5)));
                NewTetrahedron.push_back(Tetra);Tetra.clear();
                // eigthth tetrahedron
                Tetra.push_back(*(std::next(it,1)));Tetra.push_back(*(std::next(it,2)));
                Tetra.push_back(*(std::next(it,3)));Tetra.push_back(*it);
                NewTetrahedron.push_back(Tetra);Tetra.clear();

            }
        }
//        T2 k=0;
//        for(auto nt =NewTetrahedron.begin();nt!=NewTetrahedron.end();++nt){
//            for(size_t iD=0;iD<4;++iD){
//                T2 iB = (iD+1)%4, iC=(iD+2)%4, iA=(iD+3)%4;
//                T1 ang1=dot((*nt)[iB]-(*nt)[iA],(*nt)[iC]-(*nt)[iA])/(norm((*nt)[iB]-(*nt)[iA])*norm((*nt)[iC]-(*nt)[iA]));
//                ang1=acos(ang1)*180.0/(4*atan(1));
//                T1 ang2=dot((*nt)[iA]-(*nt)[iB],(*nt)[iC]-(*nt)[iB])/(norm((*nt)[iA]-(*nt)[iB])*norm((*nt)[iC]-(*nt)[iB]));
//                ang2=acos(ang2)*180.0/(4*atan(1));
//                T1 ang3=dot((*nt)[iB]-(*nt)[iC],(*nt)[iA]-(*nt)[iC])/(norm((*nt)[iB]-(*nt)[iC])*norm((*nt)[iA]-(*nt)[iC]));
//                ang3=acos(ang3)*180.0/(4*atan(1));
//                if(ang1>92.0||ang2>92.0||ang3>92.0){
//                    //std::cout<<"we have an obtuse angle in the new tetrahedron N "<<k<<": angle1  : "<<ang1<<" ; angle 2 : ";
//                    //std::cout<<ang2<<" ; angle 3 : "<< ang3<<std::endl;
//                }
//            }
//            k++;
//        }
//        if(areCoplanar(*(std::next(it,1)), *(std::next(it,2)), *(std::next(it,4)), *(std::next(it,3)))){
//            sxyz<T1> vect1(*(std::next(it,1))-*(std::next(it,2)));
//            sxyz<T1> vect2(*(std::next(it,1))-*(std::next(it,3)));
//            if(dot(vect1,vect2)/(norm(vect1)*norm(vect2))<=0){
//
//                // fifth tetrahedron
//                Tetra.push_back(*it);Tetra.push_back(*(std::next(it,2)));
//                Tetra.push_back(*(std::next(it,1)));Tetra.push_back(*(std::next(it,4)));
//                NewTetrahedron.push_back(Tetra);Tetra.clear();
//                // sixth tetrahedron
//                Tetra.push_back(*(std::next(it,2)));Tetra.push_back(*(std::next(it,5)));
//                Tetra.push_back(*(std::next(it,1)));Tetra.push_back(*(std::next(it,4)));
//                NewTetrahedron.push_back(Tetra);Tetra.clear();
//                // seventh tetrahedron
//                Tetra.push_back*((std::next(it,1)));Tetra.push_back(*(std::next(it,4)));
//                Tetra.push_back(*(std::next(it,3)));Tetra.push_back(*(std::next(it,5)));
//                NewTetrahedron.push_back(Tetra);Tetra.clear();
//                // eigthth tetrahedron
//                Tetra.push_back(*(std::next(it,1)));Tetra.push_back(*(std::next(it,4)));
//                Tetra.push_back(*(std::next(it,3)));Tetra.push_back(*it);
//                NewTetrahedron.push_back(Tetra);Tetra.clear();
//                return 0;
//            }else{              // fifth tetrahedron
//                Tetra.push_back(*it);Tetra.push_back(*(std::next(it,2)));
//                Tetra.push_back(*(std::next(it,4)));Tetra.push_back(*(std::next(it,3)));
//                NewTetrahedron.push_back(Tetra);Tetra.clear();
//                // sixth tetrahedron
//                Tetra.push_back(*it);Tetra.push_back(*(std::next(it,1)));
//                Tetra.push_back(*(std::next(it,3)));Tetra.push_back(*(std::next(it,2)));
//                NewTetrahedron.push_back(Tetra);Tetra.clear();
//                // seventh tetrahedron
//                Tetra.push_back(*(std::next(it,2)));Tetra.push_back(*(std::next(it,5)));
//                Tetra.push_back(*(std::next(it,4)));Tetra.push_back(*(std::next(it,3)));
//                NewTetrahedron.push_back(Tetra);Tetra.clear();
//                // eigthth tetrahedron
//                Tetra.push_back(*(std::next(it,2)));Tetra.push_back(*(std::next(it,3)));
//                Tetra.push_back(*(std::next(it,5)));Tetra.push_back(*(std::next(it,1)));
//                NewTetrahedron.push_back(Tetra);Tetra.clear();
//                return 0;
//
//            }
//        }
//        if(areCoplanar(*(it), *(std::next(it,1)), *(std::next(it,4)), *(std::next(it,5)))){
//            sxyz<T1> vect1((*it)-*(std::next(it,1)));
//            sxyz<T1> vect2((*it)-*(std::next(it,4)));
//            if(dot(vect1,vect2)/(norm(vect1)*norm(vect2))<=0){
//
//                // fifth tetrahedron
//                Tetra.push_back(*it);Tetra.push_back(*(std::next(it,2)));
//                Tetra.push_back(*(std::next(it,5)));Tetra.push_back(*(std::next(it,4)));
//                NewTetrahedron.push_back(Tetra);Tetra.clear();
//                // sixth tetrahedron
//                Tetra.push_back(*it);Tetra.push_back(*(std::next(it,5)));
//                Tetra.push_back(*(std::next(it,3)));Tetra.push_back(*(std::next(it,4)));
//                NewTetrahedron.push_back(Tetra);Tetra.clear();
//                // seventh tetrahedron
//                Tetra.push_back(*it);Tetra.push_back(*(std::next(it,2)));
//                Tetra.push_back(*(std::next(it,1)));Tetra.push_back(*(std::next(it,5)));
//                NewTetrahedron.push_back(Tetra);Tetra.clear();
//                // eigthth tetrahedron
//                Tetra.push_back(*(std::next(it,5)));Tetra.push_back(*(std::next(it,1)));
//                Tetra.push_back(*(std::next(it,2)));Tetra.push_back(*it);
//                NewTetrahedron.push_back(Tetra);Tetra.clear();
//                return 0;
//            }else{              // fifth tetrahedron
//                Tetra.push_back(*it);Tetra.push_back(*(std::next(it,2)));
//                Tetra.push_back(*(std::next(it,4)));Tetra.push_back(*(std::next(it,1)));
//                NewTetrahedron.push_back(Tetra);Tetra.clear();
//                // sixth tetrahedron
//                Tetra.push_back(*it);Tetra.push_back(*(std::next(it,1)));
//                Tetra.push_back(*(std::next(it,3)));Tetra.push_back(*(std::next(it,4)));
//                NewTetrahedron.push_back(Tetra);Tetra.clear();
//                // seventh tetrahedron
//                Tetra.push_back(*(std::next(it,1)));Tetra.push_back(*(std::next(it,5)));
//                Tetra.push_back(*(std::next(it,4)));Tetra.push_back(*(std::next(it,3)));
//                NewTetrahedron.push_back(Tetra);Tetra.clear();
//                // eigthth tetrahedron
//                Tetra.push_back(*(std::next(it,2)));Tetra.push_back(*(std::next(it,4)));
//                Tetra.push_back(*(std::next(it,5)));Tetra.push_back(*(std::next(it,1)));
//                NewTetrahedron.push_back(Tetra);Tetra.clear();
//                return 0;
//
//            }
//        }
        return 0;
    }
    template<typename T1,typename T2,typename NODE>
    T2 Grid3Dui<T1,T2,NODE>::addNodes(const sxyz<T1>& Pnt,std::vector<tetrahedronElem<T2>>& Delete_CeLL,const size_t & nt)const {
        T2 CellNO=getCellNo(Pnt);
        std::set<T2> Celli_indices;
        for(auto Nodeneighbor=neighbors[CellNO].begin();Nodeneighbor!=neighbors[CellNO].end();++Nodeneighbor){
            for(T2 nc0=0;nc0!=nodes[*Nodeneighbor].getOwners().size();++nc0){
                T2 NearCell=nodes[*Nodeneighbor].getOwners()[nc0];
                for (T2 nn=0;nn<4;++nn){
                    for (T2 nc=0;nc<nodes[neighbors[NearCell][nn]].getOwners().size();++nc){
                        T2 Celli=nodes[neighbors[NearCell][nn]].getOwners()[nc];
                        if(Celli_indices.find(Celli)!=Celli_indices.end()) continue;
                        Celli_indices.insert(Celli);
                        if (nodes[neighbors[Celli][0]].getprimary()!=10 && nodes[neighbors[Celli][1]].getprimary()!=10 &&
                            nodes[neighbors[Celli][2]].getprimary()!=10 && nodes[neighbors[Celli][3]].getprimary()!=10)
                            Delete_CeLL.push_back(tetrahedra[Celli]);
                        std::vector<sxyz<T1>> Tetrahedron;
                        for(size_t n=0;n<4;++n)
                            Tetrahedron.push_back(sxyz<T1>(nodes[neighbors[Celli][n]]));
                        std::vector<std::vector<sxyz<T1>>> NewCells;
                        DivideTetraHedron(Tetrahedron, NewCells);
                        for(size_t cel =0;cel<NewCells.size();++cel){
                            if(cel<4){
                                T2 Cell1[4];
                                for(size_t n=1;n<4;++n){
                                    bool found(false);
                                    for(auto no=nodes.begin()+nPrimary-2;no!=nodes.end();++no){
                                        if((*no).getX()==NewCells[cel][n].x&&(*no).getY()==NewCells[cel][n].y&&(*no).getZ()==NewCells[cel][n].z){
                                            Cell1[n]=(*no).getGridIndex();
                                            found=true;
                                            break;
                                        }
                                    }
                                    if(!found){
                                        NODE Newnode(nThreads);
                                        nodes.push_back(Newnode);
                                        nodes.back().setPrimary(10);
                                        nodes.back().setXYZindex(NewCells[cel][n],static_cast<T2>(nodes.size()-1));
                                        T1 Slow;
                                        for (size_t iD=0;iD<4;++iD){
                                            for(size_t n1=1;n1<4;++n1){
                                                if (areCollinear(NewCells[cel][n], neighbors[Celli][iD], neighbors[Celli][(iD+n1)%4])){
                                                    Slow=0.5*(nodes[neighbors[Celli][iD]].getNodeSlowness()+nodes[neighbors[Celli][(iD+n1)%4]].getNodeSlowness());
//                                                    if(nodes[neighbors[Celli][iD]].getprimary()!=10){
//                                                        nodes.back().SetPrincipals(neighbors[Celli][iD]);
//                                                    }else{
   //                                                 std::vector<T2> Principe=nodes[neighbors[Celli][iD]].getPrincipals();
                                                        for (auto P=nodes[neighbors[Celli][iD]].getPrincipals().begin();P!=nodes[neighbors[Celli][iD]].getPrincipals().end();++P)
                                                            nodes.back().SetPrincipals(*P);
//                                                    }
//                                                    
//                                                    if(nodes[neighbors[Celli][(iD+n1)%4]].getprimary()!=10){
//                                                        nodes.back().SetPrincipals(neighbors[Celli][(iD+n1)%4]);
//                                                    }else{
 //                                                       std::vector<T2> Princip=nodes[neighbors[Celli][(iD+n1)%4]].getPrincipals();
                                                        for (auto P=nodes[neighbors[Celli][(iD+n1)%4]].getPrincipals().begin();P!=nodes[neighbors[Celli][(iD+n1)%4]].getPrincipals().end();++P)
                                                            nodes.back().SetPrincipals(*P);
//                                                    }
                                                    break;
                                                }
                                            }
                                        }
                                        nodes.back().setNodeSlowness(Slow);
                                        Cell1[n]=nodes.back().getGridIndex();
                                    }
                                }
                                Cell1[0]=nodes[neighbors[Celli][cel]].getGridIndex();
                                neighbors.push_back( std::vector<T2>(Cell1,Cell1+sizeof(Cell1)/sizeof(T2)));
//                                if (Celli ==CellNO)
//                                    plotCell(neighbors.size()-1, sxyz<T1>(nodes[neighbors.back()[0]]), sxyz<T1>(0.0,0.0,0.0));
                                tetrahedra.push_back(tetrahedronElem<T2>(Cell1[0],Cell1[1],Cell1[2],Cell1[3],tetrahedra[Celli].physical_entity));
                            }
                            if(cel>=4){
                                T2 Cell1[4];
                                for(size_t n=0;n<4;++n){
                                    for(auto no=nodes.begin()+nPrimary-1;no!=nodes.end();++no){
                                        if((*no).getX()==NewCells[cel][n].x&&(*no).getY()==NewCells[cel][n].y&&(*no).getZ()==NewCells[cel][n].z){
                                            Cell1[n]=(*no).getGridIndex();
                                            break;
                                        }
                                    }
                                }
                                neighbors.push_back( std::vector<T2>(Cell1,Cell1+sizeof(Cell1)/sizeof(T2)));
//                                if (Celli ==CellNO)
//                                    plotCell(neighbors.size()-1, sxyz<T1>(nodes[neighbors.back()[0]]), sxyz<T1>(0.0,0.0,0.0));
                                tetrahedra.push_back(tetrahedronElem<T2>(Cell1[0],Cell1[1],Cell1[2],Cell1[3],tetrahedra[Celli].physical_entity));
                            }
                        }
                    }
                }
            }
        }
        for(auto dc=Celli_indices.crbegin();dc!=Celli_indices.crend();++dc){
            neighbors.erase(neighbors.begin()+(*dc));
            tetrahedra.erase(tetrahedra.begin()+(*dc));
        }
        for( auto No=nodes.begin();No!=nodes.end();++No){
            (*No).ReinitialOwner();
        }
        for(size_t nc=0;nc<neighbors.size();++nc){
            for(size_t nn=0;nn<neighbors[nc].size();++nn){
              nodes[neighbors[nc][nn]].pushOwner(nc);
            }
        }
        return 0;
    }
    template<typename T1,typename  T2,typename NODE>
    void Grid3Dui<T1,T2,NODE>:: delateCells(const std::vector<tetrahedronElem<T2>> & delatedCells) const{
        std::set<T2> dalatedNodes;
        std::set<T2> dalatedcells;
        for (size_t n=0;n<neighbors.size();++n){
            bool found=false;
            for(size_t nn=0;nn<4;++nn){
                if(nodes[neighbors[n][nn]].getprimary()==10){
                    found=true;
                    dalatedNodes.insert(neighbors[n][nn]);
                }
            }
            if (found){
                dalatedcells.insert(n);
            }

        }
        for(auto dc=dalatedcells.crbegin();dc!=dalatedcells.crend();++dc){
            neighbors.erase(neighbors.begin()+(*dc));
            tetrahedra.erase(tetrahedra.begin()+(*dc));
        }
        for(auto N=dalatedNodes.crbegin();N!=dalatedNodes.crend();++N){
            nodes.erase(nodes.begin()+(*N));
        }
        for(size_t nc=0;nc<delatedCells.size();++nc){
            neighbors.push_back({delatedCells[nc].i[0],delatedCells[nc].i[1],delatedCells[nc].i[2],delatedCells[nc].i[3]});
            tetrahedra.push_back(delatedCells[nc]);
        }
        for( auto No=nodes.begin();No!=nodes.end();++No){
            (*No).ReinitialOwner();
        }
        for(size_t nc=0;nc<neighbors.size();++nc){
            for(size_t nn=0;nn<neighbors[nc].size();++nn){
                nodes[neighbors[nc][nn]].pushOwner(nc);
            }
        }
    }
    template<typename T1, typename T2, typename NODE>
    sxyz<T1> Grid3Dui<T1,T2,NODE>::reflectedGradient(const sxyz<T1>& g1,const std::array<T2,3> & faces,const T2 & cellNo)const{
        // calculate normal to the face plan
        sxyz<T1> Vect1={nodes[faces[0]].getX()-nodes[faces[1]].getX(),
                        nodes[faces[0]].getY()-nodes[faces[1]].getY(),
                        nodes[faces[0]].getZ()-nodes[faces[1]].getZ()};
        sxyz<T1> Vect2={nodes[faces[0]].getX()-nodes[faces[2]].getX(),
                        nodes[faces[0]].getY()-nodes[faces[2]].getY(),
                        nodes[faces[0]].getZ()-nodes[faces[2]].getZ()};
        sxyz<T1> Normal=cross(Vect1,Vect2);
        sxyz<T1> g1Normal=dot(Normal,g1)*Normal/(norm(Normal)*norm(Normal));
        g1Normal*=1.05;
        return (g1-g1Normal);
    }
    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dui<T1,T2,NODE>::getTraveltime(const sxyz<T1>& Rx,
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
    bool Grid3Dui<T1,T2,NODE>:: areInLimits(const T2 & Node) const {
        if (1.0/nodes[Node].getNodeSlowness()<=Vair)
            return true;
        for(auto nc=nodes[Node].getOwners().begin();nc!=nodes[Node].getOwners().end();++nc){
            for (T2 n=0;n<4;++n){
                if (1.0/nodes[neighbors[*nc][n]].getNodeSlowness()<=Vair)
                    return true;
            }
        }
        return false;
    }

    template<typename T1, typename T2, typename NODE>
    int Grid3Dui<T1,T2,NODE>::checkPts(const std::vector<sxyz<T1>>& pts) const {

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
                std::cerr << "Error: point no " << (n+1)
                << " outside the mesh (" << pts[n].x << "," << pts[n].y
                << "," << pts[n].z << ").\n";
                return 1;
            }
        }
        return 0;
    }



    template<typename T1, typename T2, typename NODE>
    bool Grid3Dui<T1,T2,NODE>::insideTetrahedron(const sxyz<T1>& v, const T2 nt) const {


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
//        if (fabs(D0)-(fabs(D1)+fabs(D2)+fabs(D3)+fabs(D4))<=small)
//            return 1;
        if ( fabs(D1)<small*small) {
            // points are coplanar, check if pt is inside triangle
            if ( testInTriangle(&(nodes[ tetrahedra[nt].i[1] ]),
                                &(nodes[ tetrahedra[nt].i[2] ]),
                                &(nodes[ tetrahedra[nt].i[3] ]), v))
                return 1;
            else
                t1 = (signum(D0)==signum(D1));

        } else {
            t1 = (signum(D0)==signum(D1));
        }

        if ( fabs(D2)<small*small ) {
            if ( testInTriangle(&(nodes[ tetrahedra[nt].i[0] ]),
                                &(nodes[ tetrahedra[nt].i[2] ]),
                                &(nodes[ tetrahedra[nt].i[3] ]), v))
                return 1;
            else
                 t2 = (signum(D0)==signum(D2));
        } else {
            t2 = (signum(D0)==signum(D2));
        }

        if ( fabs(D3)<small*small ) {
            if ( testInTriangle(&(nodes[ tetrahedra[nt].i[0] ]),
                                &(nodes[ tetrahedra[nt].i[1] ]),
                                &(nodes[ tetrahedra[nt].i[3] ]), v))
                return 1;
            else
                t3 = (signum(D0)==signum(D3));
        } else {
            t3 = (signum(D0)==signum(D3));
        }

        if ( fabs(D4)<small*small ) {
            if ( testInTriangle(&(nodes[ tetrahedra[nt].i[0] ]),
                                &(nodes[ tetrahedra[nt].i[1] ]),
                                &(nodes[ tetrahedra[nt].i[2] ]), v))
                return 1;
            else
                 t4 = (signum(D0)==signum(D4));
        } else {
            t4 = (signum(D0)==signum(D4));
        }

        return t1 && t2 && t3 && t4;
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dui<T1,T2,NODE>::saveTT(const std::string &fname, const int all,
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
                fout<< nodes[n].getX() << '\t'
                << nodes[n].getY() << '\t'
                << nodes[n].getZ() << '\t'
                << nodes[n].getTT(nt) << '\n';
            }
            fout.close();
        }
    }

#ifdef VTK
    template<typename T1, typename T2, typename NODE>
    void Grid3Dui<T1,T2,NODE>::saveModelVTU(const std::string &fname,
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
    void Grid3Dui<T1,T2,NODE>::localUpdate3D(NODE *vertexD,
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
            
            if ( vertexA->getTT(threadNo) == std::numeric_limits<T1>::max() ) {
                continue;
            }
            
            T1 s=0.25*(vertexD->getNodeSlowness()+vertexA->getNodeSlowness()
                       +vertexB->getNodeSlowness()+vertexC->getNodeSlowness());
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
                
                sxyz<T1> v_n = cross(v_b, v_c);
                
                T1 b = norm( v_b );
                T1 c = norm( v_c );
                T1 d2 = dot(v_b, v_c);
                
                //T1 alpha = acos( d2 / (b*c) );
                
                T1 phi = norm(v_n);//
                T1 w_tilde2=s*s*phi*phi -u*u*b*b - v*v*c*c + 2.0*u*v*d2;
                if(w_tilde2<0.0){
                    T1 t = localUpdate2D(vertexA, vertexB, vertexD, tetNo, threadNo);
                    if ( t < tABC ) tABC = t;
                    t = localUpdate2D(vertexA, vertexC, vertexD, tetNo, threadNo);
                    if ( t < tABC ) tABC = t;
                    
                    if ( vertexB->getTT(threadNo) > vertexC->getTT(threadNo) ) {
                        t = localUpdate2D(vertexC, vertexB, vertexD, tetNo, threadNo);
                    }else{
                        t = localUpdate2D(vertexB, vertexC, vertexD, tetNo, threadNo);
                    }
                    if ( t < tABC ) tABC = t;
                    if ( tABC<vertexD->getTT(threadNo) )
                        vertexD->setTT(tABC, threadNo);
                    continue;
                }
                T1 w_tilde = sqrt( w_tilde2 );
                
                // project D on plane
                
                T1 d_tmp = -vertexA->getX()*v_n.x - vertexA->getY()*v_n.y - vertexA->getZ()*v_n.z;
                
                T1 k = -(d_tmp + v_n.x*vertexD->getX() + v_n.y*vertexD->getY() + v_n.z*vertexD->getZ())/
                norm2(v_n);
                
                sxyz<T1> pt;
                pt.x = vertexD->getX() + k*v_n.x;
                pt.y = vertexD->getY() + k*v_n.y;
                pt.z = vertexD->getZ() + k*v_n.z;
                
                T1 rho0 = vertexD->getDistance( pt );
                
                sxyz<T1> v_pt = {pt.x-vertexA->getX(), pt.y-vertexA->getY(), pt.z-vertexA->getZ()};
                //// decomposition of Ap
                sxz<T1> AtA_Vect1={b*b,d2};
                sxz<T1> AtA_Vect2={d2,c*c};
                sxz<T1> Atb={dot(v_b,v_pt),dot(v_c,v_pt)};
                T1 DeT=det(AtA_Vect1,AtA_Vect2);
                T1 x=det(AtA_Vect1,Atb)/DeT;
                T1 y=det(Atb,AtA_Vect2)/DeT;
                //               sxyz<T1> Diff=x*v_c+y*v_b-v_pt;
                // project point on AB
                // T1 xi0 = dot(v_pt,v_c)/dot(v_c,v_c);
                T1 xi0=x;
                // project point on AC
                //  T1 zeta0 = dot(v_pt,v_b)/dot(v_b,v_b);
                T1 zeta0 = y;
                
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
            
            //            T1 t = vertexA->getTT(threadNo) + vertexD->getNodeSlowness() * vertexD->getDistance( *vertexA );
            //            if ( t < tABC ) tABC = t;
            //            t = vertexB->getTT(threadNo) + vertexD->getNodeSlowness() * vertexD->getDistance( *vertexB );
            //            if ( t < tABC ) tABC = t;
            //            t = vertexC->getTT(threadNo) + vertexD->getNodeSlowness() * vertexD->getDistance( *vertexC );
            //            if ( t < tABC ) tABC = t;
            
            T1 t = localUpdate2D(vertexA, vertexB, vertexD, tetNo, threadNo);
            if ( t < tABC ) tABC = t;
            t = localUpdate2D(vertexA, vertexC, vertexD, tetNo, threadNo);
            if ( t < tABC ) tABC = t;
            
            if ( vertexB->getTT(threadNo) > vertexC->getTT(threadNo) ) {
                t = localUpdate2D(vertexC, vertexB, vertexD, tetNo, threadNo);
            }else{
                t = localUpdate2D(vertexB, vertexC, vertexD, tetNo, threadNo);
            }
            if ( t < tABC ) tABC = t;
            if ( tABC<vertexD->getTT(threadNo) )
                vertexD->setTT(tABC, threadNo);
            
        }
    }
    
    
    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dui<T1,T2,NODE>::localUpdate2D(const NODE *vertexA,
                                           const NODE *vertexB,
                                           const NODE *vertexC,
                                           const T2 tetNo,
                                           const size_t threadNo) const {
        // method of Lelievre et al. 2011
        
        if ( vertexB->getTT(threadNo)==std::numeric_limits<T1>::max() &&
            vertexA->getTT(threadNo)==std::numeric_limits<T1>::max() ) {
            return std::numeric_limits<T1>::max();
        }
        T1 s=(vertexC->getNodeSlowness()+vertexA->getNodeSlowness()+vertexB->getNodeSlowness())/3.0;
        T1 t=std::numeric_limits<T1>::max();
        
        T1 u = vertexB->getTT(threadNo) - vertexA->getTT(threadNo);
        
        sxyz<T1> v_b = { vertexC->getX() - vertexA->getX(),
            vertexC->getY() - vertexA->getY(),
            vertexC->getZ() - vertexA->getZ() };
        sxyz<T1> v_c = { vertexB->getX() - vertexA->getX(),
            vertexB->getY() - vertexA->getY(),
            vertexB->getZ() - vertexA->getZ() };
        
        T1 c = norm( v_c );
        
        T1 w2 = s*s*c*c - u*u;
        if ( w2 < 0.0 ){
            T1 b=norm(v_b);
            T1 a=vertexC->getDistance(*vertexB);
            T1 tAC=vertexA->getTT(threadNo)+vertexC->getNodeSlowness()*b;
            T1 tBC=vertexB->getTT(threadNo)+vertexC->getNodeSlowness()*a;
            return(tAC<tBC?tAC:tBC);
        }
        
        T1 w = sqrt( w2 );
        
        T1 k = dot(v_b,v_c)/dot(v_c,v_c);
        sxyz<T1> pt;
        pt.x = vertexA->getX() + k*v_c.x;
        pt.y = vertexA->getY() + k*v_c.y;
        pt.z = vertexA->getZ() + k*v_c.z;
        //
        //T1 D=vertexA->getDistance(pt)+vertexB->getDistance(pt)-vertexA->getDistance(*vertexB);
        T1 rho0 = vertexC->getDistance( pt );
        //        T1 xi0 = vertexA->getDistance( pt )/c;
        T1 xi0 = k;
        
        T1 xi = xi0 - u*rho0/(w*c);
        
        if ( 0.0<xi && xi<1.0 ) {
            t = vertexA->getTT(threadNo) + u*xi0 + w*rho0/c;
        } else {
            t = std::numeric_limits<T1>::max();
        }
        T1 b=norm(v_b);
        T1 a=vertexC->getDistance(*vertexB);
        T1 tAC=vertexA->getTT(threadNo)+0.5*(vertexC->getNodeSlowness()+vertexA->getNodeSlowness())*b;
        T1 tBC=vertexB->getTT(threadNo)+0.5*(vertexC->getNodeSlowness()+vertexB->getNodeSlowness())*a;
        t=t<tBC?t:tBC;
        t=t<tAC?t:tAC;
        return t;
    }


    template<typename T1, typename T2, typename NODE>
    void Grid3Dui<T1,T2,NODE>::local3Dsolver(NODE *vertexD,
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
                    // check solution
//                    sxyz<T1> N1 (n[0][0],n[0][1],n[0][2]);
//                    sxyz<T1> N2 (n[1][0],n[1][1],n[1][2]);
//                    sxyz<T1> AB (ab[0],ab[1],ab[2]);
//                    sxyz<T1> AC (ac[0],ac[1],ac[2]);
//                    if (fabs(norm(N1)-1)<small && fabs(dot(AB,N1)-ab[3])<small && fabs(dot(AC,N1)-ac[3])<small){
//                        std::cout<<"norm N1 "<<norm(N1)<<std::endl;
//                    }else{
//                        std::cout<<" bed solution N1"<<std::endl;
//                        std::cout<<"norm N1 "<<norm(N1)<<std::endl;
//                        std::cout<<"AB.N1/c1 "<<fabs(dot(AB,N1)-ab[3])<<std::endl;
//                        std::cout<<"AB.N1/c2 "<<fabs(dot(AC,N1)-ac[3])<<std::endl;
//                    }
//                    if (fabs(norm(N2)-1)<small && fabs(dot(AB,N2)-ab[3])<small && fabs(dot(AC,N2)-ac[3])<small){
//                        std::cout<<"norm N2 "<<norm(N1)<<std::endl;
//                    }else{
//                        std::cout<<" bed solution N2"<<std::endl;
//                        std::cout<<"norm N2 "<<norm(N2)<<std::endl;
//                        std::cout<<"AB.N2/c1 "<<fabs(dot(AB,N2)-ab[3])<<std::endl;
//                        std::cout<<"AB.N2/c2 "<<fabs(dot(AC,N2)-ac[3])<<std::endl;
//                    }

                    for ( size_t ns=0; ns<2; ++ns ) {
                        sxyz<T1> AD ={vertexD->getX()-vertexA->getX(),vertexD->getY()-vertexA->getY(),
                            vertexD->getZ()-vertexA->getZ()};
                        if (dot(sxyz<T1>(n[ns][0],n[ns][1],n[ns][2]),AD)<0){
                            continue;
                        }
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

//                        if ( fabs(AreaTriangle(vertexA,vertexB,E)+AreaTriangle(vertexA,vertexC,E)+
//                            AreaTriangle(vertexC,vertexB,E)-AreaTriangle(vertexA,vertexB,vertexC)<small*small)) {
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

                            if ( fabs(d3-d4)>small*small ) {
                                std::cout << " d3 ne d4: " << d3 << '\t' << d4 << '\t' << d2 << '\n';
                                //plotCell(tetNo, pt, sxyz<T1>(n[ns][0],n[ns][1],n[ns][2]));
                            }

                            T1 t = vertexA->getTT(threadNo) +
                            d3*vertexD->getNodeSlowness();// d4 because (AD.n)*sd d3 in the beginer

                            if ( t<vertexD->getTT(threadNo) )
                                vertexD->setTT(t, threadNo);

                            apply2Dsolvers = false;
                            break;
                        }
                    }
                }
            }

            if ( apply2Dsolvers ) {
                T1 tABD = local2Dsolver(vertexA, vertexB, vertexD, threadNo);
                T1 tACD = local2Dsolver(vertexA, vertexC, vertexD, threadNo);
                
                if ( vertexB->getTT(threadNo) > vertexC->getTT(threadNo) ) {
                    std::swap(iB, iC);
                    std::swap(vertexB, vertexC);
                }
                T1 tBCD = local2Dsolver(vertexB, vertexC, vertexD, threadNo);

                T1 t = tABD < tACD ? tABD : tACD;
                t = t < tBCD ? t : tBCD;

                if ( t<vertexD->getTT(threadNo) )
                    vertexD->setTT(t, threadNo);
            }

        }
    }

    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dui<T1,T2,NODE>::local2Dsolver(const NODE *vertexA,
                                           const NODE *vertexB,
                                           const NODE *vertexC,
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
            
//            T1 gamma = acos((a*a + b*b - c*c)/(2.*a*b));
//            
//            if ( gamma > pi2 ) {
//                std::cout << "*** Obtuse angle: " << gamma*57.2957795 << " ***\n";
//            } else {
//                std::cout << "Accute angle: " << gamma*57.2957795 << " \n";
//            }
            
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
    void Grid3Dui<T1,T2,NODE>::ObtusAngle2Dsolver(const NODE *vertexA,
                                                 const NODE *vertexB,
                                                 const NODE *vertexC,
                                                std::vector<std::vector<T2>> &Vertexes) const{
        Vertexes.push_back({vertexC->getGridIndex(),vertexA->getGridIndex(),vertexB->getGridIndex()});
        bool obtusAngle;
        do{
            size_t Dimension=Vertexes.size();
            for (size_t tri=0;tri<Dimension;++tri){
                T1 a = nodes[Vertexes[tri][0]].getDistance(nodes[Vertexes[tri][2]]);
                T1 b = nodes[Vertexes[tri][1]].getDistance(nodes[Vertexes[tri][0]]);
                T1 c = nodes[Vertexes[tri][1]].getDistance(nodes[Vertexes[tri][2]]);
                T1 gamma = acos((a*a + b*b - c*c)/(2.*a*b));
                obtusAngle=false;
                if(gamma>(pi/2)){
                    //obtusAngle=true;
                    for ( size_t no=0; no<nodes[Vertexes[tri][1]].getOwners().size(); ++no ){
                        T2 tetNo = nodes[Vertexes[tri][1]].getOwners()[no];
                        if (std::find(neighbors[tetNo].begin(),neighbors[tetNo].begin(),nodes[Vertexes[tri][2]].getGridIndex())!=neighbors[tetNo].end()){
                            for(auto Ind=neighbors[tetNo].begin();Ind!=neighbors[tetNo].end();++Ind ){
                                if(nodes[Vertexes[tri][1]].getGridIndex() != *Ind && nodes[Vertexes[tri][2]].getGridIndex() != *Ind){
                                    Vertexes.push_back({Vertexes[tri][0],(*Ind),Vertexes[tri][1]});
                                    Vertexes.push_back({Vertexes[tri][0],(*Ind),Vertexes[tri][2]});
                                }
                            }
                        }
                    }
                    //std::cout<< " Vertexes size "<< Vertexes.size()<<std::endl;
                    Vertexes.erase(Vertexes.begin()+tri);
                    //std::cout<< " Vertexes size "<< Vertexes.size()<<std::endl;
                }
            }
        }while (obtusAngle);
    }

    template<typename T1, typename T2, typename NODE>
    int Grid3Dui<T1,T2,NODE>::solveEq23(const T1 a[], const T1 b[], T1 n[][3]) const {
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
            
            T1 Normn1=sqrt(n[0][0] *n[0][0] +n[0][1] *n[0][1]+n[0][2] *n[0][2]);
            n[0][0] /=Normn1;
            n[0][1] /=Normn1;
            n[0][2] /=Normn1;

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
            T1 Normn2=sqrt(n[1][0] *n[1][0] +n[1][1] *n[1][1]+n[1][2] *n[1][2]);
            n[1][0] /=Normn2;
            n[1][1] /=Normn2;
            n[1][2] /=Normn2;
        }
        return 1;
    }



    template<typename T1, typename T2, typename NODE>
    void Grid3Dui<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
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
                
                // find adjacent cells
                T2 ind[6][2] = {
                    {neighbors[txCell[nt]][0], neighbors[txCell[nt]][1]},
                    {neighbors[txCell[nt]][0], neighbors[txCell[nt]][2]},
                    {neighbors[txCell[nt]][0], neighbors[txCell[nt]][3]},
                    {neighbors[txCell[nt]][1], neighbors[txCell[nt]][2]},
                    {neighbors[txCell[nt]][1], neighbors[txCell[nt]][3]},
                    {neighbors[txCell[nt]][2], neighbors[txCell[nt]][3]} };
                
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
        
        bool InLimits=false;
        bool onNode = false;
        bool onEdge = false;
        bool onFace = false;
        bool onNodePrev = false;
        bool onEdgePrev = false;
        bool onFacePrev = false;
        std::array<T2,2> edgeNodes, edgeNodesPrev;
        std::array<T2,3> faceNodes={{0,0,0}};
        std::array<T2,3> faceNodesPrev;
        Grad3D <T1,NODE> grad3d;
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
            
            T2 ind[6][2] = {
                {neighbors[cellNo][0], neighbors[cellNo][1]},
                {neighbors[cellNo][0], neighbors[cellNo][2]},
                {neighbors[cellNo][0], neighbors[cellNo][3]},
                {neighbors[cellNo][1], neighbors[cellNo][2]},
                {neighbors[cellNo][1], neighbors[cellNo][3]},
                {neighbors[cellNo][2], neighbors[cellNo][3]} };
            
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
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][2] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][2], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][1], neighbors[cellNo][2], neighbors[cellNo][3]} } };
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                
                for ( size_t n=0; n<4; ++n ) {
                    if ( areCoplanar(curr_pt, ind[n][0], ind[n][1], ind[n][2]) ) {
                        onFace = true;
                        faceNodes[0] = ind[n][0];
                        faceNodes[1] = ind[n][1];
                        faceNodes[2] = ind[n][2];
                        break;
                    }
                }
            }
        }
        if (!onEdge && ! onNode && ! onFace){
            onFace=true;
        }
        for(auto t=0;t<txCell.size();++t){
            if (getCellNo( Rx )==txCell[t]){
                r_tmp.emplace_back(Tx[t]);
                reachedTx=true;
                break;
            }
        }
        sxyz<T1> g;
        T2 N=0;
        while ( reachedTx == false && N<150) {
            ++N;

            if ( onNode ) {
                bool foundIntersection = false;
                if (!InLimits){// if the current point is not on the limitis of mesh
                    // find cells common to edge
                    std::set<NODE*> nnodes;
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                        //getNeighborNodes(*nc, nnodes);
                        for ( size_t nn=0; nn<4; ++nn ) {
                            nnodes.insert( &(nodes[ neighbors[*nc][nn]]) );
                        }
                    }

                    // compute gradient with nodes from all common cells
                    sxyz<T1> g = grad3d.ls_grad(nnodes, threadNo,curr_pt);
                    for(size_t nt=0;nt<Tx.size();++nt){
                        for(auto n=neighbors[txCell[nt]].begin();n!=neighbors[txCell[nt]].begin()+4;++n){
                            bool Find=false;
                            for(auto nc=nodes[*n].getOwners().begin();nc!=nodes[*n].getOwners().end();++nc){
                                if(*nc==cellNo){
                                    g=Tx[nt]-curr_pt;
                                    Find=true;
                                    break;
                                }
                            }
                            if (Find)
                                break;
                        }
                        
                    }
                    // find cell for which gradient intersect opposing face
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                        
                        std::array<T2,3> nb;
                        size_t n=0;
                        for (auto nn=neighbors[*nc].begin(); nn!=neighbors[*nc].begin()+4; ++nn ) {
                            if ( *nn != nodeNo ) {
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
                        InLimits=false;
                        r_tmp.push_back( curr_pt );
                        
                        bool break_flag = false;
                        for ( n=0; n<3; ++n ) {
                            if ( nodes[ nb[n] ].getDistance( curr_pt ) < small*small ) {
                                nodeNoPrev = nodeNo;
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;
                                
                                nodeNo = nb[n];
                                onNode = true;
                                onEdge = false;
                                onFace = false;
                                
                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;
                        
                        for ( size_t n1=0; n1<3; ++n1 ) {
                            size_t n2 = (n1+1)%3;
                            if (nodes[nb[n1]].getDistance(curr_pt)+nodes[nb[n2]].getDistance(curr_pt)-nodes[nb[n2]].getDistance(nodes[nb[n1]])<minDist*minDist) {
                                edgeNodesPrev = edgeNodes;
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;
                                
                                edgeNodes[0] = nb[n1];
                                edgeNodes[1] = nb[n2];
                                onNode = false;
                                onEdge = true;
                                onFace = false;
                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;
                        
                        onNodePrev = onNode;
                        onEdgePrev = onEdge;
                        onFacePrev = onFace;
                        onNode = false;
                        onEdge = false;
                        onFace = true;
                        
                        faceNodesPrev = faceNodes;
                        faceNodes = nb;
                        // find next cell
                        T2 PrevCell=cellNo;
                        cellNo = findAdjacentCell1(faceNodes, nodeNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {// cuurent point is in the limit of mesh
                            InLimits=true;
                            cellNo=PrevCell;
                        }
                        break;
                    }
                    
                }else{// cuurent point is in the limit of mesh
                    // we project gradient onto all the bordering faces
                    std::priority_queue<sxyz<T1>,std::vector<sxyz<T1>>,Comparesxyz_norm<T1>> ProjectedGradients;
                    //plotCell(cellNo, curr_pt, g);
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ){
                        std::array<T2,3> InBorderFace;
                        for(T2 nn=0;nn<4;++nn){
                            InBorderFace[0]=neighbors[*nc][(nn+1)%4];
                            InBorderFace[1]=neighbors[*nc][(nn+2)%4];
                            InBorderFace[2]=neighbors[*nc][(nn+3)%4];
                            if(findAdjacentCell2(InBorderFace, *nc, curr_pt)==std::numeric_limits <T2>::max())
                                ProjectedGradients.push((reflectedGradient(g,InBorderFace,*nc)));
                        }
                    }
                    // take the gradient with the biggest norm
                    g=ProjectedGradients.top();
                    while (! ProjectedGradients.empty()){
                        sxyz<T1> g_reflected=ProjectedGradients.top();
                        ProjectedGradients.pop();
                        // find cell for which gradient intersect opposing face
                        for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                            //plotCell(*nc, curr_pt, g_reflected);
                            std::array<T2,3> nb;
                            size_t n=0;
                            for (auto nn=neighbors[*nc].begin(); nn!=neighbors[*nc].begin()+4; ++nn ) {
                                if ( *nn != nodeNo ) {
                                    nb[n++] = *nn;
                                }
                            }
                            std::sort(nb.begin(), nb.end());
                            
                            sxyz<T1> pt_i;
                            foundIntersection = intersectVecTriangle( nodeNo,g_reflected, nb[0], nb[1], nb[2], pt_i);
                            if ( !foundIntersection ) {
                                continue;
                            }
                            
                            prev_pt = curr_pt;
                            curr_pt = pt_i;
                            InLimits=false;
                            r_tmp.push_back( curr_pt );
                            
                            bool break_flag = false;
                            for ( n=0; n<3; ++n ) {
                                if ( nodes[ nb[n] ].getDistance( curr_pt ) < small*small ) {
                                    nodeNoPrev = nodeNo;
                                    onNodePrev = onNode;
                                    onEdgePrev = onEdge;
                                    onFacePrev = onFace;
                                    
                                    nodeNo = nb[n];
                                    onNode = true;
                                    onEdge = false;
                                    onFace = false;
                                    break_flag = true;
                                    break;
                                }
                            }
                            if ( break_flag ) break;
                            
                            for ( size_t n1=0; n1<3; ++n1 ) {
                                size_t n2 = (n1+1)%3;
                                if (nodes[nb[n1]].getDistance(curr_pt)+nodes[nb[n2]].getDistance(curr_pt)-nodes[nb[n2]].getDistance(nodes[nb[n1]])<minDist*minDist) {
                                    edgeNodesPrev = edgeNodes;
                                    onNodePrev = onNode;
                                    onEdgePrev = onEdge;
                                    onFacePrev = onFace;
                                    
                                    edgeNodes[0] = nb[n1];
                                    edgeNodes[1] = nb[n2];
                                    onNode = false;
                                    onEdge = true;
                                    onFace = false;
                                    break_flag = true;
                                    break;
                                }
                            }
                            if ( break_flag ) break;
                            
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            onNode = false;
                            onEdge = false;
                            onFace = true;
                            
                            faceNodesPrev = faceNodes;
                            faceNodes = nb;
                            // find next cell
                            T2 PrevCell=cellNo;
                            cellNo = findAdjacentCell1(faceNodes, nodeNo);
                            if ( cellNo == std::numeric_limits<T2>::max() ) {
                                InLimits=true;
                                cellNo=PrevCell;
                            }
                            break;
                        }
                        
                        if (foundIntersection)
                            break;
                    }
                    
                    if (InLimits && ! foundIntersection){
                        // if the current point is on limits but we don't find intersection, we project again the gradient.
                        foundIntersection=true;
                        
                    }
                    if ( foundIntersection == false && !InLimits) {
                        std::cout << "\n\nWarning: finding raypath on node failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    
                }
                if (!foundIntersection)//
                    InLimits=InLimits?false:true;
                
            } else if ( onEdge ) {
                
                // find cells common to edge
                std::vector<T2> cells;
                std::set<NODE*> nnodes;
                std::vector<sxyz<T1>> Gradients;
                std::array<T2,2> edgeNodes2;
                std::array<T2,2> edgeNodes1;
                for(T2 Edge=0;Edge<2;++Edge){
                    for ( auto nc0=nodes[edgeNodes[Edge]].getOwners().begin(); nc0!=nodes[edgeNodes[Edge]].getOwners().end(); ++nc0 ) {
                        T2 Celli=*nc0;
                        for(T2 n=0;n<4;++n){
                            for(auto nc=nodes[neighbors[Celli][n]].getOwners().begin();nc!=nodes[neighbors[Celli][n]].getOwners().end();++nc){
                                //get cells common to edge
                                for(T2 iD=0;iD<4;++iD){
                                    T2 iA =neighbors[(*nc)][((iD+1)%4)];
                                    T2 iB =neighbors[(*nc)][((iD+2)%4)];
                                    T2 iC =neighbors[(*nc)][((iD+3)%4)];
                                    if ((nodes[iA].getDistance(curr_pt)+nodes[iB].getDistance(curr_pt)-nodes[iA].getDistance(nodes[iB])<minDist*minDist)||
                                        (nodes[iA].getDistance(curr_pt)+nodes[iC].getDistance(curr_pt)-nodes[iA].getDistance(nodes[iC])<minDist*minDist)||
                                        (nodes[iB].getDistance(curr_pt)+nodes[iC].getDistance(curr_pt)-nodes[iB].getDistance(nodes[iC])<minDist*minDist)){
                                        if (std::find(cells.begin(),cells.end(),(*nc))==cells.end())
                                            cells.push_back( (*nc) );
                                        getNeighborNodes(*nc, nnodes);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                if (!InLimits){
                    std::array<T1, 2> Weights={curr_pt.getDistance(nodes[edgeNodes[1]])/nodes[edgeNodes[0]].getDistance(nodes[edgeNodes[1]]),curr_pt.getDistance(nodes[edgeNodes[0]])/nodes[edgeNodes[0]].getDistance(nodes[edgeNodes[1]])};
                    T1 t=nodes[edgeNodes[0]].getTT(threadNo)*Weights[0]+nodes[edgeNodes[1]].getTT(threadNo)*Weights[1];
                    //g= grad3d.RM4D_grad(nnodes, threadNo,curr_pt);
                     //g = grad3d.ls_grad(nnodes, threadNo,curr_pt);
                    g=grad3d.ls_grad(nnodes, threadNo, t, curr_pt);
                    for(size_t nt=0;nt<Tx.size();++nt){
                        for(auto n=neighbors[txCell[nt]].begin();n!=neighbors[txCell[nt]].begin()+4;++n){
                            bool Find=false;
                            for(auto nc=nodes[*n].getOwners().begin();nc!=nodes[*n].getOwners().end();++nc){
                                if(*nc==cellNo){
                                    g=Tx[nt]-curr_pt;
                                    Find=true;
                                    break;
                                }
                            }
                            if (Find)
                                break;
                        }
                        if(curr_pt.getDistance(Tx[nt])<50.0*1.e-3){
                            g=Tx[nt]-curr_pt;
                            break;
                        }
                    }
                    Gradients.push_back(g);
                }else{
                    std::array<T2,3> InBorderFace;
                    for (size_t n=0; n<cells.size(); ++n ) {
                        for ( T2 nn=0; nn<4; ++nn ){
                            InBorderFace[0]=neighbors[cells[n]][(nn+1)%4];
                            InBorderFace[1]=neighbors[cells[n]][(nn+2)%4];
                            InBorderFace[2]=neighbors[cells[n]][(nn+3)%4];
                            if (findAdjacentCell2(InBorderFace, cells[n], curr_pt)==std::numeric_limits<T2>::max() &&
                                testInTriangle(&nodes[InBorderFace[0]], &nodes[InBorderFace[1]], &nodes[InBorderFace[2]], curr_pt)){
                                Gradients.push_back(reflectedGradient(g, InBorderFace,  cells[n]));
                            }
                        }
                    }
                    if (Gradients.size()>1){// take gradient with biggest norm
                        if(norm(Gradients[0])<norm(Gradients[1])){
                            sxyz<T1> Permute=Gradients[0];
                            Gradients[0]=Gradients[1];
                            Gradients[1]=Permute;
                        }
                    }
                }
                bool foundIntersection=false;
                for(T2 ng=0; ng<Gradients.size();++ng){
                    g=Gradients[ng];
                    for (size_t n=0; n<cells.size(); ++n ) {
                        cellNo = cells[n];
                        // there are 2 faces that might be intersected
                        size_t n2=0;
                        size_t n1=0;
                        for ( auto nn=neighbors[cellNo].begin(); nn!= neighbors[cellNo].begin()+4; nn++ ) {
                            if ((*nn)==edgeNodes[0] || (*nn)==edgeNodes[1]||
                                areCollinear(sxyz<T1>(nodes[*nn]), edgeNodes[0], edgeNodes[1])){
                                edgeNodes1[n1++] = *nn;// Edge contains current point
                            }else{
                                edgeNodes2[n2++] = *nn; // opposite edge
                            }
                        }
                        
                        sxyz<T1> pt_i;
                        T2 itmpNode;
                        foundIntersection = intersectVecTriangle(curr_pt, Gradients[ng],
                                                                 edgeNodes1[0],
                                                                 edgeNodes2[0],
                                                                 edgeNodes2[1], pt_i);
                        itmpNode = edgeNodes1[0];
                        if ( !foundIntersection ) {
                            foundIntersection = intersectVecTriangle(curr_pt, Gradients[ng],
                                                                     edgeNodes1[1],
                                                                     edgeNodes2[0],
                                                                     edgeNodes2[1], pt_i);
                            itmpNode = edgeNodes1[1];
                        }
                        if ( !foundIntersection ) {
                            continue;
                        }
                        InLimits=false;
                        prev_pt = curr_pt;
                        curr_pt = pt_i;
                        r_tmp.push_back( curr_pt );
                        bool break_flag = false;
                        for ( size_t n2=0; n2<4; ++n2 ) {
                            if ( nodes[ neighbors[cellNo][n2] ].getDistance( curr_pt ) < small*small ) {
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;
                                edgeNodesPrev = edgeNodes;
                                nodeNo = neighbors[cellNo][n2];
                                onNode = true;
                                onEdge = false;
                                onFace = false;
                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;
                        
                        if ((nodes[itmpNode].getDistance(curr_pt)+nodes[edgeNodes2[0]].getDistance(curr_pt)-nodes[ itmpNode].getDistance(nodes[edgeNodes2[0]])<minDist*minDist)) {
                            edgeNodesPrev = edgeNodes;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            edgeNodes[0] = itmpNode;
                            edgeNodes[1] = edgeNodes2[0];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            
                            break_flag = true;
                            break;
                        } else if((nodes[itmpNode].getDistance(curr_pt)+nodes[edgeNodes2[1]].getDistance(curr_pt)-nodes[ itmpNode].getDistance(nodes[edgeNodes2[1]])<minDist*minDist)){
                            edgeNodesPrev = edgeNodes;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            edgeNodes[0] = itmpNode;
                            edgeNodes[1] = edgeNodes2[1];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            
                            break_flag = true;
                            break;
                        } else if ((nodes[edgeNodes2[1]].getDistance(curr_pt)+nodes[edgeNodes2[0]].getDistance(curr_pt)-nodes[edgeNodes2[1]].getDistance(nodes[edgeNodes2[0]])<minDist*minDist)){
                            edgeNodesPrev = edgeNodes;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            edgeNodes[0] = edgeNodes2[0];
                            edgeNodes[1] = edgeNodes2[1];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            break_flag = true;
                            break;
                        }
                        if ( break_flag ) break;
                        
                        onNodePrev = onNode;
                        onEdgePrev = onEdge;
                        onFacePrev = onFace;
                        onNode = false;
                        onEdge = false;
                        onFace = true;
                        
                        //faceNodesPrev = faceNodes;
                        edgeNodesPrev = edgeNodes;
                        faceNodes[0] = itmpNode;
                        faceNodes[1] = edgeNodes2[0];
                        faceNodes[2] = edgeNodes2[1];
                        std::sort(faceNodes.begin(), faceNodes.end());
                        // find next cell
                        T2 PrevCell(cellNo);
                        cellNo=findAdjacentCell2(faceNodes, cellNo,curr_pt);
                        
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            InLimits=true;
                            cellNo=PrevCell;
                        }
                        break;
                    }
                    if (foundIntersection)
                        break;
                }
                if (InLimits && ! foundIntersection){
                    foundIntersection=true;// we project again th gradient
                    
                }
                if (!foundIntersection){
                    InLimits=InLimits? false:true;
                }
                if ( foundIntersection == false && InLimits==false) {
                    std::cout << "\n\nWarning: finding raypath on edge failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    reachedTx = true;
                }
            } else{ // on Face
                bool foundIntersection=false;
                std::set<NODE*> nnodes;
                if (!InLimits){
                    getNeighborNodes(cellNo, nnodes);
                    T1 t;
                    if (r_tmp.size()<=1){
                        t=Interpolator<T1>::TrilinearTime(curr_pt, nodes[neighbors[cellNo][0]], nodes[neighbors[cellNo][1]], nodes[neighbors[cellNo][2]], nodes[neighbors[cellNo][3]], threadNo);
                    }else{
                        t=Interpolator<T1>::bilinearTime(curr_pt, nodes[faceNodes[0]], nodes[faceNodes[1]],  nodes[faceNodes[2]], threadNo);
                    }
                    //g = grad3d.ls_grad(nodes[neighbors[cellNo][0]], nodes[neighbors[cellNo][1]], nodes[neighbors[cellNo][2]], nodes[neighbors[cellNo][3]], threadNo, curr_pt);
                    //g = grad3d.ls_grad(nnodes, threadNo,curr_pt);
                    g=grad3d.ls_grad(nnodes, threadNo, t, curr_pt);
                    for(size_t nt=0;nt<Tx.size();++nt){
                        for(auto n=neighbors[txCell[nt]].begin();n!=neighbors[txCell[nt]].begin()+4;++n){
                            bool Find=false;
                            for(auto nc=nodes[*n].getOwners().begin();nc!=nodes[*n].getOwners().end();++nc){
                                if(*nc==cellNo){
                                    g=Tx[nt]-curr_pt;
                                    Find=true;
                                    break;
                                }
                            }
                            if (Find)
                                break;
                        }
                        if(curr_pt.getDistance(Tx[nt])<50.0*1.e-3){
                            g=Tx[nt]-curr_pt;
                            break;
                        }
                        
                    }
                }else{
                    g=reflectedGradient(g, faceNodes,cellNo);
                }
                //std::cout<<"nomr g : "<<norm(g)<<std::endl;
                std::array<T2,3> ind[4] = {
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][2] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][2], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][1], neighbors[cellNo][2], neighbors[cellNo][3] } }
                };
                
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                // there are 3 faces that might be intersected
                for ( size_t n=0; n<4; ++n ) {
                    if ( ind[n] == faceNodes) continue;
                    
                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);
                    
                    if ( !foundIntersection )
                        continue;
                    InLimits=false;
                    prev_pt = curr_pt;
                    curr_pt = pt_i;
                    r_tmp.push_back( curr_pt );
                    
                    bool break_flag = false;
                    for ( size_t n2=0; n2<3; ++n2 ) {
                        if ( nodes[ ind[n][n2] ].getDistance( curr_pt ) < small*small ) {
                            //nodeNoPrev = nodeNo;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            nodeNo = ind[n][n2];
                            onNode = true;
                            onEdge = false;
                            onFace = false;
                            faceNodesPrev = faceNodes;
                            
                            break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;
                    
                    for ( size_t n1=0; n1<3; ++n1 ) {
                        size_t n2 = (n1+1)%3;
                        if (nodes[ind[n][n1]].getDistance(curr_pt)+nodes[ind[n][n2]].getDistance(curr_pt)-nodes[ind[n][n1]].getDistance(nodes[ind[n][n2]])<minDist*minDist) {
                            //edgeNodesPrev = edgeNodes;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            edgeNodes[0] = ind[n][n1];
                            edgeNodes[1] = ind[n][n2];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            faceNodesPrev = faceNodes;
                            
                            break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;
                    
                    onNodePrev = onNode;
                    onEdgePrev = onEdge;
                    onFacePrev = onFace;
                    onNode = false;
                    onEdge = false;
                    onFace = true;
                    
                    faceNodesPrev = faceNodes;
                    
                    faceNodes = ind[n];
                    // find next cell
                    T2 PrevCell(cellNo);
                    cellNo=findAdjacentCell2(faceNodes, cellNo,curr_pt);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        InLimits=true;
                        cellNo=PrevCell;
                    }
                    break;
                }
                if ( foundIntersection == false ) {
                    
                    // we must be on an face with gradient pointing slightly outward tetrahedron
                    // return in other cell but keep gradient
                    T2 PrevCell(cellNo);
                    cellNo=findAdjacentCell2(faceNodes, cellNo,curr_pt);
                    if (cellNo==std::numeric_limits<T2>::max())
                        cellNo=PrevCell;
                    std::array<T2,3> ind[4];
                    ind[0] = { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][2] } };
                    ind[1] = { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][3] } };
                    ind[2] = { { neighbors[cellNo][0], neighbors[cellNo][2], neighbors[cellNo][3] } };
                    ind[3] = { { neighbors[cellNo][1], neighbors[cellNo][2], neighbors[cellNo][3] } };
                    
                    for ( size_t n=0; n<4; ++n )
                        std::sort( ind[n].begin(), ind[n].end() );
                    
                    for ( size_t n=0; n<4; ++n ) {
                        if ( ind[n] == faceNodes || AreSameFace(ind[n], faceNodes)) continue;
                        
                        sxyz<T1> pt_i;
                        foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                                 ind[n][1], ind[n][2],
                                                                 pt_i);
                        
                        if ( !foundIntersection ) {
                            continue;
                        }
                        prev_pt = curr_pt;
                        InLimits=false;
                        curr_pt = pt_i;
                        r_tmp.push_back( curr_pt );
                        bool break_flag = false;
                        for ( size_t n2=0; n2<3; ++n2 ) {
                            if ( nodes[ ind[n][n2] ].getDistance( curr_pt ) < small*small ) {
                                //nodeNoPrev = nodeNo;
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;
                                
                                nodeNo = ind[n][n2];
                                onNode = true;
                                onEdge = false;
                                onFace = false;
                                faceNodesPrev = faceNodes;
                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;
                        
                        for ( size_t n1=0; n1<3; ++n1 ) {
                            size_t n2 = (n1+1)%3;
                            if (nodes[ind[n][n1]].getDistance(curr_pt)+nodes[ind[n][n2]].getDistance(curr_pt)-nodes[ind[n][n1]].getDistance(nodes[ind[n][n2]])<minDist*minDist) {
                                //edgeNodesPrev = edgeNodes;
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;
                                
                                edgeNodes[0] = ind[n][n1];
                                edgeNodes[1] = ind[n][n2];
                                onNode = false;
                                onEdge = true;
                                onFace = false;
                                faceNodesPrev = faceNodes;
                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;
                        
                        onNodePrev = onNode;
                        onEdgePrev = onEdge;
                        onFacePrev = onFace;
                        onNode = false;
                        onEdge = false;
                        onFace = true;
                        
                        faceNodesPrev = faceNodes;
                        faceNodes = ind[n];
                        // find next cell
                        T2 PrevCell(cellNo);
                        cellNo=findAdjacentCell2(faceNodes, cellNo,curr_pt);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            InLimits=true;
                            cellNo=PrevCell;
                            
                        }
                        break;
                    }
                }
                if (!foundIntersection)
                    InLimits=InLimits? false:true;
                if ( foundIntersection == false && InLimits==false ) {
                    std::cout << "\n\nWarning: finding raypath on face failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    reachedTx = true;
                }
            }
            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( curr_pt.getDistance( Tx[nt] ) < minDist*minDist ) {
                        reachedTx = true;
                        break;
                    }
                    for(size_t n=0;n<4;++n){
                        sxyz<T1> NearTx={nodes[neighbors[txCell[nt]][n]]};
                        if (curr_pt.getDistance(NearTx)< minDist*minDist){
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
    void Grid3Dui<T1,T2,NODE>::getRaypathABM(const std::vector<sxyz<T1>>& Tx,
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
                
                // find adjacent cells
                T2 ind[6][2] = {
                    {neighbors[txCell[nt]][0], neighbors[txCell[nt]][1]},
                    {neighbors[txCell[nt]][0], neighbors[txCell[nt]][2]},
                    {neighbors[txCell[nt]][0], neighbors[txCell[nt]][3]},
                    {neighbors[txCell[nt]][1], neighbors[txCell[nt]][2]},
                    {neighbors[txCell[nt]][1], neighbors[txCell[nt]][3]},
                    {neighbors[txCell[nt]][2], neighbors[txCell[nt]][3]} };
                
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
        
        bool InLimits=false;
        bool onNode = false;
        bool onEdge = false;
        bool onFace = false;
        bool onNodePrev = false;
        bool onEdgePrev = false;
        bool onFacePrev = false;
        std::array<T2,2> edgeNodes, edgeNodesPrev;
        std::array<T2,3> faceNodes={{0,0,0}};
        std::array<T2,3> faceNodesPrev;
        Grad3D <T1,NODE> grad3d;
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
            
            T2 ind[6][2] = {
                {neighbors[cellNo][0], neighbors[cellNo][1]},
                {neighbors[cellNo][0], neighbors[cellNo][2]},
                {neighbors[cellNo][0], neighbors[cellNo][3]},
                {neighbors[cellNo][1], neighbors[cellNo][2]},
                {neighbors[cellNo][1], neighbors[cellNo][3]},
                {neighbors[cellNo][2], neighbors[cellNo][3]} };
            
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
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][2] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][2], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][1], neighbors[cellNo][2], neighbors[cellNo][3]} } };
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                
                for ( size_t n=0; n<4; ++n ) {
                    if ( areCoplanar(curr_pt, ind[n][0], ind[n][1], ind[n][2]) ) {
                        onFace = true;
                        faceNodes[0] = ind[n][0];
                        faceNodes[1] = ind[n][1];
                        faceNodes[2] = ind[n][2];
                        break;
                    }
                }
            }
        }
        if (!onEdge && ! onNode && ! onFace){
            onFace=true;
        }
        for(auto t=0;t<txCell.size();++t){
            if (getCellNo( Rx )==txCell[t]){
                r_tmp.emplace_back(Tx[t]);
                reachedTx=true;
                break;
            }
        }
        sxyz<T1> g;
        T2 N=0;
        while ( reachedTx == false && N<500) {
//           if (N==66)
//               std::cout<<"arret"<<std::endl;
            ++N;
            //            }
            if ( onNode ) {
                bool foundIntersection = false;
                if (!InLimits){// if the current point is not on the limitis of mesh
                    // compute gradient with nodes from all common cells
                    sxyz<T1> g={0.0,0.0,0.0};
                    T1 sum_wi=0;
                    for (auto Neighborcell=nodes[nodeNo].getOwners().begin();Neighborcell!=nodes[nodeNo].getOwners().end();Neighborcell++){
                        T2 iA(0);
                        T2 iB(1);
                        T2 iC(2);
                        for (T2 n=0;n<4;++n){
                            if (neighbors[(*Neighborcell)][n]==nodeNo){
                                iA=neighbors[(*Neighborcell)][(n+1)%4];
                                iB=neighbors[(*Neighborcell)][(n+2)%4];
                                iC=neighbors[(*Neighborcell)][(n+3)%4];
                                break;
                            }
                        }
                        sxyz<T1>Centroid={0.25*(nodes[nodeNo].getX()+nodes[iA].getX()+nodes[iB].getX()+nodes[iC].getX()),0.25*(nodes[nodeNo].getY()+nodes[iA].getY()+nodes[iB].getY()+nodes[iC].getY()),0.25*(nodes[nodeNo].getZ()+nodes[iA].getZ()+nodes[iB].getZ()+nodes[iC].getZ())};
                        T1 wi=1.0/Centroid.getDistance(nodes[nodeNo]);
                        sum_wi+=wi;
                        g+=wi*grad3d.ABM_grad(nodes[nodeNo], nodes[iA], nodes[iB], nodes[iC], threadNo);
                    }

                    
                    // find cell for which gradient intersect opposing face
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                        
                        std::array<T2,3> nb;
                        size_t n=0;
                        for (auto nn=neighbors[*nc].begin(); nn!=neighbors[*nc].begin()+4; ++nn ) {
                            if ( *nn != nodeNo ) {
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
                        InLimits=false;
                        r_tmp.push_back( curr_pt );
                        
                        bool break_flag = false;
                        for ( n=0; n<3; ++n ) {
                            if ( nodes[ nb[n] ].getDistance( curr_pt ) < small*small ) {
                                nodeNoPrev = nodeNo;
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;
                                
                                nodeNo = nb[n];
                                onNode = true;
                                onEdge = false;
                                onFace = false;
                                
                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;
                        
                        for ( size_t n1=0; n1<3; ++n1 ) {
                            size_t n2 = (n1+1)%3;
                            if (nodes[nb[n1]].getDistance(curr_pt)+nodes[nb[n2]].getDistance(curr_pt)-nodes[nb[n2]].getDistance(nodes[nb[n1]])<minDist*minDist) {
                                edgeNodesPrev = edgeNodes;
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;
                                
                                edgeNodes[0] = nb[n1];
                                edgeNodes[1] = nb[n2];
                                onNode = false;
                                onEdge = true;
                                onFace = false;
                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;
                        
                        onNodePrev = onNode;
                        onEdgePrev = onEdge;
                        onFacePrev = onFace;
                        onNode = false;
                        onEdge = false;
                        onFace = true;
                        
                        faceNodesPrev = faceNodes;
                        faceNodes = nb;
                        // find next cell
                        T2 PrevCell=cellNo;
                        cellNo = findAdjacentCell1(faceNodes, nodeNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {// cuurent point is in the limit of mesh
                            InLimits=true;
                            cellNo=PrevCell;
                        }
                        break;
                    }
                    
                }else{// cuurent point is in the limit of mesh
                    // we project gradient onto all the bordering faces
                    std::priority_queue<sxyz<T1>,std::vector<sxyz<T1>>,Comparesxyz_norm<T1>> ProjectedGradients;
                    //plotCell(cellNo, curr_pt, g);
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ){
                        std::array<T2,3> InBorderFace;
                        for(T2 nn=0;nn<4;++nn){
                            InBorderFace[0]=neighbors[*nc][(nn+1)%4];
                            InBorderFace[1]=neighbors[*nc][(nn+2)%4];
                            InBorderFace[2]=neighbors[*nc][(nn+3)%4];
                            if(findAdjacentCell2(InBorderFace, *nc, curr_pt)==std::numeric_limits <T2>::max())
                                ProjectedGradients.push((reflectedGradient(g,InBorderFace,*nc)));
                        }
                    }
                    // take the gradient with the biggest norm
                    g=ProjectedGradients.top();
                    while (! ProjectedGradients.empty()){
                        sxyz<T1> g_reflected=ProjectedGradients.top();
                        ProjectedGradients.pop();
                        // find cell for which gradient intersect opposing face
                        for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                            //plotCell(*nc, curr_pt, g_reflected);
                            std::array<T2,3> nb;
                            size_t n=0;
                            for (auto nn=neighbors[*nc].begin(); nn!=neighbors[*nc].begin()+4; ++nn ) {
                                if ( *nn != nodeNo ) {
                                    nb[n++] = *nn;
                                }
                            }
                            std::sort(nb.begin(), nb.end());
                            
                            sxyz<T1> pt_i;
                            foundIntersection = intersectVecTriangle( nodeNo,g_reflected, nb[0], nb[1], nb[2], pt_i);
                            if ( !foundIntersection ) {
                                continue;
                            }
                            
                            prev_pt = curr_pt;
                            curr_pt = pt_i;
                            InLimits=false;
                            r_tmp.push_back( curr_pt );
                            
                            bool break_flag = false;
                            for ( n=0; n<3; ++n ) {
                                if ( nodes[ nb[n] ].getDistance( curr_pt ) < small*small ) {
                                    nodeNoPrev = nodeNo;
                                    onNodePrev = onNode;
                                    onEdgePrev = onEdge;
                                    onFacePrev = onFace;
                                    
                                    nodeNo = nb[n];
                                    onNode = true;
                                    onEdge = false;
                                    onFace = false;
                                    break_flag = true;
                                    break;
                                }
                            }
                            if ( break_flag ) break;
                            
                            for ( size_t n1=0; n1<3; ++n1 ) {
                                size_t n2 = (n1+1)%3;
                                if (nodes[nb[n1]].getDistance(curr_pt)+nodes[nb[n2]].getDistance(curr_pt)-nodes[nb[n2]].getDistance(nodes[nb[n1]])<minDist*minDist) {
                                    edgeNodesPrev = edgeNodes;
                                    onNodePrev = onNode;
                                    onEdgePrev = onEdge;
                                    onFacePrev = onFace;
                                    
                                    edgeNodes[0] = nb[n1];
                                    edgeNodes[1] = nb[n2];
                                    onNode = false;
                                    onEdge = true;
                                    onFace = false;
                                    break_flag = true;
                                    break;
                                }
                            }
                            if ( break_flag ) break;
                            
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            onNode = false;
                            onEdge = false;
                            onFace = true;
                            
                            faceNodesPrev = faceNodes;
                            faceNodes = nb;
                            // find next cell
                            T2 PrevCell=cellNo;
                            cellNo = findAdjacentCell1(faceNodes, nodeNo);
                            if ( cellNo == std::numeric_limits<T2>::max() ) {
                                InLimits=true;
                                cellNo=PrevCell;
                            }
                            break;
                        }
                        
                        if (foundIntersection)
                            break;
                    }
                    
                    if (InLimits && ! foundIntersection){
                        // if the current point is on limits but we don't find intersection, we project again the gradient.
                        foundIntersection=true;
                        
                    }
                    if ( foundIntersection == false && !InLimits) {
                        std::cout << "\n\nWarning: finding raypath on node failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    
                }
                if (!foundIntersection)//
                    InLimits=InLimits?false:true;
                
            } else if ( onEdge ) {
                
                // find cells common to edge
                std::vector<T2> cells;
                std::set<NODE*> nnodes;
                std::vector<sxyz<T1>> Gradients;
                std::array<T2,2> edgeNodes2;
                std::array<T2,2> edgeNodes1;
                for(T2 Edge=0;Edge<2;++Edge){
                    for ( auto nc0=nodes[edgeNodes[Edge]].getOwners().begin(); nc0!=nodes[edgeNodes[Edge]].getOwners().end(); ++nc0 ) {
                        T2 Celli=*nc0;
                        for(T2 n=0;n<4;++n){
                            for(auto nc=nodes[neighbors[Celli][n]].getOwners().begin();nc!=nodes[neighbors[Celli][n]].getOwners().end();++nc){
                                //get cells common to edge
                                for(T2 iD=0;iD<4;++iD){
                                    T2 iA =neighbors[(*nc)][((iD+1)%4)];
                                    T2 iB =neighbors[(*nc)][((iD+2)%4)];
                                    T2 iC =neighbors[(*nc)][((iD+3)%4)];
                                    if ((nodes[iA].getDistance(curr_pt)+nodes[iB].getDistance(curr_pt)-nodes[iA].getDistance(nodes[iB])<minDist*minDist)||
                                        (nodes[iA].getDistance(curr_pt)+nodes[iC].getDistance(curr_pt)-nodes[iA].getDistance(nodes[iC])<minDist*minDist)||
                                        (nodes[iB].getDistance(curr_pt)+nodes[iC].getDistance(curr_pt)-nodes[iB].getDistance(nodes[iC])<minDist*minDist)){
                                        if (std::find(cells.begin(),cells.end(),(*nc))==cells.end())
                                            cells.push_back( (*nc) );
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                if (!InLimits){
                    g={0.0,0.0,0.0};
                    std::array<T1, 2> Weights={curr_pt.getDistance(nodes[edgeNodes[1]])/nodes[edgeNodes[0]].getDistance(nodes[edgeNodes[1]]),curr_pt.getDistance(nodes[edgeNodes[0]])/nodes[edgeNodes[0]].getDistance(nodes[edgeNodes[1]])};
                    for (size_t i=0; i<2;++i){
                        T2 Nodei=edgeNodes[i];
                        sxyz<T1> NodeGrad={0.0,0.0,0.0};
                        T1 sum_wi=0;
                        for (auto Neighborcell=nodes[Nodei].getOwners().begin();Neighborcell!=nodes[Nodei].getOwners().end();Neighborcell++){
                            T2 iA(0);
                            T2 iB(1);
                            T2 iC(2);
                            for (T2 n=0;n<4;++n){
                                if (neighbors[(*Neighborcell)][n]==Nodei){
                                    iA=neighbors[(*Neighborcell)][(n+1)%4];
                                    iB=neighbors[(*Neighborcell)][(n+2)%4];
                                    iC=neighbors[(*Neighborcell)][(n+3)%4];
                                    break;
                                }
                            }
                            sxyz<T1>Centroid={0.25*(nodes[Nodei].getX()+nodes[iA].getX()+nodes[iB].getX()+nodes[iC].getX()),0.25*(nodes[Nodei].getY()+nodes[iA].getY()+nodes[iB].getY()+nodes[iC].getY()),0.25*(nodes[Nodei].getZ()+nodes[iA].getZ()+nodes[iB].getZ()+nodes[iC].getZ())};
                            T1 wi=1.0/Centroid.getDistance(nodes[Nodei]);
                            sum_wi+=wi;
                            NodeGrad+=wi*grad3d.ABM_grad(nodes[Nodei], nodes[iA], nodes[iB], nodes[iC], threadNo);
                        }
                        NodeGrad/=sum_wi;
                        g+=Weights[i]*NodeGrad;
                    }
                    for(size_t nt=0;nt<Tx.size();++nt){
                        for(auto n=neighbors[txCell[nt]].begin();n!=neighbors[txCell[nt]].begin()+4;++n){
                            bool Find=false;
                            for(auto nc=nodes[*n].getOwners().begin();nc!=nodes[*n].getOwners().end();++nc){
                                if(*nc==cellNo){
                                    g=Tx[nt]-curr_pt;
                                    Find=true;
                                    break;
                                }
                            }
                            if (Find)
                                break;
                        }
                        if(curr_pt.getDistance(Tx[nt])<10.0*1.e-3){
                            g=Tx[nt]-curr_pt;
                            break;
                        }
                    }
                    Gradients.push_back(g);
                }else{
                    std::array<T2,3> InBorderFace;
                    for (size_t n=0; n<cells.size(); ++n ) {
                        for ( T2 nn=0; nn<4; ++nn ){
                            InBorderFace[0]=neighbors[cells[n]][(nn+1)%4];
                            InBorderFace[1]=neighbors[cells[n]][(nn+2)%4];
                            InBorderFace[2]=neighbors[cells[n]][(nn+3)%4];
                            if (findAdjacentCell2(InBorderFace, cells[n], curr_pt)==std::numeric_limits<T2>::max() &&
                                testInTriangle(&nodes[InBorderFace[0]], &nodes[InBorderFace[1]], &nodes[InBorderFace[2]], curr_pt)){
                                Gradients.push_back(reflectedGradient(g, InBorderFace,  cells[n]));
                            }
                        }
                    }
                    if (Gradients.size()>1){// take gradient with biggest norm
                        if(norm(Gradients[0])<norm(Gradients[1])){
                            sxyz<T1> Permute=Gradients[0];
                            Gradients[0]=Gradients[1];
                            Gradients[1]=Permute;
                        }
                    }
                }
                bool foundIntersection=false;
                for(T2 ng=0; ng<Gradients.size();++ng){
                    g=Gradients[ng];
                    for (size_t n=0; n<cells.size(); ++n ) {
                        cellNo = cells[n];
                        // there are 2 faces that might be intersected
                        size_t n2=0;
                        size_t n1=0;
                        for ( auto nn=neighbors[cellNo].begin(); nn!= neighbors[cellNo].begin()+4; nn++ ) {
                            if ((*nn)==edgeNodes[0] || (*nn)==edgeNodes[1]||
                                areCollinear(sxyz<T1>(nodes[*nn]), edgeNodes[0], edgeNodes[1])){
                                edgeNodes1[n1++] = *nn;// Edge contains current point
                            }else{
                                edgeNodes2[n2++] = *nn; // opposite edge
                            }
                        }
                        
                        sxyz<T1> pt_i;
                        T2 itmpNode;
                        foundIntersection = intersectVecTriangle(curr_pt, Gradients[ng],
                                                                 edgeNodes1[0],
                                                                 edgeNodes2[0],
                                                                 edgeNodes2[1], pt_i);
                        itmpNode = edgeNodes1[0];
                        if ( !foundIntersection ) {
                            foundIntersection = intersectVecTriangle(curr_pt, Gradients[ng],
                                                                     edgeNodes1[1],
                                                                     edgeNodes2[0],
                                                                     edgeNodes2[1], pt_i);
                            itmpNode = edgeNodes1[1];
                        }
                        if ( !foundIntersection ) {
                            continue;
                        }
                        InLimits=false;
                        prev_pt = curr_pt;
                        curr_pt = pt_i;
                        r_tmp.push_back( curr_pt );
                        bool break_flag = false;
                        for ( size_t n2=0; n2<4; ++n2 ) {
                            if ( nodes[ neighbors[cellNo][n2] ].getDistance( curr_pt ) < small*small ) {
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;
                                edgeNodesPrev = edgeNodes;
                                nodeNo = neighbors[cellNo][n2];
                                onNode = true;
                                onEdge = false;
                                onFace = false;
                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;
                        
                        if ((nodes[itmpNode].getDistance(curr_pt)+nodes[edgeNodes2[0]].getDistance(curr_pt)-nodes[ itmpNode].getDistance(nodes[edgeNodes2[0]])<minDist*minDist)) {
                            edgeNodesPrev = edgeNodes;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            edgeNodes[0] = itmpNode;
                            edgeNodes[1] = edgeNodes2[0];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            
                            break_flag = true;
                            break;
                        } else if((nodes[itmpNode].getDistance(curr_pt)+nodes[edgeNodes2[1]].getDistance(curr_pt)-nodes[ itmpNode].getDistance(nodes[edgeNodes2[1]])<minDist*minDist)){
                            edgeNodesPrev = edgeNodes;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            edgeNodes[0] = itmpNode;
                            edgeNodes[1] = edgeNodes2[1];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            
                            break_flag = true;
                            break;
                        } else if ((nodes[edgeNodes2[1]].getDistance(curr_pt)+nodes[edgeNodes2[0]].getDistance(curr_pt)-nodes[edgeNodes2[1]].getDistance(nodes[edgeNodes2[0]])<minDist*minDist)){
                            edgeNodesPrev = edgeNodes;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            edgeNodes[0] = edgeNodes2[0];
                            edgeNodes[1] = edgeNodes2[1];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            break_flag = true;
                            break;
                        }
                        if ( break_flag ) break;
                        
                        onNodePrev = onNode;
                        onEdgePrev = onEdge;
                        onFacePrev = onFace;
                        onNode = false;
                        onEdge = false;
                        onFace = true;
                        
                        //faceNodesPrev = faceNodes;
                        edgeNodesPrev = edgeNodes;
                        faceNodes[0] = itmpNode;
                        faceNodes[1] = edgeNodes2[0];
                        faceNodes[2] = edgeNodes2[1];
                        std::sort(faceNodes.begin(), faceNodes.end());
                        // find next cell
                        T2 PrevCell(cellNo);
                        cellNo=findAdjacentCell2(faceNodes, cellNo,curr_pt);
                        
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            InLimits=true;
                            cellNo=PrevCell;
                        }
                        break;
                    }
                    if (foundIntersection)
                        break;
                }
                if (InLimits && ! foundIntersection){
                    foundIntersection=true;// we project again th gradient
                    
                }
                if (!foundIntersection){
                    InLimits=InLimits? false:true;
                }
                if ( foundIntersection == false && InLimits==false) {
                    std::cout << "\n\nWarning: finding raypath on edge failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    reachedTx = true;
                }
            } else{ // on Face
                bool foundIntersection=false;
                std::set<NODE*> nnodes;
                if (!InLimits){
                    g={0.0,0.0,0.0};
                    if (r_tmp.size()==1){
                        std::array<T1, 4> Weights;
                        Interpolator<T1>::TrilinearTriangleWeights(curr_pt, nodes[neighbors[cellNo][0]], nodes[neighbors[cellNo][1]], nodes[neighbors[cellNo][2]], nodes[neighbors[cellNo][3]], Weights);
                        for (size_t i=0; i<4;++i){
                            T2 Nodei=neighbors[cellNo][i];
                            sxyz<T1> NodeGrad={0.0,0.0,0.0};
                            T1 sum_wi=0;
                            for (auto Neighborcell=nodes[Nodei].getOwners().begin();Neighborcell!=nodes[Nodei].getOwners().end();Neighborcell++){
                                T2 iA(0);
                                T2 iB(1);
                                T2 iC(2);
                                for (T2 n=0;n<4;++n){
                                    if (neighbors[(*Neighborcell)][n]==Nodei){
                                        iA=neighbors[(*Neighborcell)][(n+1)%4];
                                        iB=neighbors[(*Neighborcell)][(n+2)%4];
                                        iC=neighbors[(*Neighborcell)][(n+3)%4];
                                        break;
                                    }
                                }
                                sxyz<T1>Centroid={0.25*(nodes[Nodei].getX()+nodes[iA].getX()+nodes[iB].getX()+nodes[iC].getX()),0.25*(nodes[Nodei].getY()+nodes[iA].getY()+nodes[iB].getY()+nodes[iC].getY()),0.25*(nodes[Nodei].getZ()+nodes[iA].getZ()+nodes[iB].getZ()+nodes[iC].getZ())};
                                T1 wi=1.0/Centroid.getDistance(nodes[Nodei]);
                                sum_wi+=wi;
                                NodeGrad+=wi*grad3d.ABM_grad(nodes[Nodei], nodes[iA], nodes[iB], nodes[iC], threadNo);
                            }
                            NodeGrad/=sum_wi;
                            g+=Weights[i]*NodeGrad;
                        }
                    }else{
                        std::array<T1, 3> Weights;
                        Interpolator<T1>::bilinearTriangleWeight(curr_pt, nodes[faceNodes[0]], nodes[faceNodes[1]], nodes[faceNodes[2]], Weights);
                        for (size_t i=0; i<3;++i){
                            T2 Nodei=faceNodes[i];
                            sxyz<T1> NodeGrad={0.0,0.0,0.0};
                            T1 sum_wi=0;
                            for (auto Neighborcell=nodes[Nodei].getOwners().begin();Neighborcell!=nodes[Nodei].getOwners().end();Neighborcell++){
                                T2 iA(0);
                                T2 iB(1);
                                T2 iC(2);
                                for (T2 n=0;n<4;++n){
                                    if (neighbors[(*Neighborcell)][n]==Nodei){
                                        iA=neighbors[(*Neighborcell)][(n+1)%4];
                                        iB=neighbors[(*Neighborcell)][(n+2)%4];
                                        iC=neighbors[(*Neighborcell)][(n+3)%4];
                                        break;
                                    }
                                }
                                sxyz<T1>Centroid={0.25*(nodes[Nodei].getX()+nodes[iA].getX()+nodes[iB].getX()+nodes[iC].getX()),0.25*(nodes[Nodei].getY()+nodes[iA].getY()+nodes[iB].getY()+nodes[iC].getY()),0.25*(nodes[Nodei].getZ()+nodes[iA].getZ()+nodes[iB].getZ()+nodes[iC].getZ())};
                                T1 wi=1.0/Centroid.getDistance(nodes[Nodei]);
                                sum_wi+=wi;
                                NodeGrad+=wi*grad3d.ABM_grad(nodes[Nodei], nodes[iA], nodes[iB], nodes[iC], threadNo);
                                }
                            NodeGrad/=sum_wi;
                            g+=Weights[i]*NodeGrad;
                            }
                        }
                    for(size_t nt=0;nt<Tx.size();++nt){
                        for(auto n=neighbors[txCell[nt]].begin();n!=neighbors[txCell[nt]].begin()+4;++n){
                            bool Find=false;
                            for(auto nc=nodes[*n].getOwners().begin();nc!=nodes[*n].getOwners().end();++nc){
                                if(*nc==cellNo){
                                    g=Tx[nt]-curr_pt;
                                    Find=true;
                                    break;
                                }
                            }
                            if (Find)
                                break;
                        }
                        if(curr_pt.getDistance(Tx[nt])<35.0*1.e-3){
                            g=Tx[nt]-curr_pt;
                            break;
                        }
                        
                    }
                }else{
                    g=reflectedGradient(g, faceNodes,cellNo);
                }
                //std::cout<<"nomr g : "<<norm(g)<<std::endl;
                std::array<T2,3> ind[4] = {
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][2] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][2], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][1], neighbors[cellNo][2], neighbors[cellNo][3] } }
                };
                
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                // there are 3 faces that might be intersected
                for ( size_t n=0; n<4; ++n ) {
                    if ( ind[n] == faceNodes ) continue;
                    
                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);
                    
                    if ( !foundIntersection )
                        continue;
                    InLimits=false;
                    prev_pt = curr_pt;
                    curr_pt = pt_i;
                    r_tmp.push_back( curr_pt );
                    
                    bool break_flag = false;
                    for ( size_t n2=0; n2<3; ++n2 ) {
                        if ( nodes[ ind[n][n2] ].getDistance( curr_pt ) < small*small ) {
                            //nodeNoPrev = nodeNo;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            nodeNo = ind[n][n2];
                            onNode = true;
                            onEdge = false;
                            onFace = false;
                            faceNodesPrev = faceNodes;
                            
                            break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;
                    
                    for ( size_t n1=0; n1<3; ++n1 ) {
                        size_t n2 = (n1+1)%3;
                        if (nodes[ind[n][n1]].getDistance(curr_pt)+nodes[ind[n][n2]].getDistance(curr_pt)-nodes[ind[n][n1]].getDistance(nodes[ind[n][n2]])<minDist*minDist) {
                            //edgeNodesPrev = edgeNodes;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            edgeNodes[0] = ind[n][n1];
                            edgeNodes[1] = ind[n][n2];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            faceNodesPrev = faceNodes;
                            
                            break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;
                    
                    onNodePrev = onNode;
                    onEdgePrev = onEdge;
                    onFacePrev = onFace;
                    onNode = false;
                    onEdge = false;
                    onFace = true;
                    
                    faceNodesPrev = faceNodes;
                    
                    faceNodes = ind[n];
                    // find next cell
                    T2 PrevCell(cellNo);
                    cellNo=findAdjacentCell2(faceNodes, cellNo,curr_pt);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        InLimits=true;
                        cellNo=PrevCell;
                    }
                    break;
                }
                if ( foundIntersection == false ) {
                    
                    // we must be on an face with gradient pointing slightly outward tetrahedron
                    // return in other cell but keep gradient
                    T2 PrevCell(cellNo);
                    cellNo=findAdjacentCell2(faceNodes, cellNo,curr_pt);
                    if (cellNo==std::numeric_limits<T2>::max())
                        cellNo=PrevCell;
                    std::array<T2,3> ind[4];
                    ind[0] = { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][2] } };
                    ind[1] = { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][3] } };
                    ind[2] = { { neighbors[cellNo][0], neighbors[cellNo][2], neighbors[cellNo][3] } };
                    ind[3] = { { neighbors[cellNo][1], neighbors[cellNo][2], neighbors[cellNo][3] } };
                    
                    for ( size_t n=0; n<4; ++n )
                        std::sort( ind[n].begin(), ind[n].end() );
                    
                    for ( size_t n=0; n<4; ++n ) {
                        if ( ind[n] == faceNodes || AreSameFace(ind[n], faceNodes)) continue;
                        
                        sxyz<T1> pt_i;
                        foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                                 ind[n][1], ind[n][2],
                                                                 pt_i);
                        
                        if ( !foundIntersection ) {
                            continue;
                        }
                        prev_pt = curr_pt;
                        InLimits=false;
                        curr_pt = pt_i;
                        r_tmp.push_back( curr_pt );
                        bool break_flag = false;
                        for ( size_t n2=0; n2<3; ++n2 ) {
                            if ( nodes[ ind[n][n2] ].getDistance( curr_pt ) < small*small ) {
                                //nodeNoPrev = nodeNo;
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;
                                
                                nodeNo = ind[n][n2];
                                onNode = true;
                                onEdge = false;
                                onFace = false;
                                faceNodesPrev = faceNodes;
                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;
                        
                        for ( size_t n1=0; n1<3; ++n1 ) {
                            size_t n2 = (n1+1)%3;
                            if (nodes[ind[n][n1]].getDistance(curr_pt)+nodes[ind[n][n2]].getDistance(curr_pt)-nodes[ind[n][n1]].getDistance(nodes[ind[n][n2]])<minDist*minDist) {
                                //edgeNodesPrev = edgeNodes;
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;
                                
                                edgeNodes[0] = ind[n][n1];
                                edgeNodes[1] = ind[n][n2];
                                onNode = false;
                                onEdge = true;
                                onFace = false;
                                faceNodesPrev = faceNodes;
                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;
                        
                        onNodePrev = onNode;
                        onEdgePrev = onEdge;
                        onFacePrev = onFace;
                        onNode = false;
                        onEdge = false;
                        onFace = true;
                        
                        faceNodesPrev = faceNodes;
                        faceNodes = ind[n];
                        // find next cell
                        T2 PrevCell(cellNo);
                        cellNo=findAdjacentCell2(faceNodes, cellNo,curr_pt);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            InLimits=true;
                            cellNo=PrevCell;
                            
                        }
                        break;
                    }
                }
                if (!foundIntersection)
                    InLimits=InLimits? false:true;
                if ( foundIntersection == false && InLimits==false ) {
                    std::cout << "\n\nWarning: finding raypath on face failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    reachedTx = true;
                }
            }
            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( curr_pt.getDistance( Tx[nt] ) < minDist*minDist ) {
                        reachedTx = true;
                        break;
                    }
                    for(size_t n=0;n<4;++n){
                        sxyz<T1> NearTx={nodes[neighbors[txCell[nt]][n]]};
                        if (curr_pt.getDistance(NearTx)< minDist*minDist){
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
    void Grid3Dui<T1,T2,NODE>::getRaypathBLIT(const std::vector<sxyz<T1>>& Tx,
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
                
                // find adjacent cells
                T2 ind[6][2] = {
                    {neighbors[txCell[nt]][0], neighbors[txCell[nt]][1]},
                    {neighbors[txCell[nt]][0], neighbors[txCell[nt]][2]},
                    {neighbors[txCell[nt]][0], neighbors[txCell[nt]][3]},
                    {neighbors[txCell[nt]][1], neighbors[txCell[nt]][2]},
                    {neighbors[txCell[nt]][1], neighbors[txCell[nt]][3]},
                    {neighbors[txCell[nt]][2], neighbors[txCell[nt]][3]} };
                
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
        Grad3D <T1,NODE> grad3d;
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
            
            T2 ind[6][2] = {
                {neighbors[cellNo][0], neighbors[cellNo][1]},
                {neighbors[cellNo][0], neighbors[cellNo][2]},
                {neighbors[cellNo][0], neighbors[cellNo][3]},
                {neighbors[cellNo][1], neighbors[cellNo][2]},
                {neighbors[cellNo][1], neighbors[cellNo][3]},
                {neighbors[cellNo][2], neighbors[cellNo][3]} };
            
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
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][2] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][2], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][1], neighbors[cellNo][2], neighbors[cellNo][3]} } };
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                
                for ( size_t n=0; n<4; ++n ) {
                    if ( areCoplanar(curr_pt, ind[n][0], ind[n][1], ind[n][2]) ) {
                        onFace = true;
                        faceNodes[0] = ind[n][0];
                        faceNodes[1] = ind[n][1];
                        faceNodes[2] = ind[n][2];
                        break;
                    }
                }
            }
        }
         T1 time=std::numeric_limits<T1>::max();
        if (!onEdge && ! onNode && ! onFace){
            onFace=true;
            time=Interpolator<T1>::TrilinearTime(curr_pt, nodes[neighbors[cellNo][0]], nodes[neighbors[cellNo][1]], nodes[neighbors[cellNo][2]], nodes[neighbors[cellNo][3]], threadNo);
        }
        for(auto nt=0;nt<txCell.size();++nt){
            if (getCellNo( Rx )==txCell[nt]){
                r_tmp.emplace_back(Tx[nt]);
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
                        if (BLTISolver_ArroundSource(NodeSource, curr_pt, ind, Barycenter)==true){
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
                            if (BLTISolver_ArroundSource(NodeSource, curr_pt, ind[n], Barycenter)==true){
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
                        if (BLTISolver_ArroundSource(NodeSource, curr_pt, ind[n], Barycenter)==true){
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
                        reachedTx = true;
                        break;
                    }
                    for(size_t n=0;n<4;++n){
                        sxyz<T1> NearTx={nodes[neighbors[txCell[nt]][n]]};
                        if (curr_pt.getDistance(NearTx)< minDist){
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
    void Grid3Dui<T1,T2,NODE>::getRaypath_ho(const std::vector<sxyz<T1>>& Tx,
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
                
                // find adjacent cells
                T2 ind[6][2] = {
                    {neighbors[txCell[nt]][0], neighbors[txCell[nt]][1]},
                    {neighbors[txCell[nt]][0], neighbors[txCell[nt]][2]},
                    {neighbors[txCell[nt]][0], neighbors[txCell[nt]][3]},
                    {neighbors[txCell[nt]][1], neighbors[txCell[nt]][2]},
                    {neighbors[txCell[nt]][1], neighbors[txCell[nt]][3]},
                    {neighbors[txCell[nt]][2], neighbors[txCell[nt]][3]} };
                
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

        bool InLimits=false;
        bool onNode = false;
        bool onEdge = false;
        bool onFace = false;
        bool onNodePrev = false;
        bool onEdgePrev = false;
        bool onFacePrev = false;
        std::array<T2,2> edgeNodes, edgeNodesPrev;
        std::array<T2,3> faceNodes={{0,0,0}};
        std::array<T2,3> faceNodesPrev;
        Grad3D_ho <T1,NODE> grad3d;
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
            
            T2 ind[6][2] = {
                {neighbors[cellNo][0], neighbors[cellNo][1]},
                {neighbors[cellNo][0], neighbors[cellNo][2]},
                {neighbors[cellNo][0], neighbors[cellNo][3]},
                {neighbors[cellNo][1], neighbors[cellNo][2]},
                {neighbors[cellNo][1], neighbors[cellNo][3]},
                {neighbors[cellNo][2], neighbors[cellNo][3]} };
            
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
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][2] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][2], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][1], neighbors[cellNo][2], neighbors[cellNo][3] } }
                };
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                
                for ( size_t n=0; n<4; ++n ) {
                    if ( areCoplanar(curr_pt, ind[n][0], ind[n][1], ind[n][2]) ) {
                        onFace = true;
                        faceNodes[0] = ind[n][0];
                        faceNodes[1] = ind[n][1];
                        faceNodes[2] = ind[n][2];
                        break;
                    }
                }
            }
        }
        if ( !onNode ) {
            cellNo = getCellNo( curr_pt );
            
            T2 ind[6][2] = {
                {neighbors[cellNo][0], neighbors[cellNo][1]},
                {neighbors[cellNo][0], neighbors[cellNo][2]},
                {neighbors[cellNo][0], neighbors[cellNo][3]},
                {neighbors[cellNo][1], neighbors[cellNo][2]},
                {neighbors[cellNo][1], neighbors[cellNo][3]},
                {neighbors[cellNo][2], neighbors[cellNo][3]} };
            
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
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][2] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][2], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][1], neighbors[cellNo][2], neighbors[cellNo][3]} } };
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                
                for ( size_t n=0; n<4; ++n ) {
                    if ( areCoplanar(curr_pt, ind[n][0], ind[n][1], ind[n][2]) ) {
                        onFace = true;
                        faceNodes[0] = ind[n][0];
                        faceNodes[1] = ind[n][1];
                        faceNodes[2] = ind[n][2];
                        break;
                    }
                }
            }
        }
        if (!onEdge && ! onNode && ! onFace){
            onFace=true;
        }
        for(auto t=0;t<txCell.size();++t){
            if (getCellNo( Rx )==txCell[t]){
                r_tmp.emplace_back(Tx[t]);
                reachedTx=true;
                break;
            }
        }
        sxyz<T1> g;
        T2 N=0;
        while ( reachedTx == false && N<250) {
            ++N;
            //            }
            if ( onNode ) {
                bool foundIntersection = false;
                if (!InLimits){// if the current point is not on the limitis of mesh
                    // find cells common to edge
                    std::set<NODE*> nnodes;
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                        getNeighborNodes(*nc, nnodes);
                    }
                    // compute gradient with nodes from all common cells
                    sxyz<T1> g = grad3d.ls_grad(nnodes, threadNo,curr_pt);
                    
                    // find cell for which gradient intersect opposing face
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                        
                        std::array<T2,3> nb;
                        size_t n=0;
                        for (auto nn=neighbors[*nc].begin(); nn!=neighbors[*nc].begin()+4; ++nn ) {
                            if ( *nn != nodeNo ) {
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
                        InLimits=false;
                        r_tmp.push_back( curr_pt );
                        
                        bool break_flag = false;
                        for ( n=0; n<3; ++n ) {
                            if ( nodes[ nb[n] ].getDistance( curr_pt ) < small*small ) {
                                nodeNoPrev = nodeNo;
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;
                                
                                nodeNo = nb[n];
                                onNode = true;
                                onEdge = false;
                                onFace = false;

                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;
                        
                        for ( size_t n1=0; n1<3; ++n1 ) {
                            size_t n2 = (n1+1)%3;
                            if (nodes[nb[n1]].getDistance(curr_pt)+nodes[nb[n2]].getDistance(curr_pt)-nodes[nb[n2]].getDistance(nodes[nb[n1]])<minDist*minDist) {
                                edgeNodesPrev = edgeNodes;
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;
                                
                                edgeNodes[0] = nb[n1];
                                edgeNodes[1] = nb[n2];
                                onNode = false;
                                onEdge = true;
                                onFace = false;
                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;
                        
                        onNodePrev = onNode;
                        onEdgePrev = onEdge;
                        onFacePrev = onFace;
                        onNode = false;
                        onEdge = false;
                        onFace = true;
                        
                        faceNodesPrev = faceNodes;
                        faceNodes = nb;
                        // find next cell
                        T2 PrevCell=cellNo;
                        cellNo = findAdjacentCell1(faceNodes, nodeNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {// cuurent point is in the limit of mesh
                            InLimits=true;
                            cellNo=PrevCell;
                        }
                        break;
                    }
                    
                }else{// cuurent point is in the limit of mesh
                    // we project gradient onto all the bordering faces
                    std::priority_queue<sxyz<T1>,std::vector<sxyz<T1>>,Comparesxyz_norm<T1>> ProjectedGradients;
                    //plotCell(cellNo, curr_pt, g);
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ){
                        std::array<T2,3> InBorderFace;
                        for(T2 nn=0;nn<4;++nn){
                            InBorderFace[0]=neighbors[*nc][(nn+1)%4];
                            InBorderFace[1]=neighbors[*nc][(nn+2)%4];
                            InBorderFace[2]=neighbors[*nc][(nn+3)%4];
                            if(findAdjacentCell2(InBorderFace, *nc, curr_pt)==std::numeric_limits <T2>::max())
                                ProjectedGradients.push((reflectedGradient(g,InBorderFace,*nc)));
                        }
                    }
                    // take the gradient with the biggest norm
                    g=ProjectedGradients.top();
                    while (! ProjectedGradients.empty()){
                        sxyz<T1> g_reflected=ProjectedGradients.top();
                        ProjectedGradients.pop();
                        // find cell for which gradient intersect opposing face
                        for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                            //plotCell(*nc, curr_pt, g_reflected);
                            std::array<T2,3> nb;
                            size_t n=0;
                            for (auto nn=neighbors[*nc].begin(); nn!=neighbors[*nc].begin()+4; ++nn ) {
                                if ( *nn != nodeNo ) {
                                    nb[n++] = *nn;
                                }
                            }
                            std::sort(nb.begin(), nb.end());
                            
                            sxyz<T1> pt_i;
                            foundIntersection = intersectVecTriangle( nodeNo,g_reflected, nb[0], nb[1], nb[2], pt_i);
                            if ( !foundIntersection ) {
                                continue;
                            }
                            
                            prev_pt = curr_pt;
                            curr_pt = pt_i;
                            InLimits=false;
                            r_tmp.push_back( curr_pt );
                            
                            bool break_flag = false;
                            for ( n=0; n<3; ++n ) {
                                if ( nodes[ nb[n] ].getDistance( curr_pt ) < small*small ) {
                                    nodeNoPrev = nodeNo;
                                    onNodePrev = onNode;
                                    onEdgePrev = onEdge;
                                    onFacePrev = onFace;
                                    
                                    nodeNo = nb[n];
                                    onNode = true;
                                    onEdge = false;
                                    onFace = false;
                                    break_flag = true;
                                    break;
                                }
                            }
                            if ( break_flag ) break;
                            
                            for ( size_t n1=0; n1<3; ++n1 ) {
                                size_t n2 = (n1+1)%3;
                                if (nodes[nb[n1]].getDistance(curr_pt)+nodes[nb[n2]].getDistance(curr_pt)-nodes[nb[n2]].getDistance(nodes[nb[n1]])<minDist*minDist) {
                                    edgeNodesPrev = edgeNodes;
                                    onNodePrev = onNode;
                                    onEdgePrev = onEdge;
                                    onFacePrev = onFace;
                                    
                                    edgeNodes[0] = nb[n1];
                                    edgeNodes[1] = nb[n2];
                                    onNode = false;
                                    onEdge = true;
                                    onFace = false;
                                    break_flag = true;
                                    break;
                                }
                            }
                            if ( break_flag ) break;
                            
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            onNode = false;
                            onEdge = false;
                            onFace = true;
                            
                            faceNodesPrev = faceNodes;
                            faceNodes = nb;
                            // find next cell
                            T2 PrevCell=cellNo;
                            cellNo = findAdjacentCell1(faceNodes, nodeNo);
                            if ( cellNo == std::numeric_limits<T2>::max() ) {
                                InLimits=true;
                                cellNo=PrevCell;
                            }
                            break;
                        }
                        
                        if (foundIntersection)
                            break;
                    }
                    
                    if (InLimits && ! foundIntersection){
                        // if the current point is on limits but we don't find intersection, we project again the gradient.
                        foundIntersection=true;
                        
                    }
                    if ( foundIntersection == false && !InLimits) {
                        std::cout << "\n\nWarning: finding raypath on node failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    
                }
                if (!foundIntersection)//
                    InLimits=InLimits?false:true;
                
            } else if ( onEdge ) {
                
                // find cells common to edge
                std::vector<T2> cells;
                std::set<NODE*> nnodes;
                std::vector<sxyz<T1>> Gradients;
                std::array<T2,2> edgeNodes2;
                std::array<T2,2> edgeNodes1;
                for(T2 Edge=0;Edge<2;++Edge){
                    for ( auto nc0=nodes[edgeNodes[Edge]].getOwners().begin(); nc0!=nodes[edgeNodes[Edge]].getOwners().end(); ++nc0 ) {
                        T2 Celli=*nc0;
                        for(T2 n=0;n<4;++n){
                            for(auto nc=nodes[neighbors[Celli][n]].getOwners().begin();nc!=nodes[neighbors[Celli][n]].getOwners().end();++nc){
                                //get cells common to edge
                                for(T2 iD=0;iD<4;++iD){
                                    T2 iA =neighbors[(*nc)][((iD+1)%4)];
                                    T2 iB =neighbors[(*nc)][((iD+2)%4)];
                                    T2 iC =neighbors[(*nc)][((iD+3)%4)];
                                    if ((nodes[iA].getDistance(curr_pt)+nodes[iB].getDistance(curr_pt)-nodes[iA].getDistance(nodes[iB])<minDist*minDist)||
                                        (nodes[iA].getDistance(curr_pt)+nodes[iC].getDistance(curr_pt)-nodes[iA].getDistance(nodes[iC])<minDist*minDist)||
                                        (nodes[iB].getDistance(curr_pt)+nodes[iC].getDistance(curr_pt)-nodes[iB].getDistance(nodes[iC])<minDist*minDist)){
                                        if (std::find(cells.begin(),cells.end(),(*nc))==cells.end())
                                            cells.push_back( (*nc) );
                                        getNeighborNodes(*nc, nnodes);
                                    }
                                }
                            }
                        }
                    }
                }
                if (!InLimits){
                    g= grad3d.ls_grad(nnodes, threadNo,curr_pt);
                    for(size_t nt=0;nt<Tx.size();++nt){
                        for(auto n=neighbors[txCell[nt]].begin();n!=neighbors[txCell[nt]].begin()+4;++n){
                            bool Find=false;
                            for(auto nc=nodes[*n].getOwners().begin();nc!=nodes[*n].getOwners().end();++nc){
                                if(*nc==cellNo){
                                    g=Tx[nt]-curr_pt;
                                    Find=true;
                                    break;
                                }
                            }
                            if (Find)
                                break;
                        }
                        if(curr_pt.getDistance(Tx[nt])<50.0*1.e-3){
                            g=Tx[nt]-curr_pt;
                            break;
                        }
                    }
                    Gradients.push_back(g);
                }else{
                    std::array<T2,3> InBorderFace;
                    for (size_t n=0; n<cells.size(); ++n ) {
                        for ( T2 nn=0; nn<4; ++nn ){
                            InBorderFace[0]=neighbors[cells[n]][(nn+1)%4];
                            InBorderFace[1]=neighbors[cells[n]][(nn+2)%4];
                            InBorderFace[2]=neighbors[cells[n]][(nn+3)%4];
                            if (findAdjacentCell2(InBorderFace, cells[n], curr_pt)==std::numeric_limits<T2>::max() &&
                                testInTriangle(&nodes[InBorderFace[0]], &nodes[InBorderFace[1]], &nodes[InBorderFace[2]], curr_pt)){
                                Gradients.push_back(reflectedGradient(g, InBorderFace,  cells[n]));
                            }
                        }
                    }
                    if (Gradients.size()>1){// take gradient with biggest norm
                        if(norm(Gradients[0])<norm(Gradients[1])){
                            sxyz<T1> Permute=Gradients[0];
                            Gradients[0]=Gradients[1];
                            Gradients[1]=Permute;
                        }
                    }
                }
                bool foundIntersection=false;
                for(T2 ng=0; ng<Gradients.size();++ng){
                     g=Gradients[ng];
                    for (size_t n=0; n<cells.size(); ++n ) {
                        cellNo = cells[n];
                        // there are 2 faces that might be intersected
                        size_t n2=0;
                        size_t n1=0;
                        for ( auto nn=neighbors[cellNo].begin(); nn!= neighbors[cellNo].begin()+4; nn++ ) {
                            if ((*nn)==edgeNodes[0] || (*nn)==edgeNodes[1]||
                                areCollinear(sxyz<T1>(nodes[*nn]), edgeNodes[0], edgeNodes[1])){
                                edgeNodes1[n1++] = *nn;// Edge contains current point
                            }else{
                                edgeNodes2[n2++] = *nn; // opposite edge
                            }
                        }
                        
                        sxyz<T1> pt_i;
                        T2 itmpNode;
                        foundIntersection = intersectVecTriangle(curr_pt, Gradients[ng],
                                                                 edgeNodes1[0],
                                                                 edgeNodes2[0],
                                                                 edgeNodes2[1], pt_i);
                        itmpNode = edgeNodes1[0];
                        if ( !foundIntersection ) {
                            foundIntersection = intersectVecTriangle(curr_pt, Gradients[ng],
                                                                     edgeNodes1[1],
                                                                     edgeNodes2[0],
                                                                     edgeNodes2[1], pt_i);
                            itmpNode = edgeNodes1[1];
                        }
                        if ( !foundIntersection ) {
                            continue;
                        }
                        InLimits=false;
                        prev_pt = curr_pt;
                        curr_pt = pt_i;
                        r_tmp.push_back( curr_pt );
                        bool break_flag = false;
                        for ( size_t n2=0; n2<4; ++n2 ) {
                            if ( nodes[ neighbors[cellNo][n2] ].getDistance( curr_pt ) < small*small ) {
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;
                                edgeNodesPrev = edgeNodes;
                                nodeNo = neighbors[cellNo][n2];
                                onNode = true;
                                onEdge = false;
                                onFace = false;
                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;
                        
                        if ((nodes[itmpNode].getDistance(curr_pt)+nodes[edgeNodes2[0]].getDistance(curr_pt)-nodes[ itmpNode].getDistance(nodes[edgeNodes2[0]])<minDist*minDist)) {
                            edgeNodesPrev = edgeNodes;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            edgeNodes[0] = itmpNode;
                            edgeNodes[1] = edgeNodes2[0];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
   
                            break_flag = true;
                            break;
                        } else if((nodes[itmpNode].getDistance(curr_pt)+nodes[edgeNodes2[1]].getDistance(curr_pt)-nodes[ itmpNode].getDistance(nodes[edgeNodes2[1]])<minDist*minDist)){
                            edgeNodesPrev = edgeNodes;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            edgeNodes[0] = itmpNode;
                            edgeNodes[1] = edgeNodes2[1];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
             
                            break_flag = true;
                            break;
                        } else if ((nodes[edgeNodes2[1]].getDistance(curr_pt)+nodes[edgeNodes2[0]].getDistance(curr_pt)-nodes[edgeNodes2[1]].getDistance(nodes[edgeNodes2[0]])<minDist*minDist)){
                            edgeNodesPrev = edgeNodes;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            edgeNodes[0] = edgeNodes2[0];
                            edgeNodes[1] = edgeNodes2[1];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            break_flag = true;
                            break;
                        }
                        if ( break_flag ) break;
                        
                        onNodePrev = onNode;
                        onEdgePrev = onEdge;
                        onFacePrev = onFace;
                        onNode = false;
                        onEdge = false;
                        onFace = true;
                        
                        //faceNodesPrev = faceNodes;
                        edgeNodesPrev = edgeNodes;
                        faceNodes[0] = itmpNode;
                        faceNodes[1] = edgeNodes2[0];
                        faceNodes[2] = edgeNodes2[1];
                        std::sort(faceNodes.begin(), faceNodes.end());
                        // find next cell
                        T2 PrevCell(cellNo);
                        cellNo=findAdjacentCell2(faceNodes, cellNo,curr_pt);
                        
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            InLimits=true;
                            cellNo=PrevCell;
                        }
                        break;
                    }
                    if (foundIntersection)
                        break;
                }
                if (InLimits && ! foundIntersection){
                    foundIntersection=true;// we project again th gradient
                    
                }
                if (!foundIntersection){
                    InLimits=InLimits? false:true;
                }
                if ( foundIntersection == false && InLimits==false) {
                    std::cout << "\n\nWarning: finding raypath on edge failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    reachedTx = true;
                }
            } else{ // on Face
                bool foundIntersection=false;
                std::set<NODE*> nnodes;
                if (!InLimits){
                    getNeighborNodes(cellNo, nnodes);
                    T1 t;
                    if(r_tmp.size()<=1){
                        t=Interpolator<T1>::TrilinearTime(curr_pt, nodes[neighbors[cellNo][0]], nodes[neighbors[cellNo][1]], nodes[neighbors[cellNo][2]], nodes[neighbors[cellNo][3]], threadNo);
                    }else{
                        t=Interpolator<T1>::bilinearTime(curr_pt, nodes[faceNodes[0]],  nodes[faceNodes[1]],  nodes[faceNodes[2]], threadNo);
                    }
                    g = grad3d.ls_grad(nnodes,threadNo,t,curr_pt);
                    for(size_t nt=0;nt<Tx.size();++nt){
                        for(auto n=neighbors[txCell[nt]].begin();n!=neighbors[txCell[nt]].begin()+4;++n){
                            bool Find=false;
                            for(auto nc=nodes[*n].getOwners().begin();nc!=nodes[*n].getOwners().end();++nc){
                                if(*nc==cellNo){
                                    g=Tx[nt]-curr_pt;
                                    Find=true;
                                    break;
                                }
                            }
                            if (Find)
                                break;
                        }
                        if(curr_pt.getDistance(Tx[nt])<50.0*1.e-3){
                            g=Tx[nt]-curr_pt;
                            break;
                        }

                    }
                }else{
                    g=reflectedGradient(g, faceNodes,cellNo);
                }
                //std::cout<<"nomr g : "<<norm(g)<<std::endl;
                std::array<T2,3> ind[4] = {
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][2] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][2], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][1], neighbors[cellNo][2], neighbors[cellNo][3] } }
                };
                
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                // there are 3 faces that might be intersected
                for ( size_t n=0; n<4; ++n ) {
                    if ( ind[n] == faceNodes) continue;
                    
                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);
                    
                    if ( !foundIntersection )
                        continue;
                    InLimits=false;
                    prev_pt = curr_pt;
                    curr_pt = pt_i;
                    r_tmp.push_back( curr_pt );
                    
                    bool break_flag = false;
                    for ( size_t n2=0; n2<3; ++n2 ) {
                        if ( nodes[ ind[n][n2] ].getDistance( curr_pt ) < small*small ) {
                            //nodeNoPrev = nodeNo;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            nodeNo = ind[n][n2];
                            onNode = true;
                            onEdge = false;
                            onFace = false;
                            faceNodesPrev = faceNodes;

                            break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;
                    
                    for ( size_t n1=0; n1<3; ++n1 ) {
                        size_t n2 = (n1+1)%3;
                        if (nodes[ind[n][n1]].getDistance(curr_pt)+nodes[ind[n][n2]].getDistance(curr_pt)-nodes[ind[n][n1]].getDistance(nodes[ind[n][n2]])<minDist*minDist) {
                            //edgeNodesPrev = edgeNodes;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            edgeNodes[0] = ind[n][n1];
                            edgeNodes[1] = ind[n][n2];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            faceNodesPrev = faceNodes;
 
                            break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;
                    
                    onNodePrev = onNode;
                    onEdgePrev = onEdge;
                    onFacePrev = onFace;
                    onNode = false;
                    onEdge = false;
                    onFace = true;
                    
                    faceNodesPrev = faceNodes;
                    
                    faceNodes = ind[n];
                    // find next cell
                    T2 PrevCell(cellNo);
                    cellNo=findAdjacentCell2(faceNodes, cellNo,curr_pt);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        InLimits=true;
                        cellNo=PrevCell;
                    }
                    break;
                }
                if ( foundIntersection == false ) {
                    
                    // we must be on an face with gradient pointing slightly outward tetrahedron
                    // return in other cell but keep gradient
                    T2 PrevCell(cellNo);
                    cellNo=findAdjacentCell2(faceNodes, cellNo,curr_pt);
                    if (cellNo==std::numeric_limits<T2>::max())
                        cellNo=PrevCell;
                    std::array<T2,3> ind[4];
                    ind[0] = { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][2] } };
                    ind[1] = { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][3] } };
                    ind[2] = { { neighbors[cellNo][0], neighbors[cellNo][2], neighbors[cellNo][3] } };
                    ind[3] = { { neighbors[cellNo][1], neighbors[cellNo][2], neighbors[cellNo][3] } };
                    
                    for ( size_t n=0; n<4; ++n )
                        std::sort( ind[n].begin(), ind[n].end() );
                    
                    for ( size_t n=0; n<4; ++n ) {
                        if ( ind[n] == faceNodes || AreSameFace(ind[n], faceNodes)) continue;
                        
                        sxyz<T1> pt_i;
                        foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                                 ind[n][1], ind[n][2],
                                                                 pt_i);
                        
                        if ( !foundIntersection ) {
                            continue;
                        }
                        prev_pt = curr_pt;
                        InLimits=false;
                        curr_pt = pt_i;
                        r_tmp.push_back( curr_pt );
                        bool break_flag = false;
                        for ( size_t n2=0; n2<3; ++n2 ) {
                            if ( nodes[ ind[n][n2] ].getDistance( curr_pt ) < small*small ) {
                                //nodeNoPrev = nodeNo;
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;
                                
                                nodeNo = ind[n][n2];
                                onNode = true;
                                onEdge = false;
                                onFace = false;
                                faceNodesPrev = faceNodes;
                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;
                        
                        for ( size_t n1=0; n1<3; ++n1 ) {
                            size_t n2 = (n1+1)%3;
                            if (nodes[ind[n][n1]].getDistance(curr_pt)+nodes[ind[n][n2]].getDistance(curr_pt)-nodes[ind[n][n1]].getDistance(nodes[ind[n][n2]])<minDist*minDist) {
                                //edgeNodesPrev = edgeNodes;
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;
                                
                                edgeNodes[0] = ind[n][n1];
                                edgeNodes[1] = ind[n][n2];
                                onNode = false;
                                onEdge = true;
                                onFace = false;
                                faceNodesPrev = faceNodes;
                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;
                        
                        onNodePrev = onNode;
                        onEdgePrev = onEdge;
                        onFacePrev = onFace;
                        onNode = false;
                        onEdge = false;
                        onFace = true;
                        
                        faceNodesPrev = faceNodes;
                        faceNodes = ind[n];
                        // find next cell
                        T2 PrevCell(cellNo);
                        cellNo=findAdjacentCell2(faceNodes, cellNo,curr_pt);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            InLimits=true;
                            cellNo=PrevCell;
                            
                        }
                        break;
                    }
                }
                if (!foundIntersection)
                    InLimits=InLimits? false:true;
                if ( foundIntersection == false && InLimits==false ) {
                    std::cout << "\n\nWarning: finding raypath on face failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    reachedTx = true;
                }
            }
            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( curr_pt.getDistance( Tx[nt] ) < minDist*minDist ) {
                        reachedTx = true;
                        break;
                    }
                    for(size_t n=0;n<4;++n){
                        sxyz<T1> NearTx={nodes[neighbors[txCell[nt]][n]]};
                        if (curr_pt.getDistance(NearTx)< minDist*minDist){
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
    void Grid3Dui<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
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

                // find adjacent cells
                T2 ind[6][2] = {
                    {neighbors[txCell[nt]][0], neighbors[txCell[nt]][1]},
                    {neighbors[txCell[nt]][0], neighbors[txCell[nt]][2]},
                    {neighbors[txCell[nt]][0], neighbors[txCell[nt]][3]},
                    {neighbors[txCell[nt]][1], neighbors[txCell[nt]][2]},
                    {neighbors[txCell[nt]][1], neighbors[txCell[nt]][3]},
                    {neighbors[txCell[nt]][2], neighbors[txCell[nt]][3]} };

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
        Grad3D<T1,NODE> grad3d;
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

            T2 ind[6][2] = {
                {neighbors[cellNo][0], neighbors[cellNo][1]},
                {neighbors[cellNo][0], neighbors[cellNo][2]},
                {neighbors[cellNo][0], neighbors[cellNo][3]},
                {neighbors[cellNo][1], neighbors[cellNo][2]},
                {neighbors[cellNo][1], neighbors[cellNo][3]},
                {neighbors[cellNo][2], neighbors[cellNo][3]} };

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
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][2] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][2], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][1], neighbors[cellNo][2], neighbors[cellNo][3] } }
                };
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );

                for ( size_t n=0; n<4; ++n ) {
                    if ( areCoplanar(curr_pt, ind[n][0], ind[n][1], ind[n][2]) ) {
                        onFace = true;
                        faceNodes[0] = ind[n][0];
                        faceNodes[1] = ind[n][1];
                        faceNodes[2] = ind[n][2];
                        break;
                    }
                }
            }
        }
		T1 s, ds;
        while ( reachedTx == false ) {

            if ( onNode ) {

                // find cells common to edge
                std::vector<T2> cells;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                    cells.push_back( *nc );
                }

                // compute gradient with nodes from all common cells
                std::set<NODE*> nnodes;
                for (size_t n=0; n<cells.size(); ++n ) {
                    for ( size_t no=0; no<4; ++no ) {
                        nnodes.insert( &(nodes[ neighbors[cells[n]][no] ]) );
                    }
                }
                sxyz<T1> g = grad3d.ls_grad(nnodes, threadNo,curr_pt);

                // find cell for which gradient intersect opposing face
                bool foundIntersection = false;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                    std::array<T2,3> nb;
                    size_t n=0;
                    for (auto nn=neighbors[*nc].begin(); nn!=neighbors[*nc].end(); ++nn ) {
                        if ( *nn != nodeNo ) {
                            nb[n++] = *nn;
                        }
                    }
                    std::sort(nb.begin(), nb.end());

                    foundIntersection = intersectVecTriangle( nodeNo, g, nb[0], nb[1], nb[2], curr_pt);
                    if ( !foundIntersection ) {
                        continue;
                    }

                    prev_pt = curr_pt;
                    r_tmp.push_back( curr_pt );

                    if ( r_tmp.size() > 1 ) {
                        // compute terms of matrix M
                        mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                        s = computeSlowness(mid_pt);
                        s *= s;
                        ds = curr_pt.getDistance( prev_pt );
					}

                    bool break_flag = false;
                    for ( n=0; n<3; ++n ) {
                        if ( nodes[ nb[n] ].getDistance( curr_pt ) < small*small ) {
							nodeNoPrev = nodeNo;
							onNodePrev = onNode;
							onEdgePrev = onEdge;
							onFacePrev = onFace;

                            nodeNo = nb[n];
                            onNode = true;
                            onEdge = false;
                            onFace = false;

							if ( r_tmp.size() > 1) {
								std::set<T2> allNodes;
								if (onNodePrev) allNodes.insert(nodeNoPrev);
								if (onNode) allNodes.insert(nodeNo);
								if (onEdgePrev) {
									allNodes.insert( edgeNodesPrev[0] );
									allNodes.insert( edgeNodesPrev[1] );
								}
								if (onEdge) {
									allNodes.insert( edgeNodes[0] );
									allNodes.insert( edgeNodes[1] );
								}
								if (onFacePrev) {
									allNodes.insert( faceNodesPrev[0] );
									allNodes.insert( faceNodesPrev[1] );
									allNodes.insert( faceNodesPrev[2] );
								}
								if (onFace) {
									allNodes.insert( faceNodes[0] );
									allNodes.insert( faceNodes[1] );
									allNodes.insert( faceNodes[2] );
								}

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

							break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;

                    for ( size_t n1=0; n1<3; ++n1 ) {   // changed n1<2  -> n1<3
                        size_t n2 = (n1+1)%3;
                        if ( areCollinear(curr_pt, nb[n1], nb[n2]) ) {
							edgeNodesPrev = edgeNodes;
							onNodePrev = onNode;
							onEdgePrev = onEdge;
							onFacePrev = onFace;

							edgeNodes[0] = nb[n1];
                            edgeNodes[1] = nb[n2];
                            onNode = false;
                            onEdge = true;
                            onFace = false;

							if ( r_tmp.size() > 1) {
								std::set<T2> allNodes;
								if (onNodePrev) allNodes.insert(nodeNoPrev);
								if (onNode) allNodes.insert(nodeNo);
								if (onEdgePrev) {
									allNodes.insert( edgeNodesPrev[0] );
									allNodes.insert( edgeNodesPrev[1] );
								}
								if (onEdge) {
									allNodes.insert( edgeNodes[0] );
									allNodes.insert( edgeNodes[1] );
								}
								if (onFacePrev) {
									allNodes.insert( faceNodesPrev[0] );
									allNodes.insert( faceNodesPrev[1] );
									allNodes.insert( faceNodesPrev[2] );
								}
								if (onFace) {
									allNodes.insert( faceNodes[0] );
									allNodes.insert( faceNodes[1] );
									allNodes.insert( faceNodes[2] );
								}

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

                            break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;

					onNodePrev = onNode;
					onEdgePrev = onEdge;
					onFacePrev = onFace;
					onNode = false;
                    onEdge = false;
                    onFace = true;

					faceNodesPrev = faceNodes;
                    faceNodes = nb;

					if ( r_tmp.size() > 1) {
						std::set<T2> allNodes;
						if (onNodePrev) allNodes.insert(nodeNoPrev);
						if (onNode) allNodes.insert(nodeNo);
						if (onEdgePrev) {
							allNodes.insert( edgeNodesPrev[0] );
							allNodes.insert( edgeNodesPrev[1] );
						}
						if (onEdge) {
							allNodes.insert( edgeNodes[0] );
							allNodes.insert( edgeNodes[1] );
						}
						if (onFacePrev) {
							allNodes.insert( faceNodesPrev[0] );
							allNodes.insert( faceNodesPrev[1] );
							allNodes.insert( faceNodesPrev[2] );
						}
						if (onFace) {
							allNodes.insert( faceNodes[0] );
							allNodes.insert( faceNodes[1] );
							allNodes.insert( faceNodes[2] );
						}

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
                    if ( std::find(nodes[edgeNodes[1]].getOwners().begin(), nodes[edgeNodes[1]].getOwners().end(), *nc0)!=nodes[edgeNodes[1]].getOwners().end() )
                        cells.push_back( *nc0 );
                }
                // compute gradient with nodes from all common cells
                std::set<NODE*> nnodes;
                for (size_t n=0; n<cells.size(); ++n ) {
                    for ( size_t no=0; no<4; ++no ) {
                        nnodes.insert( &(nodes[ neighbors[cells[n]][no] ]) );
                    }
                }
                sxyz<T1> g = grad3d.ls_grad(nnodes, threadNo,curr_pt);

                bool foundIntersection=false;
                for (size_t n=0; n<cells.size(); ++n ) {

                    cellNo = cells[n];

                    // there are 2 faces that might be intersected
                    std::array<T2,2> edgeNodes2;
                    size_t n2=0;
                    for ( auto nn=neighbors[cellNo].begin(); nn!= neighbors[cellNo].end(); ++nn ) {
                        if ( *nn!=edgeNodes[0] && *nn!=edgeNodes[1] ) {
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
                        if ( nodes[ neighbors[cellNo][n2] ].getDistance( curr_pt ) < small*small ) {
							onNodePrev = onNode;
							onEdgePrev = onEdge;
							onFacePrev = onFace;
							nodeNo = neighbors[cellNo][n2];

                            onNode = true;
                            onEdge = false;
                            onFace = false;
                            break_flag = true;

							if ( r_tmp.size() > 1) {
								std::set<T2> allNodes;
								if (onNodePrev) allNodes.insert(nodeNoPrev);
								if (onNode) allNodes.insert(nodeNo);
								if (onEdgePrev) {
									allNodes.insert( edgeNodesPrev[0] );
									allNodes.insert( edgeNodesPrev[1] );
								}
								if (onEdge) {
									allNodes.insert( edgeNodes[0] );
									allNodes.insert( edgeNodes[1] );
								}
								if (onFacePrev) {
									allNodes.insert( faceNodesPrev[0] );
									allNodes.insert( faceNodesPrev[1] );
									allNodes.insert( faceNodesPrev[2] );
								}
								if (onFace) {
									allNodes.insert( faceNodes[0] );
									allNodes.insert( faceNodes[1] );
									allNodes.insert( faceNodes[2] );
								}

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

							break;
                        }
                    }
                    if ( break_flag ) break;

                    if ( areCollinear(curr_pt, itmpNode, edgeNodes2[0]) ) {
						edgeNodesPrev = edgeNodes;
						onNodePrev = onNode;
						onEdgePrev = onEdge;
						onFacePrev = onFace;

						edgeNodes[0] = itmpNode;
                        edgeNodes[1] = edgeNodes2[0];
                        onNode = false;
                        onEdge = true;
                        onFace = false;

						if ( r_tmp.size() > 1) {
							std::set<T2> allNodes;
							if (onNodePrev) allNodes.insert(nodeNoPrev);
							if (onNode) allNodes.insert(nodeNo);
							if (onEdgePrev) {
								allNodes.insert( edgeNodesPrev[0] );
								allNodes.insert( edgeNodesPrev[1] );
							}
							if (onEdge) {
								allNodes.insert( edgeNodes[0] );
								allNodes.insert( edgeNodes[1] );
							}
							if (onFacePrev) {
								allNodes.insert( faceNodesPrev[0] );
								allNodes.insert( faceNodesPrev[1] );
								allNodes.insert( faceNodesPrev[2] );
							}
							if (onFace) {
								allNodes.insert( faceNodes[0] );
								allNodes.insert( faceNodes[1] );
								allNodes.insert( faceNodes[2] );
							}

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

						break_flag = true;
                        break;
                    } else if ( areCollinear(curr_pt, itmpNode, edgeNodes2[1]) ) {
						edgeNodesPrev = edgeNodes;
						onNodePrev = onNode;
						onEdgePrev = onEdge;
						onFacePrev = onFace;

						edgeNodes[0] = itmpNode;
                        edgeNodes[1] = edgeNodes2[1];
                        onNode = false;
                        onEdge = true;
                        onFace = false;

						if ( r_tmp.size() > 1) {
							std::set<T2> allNodes;
							if (onNodePrev) allNodes.insert(nodeNoPrev);
							if (onNode) allNodes.insert(nodeNo);
							if (onEdgePrev) {
								allNodes.insert( edgeNodesPrev[0] );
								allNodes.insert( edgeNodesPrev[1] );
							}
							if (onEdge) {
								allNodes.insert( edgeNodes[0] );
								allNodes.insert( edgeNodes[1] );
							}
							if (onFacePrev) {
								allNodes.insert( faceNodesPrev[0] );
								allNodes.insert( faceNodesPrev[1] );
								allNodes.insert( faceNodesPrev[2] );
							}
							if (onFace) {
								allNodes.insert( faceNodes[0] );
								allNodes.insert( faceNodes[1] );
								allNodes.insert( faceNodes[2] );
							}

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

						break_flag = true;
                        break;
                    } else if ( areCollinear(curr_pt, edgeNodes2[0], edgeNodes2[1]) ) {
						edgeNodesPrev = edgeNodes;
						onNodePrev = onNode;
						onEdgePrev = onEdge;
						onFacePrev = onFace;

						edgeNodes[0] = edgeNodes2[0];
                        edgeNodes[1] = edgeNodes2[1];
                        onNode = false;
                        onEdge = true;
                        onFace = false;

						if ( r_tmp.size() > 1) {
							std::set<T2> allNodes;
							if (onNodePrev) allNodes.insert(nodeNoPrev);
							if (onNode) allNodes.insert(nodeNo);
							if (onEdgePrev) {
								allNodes.insert( edgeNodesPrev[0] );
								allNodes.insert( edgeNodesPrev[1] );
							}
							if (onEdge) {
								allNodes.insert( edgeNodes[0] );
								allNodes.insert( edgeNodes[1] );
							}
							if (onFacePrev) {
								allNodes.insert( faceNodesPrev[0] );
								allNodes.insert( faceNodesPrev[1] );
								allNodes.insert( faceNodesPrev[2] );
							}
							if (onFace) {
								allNodes.insert( faceNodes[0] );
								allNodes.insert( faceNodes[1] );
								allNodes.insert( faceNodes[2] );
							}

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

						break_flag = true;
                        break;
                    }
                    if ( break_flag ) break;

					onNodePrev = onNode;
					onEdgePrev = onEdge;
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
						if (onNodePrev) allNodes.insert(nodeNoPrev);
						if (onNode) allNodes.insert(nodeNo);
						if (onEdgePrev) {
							allNodes.insert( edgeNodesPrev[0] );
							allNodes.insert( edgeNodesPrev[1] );
						}
						if (onEdge) {
							allNodes.insert( edgeNodes[0] );
							allNodes.insert( edgeNodes[1] );
						}
						if (onFacePrev) {
							allNodes.insert( faceNodesPrev[0] );
							allNodes.insert( faceNodesPrev[1] );
							allNodes.insert( faceNodesPrev[2] );
						}
						if (onFace) {
							allNodes.insert( faceNodes[0] );
							allNodes.insert( faceNodes[1] );
							allNodes.insert( faceNodes[2] );
						}

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

                assert(neighbors[cellNo].size()==4);
                sxyz<T1> g = grad3d.ls_grad(nodes[ neighbors[cellNo][0] ],
                                            nodes[ neighbors[cellNo][1] ],
                                            nodes[ neighbors[cellNo][2] ],
                                            nodes[ neighbors[cellNo][3] ], threadNo);

                std::array<T2,3> ind[4] = {
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][2] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][2], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][1], neighbors[cellNo][2], neighbors[cellNo][3] } }
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
                        if ( nodes[ ind[n][n2] ].getDistance( curr_pt ) < small*small ) {
							nodeNoPrev = nodeNo;
							onNodePrev = onNode;
							onEdgePrev = onEdge;
							onFacePrev = onFace;

							nodeNo = ind[n][n2];
                            onNode = true;
                            onEdge = false;
                            onFace = false;

							if ( r_tmp.size() > 1) {
								std::set<T2> allNodes;
								if (onNodePrev) allNodes.insert(nodeNoPrev);
								if (onNode) allNodes.insert(nodeNo);
								if (onEdgePrev) {
									allNodes.insert( edgeNodesPrev[0] );
									allNodes.insert( edgeNodesPrev[1] );
								}
								if (onEdge) {
									allNodes.insert( edgeNodes[0] );
									allNodes.insert( edgeNodes[1] );
								}
								if (onFacePrev) {
									allNodes.insert( faceNodesPrev[0] );
									allNodes.insert( faceNodesPrev[1] );
									allNodes.insert( faceNodesPrev[2] );
								}
								if (onFace) {
									allNodes.insert( faceNodes[0] );
									allNodes.insert( faceNodes[1] );
									allNodes.insert( faceNodes[2] );
								}

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
							onFacePrev = onFace;

							edgeNodes[0] = ind[n][n1];
                            edgeNodes[1] = ind[n][n2];
                            onNode = false;
                            onEdge = true;
                            onFace = false;

							if ( r_tmp.size() > 1) {
								std::set<T2> allNodes;
								if (onNodePrev) allNodes.insert(nodeNoPrev);
								if (onNode) allNodes.insert(nodeNo);
								if (onEdgePrev) {
									allNodes.insert( edgeNodesPrev[0] );
									allNodes.insert( edgeNodesPrev[1] );
								}
								if (onEdge) {
									allNodes.insert( edgeNodes[0] );
									allNodes.insert( edgeNodes[1] );
								}
								if (onFacePrev) {
									allNodes.insert( faceNodesPrev[0] );
									allNodes.insert( faceNodesPrev[1] );
									allNodes.insert( faceNodesPrev[2] );
								}
								if (onFace) {
									allNodes.insert( faceNodes[0] );
									allNodes.insert( faceNodes[1] );
									allNodes.insert( faceNodes[2] );
								}

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

							break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;

					onNodePrev = onNode;
					onEdgePrev = onEdge;
					onFacePrev = onFace;
					onNode = false;
                    onEdge = false;
                    onFace = true;

					faceNodesPrev = faceNodes;
                    faceNodes = ind[n];

					if ( r_tmp.size() > 1) {
						std::set<T2> allNodes;
						if (onNodePrev) allNodes.insert(nodeNoPrev);
						if (onNode) allNodes.insert(nodeNo);
						if (onEdgePrev) {
							allNodes.insert( edgeNodesPrev[0] );
							allNodes.insert( edgeNodesPrev[1] );
						}
						if (onEdge) {
							allNodes.insert( edgeNodes[0] );
							allNodes.insert( edgeNodes[1] );
						}
						if (onFacePrev) {
							allNodes.insert( faceNodesPrev[0] );
							allNodes.insert( faceNodesPrev[1] );
							allNodes.insert( faceNodesPrev[2] );
						}
						if (onFace) {
							allNodes.insert( faceNodes[0] );
							allNodes.insert( faceNodes[1] );
							allNodes.insert( faceNodes[2] );
						}

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

                    ind[0] = { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][2] } };
                    ind[1] = { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][3] } };
                    ind[2] = { { neighbors[cellNo][0], neighbors[cellNo][2], neighbors[cellNo][3] } };
                    ind[3] = { { neighbors[cellNo][1], neighbors[cellNo][2], neighbors[cellNo][3] } };

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
                            if ( nodes[ ind[n][n2] ].getDistance( curr_pt ) < small*small ) {
								nodeNoPrev = nodeNo;
								onNodePrev = onNode;
								onEdgePrev = onEdge;
								onFacePrev = onFace;

								nodeNo = ind[n][n2];
                                onNode = true;
                                onEdge = false;
                                onFace = false;

								if ( r_tmp.size() > 1) {
									std::set<T2> allNodes;
									if (onNodePrev) allNodes.insert(nodeNoPrev);
									if (onNode) allNodes.insert(nodeNo);
									if (onEdgePrev) {
										allNodes.insert( edgeNodesPrev[0] );
										allNodes.insert( edgeNodesPrev[1] );
									}
									if (onEdge) {
										allNodes.insert( edgeNodes[0] );
										allNodes.insert( edgeNodes[1] );
									}
									if (onFacePrev) {
										allNodes.insert( faceNodesPrev[0] );
										allNodes.insert( faceNodesPrev[1] );
										allNodes.insert( faceNodesPrev[2] );
									}
									if (onFace) {
										allNodes.insert( faceNodes[0] );
										allNodes.insert( faceNodes[1] );
										allNodes.insert( faceNodes[2] );
									}

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
								onFacePrev = onFace;

								edgeNodes[0] = ind[n][n1];
                                edgeNodes[1] = ind[n][n2];
                                onNode = false;
                                onEdge = true;
                                onFace = false;

								if ( r_tmp.size() > 1) {
									std::set<T2> allNodes;
									if (onNodePrev) allNodes.insert(nodeNoPrev);
									if (onNode) allNodes.insert(nodeNo);
									if (onEdgePrev) {
										allNodes.insert( edgeNodesPrev[0] );
										allNodes.insert( edgeNodesPrev[1] );
									}
									if (onEdge) {
										allNodes.insert( edgeNodes[0] );
										allNodes.insert( edgeNodes[1] );
									}
									if (onFacePrev) {
										allNodes.insert( faceNodesPrev[0] );
										allNodes.insert( faceNodesPrev[1] );
										allNodes.insert( faceNodesPrev[2] );
									}
									if (onFace) {
										allNodes.insert( faceNodes[0] );
										allNodes.insert( faceNodes[1] );
										allNodes.insert( faceNodes[2] );
									}

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

								break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;

						onNodePrev = onNode;
						onEdgePrev = onEdge;
						onFacePrev = onFace;
						onNode = false;
                        onEdge = false;
                        onFace = true;

						faceNodesPrev = faceNodes;
                        faceNodes = ind[n];

						if ( r_tmp.size() > 1) {
							std::set<T2> allNodes;
							if (onNodePrev) allNodes.insert(nodeNoPrev);
							if (onNode) allNodes.insert(nodeNo);
							if (onEdgePrev) {
								allNodes.insert( edgeNodesPrev[0] );
								allNodes.insert( edgeNodesPrev[1] );
							}
							if (onEdge) {
								allNodes.insert( edgeNodes[0] );
								allNodes.insert( edgeNodes[1] );
							}
							if (onFacePrev) {
								allNodes.insert( faceNodesPrev[0] );
								allNodes.insert( faceNodesPrev[1] );
								allNodes.insert( faceNodesPrev[2] );
							}
							if (onFace) {
								allNodes.insert( faceNodes[0] );
								allNodes.insert( faceNodes[1] );
								allNodes.insert( faceNodes[2] );
							}

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
    void Grid3Dui<T1,T2,NODE>::getRaypath_ho(const std::vector<sxyz<T1>>& Tx,
                                             const sxyz<T1>& Rx,
                                             std::vector<sxyz<T1>>& r_data,
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

                // find adjacent cells
                T2 ind[6][2] = {
                    {neighbors[txCell[nt]][0], neighbors[txCell[nt]][1]},
                    {neighbors[txCell[nt]][0], neighbors[txCell[nt]][2]},
                    {neighbors[txCell[nt]][0], neighbors[txCell[nt]][3]},
                    {neighbors[txCell[nt]][1], neighbors[txCell[nt]][2]},
                    {neighbors[txCell[nt]][1], neighbors[txCell[nt]][3]},
                    {neighbors[txCell[nt]][2], neighbors[txCell[nt]][3]} };

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
        
        bool InLimits=false;
        bool onNode = false;
        bool onEdge = false;
        bool onFace = false;
        bool onNodePrev = false;
        bool onEdgePrev = false;
        bool onFacePrev = false;
        bool secondNodes=nodes.size()>nPrimary;
        std::array<T2,2> edgeNodes, edgeNodesPrev;
        std::array<T2,3> faceNodes={{0,0,0}};
        std::array<T2,3> faceNodesPrev;
        Grad3D_ho <T1,NODE> grad3d;
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

            T2 ind[6][2] = {
                {neighbors[cellNo][0], neighbors[cellNo][1]},
                {neighbors[cellNo][0], neighbors[cellNo][2]},
                {neighbors[cellNo][0], neighbors[cellNo][3]},
                {neighbors[cellNo][1], neighbors[cellNo][2]},
                {neighbors[cellNo][1], neighbors[cellNo][3]},
                {neighbors[cellNo][2], neighbors[cellNo][3]} };

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
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][2] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][2], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][1], neighbors[cellNo][2], neighbors[cellNo][3]} } };
                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );

                for ( size_t n=0; n<4; ++n ) {
                    if ( areCoplanar(curr_pt, ind[n][0], ind[n][1], ind[n][2]) ) {
                        onFace = true;
                        faceNodes[0] = ind[n][0];
                        faceNodes[1] = ind[n][1];
                        faceNodes[2] = ind[n][2];
                        break;
                    }
                }
            }
        }
        if (!onEdge && ! onNode && ! onFace){
            onFace=true;
        }
        for(auto t=0;t<txCell.size();++t){
            if (getCellNo( Rx )==txCell[t]){
                r_tmp.emplace_back(Tx[t]);
                reachedTx=true;
                break;
            }
        }
        T1 s, ds;
        sxyz<T1> g;
        T2 N=0;
        while ( reachedTx == false && N<500) {
            ++N;
            if ( onNode ) {
                bool foundIntersection = false;
                if (!InLimits){// if the current point is not on the limitis of mesh
                    // find cells common to edge
                    std::set<NODE*> nnodes;
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                        getNeighborNodes(*nc, nnodes);
                    }
                    // compute gradient with nodes from all common cells
                    sxyz<T1> g = grad3d.ls_grad(nnodes, threadNo,curr_pt);
                    
                    // find cell for which gradient intersect opposing face
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                        
                        std::array<T2,3> nb;
                        size_t n=0;
                        for (auto nn=neighbors[*nc].begin(); nn!=neighbors[*nc].begin()+4; ++nn ) {
                            if ( *nn != nodeNo ) {
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
                        InLimits=false;
                        r_tmp.push_back( curr_pt );
                        ds = curr_pt.getDistance( prev_pt );
                        if ( r_tmp.size() > 1 && ds>minDist) {
                            // compute terms of matrix M
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            s=computeSlowness(mid_pt,cellNo );
                            s *= s;
                        }
                        
                        bool break_flag = false;
                        for ( n=0; n<3; ++n ) {
                            if ( nodes[ nb[n] ].getDistance( curr_pt ) < small*small ) {
                                nodeNoPrev = nodeNo;
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;
                                
                                nodeNo = nb[n];
                                onNode = true;
                                onEdge = false;
                                onFace = false;
                                
                                if ( r_tmp.size() > 1 && ds> 0.01) {
                                    std::set<T2> allNodes;
                                    if (onNodePrev) allNodes.insert(nodeNoPrev);
                                    if (onNode) allNodes.insert(nodeNo);
                                    if (onEdgePrev) {
                                        allNodes.insert( edgeNodesPrev[0] );
                                        allNodes.insert( edgeNodesPrev[1] );
                                    }
                                    if (onEdge) {
                                        allNodes.insert( edgeNodes[0] );
                                        allNodes.insert( edgeNodes[1] );
                                    }
                                    if (onFacePrev) {
                                        allNodes.insert( faceNodesPrev[0] );
                                        allNodes.insert( faceNodesPrev[1] );
                                        allNodes.insert( faceNodesPrev[2] );
                                    }
                                    if (onFace) {
                                        allNodes.insert( faceNodes[0] );
                                        allNodes.insert( faceNodes[1] );
                                        allNodes.insert( faceNodes[2] );
                                    }
                                    if (secondNodes){
                                        auto no=--allNodes.end();
                                        while (no!=allNodes.begin()) {
                                            if(nodes[*no].getprimary()==10){
                                                std::copy(nodes[*no].getPrincipals().begin(),nodes[*no].getPrincipals().end(), std::inserter(allNodes,allNodes.end()));
                                                allNodes.erase(no);
                                            }
                                            --no;
                                        }
                                    }
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
                                
                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;
                        
                        for ( size_t n1=0; n1<3; ++n1 ) {
                            size_t n2 = (n1+1)%3;
                            if (nodes[nb[n1]].getDistance(curr_pt)+nodes[nb[n2]].getDistance(curr_pt)-nodes[nb[n2]].getDistance(nodes[nb[n1]])<minDist*minDist) {
                                nodeNoPrev = nodeNo;
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;
                                
                                edgeNodes[0] = nb[n1];
                                edgeNodes[1] = nb[n2];
                                onNode = false;
                                onEdge = true;
                                onFace = false;
                                
                                if ( r_tmp.size() > 1 && ds>minDist) {
                                    std::set<T2> allNodes;
                                    if (onNodePrev) allNodes.insert(nodeNoPrev);
                                    if (onNode) allNodes.insert(nodeNo);
                                    if (onEdgePrev) {
                                        allNodes.insert( edgeNodesPrev[0] );
                                        allNodes.insert( edgeNodesPrev[1] );
                                    }
                                    if (onEdge) {
                                        allNodes.insert( edgeNodes[0] );
                                        allNodes.insert( edgeNodes[1] );
                                    }
                                    if (onFacePrev) {
                                        allNodes.insert( faceNodesPrev[0] );
                                        allNodes.insert( faceNodesPrev[1] );
                                        allNodes.insert( faceNodesPrev[2] );
                                    }
                                    if (onFace) {
                                        allNodes.insert( faceNodes[0] );
                                        allNodes.insert( faceNodes[1] );
                                        allNodes.insert( faceNodes[2] );
                                    }
                                    if (secondNodes){
                                        auto no=--allNodes.end();
                                        while (no!=allNodes.begin()) {
                                            if(nodes[*no].getprimary()==10){
                                                std::copy(nodes[*no].getPrincipals().begin(),nodes[*no].getPrincipals().end(), std::inserter(allNodes,allNodes.end()));
                                                allNodes.erase(no);
                                            }
                                            --no;
                                        }
                                    }
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
                                
                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;
                        
                        onNodePrev = onNode;
                        onEdgePrev = onEdge;
                        onFacePrev = onFace;
                        onNode = false;
                        onEdge = false;
                        onFace = true;
                        
                        nodeNoPrev = nodeNo;
                        faceNodes = nb;
                        if ( r_tmp.size() > 1 && ds>minDist) {
                            std::set<T2> allNodes;
                            if (onNodePrev) allNodes.insert(nodeNoPrev);
                            if (onNode) allNodes.insert(nodeNo);
                            if (onEdgePrev) {
                                allNodes.insert( edgeNodesPrev[0] );
                                allNodes.insert( edgeNodesPrev[1] );
                            }
                            if (onEdge) {
                                allNodes.insert( edgeNodes[0] );
                                allNodes.insert( edgeNodes[1] );
                            }
                            if (onFacePrev) {
                                allNodes.insert( faceNodesPrev[0] );
                                allNodes.insert( faceNodesPrev[1] );
                                allNodes.insert( faceNodesPrev[2] );
                            }
                            if (onFace) {
                                allNodes.insert( faceNodes[0] );
                                allNodes.insert( faceNodes[1] );
                                allNodes.insert( faceNodes[2] );
                            }
                            if (secondNodes){
                                auto no=--allNodes.end();
                                while (no!=allNodes.begin()) {
                                    if(nodes[*no].getprimary()==10){
                                        std::copy(nodes[*no].getPrincipals().begin(),nodes[*no].getPrincipals().end(), std::inserter(allNodes,allNodes.end()));
                                        allNodes.erase(no);
                                    }
                                    --no;
                                }
                            }
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
                        
                        // find next cell
                        T2 PrevCell=cellNo;
                        cellNo = findAdjacentCell1(faceNodes, nodeNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {// cuurent point is in the limit of mesh
                            InLimits=true;
                            cellNo=PrevCell;
                        }
                        break;
                    }
                    
                }else{// cuurent point is in the limit of mesh
                    // we project gradient onto all the bordering faces
                    std::priority_queue<sxyz<T1>,std::vector<sxyz<T1>>,Comparesxyz_norm<T1>> ProjectedGradients;
                    //plotCell(cellNo, curr_pt, g);
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ){
                        std::array<T2,3> InBorderFace;
                        for(T2 nn=0;nn<4;++nn){
                            InBorderFace[0]=neighbors[*nc][(nn+1)%4];
                            InBorderFace[1]=neighbors[*nc][(nn+2)%4];
                            InBorderFace[2]=neighbors[*nc][(nn+3)%4];
                            if(findAdjacentCell2(InBorderFace, *nc, curr_pt)==std::numeric_limits <T2>::max())
                                ProjectedGradients.push((reflectedGradient(g,InBorderFace,*nc)));
                        }
                    }
                    // take the gradient with the biggest norm
                    g=ProjectedGradients.top();
                    while (! ProjectedGradients.empty()){
                        sxyz<T1> g_reflected=ProjectedGradients.top();
                        ProjectedGradients.pop();
                        // find cell for which gradient intersect opposing face
                        for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                            //plotCell(*nc, curr_pt, g_reflected);
                            std::array<T2,3> nb;
                            size_t n=0;
                            for (auto nn=neighbors[*nc].begin(); nn!=neighbors[*nc].begin()+4; ++nn ) {
                                if ( *nn != nodeNo ) {
                                    nb[n++] = *nn;
                                }
                            }
                            std::sort(nb.begin(), nb.end());
                            
                            sxyz<T1> pt_i;
                            foundIntersection = intersectVecTriangle( nodeNo,g_reflected, nb[0], nb[1], nb[2], pt_i);
                            if ( !foundIntersection ) {
                                continue;
                            }
                            
                            prev_pt = curr_pt;
                            curr_pt = pt_i;
                            InLimits=false;
                            r_tmp.push_back( curr_pt );
                            ds = curr_pt.getDistance( prev_pt );
                            if ( r_tmp.size() > 1 && ds>minDist) {
                                // compute terms of matrix M
                                mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                                s=computeSlowness(mid_pt,cellNo );
                                s *= s;
                            }
                            
                            bool break_flag = false;
                            for ( n=0; n<3; ++n ) {
                                if ( nodes[ nb[n] ].getDistance( curr_pt ) < small*small ) {
                                    nodeNoPrev = nodeNo;
                                    onNodePrev = onNode;
                                    onEdgePrev = onEdge;
                                    onFacePrev = onFace;
                                    
                                    nodeNo = nb[n];
                                    onNode = true;
                                    onEdge = false;
                                    onFace = false;
                                    
                                    if ( r_tmp.size() > 1 && ds>minDist) {
                                        std::set<T2> allNodes;
                                        if (onNodePrev) allNodes.insert(nodeNoPrev);
                                        if (onNode) allNodes.insert(nodeNo);
                                        if (onEdgePrev) {
                                            allNodes.insert( edgeNodesPrev[0] );
                                            allNodes.insert( edgeNodesPrev[1] );
                                        }
                                        if (onEdge) {
                                            allNodes.insert( edgeNodes[0] );
                                            allNodes.insert( edgeNodes[1] );
                                        }
                                        if (onFacePrev) {
                                            allNodes.insert( faceNodesPrev[0] );
                                            allNodes.insert( faceNodesPrev[1] );
                                            allNodes.insert( faceNodesPrev[2] );
                                        }
                                        if (onFace) {
                                            allNodes.insert( faceNodes[0] );
                                            allNodes.insert( faceNodes[1] );
                                            allNodes.insert( faceNodes[2] );
                                        }
                                        if (secondNodes){
                                            auto no=--allNodes.end();
                                            while (no!=allNodes.begin()) {
                                                if(nodes[*no].getprimary()==10){
                                                    std::copy(nodes[*no].getPrincipals().begin(),nodes[*no].getPrincipals().end(), std::inserter(allNodes,allNodes.end()));
                                                    allNodes.erase(no);
                                                }
                                                --no;
                                            }
                                        }
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
                                    
                                    break_flag = true;
                                    break;
                                }
                            }
                            if ( break_flag ) break;
                            
                            for ( size_t n1=0; n1<3; ++n1 ) {
                                size_t n2 = (n1+1)%3;
                                if (nodes[nb[n1]].getDistance(curr_pt)+nodes[nb[n2]].getDistance(curr_pt)-nodes[nb[n2]].getDistance(nodes[nb[n1]])<minDist*minDist) {
                                    nodeNoPrev = nodeNo;
                                    onNodePrev = onNode;
                                    onEdgePrev = onEdge;
                                    onFacePrev = onFace;
                                    
                                    edgeNodes[0] = nb[n1];
                                    edgeNodes[1] = nb[n2];
                                    onNode = false;
                                    onEdge = true;
                                    onFace = false;
                                    
                                    if ( r_tmp.size() > 1 && ds>minDist) {
                                        std::set<T2> allNodes;
                                        if (onNodePrev) allNodes.insert(nodeNoPrev);
                                        if (onNode) allNodes.insert(nodeNo);
                                        if (onEdgePrev) {
                                            allNodes.insert( edgeNodesPrev[0] );
                                            allNodes.insert( edgeNodesPrev[1] );
                                        }
                                        if (onEdge) {
                                            allNodes.insert( edgeNodes[0] );
                                            allNodes.insert( edgeNodes[1] );
                                        }
                                        if (onFacePrev) {
                                            allNodes.insert( faceNodesPrev[0] );
                                            allNodes.insert( faceNodesPrev[1] );
                                            allNodes.insert( faceNodesPrev[2] );
                                        }
                                        if (onFace) {
                                            allNodes.insert( faceNodes[0] );
                                            allNodes.insert( faceNodes[1] );
                                            allNodes.insert( faceNodes[2] );
                                        }
                                        if (secondNodes){
                                            auto no=--allNodes.end();
                                            while (no!=allNodes.begin()) {
                                                if(nodes[*no].getprimary()==10){
                                                    std::copy(nodes[*no].getPrincipals().begin(),nodes[*no].getPrincipals().end(), std::inserter(allNodes,allNodes.end()));
                                                    allNodes.erase(no);
                                                }
                                                --no;
                                            }
                                        }
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
                                    
                                    break_flag = true;
                                    break;
                                }
                            }
                            if ( break_flag ) break;
                            
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            onNode = false;
                            onEdge = false;
                            onFace = true;
                            faceNodes = nb;
                            if ( r_tmp.size() > 1 && ds>minDist) {
                                std::set<T2> allNodes;
                                if (onNodePrev) allNodes.insert(nodeNoPrev);
                                if (onNode) allNodes.insert(nodeNo);
                                if (onEdgePrev) {
                                    allNodes.insert( edgeNodesPrev[0] );
                                    allNodes.insert( edgeNodesPrev[1] );
                                }
                                if (onEdge) {
                                    allNodes.insert( edgeNodes[0] );
                                    allNodes.insert( edgeNodes[1] );
                                }
                                if (onFacePrev) {
                                    allNodes.insert( faceNodesPrev[0] );
                                    allNodes.insert( faceNodesPrev[1] );
                                    allNodes.insert( faceNodesPrev[2] );
                                }
                                if (onFace) {
                                    allNodes.insert( faceNodes[0] );
                                    allNodes.insert( faceNodes[1] );
                                    allNodes.insert( faceNodes[2] );
                                }
                                if (secondNodes){
                                    auto no=--allNodes.end();
                                    while (no!=allNodes.begin()) {
                                        if(nodes[*no].getprimary()==10){
                                            std::copy(nodes[*no].getPrincipals().begin(),nodes[*no].getPrincipals().end(), std::inserter(allNodes,allNodes.end()));
                                            allNodes.erase(no);
                                        }
                                        --no;
                                    }
                                }
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
                            
                            // find next cell
                            T2 PrevCell=cellNo;
                            cellNo = findAdjacentCell1(faceNodes, nodeNo);
                            if ( cellNo == std::numeric_limits<T2>::max() ) {
                                InLimits=true;
                                cellNo=PrevCell;
                            }
                            break;
                        }
                        
                        if (foundIntersection)
                            break;
                    }
                    
                    if (InLimits && ! foundIntersection){
                        // if the current point is on limits but we don't find intersection, we project again the gradient.
                        foundIntersection=true;

                    }
                    if ( foundIntersection == false && !InLimits) {
                        std::cout << "\n\nWarning: finding raypath on node failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }

                }
                if (!foundIntersection)//
                    InLimits=InLimits?false:true;

            } else if ( onEdge ) {

                // find cells common to edge
                std::vector<T2> cells;
                std::set<NODE*> nnodes;
                std::vector<sxyz<T1>> Gradients;
                std::array<T2,2> edgeNodes2;
                std::array<T2,2> edgeNodes1;
                for(T2 Edge=0;Edge<2;++Edge){
                    for ( auto nc0=nodes[edgeNodes[Edge]].getOwners().begin(); nc0!=nodes[edgeNodes[Edge]].getOwners().end(); ++nc0 ) {
                        T2 Celli=*nc0;
                        for(T2 n=0;n<4;++n){
                            for(auto nc=nodes[neighbors[Celli][n]].getOwners().begin();nc!=nodes[neighbors[Celli][n]].getOwners().end();++nc){
                                //get cells common to edge
                                for(T2 iD=0;iD<4;++iD){
                                    T2 iA =neighbors[(*nc)][((iD+1)%4)];
                                    T2 iB =neighbors[(*nc)][((iD+2)%4)];
                                    T2 iC =neighbors[(*nc)][((iD+3)%4)];
                                    if ((nodes[iA].getDistance(curr_pt)+nodes[iB].getDistance(curr_pt)-nodes[iA].getDistance(nodes[iB])<minDist*minDist)||
                                        (nodes[iA].getDistance(curr_pt)+nodes[iC].getDistance(curr_pt)-nodes[iA].getDistance(nodes[iC])<minDist*minDist)||
                                        (nodes[iB].getDistance(curr_pt)+nodes[iC].getDistance(curr_pt)-nodes[iB].getDistance(nodes[iC])<minDist*minDist)){
                                        if (std::find(cells.begin(),cells.end(),(*nc))==cells.end())
                                            cells.push_back( (*nc) );
                                        getNeighborNodes(*nc, nnodes);
                                    }
                                }
                            }
                        }
                    }
                }
                if (!InLimits){
                    std::array<T1, 2> Weights={curr_pt.getDistance(nodes[edgeNodes[1]])/nodes[edgeNodes[0]].getDistance(nodes[edgeNodes[1]]),curr_pt.getDistance(nodes[edgeNodes[0]])/nodes[edgeNodes[0]].getDistance(nodes[edgeNodes[1]])};
                    T1 t=nodes[edgeNodes[0]].getTT(threadNo)*Weights[0]+nodes[edgeNodes[1]].getTT(threadNo)*Weights[1];
                    g= grad3d.ls_grad(nnodes, threadNo,t,curr_pt);
                    for(size_t nt=0;nt<Tx.size();++nt){
                        for(auto n=neighbors[txCell[nt]].begin();n!=neighbors[txCell[nt]].begin()+4;++n){
                            bool Find=false;
                            for(auto nc=nodes[*n].getOwners().begin();nc!=nodes[*n].getOwners().end();++nc){
                                if(*nc==cellNo){
                                    g=Tx[nt]-curr_pt;
                                    Find=true;
                                    break;
                                }
                            }
                            if (Find)
                                break;
                        }
                        
                    }
                    Gradients.push_back(g);
                }else{
                    std::array<T2,3> InBorderFace;
                    for (size_t n=0; n<cells.size(); ++n ) {
                        for ( T2 nn=0; nn<4; ++nn ){
                            InBorderFace[0]=neighbors[cells[n]][(nn+1)%4];
                            InBorderFace[1]=neighbors[cells[n]][(nn+2)%4];
                            InBorderFace[2]=neighbors[cells[n]][(nn+3)%4];
                            if (findAdjacentCell2(InBorderFace, cells[n], curr_pt)==std::numeric_limits<T2>::max() &&
                                testInTriangle(&nodes[InBorderFace[0]], &nodes[InBorderFace[1]], &nodes[InBorderFace[2]], curr_pt)){
                                Gradients.push_back(reflectedGradient(g, InBorderFace,  cells[n]));
                            }
                        }
                    }
                    if (Gradients.size()>1){// take gradient with biggest norm
                        if(norm(Gradients[0])<norm(Gradients[1])){
                            sxyz<T1> Permute=Gradients[0];
                            Gradients[0]=Gradients[1];
                            Gradients[1]=Permute;
                        }
                    }
                }
                bool foundIntersection=false;
                for(T2 ng=0; ng<Gradients.size();++ng){
                    g=Gradients[ng];
                    for (size_t n=0; n<cells.size(); ++n ) {
                        cellNo = cells[n];
                        // there are 2 faces that might be intersected
                        size_t n2=0;
                        size_t n1=0;
                        for ( auto nn=neighbors[cellNo].begin(); nn!= neighbors[cellNo].begin()+4; nn++ ) {
                            if ((*nn)==edgeNodes[0] || (*nn)==edgeNodes[1]||nodes[edgeNodes[0]].getDistance(nodes[*nn])+
                                nodes[edgeNodes[1]].getDistance(nodes[*nn])-nodes[edgeNodes[0]].getDistance(nodes[edgeNodes[1]])<minDist*minDist){
                                edgeNodes1[n1++] = *nn;// Edge contains current point
                            }else{
                                edgeNodes2[n2++] = *nn; // opposite edge
                            }
                        }
                        
                        sxyz<T1> pt_i;
                        T2 itmpNode;
                            foundIntersection = intersectVecTriangle(curr_pt, Gradients[ng],
                                                                     edgeNodes1[0],
                                                                     edgeNodes2[0],
                                                                     edgeNodes2[1], pt_i);
                            itmpNode = edgeNodes1[0];
                            if ( !foundIntersection ) {
                                foundIntersection = intersectVecTriangle(curr_pt, Gradients[ng],
                                                                         edgeNodes1[1],
                                                                         edgeNodes2[0],
                                                                         edgeNodes2[1], pt_i);
                                itmpNode = edgeNodes1[1];
                            }
                        if ( !foundIntersection ) {
                            continue;
                        }
                        InLimits=false;
                        prev_pt = curr_pt;
                        curr_pt = pt_i;
                        ds = curr_pt.getDistance( prev_pt );
                        r_tmp.push_back( curr_pt );
                        if (r_tmp.size() > 1 && ds >minDist) {
                            // compute terms of matrix M
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            s=computeSlowness(mid_pt,cellNo );
                            s *= s;
                        }
                        
                        bool break_flag = false;
                        for ( size_t n2=0; n2<4; ++n2 ) {
                            if ( nodes[ neighbors[cellNo][n2] ].getDistance( curr_pt ) < small*small ) {
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;
                                edgeNodesPrev = edgeNodes;
                                nodeNo = neighbors[cellNo][n2];
                                onNode = true;
                                onEdge = false;
                                onFace = false;
                                break_flag = true;
                                if ( r_tmp.size() > 1 && ds>minDist) {
                                    std::set<T2> allNodes;
                                    allNodes.insert(nodeNo);
                                    
                                    allNodes.insert( edgeNodesPrev[0] );
                                    allNodes.insert( edgeNodesPrev[1] );
                                    if (secondNodes){
                                        auto no=--allNodes.end();
                                        while (no!=allNodes.begin()) {
                                            if(nodes[*no].getprimary()==10){
                                                std::copy(nodes[*no].getPrincipals().begin(),nodes[*no].getPrincipals().end(), std::inserter(allNodes,allNodes.end()));
                                                allNodes.erase(no);
                                            }
                                            --no;
                                        }
                                    }
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
                                
                                break;
                            }
                        }
                        if ( break_flag ) break;
                        
                        if ((nodes[itmpNode].getDistance(curr_pt)+nodes[edgeNodes2[0]].getDistance(curr_pt)-nodes[ itmpNode].getDistance(nodes[edgeNodes2[0]])<minDist*minDist)) {
                            edgeNodesPrev = edgeNodes;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            edgeNodes[0] = itmpNode;
                            edgeNodes[1] = edgeNodes2[0];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            
                            if ( r_tmp.size() > 1 && ds>minDist) {
                                std::set<T2> allNodes;
                                allNodes.insert( edgeNodesPrev[0] );
                                allNodes.insert( edgeNodesPrev[1] );
                                
                                allNodes.insert( edgeNodes[0] );
                                allNodes.insert( edgeNodes[1] );
                                if (secondNodes){
                                    auto no=--allNodes.end();
                                    while (no!=allNodes.begin()) {
                                        if(nodes[*no].getprimary()==10){
                                            std::copy(nodes[*no].getPrincipals().begin(),nodes[*no].getPrincipals().end(), std::inserter(allNodes,allNodes.end()));
                                            allNodes.erase(no);
                                        }
                                        --no;
                                    }
                                }
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
                            
                            break_flag = true;
                            break;
                        } else if((nodes[itmpNode].getDistance(curr_pt)+nodes[edgeNodes2[1]].getDistance(curr_pt)-nodes[ itmpNode].getDistance(nodes[edgeNodes2[1]])<minDist*minDist)){
                            edgeNodesPrev = edgeNodes;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            edgeNodes[0] = itmpNode;
                            edgeNodes[1] = edgeNodes2[1];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            
                            if ( r_tmp.size() > 1  && ds>minDist) {
                                std::set<T2> allNodes;
                                allNodes.insert( edgeNodesPrev[0] );
                                allNodes.insert( edgeNodesPrev[1] );
                                
                                allNodes.insert( edgeNodes[0] );
                                allNodes.insert( edgeNodes[1] );
                                if (secondNodes){
                                    auto no=--allNodes.end();
                                    while (no!=allNodes.begin()) {
                                        if(nodes[*no].getprimary()==10){
                                            std::copy(nodes[*no].getPrincipals().begin(),nodes[*no].getPrincipals().end(), std::inserter(allNodes,allNodes.end()));
                                            allNodes.erase(no);
                                        }
                                        --no;
                                    }
                                }
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
                            
                            break_flag = true;
                            break;
                        } else if ((nodes[edgeNodes2[1]].getDistance(curr_pt)+nodes[edgeNodes2[0]].getDistance(curr_pt)-nodes[edgeNodes2[1]].getDistance(nodes[edgeNodes2[0]])<minDist*minDist)){
                            edgeNodesPrev = edgeNodes;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;
                            
                            edgeNodes[0] = edgeNodes2[0];
                            edgeNodes[1] = edgeNodes2[1];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            
                            if ( r_tmp.size() > 1 && ds>minDist) {
                                std::set<T2> allNodes;
                                allNodes.insert( edgeNodesPrev[0] );
                                allNodes.insert( edgeNodesPrev[1] );
                                
                                allNodes.insert( edgeNodes[0] );
                                allNodes.insert( edgeNodes[1] );
                                if (secondNodes){
                                    auto no=--allNodes.end();
                                    while (no!=allNodes.begin()) {
                                        if(nodes[*no].getprimary()==10){
                                            std::copy(nodes[*no].getPrincipals().begin(),nodes[*no].getPrincipals().end(), std::inserter(allNodes,allNodes.end()));
                                            allNodes.erase(no);
                                        }
                                        --no;
                                    }
                                }
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
                            
                            break_flag = true;
                            break;
                        }
                        if ( break_flag ) break;
                        
                        onNodePrev = onNode;
                        onEdgePrev = onEdge;
                        onFacePrev = onFace;
                        onNode = false;
                        onEdge = false;
                        onFace = true;
                        
                        //faceNodesPrev = faceNodes;
                        edgeNodesPrev = edgeNodes;
                        faceNodes[0] = itmpNode;
                        faceNodes[1] = edgeNodes2[0];
                        faceNodes[2] = edgeNodes2[1];
                        std::sort(faceNodes.begin(), faceNodes.end());
                        
                        if ( r_tmp.size() > 1 && ds>minDist) {
                            std::set<T2> allNodes;
                            allNodes.insert( edgeNodesPrev[0] );
                            allNodes.insert( edgeNodesPrev[1] );
                            
                            allNodes.insert( faceNodes[0] );
                            allNodes.insert( faceNodes[1] );
                            allNodes.insert( faceNodes[2] );
                            if (secondNodes){
                                auto no=--allNodes.end();
                                while (no!=allNodes.begin()) {
                                    if(nodes[*no].getGridIndex()>nPrimary){
                                        std::copy(nodes[*no].getPrincipals().begin(),nodes[*no].getPrincipals().end(),std::inserter(allNodes,allNodes.end()));
                                        allNodes.erase(no);
                                    }
                                    --no;
                                }
                            }
                            std::vector<T1> w;
                            T1 sum_w = 0.0;
                            for ( auto it=allNodes.begin(); it!=allNodes.end(); ++it ) {
                                if (*it<nodes.size()){
                                    w.push_back( 1./nodes[*it].getDistance( mid_pt ) );
                                    sum_w += w.back();
                                }
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
                        
                        // find next cell
                        T2 PrevCell(cellNo);
                        cellNo=findAdjacentCell2(faceNodes, cellNo,curr_pt);
                        
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            InLimits=true;
                            cellNo=PrevCell;
                        }
                        break;
                    }
                    if (foundIntersection)
                        break;
                }
                if (InLimits && ! foundIntersection){
                    foundIntersection=true;// we project again th gradient

                }
                if (!foundIntersection){
                    InLimits=InLimits? false:true;
                }
                if ( foundIntersection == false && InLimits==false) {
                    std::cout << "\n\nWarning: finding raypath on edge failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    reachedTx = true;
                }
            } else{ // on Face
                bool foundIntersection=false;
                std::set<NODE*> nnodes;
                if (!InLimits){
                    getNeighborNodes(cellNo, nnodes);
                    T1 t;
                    if (r_tmp.size()<=1){
                        t=Interpolator<T1>::TrilinearTime(curr_pt, nodes[neighbors[cellNo][0]], nodes[neighbors[cellNo][1]], nodes[neighbors[cellNo][2]], nodes[neighbors[cellNo][3]], threadNo);
                        
                    }else{
                        t=Interpolator<T1>::bilinearTime(curr_pt, nodes[faceNodes[0]], nodes[faceNodes[1]],  nodes[faceNodes[2]], threadNo);
                    }
                    g = grad3d.ls_grad(nnodes,threadNo,t,curr_pt);
                    for(size_t nt=0;nt<Tx.size();++nt){
                        for(auto n=neighbors[txCell[nt]].begin();n!=neighbors[txCell[nt]].begin()+4;++n){
                            bool Find=false;
                            for(auto nc=nodes[*n].getOwners().begin();nc!=nodes[*n].getOwners().end();++nc){
                                if(*nc==cellNo){
                                    g=Tx[nt]-curr_pt;
                                    Find=true;
                                    break;
                                }
                            }
                            if (Find)
                                break;
                        }

                    }
                }else{
                    g=reflectedGradient(g, faceNodes,cellNo);
                }
                //std::cout<<"nomr g : "<<norm(g)<<std::endl;

                std::array<T2,3> ind[4] = {
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][2] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][0], neighbors[cellNo][2], neighbors[cellNo][3] } },
                    { { neighbors[cellNo][1], neighbors[cellNo][2], neighbors[cellNo][3] } }
                };

                for ( size_t n=0; n<4; ++n )
                    std::sort( ind[n].begin(), ind[n].end() );
                // there are 3 faces that might be intersected
                for ( size_t n=0; n<4; ++n ) {
                    if ( ind[n] == faceNodes || AreSameFace(ind[n], faceNodes)) continue;

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);
                    
                    if ( !foundIntersection )
                        continue;
                    InLimits=false;
                    prev_pt = curr_pt;
                    curr_pt = pt_i;
                    ds = curr_pt.getDistance( prev_pt );
                    r_tmp.push_back( curr_pt );
                    if (r_tmp.size() > 1 && ds>minDist) {
                        // compute terms of matrix M
                        mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                        s=computeSlowness(mid_pt,cellNo );
                        s *= s;

                    }

                    bool break_flag = false;
                    for ( size_t n2=0; n2<3; ++n2 ) {
                        if ( nodes[ ind[n][n2] ].getDistance( curr_pt ) < small*small ) {
                            //nodeNoPrev = nodeNo;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;

                            nodeNo = ind[n][n2];
                            onNode = true;
                            onEdge = false;
                            onFace = false;
                            faceNodesPrev = faceNodes;
                            if ( r_tmp.size() > 1 && ds>minDist) {
                                std::set<T2> allNodes;
                                allNodes.insert(nodeNo);
                                allNodes.insert( faceNodesPrev[0] );
                                allNodes.insert( faceNodesPrev[1] );
                                allNodes.insert( faceNodesPrev[2] );
                                if (secondNodes){
                                    auto no=--allNodes.end();
                                    while (no!=allNodes.begin()) {
                                        if(nodes[*no].getprimary()==10){
                                            std::copy(nodes[*no].getPrincipals().begin(),nodes[*no].getPrincipals().end(), std::inserter(allNodes,allNodes.end()));
                                            allNodes.erase(no);
                                        }
                                        --no;
                                    }
                                }
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

                            break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;

                    for ( size_t n1=0; n1<3; ++n1 ) {
                        size_t n2 = (n1+1)%3;
                        if (nodes[ind[n][n1]].getDistance(curr_pt)+nodes[ind[n][n2]].getDistance(curr_pt)-nodes[ind[n][n1]].getDistance(nodes[ind[n][n2]])<minDist*minDist) {
                            //edgeNodesPrev = edgeNodes;
                            onNodePrev = onNode;
                            onEdgePrev = onEdge;
                            onFacePrev = onFace;

                            edgeNodes[0] = ind[n][n1];
                            edgeNodes[1] = ind[n][n2];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            faceNodesPrev = faceNodes;
                            if ( r_tmp.size() > 1 && ds>minDist) {
                                std::set<T2> allNodes;
                                allNodes.insert( edgeNodes[0] );
                                allNodes.insert( edgeNodes[1] );
                                allNodes.insert( faceNodesPrev[0] );
                                allNodes.insert( faceNodesPrev[1] );
                                allNodes.insert( faceNodesPrev[2] );
                                if (secondNodes){
                                    auto no=--allNodes.end();
                                    while (no!=allNodes.begin()) {
                                        if(nodes[*no].getprimary()==10){
                                            std::copy(nodes[*no].getPrincipals().begin(),nodes[*no].getPrincipals().end(), std::inserter(allNodes,allNodes.end()));
                                            allNodes.erase(no);
                                        }
                                        --no;
                                    }
                                }
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

                            break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;

                    onNodePrev = onNode;
                    onEdgePrev = onEdge;
                    onFacePrev = onFace;
                    onNode = false;
                    onEdge = false;
                    onFace = true;

                    faceNodesPrev = faceNodes;
                    
                    faceNodes = ind[n];
                    if ( r_tmp.size() > 1 && ds>minDist) {
                        std::set<T2> allNodes;
                        if (r_tmp.size()>2){
                            allNodes.insert( faceNodesPrev[0] );
                            allNodes.insert( faceNodesPrev[1] );
                            allNodes.insert( faceNodesPrev[2] );
                            allNodes.insert( faceNodes[0] );
                            allNodes.insert( faceNodes[1] );
                            allNodes.insert( faceNodes[2] );
                        }else{
                            for (auto nc=neighbors[cellNo].begin();nc!=neighbors[cellNo].begin()+4;++nc){
                                allNodes.insert(*nc);
                            }
                        }
                        if (secondNodes){
                            auto no=--allNodes.end();
                            while (no!=allNodes.begin()) {
                                if(nodes[*no].getprimary()==10){
                                    std::copy(nodes[*no].getPrincipals().begin(),nodes[*no].getPrincipals().end(), std::inserter(allNodes,allNodes.end()));
                                    allNodes.erase(no);
                                }
                                --no;
                            }
                        }
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

                    // find next cell
                    T2 PrevCell(cellNo);
                    cellNo=findAdjacentCell2(faceNodes, cellNo,curr_pt);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        InLimits=true;
                        cellNo=PrevCell;
                    }
                    break;
                }
                if ( foundIntersection == false ) {

                    // we must be on an face with gradient pointing slightly outward tetrahedron
                    // return in other cell but keep gradient
                    T2 PrevCell(cellNo);
                    cellNo=findAdjacentCell2(faceNodes, cellNo,curr_pt);
                    if (cellNo==std::numeric_limits<T2>::max())
                        cellNo=PrevCell;
                    std::array<T2,3> ind[4];
                    ind[0] = { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][2] } };
                    ind[1] = { { neighbors[cellNo][0], neighbors[cellNo][1], neighbors[cellNo][3] } };
                    ind[2] = { { neighbors[cellNo][0], neighbors[cellNo][2], neighbors[cellNo][3] } };
                    ind[3] = { { neighbors[cellNo][1], neighbors[cellNo][2], neighbors[cellNo][3] } };

                    for ( size_t n=0; n<4; ++n )
                        std::sort( ind[n].begin(), ind[n].end() );

                    for ( size_t n=0; n<4; ++n ) {
                        if ( ind[n] == faceNodes || AreSameFace(ind[n], faceNodes)) continue;

                        sxyz<T1> pt_i;
                        foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                                 ind[n][1], ind[n][2],
                                                                 pt_i);

                        if ( !foundIntersection ) {
                            continue;
                        }
                        prev_pt = curr_pt;
                        InLimits=false;
                        curr_pt = pt_i;
                        ds = curr_pt.getDistance( prev_pt );
                        r_tmp.push_back( curr_pt );
                        if (r_tmp.size() > 1 && ds>minDist) {
                            // compute terms of matrix M
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            s = computeSlowness(mid_pt,cellNo);
                            s *= s;
                        }

                        bool break_flag = false;
                        for ( size_t n2=0; n2<3; ++n2 ) {
                            if ( nodes[ ind[n][n2] ].getDistance( curr_pt ) < small*small ) {
                                //nodeNoPrev = nodeNo;
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;

                                nodeNo = ind[n][n2];
                                onNode = true;
                                onEdge = false;
                                onFace = false;
                                faceNodesPrev = faceNodes;
                                if ( r_tmp.size() > 1 && ds>minDist) {
                                    std::set<T2> allNodes;
                                    allNodes.insert(nodeNo);
                                    allNodes.insert( faceNodesPrev[0] );
                                    allNodes.insert( faceNodesPrev[1] );
                                    allNodes.insert( faceNodesPrev[2] );
                                    if (secondNodes){
                                        auto no=--allNodes.end();
                                        while (no!=allNodes.begin()) {
                                            if(nodes[*no].getprimary()==10){
                                                std::copy(nodes[*no].getPrincipals().begin(),nodes[*no].getPrincipals().end(), std::inserter(allNodes,allNodes.end()));
                                                allNodes.erase(no);
                                            }
                                            --no;
                                        }
                                    }
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

                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;

                        for ( size_t n1=0; n1<3; ++n1 ) {
                            size_t n2 = (n1+1)%3;
                            if (nodes[ind[n][n1]].getDistance(curr_pt)+nodes[ind[n][n2]].getDistance(curr_pt)-nodes[ind[n][n1]].getDistance(nodes[ind[n][n2]])<minDist*minDist) {
                                //edgeNodesPrev = edgeNodes;
                                onNodePrev = onNode;
                                onEdgePrev = onEdge;
                                onFacePrev = onFace;

                                edgeNodes[0] = ind[n][n1];
                                edgeNodes[1] = ind[n][n2];
                                onNode = false;
                                onEdge = true;
                                onFace = false;
                                faceNodesPrev = faceNodes;
                                if ( r_tmp.size() > 1 && ds>minDist) {
                                    std::set<T2> allNodes;
                                    allNodes.insert( edgeNodes[0] );
                                    allNodes.insert( edgeNodes[1] );
                                    
                                    allNodes.insert( faceNodesPrev[0] );
                                    allNodes.insert( faceNodesPrev[1] );
                                    allNodes.insert( faceNodesPrev[2] );
                                    if (secondNodes){
                                        auto no=--allNodes.end();
                                        while (no!=allNodes.begin()) {
                                            if(nodes[*no].getprimary()==10){
                                                std::copy(nodes[*no].getPrincipals().begin(),nodes[*no].getPrincipals().end(), std::inserter(allNodes,allNodes.end()));
                                                allNodes.erase(no);
                                            }
                                            --no;
                                        }
                                    }
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

                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;

                        onNodePrev = onNode;
                        onEdgePrev = onEdge;
                        onFacePrev = onFace;
                        onNode = false;
                        onEdge = false;
                        onFace = true;

                        faceNodesPrev = faceNodes;
                        faceNodes = ind[n];
                        if ( r_tmp.size() > 1 && ds>minDist) {
                            std::set<T2> allNodes;
                            allNodes.insert( faceNodesPrev[0] );
                            allNodes.insert( faceNodesPrev[1] );
                            allNodes.insert( faceNodesPrev[2] );
                            
                            allNodes.insert( faceNodes[0] );
                            allNodes.insert( faceNodes[1] );
                            allNodes.insert( faceNodes[2] );
                            if (secondNodes){
                                auto no=--allNodes.end();
                                while (no!=allNodes.begin()) {
                                    if(nodes[*no].getprimary()==10){
                                        std::copy(nodes[*no].getPrincipals().begin(),nodes[*no].getPrincipals().end(), std::inserter(allNodes,allNodes.end()));
                                        allNodes.erase(no);
                                    }
                                    --no;
                                }
                            }
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

                        // find next cell
                        T2 PrevCell(cellNo);
                        cellNo=findAdjacentCell2(faceNodes, cellNo,curr_pt);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            InLimits=true;
                            cellNo=PrevCell;
                            
                        }
                        break;
                    }
            }
                if (!foundIntersection)
                    InLimits=InLimits? false:true;
                if ( foundIntersection == false && InLimits==false ) {
                    std::cout << "\n\nWarning: finding raypath on face failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    reachedTx = true;
                }
            }
            for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                if ( curr_pt.getDistance( Tx[nt] ) < minDist*minDist) {
                    reachedTx = true;
                    break;
                }
            }
            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( curr_pt.getDistance( Tx[nt] ) < minDist*minDist ) {
                        reachedTx = true;
                        break;
                    }
                    for(size_t n=0;n<4;++n){
                        sxyz<T1> NearTx={nodes[neighbors[txCell[nt]][n]]};
                        if (curr_pt.getDistance(NearTx)< minDist*minDist){
                            reachedTx = true;
                            r_tmp.push_back( Tx[nt] );
                            prev_pt=curr_pt;
                            curr_pt=Tx[nt];
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
                                r_tmp.push_back( Tx[nt] );
                                reachedTx = true;
                                prev_pt=curr_pt;
                                curr_pt=Tx[nt];
                                break;
                            }
                        }
                    } else {
                        if ( cellNo == txCell[nt] ) {
                            r_tmp.push_back( Tx[nt] );
                            reachedTx = true;
                            prev_pt=curr_pt;
                            curr_pt=Tx[nt];
                        }
                        
                        else {
                            for ( size_t nn=0; nn<txNeighborCells[nt].size(); ++nn ) {
                                if ( cellNo == txNeighborCells[nt][nn] ) {
                                    bool find=false;
                                    for (T2 i=0;i<4;++i){
                                        T2 iA=neighbors[txCell[nt]][(i+1)%4];
                                        T2 iB=neighbors[txCell[nt]][(i+2)%4];
                                        T2 iC=neighbors[txCell[nt]][(i+3)%4];
                                        if (nodes[iA].getDistance(curr_pt)+nodes[iB].getDistance(curr_pt)-nodes[iB].getDistance(nodes[iA])<minDist*minDist||nodes[iA].getDistance(curr_pt)+nodes[iC].getDistance(curr_pt)-nodes[iC].getDistance(nodes[iA])<minDist*minDist||nodes[iC].getDistance(curr_pt)+nodes[iB].getDistance(curr_pt)-nodes[iB].getDistance(nodes[iC])<minDist*minDist){
                                            find=true;
                                            break;
                                        }
                                    }
                                    if (find){
                                        r_tmp.push_back( Tx[nt] );
                                        prev_pt=curr_pt;
                                        curr_pt=Tx[nt];
                                        reachedTx = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if ( reachedTx ){
                        if (r_tmp.size() > 1 && ds>minDist) {
                            // compute terms of matrix M
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            s = computeSlowness(mid_pt,cellNo);
                            s *= s;
                        }
                        if ( r_tmp.size() > 1 && ds>minDist) {
                            std::set<T2> allNodes;
                            for (T2 i=0;i<4;++i)
                                allNodes.insert(neighbors[cellNo][i]);
                            if (secondNodes){
                                auto no=--allNodes.end();
                                while (no!=allNodes.begin()) {
                                    if(nodes[*no].getprimary()==10){
                                        std::copy(nodes[*no].getPrincipals().begin(),nodes[*no].getPrincipals().end(), std::inserter(allNodes,allNodes.end()));
                                        allNodes.erase(no);
                                    }
                                    --no;
                                }
                            }
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
                        break;}
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
    bool Grid3Dui<T1,T2,NODE>::intersectVecTriangle(const T2 iO, const sxyz<T1> &vec,
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
    bool Grid3Dui<T1,T2,NODE>::intersectVecTriangle(const sxyz<T1> &O, const sxyz<T1> &vect,
                                                    const T2 iA, T2 iB, T2 iC,
                                                    sxyz<T1> &pt_i) const {

        sxyz<T1> vec(vect);
        vec/=norm(vec);
        vec*=10;
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
        if ( -u > small*small) return false;
        v = tripleScalar(vec, OA, OC);
        if ( -v > small*small) return false;
        w = tripleScalar(vec, OB, OA);
        if ( -w > small*small) return false;

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
    bool Grid3Dui<T1,T2,NODE>::areCollinear(const sxyz<T1> &pt, const T2 i0, const T2 i1) const {

        // http://mathworld.wolfram.com/Collinear.html
        //
        sxyz<T1> v1 = {pt.x-nodes[i0].getX(), pt.y-nodes[i0].getY(), pt.z-nodes[i0].getZ()};
        sxyz<T1> v2 = {pt.x-nodes[i1].getX(), pt.y-nodes[i1].getY(), pt.z-nodes[i1].getZ()};
        sxyz<T1> v3 = cross(v1, v2);
        return norm(v3)<small*small;

    }

    template<typename T1, typename T2, typename NODE>
    bool Grid3Dui<T1,T2,NODE>::areCoplanar(const sxyz<T1> &x1, const T2 i0, const T2 i1, const T2 i2) const {

        // http://mathworld.wolfram.com/Coplanar.html
        //
        sxyz<T1> x2 = {nodes[i0].getX(), nodes[i0].getY(), nodes[i0].getZ()};
        sxyz<T1> x3 = {nodes[i1].getX(), nodes[i1].getY(), nodes[i1].getZ()};
        sxyz<T1> x4 = {nodes[i2].getX(), nodes[i2].getY(), nodes[i2].getZ()};

        return fabs( dot( x3-x1, cross(x2-x1, x4-x3) ) )<small*small*small;
    }
    template<typename T1, typename T2, typename NODE>
    bool Grid3Dui<T1,T2,NODE>::areCoplanar(const sxyz<T1> &x1, const sxyz<T1> &x2
                                           , const sxyz<T1> &x3, const sxyz<T1> &x4) const {
        return (fabs( dot( x3-x1, cross(x2-x1, x4-x3) ) )<small*small);
    }


    template<typename T1, typename T2, typename NODE>
    T2 Grid3Dui<T1,T2,NODE>::findAdjacentCell1(const std::array<T2,3> &faceNodes,
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
    T2 Grid3Dui<T1,T2,NODE>::findAdjacentCell2(const std::array<T2,3> &faceNodes,
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
        if(cells.size()==2){
            if ( cellNo == cells[0] ) {
                return cells[1];
            } else if ( cellNo == cells[1] ) {
                return cells[0];
            }
        }
        return std::numeric_limits<T2>::max();
    }
    template<typename T1, typename T2, typename NODE>
    void Grid3Dui<T1,T2,NODE>::findAdjacentCell2(const std::array<T2,3> &faceNodes,
                                                 const T2 & cellNo, std::set<T2> & AdjacentCells) const {
        
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
    }
    template<typename T1, typename T2, typename NODE>
    T2 Grid3Dui<T1,T2,NODE>::findAdjacentCell2(const std::array<T2,3> &faceNodes,
                                                 const T2 & cellNo, const sxyz<T1> & curr_pt )const {
        
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
    void Grid3Dui<T1,T2,NODE>::plotCell(const T2 cellNo, const sxyz<T1> &pt, const sxyz<T1> &g) const {


        if ( true ) {
            T2 i0 = neighbors[cellNo][0];
            T2 i1 = neighbors[cellNo][1];
            T2 i2 = neighbors[cellNo][2];
            T2 i3 = neighbors[cellNo][3];

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
            std::cout << "plot3(["<<pt.x<< ' ' << pt.x+g.x<<"],["<<pt.y<< ' ' << pt.y+g.y<<"],["<<pt.z<< ' ' << pt.z+g.z<<"],'*-r')\naxis equal\n\n";
        }
    }
    template<typename T1, typename T2, typename NODE>
    void Grid3Dui<T1,T2,NODE>::getNeighborNodes(const T2 cellNo,
                                                std::set<NODE*> &nnodes) const {

        for ( size_t n=0; n<4; ++n ) {
            T2 nodeNo = neighbors[cellNo][n];
            nnodes.insert( &(nodes[nodeNo]) );

            for ( auto nc=nodes[nodeNo].getOwners().cbegin(); nc!=nodes[nodeNo].getOwners().cend(); ++nc ) {
                for ( size_t nn=0; nn<4; ++nn ) {
                    //if (nodes[ neighbors[*nc][nn]].getGridIndex()<nPrimary)
                        nnodes.insert( &(nodes[ neighbors[*nc][nn]]) );
                }
            }
        }
    }
    template<typename T1, typename T2, typename NODE>
    void Grid3Dui<T1,T2,NODE>::getNeighborNodes2(const T2 cellNo,
                                                std::set<NODE*> &nnodes) const {

        for ( size_t n=0; n<4; ++n ) {
            T2 nodeNo = neighbors[cellNo][n];
            nnodes.insert( &(nodes[nodeNo]) );

            for ( auto nc=nodes[nodeNo].getOwners().cbegin(); nc!=nodes[nodeNo].getOwners().cend(); ++nc ) {
                for ( T2 nn=0; nn<4; ++nn ) {
                    NODE Nod2(nodes[ neighbors[(*nc)][nn] ]);
                    for ( auto nnc=Nod2.getOwners().cbegin(); nnc!=Nod2.getOwners().cend(); ++nnc ) {
                        for(size_t N=0;N<4;++N){
                            if (nodes[ neighbors[*nc][nn]].getGridIndex()<nPrimary)
                                nnodes.insert( &(nodes[ neighbors[*nc][N]]));
                        }
                    }
                }
            }
        }
    }


    template<typename T1, typename T2, typename NODE>
    bool Grid3Dui<T1,T2,NODE>::testInTriangle(const NODE *vertexA,
                                              const NODE *vertexB,
                                              const NODE *vertexC,
                                              const sxyz<T1> &E) const {
        if (vertexA->getDistance(E)+vertexB->getDistance(E)-vertexA->getDistance(*vertexB)<small*small)
            return true;
        if (vertexA->getDistance(E)+vertexC->getDistance(E)-vertexA->getDistance(*vertexC)<small*small)
            return true;
        if (vertexC->getDistance(E)+vertexB->getDistance(E)-vertexC->getDistance(*vertexB)<small*small)
            return true;
        if (!areCoplanar(sxyz<T1>((*vertexA)), sxyz<T1>((*vertexB)), sxyz<T1>((*vertexC)), E))
            return false;
        T1 u, v, w;
        barycentric(vertexA, vertexB, vertexC, E, u, v, w);
        return (v) >=0.0 && (w )>= 0.0 && (v + w) <= 1.0;
    }
    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dui<T1,T2,NODE>::AreaTriangle(const NODE *vertexA,
                                              const NODE *vertexB,
                                              const NODE *vertexC) const {
        T1 a=vertexB->getDistance(*vertexC);
        T1 b=vertexA->getDistance(*vertexC);
        T1 c=vertexA->getDistance(*vertexB);
        T1 Area=0.25*sqrt((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c));
        return Area;
    }
    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dui<T1,T2,NODE>::AreaTriangle(const NODE *vertexA,
                                          const NODE *vertexB,
                                          const sxyz<T1>& E) const {
        T1 a=vertexB->getDistance(E);
        T1 b=vertexA->getDistance(E);
        T1 c=vertexA->getDistance(*vertexB);
        T1 Area=0.25*sqrt((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c));
        return Area;
    }
    template<typename T1, typename T2, typename NODE>
    void Grid3Dui<T1,T2,NODE>::barycentric(const NODE *a,
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
    T1 Grid3Dui<T1,T2,NODE>::computeSlowness( const sxyz<T1>& Rx ) const {

        //Calculate the slowness of any point that is not on a node

        T2 cellNo = this->getCellNo( Rx );

        return(Interpolator<T1>::TrilinearTriangleVel(Rx,nodes[neighbors[cellNo][0]] , nodes[neighbors[cellNo][1]], nodes[neighbors[cellNo][2]], nodes[neighbors[cellNo][3]]));


    }
    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dui<T1,T2,NODE>::Highorder_computeSlowness( const sxyz<T1>& Rx ) const {
        
        //Calculate the slowness of any point that is not on a node
        
        T2 cellNo = this->getCellNo( Rx );
        std::set<T2> nnodes;
        for (size_t i=0;i<4;++i){
            nnodes.insert(neighbors[cellNo][i]);
            for (auto nt=nodes[neighbors[cellNo][i]].getOwners().begin();nt!=nodes[neighbors[cellNo][i]].getOwners().end();++nt){
                nnodes.insert(neighbors[*nt][0]);
                nnodes.insert(neighbors[*nt][1]);
                nnodes.insert(neighbors[*nt][2]);
                nnodes.insert(neighbors[*nt][3]);
            }
        }
        Eigen::Matrix<T1, Eigen::Dynamic, 10> A;
        Eigen::Matrix<T1, 10, 1> x;
        Eigen::Matrix<T1, Eigen::Dynamic, 1> b;
        A.resize( nnodes.size(), 10);
        b.resize( nnodes.size(), 1 );
        size_t i=0;
        for ( auto n=nnodes.cbegin(); n!=nnodes.cend(); ++n ) {
            A(i,0)=1.0;
            A(i,1)=nodes[*n].getX();
            A(i,2)=nodes[*n].getY();
            A(i,3)=nodes[*n].getZ();
            A(i,4)=nodes[*n].getX()*nodes[*n].getY();
            A(i,5)=nodes[*n].getX()*nodes[*n].getZ();
            A(i,6)=nodes[*n].getZ()*nodes[*n].getY();
            //A(i,7)=nodes[*n].getX()*nodes[*n].getY()*nodes[*n].getZ();
            A(i,7)=nodes[*n].getX()*nodes[*n].getX();
            A(i,8)=nodes[*n].getY()*nodes[*n].getY();
            A(i,9)=nodes[*n].getZ()*nodes[*n].getZ();
            b(i,0)=nodes[*n].getNodeSlowness();
            i++;
        }
         x = A.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(b);
        return(x(0)+x(1)*Rx.x+x(2)*Rx.y+x(3)*Rx.z+x(4)*Rx.x*Rx.y+x(5)*Rx.x*Rx.z+x(6)*Rx.y*Rx.z+x(7)*Rx.x*Rx.x+x(8)*Rx.y*Rx.y+x(9)*Rx.z*Rx.z);
    }
    template <typename T1,typename T2,typename NODE>
    bool Grid3Dui<T1,T2,NODE>::AreSameFace(const std::array<T2, 3>& Face1,const std::array<T2, 3>& Face2) const{
        if (areCoplanar(sxyz<T1> (nodes[Face2[0]]), Face1[0], Face1[1], Face1[2])&&
            areCoplanar(sxyz<T1> (nodes[Face2[1]]), Face1[0], Face1[1], Face1[2])&&
            areCoplanar(sxyz<T1> (nodes[Face2[2]]), Face1[0], Face1[1], Face1[2]))
            return true;
        return false;
    }
    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dui<T1,T2,NODE>::computeDt(const NODE& source, const NODE& node) const{
        T1 slo=node.getNodeSlowness();
        return (computeDt(source, node,slo));
    }
    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dui<T1,T2,NODE>::computeDt(const NODE& source, const sxyz<T1>& node, T1 slo) const{
        
        return (source.getDistance(node)*0.5*(slo+source.getNodeSlowness()));
    }
    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dui<T1,T2,NODE>::computeSlowness( const sxyz<T1>& Rx,const T2& cellNo ) const {
        
        std::vector<T2> list;
        for (size_t n3=0; n3 < neighbors[ cellNo ].size(); n3++){
                list.push_back(neighbors[ cellNo ][n3]);
        }
        
        std::vector<NODE*> interpNodes;
        
        for ( size_t nn=0; nn<list.size(); ++nn )
            interpNodes.push_back( &(nodes[list[nn] ]) );
        
        return Interpolator<T1>::inverseDistance( Rx, interpNodes );
        
    }
     template<typename T1, typename T2, typename NODE>
    bool Grid3Dui<T1,T2,NODE>::blti_raytrace(const sxyz<T1> & curr_pt,const std::array<T2,3> &face, sxyz<T1> & next_pt,const size_t threadNo,const T1 & s)const{
        
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
    bool Grid3Dui<T1,T2,NODE>::blti2D_raytrace(const sxyz<T1> & curr_pt,const T2 & node1,const T2 &node2, sxyz<T1> & next_pt,const size_t threadNo, const T1 & s) const{
        
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
        //
        //T1 D=vertexA->getDistance(pt)+vertexB->getDistance(pt)-vertexA->getDistance(*vertexB);
        T1 rho0 = curr_pt.getDistance( pt );
        //        T1 xi0 = vertexA->getDistance( pt )/c;
        T1 xi0 = k;
        T1 xi = xi0 -u*rho0/(w*c);
        
        if ( 0.0<xi && xi<1.0 ) {
            next_pt.x=(1-xi)*vertexA->getX()+xi*vertexB->getX();
            next_pt.y=(1-xi)*vertexA->getY()+xi*vertexB->getY();
            next_pt.z=(1-xi)*vertexA->getZ()+xi*vertexB->getZ();
//            T1 d=(pt.getDistance(*vertexA)+pt.getDistance(*vertexB))- vertexA->getDistance(*vertexB);
            return true;
        }
        return false;
    }
    template<typename T1, typename T2, typename NODE>
    bool Grid3Dui<T1,T2,NODE>::BLTISolver_ArroundSource(const sxyz<T1>& Source,const sxyz<T1>& curr_pt,const std::array<T2, 3>& face, std::array<T1, 3>& barycenters) const{
        T2 node1=face[0];
        T2 node2=face[1];
        T2 node3=face[2];
        sxyz<T1> PQ=Source-curr_pt;
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
}

#endif

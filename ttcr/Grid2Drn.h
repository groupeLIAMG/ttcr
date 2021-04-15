//
//  Grid2Drn.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-04-22.
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

#ifndef ttcr_Grid2Drn_h
#define ttcr_Grid2Drn_h

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifdef VTK
#include "vtkDoubleArray.h"
#include "vtkPolyData.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLRectilinearGridWriter.h"
#endif

#include <boost/math/special_functions/sign.hpp>

#include "Grid2D.h"

namespace ttcr {

    template<typename T1, typename T2, typename S, typename NODE>
    class Grid2Drn : public Grid2D<T1,T2,S> {
    public:
        Grid2Drn(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
                 const T1 minx, const T1 minz, const bool ttrp,
                 const size_t nt=1) :
        Grid2D<T1,T2,S>(nx*nz, ttrp, nt),
        dx(ddx), dz(ddz), xmin(minx), zmin(minz),
        xmax(minx+nx*ddx), zmax(minz+nz*ddz),
        ncx(nx), ncz(nz),
        nodes(std::vector<NODE>( (ncx+1) * (ncz+1), NODE(nt) ))
        {
        }

        virtual ~Grid2Drn() {
        }

        void setSlowness(const std::vector<T1>& s) {
            if ( nodes.size() != s.size() ) {
                throw std::length_error("Error: slowness vectors of incompatible size.");
            }
            for ( size_t n=0; n<nodes.size(); ++n ) {
                nodes[n].setNodeSlowness( s[n] );
            }
        }

        void getSlowness(std::vector<T1>& slowness) const {
            if (slowness.size() != (ncx+1) * (ncz+1)) {
                slowness.resize((ncx+1) * (ncz+1));
            }
            for ( size_t n=0; n<slowness.size(); ++n ) {
                slowness[n] = nodes[n].getNodeSlowness();
            }
        }

        size_t getNumberOfNodes(const bool primary=false) const {
            if ( primary ) {
                return (ncx+1) * (ncz+1);
            } else {
                return nodes.size();
            }
        }

        size_t getNumberOfCells() const { return ncx*ncz; }

        void getTT(std::vector<T1>& tt, const size_t threadNo=0) const {
            size_t nPrimary = (ncx+1) * (ncz+1);
            tt.resize(nPrimary);
            size_t n = 0;
            for ( size_t nn=0; nn<nodes.size(); ++nn ) {
                if ( nodes[nn].isPrimary() ) {
                    tt[n++] = nodes[nn].getTT(threadNo);
                }
            }
        }

        void saveSlownessXYZ(const char filename[]) const {
            std::ofstream fout( filename );

            for ( size_t n=0; n<nodes.size(); ++n ) {
                fout << nodes[n].getX() << "   " << nodes[n].getZ() << "   " << nodes[n].getNodeSlowness() << '\n';
            }

            fout.close();
        }

        void saveTT(const std::string &, const int, const size_t nt=0,
                    const int format=1) const;

        void saveTTgrad(const std::string &, const size_t nt=0,
                        const bool vtkFormat=0) const;

        const T1 getXmin() const { return xmin; }
        const T1 getZmin() const { return zmin; }
        const T1 getDx() const { return dx; }
        const T1 getDz() const { return dz; }
        const T2 getNcx() const { return ncx; }
        const T2 getNcz() const { return ncz; }

    protected:
        T1 dx;           // cell size in x
        T1 dz;           // cell size in z
        T1 xmin;         // x origin of the grid
        T1 zmin;         // z origin of the grid
        T1 xmax;         // x end of the grid
        T1 zmax;         // z end of the grid
        T2 ncx;          // number of cells in x
        T2 ncz;          // number of cells in x

        mutable std::vector<NODE> nodes;

        T1 computeDt(const NODE& source, const NODE& node) const {
            return (node.getNodeSlowness()+source.getNodeSlowness())/2 * source.getDistance( node );
        }

        T1 computeDt(const NODE& source, const S& node, T1 slo) const {
            return (slo+source.getNodeSlowness())/2 * source.getDistance( node );
        }

        void checkPts(const std::vector<S>&) const;

        bool inPolygon(const S& p, const S poly[], const size_t N) const;

        void grad(S &g, const size_t i, const size_t j, const size_t nt=0) const;

        void grad(S &g, const S &pt, const size_t nt=0) const;

        void getRaypath(const std::vector<S>& Tx,
                        const S &Rx,
                        std::vector<S> &r_data,
                        const size_t threadNo=0) const;
        void getRaypath_old(const std::vector<S>& Tx,
                            const S &Rx,
                            std::vector<S> &r_data,
                            const size_t threadNo=0) const;
        void getRaypath_old2(const std::vector<S>& Tx,
                             const S &Rx,
                             std::vector<S> &r_data,
                             const size_t threadNo=0) const;

        T2 getCellNo(const S& pt) const {
            T1 x = xmax-pt.x < small ? xmax-.5*dx : pt.x;
            T1 z = zmax-pt.z < small ? zmax-.5*dz : pt.z;
            T2 nx = static_cast<T2>( small + (x-xmin)/dx );
            T2 nz = static_cast<T2>( small + (z-zmin)/dz );
            return nx*ncz + nz;
        }

        void getCellIJ(const T2 cellNo, sij<T2>& ind) const {
            ind.i = cellNo / ncz;
            ind.j = cellNo - ncz * ind.i;
        }

        void getIJ(const S& pt, T2& i, T2& j) const {
            i = static_cast<T2>( small + (pt.x-xmin)/dx );
            j = static_cast<T2>( small + (pt.z-zmin)/dz );
        }

        void getIJ(const S& pt, long long& i, long long& j) const {
            i = static_cast<long long>( small + (pt.x-xmin)/dx );
            j = static_cast<long long>( small + (pt.z-zmin)/dz );
        }

        void sweep(const std::vector<bool>& frozen,
                   const size_t threadNo) const;
        void sweep45(const std::vector<bool>& frozen,
                     const size_t threadNo) const;
        void sweep_xz(const std::vector<bool>& frozen,
                      const size_t threadNo) const;
        void sweep_weno3(const std::vector<bool>& frozen,
                         const size_t threadNo) const;
        void sweep_weno3_xz(const std::vector<bool>& frozen,
                            const size_t threadNo) const;

        void update_node(const size_t, const size_t, const size_t=0) const;
        void update_node45(const size_t, const size_t, const size_t=0) const;
        void update_node_xz(const size_t, const size_t, const size_t=0) const;
        void update_node_weno3(const size_t, const size_t, const size_t=0) const;
        void update_node_weno3_xz(const size_t, const size_t, const size_t=0) const;

        void initFSM(const std::vector<S>& Tx,
                     const std::vector<T1>& t0, std::vector<bool>& frozen,
                     const int npts, const size_t threadNo) const;

        T1 getSlowness(const S& Rx) const;

        void dump_secondary(std::ofstream& os) const {
            size_t nPrimary = (ncx+1) * (ncz+1);
            for ( size_t n=nPrimary; n<nodes.size(); ++n ) {
                os << nodes[n].getX() << ' ' << nodes[n].getZ() << '\n';
            }
        }

    private:
        Grid2Drn() {}
        Grid2Drn(const Grid2Drn<T1,T2,S,NODE>& g) {}
        Grid2Drn<T1,T2,S,NODE>& operator=(const Grid2Drn<T1,T2,S,NODE>& g) {}

        T1 getTraveltime(const S& Rx, const size_t threadNo) const;

        T1 getTraveltimeFromRaypath(const std::vector<S>& Tx,
                                    const std::vector<T1>& t0,
                                    const S& Rx,
                                    const size_t threadNo) const final;

        void getRaypath(const std::vector<S>& Tx,
                        const std::vector<T1>& t0,
                        const S& Rx,
                        std::vector<S>& r_data,
                        T1 &tt,
                        const size_t threadNo) const final;

        void getRaypath(const std::vector<S>& Tx,
                        const std::vector<T1>& t0,
                        const S& Rx,
                        std::vector<S>& r_data,
                        std::vector<siv<T1>> &l_data,
                        T1 &tt,
                        const size_t threadNo) const final;
    };

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::checkPts(const std::vector<S>& pts) const {
        for (size_t n=0; n<pts.size(); ++n) {
            if ( pts[n].x < xmin || pts[n].x > xmax ||
                pts[n].z < zmin || pts[n].z > zmax ) {
                std::ostringstream msg;
                msg << "Error: Point (" << pts[n].x << ", "<< pts[n] .z << ") outside grid.";
                throw std::runtime_error(msg.str());
            }
        }
    }



    template<typename T1, typename T2, typename S, typename NODE>
    bool Grid2Drn<T1,T2,S,NODE>::inPolygon(const S& p, const S poly[], const size_t N) const {
        bool c = false;
        for (size_t i = 0, j = N-1; i < N; j = i++) {
            if ((((poly[i].z <= p.z) && (p.z < poly[j].z)) ||
                 ((poly[j].z <= p.z) && (p.z < poly[i].z))) &&
                (p.x < (poly[j].x - poly[i].x) * (p.z - poly[i].z) / (poly[j].z - poly[i].z) + poly[i].x))
                c = !c;
        }
        return c;
    }

    template<typename T1, typename T2, typename S, typename NODE>
    T1 Grid2Drn<T1,T2,S,NODE>::getTraveltime(const S &pt, const size_t nt) const {

        const size_t nnz = ncz+1;

        // bilinear interpolation if not on node

        T1 tt;
        T2 i, j;

        getIJ(pt, i, j);

        if ( std::abs(pt.x - (xmin+i*dx))<small && std::abs(pt.z - (zmin+j*dz))<small ) {
            // on node
            return nodes[i*nnz+j].getTT(nt);
        } else if ( std::abs(pt.x - (xmin+i*dx))<small ) {

            // on edge
            T1 t1 = nodes[i*nnz+j].getTT(nt);
            T1 t2 = nodes[i*nnz+j+1].getTT(nt);

            T1 w1 = (zmin+(j+1)*dz - pt.z)/dz;
            T1 w2 = (pt.z - (zmin+j*dz))/dz;

            tt = t1*w1 + t2*w2;

        } else if ( std::abs(pt.z - (zmin+j*dz))<small ) {

            // on edge
            T1 t1 = nodes[i*nnz+j].getTT(nt);
            T1 t2 = nodes[(i+1)*nnz+j].getTT(nt);

            T1 w1 = (xmin+(i+1)*dx - pt.x)/dx;
            T1 w2 = (pt.x - (xmin+i*dx))/dx;

            tt = t1*w1 + t2*w2;

        } else {

            T1 t1 = nodes[    i*nnz+j  ].getTT(nt);
            T1 t2 = nodes[(i+1)*nnz+j  ].getTT(nt);
            T1 t3 = nodes[    i*nnz+j+1].getTT(nt);
            T1 t4 = nodes[(i+1)*nnz+j+1].getTT(nt);

            T1 w1 = (xmin+(i+1)*dx - pt.x)/dx;
            T1 w2 = (pt.x - (xmin+i*dx))/dx;

            t1 = t1*w1 + t2*w2;
            t2 = t3*w1 + t4*w2;

            w1 = (zmin+(j+1)*dz - pt.z)/dz;
            w2 = (pt.z - (zmin+j*dz))/dz;

            tt = t1*w1 + t2*w2;
        }

        return tt;
    }


    //template<typename T1, typename T2, typename S, typename NODE>
    //T1 Grid2Drn<T1,T2,S,NODE>::getTraveltime(const S& Rx,
    //                                       const std::vector<NODE>& nodes,
    //                                       T2& nodeParentRx, T2& cellParentRx,
    //                                       const size_t threadNo) const {
    //
    //    for ( size_t nn=0; nn<nodes.size(); ++nn ) {
    //        if ( nodes[nn] == Rx ) {
    //            nodeParentRx = nodes[nn].getNodeParent(threadNo);
    //            cellParentRx = nodes[nn].getCellParent(threadNo);
    //            return nodes[nn].getTT(threadNo);
    //        }
    //    }
    //
    //    T1 slownessRx = getSlowness( Rx );
    //
    //    T2 cellNo = getCellNo( Rx );
    //    T2 neibNo = this->neighbors[cellNo][0];
    //    T1 dt = computeDt(nodes[neibNo], Rx, slownessRx);
    //
    //    T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
    //    nodeParentRx = neibNo;
    //    cellParentRx = cellNo;
    //    for ( size_t k=1; k< this->neighbors[cellNo].size(); ++k ) {
    //        neibNo = this->neighbors[cellNo][k];
    //        dt = computeDt(nodes[neibNo], Rx, slownessRx);
    //        if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
    //            traveltime =  nodes[neibNo].getTT(threadNo)+dt;
    //            nodeParentRx = neibNo;
    //        }
    //    }
    //    return traveltime;
    //}

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::saveTT(const std::string& fname, const int all,
                                        const size_t nt,
                                        const int format) const {
        if ( format == 1 ) {
            std::string filename = fname+".dat";
            std::ofstream fout(filename.c_str());
            fout.precision(12);
            for ( T2 n=0; n<nodes.size(); ++n ) {
                if ( nodes[n].isPrimary() || all==1 ) {

                    fout << nodes[n].getX() << '\t'
                    << nodes[n].getZ() << '\t'
                    << nodes[n].getTT(nt) << '\n';
                }
            }
            fout.close();
        } else  if ( format == 2 ) {
#ifdef VTK

            std::string filename = fname+".vtr";
            int nn[3] = {static_cast<int>(ncx+1), 1, static_cast<int>(ncz+1)};

            vtkSmartPointer<vtkDoubleArray> xCoords = vtkSmartPointer<vtkDoubleArray>::New();
            for (size_t n=0; n<nn[0]; ++n) {
                xCoords->InsertNextValue( xmin + n*dx );
            }
            vtkSmartPointer<vtkDoubleArray> yCoords = vtkSmartPointer<vtkDoubleArray>::New();
            yCoords->InsertNextValue( 0.0 );
            vtkSmartPointer<vtkDoubleArray> zCoords = vtkSmartPointer<vtkDoubleArray>::New();
            for (size_t n=0; n<nn[2]; ++n) {
                zCoords->InsertNextValue( zmin + n*dz );
            }

            vtkSmartPointer<vtkRectilinearGrid> rgrid = vtkSmartPointer<vtkRectilinearGrid>::New();
            rgrid->SetDimensions( nn );
            rgrid->SetXCoordinates(xCoords);
            rgrid->SetYCoordinates(yCoords);
            rgrid->SetZCoordinates(zCoords);


            vtkSmartPointer<vtkDoubleArray> newScalars =
            vtkSmartPointer<vtkDoubleArray>::New();

            newScalars->SetName("Travel time");
            newScalars->SetNumberOfComponents(1);
            newScalars->SetNumberOfTuples( rgrid->GetNumberOfPoints() );


            for ( size_t n=0; n<nodes.size(); ++n ) {
                if ( nodes[n].isPrimary() ) {
                    vtkIdType id = rgrid->FindPoint(nodes[n].getX(), 0.0, nodes[n].getZ());
                    newScalars->SetTuple1(id, nodes[n].getTT(nt) );
                }
            }
            rgrid->GetPointData()->SetScalars(newScalars);

            vtkSmartPointer<vtkXMLRectilinearGridWriter> writer =
            vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();

            writer->SetFileName( filename.c_str() );
            writer->SetInputData( rgrid );
            writer->SetDataModeToBinary();
            writer->Update();
#else
            std::cerr << "VTK not included during compilation.\nNothing saved.\n";
#endif
        } else if ( format == 3 ){
            std::string filename = fname+".bin";
            std::ofstream fout(filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
            for ( T2 n=0; n<nodes.size(); ++n ) {
                if ( nodes[n].isPrimary() || all==1 ) {
                    T1 tmp[] = { nodes[n].getX(), nodes[n].getZ(), nodes[n].getTT(nt) };
                    fout.write( (char*)tmp, 3*sizeof(T1) );
                }
            }
            fout.close();
        } else {
            throw std::runtime_error("Unsupported format for saving traveltimes");
        }
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::saveTTgrad(const std::string& fname,
                                            const size_t nt,
                                            const bool vtkFormat) const {

        if (vtkFormat) {
#ifdef VTK

            std::string filename = fname+".vtr";
            int nn[3] = {static_cast<int>(ncx), 1, static_cast<int>(ncz)};

            vtkSmartPointer<vtkDoubleArray> xCoords = vtkSmartPointer<vtkDoubleArray>::New();
            for (size_t n=0; n<nn[0]; ++n) {
                xCoords->InsertNextValue( xmin + (0.5+n)*dx );
            }
            vtkSmartPointer<vtkDoubleArray> yCoords = vtkSmartPointer<vtkDoubleArray>::New();
            yCoords->InsertNextValue( 0.0 );
            vtkSmartPointer<vtkDoubleArray> zCoords = vtkSmartPointer<vtkDoubleArray>::New();
            for (size_t n=0; n<nn[2]; ++n) {
                zCoords->InsertNextValue( zmin + (0.5+n)*dz );
            }

            vtkSmartPointer<vtkRectilinearGrid> rgrid = vtkSmartPointer<vtkRectilinearGrid>::New();
            rgrid->SetDimensions( nn );
            rgrid->SetXCoordinates(xCoords);
            rgrid->SetYCoordinates(yCoords);
            rgrid->SetZCoordinates(zCoords);


            vtkSmartPointer<vtkDoubleArray> grad_tt =
            vtkSmartPointer<vtkDoubleArray>::New();

            grad_tt->SetName("grad tt");
            grad_tt->SetNumberOfComponents(3);
            grad_tt->SetComponentName(0, "x");
            grad_tt->SetComponentName(1, "y");
            grad_tt->SetComponentName(2, "z");
            grad_tt->SetNumberOfTuples( rgrid->GetNumberOfPoints() );


            double x[3];
            x[1] = 0.0;
            for ( size_t i=0; i<ncx; ++i ) {
                for ( size_t j=0; j<ncz; ++j ) {
                    S g;
                    grad(g, i, j, nt);

                    x[0] = xmin + (i+0.5)*dx;
                    x[2] = zmin + (j+0.5)*dz;

                    vtkIdType id = rgrid->FindPoint(x);
                    grad_tt->SetTuple3(id, g.x, 0.0, g.z );
                }
            }
            rgrid->GetPointData()->SetVectors( grad_tt );


            vtkSmartPointer<vtkXMLRectilinearGridWriter> writer =
            vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();

            writer->SetFileName( filename.c_str() );
            writer->SetInputData( rgrid );
            writer->SetDataModeToBinary();
            writer->Update();
#else
            std::cerr << "VTK not included during compilation.\nNothing saved.\n";
#endif
        } else {
            std::string filename = fname+".dat";
            std::ofstream fout(filename.c_str());
            fout.precision(12);
            for ( size_t i=0; i<ncx; ++i ) {
                for ( size_t j=0; j<ncz; ++j ) {
                    S g;
                    grad(g, i, j, nt);

                    T1 x = xmin + (i+0.5)*dx;
                    T1 z = zmin + (j+0.5)*dz;

                    fout << x << ' ' << z << ' ' << g.x << ' ' << g.z << '\n';
                }
            }
            fout.close();
        }
    }


    //template<typename T1, typename T2, typename S, typename NODE>
    //void Grid2Drn<T1,T2,S,NODE>::saveTTgrad2(const std::string& fname,
    //                                      const size_t nt,
    //                                      const bool vtkFormat) const {
    //
    //    if (vtkFormat) {
    //#ifdef VTK
    //
    //        std::string filename = fname+".vtp";
    //
    //
    //
    //        vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    //        vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    //        vtkSmartPointer<vtkDoubleArray> grad_tt =
    //        vtkSmartPointer<vtkDoubleArray>::New();
    //
    //        grad_tt->SetName("grad tt");
    //        grad_tt->SetNumberOfComponents(3);
    //        grad_tt->SetComponentName(0, "x");
    //        grad_tt->SetComponentName(1, "y");
    //        grad_tt->SetComponentName(2, "z");
    //
    //        size_t npts=12;
    //        pts->SetNumberOfPoints(npts);
    //        grad_tt->SetNumberOfTuples( npts );
    //
    //        S g, pt;
    //
    //        pt.x = 50.23;
    //        pt.z = 100.1;
    //        grad(g, pt, nt);
    //        pts->InsertNextPoint(pt.x, 0.0, pt.z);
    //        grad_tt->InsertNextTuple3(g.x, 0.0, g.z);
    //
    //        pt.x = 49.23;
    //        pt.z = 100.1;
    //        grad(g, pt, nt);
    //        pts->InsertNextPoint(pt.x, 0.0, pt.z);
    //        grad_tt->InsertNextTuple3(g.x, 0.0, g.z);
    //
    //        pt.x = 50.23;
    //        pt.z = 99.9;
    //        grad(g, pt, nt);
    //        pts->InsertNextPoint(pt.x, 0.0, pt.z);
    //        grad_tt->InsertNextTuple3(g.x, 0.0, g.z);
    //
    //        pt.x = 50.23;
    //        pt.z = 100.9;
    //        grad(g, pt, nt);
    //        pts->InsertNextPoint(pt.x, 0.0, pt.z);
    //        grad_tt->InsertNextTuple3(g.x, 0.0, g.z);
    //
    //        pt.x = 49.23;
    //        pt.z = 99.9;
    //        grad(g, pt, nt);
    //        pts->InsertNextPoint(pt.x, 0.0, pt.z);
    //        grad_tt->InsertNextTuple3(g.x, 0.0, g.z);
    //
    //        pt.x = 50.7;
    //        pt.z = 99.8;
    //        grad(g, pt, nt);
    //        pts->InsertNextPoint(pt.x, 0.0, pt.z);
    //        grad_tt->InsertNextTuple3(g.x, 0.0, g.z);
    //
    //        pt.x = 50.0;
    //        pt.z = 100.0;
    //        grad(g, pt, nt);
    //        pts->InsertNextPoint(pt.x, 0.0, pt.z);
    //        grad_tt->InsertNextTuple3(g.x, 0.0, g.z);
    //
    //        pt.x = 50.0;
    //        pt.z = 100.2;
    //        grad(g, pt, nt);
    //        pts->InsertNextPoint(pt.x, 0.0, pt.z);
    //        grad_tt->InsertNextTuple3(g.x, 0.0, g.z);
    //
    //        pt.x = 50.3;
    //        pt.z = 100.0;
    //        grad(g, pt, nt);
    //        pts->InsertNextPoint(pt.x, 0.0, pt.z);
    //        grad_tt->InsertNextTuple3(g.x, 0.0, g.z);
    //
    //        pt.x = 50.7;
    //        pt.z = 100.7;
    //        grad(g, pt, nt);
    //        pts->InsertNextPoint(pt.x, 0.0, pt.z);
    //        grad_tt->InsertNextTuple3(g.x, 0.0, g.z);
    //
    //        pt.x = 50.7;
    //        pt.z = 100.7;
    //        grad(g, pt, nt);
    //        pts->InsertNextPoint(pt.x, 0.0, pt.z);
    //        grad_tt->InsertNextTuple3(g.x, 0.0, g.z);
    //
    //        pt.x = 49.7;
    //        pt.z = 98.0;
    //        grad(g, pt, nt);
    //        pts->InsertNextPoint(pt.x, 0.0, pt.z);
    //        grad_tt->InsertNextTuple3(g.x, 0.0, g.z);
    //
    //        pt.x = 49.2;
    //        pt.z = 98.0;
    //        grad(g, pt, nt);
    //        pts->InsertNextPoint(pt.x, 0.0, pt.z);
    //        grad_tt->InsertNextTuple3(g.x, 0.0, g.z);
    //
    //
    //
    //        polydata->SetPoints(pts);
    //        polydata->GetPointData()->SetVectors(grad_tt);
    //
    //        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    //        writer->SetFileName( filename.c_str() );
    //        writer->SetInputData( polydata );
    //        writer->SetDataModeToBinary();
    //        writer->Update();
    //
    //
    //        #else
    //        std::cerr << "VTK not included during compilation.\nNothing saved.\n";
    //#endif
    //    } else {
    //        std::ofstream fout(fname.c_str());
    //        fout.precision(12);
    //        for ( size_t i=0; i<ncx; ++i ) {
    //            for ( size_t j=0; j<ncz; ++j ) {
    //                S g;
    //                grad(g, i, j, nt);
    //
    //                T1 x = xmin + (i+0.5)*dx;
    //                T1 z = zmin + (j+0.5)*dz;
    //
    //                fout << x << ' ' << z << ' ' << g.x << ' ' << g.z << '\n';
    //            }
    //        }
    //        fout.close();
    //    }
    //}


    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::grad(S& g, const size_t i, const size_t j,
                                      const size_t nt) const {

        // compute average gradient for cell (i,j)

        const size_t nnz = ncz+1;

        g.x = 0.5*(( nodes[(i+1)*nnz+j].getTT(nt)+nodes[(i+1)*nnz+j+1].getTT(nt) ) -
                   ( nodes[    i*nnz+j].getTT(nt)+nodes[    i*nnz+j+1].getTT(nt) ))/dx;
        g.z = 0.5*(( nodes[i*nnz+j+1].getTT(nt)+nodes[(i+1)*nnz+j+1].getTT(nt) ) -
                   ( nodes[i*nnz+j  ].getTT(nt)+nodes[(i+1)*nnz+j  ].getTT(nt) ))/dz;
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::grad(S &g, const S &pt,
                                      const size_t nt) const {

        // compute travel time gradient at point pt

        T1 p1 = pt.x - dx/2.0;
        if (p1 < xmin) {
            p1 = xmin;
        }
        T1 p2 = p1 + dx;
        if (p2 > xmax) {
            p2 = xmax;
            p1 = xmax - dx;
        }
        //        if (p1 == xmin) {
        //            // use 3-pt operator
        //            T1 p3 = p2 + dx;
        //            g.x = (-1.5*getTraveltime({p1, pt.z}, nt) +
        //                   2.0*getTraveltime({p2, pt.z}, nt) -
        //                   0.5*getTraveltime({p3, pt.z}, nt)) / dx;
        //        } else if (p2 == xmax) {
        //            T1 p3 = xmax;
        //            p2 = p1;
        //            p1 = p2 - dx;
        //            g.x = (1.5*getTraveltime({p3, pt.z}, nt) -
        //                   2.0*getTraveltime({p2, pt.z}, nt) +
        //                   0.5*getTraveltime({p1, pt.z}, nt)) / dx;
        //        } else {
        g.x = (getTraveltime({p2, pt.z}, nt) - getTraveltime({p1, pt.z}, nt)) / dx;
        //        }

        p1 = pt.z - dz/2.0;
        if (p1 < zmin) {
            p1 = zmin;
        }
        p2 = p1 + dz;
        if (p2 > zmax) {
            p2 = zmax;
            p1 = zmax - dz;
        }
        //        if (p1 == zmin) {
        //            T1 p3 = p2 + dz;
        //            g.z = (-1.5*getTraveltime({pt.x, p1}, nt) +
        //                   2.0*getTraveltime({pt.x, p2}, nt) -
        //                   0.5*getTraveltime({pt.x, p3}, nt)) / dz;
        //        } else if (p2 == zmax) {
        //            T1 p3 = zmax;
        //            p2 = p1;
        //            p1 = p2 - dx;
        //            g.z = (1.5*getTraveltime({pt.x, p3}, nt) -
        //                   2.0*getTraveltime({pt.x, p2}, nt) +
        //                   0.5*getTraveltime({pt.x, p1}, nt)) / dz;
        //        } else {
        g.z = (getTraveltime({pt.x, p2}, nt) - getTraveltime({pt.x, p1}, nt)) / dz;
        //        }
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::getRaypath(const std::vector<S>& Tx,
                                            const S &Rx,
                                            std::vector<S> &r_data,
                                            const size_t threadNo) const {

        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                return;
            }
        }

        S curr_pt( Rx );
        // distance between opposite nodes of a voxel
        const T1 maxDist = sqrt( dx*dx + dz*dz );
        S g;

        bool reachedTx = false;
        while ( reachedTx == false ) {

            grad(g, curr_pt, threadNo);
            g *= -1.0;

            long long i, k;
            getIJ(curr_pt, i, k);

            // planes we will intersect
            T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
            T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

            if ( std::abs(xp-curr_pt.x)<small) {
                xp += dx*boost::math::sign(g.x);
            }
            if ( std::abs(zp-curr_pt.z)<small) {
                zp += dz*boost::math::sign(g.z);
            }

            // dist to planes
            T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
            T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

            if ( tx<tz ) { // closer to xp
                curr_pt += tx*g;
                curr_pt.x = xp;     // make sure we don't accumulate rounding errors
            } else {
                curr_pt += tz*g;
                curr_pt.z = zp;
            }

            if ( curr_pt.x < xmin || curr_pt.x > xmax ||
                curr_pt.z < zmin || curr_pt.z > zmax ) {
                //  we are going oustide the grid!
                std::ostringstream msg;
                msg << "Error while computing raypaths: going outside grid \n\
                Rx: " << Rx << "\n\
                Tx: " << Tx[0] << "\n";
                for ( size_t ns=1; ns<Tx.size(); ++ns ) {
                    msg << "\
                    " << Tx[ns] << "\n";
                }
                throw std::runtime_error(msg.str());
            }

            r_data.push_back( curr_pt );

            // are we close enough to one the Tx nodes ?
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                if ( curr_pt.getDistance( Tx[ns] ) < maxDist ) {
                    r_data.push_back( Tx[ns] );
                    reachedTx = true;
                }
            }

        }
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::getRaypath_old(const std::vector<S>& Tx,
                                                const S &Rx,
                                                std::vector<S> &r_data,
                                                const size_t threadNo) const {

        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                return;
            }
        }


        long long int iIn, kIn, iOut=-1, kOut=-1; // Out for cell we are exiting; In for cell we are entering
        S curr_pt( Rx );
        S gOut = {0.0, 0.0};

        // distance between opposite nodes of a voxel
        const T1 maxDist = sqrt( dx*dx + dz*dz );

        getIJ(curr_pt, iIn, kIn);

        bool reachedTx = false;
        while ( reachedTx == false ) {

            bool onNode=false;

            if ( std::abs(remainder(curr_pt.x,dx))<small &&
                std::abs(remainder(curr_pt.z,dz))<small ) {
                onNode = true;
            }

            if ( onNode ) {

                T2 i, k;
                getIJ(curr_pt, i, k);
                std::vector<sij<T2>> cells;

                // find voxels touching node
                if ( i<=ncx && k<=ncz )
                    cells.push_back( {i,k} );
                if ( i<=ncx && k>0 )
                    cells.push_back( {i,k-1} );
                if ( i>0 && k<=ncz )
                    cells.push_back( {i-1,k} );
                if ( i>0 && k>0 )
                    cells.push_back( {i-1,k-1} );

                gOut = static_cast<T1>(0.0);
                S g;
                T2 nc;
                for ( nc=0; nc<cells.size(); ++nc ) {
                    grad( g, cells[nc].i, cells[nc].j, threadNo );
                    g *= -1.0;
                    gOut += g;
                }
                gOut /= nc;  // gOut holds average grad

                if ((gOut.x<0.0 && i==0) || (gOut.x>0.0 && i==ncx+1) ||
                    (gOut.z<0.0 && k==0) || (gOut.z>0.0 && k==ncz+1)) {
                    //  we are going oustide the grid!
                    std::ostringstream msg;
                    msg << "Error while computing raypaths: going outside grid \n\
                    Rx: " << Rx << "\n\
                    Tx: " << Tx[0] << "\n";
                    for ( size_t ns=1; ns<Tx.size(); ++ns ) {
                        msg << "\
                        " << Tx[ns] << "\n";
                    }
                    throw std::runtime_error(msg.str());
                }

                iOut = boost::math::sign(gOut.x)<0.0 ? i-1 : i;
                kOut = boost::math::sign(gOut.z)<0.0 ? k-1 : k;

                // planes we will intersect
                T1 xp = xmin + dx*(i + boost::math::sign(gOut.x));
                T1 zp = zmin + dz*(k + boost::math::sign(gOut.z));


                // dist to planes
                T1 tx = (xp - curr_pt.x)/gOut.x;
                T1 tz = (zp - curr_pt.z)/gOut.z;

                if ( tx<tz ) { // closer to xp
                    curr_pt += tx*gOut;
                    curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                    iIn = iOut + boost::math::sign(gOut.x);
                    kIn = kOut;
                } else {
                    curr_pt += tz*gOut;
                    curr_pt.z = zp;
                    iIn = iOut;
                    kIn = kOut + boost::math::sign(gOut.z);
                }

            } else {

                S gIn;
                grad( gIn, iIn, kIn, threadNo );
                gIn *= -1.0;

                if ( iIn == iOut && kIn == kOut) {
                    // ray is returning to cell it was exiting
                    // we might be at grazing incidence
                    // check if gIn is significantly different from gOut

                    S diff = normalize(gOut)-normalize(gIn);
                    if ( norm(diff) > small ) {
                        throw std::runtime_error("Error while computing raypaths: raypath not converging!");
                    }
                }

                gOut = gIn;
                iOut = iIn;
                kOut = kIn;

                if ((gOut.x<0.0 && iOut==0) || (gOut.x>0.0 && iOut==ncx+1) ||
                    (gOut.z<0.0 && kOut==0) || (gOut.z>0.0 && kOut==ncz+1)) {
                    //  we are going oustide the grid!
                    std::ostringstream msg;
                    msg << "Error while computing raypaths: going outside grid \n\
                    Rx: " << Rx << "\n\
                    Tx: " << Tx[0] << "\n";
                    for ( size_t ns=1; ns<Tx.size(); ++ns ) {
                        msg << "\
                        " << Tx[ns] << "\n";
                    }
                    throw std::runtime_error(msg.str());
                }

                // planes we will intersect
                T1 xp = xmin + dx*(iIn + (boost::math::sign(gOut.x)>0.0 ? 1.0 : 0.0));
                T1 zp = zmin + dz*(kIn + (boost::math::sign(gOut.z)>0.0 ? 1.0 : 0.0));

                if ( std::abs(xp-curr_pt.x)<small) {
                    xp += dx*boost::math::sign(gOut.x);
                    iOut += boost::math::sign(gOut.x);
                }
                if ( std::abs(zp-curr_pt.z)<small) {
                    zp += dz*boost::math::sign(gOut.z);
                    kOut += boost::math::sign(gOut.z);
                }

                // dist to planes
                T1 tx = (xp - curr_pt.x)/gOut.x;
                T1 tz = (zp - curr_pt.z)/gOut.z;

                if ( tx<tz ) { // closer to xp
                    curr_pt += tx*gOut;
                    curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                    iIn = iOut + boost::math::sign(gOut.x);
                    kIn = kOut;
                } else {
                    curr_pt += tz*gOut;
                    curr_pt.z = zp;
                    iIn = iOut;
                    kIn = kOut + boost::math::sign(gOut.z);
                }

            }


            r_data.push_back( curr_pt );

            // are we close enough to one the of Tx nodes ?
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                if ( curr_pt.getDistance( Tx[ns] ) < maxDist ) {
                    r_data.push_back( Tx[ns] );
                    reachedTx = true;
                }
            }
        }

    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::getRaypath_old2(const std::vector<S>& Tx,
                                                 const S &Rx,
                                                 std::vector<S> &r_data,
                                                 const size_t threadNo) const {

        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                return;
            }
        }

        long long int iIn, jIn, iOut, jOut; // iOut, jOut for cell we are exiting; iIn, jIn for cell we are entering
        S curr_pt( Rx );
        S gOut;

        // distance between opposite nodes of a cell
        const T1 maxDist = sqrt( dx*dx + dz*dz );

        bool onNode=false;
        if ( std::abs(remainder(curr_pt.x,dx))<small && std::abs(remainder(curr_pt.z,dz))<small ) {
            onNode = true;
        }

        if ( !onNode ) {

            // go to first edge
            T2 i, j;
            getIJ(curr_pt, i, j);
            iOut = i;
            jOut = j;
            grad( gOut, iOut, jOut, threadNo );
            gOut *= -1.0;

            T1 theta = atan2( gOut.z, gOut.x );

            if ( theta > pi/2. ) {
                T1 x1 = xmin+iOut*dx;
                T1 ddx = curr_pt.x - x1;
                if ( ddx == 0.0 ) { x1 -= dx; ddx = curr_pt.x - x1; } // we are on the rightmost edge
                T1 ddz = ddx*gOut.z/gOut.x;    // ratio is negative
                T1 z1 = curr_pt.z - ddz;
                T1 d1 = sqrt(ddx*ddx + ddz*ddz);

                T1 z2 = zmin+(jOut+1)*dz;
                ddz = z2 - curr_pt.z;
                if ( ddz == 0.0 ) { z2 += dz; ddz = z2 - curr_pt.z; }
                ddx = ddz*gOut.x/gOut.z;
                T1 x2 = curr_pt.x + ddx;
                T1 d2 = sqrt(ddx*ddx + ddz*ddz);

                if ( d1 <= d2 ) {
                    curr_pt.x = x1;
                    curr_pt.z = z1;
                    iIn = iOut-1;
                    jIn = jOut;
                } else {
                    curr_pt.x = x2;
                    curr_pt.z = z2;
                    iIn = iOut;
                    jIn = jOut+1;
                }
            } else if ( theta > 0. ) {
                T1 x1 = xmin+(iOut+1)*dx;
                T1 ddx = x1 - curr_pt.x;
                if ( ddx == 0.0 ) { x1 += dx; ddx = x1 - curr_pt.x; }
                T1 ddz = ddx*gOut.z/gOut.x;
                T1 z1 = curr_pt.z + ddz;
                T1 d1 = sqrt(ddx*ddx + ddz*ddz);

                T1 z2 = zmin+(jOut+1)*dz;
                ddz = z2 - curr_pt.z;
                if ( ddz == 0.0 ) { z2 += dz; ddz = z2 - curr_pt.z; }
                ddx = ddz*gOut.x/gOut.z;
                T1 x2 = curr_pt.x + ddx;
                T1 d2 = sqrt(ddx*ddx + ddz*ddz);

                if ( d1 <= d2 ) {
                    curr_pt.x = x1;
                    curr_pt.z = z1;
                    iIn = iOut+1;
                    jIn = jOut;
                } else {
                    curr_pt.x = x2;
                    curr_pt.z = z2;
                    iIn = iOut;
                    jIn = jOut+1;
                }
            } else if ( theta > -pi/2. ) {
                T1 x1 = xmin+(iOut+1)*dx;
                T1 ddx = x1 - curr_pt.x;
                if ( ddx == 0.0 ) { x1 += dx; ddx = x1 = curr_pt.x; }
                T1 ddz = ddx*gOut.z/gOut.x;  // ratio is negative
                T1 z1 = curr_pt.z + ddz;
                T1 d1 = sqrt(ddx*ddx + ddz*ddz);

                T1 z2 = zmin+jOut*dz;
                ddz = curr_pt.z - z2;
                if ( ddz == 0.0 ) { z2 -= dz; ddz = curr_pt.z - z2; }
                ddx = ddz*gOut.x/gOut.z;
                T1 x2 = curr_pt.x - ddx;
                T1 d2 = sqrt(ddx*ddx + ddz*ddz);

                if ( d1 <= d2 ) {
                    curr_pt.x = x1;
                    curr_pt.z = z1;
                    iIn = iOut+1;
                    jIn = jOut;
                } else {
                    curr_pt.x = x2;
                    curr_pt.z = z2;
                    iIn = iOut;
                    jIn = jOut-1;
                }

            } else {
                T1 x1 = xmin+iOut*dx;
                T1 ddx = curr_pt.x - x1;
                if ( ddx == 0.0 ) { x1 -= dx; ddx = curr_pt.x - x1; }
                T1 ddz = ddx*gOut.z/gOut.x;
                T1 z1 = curr_pt.z - ddz;
                T1 d1 = sqrt(ddx*ddx + ddz*ddz);

                T1 z2 = zmin+jOut*dz;
                ddz = curr_pt.z - z2;
                if ( ddz == 0.0 ) { z2 -= dz; ddz = curr_pt.z - z2; }
                ddx = ddz*gOut.x/gOut.z;
                T1 x2 = curr_pt.x - ddx;
                T1 d2 = sqrt(ddx*ddx + ddz*ddz);

                if ( d1 <= d2 ) {
                    curr_pt.x = x1;
                    curr_pt.z = z1;
                    iIn = iOut-1;
                    jIn = jOut;
                } else {
                    curr_pt.x = x2;
                    curr_pt.z = z2;
                    iIn = iOut;
                    jIn = jOut-1;
                }
            }

            if ( iIn<0 || iIn>ncx || jIn<0 || jIn>ncz ) {
                //  we are going oustide the grid!
                std::ostringstream msg;
                msg << "Error while computing raypaths: going outside grid \n\
                Rx: " << Rx << "\n\
                Tx: " << Tx[0] << "\n";
                for ( size_t ns=1; ns<Tx.size(); ++ns ) {
                    msg << "\
                    " << Tx[ns] << "\n";
                }
                throw std::runtime_error(msg.str());
            }

            r_data.push_back( curr_pt );

            // are we close enough to one the of Tx nodes ?
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                if ( curr_pt.getDistance( Tx[ns] ) < maxDist ) {
                    r_data.push_back( Tx[ns] );
                    return;
                }
            }

            onNode = false;
            if ( std::abs(remainder(curr_pt.x,dx))<small && std::abs(remainder(curr_pt.z,dz))<small ) {
                onNode = true;
            }

        }

        bool reachedTx = false;
        while ( reachedTx == false ) {

            if ( onNode ) {

                // compute average gradient
                T2 i, j;
                getIJ(curr_pt, i, j);
                std::vector<sij<T2>> cells;

                // find cells touching node
                if ( i<=ncx && j<=ncz )
                    cells.push_back( {i, j} );
                if ( i>0 && j<=ncz )
                    cells.push_back( {i-1, j} );
                if ( i<=ncx && j>0 )
                    cells.push_back( {i, j-1} );
                if ( i>0 && j>0 )
                    cells.push_back( {i-1, j-1} );

                gOut = static_cast<T1>(0.0);
                S g;
                T2 nc;
                for ( nc=0; nc<cells.size(); ++nc ) {
                    grad( g, cells[nc].i, cells[nc].j, threadNo );
                    g *= -1.0;
                    gOut += g;
                }
                gOut /= nc;  // gOut holds average grad


                if ((gOut.x<0.0 && i==0) || (gOut.x>0.0 && i==ncx+1) ||
                    (gOut.z<0.0 && j==0) || (gOut.z>0.0 && j==ncz+1)) {
                    //  we are going oustide the grid!
                    std::ostringstream msg;
                    msg << "Error while computing raypaths: going outside grid \n\
                    Rx: " << Rx << "\n\
                    Tx: " << Tx[0] << "\n";
                    for ( size_t ns=1; ns<Tx.size(); ++ns ) {
                        msg << "\
                        " << Tx[ns] << "\n";
                    }
                    throw std::runtime_error(msg.str());
                }


                T1 theta = atan2( gOut.z, gOut.x );

                if ( gOut.z==0.0 && gOut.x<0.0 ) {          // toward node i-1, j
                    curr_pt.x -= dx;
                    r_data.push_back( curr_pt );
                    onNode = true;

                } else if ( theta>pi/2.0 ) {                // toward cell i-1, j
                    iOut = i-1;
                    jOut = j;
                    if ( gOut.z/gOut.x > -dz/dx ) {  // ratio is negative
                        curr_pt.x -= dx;
                        curr_pt.z += dx*sin(pi-theta);
                        r_data.push_back( curr_pt );
                        iIn = iOut - 1;
                        jIn = jOut;
                        onNode = false;
                    } else if ( gOut.z/gOut.x < -dz/dx ) {
                        curr_pt.x -= dz*sin(theta-pi/2.0);
                        curr_pt.z += dz;
                        iIn = iOut;
                        jIn = jOut + 1;
                        onNode = false;
                    } else { // gOut.z/gOut.x == dz/dx   ->    toward node i-1, j+1
                        curr_pt.x -= dx;
                        curr_pt.z += dz;
                        r_data.push_back( curr_pt );
                        onNode = true;
                    }

                } else if ( gOut.z>0.0 && gOut.x==0.0 ) {   // toward node i, j+1
                    curr_pt.z += dz;
                    r_data.push_back( curr_pt );
                    onNode = true;

                } else if ( theta>0.0 ) {                   // toward cell i, j
                    iOut = i;
                    jOut = j;
                    if ( gOut.z/gOut.x < dz/dx ) {  // ratio is positive
                        curr_pt.x += dx;
                        curr_pt.z += dx*sin(theta);
                        r_data.push_back( curr_pt );
                        iIn = iOut + 1;
                        jIn = jOut;
                        onNode = false;
                    } else if ( gOut.z/gOut.x > dz/dx ) {
                        curr_pt.x += dz*sin(pi/2.-theta);
                        curr_pt.z += dz;
                        r_data.push_back( curr_pt );
                        iIn = iOut;
                        jIn = jOut + 1;
                        onNode = false;
                    } else { // gOut.z/gOut.x == dz/dx   ->    toward node i+1, j+1
                        curr_pt.x += dx;
                        curr_pt.z += dz;
                        r_data.push_back( curr_pt );
                        onNode = true;
                    }

                } else if ( gOut.z==0.0 && gOut.x>0.0 ) {   // toward node i+1, j
                    curr_pt.x += dx;
                    r_data.push_back( curr_pt );
                    onNode = true;

                } else if ( theta>-pi/2.0 ) {               // toward cell i, j-1
                    iOut = i;
                    jOut = j-1;
                    if ( gOut.z/gOut.x > -dz/dx ) {  // ratio is negative
                        curr_pt.x += dx;
                        curr_pt.z -= dx*sin(-theta);
                        r_data.push_back( curr_pt );
                        iIn = iOut + 1;
                        jIn = jOut;
                        onNode = false;
                    } else if ( gOut.z/gOut.x < -dz/dx ) {
                        curr_pt.x += dz*sin(pi/2.+theta);
                        curr_pt.z -= dz;
                        r_data.push_back( curr_pt );
                        iIn = iOut;
                        jIn = jOut - 1;
                        onNode = false;
                    } else { // gOut.z/gOut.x == dz/dx   ->    toward node i+1, j-1
                        curr_pt.x += dx;
                        curr_pt.z -= dz;
                        r_data.push_back( curr_pt );
                        onNode = true;
                    }

                } else if ( gOut.z<0.0 && gOut.x==0.0 ) {   // toward node i, j-1
                    curr_pt.z -= dz;
                    r_data.push_back( curr_pt );
                    onNode = true;

                } else {                                    // toward cell i-1, j-1
                    iOut = i-1;
                    jOut = j-1;
                    if ( gOut.z/gOut.x < dz/dx ) {  // ratio is positive
                        curr_pt.x -= dx;
                        curr_pt.z -= dx*sin(pi+theta);
                        r_data.push_back( curr_pt );
                        iIn = iOut - 1;
                        jIn = jOut;
                        onNode = false;
                    } else if ( gOut.z/gOut.x > dz/dx ) {
                        curr_pt.x -= dz*sin(-pi/2.0-theta);
                        curr_pt.z -= dz;
                        r_data.push_back( curr_pt );
                        iIn = iOut;
                        jIn = jOut - 1;
                        onNode = false;
                    } else { // gOut.z/gOut.x == dz/dx   ->    toward node i-1, j-1
                        curr_pt.x -= dx;
                        curr_pt.z -= dz;
                        r_data.push_back( curr_pt );
                        onNode = true;
                    }
                }

            } else { // not on node, must be on an edge

                S gIn;
                grad( gIn, iIn, jIn, threadNo );
                gIn *= -1.0;
                T1 theta = atan2( gIn.z, gIn.x );

                if ( iIn == iOut && jIn == jOut ) {
                    // ray is returning to cell it was exiting
                    // we might be at grazing incidence
                    // check if gIn is significantly different from gOut

                    T1 thetaOut = atan2( gOut.z, gOut.x );
                    if ( std::abs(theta-thetaOut) < pi/180. ) {
                        // less that 1 degree difference, OK
                        gOut = gIn;
                    } else {
                        throw std::runtime_error("Error while computing raypaths: raypath not converging!");
                    }
                }

                if ( theta > pi/2.0 ) {
                    T1 x1 = xmin+iIn*dx;
                    T1 ddx = curr_pt.x - x1;
                    if ( ddx == 0 ) { x1 -= dx; ddx = curr_pt.x - x1; iIn--; }
                    T1 ddz = ddx*gIn.z/gIn.x;    // ratio is negative
                    T1 z1 = curr_pt.z - ddz;
                    T1 d1 = sqrt(ddx*ddx + ddz*ddz);

                    T1 z2 = zmin+(jIn+1)*dz;
                    ddz = z2 - curr_pt.z;
                    if ( ddz == 0.0 ) { z2 += dz; ddz = z2 - curr_pt.z; jIn++; }
                    ddx = ddz*gIn.x/gIn.z;
                    T1 x2 = curr_pt.x + ddx;
                    T1 d2 = sqrt(ddx*ddx + ddz*ddz);

                    if ( d1 <= d2 ) {
                        curr_pt.x = x1;
                        curr_pt.z = z1;
                        iOut = iIn;
                        jOut = jIn;
                        iIn--;
                    } else {
                        curr_pt.x = x2;
                        curr_pt.z = z2;
                        iOut = iIn;
                        jOut = jIn;
                        jIn++;
                    }

                } else if ( theta > 0. ) {
                    T1 x1 = xmin+(iIn+1)*dx;
                    T1 ddx = x1 - curr_pt.x;
                    if ( ddx == 0.0 ) { x1 += dx; ddx = x1 - curr_pt.x; iIn++; }
                    T1 ddz = ddx*gIn.z/gIn.x;
                    T1 z1 = curr_pt.z + ddz;
                    T1 d1 = sqrt(ddx*ddx + ddz*ddz);

                    T1 z2 = zmin+(jIn+1)*dz;
                    ddz = z2 - curr_pt.z;
                    if ( ddz == 0.0 ) { z2 += dz; ddz = z2 - curr_pt.z; jIn++; }
                    ddx = ddz*gIn.x/gIn.z;
                    T1 x2 = curr_pt.x + ddx;
                    T1 d2 = sqrt(ddx*ddx + ddz*ddz);

                    if ( d1 <= d2 ) {
                        curr_pt.x = x1;
                        curr_pt.z = z1;
                        iOut = iIn;
                        jOut = jIn;
                        iIn++;
                    } else {
                        curr_pt.x = x2;
                        curr_pt.z = z2;
                        iOut = iIn;
                        jOut = jIn;
                        jIn++;
                    }

                } else if ( theta > -pi/2. ) {
                    T1 x1 = xmin+(iIn+1)*dx;
                    T1 ddx = x1 - curr_pt.x;
                    if ( ddx == 0.0 ) { x1 += dx; ddx = x1 - curr_pt.x; iIn++; }
                    T1 ddz = ddx*gIn.z/gIn.x;  // ratio is negative
                    T1 z1 = curr_pt.z + ddz;
                    T1 d1 = sqrt(ddx*ddx + ddz*ddz);

                    T1 z2 = zmin+jIn*dz;
                    ddz = curr_pt.z - z2;
                    if ( ddz == 0.0 ) { z2 -= dz; ddz = curr_pt.z - z2; jIn--; }
                    ddx = ddz*gIn.x/gIn.z;
                    T1 x2 = curr_pt.x - ddx;
                    T1 d2 = sqrt(ddx*ddx + ddz*ddz);

                    if ( d1 <= d2 ) {
                        curr_pt.x = x1;
                        curr_pt.z = z1;
                        iOut = iIn;
                        jOut = jIn;
                        iIn++;
                    } else {
                        curr_pt.x = x2;
                        curr_pt.z = z2;
                        iOut = iIn;
                        jOut = jIn;
                        jIn--;
                    }

                } else {
                    T1 x1 = xmin+iIn*dx;
                    T1 ddx = curr_pt.x - x1;
                    if ( ddx == 0.0 ) { x1 -= dx; ddx = curr_pt.x - x1; iIn--; }
                    T1 ddz = ddx*gIn.z/gIn.x;
                    T1 z1 = curr_pt.z - ddz;
                    T1 d1 = sqrt(ddx*ddx + ddz*ddz);

                    T1 z2 = zmin+jIn*dz;
                    ddz = curr_pt.z - z2;
                    if ( ddz == 0.0 ) { z2 -= dz; ddz = curr_pt.z - z2; jIn--; }
                    ddx = ddz*gIn.x/gIn.z;
                    T1 x2 = curr_pt.x - ddx;
                    T1 d2 = sqrt(ddx*ddx + ddz*ddz);

                    if ( d1 <= d2 ) {
                        curr_pt.x = x1;
                        curr_pt.z = z1;
                        iOut = iIn;
                        jOut = jIn;
                        iIn--;
                    } else {
                        curr_pt.x = x2;
                        curr_pt.z = z2;
                        iOut = iIn;
                        jOut = jIn;
                        jIn--;
                    }
                }

                onNode = false;
                if ( std::abs(remainder(curr_pt.x,dx))<small && std::abs(remainder(curr_pt.z,dz))<small ) {
                    onNode = true;
                }

                gOut = gIn;

            }

            r_data.push_back( curr_pt );

            //        std::cout << curr_pt.x << '\t' << curr_pt.z << '\t' << iIn << '\t' << jIn << '\t' << iOut << '\t' << jOut << '\n';

            // are we close enough to one the of Tx nodes ?
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                if ( curr_pt.getDistance( Tx[ns] ) < maxDist ) {
                    r_data.push_back( Tx[ns] );
                    reachedTx = true;
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::sweep(const std::vector<bool>& frozen,
                                       const size_t threadNo) const {

        //    std::cout << '\n';
        //    for ( int j=ncz; j>=0; --j ) {
        //        for ( size_t i=0; i<=ncx; ++i ) {
        //            std::cout << nodes[i*(ncz+1)+j].getTT(threadNo) << ' ';
        //        }
        //        std::cout << '\n';
        //    }

        // sweep first direction
        for ( size_t i=0; i<=ncx; ++i ) {
            for ( size_t j=0; j<=ncz; ++j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node(i, j, threadNo);
                }
            }
        }
        //    std::cout << '\n';
        //    for ( int j=ncz; j>=0; --j ) {
        //        for ( size_t i=0; i<=ncx; ++i ) {
        //            std::cout << nodes[i*(ncz+1)+j].getTT(threadNo) << ' ';
        //        }
        //        std::cout << '\n';
        //    }

        // sweep second direction
        for ( long int i=ncx; i>=0; --i ) {
            for ( size_t j=0; j<=ncz; ++j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node(i, j, threadNo);
                }
            }
        }
        //    std::cout << '\n';
        //    for ( int j=ncz; j>=0; --j ) {
        //        for ( size_t i=0; i<=ncx; ++i ) {
        //            std::cout << nodes[i*(ncz+1)+j].getTT(threadNo) << ' ';
        //        }
        //        std::cout << '\n';
        //    }

        // sweep third direction
        for ( long int i=ncx; i>=0; --i ) {
            for ( long int j=ncz; j>=0; --j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node(i, j, threadNo);
                }
            }
        }
        //    std::cout << '\n';
        //    for ( int j=ncz; j>=0; --j ) {
        //        for ( size_t i=0; i<=ncx; ++i ) {
        //            std::cout << nodes[i*(ncz+1)+j].getTT(threadNo) << ' ';
        //        }
        //        std::cout << '\n';
        //    }

        // sweep fourth direction
        for ( size_t i=0; i<=ncx; ++i ) {
            for ( long int j=ncz; j>=0; --j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node(i, j, threadNo);
                }
            }
        }
        //    std::cout << '\n';
        //    for ( int j=ncz; j>=0; --j ) {
        //        for ( size_t i=0; i<=ncx; ++i ) {
        //            std::cout << nodes[i*(ncz+1)+j].getTT(threadNo) << ' ';
        //        }
        //        std::cout << '\n';
        //    }
    }


    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::sweep45(const std::vector<bool>& frozen,
                                         const size_t threadNo) const {

        // sweep first direction
        for ( size_t i=0; i<=ncx; ++i ) {
            for ( size_t j=0; j<=ncz; ++j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node45(i, j, threadNo);
                }
            }
        }

        // sweep second direction
        for ( long int i=ncx; i>=0; --i ) {
            for ( size_t j=0; j<=ncz; ++j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node45(i, j, threadNo);
                }
            }
        }

        // sweep third direction
        for ( long int i=ncx; i>=0; --i ) {
            for ( long int j=ncz; j>=0; --j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node45(i, j, threadNo);
                }
            }
        }

        // sweep fourth direction
        for ( size_t i=0; i<=ncx; ++i ) {
            for ( long int j=ncz; j>=0; --j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node45(i, j, threadNo);
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::sweep_xz(const std::vector<bool>& frozen,
                                          const size_t threadNo) const {

        // sweep first direction
        for ( size_t i=0; i<=ncx; ++i ) {
            for ( size_t j=0; j<=ncz; ++j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_xz(i, j, threadNo);
                }
            }
        }

        // sweep second direction
        for ( long int i=ncx; i>=0; --i ) {
            for ( size_t j=0; j<=ncz; ++j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_xz(i, j, threadNo);
                }
            }
        }

        // sweep third direction
        for ( long int i=ncx; i>=0; --i ) {
            for ( long int j=ncz; j>=0; --j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_xz(i, j, threadNo);
                }
            }
        }

        // sweep fourth direction
        for ( size_t i=0; i<=ncx; ++i ) {
            for ( long int j=ncz; j>=0; --j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_xz(i, j, threadNo);
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::sweep_weno3(const std::vector<bool>& frozen,
                                             const size_t threadNo) const {

        // sweep first direction
        for ( size_t i=0; i<=ncx; ++i ) {
            for ( size_t j=0; j<=ncz; ++j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_weno3(i, j, threadNo);
                }
            }
        }

        // sweep second direction
        for ( long int i=ncx; i>=0; --i ) {
            for ( size_t j=0; j<=ncz; ++j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_weno3(i, j, threadNo);
                }
            }
        }

        // sweep third direction
        for ( long int i=ncx; i>=0; --i ) {
            for ( long int j=ncz; j>=0; --j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_weno3(i, j, threadNo);
                }
            }
        }

        // sweep fourth direction
        for ( size_t i=0; i<=ncx; ++i ) {
            for ( long int j=ncz; j>=0; --j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_weno3(i, j, threadNo);
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::sweep_weno3_xz(const std::vector<bool>& frozen,
                                                const size_t threadNo) const {

        // sweep first direction
        for ( size_t i=0; i<=ncx; ++i ) {
            for ( size_t j=0; j<=ncz; ++j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_weno3_xz(i, j, threadNo);
                }
            }
        }

        // sweep second direction
        for ( long int i=ncx; i>=0; --i ) {
            for ( size_t j=0; j<=ncz; ++j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_weno3_xz(i, j, threadNo);
                }
            }
        }

        // sweep third direction
        for ( long int i=ncx; i>=0; --i ) {
            for ( long int j=ncz; j>=0; --j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_weno3_xz(i, j, threadNo);
                }
            }
        }

        // sweep fourth direction
        for ( size_t i=0; i<=ncx; ++i ) {
            for ( long int j=ncz; j>=0; --j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_weno3_xz(i, j, threadNo);
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::update_node(const size_t i, const size_t j,
                                             const size_t threadNo) const {

        T1 a, b, t;
        if (i==0)
            a = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo);
        else if (i==ncx)
            a = nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
        else {
            a = nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
            t = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo);
            a = a<t ? a : t;
        }

        if (j==0)
            b = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo);
        else if (j==ncz)
            b = nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
        else {
            b = nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
            t = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo);
            b = b<t ? b : t;
        }

        T1 fh = nodes[i*(ncz+1)+j].getNodeSlowness() * dx;

        if ( std::abs(a-b) >= fh )
            t = (a<b ? a : b) + fh;
        else
            t = 0.5*( a+b + sqrt(2.*fh*fh - (a-b)*(a-b) ) );

        if ( t<nodes[i*(ncz+1)+j].getTT(threadNo) )
            nodes[i*(ncz+1)+j].setTT(t,threadNo);

    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::update_node45(const size_t i, const size_t j,
                                               const size_t threadNo) const {
        // stencil rotated pi/4

        T1 a, b, t;
        if (i==0) {
            if (j!=ncz)
                a = nodes[ (i+1)*(ncz+1)+j+1 ].getTT(threadNo);
            else
                a = std::numeric_limits<T1>::max();
        } else if (i==ncx) {
            if (j!=0)
                a = nodes[ (i-1)*(ncz+1)+j-1 ].getTT(threadNo);
            else
                a = std::numeric_limits<T1>::max();
        } else {
            if (j!=ncz)
                a = nodes[ (i+1)*(ncz+1)+j+1 ].getTT(threadNo);
            else
                a = std::numeric_limits<T1>::max();
            if (j!=0)
                t = nodes[ (i-1)*(ncz+1)+j-1 ].getTT(threadNo);
            else
                t = std::numeric_limits<T1>::max();
            a = a<t ? a : t;
        }

        if (i==0) {
            if (j!=0)
                b = nodes[ (i+1)*(ncz+1)+j-1 ].getTT(threadNo);
            else
                b = std::numeric_limits<T1>::max();
        } else if (i==ncx) {
            if (j!=ncz)
                b = nodes[ (i-1)*(ncz+1)+j+1 ].getTT(threadNo);
            else
                b = std::numeric_limits<T1>::max();
        } else {
            if (j!=0)
                b = nodes[ (i+1)*(ncz+1)+j-1 ].getTT(threadNo);
            else
                b = std::numeric_limits<T1>::max();
            if (j!=ncz)
                t = nodes[ (i-1)*(ncz+1)+j+1 ].getTT(threadNo);
            else
                t = std::numeric_limits<T1>::max();
            b = b<t ? b : t;
        }

        T1 fh = 1.414213562373095 * nodes[i*(ncz+1)+j].getNodeSlowness() *
        dx;
        if ( std::abs(a-b) >= fh )
            t = (a<b ? a : b) + fh;
        else
            t = 0.5*( a+b + sqrt(2.*fh*fh - (a-b)*(a-b) ) );

        if ( t<nodes[i*(ncz+1)+j].getTT(threadNo) )
            nodes[i*(ncz+1)+j].setTT(t,threadNo);
    }


    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::update_node_xz(const size_t i, const size_t j,
                                                const size_t threadNo) const {

        T1 a, b, t;
        if (i==0)
            a = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo);
        else if (i==ncx)
            a = nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
        else {
            a = nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
            t = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo);
            a = a<t ? a : t;
        }

        if (j==0)
            b = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo);
        else if (j==ncz)
            b = nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
        else {
            b = nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
            t = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo);
            b = b<t ? b : t;
        }

        if ( a<b && ((b-a)/dx)>nodes[i*(ncz+1)+j].getNodeSlowness() ) {
            t = a + nodes[i*(ncz+1)+j].getNodeSlowness()*dx;
        } else if ( a>b && ((a-b)/dz)>nodes[i*(ncz+1)+j].getNodeSlowness() ) {
            t = b + nodes[i*(ncz+1)+j].getNodeSlowness()*dz;
        } else {
            T1 dx2 = dx*dx;
            T1 dz2 = dz*dz;
            T1 s2 = nodes[i*(ncz+1)+j].getNodeSlowness()*nodes[i*(ncz+1)+j].getNodeSlowness();
            t = (b*dx2 + a*dz2)/(dx2 + dz2) + sqrt((2.0*a*b*dx2*dz2 - a*a*dx2*dz2 -
                                                    b*b*dx2*dz2 + dx2*dx2*dz2*s2 +
                                                    dx2*dz2*dz2*s2)/((dx2 + dz2)*(dx2 + dz2)));
        }

        if ( t<nodes[i*(ncz+1)+j].getTT(threadNo) )
            nodes[i*(ncz+1)+j].setTT(t,threadNo);
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::update_node_weno3(const size_t i, const size_t j,
                                                   const size_t threadNo) const {

        // not valid if dx != dz

        //    @Article{zhang06,
        //        Title                    = {High Order Fast Sweeping Methods for Static {H}amilton–{J}acobi Equations},
        //        Author                   = {Yong-Tao Zhang and Hong-Kai Zhao and Jianliang Qian},
        //        Journal                  = {Journal of Scientific Computing},
        //        Year                     = {2006},
        //        Number                   = {1},
        //        Pages                    = {25--56},
        //        Volume                   = {29},
        //        DOI                      = {10.1007/s10915-005-9014-3},
        //        URL                      = {http://dx.doi.org/10.1007/s10915-005-9014-3}
        //        }


        T1 a, b, t;
        if (i==0) {
            a = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo);  // fist order
        } else if (i==1) {
            T1 num = nodes[ (i+2)*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) + nodes[ i*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 ap = (1.-w)*(nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo)-nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo))/(2.*dx) +
            w*(-nodes[ (i+2)*(ncz+1)+j ].getTT(threadNo) +4.*nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) -3.*nodes[ i*(ncz+1)+j ].getTT(threadNo))/(2.*dx);

            a = nodes[ i*(ncz+1)+j ].getTT(threadNo) + dx*ap;

            t = nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo); // first order for left
            a = a<t ? a : t;

        } else if (i==ncx) {
            a = nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
        } else if (i==ncx-1) {
            T1 num = nodes[ i*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-2)*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 am = (1.-w)*(nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo)-nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo))/(2.*dx) +
            w*(3.*nodes[ i*(ncz+1)+j ].getTT(threadNo) -4.*nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-2)*(ncz+1)+j ].getTT(threadNo))/(2.*dx);

            a = nodes[ i*(ncz+1)+j ].getTT(threadNo) - dx*am;

            t = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo); // first order for right
            a = a<t ? a : t;

        } else {
            T1 num = nodes[ (i+2)*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) + nodes[ i*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 ap = (1.-w)*(nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo)-nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo))/(2.*dx) +
            w*(-nodes[ (i+2)*(ncz+1)+j ].getTT(threadNo) +4.*nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) -3.*nodes[ i*(ncz+1)+j ].getTT(threadNo))/(2.*dx);

            num = nodes[ i*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-2)*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            w = 1./(1.+2.*r*r);

            T1 am = (1.-w)*(nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo)-nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo))/(2.*dx) +
            w*(3.*nodes[ i*(ncz+1)+j ].getTT(threadNo) -4.*nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-2)*(ncz+1)+j ].getTT(threadNo))/(2.*dx);

            a = nodes[ i*(ncz+1)+j ].getTT(threadNo) - dx*am < nodes[ i*(ncz+1)+j ].getTT(threadNo) + dx*ap ?
            nodes[ i*(ncz+1)+j ].getTT(threadNo) - dx*am : nodes[ i*(ncz+1)+j ].getTT(threadNo) + dx*ap;

        }

        if (j==0) {
            b = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo);
        } else if (j==1) {
            T1 num = nodes[ i*(ncz+1)+j+2 ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 bp = (1.-w)*(nodes[ i*(ncz+1)+j+1 ].getTT(threadNo)-nodes[ i*(ncz+1)+j-1 ].getTT(threadNo))/(2.*dx) +
            w*(-nodes[ i*(ncz+1)+j+2 ].getTT(threadNo) +4.*nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) -3.*nodes[ i*(ncz+1)+j ].getTT(threadNo))/(2.*dx);

            b = nodes[ i*(ncz+1)+j ].getTT(threadNo) + dx*bp;

            t = nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
            b = b<t ? b : t;

        } else if (j==ncz) {
            b = nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
        } else if (j==ncz-1) {
            T1 num = nodes[ i*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j-1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j-2 ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 bm = (1.-w)*(nodes[ i*(ncz+1)+j+1 ].getTT(threadNo)-nodes[ i*(ncz+1)+j-1 ].getTT(threadNo))/(2.*dx) +
            w*(3.*nodes[ i*(ncz+1)+j ].getTT(threadNo) -4.*nodes[ i*(ncz+1)+j-1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j-2 ].getTT(threadNo))/(2.*dx);

            b = nodes[ i*(ncz+1)+j ].getTT(threadNo) - dx*bm;

            t = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo); // first order for right
            b = b<t ? b : t;

        } else {
            T1 num = nodes[ i*(ncz+1)+j+2 ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 bp = (1.-w)*(nodes[ i*(ncz+1)+j+1 ].getTT(threadNo)-nodes[ i*(ncz+1)+j-1 ].getTT(threadNo))/(2.*dx) +
            w*(-nodes[ i*(ncz+1)+j+2 ].getTT(threadNo) +4.*nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) -3.*nodes[ i*(ncz+1)+j ].getTT(threadNo))/(2.*dx);

            num = nodes[ i*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j-1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j-2 ].getTT(threadNo);
            num *= num;
            r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            w = 1./(1.+2.*r*r);

            T1 bm = (1.-w)*(nodes[ i*(ncz+1)+j+1 ].getTT(threadNo)-nodes[ i*(ncz+1)+j-1 ].getTT(threadNo))/(2.*dx) +
            w*(3.*nodes[ i*(ncz+1)+j ].getTT(threadNo) -4.*nodes[ i*(ncz+1)+j-1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j-2 ].getTT(threadNo))/(2.*dx);

            b = nodes[ i*(ncz+1)+j ].getTT(threadNo) - dx*bm < nodes[ i*(ncz+1)+j ].getTT(threadNo) + dx*bp ?
            nodes[ i*(ncz+1)+j ].getTT(threadNo) - dx*bm : nodes[ i*(ncz+1)+j ].getTT(threadNo) + dx*bp;

        }

        T1 fh = nodes[i*(ncz+1)+j].getNodeSlowness() * dx;

        if ( std::abs(a-b) >= fh )
            t = (a<b ? a : b) + fh;
        else
            t = 0.5*( a+b + sqrt(2.*fh*fh - (a-b)*(a-b) ) );

        if ( t<nodes[i*(ncz+1)+j].getTT(threadNo) )
            nodes[i*(ncz+1)+j].setTT(t,threadNo);

        if ( t<nodes[i*(ncz+1)+j].getTT(threadNo) )
            nodes[i*(ncz+1)+j].setTT(t,threadNo);

    }


    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::update_node_weno3_xz(const size_t i, const size_t j,
                                                      const size_t threadNo) const {

        T1 a, b, t;
        if (i==0) {
            a = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo);  // fist order
        } else if (i==1) {
            T1 num = nodes[ (i+2)*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) + nodes[ i*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 ap = (1.-w)*(nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo)-nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo))/(2.*dx) +
            w*(-nodes[ (i+2)*(ncz+1)+j ].getTT(threadNo) +4.*nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) -3.*nodes[ i*(ncz+1)+j ].getTT(threadNo))/(2.*dx);

            a = nodes[ i*(ncz+1)+j ].getTT(threadNo) + dx*ap;

            t = nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo); // first order for left
            a = a<t ? a : t;

        } else if (i==ncx) {
            a = nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
        } else if (i==ncx-1) {
            T1 num = nodes[ i*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-2)*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 am = (1.-w)*(nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo)-nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo))/(2.*dx) +
            w*(3.*nodes[ i*(ncz+1)+j ].getTT(threadNo) -4.*nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-2)*(ncz+1)+j ].getTT(threadNo))/(2.*dx);

            a = nodes[ i*(ncz+1)+j ].getTT(threadNo) - dx*am;

            t = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo); // first order for right
            a = a<t ? a : t;

        } else {
            T1 num = nodes[ (i+2)*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) + nodes[ i*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 ap = (1.-w)*(nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo)-nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo))/(2.*dx) +
            w*(-nodes[ (i+2)*(ncz+1)+j ].getTT(threadNo) +4.*nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) -3.*nodes[ i*(ncz+1)+j ].getTT(threadNo))/(2.*dx);

            num = nodes[ i*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-2)*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            w = 1./(1.+2.*r*r);

            T1 am = (1.-w)*(nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo)-nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo))/(2.*dx) +
            w*(3.*nodes[ i*(ncz+1)+j ].getTT(threadNo) -4.*nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-2)*(ncz+1)+j ].getTT(threadNo))/(2.*dx);

            a = nodes[ i*(ncz+1)+j ].getTT(threadNo) - dx*am < nodes[ i*(ncz+1)+j ].getTT(threadNo) + dx*ap ?
            nodes[ i*(ncz+1)+j ].getTT(threadNo) - dx*am : nodes[ i*(ncz+1)+j ].getTT(threadNo) + dx*ap;

        }

        if (j==0) {
            b = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo);
        } else if (j==1) {
            T1 num = nodes[ i*(ncz+1)+j+2 ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 bp = (1.-w)*(nodes[ i*(ncz+1)+j+1 ].getTT(threadNo)-nodes[ i*(ncz+1)+j-1 ].getTT(threadNo))/(2.*dz) +
            w*(-nodes[ i*(ncz+1)+j+2 ].getTT(threadNo) +4.*nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) -3.*nodes[ i*(ncz+1)+j ].getTT(threadNo))/(2.*dz);

            b = nodes[ i*(ncz+1)+j ].getTT(threadNo) + dz*bp;

            t = nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
            b = b<t ? b : t;

        } else if (j==ncz) {
            b = nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
        } else if (j==ncz-1) {
            T1 num = nodes[ i*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j-1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j-2 ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 bm = (1.-w)*(nodes[ i*(ncz+1)+j+1 ].getTT(threadNo)-nodes[ i*(ncz+1)+j-1 ].getTT(threadNo))/(2.*dz) +
            w*(3.*nodes[ i*(ncz+1)+j ].getTT(threadNo) -4.*nodes[ i*(ncz+1)+j-1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j-2 ].getTT(threadNo))/(2.*dz);

            b = nodes[ i*(ncz+1)+j ].getTT(threadNo) - dz*bm;

            t = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo); // first order for right
            b = b<t ? b : t;

        } else {
            T1 num = nodes[ i*(ncz+1)+j+2 ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 bp = (1.-w)*(nodes[ i*(ncz+1)+j+1 ].getTT(threadNo)-nodes[ i*(ncz+1)+j-1 ].getTT(threadNo))/(2.*dz) +
            w*(-nodes[ i*(ncz+1)+j+2 ].getTT(threadNo) +4.*nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) -3.*nodes[ i*(ncz+1)+j ].getTT(threadNo))/(2.*dz);

            num = nodes[ i*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j-1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j-2 ].getTT(threadNo);
            num *= num;
            r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            w = 1./(1.+2.*r*r);

            T1 bm = (1.-w)*(nodes[ i*(ncz+1)+j+1 ].getTT(threadNo)-nodes[ i*(ncz+1)+j-1 ].getTT(threadNo))/(2.*dz) +
            w*(3.*nodes[ i*(ncz+1)+j ].getTT(threadNo) -4.*nodes[ i*(ncz+1)+j-1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j-2 ].getTT(threadNo))/(2.*dz);

            b = nodes[ i*(ncz+1)+j ].getTT(threadNo) - dz*bm < nodes[ i*(ncz+1)+j ].getTT(threadNo) + dz*bp ?
            nodes[ i*(ncz+1)+j ].getTT(threadNo) - dz*bm : nodes[ i*(ncz+1)+j ].getTT(threadNo) + dz*bp;

        }

        if ( a<b && ((b-a)/dx)>nodes[i*(ncz+1)+j].getNodeSlowness() ) {
            t = a + nodes[i*(ncz+1)+j].getNodeSlowness()*dx;
        } else if ( a>b && ((a-b)/dz)>nodes[i*(ncz+1)+j].getNodeSlowness() ) {
            t = b + nodes[i*(ncz+1)+j].getNodeSlowness()*dz;
        } else {
            T1 dx2 = dx*dx;
            T1 dz2 = dz*dz;
            T1 s2 = nodes[i*(ncz+1)+j].getNodeSlowness()*nodes[i*(ncz+1)+j].getNodeSlowness();
            t = (b*dx2 + a*dz2)/(dx2 + dz2) + sqrt((2.0*a*b*dx2*dz2 - a*a*dx2*dz2 -
                                                    b*b*dx2*dz2 + dx2*dx2*dz2*s2 +
                                                    dx2*dz2*dz2*s2)/((dx2 + dz2)*(dx2 + dz2)));
        }

        if ( t<nodes[i*(ncz+1)+j].getTT(threadNo) )
            nodes[i*(ncz+1)+j].setTT(t,threadNo);
    }


    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::initFSM(const std::vector<S>& Tx,
                                         const std::vector<T1>& t0,
                                         std::vector<bool>& frozen,
                                         const int npts,
                                         const size_t threadNo) const {

        for (size_t n=0; n<Tx.size(); ++n) {
            bool found = false;
            for ( long long nn=0; nn<nodes.size(); ++nn ) {
                if ( nodes[nn] == Tx[n] ) {
                    found = true;
                    nodes[nn].setTT( t0[n], threadNo );
                    frozen[nn] = true;

                    long long i = nn/(ncz+1);
                    long long j = nn - i*(ncz+1);

                    for ( long long ii=i-npts; ii<=i+npts; ++ii ) {
                        if ( ii>=0 && ii<=ncx ) {
                            for ( long long jj=j-npts; jj<=j+npts; ++jj ) {
                                if ( jj>=0 && jj<=ncz && !(ii==i && jj==j) ) {

                                    size_t nnn = ii*(ncz+1) + jj;

                                    T1 tt = t0[n] + nodes[nnn].getDistance(Tx[n]) * 0.5*(nodes[nnn].getNodeSlowness() + nodes[nn].getNodeSlowness());
                                    nodes[nnn].setTT( tt, threadNo );
                                    frozen[nnn] = true;
                                }
                            }
                        }
                    }

                    break;
                }
            }
            if ( found==false ) {

                // find cell where Tx resides
                long long cellNo = getCellNo(Tx[n]);

                long long i = cellNo/ncz;
                long long j = cellNo - i*ncz;

                for ( long long ii=i-(npts-1); ii<=i+npts; ++ii ) {
                    if ( ii>=0 && ii<=ncx ) {
                        for ( long long jj=j-(npts-1); jj<=j+npts; ++jj ) {
                            if ( jj>=0 && jj<=ncz ) {
                                size_t nnn = ii*(ncz+1) + jj;

                                T1 tt = t0[n] + nodes[nnn].getDistance(Tx[n]) * nodes[nnn].getNodeSlowness();
                                nodes[nnn].setTT( tt, threadNo );
                                frozen[nnn] = true;
                            }
                        }
                    }
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE>
    T1 Grid2Drn<T1,T2,S,NODE>::getSlowness(const S& pt) const {

        // bilinear interpolation if not on node

        T1 s;
        T2 i, j;

        getIJ(pt, i, j);

        if ( std::abs(pt.x - (xmin+i*dx))<small && std::abs(pt.z - (zmin+j*dz))<small ) {
            // on node
            return nodes[i*(ncz+1)+j].getNodeSlowness();
        } else if ( std::abs(pt.x - (xmin+i*dx))<small ) {

            // on edge
            T1 t1 = nodes[i*(ncz+1)+j].getNodeSlowness();
            T1 t2 = nodes[i*(ncz+1)+j+1].getNodeSlowness();

            T1 w1 = (zmin+(j+1)*dz - pt.z)/dz;
            T1 w2 = (pt.z - (zmin+j*dz))/dz;

            s = t1*w1 + t2*w2;

        } else if ( std::abs(pt.z - (zmin+j*dz))<small ) {

            // on edge
            T1 t1 = nodes[i*(ncz+1)+j].getNodeSlowness();
            T1 t2 = nodes[(i+1)*(ncz+1)+j].getNodeSlowness();

            T1 w1 = (xmin+(i+1)*dx - pt.x)/dx;
            T1 w2 = (pt.x - (xmin+i*dx))/dx;

            s = t1*w1 + t2*w2;

        } else {

            T1 t1 = nodes[    i*(ncz+1)+j  ].getNodeSlowness();
            T1 t2 = nodes[(i+1)*(ncz+1)+j  ].getNodeSlowness();
            T1 t3 = nodes[    i*(ncz+1)+j+1].getNodeSlowness();
            T1 t4 = nodes[(i+1)*(ncz+1)+j+1].getNodeSlowness();

            T1 w1 = (xmin+(i+1)*dx - pt.x)/dx;
            T1 w2 = (pt.x - (xmin+i*dx))/dx;

            t1 = t1*w1 + t2*w2;
            t2 = t3*w1 + t4*w2;

            w1 = (zmin+(j+1)*dz - pt.z)/dz;
            w2 = (pt.z - (zmin+j*dz))/dz;

            s = t1*w1 + t2*w2;
        }

        return s;
    }

    template<typename T1, typename T2, typename S, typename NODE>
    T1 Grid2Drn<T1,T2,S,NODE>::getTraveltimeFromRaypath(const std::vector<S>& Tx,
                                                        const std::vector<T1>& t0,
                                                        const S& Rx,
                                                        const size_t threadNo) const {

        T1 tt = 0.0;
        T1 s1, s2;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return tt;
            }
        }

        S prev_pt( Rx );
        S curr_pt( Rx );
        s1 = getSlowness( curr_pt );
        // distance between opposite nodes of a voxel
        const T1 maxDist = sqrt( dx*dx + dz*dz );
        S g;

        bool reachedTx = false;
        while ( reachedTx == false ) {

            grad(g, curr_pt, threadNo);
            g *= -1.0;

            long long i, k;
            getIJ(curr_pt, i, k);

            // planes we will intersect
            T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
            T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

            if ( std::abs(xp-curr_pt.x)<small) {
                xp += dx*boost::math::sign(g.x);
            }
            if ( std::abs(zp-curr_pt.z)<small) {
                zp += dz*boost::math::sign(g.z);
            }

            // dist to planes
            T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
            T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

            if ( tx<tz ) { // closer to xp
                curr_pt += tx*g;
                curr_pt.x = xp;     // make sure we don't accumulate rounding errors
            } else {
                curr_pt += tz*g;
                curr_pt.z = zp;
            }

            if ( curr_pt.x < xmin || curr_pt.x > xmax ||
                curr_pt.z < zmin || curr_pt.z > zmax ) {
                //  we are going oustide the grid!
                // make gardient point along outside face
                if ( abs(g.x) > abs(g.z) ) {
                    g.x = boost::math::sign(g.x);
                    g.z = 0.0;
                } else {
                    g.x = 0.0;
                    g.z = boost::math::sign(g.z);
                }

                // put back previous coordinates
                curr_pt = prev_pt;

                // planes we will intersect
                T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
                T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

                if ( std::abs(xp-curr_pt.x)<small) {
                    xp += dx*boost::math::sign(g.x);
                }
                if ( std::abs(zp-curr_pt.z)<small) {
                    zp += dz*boost::math::sign(g.z);
                }

                // dist to planes
                T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
                T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

                if ( tx<tz ) { // closer to xp
                    curr_pt += tx*g;
                    curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                } else {
                    curr_pt += tz*g;
                    curr_pt.z = zp;
                }

                if ( curr_pt.x < xmin || curr_pt.x > xmax ||
                    curr_pt.z < zmin || curr_pt.z > zmax ) {
                    std::ostringstream msg;
                    msg << "Error while computing raypaths: going outside grid \n\
                    Rx: " << Rx << "\n\
                    Tx: " << Tx[0] << "\n";
                    for ( size_t ns=1; ns<Tx.size(); ++ns ) {
                        msg << "\
                        " << Tx[ns] << "\n";
                    }
                    throw std::runtime_error(msg.str());
                }
            }

            s2 = getSlowness( curr_pt );
            tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
            s1 = s2;
            prev_pt = curr_pt;

            // are we close enough to one the Tx nodes ?
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                T1 dist = curr_pt.getDistance( Tx[ns] );
                if ( dist < maxDist ) {

                    g = Tx[ns] - curr_pt;

                    T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
                    T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

                    if ( std::abs(xp-curr_pt.x)<small) {
                        xp += dx*boost::math::sign(g.x);
                    }
                    if ( std::abs(zp-curr_pt.z)<small) {
                        zp += dz*boost::math::sign(g.z);
                    }

                    // dist to planes
                    T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
                    T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

                    if ( tx<tz ) { // closer to xp
                        curr_pt += tx*g;
                        curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                    } else {
                        curr_pt += tz*g;
                        curr_pt.z = zp;
                    }

                    if ( curr_pt.getDistance(prev_pt) > dist ||  // we do not intersect
                        curr_pt == Tx[ns] ) {  // we have arrived
                        s2 = getSlowness( Tx[ns] );
                        tt += 0.5*(s1 + s2) * prev_pt.getDistance( Tx[ns] );
                    } else {
                        // to intersection
                        s2 = getSlowness( curr_pt );
                        tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
                        s1 = s2;
                        // to Tx
                        s2 = getSlowness( Tx[ns] );
                        tt += 0.5*(s1 + s2) * curr_pt.getDistance( Tx[ns] );
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
        return tt;
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::getRaypath(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            const S& Rx,
                                            std::vector<S>& r_data,
                                            T1 &tt,
                                            const size_t threadNo) const {

        tt = 0.0;
        T1 s1, s2;

        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        S curr_pt( Rx );
        s1 = getSlowness( curr_pt );
        // distance between opposite nodes of a voxel
        const T1 maxDist = sqrt( dx*dx + dz*dz );
        S g;

        bool reachedTx = false;
        while ( reachedTx == false ) {

            grad(g, curr_pt, threadNo);
            g *= -1.0;

            long long i, k;
            getIJ(curr_pt, i, k);

            // planes we will intersect
            T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
            T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

            if ( std::abs(xp-curr_pt.x)<small) {
                xp += dx*boost::math::sign(g.x);
            }
            if ( std::abs(zp-curr_pt.z)<small) {
                zp += dz*boost::math::sign(g.z);
            }

            // dist to planes
            T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
            T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

            if ( tx<tz ) { // closer to xp
                curr_pt += tx*g;
                curr_pt.x = xp;     // make sure we don't accumulate rounding errors
            } else {
                curr_pt += tz*g;
                curr_pt.z = zp;
            }

            if ( curr_pt.x < xmin || curr_pt.x > xmax ||
                curr_pt.z < zmin || curr_pt.z > zmax ) {
                //  we are going oustide the grid!
                // make gardient point along outside face
                if ( abs(g.x) > abs(g.z) ) {
                    g.x = boost::math::sign(g.x);
                    g.z = 0.0;
                } else {
                    g.x = 0.0;
                    g.z = boost::math::sign(g.z);
                }

                // put back previous coordinates
                curr_pt = r_data.back();

                // planes we will intersect
                T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
                T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

                if ( std::abs(xp-curr_pt.x)<small) {
                    xp += dx*boost::math::sign(g.x);
                }
                if ( std::abs(zp-curr_pt.z)<small) {
                    zp += dz*boost::math::sign(g.z);
                }

                // dist to planes
                T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
                T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

                if ( tx<tz ) { // closer to xp
                    curr_pt += tx*g;
                    curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                } else {
                    curr_pt += tz*g;
                    curr_pt.z = zp;
                }

                if ( curr_pt.x < xmin || curr_pt.x > xmax ||
                    curr_pt.z < zmin || curr_pt.z > zmax ) {
                    std::ostringstream msg;
                    msg << "Error while computing raypaths: going outside grid \n\
                    Rx: " << Rx << "\n\
                    Tx: " << Tx[0] << "\n";
                    for ( size_t ns=1; ns<Tx.size(); ++ns ) {
                        msg << "\
                        " << Tx[ns] << "\n";
                    }
                    throw std::runtime_error(msg.str());
                }
            }

            s2 = getSlowness( curr_pt );
            tt += 0.5*(s1 + s2) * r_data.back().getDistance( curr_pt );
            s1 = s2;
            r_data.push_back( curr_pt );

            // are we close enough to one the Tx nodes ?
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                T1 dist = curr_pt.getDistance( Tx[ns] );
                if ( dist < maxDist ) {

                    g = Tx[ns] - curr_pt;

                    T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
                    T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

                    if ( std::abs(xp-curr_pt.x)<small) {
                        xp += dx*boost::math::sign(g.x);
                    }
                    if ( std::abs(zp-curr_pt.z)<small) {
                        zp += dz*boost::math::sign(g.z);
                    }

                    // dist to planes
                    T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
                    T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

                    if ( tx<tz ) { // closer to xp
                        curr_pt += tx*g;
                        curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                    } else {
                        curr_pt += tz*g;
                        curr_pt.z = zp;
                    }

                    if ( curr_pt.getDistance(r_data.back()) > dist ||  // we do not intersect
                        curr_pt == Tx[ns] ) {  // we have arrived
                        s2 = getSlowness( Tx[ns] );
                        tt += 0.5*(s1 + s2) * r_data.back().getDistance( Tx[ns] );
                        r_data.push_back( Tx[ns] );
                    } else {
                        // to intersection
                        s2 = getSlowness( curr_pt );
                        tt += 0.5*(s1 + s2) * r_data.back().getDistance( curr_pt );
                        r_data.push_back( curr_pt );
                        s1 = s2;
                        // to Tx
                        s2 = getSlowness( Tx[ns] );
                        tt += 0.5*(s1 + s2) * curr_pt.getDistance( Tx[ns] );
                        r_data.push_back( Tx[ns] );
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::getRaypath(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            const S& Rx,
                                            std::vector<S>& r_data,
                                            std::vector<siv<T1>> &l_data,
                                            T1 &tt,
                                            const size_t threadNo) const {

        tt = 0.0;
        T1 s1, s2;

        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        S curr_pt( Rx );
        s1 = getSlowness( curr_pt );
        // distance between opposite nodes of a voxel
        const T1 maxDist = sqrt( dx*dx + dz*dz );
        S g;

        siv<T1> cell;
        bool reachedTx = false;
        while ( reachedTx == false ) {

            grad(g, curr_pt, threadNo);
            g *= -1.0;

            long long i, k;
            getIJ(curr_pt, i, k);

            // planes we will intersect
            T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
            T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

            if ( std::abs(xp-curr_pt.x)<small) {
                xp += dx*boost::math::sign(g.x);
            }
            if ( std::abs(zp-curr_pt.z)<small) {
                zp += dz*boost::math::sign(g.z);
            }

            // dist to planes
            T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
            T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

            if ( tx<tz ) { // closer to xp
                curr_pt += tx*g;
                curr_pt.x = xp;     // make sure we don't accumulate rounding errors
            } else {
                curr_pt += tz*g;
                curr_pt.z = zp;
            }

            if ( curr_pt.x < xmin || curr_pt.x > xmax ||
                curr_pt.z < zmin || curr_pt.z > zmax ) {
                //  we are going oustide the grid!
                std::ostringstream msg;
                msg << "Error while computing raypaths: going outside grid \n\
                Rx: " << Rx << "\n\
                Tx: " << Tx[0] << "\n";
                for ( size_t ns=1; ns<Tx.size(); ++ns ) {
                    msg << "\
                    " << Tx[ns] << "\n";
                }
                throw std::runtime_error(msg.str());
            }
            S mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
            cell.i = getCellNo(mid_pt);
            cell.v = curr_pt.getDistance(r_data.back());
            l_data.push_back(cell);
            s2 = getSlowness( curr_pt );
            tt += 0.5*(s1 + s2) * r_data.back().getDistance( curr_pt );
            s1 = s2;
            r_data.push_back( curr_pt );

            // are we close enough to one the Tx nodes ?
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                T1 dist = curr_pt.getDistance( Tx[ns] );
                if ( dist < maxDist ) {

                    g = Tx[ns] - curr_pt;

                    T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
                    T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

                    if ( std::abs(xp-curr_pt.x)<small) {
                        xp += dx*boost::math::sign(g.x);
                    }
                    if ( std::abs(zp-curr_pt.z)<small) {
                        zp += dz*boost::math::sign(g.z);
                    }

                    // dist to planes
                    T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
                    T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

                    if ( tx<tz ) { // closer to xp
                        curr_pt += tx*g;
                        curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                    } else {
                        curr_pt += tz*g;
                        curr_pt.z = zp;
                    }

                    if ( curr_pt.getDistance(r_data.back()) > dist ||  // we do not intersect
                        curr_pt == Tx[ns] ) {  // we have arrived
                        cell.i = getCellNo(Tx[ns]);
                        cell.v = Tx[ns].getDistance(r_data.back());
                        l_data.push_back(cell);s2 = getSlowness( Tx[ns] );
                        tt += 0.5*(s1 + s2) * r_data.back().getDistance( Tx[ns] );
                        r_data.push_back( Tx[ns] );
                    } else {
                        // to intersection
                        mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
                        cell.i = getCellNo(mid_pt);
                        cell.v = curr_pt.getDistance(r_data.back());
                        l_data.push_back(cell);s2 = getSlowness( curr_pt );
                        tt += 0.5*(s1 + s2) * r_data.back().getDistance( curr_pt );
                        r_data.push_back( curr_pt );
                        s1 = s2;
                        // to Tx
                        cell.i = getCellNo(Tx[ns]);
                        cell.v = Tx[ns].getDistance(r_data.back());
                        l_data.push_back(cell);s2 = getSlowness( Tx[ns] );
                        tt += 0.5*(s1 + s2) * curr_pt.getDistance( Tx[ns] );
                        r_data.push_back( Tx[ns] );
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
    }
}

#endif

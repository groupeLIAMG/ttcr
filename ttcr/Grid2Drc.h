//
//  Grid2Drc.h
//  ttcr
//
//  Created by Bernard Giroux on 15-12-23.
//  Copyright © 2015 Bernard Giroux. All rights reserved.
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

#ifndef ttcr_Grid2Drc_h
#define ttcr_Grid2Drc_h

#include <algorithm>
#include <cstring>
#include <exception>
#include <iostream>
#include <fstream>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifdef VTK
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkXMLRectilinearGridWriter.h"
#endif

#include <boost/math/special_functions/sign.hpp>

#include "Grid2D.h"

namespace ttcr {

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    class Grid2Drc : public Grid2D<T1,T2,S> {
    public:
        Grid2Drc(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
                 const T1 minx, const T1 minz, const bool ttrp,
                 const size_t nt=1) :
        Grid2D<T1,T2,S>(nx*nz, ttrp, nt),
        dx(ddx), dz(ddz), xmin(minx), zmin(minz),
        xmax(minx+nx*ddx), zmax(minz+nz*ddz),
        ncx(nx), ncz(nz),
        nodes(std::vector<NODE>( (ncx+1) * (ncz+1), NODE(nt) )),
        cells(ncx*ncz)
        {
        }

        virtual ~Grid2Drc() {
        }

        void getSlowness(std::vector<T1>& slowness) const {
            if (slowness.size() != ncx*ncz) {
                slowness.resize(ncx*ncz);
            }
            for (size_t n=0; n<slowness.size(); ++n) {
                slowness[n] = cells.getSlowness(n);
            }
        }

        void setSlowness(const std::vector<T1>& s) {
            try {
                cells.setSlowness( s );
            } catch (std::exception& e) {
                throw;
            }
        }
        void setXi(const std::vector<T1>& x) {
            try {
                cells.setXi( x );
            } catch (std::exception& e) {
                throw;
            }
        }
        void setTiltAngle(const std::vector<T1>& t) {
            try {
                cells.setTiltAngle( t );
            } catch (std::exception& e) {
                throw;
            }
        }
        void setVp0(const std::vector<T1>& s) {
            try {
                cells.setVp0(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        void setVs0(const std::vector<T1>& s) {
            try {
                cells.setVs0(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        void setDelta(const std::vector<T1>& s) {
            try {
                cells.setDelta(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        void setEpsilon(const std::vector<T1>& s) {
            try {
                cells.setEpsilon(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        void setGamma(const std::vector<T1>& s) {
            try {
                cells.setGamma(s);
            } catch (std::exception& e) {
                throw;
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

        T2 getCellNo(const S& pt) const {
            T1 x = xmax-pt.x < small ? xmax-.5*dx : pt.x;
            T1 z = zmax-pt.z < small ? zmax-.5*dz : pt.z;
            T2 nx = static_cast<T2>( small + (x-xmin)/dx );
            T2 nz = static_cast<T2>( small + (z-zmin)/dz );
            return nx*ncz + nz;
        }

        void saveSlownessXYZ(const char filename[]) const {
            std::ofstream fout( filename );

            for ( T2 j=0, n=0; j<ncx; ++j ) {
                T1 x = xmin + (0.5+j)*dx;
                for ( T2 i=0; i<ncz; ++i, ++n ) {
                    T1 z = zmin + (0.5+i)*dz;
                    fout << x << "   " << z << "   " << cells.getSlowness(n) << '\n';
                }
            }

            fout.close();
        }

        void saveTT(const std::string &, const int, const size_t nt=0,
                    const int format=1) const;

        void saveTTgrad(const std::string &, const size_t nt=0,
                        const bool vtkFormat=0) const {}

        const T1 getXmin() const { return xmin; }
        const T1 getZmin() const { return zmin; }
        const T1 getDx() const { return dx; }
        const T1 getDz() const { return dz; }
        const T2 getNcx() const { return ncx; }
        const T2 getNcz() const { return ncz; }

        void dump_secondary(std::ofstream& os) const {
            size_t nPrimary = (ncx+1) * (ncz+1);
            for ( size_t n=nPrimary; n<nodes.size(); ++n ) {
                os << nodes[n].getX() << ' ' << nodes[n].getZ() << '\n';
            }
        }

    protected:
        T1 dx;           // cell size in x
        T1 dz;           // cell size in z
        T1 xmin;         // x origin of the grid
        T1 zmin;         // z origin of the grid
        T1 xmax;         // x end of the grid
        T1 zmax;         // z end of the grid
        T2 ncx;          // number of cells in x
        T2 ncz;          // number of cells in z

        mutable std::vector<NODE> nodes;

        CELL cells;   // column-wise (z axis) slowness vector of the cells

        void checkPts(const std::vector<S>&) const;

        bool inPolygon(const S& p, const S poly[], const size_t N) const;

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

        T1 getTraveltime(const S& Rx, const size_t threadNo) const;

    private:
        void grad(S &g, const S &pt, const size_t nt) const;
        T1 interpolateTraveltime(const S& pt, const size_t threadNo) const;

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

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Drc<T1,T2,S,NODE,CELL>::saveTT(const std::string& fname,
                                             const int all,
                                             const size_t nt,
                                             const int format) const {

        if ( format == 1 ) {
            std::string filename = fname+".dat";
            std::ofstream fout(filename.c_str());
            fout.precision(12);
            for ( T2 n=0; n<nodes.size(); ++n ) {
                fout << nodes[n].getX() << '\t'
                << nodes[n].getZ() << '\t'
                << nodes[n].getTT(nt) << '\n';
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
        } else if ( format == 3 ) {
            std::string filename = fname+".bin";
            std::ofstream fout(filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
            for ( T2 n=0; n<nodes.size(); ++n ) {
                T1 tmp[] = { nodes[n].getX(), nodes[n].getZ(), nodes[n].getTT(nt) };
                fout.write( (char*)tmp, 3*sizeof(T1) );
            }
            fout.close();
        } else {
            throw std::runtime_error("Unsupported format for saving traveltimes");
        }
    }


    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Drc<T1,T2,S,NODE,CELL>::checkPts(const std::vector<S>& pts) const {
        for (size_t n=0; n<pts.size(); ++n) {
            if ( pts[n].x < xmin || pts[n].x > xmax ||
                pts[n].z < zmin || pts[n].z > zmax ) {
                std::ostringstream msg;
                msg << "Error: Point (" << pts[n].x << ", "<< pts[n] .z << ") outside grid.";
                throw std::runtime_error(msg.str());
            }
        }
    }


    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    bool Grid2Drc<T1,T2,S,NODE,CELL>::inPolygon(const S& p, const S poly[],
                                                const size_t N) const {
        bool c = false;
        for (size_t i = 0, j = N-1; i < N; j = i++) {
            if ((((poly[i].z <= p.z) && (p.z < poly[j].z)) ||
                 ((poly[j].z <= p.z) && (p.z < poly[i].z))) &&
                (p.x < (poly[j].x - poly[i].x) * (p.z - poly[i].z) / (poly[j].z - poly[i].z) + poly[i].x))
                c = !c;
        }
        return c;
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    T1 Grid2Drc<T1,T2,S,NODE,CELL>::interpolateTraveltime(const S& pt,
                                                          const size_t nt) const {

        const size_t nnz = ncz+1;

        // bilinear interpolation if not on node

        T1 tt;
        T2 i, k;

        getIJ(pt, i, k);
        if ( std::abs(pt.x - (xmin+i*dx))<small &&
            std::abs(pt.z - (zmin+k*dz))<small ) {
            // on node
            return nodes[i*nnz+k].getTT(nt);
        } else if ( std::abs(pt.x - (xmin+i*dx))<small ) {
            // on edge
            T1 t1 = nodes[i*nnz+k  ].getTT(nt);
            T1 t2 = nodes[i*nnz+k+1].getTT(nt);

            T1 w1 = (zmin+(k+1)*dz - pt.z)/dz;
            T1 w2 = (pt.z - (zmin+k*dz))/dz;

            tt = t1*w1 + t2*w2;
        } else if ( std::abs(pt.z - (zmin+k*dz))<small ) {
            // on edge
            T1 t1 = nodes[    i*nnz+k].getTT(nt);
            T1 t2 = nodes[(i+1)*nnz+k].getTT(nt);

            T1 w1 = (xmin+(i+1)*dx - pt.x)/dx;
            T1 w2 = (pt.x - (xmin+i*dx))/dx;

            tt = t1*w1 + t2*w2;
        } else {
            T1 t1 = nodes[    i*nnz+k  ].getTT(nt);
            T1 t2 = nodes[    i*nnz+k+1].getTT(nt);
            T1 t3 = nodes[(i+1)*nnz+k  ].getTT(nt);
            T1 t4 = nodes[(i+1)*nnz+k+1].getTT(nt);

            T1 w1 = (zmin+(k+1)*dz - pt.z)/dz;
            T1 w2 = (pt.z - (zmin+k*dz))/dz;

            t1 = t1*w1 + t2*w2;
            t2 = t3*w1 + t4*w2;

            w1 = (xmin+(i+1)*dx - pt.x)/dx;
            w2 = (pt.x - (xmin+i*dx))/dx;

            tt = t1*w1 + t2*w2;
        }
        return tt;
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Drc<T1,T2,S,NODE,CELL>::grad(S &g, const S &pt,
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
        g.x = (interpolateTraveltime({p2, pt.z}, nt) - interpolateTraveltime({p1, pt.z}, nt)) / dx;

        p1 = pt.z - dz/2.0;
        if (p1 < zmin) {
            p1 = zmin;
        }
        p2 = p1 + dz;
        if (p2 > zmax) {
            p2 = zmax;
            p1 = zmax - dz;
        }
        g.z = (interpolateTraveltime({pt.x, p2}, nt) - interpolateTraveltime({pt.x, p1}, nt)) / dz;
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    T1 Grid2Drc<T1,T2,S,NODE,CELL>::getTraveltime(const S& Rx, const size_t threadNo) const {

        for ( size_t nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == Rx ) {
                return nodes[nn].getTT(threadNo);
            }
        }

        T2 cellNo = getCellNo( Rx );
        T2 neibNo = this->neighbors[cellNo][0];
        T1 dt = cells.computeDt(nodes[neibNo], Rx, cellNo);

        T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
        for ( size_t k=1; k< this->neighbors[cellNo].size(); ++k ) {
            neibNo = this->neighbors[cellNo][k];
            dt = cells.computeDt(nodes[neibNo], Rx, cellNo);
            if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
                traveltime =  nodes[neibNo].getTT(threadNo)+dt;
            }
        }
        return traveltime;
    }


    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    T1 Grid2Drc<T1,T2,S,NODE,CELL>::getTraveltimeFromRaypath(const std::vector<S>& Tx,
                                                             const std::vector<T1>& t0,
                                                             const S& Rx,
                                                             const size_t threadNo) const {

        T1 tt = 0.0;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return tt;
            }
        }

        S prev_pt( Rx );
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

            S mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
            T2 cellNo = getCellNo(mid_pt);
            tt += cells.computeDt(prev_pt, curr_pt, cellNo);
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
                        cellNo = getCellNo(Tx[ns]);
                        tt += cells.computeDt(Tx[ns], prev_pt, cellNo);
                    } else {
                        // to intersection
                        mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
                        cellNo = getCellNo(mid_pt);
                        tt += cells.computeDt(prev_pt, curr_pt, cellNo);
                        // to Tx
                        cellNo = getCellNo(Tx[ns]);
                        tt += cells.computeDt(Tx[ns], curr_pt, cellNo);
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
        return tt;
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Drc<T1,T2,S,NODE,CELL>::getRaypath(const std::vector<S>& Tx,
                                                 const std::vector<T1>& t0,
                                                 const S& Rx,
                                                 std::vector<S>& r_data,
                                                 T1 &tt,
                                                 const size_t threadNo) const {
        tt = 0.0;
        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
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

            S mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
            T2 cellNo = getCellNo(mid_pt);
            tt += cells.computeDt(r_data.back(), curr_pt, cellNo);
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
                        cellNo = getCellNo(Tx[ns]);
                        tt += cells.computeDt(Tx[ns], r_data.back(), cellNo);
                        r_data.push_back( Tx[ns] );
                    } else {
                        // to intersection
                        mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
                        cellNo = getCellNo(mid_pt);
                        tt += cells.computeDt(r_data.back(), curr_pt, cellNo);
                        r_data.push_back( curr_pt );
                        // to Tx
                        cellNo = getCellNo(Tx[ns]);
                        tt += cells.computeDt(Tx[ns], curr_pt, cellNo);
                        r_data.push_back( Tx[ns] );
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Drc<T1,T2,S,NODE,CELL>::getRaypath(const std::vector<S>& Tx,
                                                 const std::vector<T1>& t0,
                                                 const S& Rx,
                                                 std::vector<S>& r_data,
                                                 std::vector<siv<T1>> &l_data,
                                                 T1 &tt,
                                                 const size_t threadNo) const {
        tt = 0.0;
        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        S curr_pt( Rx );
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
            tt += cells.computeDt(r_data.back(), curr_pt, cell.i);
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
                        l_data.push_back(cell);
                        tt += cells.computeDt(Tx[ns], r_data.back(), cell.i);
                        r_data.push_back( Tx[ns] );
                    } else {
                        // to intersection
                        mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
                        cell.i = getCellNo(mid_pt);
                        cell.v = curr_pt.getDistance(r_data.back());
                        l_data.push_back(cell);
                        tt += cells.computeDt(r_data.back(), curr_pt, cell.i);
                        r_data.push_back( curr_pt );
                        // to Tx
                        cell.i = getCellNo(Tx[ns]);
                        cell.v = Tx[ns].getDistance(r_data.back());
                        l_data.push_back(cell);
                        tt += cells.computeDt(Tx[ns], curr_pt, cell.i);
                        r_data.push_back( Tx[ns] );
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
    }
}

#endif /* Grid2Drc_h */

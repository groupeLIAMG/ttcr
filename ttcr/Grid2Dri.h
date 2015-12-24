//
//  Grid2Dri.h
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

#ifndef ttcr_Grid2Dri_h
#define ttcr_Grid2Dri_h

#include <algorithm>
#include <cstring>
#include <iostream>
#include <fstream>
#include <queue>
#include <vector>

#ifdef VTK
#include "vtkDoubleArray.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkXMLRectilinearGridWriter.h"
#endif

#include "Grid2D.h"

template<typename T1, typename T2, typename NODE>
class Grid2Dri : public Grid2D<T1,T2,sxz<T1>> {
public:
    Grid2Dri(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
			 const T1 minx, const T1 minz, const size_t nt=1);
    
    virtual ~Grid2Dri() {
    }
    
    T1 getDx() const { return dx; }
    T1 getDz() const { return dz; }
    T1 getXmin() const { return xmin; }
    T1 getXmax() const { return xmax; }
    T1 getZmin() const { return zmin; }
    T1 getZmax() const { return zmax; }
    T2 getNcellx() const { return nCellx; }
    T2 getNcellz() const { return nCellz; }
    
    void setSlowness(const T1 s) {
        for ( size_t n=0; n<nodes.size(); ++n ) {
        	nodes[n].setNodeSlowness( s );
        }
    }
    
    virtual int setSlowness(const std::vector<T1>& s) {
        if ( nodes.size() != s.size() ) {
            std::cerr << "Error: slowness vectors of incompatible size.";
            return 1;
        }
        for ( size_t n=0; n<nodes.size(); ++n ) {
            nodes[n].setNodeSlowness( s[n] );
        }
        return 0;
    }
    
    virtual int raytrace(const std::vector<sxz<T1>>& Tx,
                         const std::vector<T1>& t0,
                         const std::vector<sxz<T1>>& Rx,
                         std::vector<T1>& traveltimes,
                         const size_t threadNo=0) const { return 0; }
    
    virtual int raytrace(const std::vector<sxz<T1>>&,
                         const std::vector<T1>&,
                         const std::vector<const std::vector<sxz<T1>>*>&,
                         std::vector<std::vector<T1>*>&,
                         const size_t=0) const { return 0; }
    
    virtual int raytrace(const std::vector<sxz<T1>>& Tx,
                         const std::vector<T1>& t0,
                         const std::vector<sxz<T1>>& Rx,
                         std::vector<T1>& traveltimes,
                         std::vector<std::vector<sxz<T1>>>& r_data,
                         const size_t threadNo=0) const { return 0; }
    
    virtual int raytrace(const std::vector<sxz<T1>>&,
                         const std::vector<T1>&,
                         const std::vector<const std::vector<sxz<T1>>*>&,
                         std::vector<std::vector<T1>*>&,
                         std::vector<std::vector<std::vector<sxz<T1>>>*>&,
                         const size_t=0) const { return 0; }
	
    size_t getNumberOfNodes() const { return nodes.size(); }
    
    void saveSlownessXYZ(const char filename[]) const {
        std::ofstream fout( filename );
        
		for ( size_t n=0; n<nodes.size(); ++n ) {
			fout << nodes[n].getX() << "   " << nodes[n].getZ() << "   " << nodes[n].getNodeSlowness() << '\n';
        }
        
        fout.close();
    }
    
    T2 getCellNo(const sxz<T1>& pt) const {
        T1 x = xmax-pt.x < small ? xmax-.5*dx : pt.x;
        T1 z = zmax-pt.z < small ? zmax-.5*dz : pt.z;
        T2 nx = static_cast<T2>( small + (x-xmin)/dx );
        T2 nz = static_cast<T2>( small + (z-zmin)/dz );
        return nx*nCellz + nz;
    }
    
    void getIJ(const sxz<T1>& pt, T2& i, T2& j) const {
        i = static_cast<T2>( small + (pt.x-xmin)/dx );
        j = static_cast<T2>( small + (pt.z-zmin)/dz );
    }
    
    void getIJ(const sxz<T1>& pt, long long& i, long long& j) const {
        i = static_cast<long long>( small + (pt.x-xmin)/dx );
        j = static_cast<long long>( small + (pt.z-zmin)/dz );
    }
    
    void saveTT(const std::string &, const int, const size_t nt=0,
                const bool vtkFormat=0) const;
    
protected:
	size_t nThreads;
    T1 dx;           // cell size in x
    T1 dz;           // cell size in z
    T1 xmin;         // x origin of the grid
    T1 zmin;         // z origin of the grid
    T1 xmax;         // x end of the grid
    T1 zmax;         // z end of the grid
    T2 nCellx;       // number of cells in x
    T2 nCellz;       // number of cells in x
    
    mutable std::vector<NODE> nodes;
    
    std::vector<std::vector<T2>> neighbors;  // nodes common to a cell
    void buildGridNeighbors();
	
	T1 computeDt(const NODE& source, const NODE& node) const {
		return (node.getNodeSlowness()+source.getNodeSlowness())/2 * source.getDistance( node );
	}
    
	T1 computeDt(const NODE& source, const sxz<T1>& node, T1 slo) const {
		return (slo+source.getNodeSlowness())/2 * source.getDistance( node );
	}
    
    int check_pts(const std::vector<sxz<T1>>&) const;
    
    bool inPolygon(const sxz<T1>& p, const sxz<T1> poly[], const size_t N) const;
    
    T1 getTraveltime(const sxz<T1>& Rx, const std::vector<NODE>& nodes,
					 const size_t threadNo) const;
    
    T1 getTraveltime(const sxz<T1>& Rx, const std::vector<NODE>& nodes,
					 T2& nodeParentRx, T2& cellParentRx,
					 const size_t threadNo) const;
	
    void grad(sxz<T1> &g, const size_t i, const size_t j, const size_t nt=0) const;
    
    void getRaypath(const std::vector<sxz<T1>>& Tx,
                    const sxz<T1> &Rx,
                    std::vector<sxz<T1>> &r_data,
                    const size_t threadNo=0) const;
    
private:
    Grid2Dri() {}
    Grid2Dri(const Grid2Dri<T1,T2,NODE>& g) {}
    Grid2Dri<T1,T2,NODE>& operator=(const Grid2Dri<T1,T2,NODE>& g) {}
    
};

template<typename T1, typename T2, typename NODE>
Grid2Dri<T1,T2,NODE>::Grid2Dri(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
						  const T1 minx, const T1 minz, const size_t nt) : nThreads(nt),
dx(ddx), dz(ddz), xmin(minx), zmin(minz), xmax(minx+nx*ddx), zmax(minz+nz*ddz),
nCellx(nx), nCellz(nz),
nodes(std::vector<NODE>( (nCellx+1) * (nCellz+1), NODE(nt) )),
neighbors(std::vector<std::vector<T2>>(nCellx*nCellz))
{ }

template<typename T1, typename T2, typename NODE>
void Grid2Dri<T1,T2,NODE>::buildGridNeighbors() {
    for ( T2 n=0; n<nodes.size(); ++n ) {
        for ( size_t n2=0; n2<nodes[n].getOwners().size(); ++n2) {
            neighbors[ nodes[n].getOwners()[n2] ].push_back(n);
        }
    }
}

template<typename T1, typename T2, typename NODE>
int Grid2Dri<T1,T2,NODE>::check_pts(const std::vector<sxz<T1>>& pts) const {
    for (size_t n=0; n<pts.size(); ++n) {
        if ( pts[n].x < xmin || pts[n].x > xmax ||
            pts[n].z < zmin || pts[n].z > zmax ) {
            std::cerr << "Error: point no " << (n+1)
            << " outside the grid.\n";
            return 1;
        }
    }
    return 0;
}



template<typename T1, typename T2, typename NODE>
bool Grid2Dri<T1,T2,NODE>::inPolygon(const sxz<T1>& p, const sxz<T1> poly[], const size_t N) const {
    bool c = false;
    for (size_t i = 0, j = N-1; i < N; j = i++) {
        if ((((poly[i].z <= p.z) && (p.z < poly[j].z)) ||
             ((poly[j].z <= p.z) && (p.z < poly[i].z))) &&
            (p.x < (poly[j].x - poly[i].x) * (p.z - poly[i].z) / (poly[j].z - poly[i].z) + poly[i].x))
            c = !c;
    }
    return c;
}


template<typename T1, typename T2, typename NODE>
T1 Grid2Dri<T1,T2,NODE>::getTraveltime(const sxz<T1>& Rx,
								  const std::vector<NODE>& nodes,
								  const size_t threadNo) const {
    
    for ( size_t nn=0; nn<nodes.size(); ++nn ) {
        if ( nodes[nn] == Rx ) {
            return nodes[nn].getTT(threadNo);
        }
    }
    
    T2 cellNo = getCellNo( Rx );
    T2 neibNo = neighbors[cellNo][0];
    T1 dt = computeDt(nodes[neibNo], Rx, cellNo);
    
    T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
    for ( size_t k=1; k< neighbors[cellNo].size(); ++k ) {
        neibNo = neighbors[cellNo][k];
        dt = computeDt(nodes[neibNo], Rx, cellNo);
        if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
            traveltime =  nodes[neibNo].getTT(threadNo)+dt;
        }
    }
    return traveltime;
}


template<typename T1, typename T2, typename NODE>
T1 Grid2Dri<T1,T2,NODE>::getTraveltime(const sxz<T1>& Rx,
								  const std::vector<NODE>& nodes,
								  T2& nodeParentRx, T2& cellParentRx,
								  const size_t threadNo) const {
    
    for ( size_t nn=0; nn<nodes.size(); ++nn ) {
        if ( nodes[nn] == Rx ) {
            nodeParentRx = nodes[nn].getNodeParent(threadNo);
            cellParentRx = nodes[nn].getCellParent(threadNo);
            return nodes[nn].getTT(threadNo);
        }
    }
    
    T2 cellNo = getCellNo( Rx );
    T2 neibNo = neighbors[cellNo][0];
    T1 dt = computeDt(nodes[neibNo], Rx, cellNo);
    
    T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
    nodeParentRx = neibNo;
    cellParentRx = cellNo;
    for ( size_t k=1; k< neighbors[cellNo].size(); ++k ) {
        neibNo = neighbors[cellNo][k];
        dt = computeDt(nodes[neibNo], Rx, cellNo);
        if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
            traveltime =  nodes[neibNo].getTT(threadNo)+dt;
            nodeParentRx = neibNo;
        }
    }
    return traveltime;
}

template<typename T1, typename T2, typename NODE>
void Grid2Dri<T1,T2,NODE>::saveTT(const std::string& fname, const int all,
                                  const size_t nt,
                                  const bool vtkFormat) const {
    
    if (vtkFormat) {
#ifdef VTK

        std::string filename = fname+".vtr";
        int nn[3] = {static_cast<int>(nCellx+1), 1, static_cast<int>(nCellz+1)};
        
        vtkSmartPointer<vtkDoubleArray> xCoords = vtkSmartPointer<vtkDoubleArray>::New();
        for (size_t n=0; n<nn[0]; ++n)
            xCoords->InsertNextValue( xmin + n*dx );
        vtkSmartPointer<vtkDoubleArray> yCoords = vtkSmartPointer<vtkDoubleArray>::New();
        yCoords->InsertNextValue( 0.0 );
        vtkSmartPointer<vtkDoubleArray> zCoords = vtkSmartPointer<vtkDoubleArray>::New();
        for (size_t n=0; n<nn[2]; ++n)
            zCoords->InsertNextValue( zmin + n*dz );
        
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
                vtkIdType id = rgrid->FindPoint(nodes[n].getX(), 0.0, nodes[n].getZ());
                newScalars->SetTuple1(id, nodes[n].getTT(nt) );
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
    } else {
        std::ofstream fout(fname.c_str());
        fout.precision(12);
        for ( T2 n=0; n<nodes.size(); ++n ) {
            fout << nodes[n].getX() << '\t'
            << nodes[n].getZ() << '\t'
            << nodes[n].getTT(nt) << '\n';
        }
        fout.close();
    }
}

template<typename T1, typename T2, typename NODE>
void Grid2Dri<T1,T2,NODE>::grad(sxz<T1>& g, const size_t i, const size_t j,
                                const size_t nt) const {
    // compute average gradient for cell (i,j)
    
    static const size_t nnz = nCellz+1;
    
    g.x = 0.5*(( nodes[(i+1)*nnz+j].getTT(nt)+nodes[(i+1)*nnz+j+1].getTT(nt) ) -
               ( nodes[    i*nnz+j].getTT(nt)+nodes[    i*nnz+j+1].getTT(nt) ))/dx;
    g.z = 0.5*(( nodes[i*nnz+j+1].getTT(nt)+nodes[(i+1)*nnz+j+1].getTT(nt) ) -
               ( nodes[i*nnz+j  ].getTT(nt)+nodes[(i+1)*nnz+j  ].getTT(nt) ))/dz;
}

template<typename T1, typename T2, typename NODE>
void Grid2Dri<T1,T2,NODE>::getRaypath(const std::vector<sxz<T1>>& Tx,
                                      const sxz<T1> &Rx,
                                      std::vector<sxz<T1>> &r_data,
                                      const size_t threadNo) const {
    
    r_data.push_back( Rx );
    
    for ( size_t ns=0; ns<Tx.size(); ++ns ) {
        if ( Rx == Tx[ns] ) {
            return;
        }
    }
    
    std::vector<bool> txOnNode( Tx.size(), false );
    std::vector<T2> txNode( Tx.size() );
    std::vector<T2> txCell( Tx.size() );
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
        }
    }

	long long int iIn, jIn, iOut, jOut; // iOut, jOut for cell we are exiting; iIn, jIn for cell we are entering
    sxz<T1> curr_pt( Rx );
    sxz<T1> gOut;
    
    // distance between opposite nodes of a cell
    static const T1 maxDist = sqrt( dx*dx + dz*dz );
    
    bool onNode=false;
    if ( fabs(remainder(curr_pt.x,dx))<small && fabs(remainder(curr_pt.z,dz))<small ) {
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
		
		if ( iIn<0 || iIn>nCellx || jIn<0 || jIn>nCellz ) {
			//  we are going oustide the grid!
			std::cerr << "Error while computing raypaths: going outside grid!\n"
			<< "  Stopping calculations, raypaths will be incomplete.\n" << std::endl;
			return;
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
        if ( fabs(remainder(curr_pt.x,dx))<small && fabs(remainder(curr_pt.z,dz))<small ) {
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
            if ( i<=nCellx && j<=nCellz )
                cells.push_back( {i, j} );
            if ( i>0 && j<=nCellz )
                cells.push_back( {i-1, j} );
            if ( i<=nCellx && j>0 )
                cells.push_back( {i, j-1} );
            if ( i>0 && j>0 )
                cells.push_back( {i-1, j-1} );
            
            gOut = static_cast<T1>(0.0);
            sxz<T1> g;
            T2 nc;
            for ( nc=0; nc<cells.size(); ++nc ) {
                grad( g, cells[nc].i, cells[nc].j, threadNo );
                g *= -1.0;
                gOut += g;
            }
            gOut /= nc;  // ag holds average grad
            
            
            if ((gOut.x<0.0 && i==0) || (gOut.x>0.0 && i==nCellx+1) ||
                (gOut.z<0.0 && j==0) || (gOut.z>0.0 && j==nCellz+1)) {
                //  we are going oustide the grid!
                std::cerr << "Error while computing raypaths: going outside grid!\n"
                          << "  Stopping calculations, raypaths will be incomplete.\n" << std::endl;
                return;
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
            
            sxz<T1> gIn;
            grad( gIn, iIn, jIn, threadNo );
            gIn *= -1.0;
            T1 theta = atan2( gIn.z, gIn.x );

            if ( iIn == iOut && jIn == jOut ) {
                // ray is returning to cell it was exiting
                // we might be at grazing incidence
                // check if gIn is significantly different from gOut
                
                T1 thetaOut = atan2( gOut.z, gOut.x );
                if ( fabs(theta-thetaOut) < pi/180. ) {
                    // less that 1 degree difference, OK
                    gOut = gIn;
                } else {
                    std::cerr << "Error while computing raypaths: raypath not converging!\n"
                    << "  Stopping calculations, raypaths will be incomplete.\n" << std::endl;
                    return;
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
            if ( fabs(remainder(curr_pt.x,dx))<small && fabs(remainder(curr_pt.z,dz))<small ) {
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

#endif

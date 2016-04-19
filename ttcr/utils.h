//
//  utils.h
//  ttcr_u
//
//  Created by Bernard Giroux on 2014-02-15.
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

#ifndef ttcr_u_utils_h
#define ttcr_u_utils_h

#include <iostream>
#include <set>
#include <vector>

#ifdef VTK
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyLine.h"
#include "vtkSmartPointer.h"
#include "vtkXMLPolyDataWriter.h"
#endif

#include "MSHReader.h"
#include "Rcv.h"

template<typename T>
constexpr T factorial(T n)
{
    return n <= 1 ? 1 : (n * factorial(n-1));
}

template<typename T>
T triangleArea2D(T x1, T y1, T x2, T y2, T x3, T y3) {
	return (x1-x2)*(y2-y3) - (x2-x3)*(y1-y2);
}

template<typename T>
void buildReflectors(const MSHReader &reader,
					 const std::vector<sxyz<T>> &nodes,
					 const size_t nsrc,
					 const int nsecondary,
					 std::vector<Rcv<T>> &reflectors) {
	
	std::vector<std::string> reflector_names = reader.getPhysicalNames(2);
	std::vector<int> indices = reader.getPhysicalIndices(2);
	
	if ( reflector_names.size() != indices.size() ) {
        std::cerr << "Error - definition of reflectors\n";
		exit(1);
	}
	
	std::vector<triangleElem<uint32_t>> triangles;
	reader.readTriangleElements(triangles);
	
	sxyz<T> pt1, pt2, pt3, pt4, pt5, pt6, d1, d2;
	int ncut = nsecondary - 1;
	
	for ( size_t ni=0; ni<indices.size(); ++ni ) {
		
		reflectors.push_back( reflector_names[ni] );
		
		std::set<sxyz<T>> refl_pts;  // use set to avoid duplicate points
		typename std::set<sxyz<T>>::iterator it;
		
		for ( size_t nt=0; nt<triangles.size(); ++nt ) {
			if ( indices[ni] == triangles[nt].physical_entity ) {
				pt1 = nodes[ triangles[nt].i[0] ];
				pt2 = nodes[ triangles[nt].i[1] ];
				pt3 = nodes[ triangles[nt].i[2] ];
				
				// edge nodes
				// 1st edge
				d1.x = (pt2.x-pt1.x)/(nsecondary+1);
				d1.y = (pt2.y-pt1.y)/(nsecondary+1);
				d1.z = (pt2.z-pt1.z)/(nsecondary+1);
				
				refl_pts.insert( pt1 );
				for ( size_t n2=0; n2<nsecondary; ++n2 ) {
					pt4.x = pt1.x+(1+n2)*d1.x;
					pt4.y = pt1.y+(1+n2)*d1.y;
					pt4.z = pt1.z+(1+n2)*d1.z;
					refl_pts.insert( pt4 );
				}
				refl_pts.insert( pt2 );
				
				// 2nd edge
				d1.x = (pt3.x-pt2.x)/(nsecondary+1);
				d1.y = (pt3.y-pt2.y)/(nsecondary+1);
				d1.z = (pt3.z-pt2.z)/(nsecondary+1);
				
				for ( size_t n2=0; n2<nsecondary; ++n2 ) {
					pt4.x = pt2.x+(1+n2)*d1.x;
					pt4.y = pt2.y+(1+n2)*d1.y;
					pt4.z = pt2.z+(1+n2)*d1.z;
					refl_pts.insert( pt4 );
				}
				refl_pts.insert( pt3 );
				
				// 3rd edge
				d1.x = (pt1.x-pt3.x)/(nsecondary+1);
				d1.y = (pt1.y-pt3.y)/(nsecondary+1);
				d1.z = (pt1.z-pt3.z)/(nsecondary+1);
				
				for ( size_t n2=0; n2<nsecondary; ++n2 ) {
					pt4.x = pt3.x+(1+n2)*d1.x;
					pt4.y = pt3.y+(1+n2)*d1.y;
					pt4.z = pt3.z+(1+n2)*d1.z;
					refl_pts.insert( pt4 );
				}
				
				// face nodes
				d2.x = (pt1.x-pt2.x)/(nsecondary+1);
				d2.y = (pt1.y-pt2.y)/(nsecondary+1);
				d2.z = (pt1.z-pt2.z)/(nsecondary+1);
				
				for ( size_t n=0; n<ncut; ++n ) {
					
					pt4.x = pt3.x+(1+n)*d1.x;
					pt4.y = pt3.y+(1+n)*d1.y;
					pt4.z = pt3.z+(1+n)*d1.z;
					
					pt5.x = pt2.x+(1+n)*d2.x;
					pt5.y = pt2.y+(1+n)*d2.y;
					pt5.z = pt2.z+(1+n)*d2.z;
					
					size_t nseg = ncut+1-n;
					
					sxyz<T> d = { (pt5.x-pt4.x)/nseg,
						(pt5.y-pt4.y)/nseg,
						(pt5.z-pt4.z)/nseg };
					
					for ( size_t n2=0; n2<nseg-1; ++n2 ) {
						pt6.x = pt1.x+(1+n2)*d.x;
						pt6.y = pt1.y+(1+n2)*d.y;
						pt6.z = pt1.z+(1+n2)*d.z;
						refl_pts.insert( pt6 );
					}
				}
			}
		}
		for (it=refl_pts.begin(); it!=refl_pts.end(); ++it) {
			reflectors.back().add_coord( *it );
		}
		reflectors.back().init_tt( nsrc );
	}
}

template<typename T>
void saveRayPaths(const std::string &fname,
				  const std::vector<std::vector<sxyz<T>>> &r_data) {
	
#ifdef VTK
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	
	size_t npts=0;
	for ( size_t n=0; n<r_data.size(); ++n ) {
		npts += r_data[n].size();
	}
	pts->SetNumberOfPoints(npts);
	for ( size_t n=0, npts=0; n<r_data.size(); ++n ) {
		for ( size_t np=0; np<r_data[n].size(); ++np, ++npts ) {
			pts->InsertPoint(npts, r_data[n][np].x, r_data[n][np].y, r_data[n][np].z);
		}
	}
	polydata->SetPoints(pts);
    
	for ( size_t n=0, npts=0; n<r_data.size(); ++n ) {
		vtkSmartPointer<vtkPolyLine> line = vtkSmartPointer<vtkPolyLine>::New();
		line->GetPointIds()->SetNumberOfIds( r_data[n].size() );
		for ( size_t np=0; np<r_data[n].size(); ++np, ++npts ) {
			line->GetPointIds()->SetId(np, npts);
		}
		cellarray->InsertNextCell(line);
	}
	polydata->SetLines(cellarray);
    
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName( fname.c_str() );
    writer->SetInputData( polydata );
    writer->SetDataModeToBinary();
    writer->Update();
#endif
    
}

template<typename T>
void saveRayPaths(const std::string &fname,
				  const std::vector<std::vector<sxz<T>>> &r_data) {
	
#ifdef VTK
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	
	size_t npts=0;
	for ( size_t n=0; n<r_data.size(); ++n ) {
		npts += r_data[n].size();
	}
	pts->SetNumberOfPoints(npts);
	for ( size_t n=0, npts=0; n<r_data.size(); ++n ) {
		for ( size_t np=0; np<r_data[n].size(); ++np, ++npts ) {
			pts->InsertPoint(npts, r_data[n][np].x, 0.0, r_data[n][np].z);
		}
	}
	polydata->SetPoints(pts);
    
	for ( size_t n=0, npts=0; n<r_data.size(); ++n ) {
		vtkSmartPointer<vtkPolyLine> line = vtkSmartPointer<vtkPolyLine>::New();
		line->GetPointIds()->SetNumberOfIds( r_data[n].size() );
		for ( size_t np=0; np<r_data[n].size(); ++np, ++npts ) {
			line->GetPointIds()->SetId(np, npts);
		}
		cellarray->InsertNextCell(line);
	}
	polydata->SetLines(cellarray);
    
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName( fname.c_str() );
    writer->SetInputData( polydata );
    writer->SetDataModeToBinary();
    writer->Update();
#endif
    
}

template<typename T>
std::string to_string( const T & value )
{
    // utiliser un flux de sortie pour créer la chaîne
    std::ostringstream oss;
    // écrire la valeur dans le flux
    oss << value;
    // renvoyer une string
    return oss.str();
}


#endif

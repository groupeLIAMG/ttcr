//
//  utils.h
//  ttcr
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

#ifndef ttcr_utils_h
#define ttcr_utils_h

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

namespace ttcr {
    
    template<typename T>
    constexpr T factorial(T n)
    {
        return n <= 1 ? 1 : (n * factorial(n-1));
    }
    
    template<typename T>
    bool areCoplanar(const sxyz<T> &x1, const sxyz<T> &x2,
                     const sxyz<T> &x3, const sxyz<T> &x4) {
        return (std::abs( dot( x3-x1, cross(x2-x1, x4-x3) ) )<small2);
    }
    
    template<typename T, typename NODE>
    bool areCoplanar(const sxyz<T> &x1, const NODE &x2,
                     const NODE &x3, const NODE &x4) {
        return (std::abs( dot( x3-x1, cross(x2-x1, x4-x3) ) )<small2);
    }
    
    template<typename T, typename NODE>
    bool areCollinear(const sxyz<T> &pt, const NODE &n0, const NODE &n1) {
        
        // http://mathworld.wolfram.com/Collinear.html
        //
        sxyz<T> v = cross(pt-n0, pt-n1);
        return norm(v)<small2;
        
    }
    
    template<typename T>
    T triangleArea2D(T x1, T y1, T x2, T y2, T x3, T y3) {
        return (x1-x2)*(y2-y3) - (x2-x3)*(y1-y2);
    }
    
    template<typename T, typename NODE>
    void barycentric(const NODE *a,
                     const NODE *b,
                     const NODE *c,
                     const sxyz<T> &p,
                     T &u, T &v, T &w) {
        
        sxyz<T> ab = {b->getX()-a->getX(), b->getY()-a->getY(), b->getZ()-a->getZ()};
        sxyz<T> ac = {c->getX()-a->getX(), c->getY()-a->getY(), c->getZ()-a->getZ()};
        
        // Unnormalized triangle normal
        sxyz<T> m = cross(ab, ac);
        
        // Nominators and one-over-denominator for u and v ratios
        T nu, nv, ood;
        
        // Absolute components for determining projection plane
        T x = std::abs(m.x), y = std::abs(m.y), z = std::abs(m.z);
        
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

    template<typename T, typename NODE>
    bool testInTriangleBoundingBox(const NODE *vertexA,
                                   const NODE *vertexB,
                                   const NODE *vertexC,
                                   const sxyz<T> &E) {
        T xMin = vertexA->getX() < vertexB->getX() ? vertexA->getX() : vertexB->getX();
        xMin = xMin < vertexC->getX() ? xMin : vertexC->getX();
        T xMax = vertexA->getX() > vertexB->getX() ? vertexA->getX() : vertexB->getX();
        xMax = xMax > vertexC->getX() ? xMax : vertexC->getX();
        
        T yMin = vertexA->getY() < vertexB->getY() ? vertexA->getY() : vertexB->getY();
        yMin = yMin < vertexC->getY() ? yMin : vertexC->getY();
        T yMax = vertexA->getY() > vertexB->getY() ? vertexA->getY() : vertexB->getY();
        yMax = yMax > vertexC->getY() ? yMax : vertexC->getY();
        
        T zMin = vertexA->getZ() < vertexB->getZ() ? vertexA->getZ() : vertexB->getZ();
        zMin = zMin < vertexC->getZ() ? zMin : vertexC->getZ();
        T zMax = vertexA->getZ() > vertexB->getZ() ? vertexA->getZ() : vertexB->getZ();
        zMax = zMax > vertexC->getZ() ? zMax : vertexC->getZ();
        
        if ( E.x < xMin || xMax < E.x || E.y < yMin || yMax < E.y || E.z < zMin || zMax < E.z )
            return false;
        else
            return true;
    }
    
    template<typename T, typename NODE>
    T distSqPointToSegment(const NODE *vertexA,
                           const NODE *vertexB,
                           const sxyz<T> &E) {
        
        T p1_p2_squareLength = (vertexB->getX() - vertexA->getX())*(vertexB->getX() - vertexA->getX()) +
        (vertexB->getY() - vertexA->getY())*(vertexB->getY() - vertexA->getY()) +
        (vertexB->getZ() - vertexA->getZ())*(vertexB->getZ() - vertexA->getZ());
        
        T dotProd = ((E.x - vertexA->getX())*(vertexB->getX() - vertexA->getX()) +
                     (E.y - vertexA->getY())*(vertexB->getY() - vertexA->getY()) +
                     (E.z - vertexA->getZ())*(vertexB->getZ() - vertexA->getZ())) / p1_p2_squareLength;
        
        if ( dotProd < 0.0 ) {
            return (E.x-vertexA->getX())*(E.x-vertexA->getX()) +
            (E.y-vertexA->getY())*(E.y-vertexA->getY()) +
            (E.z-vertexA->getZ())*(E.z-vertexA->getZ());
        } else if ( dotProd > 1.0 ) {
            return (E.x-vertexB->getX())*(E.x-vertexB->getX()) +
            (E.y-vertexB->getY())*(E.y-vertexB->getY()) +
            (E.z-vertexB->getZ())*(E.z-vertexB->getZ());
        } else {
            T p_p1_squareLength = (vertexA->getX() - E.x)*(vertexA->getX() - E.x) +
            (vertexA->getY() - E.y)*(vertexA->getY() - E.y) +
            (vertexA->getZ() - E.z)*(vertexA->getZ() - E.z);
            return p_p1_squareLength - dotProd * dotProd * p1_p2_squareLength;
        }
    }
    
    template<typename T, typename NODE>
    bool testInTriangle(const NODE *vertexA,
                        const NODE *vertexB,
                        const NODE *vertexC,
                        const sxyz<T> &E) {
        
        if ( ! testInTriangleBoundingBox(vertexA, vertexB, vertexC, E) )
            return false;
        
        T u, v, w;
        barycentric(vertexA, vertexB, vertexC, E, u, v, w);
        if ( v >= 0.0 && w >= 0.0 && (v + w) <= 1.0 )
            return true;
        
        if (distSqPointToSegment(vertexA, vertexB, E) <= small2)
            return true;
        if (distSqPointToSegment(vertexA, vertexC, E) <= small2)
            return true;
        if (distSqPointToSegment(vertexB, vertexC, E) <= small2)
            return true;
        
        return false;
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
    
}

#endif

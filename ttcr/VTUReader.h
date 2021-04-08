//
//  VTUReader.h
//  ttcr
//
//  Created by Bernard Giroux on 2013-01-10.
//  Copyright (c) 2013 Bernard Giroux. All rights reserved.
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

#ifndef ttcr_VTUReader_h
#define ttcr_VTUReader_h

#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkIdList.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSmartPointer.h"
#include "vtkTetra.h"
#include "vtkTriangle.h"
#include "vtkXMLUnstructuredGridReader.h"

#include "ttcr_t.h"

namespace ttcr {

    class VTUReader {
    public:
        VTUReader(const char *fname) : filename(fname), valid(false), nNodes(0),
        nElements(0) {
            valid = check_format();
        }

        bool isValid() const { return valid; }

        int get2Ddim() const {
            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(filename.c_str());
            reader->Update();
            double x[3];
            double xmin=0.0;
            double xmax=0.0;
            double ymin=0.0;
            double ymax=0.0;
            double zmin=0.0;
            double zmax=0.0;

            reader->GetOutput()->GetPoint(0, x);
            xmin = xmax = x[0];
            ymin = ymax = x[1];
            zmin = zmax = x[2];
            for ( size_t n=1; n<reader->GetOutput()->GetNumberOfPoints(); ++n ) {
                reader->GetOutput()->GetPoint(n, x);
                xmin = xmin<x[0] ? xmin : x[0];
                xmax = xmax>x[0] ? xmax : x[0];
                ymin = ymin<x[1] ? ymin : x[1];
                ymax = ymax>x[1] ? ymax : x[1];
                zmin = zmin<x[2] ? zmin : x[2];
                zmax = zmax>x[2] ? zmax : x[2];
            }
            if ( xmin == xmax ) {
                throw std::runtime_error("Error: mesh should vary in X");
            }
            if ( ymin == ymax && zmin == zmax) {
                throw std::runtime_error("Error: mesh is 1D");
            }
            if ( ymin == ymax ) {
                return 2;
            } else if ( zmin == zmax ) {
                return 1;
            }
            return 0;
        }


        size_t getNumberOfElements() {
            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(filename.c_str());
            reader->Update();
            return reader->GetOutput()->GetNumberOfCells();
        }

        size_t getNumberOfNodes() {
            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(filename.c_str());
            reader->Update();
            return reader->GetOutput()->GetNumberOfPoints();
        }

        template<typename T>
        void readNodes2D(std::vector<sxz<T>>& nodes, const int d) {
            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(filename.c_str());
            reader->Update();

            double x[3];
            nodes.resize( reader->GetOutput()->GetNumberOfPoints() );
            for ( size_t n=0; n<reader->GetOutput()->GetNumberOfPoints(); ++n ) {
                reader->GetOutput()->GetPoint(n, x);
                nodes[n].x = x[0];
                nodes[n].z = x[d];
            }
        }

        template<typename T>
        void readNodes3D(std::vector<sxyz<T>>& nodes) {
            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(filename.c_str());
            reader->Update();

            double x[3];
            nodes.resize( reader->GetOutput()->GetNumberOfPoints() );
            for ( size_t n=0; n<reader->GetOutput()->GetNumberOfPoints(); ++n ) {
                reader->GetOutput()->GetPoint(n, x);
                nodes[n].x = x[0];
                nodes[n].y = x[1];
                nodes[n].z = x[2];
            }
        }

        template<typename T>
        void readTriangleElements(std::vector<triangleElem<T>>& tri) {
            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(filename.c_str());
            reader->Update();

            vtkSmartPointer<vtkIdList> list = vtkSmartPointer<vtkIdList>::New();
            tri.resize( reader->GetOutput()->GetNumberOfCells() );
            for ( size_t n=0; n<reader->GetOutput()->GetNumberOfCells(); ++n ) {
                if ( reader->GetOutput()->GetCell(n)->GetCellType() != VTK_TRIANGLE ) {
                    std::cerr << "Error: VTK file should only contain cells of type triangle\n";
                    std::abort();
                }
                reader->GetOutput()->GetCellPoints(n, list);
                tri[n].i[0] = static_cast<T>( list->GetId( 0 ) );
                tri[n].i[1] = static_cast<T>( list->GetId( 1 ) );
                tri[n].i[2] = static_cast<T>( list->GetId( 2 ) );
            }
        }

        template<typename T>
        void readTetrahedronElements(std::vector<tetrahedronElem<T>>& tet) {
            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(filename.c_str());
            reader->Update();

            vtkSmartPointer<vtkIdList> list = vtkSmartPointer<vtkIdList>::New();
            tet.resize( reader->GetOutput()->GetNumberOfCells() );
            for ( size_t n=0; n<reader->GetOutput()->GetNumberOfCells(); ++n ) {
                if ( reader->GetOutput()->GetCell(n)->GetCellType() != VTK_TETRA ) {
                    std::cerr << "Error: VTK file should only contain cells of type tetrahedron\n";
                    std::abort();
                }
                reader->GetOutput()->GetCellPoints(n, list);
                tet[n].i[0] = static_cast<T>( list->GetId( 0 ) );
                tet[n].i[1] = static_cast<T>( list->GetId( 1 ) );
                tet[n].i[2] = static_cast<T>( list->GetId( 2 ) );
                tet[n].i[3] = static_cast<T>( list->GetId( 3 ) );
            }
        }

        template<typename T>
        int readSlowness(std::vector<T>& slowness, const bool constCells=true) {
            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(filename.c_str());
            reader->Update();

            if ( constCells ) {
                if ( reader->GetOutput()->GetCellData()->HasArray("Slowness") == 0 &&
                    reader->GetOutput()->GetCellData()->HasArray("Velocity") == 0 ) {
                    std::cerr << "No Slowness data in file " << filename << std::endl;
                    return 0;
                }

                if ( reader->GetOutput()->GetCellData()->HasArray("Slowness") == 1 ) {

                    vtkSmartPointer<vtkCellData> cd = vtkSmartPointer<vtkCellData>::New();
                    cd = reader->GetOutput()->GetCellData();
                    vtkSmartPointer<vtkDoubleArray> slo = vtkSmartPointer<vtkDoubleArray>::New();
                    slo = vtkDoubleArray::SafeDownCast( cd->GetArray("Slowness") );

                    if ( slo->GetSize() != reader->GetOutput()->GetNumberOfCells() ) {
                        std::cerr << "Problem with Slowness data (wrong size)" << std::endl;
                        return 0;
                    }

                    slowness.resize( slo->GetSize() );
                    for ( size_t n=0; n<slo->GetSize(); ++n ) {
                        slowness[n] = slo->GetComponent(n, 0);
                    }
                } else {
                    vtkSmartPointer<vtkCellData> cd = vtkSmartPointer<vtkCellData>::New();
                    cd = reader->GetOutput()->GetCellData();
                    vtkSmartPointer<vtkDoubleArray> vel = vtkSmartPointer<vtkDoubleArray>::New();
                    vel = vtkDoubleArray::SafeDownCast( cd->GetArray("Velocity") );

                    if ( vel->GetSize() != reader->GetOutput()->GetNumberOfCells() ) {
                        std::cerr << "Problem with Slowness data (wrong size)" << std::endl;
                        return 0;
                    }

                    slowness.resize( vel->GetSize() );
                    for ( size_t n=0; n<vel->GetSize(); ++n ) {
                        slowness[n] = 1./vel->GetComponent(n, 0);
                    }
                }
            } else {
                if ( reader->GetOutput()->GetPointData()->HasArray( "Slowness") == 0 &&
                    reader->GetOutput()->GetPointData()->HasArray("Velocity") == 0 ) {
                    std::cerr << "No Slowness data in file " << filename << std::endl;
                    return 0;
                }

                if ( reader->GetOutput()->GetPointData()->HasArray( "Slowness") == 1 ) {
                    vtkSmartPointer<vtkPointData> pd = vtkSmartPointer<vtkPointData>::New();
                    pd = reader->GetOutput()->GetPointData();
                    vtkSmartPointer<vtkDoubleArray> slo = vtkSmartPointer<vtkDoubleArray>::New();
                    slo = vtkDoubleArray::SafeDownCast( pd->GetArray("Slowness") );

                    slowness.resize( slo->GetSize() );
                    for ( size_t n=0; n<slo->GetSize(); ++n ) {
                        slowness[n] = slo->GetComponent(n, 0);
                    }
                } else {
                    vtkSmartPointer<vtkPointData> pd = vtkSmartPointer<vtkPointData>::New();
                    pd = reader->GetOutput()->GetPointData();
                    vtkSmartPointer<vtkDoubleArray> vel = vtkSmartPointer<vtkDoubleArray>::New();
                    vel = vtkDoubleArray::SafeDownCast( pd->GetArray("Velocity") );

                    slowness.resize( vel->GetSize() );
                    for ( size_t n=0; n<vel->GetSize(); ++n ) {
                        slowness[n] = 1./vel->GetComponent(n, 0);
                    }

                }
            }

            return 1;
        }

        bool isConstCell() const {
            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(filename.c_str());
            reader->Update();
            return reader->GetOutput()->GetCellData()->HasArray("Slowness") == 1 ||
            reader->GetOutput()->GetCellData()->HasArray("Velocity") == 1;
        }


    private:
        std::string filename;
        bool valid;
        size_t nNodes;
        size_t nElements;

        bool check_format() const {

            bool constCells = isConstCell();

            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(filename.c_str());
            reader->Update();

            if ( reader->GetOutput() ) {

                if ( constCells ) { // slowness defined at cells

                    if ( reader->GetOutput()->GetCellData()->HasArray("Slowness") == 0 &&
                        reader->GetOutput()->GetCellData()->HasArray("Velocity") == 0 ) {
                        std::cerr << "No Slowness data in file " << filename << std::endl;
                        return false;
                    }

                    if ( reader->GetOutput()->GetCellData()->HasArray("Slowness") == 1 ) {
                        if ( reader->GetOutput()->GetCellData()->GetArray("Slowness")->GetSize() != reader->GetOutput()->GetNumberOfCells() ) {
                            std::cerr << "Problem with Slowness data (wrong size)" << std::endl;
                            return false;
                        }
                    } else {
                        if ( reader->GetOutput()->GetCellData()->GetArray("Velocity")->GetSize() != reader->GetOutput()->GetNumberOfCells() ) {
                            std::cerr << "Problem with Slowness data (wrong size)" << std::endl;
                            return false;
                        }
                    }

                    return true;

                } else {  // slowness defined at grid nodes

                    if ( reader->GetOutput()->GetPointData()->HasArray( "Slowness") == 0 &&
                        reader->GetOutput()->GetPointData()->HasArray( "Velocity") == 0 ) {
                        std::cerr << "No Slowness data in file " << filename << std::endl;
                        return false;
                    }

                    if ( reader->GetOutput()->GetPointData()->HasArray("Slowness") == 1 ) {
                        if ( reader->GetOutput()->GetPointData()->GetArray("Slowness")->GetSize() != reader->GetOutput()->GetNumberOfPoints() ) {
                            std::cerr << "Problem with Slowness data (wrong size)" << std::endl;
                            return false;
                        }
                    } else {
                        if ( reader->GetOutput()->GetPointData()->GetArray("Velocity")->GetSize() != reader->GetOutput()->GetNumberOfPoints() ) {
                            std::cerr << "Problem with Slowness data (wrong size)" << std::endl;
                            return false;
                        }
                    }

                    return true;
                }
            }
            return false;
        }
    };

}

#endif

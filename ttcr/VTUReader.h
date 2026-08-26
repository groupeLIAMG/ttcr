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

/**
 * @file VTUReader.h
 * @brief Reader for VTK XML unstructured grid (@c .vtu) mesh files.
 *
 * Declares ttcr::VTUReader, which extracts nodes, elements and model properties
 * from a @c .vtu file so an unstructured ttcr grid can be built from it. It is
 * the VTK counterpart of MSHReader.h, and carries more than geometry: a @c .vtu
 * also supplies the slowness or velocity field.
 *
 * @section vtu_scanning Access model
 * As with ttcr::MSHReader, the object holds only a filename and **every
 * accessor constructs a fresh reader and re-reads the file**. Read each
 * quantity once rather than calling these in a loop.
 *
 * @section vtu_slowness Slowness and velocity
 * The model property may be attached either to cells (constant slowness per
 * cell, the @c uc grid classes) or to points (slowness at the nodes, the @c un
 * classes); ttcr::VTUReader::isConstCell reports which. Either way it may be
 * stored as an array named @c "Slowness" or one named @c "Velocity" — when it
 * is @c "Velocity", ttcr::VTUReader::readSlowness inverts it, so that accessor
 * always yields slowness.
 *
 * @note This header requires VTK unconditionally — it includes the VTK headers
 *       outside any @c \#ifdef, unlike Src.h and Rcv.h whose VTK use is guarded.
 *       It can only be included in a build linked against VTK.
 *
 * @sa MSHReader.h, Grid3Dun.h, Grid3Duc.h
 */

#ifndef ttcr_VTUReader_h
#define ttcr_VTUReader_h

#include <iostream>

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

    /**
     * @brief Reads VTK XML unstructured grid (@c .vtu) files.
     *
     * Each accessor re-reads the file; see @ref vtu_scanning for the access
     * model and @ref vtu_slowness for how the model property is stored.
     */
    class VTUReader {
    public:
        /**
         * @brief Open a @c .vtu file and validate it.
         * @param[in] fname path to the file.
         * @note Does not throw on a missing or malformed file — check
         *       @ref isValid before using any accessor.
         */
        VTUReader(const char *fname) : filename(fname), valid(false), nNodes(0),
        nElements(0) {
            valid = check_format();
        }

        /**
         * @brief Whether the file is readable and carries a usable model property.
         * @return True if it holds a @c "Slowness" or @c "Velocity" array, on
         *         cells or on points, whose length matches the mesh.
         */
        bool isValid() const { return valid; }

        /**
         * @brief Identify which coordinate plays the role of depth in a planar mesh.
         * @return 2 if the mesh is flat in y (so z is the second coordinate),
         *         1 if it is flat in z (so y is), 0 if it is genuinely 3-D.
         * @throws std::runtime_error if the mesh does not vary in x, or if it is
         *         degenerate in both y and z.
         * @note The return value is the @c d argument expected by
         *       @ref readNodes2D.
         */
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


        /// @return Number of cells (triangles or tetrahedra) in the mesh.
        size_t getNumberOfElements() const {
            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(filename.c_str());
            reader->Update();
            return reader->GetOutput()->GetNumberOfCells();
        }

        /// @return Number of points (nodes) in the mesh.
        size_t getNumberOfNodes() const {
            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(filename.c_str());
            reader->Update();
            return reader->GetOutput()->GetNumberOfPoints();
        }

        /**
         * @brief Read the node coordinates of a planar mesh into the x-z plane.
         * @tparam T floating-point type of the output coordinates.
         * @param[out] nodes  resized to the point count and filled; @c x is taken
         *                    from the file's x, and @c z from coordinate @p d.
         * @param[in] d       which file coordinate supplies the depth: 1 for y,
         *                    2 for z. This is what @ref get2Ddim returns.
         * @pre @p d is 1 or 2.
         */
        template<typename T>
        void readNodes2D(std::vector<sxz<T>>& nodes, const int d) const {
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

        /**
         * @brief Read all node coordinates.
         * @tparam T floating-point type of the output coordinates.
         * @param[out] nodes resized to the point count and filled with (x, y, z).
         */
        template<typename T>
        void readNodes3D(std::vector<sxyz<T>>& nodes) const {
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

        /**
         * @brief Read the mesh cells as triangles.
         * @tparam T integer type used for vertex indices.
         * @param[out] tri resized to the cell count and filled with each cell's
         *                 three 0-based vertex indices.
         * @warning Every cell must be a @c VTK_TRIANGLE. On encountering any
         *          other cell type this reports to @c std::cerr and calls
         *          @c std::abort() — it does not throw, so the caller gets no
         *          chance to recover from a mixed-type mesh.
         */
        template<typename T>
        void readTriangleElements(std::vector<triangleElem<T>>& tri) const {
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

        /**
         * @brief Read the mesh cells as tetrahedra.
         * @tparam T integer type used for vertex indices.
         * @param[out] tet resized to the cell count and filled with each cell's
         *                 four 0-based vertex indices.
         * @warning Every cell must be a @c VTK_TETRA. On encountering any other
         *          cell type this reports to @c std::cerr and calls
         *          @c std::abort() — it does not throw.
         */
        template<typename T>
        void readTetrahedronElements(std::vector<tetrahedronElem<T>>& tet) const {
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
        
        /**
         * @brief Test whether a named data array is present.
         * @tparam T unused; present only to match the other accessors' shape,
         *           so callers must still write @c hasVariable<double>(...).
         * @param[in] name       array name to look for.
         * @param[in] constCells true to search the cell data, false the point data.
         * @return True if an array of that name exists in the chosen location.
         */
        template<typename T>
        bool hasVariable(const std::string& name, const bool constCells=true) const {
            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(filename.c_str());
            reader->Update();
            if ( constCells ) {
                return reader->GetOutput()->GetCellData()->HasArray(name.c_str()) == 1;
            } else {
                return reader->GetOutput()->GetPointData()->HasArray(name.c_str()) == 1;
            }
        }

        /**
         * @brief Read the model property, always returning slowness.
         *
         * Accepts either a @c "Slowness" or a @c "Velocity" array and prefers
         * @c "Slowness" when both are present; a velocity array is inverted
         * element-wise, so the output is slowness either way. See
         * @ref vtu_slowness.
         *
         * @tparam T floating-point type of the output.
         * @param[out] slowness   resized and filled with one value per cell (or
         *                        per point when @p constCells is false).
         * @param[in] constCells  true to read the cell data, false the point data.
         * @return 1 on success, 0 if neither array is present or the cell-data
         *         array has the wrong length.
         * @warning Reports failure through the **return value**, not an
         *          exception — check it. On failure @p slowness is left untouched.
         * @warning A zero velocity yields an infinite slowness; the inversion is
         *          not guarded.
         */
        template<typename T>
        int readSlowness(std::vector<T>& slowness, const bool constCells=true) const {
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

        /**
         * @brief Read an arbitrary named data array.
         *
         * The general-purpose counterpart of @ref readSlowness, used for the
         * extra fields an anisotropic model needs (Thomsen parameters, tilt
         * angles, and so on). No inversion or other interpretation is applied.
         *
         * @tparam T floating-point type of the output.
         * @param[in] name        array name to read.
         * @param[out] var        resized and filled with one value per cell (or
         *                        per point when @p constCells is false).
         * @param[in] constCells  true to read the cell data, false the point data.
         * @return 1 on success, 0 if the array is missing or, for cell data, has
         *         the wrong length.
         * @warning Reports failure through the **return value**, not an exception.
         */
        template<typename T>
        int readVariable(const std::string& name, std::vector<T>& var, const bool constCells=true) const {
            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(filename.c_str());
            reader->Update();

            if ( constCells ) {
                if ( reader->GetOutput()->GetCellData()->HasArray(name.c_str()) == 0 ) {
                    std::cerr << "No " << name << " data in file " << filename << std::endl;
                    return 0;
                }

                vtkSmartPointer<vtkCellData> cd = vtkSmartPointer<vtkCellData>::New();
                cd = reader->GetOutput()->GetCellData();
                vtkSmartPointer<vtkDoubleArray> slo = vtkSmartPointer<vtkDoubleArray>::New();
                slo = vtkDoubleArray::SafeDownCast( cd->GetArray(name.c_str()) );
                
                if ( slo->GetSize() != reader->GetOutput()->GetNumberOfCells() ) {
                    std::cerr << "Problem with data (wrong size)" << std::endl;
                    return 0;
                }
                
                var.resize( slo->GetSize() );
                for ( size_t n=0; n<slo->GetSize(); ++n ) {
                    var[n] = slo->GetComponent(n, 0);
                }

            } else {
                if ( reader->GetOutput()->GetPointData()->HasArray(name.c_str()) == 0 ) {
                    std::cerr << "No " << name << " data in file " << filename << std::endl;
                    return 0;
                }

                vtkSmartPointer<vtkPointData> pd = vtkSmartPointer<vtkPointData>::New();
                pd = reader->GetOutput()->GetPointData();
                vtkSmartPointer<vtkDoubleArray> slo = vtkSmartPointer<vtkDoubleArray>::New();
                slo = vtkDoubleArray::SafeDownCast( pd->GetArray(name.c_str()) );
                
                var.resize( slo->GetSize() );
                for ( size_t n=0; n<slo->GetSize(); ++n ) {
                    var[n] = slo->GetComponent(n, 0);
                }
            }

            return 1;
        }

        /**
         * @brief Whether the model property is attached to cells or to points.
         * @return True if a @c "Slowness" or @c "Velocity" array is present in
         *         the **cell** data, meaning slowness is constant within each
         *         cell; false otherwise, meaning it is defined at the nodes.
         * @note This is what selects between the @c uc and @c un families of
         *       grid classes, and it is also the value to pass as the
         *       @c constCells argument of @ref readSlowness.
         */
        bool isConstCell() const {
            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(filename.c_str());
            reader->Update();
            return reader->GetOutput()->GetCellData()->HasArray("Slowness") == 1 ||
            reader->GetOutput()->GetCellData()->HasArray("Velocity") == 1 ||
            reader->GetOutput()->GetCellData()->HasArray("Vp0") == 1 ||
            reader->GetOutput()->GetCellData()->HasArray("Vs0") == 1;
        }


    private:
        std::string filename;  ///< Path to the @c .vtu file; re-read by every accessor.
        bool valid;            ///< Result of the format check made at construction.
        size_t nNodes;         ///< Unused; node counts come from @ref getNumberOfNodes.
        size_t nElements;      ///< Unused; cell counts come from @ref getNumberOfElements.

        /**
         * @brief Check that the file is readable and carries a usable model property.
         *
         * Requires a @c "Slowness" or @c "Velocity" array — on cells or points,
         * whichever @ref isConstCell reports — of a length matching the mesh.
         *
         * @return True if all of that holds.
         * @note Diagnostics go to @c std::cerr; the caller only sees the bool.
         */
        bool check_format() const {

            bool constCells = isConstCell();

            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(filename.c_str());
            reader->Update();

            if ( reader->GetOutput() ) {

                if ( constCells ) { // slowness defined at cells

                    // a transversely isotropic medium carries axial velocities
                    // instead of a slowness; the cell class reads them itself
                    if ( reader->GetOutput()->GetCellData()->HasArray("Vp0") == 1 ||
                        reader->GetOutput()->GetCellData()->HasArray("Vs0") == 1 ) {
                        return true;
                    }

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

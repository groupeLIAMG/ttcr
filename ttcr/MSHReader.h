//
//  MSHReader.h
//  ttcr
//
//  Created by Bernard Giroux on 2012-09-19.
//  Copyright (c) 2012 Bernard Giroux. All rights reserved.
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
 * @file MSHReader.h
 * @brief Reader for Gmsh's native ASCII "MSH" mesh format.
 *
 * Declares ttcr::MSHReader, which extracts nodes, elements and physical-entity
 * names from a Gmsh @c .msh file so an unstructured ttcr grid can be built from
 * it. VTUReader.h is the equivalent for VTK unstructured grids.
 *
 * @section msh_scanning Access model
 * The reader holds only a filename: it keeps no parsed mesh in memory, and
 * **every accessor reopens the file and scans it from the beginning** for the
 * section it needs. That keeps the object tiny and makes each call independent,
 * but it means the cost is a full file scan per call — so read each quantity
 * once into your own container rather than calling these in a loop. The
 * exception is the physical-name lookup, which is cached on first use.
 *
 * @section msh_indexing Index base
 * Gmsh numbers nodes and physical entities from 1; this reader converts to
 * **0-based** indices throughout, so element vertex indices and physical-entity
 * indices are directly usable as C++ subscripts.
 *
 * @note Only MSH file format version 2.2, file type 0 (ASCII), is accepted;
 *       ttcr::MSHReader::isValid reports false for anything else, including the
 *       newer 4.x formats written by current Gmsh releases.
 *
 * @sa VTUReader.h, structs_msh2vtk.h
 */

#ifndef ttcr_MSHReader_h
#define ttcr_MSHReader_h

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "ttcr_t.h"

namespace ttcr {

    /**
     * @brief Reads Gmsh's native "MSH" ASCII format.
     *
     * Each accessor rescans the file; see @ref msh_scanning and
     * @ref msh_indexing for the access model and index convention.
     */
    class MSHReader {
    public:
        /**
         * @brief Open a mesh file and validate its format.
         * @param[in] fname path to the @c .msh file.
         * @note Does not throw on a missing or malformed file — check
         *       @ref isValid before using any accessor.
         */
        MSHReader(const char *fname) : filename(fname), valid(false),
        physicalNames(std::vector<std::vector<std::string>>(4)),
        physicalIndices(std::vector<std::vector<size_t>>(4)){
            valid = checkFormat();
        }

        /// @return True if the file exists and declares MSH version 2.2, ASCII.
        bool isValid() const { return valid; }

        /**
         * @brief Point the reader at a different file, discarding cached state.
         * @param[in] fname path to the new @c .msh file.
         * @post The physical-name cache is cleared and @ref isValid is
         *       recomputed for the new file.
         */
        void setFilename(const char *fname) {  // we reset the reader
            filename = fname;
            valid = checkFormat();
            for ( auto it=physicalNames.begin(); it!=physicalNames.end(); ++it ) {
                it->clear();
            }
            for ( auto it = physicalIndices.begin(); it!=physicalIndices.end(); ++it ) {
                it->clear();
            }
        }

        /**
         * @brief Test whether the mesh is planar, and hence usable as a 2-D model.
         * @return True if the nodes are constant in y or in z.
         * @throws std::runtime_error if the mesh is degenerate in both y and z,
         *         i.e. one-dimensional.
         */
        bool is2D() const {
            std::vector<sxyz<double>> nodes;
            readNodes3D(nodes);
            double ymin=0.0;
            double ymax=0.0;
            double zmin=0.0;
            double zmax=0.0;
            if ( nodes.size()> 1 ) {
                ymin = ymax = nodes[0].y;
                zmin = zmax = nodes[0].z;
            }
            for ( size_t n=1; n<nodes.size(); ++n ) {
                ymin = ymin<nodes[n].y ? ymin : nodes[n].y;
                ymax = ymax>nodes[n].y ? ymax : nodes[n].y;
                zmin = zmin<nodes[n].z ? zmin : nodes[n].z;
                zmax = zmax>nodes[n].z ? zmax : nodes[n].z;
            }
            if ( ymin == ymax && zmin == zmax) {
                throw std::runtime_error("Error: mesh is 1D");
            }
            return ymin == ymax || zmin == zmax;
        }

        /**
         * @brief Identify which coordinate plays the role of depth in a planar mesh.
         * @return 2 if the mesh is flat in y (so z is the second coordinate),
         *         1 if it is flat in z (so y is), 0 if it is genuinely 3-D.
         * @throws std::runtime_error if the mesh does not vary in x, or if it is
         *         degenerate in both y and z.
         * @note The return value is the @c d argument expected by
         *       @ref readNodes2D.
         */
        size_t get2Ddim() const {
            std::vector<sxyz<double>> nodes;
            readNodes3D(nodes);
            double xmin=0.0;
            double xmax=0.0;
            double ymin=0.0;
            double ymax=0.0;
            double zmin=0.0;
            double zmax=0.0;
            if ( nodes.size()> 1 ) {
                xmin = xmax = nodes[0].x;
                ymin = ymax = nodes[0].y;
                zmin = zmax = nodes[0].z;
            }
            for ( size_t n=1; n<nodes.size(); ++n ) {
                xmin = xmin<nodes[n].x ? xmin : nodes[n].x;
                xmax = xmax>nodes[n].x ? xmax : nodes[n].x;
                ymin = ymin<nodes[n].y ? ymin : nodes[n].y;
                ymax = ymax>nodes[n].y ? ymax : nodes[n].y;
                zmin = zmin<nodes[n].z ? zmin : nodes[n].z;
                zmax = zmax>nodes[n].z ? zmax : nodes[n].z;
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

        /**
         * @brief Total element count declared in the @c $Elements section.
         * @return Number of elements of every type combined.
         * @note This is the raw declared total. Use @ref getNumberOfLines,
         *       @ref getNumberOfTriangles or @ref getNumberOfTetra for the count
         *       of one element type.
         */
        size_t getNumberOfElements() const {
            std::ifstream fin(filename.c_str());
            std::string line;
            size_t nElements=0;
            while ( fin ) {
                getline( fin, line );
                if ( line.find("$Elements") != std::string::npos ) {
                    fin >> nElements;
                    break;
                }
            }
            fin.close();
            return nElements;
        }

        /**
         * @brief Node count declared in the @c $Nodes section.
         * @return Number of nodes in the mesh.
         */
        size_t getNumberOfNodes() const {
            std::ifstream fin(filename.c_str());
            std::string line;
            size_t nNodes=0;
            while ( fin ) {
                getline( fin, line );
                if ( line.find("$Nodes") != std::string::npos ) {
                    fin >> nNodes;
                    break;
                }
            }
            fin.close();
            return nNodes;
        }

        /**
         * @brief Names of the physical entities of a given dimension.
         *
         * Physical entities are the named groups Gmsh attaches to parts of a
         * mesh; ttcr uses them to identify reflectors and material regions.
         *
         * @param[in] i entity dimension: 0 points, 1 curves, 2 surfaces,
         *              3 volumes (the default).
         * @return Names, parallel to @ref getPhysicalIndices for the same @p i.
         * @pre @p i is at most 3 — it indexes a 4-element table unchecked.
         * @note Read on first use and cached, so repeated calls are cheap. This
         *       is the one accessor that does not rescan the file every time.
         */
        const std::vector<std::string>& getPhysicalNames(size_t i=3) const {
            if ( physicalNames[i].empty() ) {
                readPhysicalNames(i);
            }
            return physicalNames[i];
        }

        /**
         * @brief Indices of the physical entities of a given dimension.
         * @param[in] i entity dimension, as for @ref getPhysicalNames.
         * @return Indices, parallel to @ref getPhysicalNames for the same @p i.
         *         **0-based**, unlike the 1-based values in the file.
         * @pre @p i is at most 3 — it indexes a 4-element table unchecked.
         */
        const std::vector<size_t>& getPhysicalIndices(size_t i=3) const {
            if ( physicalIndices[i].empty() ) {
                readPhysicalNames(i);
            }
            return physicalIndices[i];
        }

        /// @return Number of 2-node line elements (Gmsh element type 1).
        size_t getNumberOfLines() const {
            return getNumberOfElements(1);
        }
        /// @return Number of 3-node triangles (Gmsh element type 2).
        size_t getNumberOfTriangles() const {
            return getNumberOfElements(2);
        }
        /// @return Number of 4-node tetrahedra (Gmsh element type 4).
        size_t getNumberOfTetra() const {
            return getNumberOfElements(4);
        }



        /**
         * @brief Read the node coordinates of a planar mesh into the x-z plane.
         * @tparam T floating-point type of the output coordinates.
         * @param[out] nodes  resized to the node count and filled; @c x is taken
         *                    from the file's x, and @c z from coordinate @p d.
         * @param[in] d       which file coordinate supplies the depth: 1 for y,
         *                    2 for z. This is what @ref get2Ddim returns.
         * @pre @p d is 1 or 2.
         */
        template<typename T>
        void readNodes2D(std::vector<sxz<T>>& nodes, const size_t d) const {
            std::ifstream fin(filename.c_str());
            std::string line;
            size_t nNodes;
            T x[3];
            while ( fin ) {
                getline( fin, line );
                if ( line.find("$Nodes") != std::string::npos ) {
                    fin >> nNodes;
                    if ( nodes.size() != nNodes ) {
                        nodes.resize( nNodes );
                    }
                    size_t index;
                    for ( size_t n=0; n<nNodes; ++n ) {
                        fin >> index >> nodes[n].x >> x[1] >> x[2];
                        nodes[n].z = x[d];
                    }
                    break;
                }
            }
            fin.close();
        }

        /**
         * @brief Read all node coordinates.
         * @tparam T floating-point type of the output coordinates.
         * @param[out] nodes resized to the node count and filled with (x, y, z).
         */
        template<typename T>
        void readNodes3D(std::vector<sxyz<T>>& nodes) const {
            std::ifstream fin(filename.c_str());
            std::string line;
            size_t nNodes;
            while ( fin ) {
                getline( fin, line );
                if ( line.find("$Nodes") != std::string::npos ) {
                    fin >> nNodes;
                    if ( nodes.size() != nNodes ) {
                        nodes.resize( nNodes );
                    }
                    size_t index;
                    for ( size_t n=0; n<nNodes; ++n ) {
                        fin >> index >> nodes[n].x >> nodes[n].y >> nodes[n].z;
                    }
                    break;
                }
            }
            fin.close();
        }

        /**
         * @brief Read the 2-node line elements (Gmsh element type 1).
         * @tparam T integer type used for vertex indices.
         * @param[out] lineElem  filled with one entry per line element and then
         *                       shrunk to the number actually found; each entry
         *                       carries its two **0-based** vertex indices and
         *                       its **0-based** physical entity tag.
         * @note Elements of other types in the same section are skipped.
         */
        template<typename T>
        void readLineElements(std::vector<lineElem<T>>& lineElem) const {
            std::ifstream fin(filename.c_str());
            std::string line;
            size_t nElements;
            while ( fin ) {
                getline( fin, line );
                if ( line.find("$Elements") != std::string::npos ) {
                    fin >> nElements;
                    if ( lineElem.size() != nElements ) {
                        lineElem.resize( nElements );
                    }
                    size_t index;
                    size_t elm_type, nTags;
                    std::vector<T> tags(10);
                    size_t nLines=0;
                    T ind;
                    for ( size_t n=0; n<nElements; ++n ) {
                        fin >> index >> elm_type >> nTags;
                        if ( tags.size() < nTags ) { tags.resize( nTags ); }
                        if ( elm_type == 1 ) {
                            for ( size_t n=0; n<nTags; ++n ) {
                                fin >> tags[n];
                            }
                            fin >> ind;
                            lineElem[nLines].i[0] = ind-1;
                            fin >> ind;
                            lineElem[nLines].i[1] = ind-1;
                            lineElem[nLines++].physical_entity = tags[0]-1;
                        } else {
                            getline( fin, line );
                        }
                    }
                    lineElem.resize( nLines );
                    break;
                }
            }
            fin.close();
        }

        /**
         * @brief Read the 3-node triangles (Gmsh element type 2).
         * @tparam T integer type used for vertex indices.
         * @param[out] tri  filled with one entry per triangle and then shrunk to
         *                  the number actually found; each entry carries its
         *                  three **0-based** vertex indices and its **0-based**
         *                  physical entity tag.
         * @note Elements of other types in the same section are skipped.
         */
        template<typename T>
        void readTriangleElements(std::vector<triangleElem<T>>& tri) const {
            std::ifstream fin(filename.c_str());
            std::string line;
            size_t nElements;
            while ( fin ) {
                getline( fin, line );
                if ( line.find("$Elements") != std::string::npos ) {
                    fin >> nElements;
                    if ( tri.size() != nElements ) {
                        tri.resize( nElements );
                    }
                    size_t index;
                    size_t elm_type, nTags;
                    std::vector<T> tags(10);
                    size_t nTriangles=0;
                    T ind;
                    for ( size_t n=0; n<nElements; ++n ) {
                        fin >> index >> elm_type >> nTags;
                        if ( tags.size() < nTags ) { tags.resize( nTags ); }
                        if ( elm_type == 2 ) {
                            for ( size_t n=0; n<nTags; ++n ) {
                                fin >> tags[n];
                            }
                            fin >> ind;
                            tri[nTriangles].i[0] = ind-1;
                            fin >> ind;
                            tri[nTriangles].i[1] = ind-1;
                            fin >> ind;
                            tri[nTriangles].i[2] = ind-1;
                            tri[nTriangles++].physical_entity = tags[0]-1;
                        } else {
                            getline( fin, line );
                        }
                    }
                    tri.resize( nTriangles );
                    break;
                }
            }
            fin.close();
        }

        /**
         * @brief Read the 4-node tetrahedra (Gmsh element type 4).
         * @tparam T integer type used for vertex indices.
         * @param[out] tet  filled with one entry per tetrahedron and then shrunk
         *                  to the number actually found; each entry carries its
         *                  four **0-based** vertex indices and its **0-based**
         *                  physical entity tag.
         * @note Elements of other types in the same section are skipped.
         */
        template<typename T>
        void readTetrahedronElements(std::vector<tetrahedronElem<T>>& tet) const {
            std::ifstream fin(filename.c_str());
            std::string line;
            size_t nElements;
            while ( fin ) {
                getline( fin, line );
                if ( line.find("$Elements") != std::string::npos ) {
                    fin >> nElements;
                    if ( tet.size() != nElements ) {
                        tet.resize( nElements );
                    }
                    size_t index;
                    size_t elm_type, nTags;
                    std::vector<T> tags(10);
                    size_t nTetrahedron=0;
                    T ind;
                    for ( size_t n=0; n<nElements; ++n ) {
                        fin >> index >> elm_type >> nTags;
                        if ( tags.size() < nTags ) { tags.resize( nTags ); }
                        if ( elm_type == 4 ) {
                            for ( size_t n=0; n<nTags; ++n ) {
                                fin >> tags[n];
                            }
                            fin >> ind;
                            tet[nTetrahedron].i[0] = ind-1;
                            fin >> ind;
                            tet[nTetrahedron].i[1] = ind-1;
                            fin >> ind;
                            tet[nTetrahedron].i[2] = ind-1;
                            fin >> ind;
                            tet[nTetrahedron].i[3] = ind-1;
                            tet[nTetrahedron++].physical_entity = tags[0]-1;
                        } else {
                            getline( fin, line );
                        }
                    }
                    tet.resize( nTetrahedron );
                    break;
                }
            }
            fin.close();
        }

        /**
         * @brief Read the MSH format version from the @c $MeshFormat section.
         * @return Declared version, or 0.0 if the section is absent or unreadable.
         */
        double getVersion() const {
            std::ifstream fin(filename.c_str());
            std::string line;
            double version = 0.0;
            size_t file_type;
            while ( fin ) {
                getline( fin, line );
                if ( line.find("$MeshFormat") != std::string::npos ) {
                    fin >> version >> file_type;
                    break;
                }
            }
            return version;
        }

    private:
        std::string filename;  ///< Path to the @c .msh file; reopened by every accessor.
        bool valid;            ///< Result of the format check made at construction.
        /// Cached physical-entity names, indexed by dimension 0-3. @c mutable so
        /// the const accessors can fill it lazily.
        mutable std::vector<std::vector<std::string>> physicalNames;
        /// Cached physical-entity indices, 0-based, parallel to @ref physicalNames.
        mutable std::vector<std::vector<size_t>> physicalIndices;

        /**
         * @brief Check that the file declares MSH version 2.2, file type 0.
         * @return True only for that exact combination.
         * @note Catches and reports stream failures rather than propagating
         *       them, returning false — so a missing file reads as "invalid".
         */
        bool checkFormat() const {
            bool format_ok = false;
            std::ifstream fin;
            fin.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
            try {
                fin.open(filename.c_str());
                std::string line;
                while ( fin ) {
                    getline( fin, line );
                    if ( line.find("$MeshFormat") != std::string::npos ) {
                        double version;
                        size_t file_type;
                        fin >> version >> file_type;
                        if ( version == 2.2 && file_type == 0 ) {
                            format_ok = true;
                            break;
                        }
                    }
                }
                fin.close();
            }
            catch (std::ifstream::failure &e) {
                std::cerr << "Exception opening/reading/closing file " << filename << std::endl;
                return false;
            }
            return format_ok;
        }

        /**
         * @brief Populate the physical-name cache for one dimension.
         * @param[in] dim entity dimension to read.
         * @post @ref physicalNames and @ref physicalIndices at @p dim hold the
         *       entities of that dimension, names stripped of their quotes and
         *       indices converted to 0-based.
         */
        void readPhysicalNames(const size_t dim) const {
            std::ifstream fin(filename.c_str());
            std::string line;

            while ( fin ) {
                getline( fin, line );
                if ( line.find("$PhysicalNames") != std::string::npos ) {
                    std::string name;
                    size_t np, dimension;
                    size_t index;
                    fin >> np;
                    for ( size_t n=0; n<np; ++n ) {
                        fin >> dimension >> index;
                        getline( fin, name );
                        if ( dimension==dim ) {
                            size_t p1 = name.find("\"")+1;
                            size_t p2 = name.rfind("\"");
                            physicalNames[dimension].push_back( name.substr(p1, p2-p1) );
                            physicalIndices[dimension].push_back(index-1);
                        }
                    }
                    break;
                }
            }
            fin.close();
        }

        /**
         * @brief Count the elements of one Gmsh element type.
         * @param[in] type Gmsh element type code: 1 line, 2 triangle, 4 tetrahedron.
         * @return Number of elements of that type.
         */
        size_t getNumberOfElements(const size_t type) const {

            std::ifstream fin(filename.c_str());
            std::string line;
            size_t tmp;
            size_t nElem=0;
            while ( fin ) {
                getline( fin, line );
                if ( line.find("$Elements") != std::string::npos ) {
                    fin >> tmp;
                    size_t index;
                    size_t elm_type;
                    for ( size_t n=0; n<tmp; ++n ) {
                        fin >> index >> elm_type;
                        if ( elm_type == type ) {
                            nElem++;
                        }
                        getline( fin, line );
                    }
                    break;
                }
            }
            fin.close();
            return nElem;
        }
    };

}

#endif

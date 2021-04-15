//
//  MSHReader.h
//  ttcr3du
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

#ifndef ttcr_MSHReader_h
#define ttcr_MSHReader_h

#include <fstream>
#include <string>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "ttcr_t.h"

namespace ttcr {

    // A class to read Gmsh's native "MSH" ASCII format
    class MSHReader {
    public:
        MSHReader(const char *fname) : filename(fname), valid(false),
        physicalNames(std::vector<std::vector<std::string>>(4)),
        physicalIndices(std::vector<std::vector<int>>(4)){
            valid = checkFormat();
        }

        bool isValid() const { return valid; }

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

        int get2Ddim() const {
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

        size_t getNumberOfElements() {
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

        //
        // Return names of physical entities and their corresponding indices
        //  Note: indices start at 0, not 1 like in MSH file
        //
        const std::vector<std::string>& getPhysicalNames(size_t i=3) const {
            if ( physicalNames[i].empty() ) {
                readPhysicalNames(i);
            }
            return physicalNames[i];
        }

        const std::vector<int>& getPhysicalIndices(size_t i=3) const {
            if ( physicalIndices[i].empty() ) {
                readPhysicalNames(i);
            }
            return physicalIndices[i];
        }

        size_t getNumberOfLines() const {
            return getNumberOfElements(1);
        }
        size_t getNumberOfTriangles() const {
            return getNumberOfElements(2);
        }
        size_t getNumberOfTetra() const {
            return getNumberOfElements(4);
        }



        template<typename T>
        void readNodes2D(std::vector<sxz<T>>& nodes, const int d) const {
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
                    int elm_type, nTags;
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
                    int elm_type, nTags;
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
                    int elm_type, nTags;
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

        double getVersion() const {
            std::ifstream fin(filename.c_str());
            std::string line;
            double version = 0.0;
            int file_type;
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
        std::string filename;
        bool valid;
        mutable std::vector<std::vector<std::string>> physicalNames;
        mutable std::vector<std::vector<int>> physicalIndices;

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
                        int file_type;
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

        void readPhysicalNames(const size_t dim) const {
            std::ifstream fin(filename.c_str());
            std::string line;

            while ( fin ) {
                getline( fin, line );
                if ( line.find("$PhysicalNames") != std::string::npos ) {
                    std::string name;
                    size_t np, dimension;
                    int index;
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

        size_t getNumberOfElements(const int type) const {

            std::ifstream fin(filename.c_str());
            std::string line;
            size_t tmp;
            size_t nElem=0;
            while ( fin ) {
                getline( fin, line );
                if ( line.find("$Elements") != std::string::npos ) {
                    fin >> tmp;
                    size_t index;
                    int elm_type;
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

//  Defines a class Src that contains the position of each sources
//
//  Src.h
//  ttcr
//
//  Created by Bernard Giroux on 2012-11-20.
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
 * @file Src.h
 * @brief Source positions and origin times for a 3-D run.
 *
 * Declares ttcr::Src, which reads a source file and holds the resulting
 * coordinates and origin times. One @c Src object corresponds to one entry of
 * ttcr::input_parameters::srcfiles; a run with several shots holds a vector of
 * them.
 *
 * @section src_formats Accepted file formats
 * ttcr::Src::init recognises three layouts, chosen by sniffing the first line:
 *
 * - **Legacy VTK ASCII** — first line contains @c "vtk". The third line must say
 *   @c ASCII; the reader then scans forward to the @c POINTS record and reads
 *   that many @c "x y z" triples.
 * - **CRT** — first line ends with @c '/'. Each record is
 *   @c "name x y z /", read until the stream fails.
 * - **Plain text** (the default) — a point count on the first line, then one
 *   @c "x y z t0" row per source.
 *
 * @note Only the plain-text format carries origin times. The VTK and CRT
 *       branches set every @c t0 to zero, so use plain text whenever the shots
 *       are not all fired at t = 0.
 *
 * @sa Src2D.h, Rcv.h, structs_ttcr.h
 */

#ifndef ttcr_Src_h
#define ttcr_Src_h

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef VTK
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkXMLPolyDataWriter.h"
#endif

#include "ttcr_t.h"

namespace ttcr {

    /**
     * @brief Source coordinates and origin times read from a file.
     *
     * @tparam T floating-point type of the coordinates and origin times.
     *
     * Construction only records the filename; nothing is read until
     * @ref init is called.
     */
    template<typename T>
    class Src {
    public:
        /// @param f path to the source file; not opened until @ref init.
        Src(const std::string &f) : filename(f) {
        }

        /**
         * @brief Read the source file into @ref coord and @ref t0.
         *
         * Detects the format as described in @ref src_formats and fills both
         * vectors, which end up the same length.
         *
         * @warning Prints to @c std::cerr and calls @c exit(1) if the file
         *          cannot be opened, or if a VTK file is not ASCII. It does not
         *          throw, so a caller cannot recover from a bad source file.
         */
        void init();
        /// @return Source coordinates.
        const std::vector<sxyz<T>>& get_coord() const { return coord; }
        /// @return Writable source coordinates, e.g. to reposition sources onto the grid.
        std::vector<sxyz<T>>& get_coord() { return coord; }
        /// @return Origin time of each source, parallel to @ref get_coord. Zero unless the file was plain text.
        const std::vector<T>& get_t0() const { return t0; }

        /**
         * @brief Write the source positions as a VTK PolyData file.
         * @param[in] fname output filename.
         * @note Compiled to an empty function unless @c VTK is defined; without
         *       VTK support this silently writes nothing.
         */
        void toVTK(const std::string &fname) const;
    private:
        std::string filename;         ///< Path handed to the constructor.
        std::vector<sxyz<T>> coord;   ///< Source coordinates.
        std::vector<T> t0;            ///< Origin times, parallel to @ref coord.
    };

    template<typename T>
    void Src<T>::init() {
        std::ifstream fin;
        fin.open(filename);

        if ( !fin ) {
            std::cerr << "Cannot open file " << filename << ".\n";
            exit(1);
        }

        std::string test;
        std::getline(fin, test);

        char lastChar = '\n';
        bool vtk = false;
        if (!test.empty()) {
            lastChar = *test.rbegin();
            if ( test.find("vtk") != std::string::npos ) vtk = true;
        }
        if ( vtk == true ) {
            std::getline(fin, test);  // 2nd line should be vtk output
            std::getline(fin, test);  // 3rd line should be ASCII
            if ( test.find("ASCII") == std::string::npos ) {
                std::cerr << "Error: vtk file should be ascii.\n";
                exit(1);
            }
            while ( test.find("POINTS") == std::string::npos ) {
                std::getline(fin, test);
            }
            std::istringstream sin( test );
            size_t nsrc;
            sin >> test >> nsrc;
            coord.resize( nsrc );
            t0.resize( nsrc );
            size_t nread = 0;
            while ( fin && nread<nsrc ) {
                fin >> coord[nread].x >> coord[nread].y >> coord[nread].z;
                nread++;
            }
            for ( size_t n=0; n<t0.size(); ++n ) {
                t0[n] = 0.0;
            }

        } else if ( lastChar == '/' ) {
            // CRT format
            coord.resize(0);
            t0.resize(0);
            sxyz<T> co;
            while ( fin ) {
                fin >> test >> co.x >> co.y >> co.z;
                fin >> test;  // read / at eol
                if ( !fin.fail() ) {
                    coord.push_back( co );
                    t0.push_back( 0.0 );
                }
            }
        } else {
            fin.seekg(std::ios_base::beg);
            size_t nsrc;
            fin >> nsrc;
            coord.resize( nsrc );
            t0.resize( nsrc );
            size_t nread = 0;
            while ( fin && nread<nsrc ) {
                fin >> coord[nread].x >> coord[nread].y >> coord[nread].z >> t0[nread];
                nread++;
            }
        }
        fin.close();
    }

    template<typename T>
    void Src<T>::toVTK(const std::string &fname) const {
#ifdef VTK
        vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();

        pts->SetNumberOfPoints(coord.size());
        for ( size_t n=0; n<coord.size(); ++n ) {
            pts->InsertPoint(n, coord[n].x, coord[n].y, coord[n].z);
        }
        polydata->SetPoints(pts);

        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer->SetFileName( fname.c_str() );
        writer->SetInputData( polydata );
        writer->SetDataModeToBinary();
        writer->Update();
#endif
    }
}

#endif

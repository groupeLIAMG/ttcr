//
//  Src2D.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-01-21.
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

/**
 * @file Src2D.h
 * @brief Source positions and origin times for a 2-D (x-z) run.
 *
 * The 2-D counterpart of Src.h: ttcr::Src2D reads a source file and holds the
 * resulting ttcr::sxz coordinates and origin times.
 *
 * @section src2d_formats Accepted file formats
 * ttcr::Src2D::init recognises two layouts, chosen from the first line:
 *
 * - **CRT** — first line ends with @c '/'. Each record is @c "name x z /",
 *   read until the stream fails; origin times are set to zero.
 * - **Plain text** (the default) — a source count on the first line, then one
 *   @c "x z t0" row per source.
 *
 * @note Unlike ttcr::Src, this class has **no legacy-VTK branch** — a @c .vtk
 *       source file that works in 3-D will not be parsed here; it falls through
 *       to the plain-text branch and misreads the header.
 *
 * @sa Src.h, Rcv2D.h
 */

#ifndef ttcr_Src2D_h
#define ttcr_Src2D_h

#include <fstream>
#include <iostream>
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
     * @brief Source coordinates and origin times for a 2-D run.
     *
     * @tparam T floating-point type of the coordinates and origin times.
     *
     * Construction only records the filename; nothing is read until
     * @ref init is called.
     */
    template<typename T>
    class Src2D {
    public:
        /// @param f path to the source file; not opened until @ref init.
        Src2D(const std::string &f) : filename(f) {
        }

        /**
         * @brief Read the source file into @ref coord and @ref t0.
         *
         * Detects the format as described in @ref src2d_formats and fills both
         * vectors, which end up the same length.
         *
         * @warning Prints to @c std::cerr and calls @c exit(1) if the file
         *          cannot be opened. It does not throw.
         */
        void init();
        /// @return Source coordinates in the x-z plane.
        const std::vector<sxz<T>>& get_coord() const { return coord; }
        /// @return Origin time of each source, parallel to @ref get_coord. Zero for CRT input.
        const std::vector<T>& get_t0() const { return t0; }

        /**
         * @brief Write the source positions as a VTK PolyData file.
         * @param[in] fname output filename.
         * @note Points are emitted as (x, 0, z) — the missing y is written as zero.
         * @note Compiled to an empty function unless @c VTK is defined.
         */
        void toVTK(const std::string &fname) const;
    private:
        std::string filename;        ///< Path handed to the constructor.
        std::vector<sxz<T>> coord;   ///< Source coordinates in the x-z plane.
        std::vector<T> t0;           ///< Origin times, parallel to @ref coord.
    };

    template<typename T>
    void Src2D<T>::init() {
        std::ifstream fin;
        fin.open(filename);

        if ( !fin ) {
            std::cerr << "Cannot open file " << filename << ".\n";
            exit(1);
        }

        std::string test;
        std::getline(fin, test);

        char lastChar = '\n';
        if (!test.empty()) {
            lastChar = *test.rbegin();
        }
        if ( lastChar == '/' ) {
            // CRT format
            coord.resize(0);
            t0.resize(0);
            sxz<T> co;
            while ( fin ) {
                fin >> test >> co.x >> co.z;
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
                fin >> coord[nread].x >> coord[nread].z >> t0[nread];
                nread++;
            }
        }
        fin.close();
    }

    template<typename T>
    void Src2D<T>::toVTK(const std::string &fname) const {
#ifdef VTK
        vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();

        pts->SetNumberOfPoints(coord.size());
        for ( size_t n=0; n<coord.size(); ++n ) {
            pts->InsertPoint(n, coord[n].x, 0.0, coord[n].z);
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

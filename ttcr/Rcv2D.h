//
//  Rcv2D.h
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
 * @file Rcv2D.h
 * @brief Receiver positions and computed traveltimes for a 2-D (x-z) run.
 *
 * The 2-D counterpart of Rcv.h: ttcr::Rcv2D holds ttcr::sxz receiver coordinates
 * and the traveltimes computed at them, indexed
 * @c [source][reflector][receiver] exactly as in @ref rcv_storage.
 *
 * @section rcv2d_formats Accepted file formats
 * ttcr::Rcv2D::init recognises three layouts, chosen from the first line:
 *
 * - **Legacy VTK ASCII** — first line contains @c "vtk". Points are read as
 *   @c "x y z" triples and the **y coordinate is discarded**, so a 3-D VTK
 *   receiver file can be reused directly for a 2-D run.
 * - **CRT** — first line ends with @c '/'; records of the form @c "name x z /".
 * - **Plain text** (the default) — a receiver count, then one @c "x z" row each.
 *
 * @note ttcr::Src2D, unlike this class, has no VTK branch — the 2-D source and
 *       receiver readers do not accept the same set of formats.
 *
 * @sa Rcv.h, Src2D.h
 */

#ifndef ttcr_Rcv2D_h
#define ttcr_Rcv2D_h

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef VTK
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdocumentation-deprecated-sync"
#pragma clang diagnostic ignored "-Wimplicit-int-conversion"
#pragma clang diagnostic ignored "-Wshorten-64-to-32"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkXMLPolyDataWriter.h"
#pragma clang diagnostic pop
#endif

#include "ttcr_t.h"

namespace ttcr {

    /**
     * @brief Receiver coordinates and computed traveltimes for a 2-D run.
     *
     * @tparam T floating-point type of the coordinates and traveltimes.
     *
     * Two ways to populate it: read a file with @ref init, or build the array
     * programmatically with @ref add_coord followed by @ref init_tt.
     */
    template<typename T>
    class Rcv2D {
    public:
        /// @param f path to the receiver file; not opened until @ref init.
        Rcv2D(const std::string &f) : filename(f) {
        }

        /**
         * @brief Read the receiver file and size the traveltime storage.
         *
         * Detects the format as described in @ref rcv2d_formats, fills
         * @ref coord, and allocates traveltime room for @p ns sources and
         * @p nr reflectors.
         *
         * @param[in] ns  number of sources the traveltimes will be computed for.
         * @param[in] nr  number of reflecting interfaces; the default 0 keeps
         *                direct arrivals only.
         * @warning Prints to @c std::cerr and calls @c exit(1) if the file
         *          cannot be opened, or if a VTK file is not ASCII.
         */
        void init(const size_t ns, const size_t nr=0);
        /// @return Receiver coordinates in the x-z plane.
        const std::vector<sxz<T>>& get_coord() const { return coord; }
        /**
         * @brief Writable traveltimes for one source and one arrival, for a solver to fill.
         * @param[in] ns  source number.
         * @param[in] nr  reflector number; 0 (the default) is the direct arrival.
         * @return Traveltime at each receiver, one entry per receiver.
         * @pre @ref init or @ref init_tt has been called, and both indices are in range.
         */
        std::vector<T>& get_tt(const size_t ns, const size_t nr=0) { return tt[ns][nr]; }
        /// @copydoc get_tt(const size_t, const size_t)
        const std::vector<T>& get_tt(const size_t ns, const size_t nr=0) const {
            return tt[ns][nr];
        }

        /**
         * @brief Write the traveltimes of one source to a text file.
         *
         * One line per receiver; the direct arrival first, then one
         * tab-separated column per reflector. Written with 9 significant digits.
         *
         * @param[in] f   output file.
         * @param[in] ns  source number to save.
         * @warning Calls @c exit(1) if the file cannot be opened for writing.
         */
        void save_tt( const std::string &f, const size_t ns) const;

        /**
         * @brief Append a receiver position.
         * @param[in] c coordinates of the new receiver.
         * @note Does not resize the traveltime storage — call @ref init_tt once
         *       every receiver has been added.
         */
        void add_coord(const sxz<T> &c) { coord.push_back(c); }
        /**
         * @brief Allocate traveltime storage for direct arrivals only.
         * @param[in] nsrc number of sources.
         * @note Allocates the source and reflector levels but leaves the innermost
         *       vectors empty; the solver writing into @ref get_tt sizes those.
         */
        void init_tt(const size_t nsrc) {
            tt.resize(nsrc);
            for (size_t ns=0; ns<nsrc; ++ns ) {
                tt[ns].resize( 1 );
            }
        }

        /**
         * @brief Write the receiver coordinates back to @ref filename.
         *
         * Emits the plain-text format — a count, then one @c "x z" row per
         * receiver — with 17 significant digits, so a round trip is lossless.
         *
         * @warning Overwrites the file the receivers were read from.
         */
        void save_rcvfile() const;
        /**
         * @brief Write the receiver positions as a VTK PolyData file.
         * @param[in] fname output filename.
         * @note Points are emitted as (x, 0, z) — the missing y is written as zero.
         * @note Compiled to an empty function unless @c VTK is defined.
         */
        void toVTK(const std::string &fname) const;
    private:
        std::string filename;       ///< Path handed to the constructor; also the target of @ref save_rcvfile.
        std::vector<sxz<T>> coord;  ///< Receiver coordinates in the x-z plane.
        std::vector<std::vector<std::vector<T>>> tt;  ///< Traveltimes, indexed [source][reflector][receiver].
    };

    template<typename T>
    void Rcv2D<T>::init(const size_t nsrc, const size_t nrefl) {
        std::ifstream fin;
        fin.open(filename);

        if ( !fin ) {
            std::cerr << "Cannot open file " << filename << " for reading.\n";
            exit(1);
        }
        tt.resize(nsrc);
        for (size_t ns=0; ns<nsrc; ++ns ) {
            tt[ns].resize( nrefl+1 );
        }

        std::string test;
        std::getline(fin, test);

        char lastChar = ' ';
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
            size_t nrcv;
            float tmp;
            sin >> test >> nrcv;
            coord.resize( nrcv );
            for (size_t n=0; n<nsrc; ++n ) {
                for ( size_t nr=0; nr<=nrefl; ++nr ) {
                    tt[n][nr].resize( nrcv );
                }
            }
            size_t nread = 0;
            while ( fin && nread<nrcv ) {
                fin >> coord[nread].x >> tmp >> coord[nread].z;  // discard y coord
                nread++;
            }

        } else if ( lastChar == '/' ) {
            // CRT format
            coord.resize(0);
            sxz<T> co;
            while ( fin ) {
                fin >> test >> co.x >> co.z;
                fin >> test;  // read / at eol
                if ( !fin.fail() ) {
                    coord.push_back( co );
                }
            }
            for (size_t ns=0; ns<nsrc; ++ns ) {
                for ( size_t nr=0; nr<=nrefl; ++nr ) {
                    tt[ns][nr].resize( coord.size() );
                }
            }
        } else {
            fin.seekg(std::ios_base::beg);
            size_t nrcv;
            fin >> nrcv;
            coord.resize( nrcv );
            for (size_t ns=0; ns<nsrc; ++ns ) {
                for ( size_t nr=0; nr<=nrefl; ++nr ) {
                    tt[ns][nr].resize( nrcv );
                }
            }
            size_t nread = 0;
            while ( fin && nread<nrcv ) {
                fin >> coord[nread].x >> coord[nread].z;
                nread++;
            }
        }
        fin.close();
    }

    template<typename T>
    void Rcv2D<T>::save_tt( const std::string &filename, const size_t ns) const {

        std::ofstream fout;
        fout.open( filename );
        if ( !fout ) {
            std::cerr << "Cannot open file " << filename << " for writing.\n";
            exit(1);
        }
        size_t nrefl = tt[ns].size();
        fout.precision(9);
        for ( size_t n=0; n<tt[ns][0].size(); ++n ) {
            fout << tt[ns][0][n];
            for ( size_t nr=1; nr<nrefl; ++nr ) {
                fout << '\t' << tt[ns][nr][n];
            }
            fout << '\n';
        }
        fout.close();
    }

    template<typename T>
    void Rcv2D<T>::save_rcvfile() const {

        std::ofstream fout;
        fout.open( filename );
        if ( !fout ) {
            std::cerr << "Cannot open file " << filename << " for writing.\n";
            exit(1);
        }
        fout << coord.size() << '\n';
        fout.precision(17);
        fout << std::scientific;
        for ( size_t n=0; n<coord.size(); ++n ) {
            fout << coord[n].x << '\t' << coord[n].z << '\n';
        }
        fout.close();
    }

    template<typename T>
    void Rcv2D<T>::toVTK(const std::string &fname) const {
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

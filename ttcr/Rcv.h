//
//  Rcv.h
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
 * @file Rcv.h
 * @brief Receiver positions and the traveltimes computed at them, for a 3-D run.
 *
 * Declares ttcr::Rcv, which reads a receiver file and owns the traveltime
 * results. Unlike ttcr::Src there is one @c Rcv per run, not one per shot: the
 * same receiver array is reused for every source, and the traveltimes are
 * indexed by source.
 *
 * @section rcv_storage Traveltime storage
 * Traveltimes live in a three-level structure indexed
 * @c [source][reflector][receiver]:
 * - the outer level has one entry per source,
 * - the middle level has @c nrefl+1 entries — index 0 is the direct arrival and
 *   1..@c nrefl are arrivals reflected off each interface,
 * - the inner level has one traveltime per receiver.
 *
 * @section rcv_formats Accepted file formats
 * ttcr::Rcv::init recognises the same three layouts as ttcr::Src, chosen by
 * sniffing the first line: legacy **VTK ASCII** (first line contains @c "vtk";
 * scans to the @c POINTS record), **CRT** (first line ends with @c '/'; records
 * of the form @c "name x y z /"), and **plain text** (a receiver count, then one
 * @c "x y z" row each). No format carries traveltimes — those are outputs.
 *
 * @sa Rcv2D.h, Src.h, Interface.h
 */

#ifndef ttcr_Rcv_h
#define ttcr_Rcv_h

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef VTK
#pragma clang diagnostic push
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
     * @brief Receiver coordinates and the traveltimes computed at them.
     *
     * @tparam T floating-point type of the coordinates and traveltimes.
     *
     * Two ways to populate it: read a file with @ref init, or build the array
     * programmatically with @ref add_coord followed by @ref init_tt.
     */
    template<typename T>
    class Rcv {
    public:
        /// @param f path to the receiver file; not opened until @ref init.
        Rcv(const std::string &f) : filename(f) {
        }

        /**
         * @brief Read the receiver file and size the traveltime storage.
         *
         * Detects the format as described in @ref rcv_formats, fills
         * @ref coord, and allocates traveltime room for @p nsrc sources and
         * @p nrefl reflectors as described in @ref rcv_storage.
         *
         * @param[in] nsrc   number of sources the traveltimes will be computed for.
         * @param[in] nrefl  number of reflecting interfaces; the default 0 keeps
         *                   direct arrivals only.
         * @warning Prints to @c std::cerr and calls @c exit(1) if the file
         *          cannot be opened, or if a VTK file is not ASCII.
         */
        void init(const size_t nsrc, const size_t nrefl=0);
        /// @return Receiver coordinates.
        const std::vector<sxyz<T>>& get_coord() const { return coord; }
        /// @return Writable receiver coordinates, e.g. to project receivers onto the grid.
        std::vector<sxyz<T>>& get_coord() { return coord; }
        /**
         * @brief Writable traveltimes for one source and one arrival, for a solver to fill.
         * @param[in] n   source number.
         * @param[in] nr  reflector number; 0 (the default) is the direct arrival.
         * @return Traveltime at each receiver, one entry per receiver.
         * @pre @ref init or @ref init_tt has been called, and both indices are in range.
         */
        std::vector<T>& get_tt(const size_t n, const size_t nr=0) { return tt[n][nr]; }
        /// @copydoc get_tt(const size_t, const size_t)
        const std::vector<T>& get_tt(const size_t n, const size_t nr=0) const { return tt[n][nr]; }

        /**
         * @brief Write the traveltimes of one source to a text file.
         *
         * One line per receiver; the direct arrival first, then one
         * tab-separated column per reflector. Written with 9 significant digits.
         *
         * @param[in] filename  output file.
         * @param[in] ns        source number to save.
         * @warning Calls @c exit(1) if the file cannot be opened for writing.
         */
        void save_tt( const std::string &filename, const size_t ns) const;

        /**
         * @brief Append a receiver position.
         * @param[in] c coordinates of the new receiver.
         * @note Does not resize the traveltime storage — call @ref init_tt once
         *       every receiver has been added.
         */
        void add_coord(const sxyz<T> &c) { coord.push_back(c); }
        /**
         * @brief Allocate traveltime storage for direct arrivals only.
         * @param[in] nsrc number of sources.
         * @note Allocates the source and reflector levels but leaves the innermost
         *       vectors empty; the solver writing into @ref get_tt sizes those.
         *       Use this with @ref add_coord; @ref init does the whole job when
         *       reading from a file.
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
         * Emits the plain-text format — a count, then one @c "x y z" row per
         * receiver — with 17 significant digits, so a round trip is lossless.
         *
         * @warning Overwrites the file the receivers were read from.
         */
        void save_rcvfile() const;
        /**
         * @brief Write the receiver positions as a VTK PolyData file.
         * @param[in] fname output filename.
         * @note Compiled to an empty function unless @c VTK is defined.
         */
        void toVTK(const std::string &fname) const;
    private:
        std::string filename;        ///< Path handed to the constructor; also the target of @ref save_rcvfile.
        std::vector<sxyz<T>> coord;  ///< Receiver coordinates.
        std::vector<std::vector<std::vector<T>>> tt;  ///< Traveltimes, indexed [source][reflector][receiver]. @sa @ref rcv_storage
    };

    template<typename T>
    void Rcv<T>::init(const size_t nsrc, const size_t nrefl) {
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
            sin >> test >> nrcv;
            coord.resize( nrcv );
            for (size_t n=0; n<nsrc; ++n ) {
                for ( size_t nr=0; nr<=nrefl; ++nr ) {
                    tt[n][nr].resize( nrcv );
                }
            }
            size_t nread = 0;
            while ( fin && nread<nrcv ) {
                fin >> coord[nread].x >> coord[nread].y >> coord[nread].z;
                nread++;
            }

        } else if ( lastChar == '/' ) {
            // CRT format
            coord.resize(0);
            sxyz<T> co;
            while ( fin ) {
                fin >> test >> co.x >> co.y >> co.z;
                fin >> test;  // read / at eol
                if ( !fin.fail() ) {
                    coord.push_back( co );
                }
            }
            for (size_t n=0; n<nsrc; ++n ) {
                for ( size_t nr=0; nr<=nrefl; ++nr ) {
                    tt[n][nr].resize( coord.size() );
                }
            }
        } else {
            fin.seekg(std::ios_base::beg);
            size_t nrcv;
            fin >> nrcv;
            coord.resize( nrcv );
            for (size_t n=0; n<nsrc; ++n ) {
                for ( size_t nr=0; nr<=nrefl; ++nr ) {
                    tt[n][nr].resize( nrcv );
                }
            }
            size_t nread = 0;
            while ( fin && nread<nrcv ) {
                fin >> coord[nread].x >> coord[nread].y >> coord[nread].z;
                nread++;
            }
        }
        fin.close();
    }

    template<typename T>
    void Rcv<T>::save_tt( const std::string &filename,
                         const size_t ns) const {

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
    void Rcv<T>::save_rcvfile() const {

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
            fout << coord[n].x << '\t' << coord[n].y << '\t' << coord[n].z << '\n';
        }
        fout.close();
    }

    template<typename T>
    void Rcv<T>::toVTK(const std::string &fname) const {
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

#endif /* defined(__ttcr__Rcv__) */

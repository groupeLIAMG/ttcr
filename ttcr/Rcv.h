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

    template<typename T>
    class Rcv {
    public:
        Rcv(const std::string &f) : filename(f) {
        }

        void init(const size_t, const size_t nr=0);
        const std::vector<sxyz<T>>& get_coord() const { return coord; }
        std::vector<sxyz<T>>& get_coord() { return coord; }
        std::vector<T>& get_tt(const size_t n, const size_t nr=0) { return tt[n][nr]; }
        const std::vector<T>& get_tt(const size_t n, const size_t nr=0) const { return tt[n][nr]; }

        void save_tt( const std::string &, const size_t) const;

        void add_coord(const sxyz<T> &c) { coord.push_back(c); }
        void init_tt(const size_t nsrc) {
            tt.resize(nsrc);
            for (size_t ns=0; ns<nsrc; ++ns ) {
                tt[ns].resize( 1 );
            }
        }

        void save_rcvfile() const;
        void toVTK(const std::string &) const;
    private:
        std::string filename;
        std::vector<sxyz<T>> coord;
        std::vector<std::vector<std::vector<T>>> tt;
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

//
//  main.cpp
//  ttcr3d_raypath
//
//  Created by Bernard Giroux on 2019-07-30.
//  Copyright © 2019 Bernard Giroux. All rights reserved.
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

#include <iostream>

#include "Grid3D.h"
#include "Rcv.h"
#include "Src.h"
#include "structs_ttcr.h"
#include "ttcr_io.h"
#include "grids.h"

using namespace std;
using namespace ttcr;


// Creates a template to be able to call body() for two formats: float or double
template<typename T>

// body() is the principal function for raytracing called by main()
int body(const input_parameters &par) {


    // Create the src class containing the position of the sources
    std::vector< Src<T> > src;
    for ( size_t n=0; n<par.srcfiles.size(); ++ n ) {
        src.push_back( Src<T>( par.srcfiles[n] ) );
        string end = " ... ";
        if ( n < par.srcfiles.size() - 1 ) end = " ...\n";
        if ( verbose ) cout << "Reading source file " << par.srcfiles[n] << end;
        src[n].init();
    }
    if ( verbose ) cout << "done.\n";

    size_t num_threads = 1;

    // ? Find the generic file name of the input model?
    string::size_type idx;  // can hold a string of any length
    idx = par.modelfile.rfind('.');
    string extension = "";
    if (idx != string::npos) {
        extension = par.modelfile.substr(idx);
    }

    // Intialize a Grid3D object containing the 3D grid and the functions to perform raytracing
    Grid3D<T,uint32_t> *g=nullptr;

    vector<Rcv<T>> reflectors;

    // Load the grid file into the GRID3D object g for different formats
    if (extension == ".grd") {
        g = buildRectilinear3D<T>(par, num_threads);
    } else if (extension == ".vtr") {
#ifdef VTK
        g = buildRectilinear3DfromVtr<T>(par, num_threads);
#else
        cerr << "Error: Program not compiled with VTK support" << endl;
        return 1;
#endif
    } else if (extension == ".vtu") {
#ifdef VTK
        g = buildUnstructured3DfromVtu<T>(par, num_threads);
#else
        cerr << "Error: Program not compiled with VTK support" << endl;
        return 1;
#endif
    } else if (extension == ".msh") {
        g = buildUnstructured3D<T>(par, reflectors, num_threads, src.size());
    } else {
        cerr << par.modelfile << " Unknown extenstion: " << extension << endl;
        return 1;
    }

    if ( g == nullptr ) {
        cerr << "Error: grid cannot be built\n";
        return 1;
    }

    // Load the receiver file into the Rcv object rcv
    Rcv<T> rcv( par.rcvfile );
    if ( par.rcvfile != "" ) {
        if ( verbose ) cout << "Reading receiver file " << par.rcvfile << " ... ";
        rcv.init( src.size() );
        if ( verbose ) cout << "done.\n";
    }

    if ( verbose ) {
        if ( par.singlePrecision ) {
            cout << "Calculations will be done in single precision.\n";
        } else {
            cout << "Calculations will be done in double precision.\n";
        }
    }

    chrono::high_resolution_clock::time_point begin, end;
    std::vector<std::vector<std::vector<sxyz<T>>>> r_data(src.size());
    if ( verbose ) { cout << "Computing raypaths ... "; cout.flush(); }
    if ( par.time ) { begin = chrono::high_resolution_clock::now(); }

    for ( size_t ns=0; ns<src.size(); ++ns ) {
        // load in traveltime data

        string srcname = par.srcfiles[ns];
        size_t pos = srcname.rfind("/");
        srcname.erase(0, pos+1);
        pos = srcname.rfind(".");
        size_t len = srcname.length()-pos;
        srcname.erase(pos, len);

        string filename = par.basename+"_"+srcname+"_all_tt";
        g->loadTT(filename, 0, 0, par.saveGridTT);

        r_data[ns].resize(rcv.get_coord().size());
        for ( size_t nr=0; nr<rcv.get_coord().size(); ++nr ) {
            g->getRaypath(src[ns].get_coord(), src[ns].get_t0(), rcv.get_coord()[nr], r_data[ns][nr], rcv.get_tt(ns)[nr]);
        }

    }
    if ( verbose ) { cout << "done.\n"; }
    if ( par.time ) { end = chrono::high_resolution_clock::now(); }
    if ( verbose && par.time ) {
        cout << "Computation time: "
        << chrono::duration<double>(end-begin).count() << '\n';
    }

    if ( src.size() == 1 ) {
        string filename = par.basename+"_tt2.dat";

        if ( par.rcvfile != "" ) {
            if ( verbose ) cout << "Saving traveltimes in " << filename <<  " ... ";
            rcv.save_tt(filename, 0);
            if ( verbose ) cout << "done.\n";
        }

        if ( par.saveRaypaths && par.rcvfile != "" ) {
            filename = par.basename+"_rp2.vtp";
            if ( verbose ) cout << "Saving raypaths in " << filename <<  " ... ";
            saveRayPaths(filename, r_data[0]);
            if ( verbose ) cout << "done.\n";
        }
    } else {
        for ( size_t ns=0; ns<src.size(); ++ns ) {

            string srcname = par.srcfiles[ns];
            size_t pos = srcname.rfind("/");
            srcname.erase(0, pos+1);
            pos = srcname.rfind(".");
            size_t len = srcname.length()-pos;
            srcname.erase(pos, len);

            string filename = par.basename+"_"+srcname+"_tt2.dat";

            if ( par.rcvfile != "" ) {
                if ( verbose ) cout << "Saving traveltimes in " << filename <<  " ... ";
                rcv.save_tt(filename, ns);
                if ( verbose ) cout << "done.\n";
            }

            if ( par.saveRaypaths && par.rcvfile != "" ) {
                filename = par.basename+"_"+srcname+"_rp2.vtp";
                if ( verbose ) cout << "Saving raypaths in " << filename <<  " ... ";
                saveRayPaths(filename, r_data[ns]);
                if ( verbose ) cout << "done.\n";
            }
        }
    }

    if ( verbose ) cout << "Normal termination of program.\n";
    return 0;
}
int main(int argc, char * argv[]) {

    input_parameters par;

    string fname = parse_input(argc, argv, par);

    if ( verbose ) {
        cout << "*** Program ttcr3d_raypath ***\n\n"
        << "Raytracing in 3D media.\n";
    }
    get_params(fname, par);

    if ( par.singlePrecision ) {
        return body<float>(par);
    } else {
        return body<double>(par);
    }

    return 0;
}

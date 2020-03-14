//
//  main.cpp
//  profile_threaded
//
//  Created by Bernard Giroux on 2020-03-12.
//  Copyright © 2020 Bernard Giroux. All rights reserved.
//

#include <stdlib.h>     /* strtod */

#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>

#include "Grid3D.h"
#include "Grid3Drnfs.h"

using namespace std;
using namespace ttcr;

int main(int argc, const char * argv[]) {

    vector<double> x;
    vector<double> y;
    vector<double> z;

    vector<vector<sxyz<double>>> src;
    vector<vector<sxyz<double>>> rcv(1);
    vector<vector<double>> t0;
    vector<vector<double>> tt(1);

    vector<vector<vector<sxyz<double>>>> r_data;
    vector<vector<vector<sijv<double>>>> m_data;

    chrono::high_resolution_clock::time_point begin, end;

    char data[100];

    ifstream fin("x.dat");
    fin >> data;
    while ( fin ) {
        x.push_back(strtod(data, NULL));
        fin >> data;
    }
    fin.close();
    fin.open("y.dat", ifstream::in);
    fin >> data;
    while ( fin ) {
        y.push_back(strtod(data, NULL));
        fin >> data;
    }
    fin.close();
    fin.open("z.dat", ifstream::in);
    fin >> data;
    while ( fin ) {
        z.push_back(strtod(data, NULL));
        fin >> data;
    }
    fin.close();

    double evID, ev_t, ev_x, ev_y, ev_z;
    fin.open("src.dat", ifstream::in);
    fin >> data;
    while ( fin ) {
        evID = strtod(data, NULL);
        fin >> data;
        ev_t = strtod(data, NULL);
        fin >> data;
        ev_x = strtod(data, NULL);
        fin >> data;
        ev_y = strtod(data, NULL);
        fin >> data;
        ev_z = strtod(data, NULL);
        fin >> data;
        src.push_back( vector<sxyz<double>>(1) );
        src.back().back() = {ev_x, ev_y, ev_z};
        t0.push_back( vector<double>(1) );
        t0.back().back() = ev_t;
    }
    fin.close();

    fin.open("rcv.dat", ifstream::in);
    fin >> data;
    while ( fin ) {
        ev_x = strtod(data, NULL);
        fin >> data;
        ev_y = strtod(data, NULL);
        fin >> data;
        ev_z = strtod(data, NULL);
        fin >> data;
        rcv[0].push_back( {ev_x, ev_y, ev_z} );
    }

    tt[0].resize( rcv[0].size() );
    for ( size_t n=1; n<src.size(); ++n ) {
        rcv.push_back(rcv[0]);
        tt.push_back( vector<double>(rcv[0].size()) );
    }
    fin.close();

    vector<double> s, slowness;
    fin.open("slowness.dat", ifstream::in);
    fin >> data;
    while ( fin ) {
        s.push_back(strtod(data, NULL));
        fin >> data;
    }
    fin.close();

    for ( size_t k=0; k<z.size(); ++k ) {
        for ( size_t j=0; j<y.size(); ++j ) {
            for ( size_t i=0; i<x.size(); ++i ) {
                // slowness is in 'C' order and we must pass it in 'F' order
                slowness.push_back(s[(i*y.size() + j)*z.size() + k]);
            }
        }
    }

    double dx = x[1] - x[0];
//    double dy = y[1] - y[0];
//    double dz = z[1] - z[0];
//    cout << dx << '\t' << dy << '\t' << dz << "\n\nsrc\n";
//    for ( size_t n=0; n<src.size(); ++n ) {
//        cout << src[n][0] << '\n';
//    }
//    cout << "\nrcv\n";
//    for ( size_t n=0; n<rcv[0].size(); ++n ) {
//        cout << rcv[0][n] << '\n';
//    }

    Grid3D<double, uint32_t> *grid;

    size_t nt = 8;
    bool weno = true;
    grid = new Grid3Drnfs<double, uint32_t>(x.size()-1, y.size()-1, z.size()-1,
                                      dx, x[0], y[0], z[0], 1.e-15,
                                      20, weno, 1, 0, nt);

    tt.resize(rcv.size());
    r_data.resize(src.size());
    m_data.resize(src.size());

    begin = chrono::high_resolution_clock::now();
    grid->setSlowness(slowness);
    if ( nt == 1 ) {
        for ( size_t n=0; n<src.size(); ++n ) {
            cout << "n = " << n << '\n';
            grid->raytrace(src[n], t0[n], rcv[n], tt[n], 0);
        }
    } else {
        grid->raytrace(src, t0, rcv, tt, r_data, m_data);
    }
    end = chrono::high_resolution_clock::now();
    cout << "Time to perform raytracing: "
    << chrono::duration<double>(end-begin).count() << '\n';
    delete grid;
    return 0;
}


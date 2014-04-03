//
//  Src2D.h
//  ttcr_u
//
//  Created by Bernard Giroux on 2014-01-21.
//  Copyright (c) 2014 Bernard Giroux. All rights reserved.
//

#ifndef ttcr_u_Src2D_h
#define ttcr_u_Src2D_h

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "ttcr_t.h"

template<typename T>
class Src2D {
public:
    Src2D(const std::string &f) : filename(f) {
    }
    
    void init();
    const std::vector<sxz<T>>& get_coord() const { return coord; }
    const std::vector<T>& get_t0() const { return t0; }
    
private:
    std::string filename;
    std::vector<sxz<T>> coord;
    std::vector<T> t0;
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


#endif

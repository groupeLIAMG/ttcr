//  Defines a class Src that contains the position of each sources
//
//  Src.h
//  ttcr_u
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

#ifndef ttcr_u_Src_h
#define ttcr_u_Src_h

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "ttcr_t.h"

template<typename T>
class Src {
public:
    Src(const std::string &f) : filename(f) {
    }
    
    void init();
    const std::vector<sxyz<T>>& get_coord() const { return coord; }
    const std::vector<T>& get_t0() const { return t0; }
    
private:
    std::string filename;
    std::vector<sxyz<T>> coord;
    std::vector<T> t0;
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
    if (!test.empty()) {
        lastChar = *test.rbegin();
    }
    if ( lastChar == '/' ) {
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


#endif

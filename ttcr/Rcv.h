//
//  Rcv.h
//  ttcr_u
//
//  Created by Bernard Giroux on 2012-11-20.
//  Copyright (c) 2012 Bernard Giroux. All rights reserved.
//

#ifndef __ttcr_u__Rcv__
#define __ttcr_u__Rcv__

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "ttcr_t.h"

template<typename T>
class Rcv {
public:
    Rcv(const std::string &f) : filename(f) {
    }
    
    void init(const size_t, const size_t nr=0);
    const std::vector<sxyz<T>>& get_coord() const { return coord; }
	std::vector<sxyz<T>>* get_coord_ptr() { return &coord; }
    std::vector<T>& get_tt(const size_t n, const size_t nr=0) { return tt[n][nr]; }
    
    void save_tt( const std::string &, const size_t) const;
    
    void add_coord(const sxyz<T> &c) { coord.push_back(c); }
	void init_tt(const size_t nsrc) {
		tt.resize(nsrc);
		for (size_t ns=0; ns<nsrc; ++ns ) tt[ns].resize( 1 );
	}

    void save_rcvfile() const;
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
    for (size_t ns=0; ns<nsrc; ++ns )
        tt[ns].resize( nrefl+1 );

    std::string test;
    std::getline(fin, test);
    
    char lastChar = ' ';
    if (!test.empty()) {
        lastChar = *test.rbegin();
    }
    
    if ( lastChar == '/' ) {
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
        for (size_t n=0; n<nsrc; ++n )
            for ( size_t nr=0; nr<=nrefl; ++nr )
				tt[n][nr].resize( coord.size() );
    } else {
        fin.seekg(std::ios_base::beg);
        size_t nrcv;
        fin >> nrcv;
        coord.resize( nrcv );
        for (size_t n=0; n<nsrc; ++n )
            for ( size_t nr=0; nr<=nrefl; ++nr )
				tt[n][nr].resize( nrcv );
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
        for ( size_t nr=1; nr<nrefl; ++nr )
            fout << '\t' << tt[ns][nr][n];
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
    for ( size_t n=0; n<coord.size(); ++n )
        fout << coord[n].x << '\t' << coord[n].y << '\t' << coord[n].z << '\n';
    fout.close();
}


#endif /* defined(__ttcr_u__Rcv__) */

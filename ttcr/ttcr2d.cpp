//
//  main.cpp
//  ttcr2d
//
//  Created by Bernard Giroux on 2014-04-03.
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

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <thread>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshorten-64-to-32"
#include <boost/asio/ip/host_name.hpp>
#pragma clang diagnostic pop

#include "Grid2D.h"
#include "Rcv2D.h"
#include "Src2D.h"
#include "structs_ttcr.h"
#include "ttcr_io.h"
#include "grids.h"

using namespace std;
using namespace ttcr;

template<typename T>
int body(const input_parameters &par) {

    vector< Src2D<T>> src;
    for ( size_t n=0; n<par.srcfiles.size(); ++ n ) {
        src.push_back( Src2D<T>( par.srcfiles[n] ) );
        string end = " ... ";
        if ( n < par.srcfiles.size() - 1 ) end = " ...\n";
        if ( verbose ) cout << "Reading source file " << par.srcfiles[n] << end;
        src[n].init();
    }
    if ( verbose ) cout << "done.\n";

    size_t const nTx = src.size();
    size_t num_threads = 1;

    if ( par.nt == 0 ) {
        size_t const hardware_threads = thread::hardware_concurrency();
        size_t const max_threads = (nTx+par.min_per_thread-1)/par.min_per_thread;
        num_threads = min((hardware_threads!=0?hardware_threads:2), max_threads);
    } else {
        num_threads = par.nt < nTx ? par.nt : nTx;
    }

    size_t blk_size = nTx/num_threads;
    if ( blk_size == 0 ) blk_size++;

    string::size_type idx;

    idx = par.modelfile.rfind('.');
    string extension = "";
    if (idx != string::npos) {
        extension = par.modelfile.substr(idx);
    }

    Grid2D<T,uint32_t,sxz<T>> *g=nullptr;
    vector<Rcv2D<T>> reflectors;
    if (extension == ".grd") {
        g = buildRectilinear2D<T>(par, num_threads);
    } else if (extension == ".vtr") {
#ifdef VTK
        g = buildRectilinear2DfromVtr<T>(par, num_threads);
#else
        cerr << "Error: Program not compiled with VTK support" << endl;
        return 1;
#endif
    } else if (extension == ".vtu") {
#ifdef VTK
        g = buildUnstructured2DfromVtu<T>(par, num_threads);
#else
        cerr << "Error: Program not compiled with VTK support" << endl;
        return 1;
#endif
    } else if (extension == ".msh") {
        g = buildUnstructured2D<T>(par, reflectors, num_threads, src.size());
    } else {
        cerr << par.modelfile << " Unknown extenstion: " << extension << endl;
        return 1;
    }

    if ( g == nullptr ) {
        cerr << "Error: grid cannot be built\n";
        return 1;
    }

    if ( par.method == SHORTEST_PATH and par.dump_secondary ) {
        string sec_file = par.basename+"_sec.dat";
        ofstream otmp(sec_file);
        if ( verbose ) {
            cout << "Dumping secondary node coordinates to " << sec_file << endl;
        }
        g->dump_secondary(otmp);
        otmp.close();
    }

    Rcv2D<T> rcv( par.rcvfile );
    if ( par.rcvfile != "" ) {
        if ( verbose ) cout << "Reading receiver file " << par.rcvfile << " ... ";
        rcv.init( src.size(), reflectors.size() );
        if ( verbose ) cout << "done.\n";
    }

    if ( par.saveModelVTK ) {
#ifdef VTK
        std::string filename = par.rcvfile;
        size_t i = filename.rfind(".dat");
        filename.replace(i, 4, ".vtp");

        if ( verbose ) std::cout << "Saving receiver in " << filename << " ... ";
        rcv.toVTK(filename);
        if ( verbose ) std::cout << "done.\n";

        for ( size_t n=0; n<par.srcfiles.size(); ++ n ) {
            filename = par.srcfiles[n];
            i = filename.rfind(".dat");
            filename.replace(i, 4, ".vtp");
            if ( verbose ) std::cout << "Saving source in " << filename << " ... ";
            src[n].toVTK(filename);
            if ( verbose ) std::cout << "done.\n";
        }
#else
        std::cerr << "Error: Program not compiled with VTK support" << std::endl;
        return nullptr;
#endif
    }

    if ( verbose ) {
        if ( par.singlePrecision ) {
            cout << "Calculations will be done in single precision.\n";
        } else {
            cout << "Calculations will be done in double precision.\n";
        }
    }
    if ( verbose && num_threads>1 ) {
        cout << "Calculations will be done using " << num_threads
        << " threads with " << blk_size << " shots per thread.\n";
    }

    vector<const vector<sxz<T>>*> all_rcv;
    if ( par.rcvfile != "" )
        all_rcv.push_back( &(rcv.get_coord()) );
    for ( size_t n=0; n<reflectors.size(); ++n ) {
        all_rcv.push_back( &(reflectors[n].get_coord()) );
    }


    chrono::high_resolution_clock::time_point begin, end;
    vector<vector<vector<sxz<T>>>> r_data(src.size());
    vector<vector<vector<vector<sxz<T>>>>> rfl_r_data(reflectors.size());
    for ( size_t n=0; n<reflectors.size(); ++n ) {
        rfl_r_data[n].resize( src.size() );
    }
    vector<vector<vector<vector<sxz<T>>>>> rfl2_r_data(reflectors.size());
    for ( size_t n=0; n<reflectors.size(); ++n ) {
        rfl2_r_data[n].resize( src.size() );
    }


    if ( verbose ) { cout << "Computing traveltimes ... "; cout.flush(); }
    if ( par.time ) { begin = chrono::high_resolution_clock::now(); }
    if ( par.saveRaypaths && par.rcvfile != "" ) {
        if ( num_threads == 1 ) {
            for ( size_t n=0; n<src.size(); ++n ) {

                vector<vector<T>*> all_tt;
                all_tt.push_back( &(rcv.get_tt(n)) );
                vector<vector<vector<sxz<T>>>*> all_r_data;
                all_r_data.push_back( &(r_data[n]) );

                for ( size_t nr=0; nr<reflectors.size(); ++nr ) {
                    all_tt.push_back( &(reflectors[nr].get_tt(n)) );
                    all_r_data.push_back( &(rfl_r_data[nr][n]) );
                }
                try {
                    g->raytrace(src[n].get_coord(), src[n].get_t0(), all_rcv,
                                all_tt, all_r_data);
                } catch (std::exception& e) {
                    std::cerr << e.what() << std::endl;
                    abort();
                }

                if ( par.saveGridTT>0 ) {

                    string srcname = par.srcfiles[n];
                    size_t pos = srcname.rfind("/");
                    srcname.erase(0, pos+1);
                    pos = srcname.rfind(".");
                    size_t len = srcname.length()-pos;
                    srcname.erase(pos, len);

                    string filename = par.basename+"_"+srcname+"_all_tt";
                    g->saveTT(filename, 0, 0, par.saveGridTT);
                    //                    filename = par.basename+"_"+srcname+"_tt_grad";
                    //                    g->saveTTgrad(filename, 0, true);
                    //                    filename = par.basename+"_"+srcname+"_tt_grad";
                    //                    g->saveTTgrad2(filename, 0, true);
                }

                for ( size_t nr=0; nr<reflectors.size(); ++nr ) {
                    try {
                        g->raytrace(reflectors[nr].get_coord(),
                                    reflectors[nr].get_tt(n), rcv.get_coord(),
                                    rcv.get_tt(n,nr+1), rfl2_r_data[nr][n]);
                    } catch (std::exception& e) {
                        std::cerr << e.what() << std::endl;
                        abort();
                    }
                }
            }
        } else {
            // threaded jobs

            vector<thread> threads(num_threads-1);
            size_t blk_start = 0;
            for ( size_t i=0; i<num_threads-1; ++i ) {

                size_t blk_end = blk_start + blk_size;

                threads[i]=thread( [&g,&src,&rcv,&r_data,&reflectors,&all_rcv,
                                    &rfl_r_data,&rfl2_r_data,
                                    blk_start,blk_end,i]{

                    for ( size_t n=blk_start; n<blk_end; ++n ) {

                        vector<vector<T>*> all_tt;
                        all_tt.push_back( &(rcv.get_tt(n)) );
                        vector<vector<vector<sxz<T>>>*> all_r_data;
                        all_r_data.push_back( &(r_data[n]) );
                        for ( size_t nr=0; nr<reflectors.size(); ++nr ) {
                            all_tt.push_back( &(reflectors[nr].get_tt(n)) );
                            all_r_data.push_back( &(rfl_r_data[nr][n]) );
                        }
                        try {
                            g->raytrace(src[n].get_coord(), src[n].get_t0(), all_rcv,
                                        all_tt, all_r_data, i+1);
                        } catch (std::exception& e) {
                            std::cerr << e.what() << std::endl;
                            abort();
                        }

                        for ( size_t nr=0; nr<reflectors.size(); ++nr ) {
                            try {
                                g->raytrace(reflectors[nr].get_coord(),
                                            reflectors[nr].get_tt(n), rcv.get_coord(),
                                            rcv.get_tt(n,nr+1), rfl2_r_data[nr][n], i+1);
                            } catch (std::exception& e) {
                                std::cerr << e.what() << std::endl;
                                abort();
                            }
                        }
                    }
                });

                blk_start = blk_end;
            }
            for ( size_t n=blk_start; n<nTx; ++n ) {

                vector<vector<T>*> all_tt;
                all_tt.push_back( &(rcv.get_tt(n)) );
                vector<vector<vector<sxz<T>>>*> all_r_data;
                all_r_data.push_back( &(r_data[n]) );
                for ( size_t nr=0; nr<reflectors.size(); ++nr ) {
                    all_tt.push_back( &(reflectors[nr].get_tt(n)) );
                    all_r_data.push_back( &(rfl_r_data[nr][n]) );
                }
                try {
                    g->raytrace(src[n].get_coord(), src[n].get_t0(), all_rcv,
                                all_tt, all_r_data, 0);
                } catch (std::exception& e) {
                    std::cerr << e.what() << std::endl;
                    abort();
                }

                for ( size_t nr=0; nr<reflectors.size(); ++nr ) {
                    try {
                        g->raytrace(reflectors[nr].get_coord(),
                                    reflectors[nr].get_tt(n), rcv.get_coord(),
                                    rcv.get_tt(n,nr+1), rfl2_r_data[nr][n], 0);
                    } catch (std::exception& e) {
                        std::cerr << e.what() << std::endl;
                        abort();
                    }
                }
            }

            for_each(threads.begin(),threads.end(), mem_fn(&thread::join));
        }
    } else {
        if ( num_threads == 1 ) {
            for ( size_t n=0; n<src.size(); ++n ) {

                vector<vector<T>*> all_tt;
                if ( par.rcvfile != "" )
                    all_tt.push_back( &(rcv.get_tt(n)) );
                for ( size_t nr=0; nr<reflectors.size(); ++nr ) {
                    all_tt.push_back( &(reflectors[nr].get_tt(n)) );
                }
                try {
                    g->raytrace(src[n].get_coord(), src[n].get_t0(), all_rcv,
                                all_tt);
                } catch (std::exception& e) {
                    std::cerr << e.what() << std::endl;
                    abort();
                }

                if ( par.saveGridTT>0 ) {

                    string srcname = par.srcfiles[n];
                    size_t pos = srcname.rfind("/");
                    srcname.erase(0, pos+1);
                    pos = srcname.rfind(".");
                    size_t len = srcname.length()-pos;
                    srcname.erase(pos, len);

                    string filename = par.basename+"_"+srcname+"_all_tt";
                    g->saveTT(filename, 0, 0, par.saveGridTT);
                    //                    filename = par.basename+"_"+srcname+"_tt_grad";
                    //                    g->saveTTgrad(filename, 0, true);
                    //                    filename = par.basename+"_"+srcname+"_tt_grad";
                    //                    g->saveTTgrad2(filename, 0, true);
                }

                for ( size_t nr=0; nr<reflectors.size(); ++nr ) {
                    try {
                        g->raytrace(reflectors[nr].get_coord(),
                                    reflectors[nr].get_tt(n), rcv.get_coord(),
                                    rcv.get_tt(n,nr+1));
                    } catch (std::exception& e) {
                        std::cerr << e.what() << std::endl;
                        abort();
                    }
                }
            }
        } else {
            // threaded jobs

            vector<thread> threads(num_threads-1);
            size_t blk_start = 0;
            for ( size_t i=0; i<num_threads-1; ++i ) {

                size_t blk_end = blk_start + blk_size;

                threads[i]=thread( [&par,&g,&src,&rcv,&all_rcv,&reflectors,blk_start,
                                    blk_end,i]{

                    for ( size_t n=blk_start; n<blk_end; ++n ) {

                        vector<vector<T>*> all_tt;
                        if ( par.rcvfile != "" )
                            all_tt.push_back( &(rcv.get_tt(n)) );
                        for ( size_t nr=0; nr<reflectors.size(); ++nr ) {
                            all_tt.push_back( &(reflectors[nr].get_tt(n)) );
                        }
                        try {
                            g->raytrace(src[n].get_coord(), src[n].get_t0(), all_rcv,
                                        all_tt, i+1);
                        } catch (std::exception& e) {
                            std::cerr << e.what() << std::endl;
                            abort();
                        }

                        for ( size_t nr=0; nr<reflectors.size(); ++nr ) {
                            try {
                                g->raytrace(reflectors[nr].get_coord(),
                                            reflectors[nr].get_tt(n), rcv.get_coord(),
                                            rcv.get_tt(n,nr+1), i+1);
                            } catch (std::exception& e) {
                                std::cerr << e.what() << std::endl;
                                abort();
                            }
                        }

                    }
                });

                blk_start = blk_end;
            }
            for ( size_t n=blk_start; n<nTx; ++n ) {

                vector<vector<T>*> all_tt;
                if ( par.rcvfile != "" )
                    all_tt.push_back( &(rcv.get_tt(n)) );
                for ( size_t nr=0; nr<reflectors.size(); ++nr ) {
                    all_tt.push_back( &(reflectors[nr].get_tt(n)) );
                }
                try {
                    g->raytrace(src[n].get_coord(), src[n].get_t0(), all_rcv,
                                all_tt, 0);
                } catch (std::exception& e) {
                    std::cerr << e.what() << std::endl;
                    abort();
                }

                for ( size_t nr=0; nr<reflectors.size(); ++nr ) {
                    try {
                        g->raytrace(reflectors[nr].get_coord(),
                                    reflectors[nr].get_tt(n), rcv.get_coord(),
                                    rcv.get_tt(n,nr+1), 0);
                    } catch (std::exception& e) {
                        std::cerr << e.what() << std::endl;
                        abort();
                    }
                }

            }

            for_each(threads.begin(),threads.end(), mem_fn(&thread::join));
        }
    }
    if ( par.time ) { end = chrono::high_resolution_clock::now(); }
    if ( verbose ) {
        cout << "done.\n";
        if ( par.method == FAST_SWEEPING ) {
            std::cout << g->get_niter() << " 1st order iterations ";
            if ( par.weno3==true ) std::cout << "and " << g->get_niterw() << " 3rd order iterations ";
            std::cout << "were needed with epsilon = " << par.epsilon << '\n';
        }
    }
    if ( par.time ) {
        cout.precision(12);
        cout << "Time to perform raytracing: " <<
        chrono::duration<double>(end-begin).count() << '\n';
    }

    delete g;

    if ( src.size() == 1 ) {
        string filename = par.basename+"_tt.dat";

        if ( par.rcvfile != "" ) {
            if ( verbose ) cout << "Saving traveltimes in " << filename <<  " ... ";
            rcv.save_tt(filename, 0);
            if ( verbose ) cout << "done.\n";
        }

        if ( par.saveRaypaths && par.rcvfile != "" ) {
            filename = par.basename+"_rp.vtp";
            if ( verbose ) cout << "Saving raypaths in " << filename <<  " ... ";
            saveRayPaths(filename, r_data[0]);
            if ( verbose ) cout << "done.\n";

            for ( size_t nr=0; nr<reflectors.size(); ++nr ) {

                vector<vector<sxz<T>>> r_tmp( rcv.get_coord().size() );
                for ( size_t irx=0; irx<rcv.get_coord().size(); ++irx ) {

                    sxz<T> pt1 = rfl2_r_data[nr][0][irx][0];
                    for ( size_t n=0; n<rfl_r_data[nr][0].size(); ++n ) {
                        if ( pt1 == rfl_r_data[nr][0][n].back() ) {

                            for ( size_t i=0; i<rfl_r_data[nr][0][n].size(); ++i ) {
                                r_tmp[irx].push_back( rfl_r_data[nr][0][n][i] );
                            }
                            for ( size_t i=1; i<rfl2_r_data[nr][0][irx].size(); ++i ) {
                                r_tmp[irx].push_back( rfl2_r_data[nr][0][irx][i] );
                            }
                            break;
                        }
                    }
                }
                filename = par.basename+"_rp"+to_string(nr+1)+".vtp";
                if ( verbose ) cout << "Saving raypaths of reflected waves in " << filename <<  " ... ";
                saveRayPaths(filename, r_tmp);
                if ( verbose ) cout << "done.\n";
            }

            if ( reflectors.size() > 0 ) {
                filename = par.basename+"_rp.bin";
                ofstream fout;
                //				fout.open(filename, ios::out | ios::binary);
                //				fout << r_data.size();
                //				fout.close();

                // TODO complete this...
            }
        }
    } else {
        for ( size_t ns=0; ns<src.size(); ++ns ) {

            string srcname = par.srcfiles[ns];
            size_t pos = srcname.rfind("/");
            srcname.erase(0, pos+1);
            pos = srcname.rfind(".");
            size_t len = srcname.length()-pos;
            srcname.erase(pos, len);

            string filename = par.basename+"_"+srcname+"_tt.dat";
            if ( par.rcvfile != "" ) {
                if ( verbose ) cout << "Saving traveltimes in " << filename <<  " ... ";
                rcv.save_tt(filename, ns);
                if ( verbose ) cout << "done.\n";
            }

            if ( par.saveRaypaths && par.rcvfile != "" ) {
                filename = par.basename+"_"+srcname+"_rp.vtp";
                if ( verbose ) cout << "Saving raypaths of direct waves in " << filename <<  " ... ";
                saveRayPaths(filename, r_data[ns]);
                if ( verbose ) cout << "done.\n";

                for ( size_t nr=0; nr<reflectors.size(); ++nr ) {

                    vector<vector<sxz<T>>> r_tmp( rcv.get_coord().size() );
                    for ( size_t irx=0; irx<rcv.get_coord().size(); ++irx ) {

                        sxz<T> pt1 = rfl2_r_data[nr][ns][irx][0];
                        for ( size_t n=0; n<rfl_r_data[nr][ns].size(); ++n ) {
                            if ( pt1 == rfl_r_data[nr][ns][n].back() ) {

                                for ( size_t i=0; i<rfl_r_data[nr][ns][n].size(); ++i ) {
                                    r_tmp[irx].push_back( rfl_r_data[nr][ns][n][i] );
                                }
                                for ( size_t i=1; i<rfl2_r_data[nr][ns][irx].size(); ++i ) {
                                    r_tmp[irx].push_back( rfl2_r_data[nr][ns][irx][i] );
                                }
                                break;
                            }
                        }
                    }
                    filename = par.basename+"_"+srcname+"_rp"+to_string(nr+1)+".vtp";
                    if ( verbose ) cout << "Saving raypaths of reflected waves in " << filename <<  " ... ";
                    saveRayPaths(filename, r_tmp);
                    if ( verbose ) cout << "done.\n";
                }
            }
            if ( verbose ) cout << '\n';
        }

        if ( par.saveRaypaths && reflectors.size() > 0 ) {
            string filename = par.basename+"_rp.bin";
            ofstream fout;
            fout.open(filename, ios::out | ios::binary);
            if ( !fout ) {
                std::cerr << "Cannot open file " << filename << " for writing.\n";
                exit(1);
            }

            if ( verbose ) cout << "Saving global raypath data in " << filename << " ... ";
            size_t size = r_data.size();
            fout.write((char*)&size,sizeof(size_t));
            for ( size_t n=0; n<r_data.size(); ++n ) {
                size = r_data[n].size();
                fout.write((char*)&size,sizeof(size_t));
                for ( size_t nr=0; nr<r_data[n].size(); ++nr ) {
                    size = r_data[n][nr].size();
                    fout.write((char*)&size,sizeof(size_t));
                    fout.write((char*)r_data[n][nr].data(),
                               sizeof(sxz<T>) * size);
                }
            }

            size = rfl_r_data.size();
            fout.write((char*)&size,sizeof(size_t));
            for ( size_t r=0; r<rfl_r_data.size(); ++r ) {
                size = rfl_r_data[r].size();
                fout.write((char*)&size,sizeof(size_t));
                for ( size_t n=0; n<rfl_r_data[r].size(); ++n ) {
                    size = rfl_r_data[r][n].size();
                    fout.write((char*)&size,sizeof(size_t));
                    for ( size_t nr=0; nr<rfl_r_data[r][n].size(); ++nr ) {
                        size = rfl_r_data[r][n][nr].size();
                        fout.write((char*)&size,sizeof(size_t));
                        fout.write((char*)rfl_r_data[r][n][nr].data(),
                                   sizeof(sxz<T>) * size);
                    }
                }
            }

            size = rfl2_r_data.size();
            fout.write((char*)&size,sizeof(size_t));
            for ( size_t r=0; r<rfl2_r_data.size(); ++r ) {
                size = rfl2_r_data[r].size();
                fout.write((char*)&size,sizeof(size_t));
                for ( size_t n=0; n<rfl2_r_data[r].size(); ++n ) {
                    size = rfl2_r_data[r][n].size();
                    fout.write((char*)&size,sizeof(size_t));
                    for ( size_t nr=0; nr<rfl2_r_data[r][n].size(); ++nr ) {
                        size = rfl2_r_data[r][n][nr].size();
                        fout.write((char*)&size,sizeof(size_t));
                        fout.write((char*)rfl2_r_data[r][n][nr].data(),
                                   sizeof(sxz<T>) * size);
                    }
                }
            }

            fout.close();
            if ( verbose ) cout << "done.\n";
        }
    }

    if ( verbose ) cout << "Normal termination of program.\n";
    return 0;
}


int main(int argc, char * argv[]) {

    input_parameters par;

    string fname = parse_input(argc, argv, par);
    if ( verbose ) {
        auto host_name = boost::asio::ip::host_name();
        char * user_name = getenv("USER");
        if (!user_name) {
            user_name = getenv("USERNAME");
        }
        cout << "\n*** Program ttcr2d ***\n\n"
        << "Raytracing in 2D media\n\n"
        << "Calculations performed on " << host_name
        << " by " << user_name << "\n";
    }

    get_params(fname, par);

    if ( verbose ) {
        switch (par.method) {
            case SHORTEST_PATH:
                cout << "Shortest path method selected.\n";
                break;
            case FAST_SWEEPING:
                cout << "Fast sweeping method selected.\n";
                break;
            case FAST_MARCHING:
                cout << "Fast marching method selected.\n";
                break;
            default:
                break;
        }
    }

    if ( par.singlePrecision ) {
        return body<float>(par);
    } else {
        return body<double>(par);
    }

}

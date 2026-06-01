//
//  ttcr_io.cpp
//  ttcr
//
//  Created by Bernard Giroux on 2012-11-19.
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

#include <cstring>
#include <fstream>
#include <sstream>

extern "C" {
#include <unistd.h>         // for getopt
}

#include "ttcr_io.h"

using namespace std;

namespace ttcr {

    int verbose = 0;

    void print_usage (std::ostream& stream, char *progname, int exit_code)
    {
        stream << "\n *** " << progname << " - Ray tracing by grid methods ***\n\n";
        stream << "Usage: " << progname << " [options] -p parameters.in\n";
        stream << "  -h  print this message\n"
        << "  -p  Specify parameter file (mandatory)\n"
        << "  -k  Save model in VTK format\n"
        << "  -v  Verbose mode\n"
        << "  -t  Measure time to build grid and perform raytracing\n"
        << "  -s  Dump secondary nodes to ascii file\n"
        << std::endl;
        exit (exit_code);
    }


    string parse_input(int argc, char * argv[], input_parameters &ip ) {

        string param_file;

        int next_option;
        const char* const short_options = "hk†p:vts";
        bool no_option = true;

        do {
            next_option = getopt (argc, argv, short_options);
            switch (next_option)
            {

                case  'h' : // -h or --help
                    // User has requested usage information. Print it to standard
                    //   output, and exit with exit code zero (normal termination).
                    no_option = false;
                    print_usage (cout, argv[0], 0);
                    break;

                case  'k' :
                    no_option = false;
                    ip.saveModelVTK = true;
                    break;

                case  'p' :
                    // This option takes an argument, the name of the input file.
                    no_option = false;
                    param_file = optarg;
                    break;

                case  'v' : // -v or --verbose
                    no_option = false;
                    verbose++;
                    break;

                case  't' :
                    no_option = false;
                    ip.time = true;
                    break;

                case  's' :
                    no_option = false;
                    ip.dump_secondary = true;
                    break;

                case  '?' : // The user specified an invalid option.
                    // Print usage information to standard error, and exit with exit
                    //code one (indicating abnormal termination).
                    print_usage (cerr, argv[0], 1);
                    break;

                case -1: // Done with options.
                    break;

                default: // Something else: unexpected.
                    abort ();
            }
        } while (next_option != -1);
        if ( no_option || param_file == "" ) print_usage (cout, argv[0], 0);

        return param_file;
    }

    void get_params(const std::string &filename, input_parameters &ip) {

        ifstream fin;
        fin.open( filename.c_str() );

        if ( !fin.is_open() ) {
            std::cerr << "Cannot open " << filename << std::endl;
            exit(1);
        }

        if ( verbose )
            std::cout << "\nReading file " << filename << std::endl;

        char value[100];
        char keyword[200];
        std::string kw;
        std::istringstream sin( value );

        while (!fin.eof()) {
            fin.get(value, 100, '#');
            if (strlen(value) == 0) {
                fin.clear();
                fin.getline(keyword, 200);
                continue;
            }
            fin.get(keyword, 200, ',');
            kw = keyword;

            if (kw.find("basename") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.basename;
            }
            else if (kw.find("modelfile") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.modelfile;
            }
            else if (kw.find("velfile") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.velfile;
            }
            else if (kw.find("slofile") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.slofile;
            }
            else if (kw.find("srcfile") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                string stmp;
                sin >> stmp;
                ip.srcfiles.push_back( stmp );
            }
            else if (kw.find("rcvfile") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.rcvfile;
            }
            else if (kw.find("secondary nodes") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                uint32_t val;
                size_t n=0;
                while ( sin >> val && n<3 ) {
                    ip.nn[n++] = val;
                }
                if ( n == 1 ) ip.nn[1] = ip.nn[2] = ip.nn[0];
            }
            else if (kw.find("number of threads") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.nt;
            }
            else if (kw.find("min nb Tx per thread") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.min_per_thread;
            }
            else if (kw.find("number of dynamic nodes") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.nTertiary;
            }
            else if (kw.find("tertiary nodes") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.nTertiary;
            }
            else if (kw.find("inverse distance") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.inverseDistance;
            }
            else if (kw.find("metric order") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.order;
            }
            else if (kw.find("epsilon") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.epsilon;
            }
            else if (kw.find("max number of iteration") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.nitermax;
            }
            else if (kw.find("saveGridTT") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.saveGridTT;
            }
            else if (kw.find("process reflectors") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.processReflectors;
            }
            else if (kw.find("single precision") < 200 ) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.singlePrecision;
            }
            else if (kw.find("saveRayPaths") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.saveRaypaths;
            }
            else if (kw.find("save M") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.saveM;
            }
            else if (kw.find("project Tx Rx") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.projectTxRx;
            }
            else if (kw.find("raypath high order") < 200 ||
                     kw.find("gradient method") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.raypath_method;
            }
            else if (kw.find("fast marching") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                int test;
                sin >> test;
                if ( test == 1 ) ip.method = FAST_MARCHING;
            }
            else if (kw.find("fast sweeping") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                int test;
                sin >> test;
                if ( test == 1 ) ip.method = FAST_SWEEPING;
                else if ( test == 2 ) ip.method = FAST_SWEEPING_OPENCL;
            }
            else if (kw.find("dynamic shortest path") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                int test;
                sin >> test;
                if ( test == 1 ) ip.method = DYNAMIC_SHORTEST_PATH;
            }
            else if (kw.find("source radius") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.source_radius;
            }
            else if (kw.find("src radius tertiary") < 200 ||
                     kw.find("radius dynamic nodes") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.radius_tertiary_nodes;
            }
            else if (kw.find("rotated template") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                int test;
                sin >> test;
                ip.rotated_template = (test == 1);
            }
            else if (kw.find("fsm high order") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                int test;
                sin >> test;
                ip.weno3 = (test == 1);
            }
            else if (kw.find("interpolate velocity") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                int test;
                sin >> test;
                ip.processVel = (test == 1);
            }
            else if (kw.find("traveltime from raypath") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                int test;
                sin >> test;
                ip.tt_from_rp = (test == 1);
            }
            else if (kw.find("raypath minimum distance") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                sin >> ip.min_distance_rp;
            }
            else if (kw.find("translate grid origin") < 200) {
                sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
                int test;
                sin >> test;
                ip.translateOrigin = (test == 1);
            }
            fin.getline(keyword, 200);
        }
        fin.close();
    }
}

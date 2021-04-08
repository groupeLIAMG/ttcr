//
//  msh2vtk_io.cpp
//  ttcr
//
//  Created by Bernard Giroux on 2014-10-18.
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

#include "msh2vtk_io.h"

#include <fstream>
#include <iostream>
#include <sstream>

extern "C" {
#include <unistd.h>         // for getopt
}

#include "msh2vtk_io.h"

using namespace std;

namespace ttcr {
    int verbose = 0;
}

void print_usage (std::ostream& stream, int exit_code)
{
    stream << "\n *** msh2vtk - Convert Gmsh format to VTK format ***\n\n";
    stream << "Usage: msh2vtk [options] -m mshFile -o vtkFile\n";
    stream << "  -h        print this message\n"
    << "  -d size   Cell size for rectilinear grids\n"
    << "  -f        Save physical reflectors\n"
    << "  -r        Save as rectilinear grid\n"
    << "  -c file   File with velocity values for the physical entities in the mesh\n"
    << "  -l file   File with slowness value at mesh nodes\n"
    << "  -s        Save slowness rather than velocity\n"
    << "  -t        Save CRT velocity file format\n"
    << "  -v        Verbose mode\n"
    << std::endl;
    exit (exit_code);
}


void parse_input(int argc, char * argv[], input_parameters &ip ) {

    int next_option;
    const char* const short_options = "hd:m:o:frc:l:stz:v";
    bool no_option = true;
    do {
        next_option = getopt (argc, argv, short_options);
        switch (next_option)
        {

            case  'h' : // -h or --help
                // User has requested usage information. Print it to standard
                //   output, and exit with exit code zero (normal termination).
                no_option = false;
                print_usage (cout, 0);

            case  'd' :
                no_option = false;
                ip.d = atof(optarg);
                break;

            case  'm' :
                no_option = false;
                ip.mshFile = optarg;
                break;

            case  'o' :
                no_option = false;
                ip.vtkFile = optarg;
                break;

            case  'r' :
                no_option = false;
                ip.rectilinear = true;
                break;

            case  'f' :
                no_option = false;
                ip.saveReflectors = true;
                break;

            case  'c' :
                no_option = false;
                ip.velFile = optarg;
                break;

            case  'l' :
                no_option = false;
                ip.sloFile = optarg;
                break;

            case  's' :
                no_option = false;
                ip.saveSlowness = true;
                break;

            case  't' :
                no_option = false;
                ip.crt = true;
                break;

            case  'v' : // -v or --verbose
                no_option = false;
                ttcr::verbose++;
                break;

            case  '?' : // The user specified an invalid option.
                // Print usage information to standard error, and exit with exit
                //code one (indicating abnormal termination).
                print_usage (cerr, 1);

            case -1: // Done with options.
                break;

            default: // Something else: unexpected.
                abort ();
        }
    } while (next_option != -1);
    if ( no_option || ip.mshFile == "" || ip.vtkFile == "" ) print_usage (cout, 0);

}

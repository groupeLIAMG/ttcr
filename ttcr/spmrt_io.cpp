//
//  spmrt_io.cpp
//  ttcr_u
//
//  Created by Bernard Giroux on 2012-11-19.
//  Copyright (c) 2012 Bernard Giroux. All rights reserved.
//

#include <cstring>
#include <fstream>
#include <sstream>

extern "C" {
#include <unistd.h>         // for getopt
}

#include "spmrt_io.h"


using namespace std;

void print_usage (std::ostream& stream, char *progname, int exit_code)
{
    stream << "\n *** " << progname << " - Ray tracing in 3D ***\n\n";
    stream << "Usage: " << progname << " [options] -p paramters.in\n";
    stream << "  -h  print this message\n"
    << "  -p  Specify parameter file (mandatory)\n"
	<< "  -r  Save Ray paths\n"
    << "  -s  Do computation in single precision\n"
    << "  -m  Use fast-marching method\n"
    << "  -w  Use fast-sweeping method\n"
	<< "  -k  Save model in VTK format\n"
    << "  -v  Verbose mode\n"
	<< "  -t  Measure time to build grid and perform raytracing\n"
    << std::endl;
    exit (exit_code);
}


string parse_input(int argc, char * argv[], input_parameters &ip ) {
	
	string param_file;
	
	int next_option;
    const char* const short_options = "hksp:rvtmw";
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
				
			case  'k' :
				no_option = false;
                ip.saveModelVTK = true;
                break;
                
			case  's' :
				no_option = false;
                ip.singlePrecision = true;
                break;
                
            case  'p' :
                // This option takes an argument, the name of the input file.
                no_option = false;
                param_file = optarg;
                break;
                
			case  'r' :
				no_option = false;
                ip.saveRaypaths = true;
                break;
                
            case  'v' : // -v or --verbose
                no_option = false;
                ip.verbose++;
                break;
                
            case  't' :
                no_option = false;
                ip.time = true;
                break;
                
            case  'm' :
                no_option = false;
                ip.method = FAST_MARCHING;
                break;
                
            case  'w' :
                no_option = false;
                ip.method = FAST_SWEEPING;
                break;
                
            case  '?' : // The user specified an invalid option.
                // Print usage information to standard error, and exit with exit
                //code one (indicating abnormal termination).
                print_usage (cerr, argv[0], 1);
                
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
	
	if ( ip.verbose )
		std::cout << "\nReading file " << filename << std::endl;
	
	char value[100];
	char parameter[200];
	std::string par;
	std::istringstream sin( value );
	
	while (!fin.eof()) {
		fin.get(value, 100, '#');
		if (strlen(value) == 0) {
			fin.clear();
			fin.getline(parameter, 200);
			continue;
		}
		fin.get(parameter, 200, ',');
		par = parameter;
		
		if (par.find("basename") < 200) {
			sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
			sin >> ip.basename;
		}
		else if (par.find("modelfile") < 200) {
			sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
			sin >> ip.modelfile;
		}
		else if (par.find("velfile") < 200) {
			sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
			sin >> ip.velfile;
		}
		else if (par.find("slofile") < 200) {
			sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
			sin >> ip.slofile;
		}
		else if (par.find("srcfile") < 200) {
			sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
            string stmp;
			sin >> stmp;
            ip.srcfiles.push_back( stmp );
		}
		else if (par.find("rcvfile") < 200) {
			sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
			sin >> ip.rcvfile;
		}
		else if (par.find("secondary nodes") < 200) {
			sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
            uint32_t val;
            size_t n=0;
            while ( sin >> val && n<3 ) {
                ip.nn[n++] = val;
            }
            if ( n == 1 ) ip.nn[1] = ip.nn[2] = ip.nn[0];
		}
		else if (par.find("number of threads") < 200) {
			sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
			sin >> ip.nt;
		}
        else if (par.find("inverse distance") < 200) {
            sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
			sin >> ip.inverseDistance;
        }
		else if (par.find("metric order") < 200) {
			sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
			sin >> ip.order;
		}
		else if (par.find("epsilon") < 200) {
			sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
			sin >> ip.epsilon;
		}
		else if (par.find("max number if iteration") < 200) {
			sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
			sin >> ip.nitermax;
		}
        else if (par.find("saveGridTT") < 200) {
			sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
			sin >> ip.saveGridTT;
        }
		else if (par.find("process reflectors") < 200) {
			sin.str( value ); sin.seekg(0, std::ios_base::beg); sin.clear();
			sin >> ip.processReflectors;
        }
		fin.getline(parameter, 200);
	}
	fin.close();

}
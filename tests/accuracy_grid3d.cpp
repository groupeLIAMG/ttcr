//
//  accuracy_grid3d.cpp
//
//  Assess the accuracy of the 3D rectilinear grid methods (non OpenCL)
//  for both float and double.
//
//  Two studies are performed:
//
//   1. Convergence study (layers, gradient) -- medium and fine resolution.
//      Single source / fixed receivers (as in test_grid3d.cpp); the mean
//      relative error is computed against the analytic reference solutions
//      stored in the ./files/*_tt.vtr files (which are sampled at the
//      receiver locations and therefore independent of grid resolution).
//
//   2. Constant-velocity study -- medium and fine resolution.
//      A configurable number of sources (default 100) is placed at random
//      inside the domain.  The velocity is constant, so the analytic travel
//      time is exactly  t = slowness * euclidean_distance ; no reference
//      file is needed.  The sources are run with the threaded raytracer so
//      the many shots are spread across cores.
//
//  Output: a CSV file (accuracy_grid3d.csv by default) with one row per
//  (precision, model, method, resolution) combination, holding the mean
//  relative error and the total raytracing wall time.
//
//  Run from the tests/ folder so the ./files/... paths resolve.
//
//  Usage:
//      ./accuracy_grid3d [options]
//        --nsrc N         number of random sources, constant study (default 100)
//        --seed S         RNG seed for the random sources (default 12345)
//        --threads N      threads for the constant study (default: hw concurrency)
//        --methods LIST   comma list among fsm,sp,dsp (default: all)
//        --precision LIST comma list among float,double (default: both)
//        --no-medium      skip medium-resolution models
//        --no-fine        skip fine-resolution models
//        --no-ref         skip the convergence study (layers/gradient)
//        --no-const       skip the constant-velocity study
//        -o FILE          output CSV file (default accuracy_grid3d.csv)
//

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <typeinfo>
#include <vector>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdocumentation-deprecated-sync"
#pragma clang diagnostic ignored "-Wimplicit-int-conversion"
#pragma clang diagnostic ignored "-Wshorten-64-to-32"
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#include <vtkPointData.h>
#include <vtkRectilinearGrid.h>
#include <vtkSmartPointer.h>
#include <vtkXMLRectilinearGridReader.h>
#pragma clang diagnostic pop

#include "Grid3D.h"
#include "Rcv.h"
#include "Src.h"
#include "structs_ttcr.h"
#include "grids.h"
#include "utils.h"

namespace ttcr {
    int verbose = 0;
    int gpu_profile = 0;
}

using namespace std;
using namespace ttcr;

// ---------------------------------------------------------------------------
// Rectilinear, non-OpenCL methods only
// ---------------------------------------------------------------------------
struct method_info {
    raytracing_method method;
    string key;     // short selector used by --methods
    string name;    // pretty name in the CSV / console
};
const vector<method_info> all_methods = {
    { FAST_SWEEPING,         "fsm", "FAST_SWEEPING"         },  // -> Grid3Dr?fs
    { SHORTEST_PATH,         "sp",  "SHORTEST_PATH"         },  // -> Grid3Dr?sp
    { DYNAMIC_SHORTEST_PATH, "dsp", "DYNAMIC_SHORTEST_PATH" },  // -> Grid3Dr?dsp
};

// Models for the convergence study (analytic reference stored in a .vtr)
struct ref_model {
    string model;
    string reference;
    string name;
    string resolution;
};
const vector<ref_model> ref_models = {
    { "./files/layers_medium.vtr",   "./files/sol_analytique_couches_tt.vtr",   "layers",   "medium" },
    { "./files/layers_fine.vtr",     "./files/sol_analytique_couches_tt.vtr",   "layers",   "fine"   },
    { "./files/gradient_medium.vtr", "./files/sol_analytique_gradient_tt.vtr",  "gradient", "medium" },
    { "./files/gradient_fine.vtr",   "./files/sol_analytique_gradient_tt.vtr",  "gradient", "fine"   },
};

// Constant-velocity models (analytic solution computed directly)
struct const_model {
    string model;
    string name;
    string resolution;
};
const vector<const_model> const_models = {
    { "./files/constant_medium.vtr", "constant", "medium" },
    { "./files/constant_fine.vtr",   "constant", "fine"   },
};

// ---------------------------------------------------------------------------
// Run-time options
// ---------------------------------------------------------------------------
struct options {
    int nsrc = 100;
    unsigned seed = 12345;
    int threads = 0;            // 0 -> hardware_concurrency
    bool do_medium = true;
    bool do_fine = true;
    bool do_ref = true;
    bool do_const = true;
    vector<string> method_keys = { "fsm", "sp", "dsp" };
    vector<string> precisions = { "double", "float" };
    string out = "accuracy_grid3d.csv";
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
template<typename T>
string get_class_name(const Grid3D<T,uint32_t> *g) {
    string name = typeid(*g).name();
    size_t start = name.find("Grid3D");
    name = name.substr(start, 11);
    if ((name[8] == 'f' && name[9] == 's') ||
        (name[8] == 's' && name[9] == 'p') ) {
        name.pop_back();
    }
    return name;
}

template<typename T>
double get_rel_error(const string& filename, const Rcv<T>& rcv) {
    // mean relative error against the analytic solution stored in filename
    vector<double> ref_tt(rcv.get_coord().size());
    vtkSmartPointer<vtkXMLRectilinearGridReader> reader =
    vtkSmartPointer<vtkXMLRectilinearGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    reader->GetOutput()->Register(reader);
    vtkRectilinearGrid *dataSet = reader->GetOutput();

    vtkPointData *pd = dataSet->GetPointData();
    for ( size_t n=0; n<rcv.get_coord().size(); ++n ) {
        double x[3] = {static_cast<double>(rcv.get_coord()[n].x),
            static_cast<double>(rcv.get_coord()[n].y),
            static_cast<double>(rcv.get_coord()[n].z)};
        vtkIdType ii = dataSet->FindPoint(x);
        ref_tt[n] = static_cast<double>(pd->GetArray(0)->GetTuple1(ii));
    }
    double error = 0.0;
    for ( size_t n=1; n<ref_tt.size(); ++n ) {
        // start at 1 to skip the node at the source location where tt = 0
        error += abs( (ref_tt[n] - rcv.get_tt(0)[n])/ref_tt[n] );
    }
    error /= (ref_tt.size() - 1);
    return error;
}

void set_method_params(input_parameters& par, raytracing_method method) {
    switch ( method ) {
        case FAST_SWEEPING:
            par.weno3 = 1;
            break;
        case SHORTEST_PATH:
            par.nn[0] = 5; par.nn[1] = 5; par.nn[2] = 5;
            break;
        case DYNAMIC_SHORTEST_PATH:
            par.radius_tertiary_nodes = 3.0;
            par.nn[0] = 2; par.nn[1] = 2; par.nn[2] = 2;
            break;
        default:
            break;
    }
}

bool method_selected(const options& opt, const string& key) {
    return find(opt.method_keys.begin(), opt.method_keys.end(), key)
           != opt.method_keys.end();
}

// ---------------------------------------------------------------------------
// Study 1 : convergence (layers / gradient) -- single source
// ---------------------------------------------------------------------------
template<typename T>
void run_reference(ofstream& csv, const string& precision_name,
                   const options& opt) {
    for ( const auto& mod : ref_models ) {
        if ( mod.resolution == "medium" && !opt.do_medium ) continue;
        if ( mod.resolution == "fine"   && !opt.do_fine )   continue;

        for ( const auto& mi : all_methods ) {
            if ( !method_selected(opt, mi.key) ) continue;

            Src<T> src("./files/src.dat");
            src.init();
            Rcv<T> rcv("./files/rcv.dat");
            rcv.init(1);

            input_parameters par;
            par.method = mi.method;
            set_method_params(par, mi.method);
            par.modelfile = mod.model;

            Grid3D<T,uint32_t> *g = buildRectilinear3DfromVtr<T>(par, 1);

            auto t0 = chrono::high_resolution_clock::now();
            try {
                g->raytrace(src.get_coord(), src.get_t0(),
                            rcv.get_coord(), rcv.get_tt(0));
            } catch (std::exception& e) {
                std::cerr << e.what() << std::endl;
                delete g;
                continue;
            }
            auto t1 = chrono::high_resolution_clock::now();
            double dt = chrono::duration<double>(t1 - t0).count();

            double error = get_rel_error(mod.reference, rcv);
            string cls = get_class_name(g);

            std::cout << "  " << precision_name << "  " << mod.name << " ("
                      << mod.resolution << ")  " << mi.name << "  [" << cls << "]"
                      << "  error = " << error << "  time = " << dt << " s\n";
            std::cout.flush();

            csv << precision_name << ',' << mod.name << ',' << mi.name << ','
                << cls << ',' << mod.resolution << ',' << error << ',' << dt
                << '\n';
            csv.flush();

            delete g;
        }
    }
}

// ---------------------------------------------------------------------------
// Study 2 : constant velocity -- many random sources, analytic solution
// ---------------------------------------------------------------------------
template<typename T>
void run_constant(ofstream& csv, const string& precision_name,
                  const options& opt, const vector<array<double,3>>& src_xyz) {

    size_t nthreads = opt.threads > 0 ? static_cast<size_t>(opt.threads)
                                      : std::max(1u, std::thread::hardware_concurrency());

    for ( const auto& mod : const_models ) {
        if ( mod.resolution == "medium" && !opt.do_medium ) continue;
        if ( mod.resolution == "fine"   && !opt.do_fine )   continue;

        // receivers (shared by every shot)
        Rcv<T> rcv("./files/rcv.dat");
        rcv.init(1);
        const vector<sxyz<T>>& rx = rcv.get_coord();

        // shots: one source per shot
        vector<vector<sxyz<T>>> Tx(src_xyz.size());
        vector<vector<T>> t0(src_xyz.size());
        vector<vector<sxyz<T>>> Rx(src_xyz.size(), rx);
        vector<vector<T>> tt(src_xyz.size());
        for ( size_t n=0; n<src_xyz.size(); ++n ) {
            Tx[n].push_back(sxyz<T>(static_cast<T>(src_xyz[n][0]),
                                    static_cast<T>(src_xyz[n][1]),
                                    static_cast<T>(src_xyz[n][2])));
            t0[n].push_back(static_cast<T>(0));
            tt[n].resize(rx.size());
        }

        for ( const auto& mi : all_methods ) {
            if ( !method_selected(opt, mi.key) ) continue;

            input_parameters par;
            par.method = mi.method;
            set_method_params(par, mi.method);
            par.modelfile = mod.model;

            Grid3D<T,uint32_t> *g = buildRectilinear3DfromVtr<T>(par, nthreads);

            // constant slowness read back from the grid (avoids hardcoding)
            vector<T> slo;
            g->getSlowness(slo);
            double s0 = static_cast<double>(slo[0]);

            auto t_start = chrono::high_resolution_clock::now();
            try {
                g->raytrace(Tx, t0, Rx, tt);
            } catch (std::exception& e) {
                std::cerr << e.what() << std::endl;
                delete g;
                continue;
            }
            auto t_end = chrono::high_resolution_clock::now();
            double dt = chrono::duration<double>(t_end - t_start).count();

            // mean relative error over all source-receiver pairs
            double error = 0.0;
            size_t npair = 0;
            for ( size_t n=0; n<src_xyz.size(); ++n ) {
                for ( size_t m=0; m<rx.size(); ++m ) {
                    double dx = static_cast<double>(rx[m].x) - src_xyz[n][0];
                    double dy = static_cast<double>(rx[m].y) - src_xyz[n][1];
                    double dz = static_cast<double>(rx[m].z) - src_xyz[n][2];
                    double dist = sqrt(dx*dx + dy*dy + dz*dz);
                    double ref = s0 * dist;
                    if ( ref == 0.0 ) continue;   // source == receiver
                    error += abs((ref - static_cast<double>(tt[n][m])) / ref);
                    ++npair;
                }
            }
            error /= npair;
            string cls = get_class_name(g);

            std::cout << "  " << precision_name << "  " << mod.name << " ("
                      << mod.resolution << ")  " << mi.name << "  [" << cls << "]"
                      << "  error = " << error << "  time = " << dt << " s  ("
                      << src_xyz.size() << " src)\n";
            std::cout.flush();

            csv << precision_name << ',' << mod.name << ',' << mi.name << ','
                << cls << ',' << mod.resolution << ',' << error << ',' << dt
                << '\n';
            csv.flush();

            delete g;
        }
    }
}

// ---------------------------------------------------------------------------
// Random source positions inside the [margin, 20-margin]^3 domain.
// Generated from doubles with a fixed seed so float and double runs use the
// exact same source set.
// ---------------------------------------------------------------------------
vector<array<double,3>> make_sources(int n, unsigned seed) {
    const double lo = 0.5, hi = 19.5;
    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<double> dist(lo, hi);
    vector<array<double,3>> src(n);
    for ( int i=0; i<n; ++i )
        src[i] = { dist(gen), dist(gen), dist(gen) };
    return src;
}

// ---------------------------------------------------------------------------
// Minimal CLI parsing
// ---------------------------------------------------------------------------
vector<string> split_csv(const string& s) {
    vector<string> out;
    size_t start = 0;
    while ( start < s.size() ) {
        size_t comma = s.find(',', start);
        if ( comma == string::npos ) { out.push_back(s.substr(start)); break; }
        out.push_back(s.substr(start, comma-start));
        start = comma + 1;
    }
    return out;
}

bool parse_args(int argc, const char* argv[], options& opt) {
    for ( int i=1; i<argc; ++i ) {
        string a = argv[i];
        auto next = [&]() -> string {
            if ( i+1 >= argc ) { throw runtime_error("missing value for " + a); }
            return argv[++i];
        };
        if      ( a == "--nsrc" )      opt.nsrc = stoi(next());
        else if ( a == "--seed" )      opt.seed = static_cast<unsigned>(stoul(next()));
        else if ( a == "--threads" )   opt.threads = stoi(next());
        else if ( a == "--methods" )   opt.method_keys = split_csv(next());
        else if ( a == "--precision" ) opt.precisions = split_csv(next());
        else if ( a == "--no-medium" ) opt.do_medium = false;
        else if ( a == "--no-fine" )   opt.do_fine = false;
        else if ( a == "--no-ref" )    opt.do_ref = false;
        else if ( a == "--no-const" )  opt.do_const = false;
        else if ( a == "-o" )          opt.out = next();
        else if ( a == "-h" || a == "--help" ) {
            std::cout <<
                "Usage: ./accuracy_grid3d [options]\n"
                "  --nsrc N         random sources, constant study (default 100)\n"
                "  --seed S         RNG seed (default 12345)\n"
                "  --threads N      threads for constant study (default: hw)\n"
                "  --methods LIST   comma list among fsm,sp,dsp (default all)\n"
                "  --precision LIST comma list among float,double (default both)\n"
                "  --no-medium      skip medium-resolution models\n"
                "  --no-fine        skip fine-resolution models\n"
                "  --no-ref         skip the convergence study\n"
                "  --no-const       skip the constant-velocity study\n"
                "  -o FILE          output CSV (default accuracy_grid3d.csv)\n";
            return false;
        }
        else {
            std::cerr << "unknown option: " << a << " (try --help)\n";
            return false;
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
int main(int argc, const char* argv[]) {
    options opt;
    try {
        if ( !parse_args(argc, argv, opt) ) return 0;
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    ofstream csv(opt.out);
    csv << "precision,model,method,class,resolution,error,time\n";

    auto src_xyz = make_sources(opt.nsrc, opt.seed);

    std::cout << "Assessing accuracy of 3D rectilinear grid methods "
                 "(non OpenCL)\n";

    if ( opt.do_ref ) {
        std::cout << "\n--- Convergence study (layers / gradient) ---\n";
        for ( const auto& prec : opt.precisions ) {
            if ( prec == "double" ) run_reference<double>(csv, prec, opt);
            else if ( prec == "float" ) run_reference<float>(csv, prec, opt);
        }
    }

    if ( opt.do_const ) {
        std::cout << "\n--- Constant-velocity study (" << opt.nsrc
                  << " random sources) ---\n";
        for ( const auto& prec : opt.precisions ) {
            if ( prec == "double" ) run_constant<double>(csv, prec, opt, src_xyz);
            else if ( prec == "float" ) run_constant<float>(csv, prec, opt, src_xyz);
        }
    }

    csv.close();
    std::cout << "\nResults written to " << opt.out << '\n';
    return 0;
}

//
//  main.cpp
//  test_grid3d
//
//  Created by Bernard Giroux on 2021-03-04.
//  Copyright © 2021 Bernard Giroux. All rights reserved.
//

#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <vector>

#define BOOST_TEST_MODULE Test Grid2D
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdocumentation-deprecated-sync"
#pragma clang diagnostic ignored "-Wimplicit-int-conversion"
#pragma clang diagnostic ignored "-Wshorten-64-to-32"
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

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
}

namespace bdata = boost::unit_test::data;

using namespace std;
using namespace ttcr;

string get_class_name(const Grid3D<double,uint32_t> *g) {
    string name = typeid(*g).name();
    size_t start = name.find("Grid3D");
    name = name.substr(start, 11);
    if ((name[8] == 'f' && name[9] == 's') ||
        (name[8] == 's' && name[9] == 'p') ) {
        name.pop_back();
    }
    return name;
}

double get_rel_error(const string& filename, const Rcv<double>& rcv) {
    // compute relative error
    
    vector<double> ref_tt(rcv.get_coord().size());
    vtkSmartPointer<vtkXMLRectilinearGridReader> reader =
    vtkSmartPointer<vtkXMLRectilinearGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    reader->GetOutput()->Register(reader);
    vtkRectilinearGrid *dataSet = reader->GetOutput();
    
    vtkPointData *pd = dataSet->GetPointData();
    for ( size_t n=0; n<rcv.get_coord().size(); ++n ) {
        vtkIdType ii = dataSet->FindPoint(rcv.get_coord()[n].x,
                                          rcv.get_coord()[n].y,
                                          rcv.get_coord()[n].z);
        ref_tt[n] = static_cast<double>(pd->GetArray(0)->GetTuple1(ii));
    }
    // compute relative error
    double error = 0.0;
    for ( size_t n=1; n<ref_tt.size(); ++n ) {
        // start at 1 to avoid node at source location where tt = 0
        error += abs( (ref_tt[n] - rcv.get_tt(0)[n])/ref_tt[n] );
    }
    error /= (ref_tt.size() - 1);
    return error;
}

const char* models[] = {
    "./files/layers_medium.vtr",
    "./files/gradient_medium.vtr",
    "./files/layers_medium.vtu",
    "./files/gradient_medium.vtu"
};
const char* references[] = {
    "./files/sol_analytique_couches_tt.vtr",
    "./files/sol_analytique_gradient_tt.vtr",
    "./files/sol_analytique_couches_tt.vtr",
    "./files/sol_analytique_gradient_tt.vtr"
};
raytracing_method methods[] = {
    FAST_SWEEPING,
    SHORTEST_PATH,
    DYNAMIC_SHORTEST_PATH
};

const char* models_tr[] = {
    "./files/layers_medium2.vtr",
    "./files/gradient_medium2.vtr",
    "./files/layers_medium2.vtu",
    "./files/gradient_medium2.vtu"
};
raytracing_method methods_tr[] = {
    DYNAMIC_SHORTEST_PATH,
    SHORTEST_PATH
};

BOOST_DATA_TEST_CASE(
                     testGrid3D,
                     (bdata::make(models) ^ bdata::make(references)) * bdata::make(methods),
                     model, ref, method) {
    Src<double> src("./files/src.dat");
    src.init();
    Rcv<double> rcv("./files/rcv.dat");
    rcv.init(1);

    input_parameters par;
    par.method = method;
    switch(method) {
        case FAST_SWEEPING:
            par.weno3 = 1;
            break;
        case SHORTEST_PATH:
            par.nn[0] = 5;
            par.nn[1] = 5;
            par.nn[2] = 5;
            break;
        case DYNAMIC_SHORTEST_PATH:
            par.radius_tertiary_nodes = 3.0;
            par.nn[0] = 2;
            par.nn[1] = 2;
            par.nn[2] = 2;
            break;
        default:
            // do nothing
            break;
    }
    par.modelfile = model;
    Grid3D<double,uint32_t> *g;
    double error_threshold = .01;
    if (string(model).find("vtr") != string::npos) {
        g = buildRectilinear3DfromVtr<double>(par, 1);
    } else {
        error_threshold = .07;
        g = buildUnstructured3DfromVtu<double>(par, 1);
    }
    try {
        g->raytrace(src.get_coord(), src.get_t0(), rcv.get_coord(), rcv.get_tt(0));
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        abort();
    }
    string filename = "./files/" + get_class_name(g) + "_tt_grid";
    g->saveTT(filename, 0, 0, 2);
    double error = get_rel_error(ref, rcv);
    BOOST_TEST_MESSAGE( "\t\t" << get_class_name(g) << " - error = " << error );

    BOOST_TEST(error < error_threshold);
}

BOOST_DATA_TEST_CASE(
                     testGrid3D_rp,
                     (bdata::make(models) ^ bdata::make(references)) * bdata::make(methods),
                     model, ref, method) {
    Src<double> src("./files/src3d_in.dat");
    src.init();
    Rcv<double> rcv("./files/rcv3d_in.dat");
    rcv.init(1);

    input_parameters par;
    par.method = method;
    switch(method) {
        case FAST_SWEEPING:
            par.weno3 = 1;
            break;
        case SHORTEST_PATH:
            par.nn[0] = 5;
            par.nn[1] = 5;
            par.nn[2] = 5;
            break;
        case DYNAMIC_SHORTEST_PATH:
            par.radius_tertiary_nodes = 3.0;
            par.nn[0] = 2;
            par.nn[1] = 2;
            par.nn[2] = 2;
            break;
        default:
            // do nothing
            break;
    }
    par.modelfile = model;
    Grid3D<double,uint32_t> *g;
    if (string(model).find("vtr") != string::npos) {
        g = buildRectilinear3DfromVtr<double>(par, 1);
    } else {
        g = buildUnstructured3DfromVtu<double>(par, 1);
    }
    vector<vector<sxyz<double>>> r_data;
    try {
        g->raytrace(src.get_coord(), src.get_t0(), rcv.get_coord(), rcv.get_tt(0), r_data);
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        abort();
    }
    string filename = "./files/" + get_class_name(g) + "_rp.vtp";
    saveRayPaths(filename, r_data);

    double error = get_rel_error(ref, rcv);
    BOOST_TEST_MESSAGE( "\t\t" << get_class_name(g) << ", r_data - error = " << error );

    BOOST_TEST(error < 0.15);
}

BOOST_DATA_TEST_CASE(
                     translate,
                     (bdata::make(models) ^ bdata::make(models_tr)) * bdata::make(methods_tr),
                     model, model2, method) {
    Src<double> src("./files/src3d_in.dat");
    src.init();
    Rcv<double> rcv("./files/rcv3d_in.dat");
    rcv.init(1);

    Src<double> src2("./files/src3d_in2.dat");
    src2.init();
    Rcv<double> rcv2("./files/rcv3d_in2.dat");
    rcv2.init(1);
    Rcv<double> rcv3("./files/rcv3d_in2.dat");
    rcv3.init(1);

    input_parameters par;
    par.method = method;
    switch(method) {
        case SHORTEST_PATH:
            par.nn[0] = 5;
            par.nn[1] = 5;
            par.nn[2] = 5;
            break;
        case DYNAMIC_SHORTEST_PATH:
            par.radius_tertiary_nodes = 3.0;
            par.nn[0] = 2;
            par.nn[1] = 2;
            par.nn[2] = 2;
            break;
        default:
            // do nothing
            break;
    }
    par.modelfile = model;
    Grid3D<double,uint32_t> *g;
    if (string(model).find("vtr") != string::npos) {
        g = buildRectilinear3DfromVtr<double>(par, 1);
    } else {
        g = buildUnstructured3DfromVtu<double>(par, 1);
    }

    par.modelfile = model2;
    Grid3D<double,uint32_t> *g2;
    if (string(model2).find("vtr") != string::npos) {
        g2 = buildRectilinear3DfromVtr<double>(par, 1);
    } else {
        g2 = buildUnstructured3DfromVtu<double>(par, 1);
    }

    par.translateOrigin = true;
    Grid3D<double,uint32_t> *g3;
    if (string(model2).find("vtr") != string::npos) {
        g3 = buildRectilinear3DfromVtr<double>(par, 1);
    } else {
        g3 = buildUnstructured3DfromVtu<double>(par, 1);
    }

    vector<vector<sxyz<double>>> r_data;
    vector<vector<sxyz<double>>> r_data2;
    vector<vector<sxyz<double>>> r_data3;

    bool ignore2 = false;
    try {
        g->raytrace(src.get_coord(), src.get_t0(), rcv.get_coord(), rcv.get_tt(0), r_data);
    } catch (std::exception& e) {
        BOOST_TEST_MESSAGE( "\n\n\t\tRaytracing failed for model: " << model
                           << ", with translateGrid false\n\t\t" << e.what() << "\n\n");
    }

    try {
        g2->raytrace(src2.get_coord(), src2.get_t0(), rcv2.get_coord(), rcv2.get_tt(0), r_data2);
    } catch (std::exception& e) {
        BOOST_TEST_MESSAGE( "\n\n\t\tRaytracing failed for model: " << model2
                           << ", with translateGrid false\n\t\t" << e.what()
                           << "\n\t\tThis is expected.\n\n");
        // this is expected
        ignore2 = true;
    }

    try {
        g3->raytrace(src2.get_coord(), src2.get_t0(), rcv3.get_coord(), rcv3.get_tt(0), r_data3);
    } catch (std::exception& e) {
        BOOST_TEST_MESSAGE( "\n\n\t\tRaytracing failed for model: " << model2
                           << ", with translateGrid true\n\t\t" << e.what() << "\n\n");
    }

    string filename = "./files/" + get_class_name(g) + "_rp_tr1.vtp";
    saveRayPaths(filename, r_data);
    filename = "./files/" + get_class_name(g) + "_rp_tr2.vtp";
    saveRayPaths(filename, r_data2);
    filename = "./files/" + get_class_name(g) + "_rp_tr3.vtp";
    saveRayPaths(filename, r_data3);

    double diff = 0.0;
    double diff2 = 0.0;
    vector<double> tt = rcv.get_tt(0);
    vector<double> tt2 = rcv2.get_tt(0);
    vector<double> tt3 = rcv3.get_tt(0);
    for (size_t n=1; n<tt.size(); ++n) {
        diff += abs((tt[n] - tt2[n]) / tt[n]);
        diff2 += abs((tt[n] - tt3[n]) / tt[n]);
    }
    diff /= (tt.size() - 1);
    diff2 /= (tt.size() - 1);
    BOOST_TEST_MESSAGE( "\t\t" << get_class_name(g) << ",\t translate - diff = "
                       << diff << ",\t diff2 = " << diff2);
    if (!ignore2 ) {
        BOOST_TEST(diff < 0.01);
    }
    BOOST_TEST(diff2 < 0.01);
}

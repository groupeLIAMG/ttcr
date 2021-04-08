//
//  main.cpp
//  test_grid2d
//
//  Created by Bernard Giroux on 2021-02-24.
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

#include "Grid2D.h"
#include "Rcv2D.h"
#include "Src2D.h"
#include "structs_ttcr.h"
#include "grids.h"
#include "utils.h"

namespace ttcr {
    int verbose = 0;
}

namespace bdata = boost::unit_test::data;

using namespace std;
using namespace ttcr;

string get_class_name(const Grid2D<double,uint32_t,sxz<double>> *g) {
    string name = typeid(*g).name();
    size_t start = name.find("Grid2D");
    name = name.substr(start, 11);
    if ((name[8] == 'f' && name[9] == 's') ||
        (name[8] == 's' && name[9] == 'p') ) {
        name.pop_back();
    }
    return name;
}

double get_rel_error(const string& filename, const Rcv2D<double>& rcv) {
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
        vtkIdType ii = dataSet->FindPoint(rcv.get_coord()[n].x, 0.0, rcv.get_coord()[n].z);
        ref_tt[n] = static_cast<double>(pd->GetArray(0)->GetTuple1(ii));
    }
    // compute relative error
    double error = 0.0;
    for ( size_t n=1; n<ref_tt.size(); ++n ) {
        // start at 1 to avoid node at source location where tt = 0
        error += abs( (ref_tt[n] - rcv.get_tt(0)[n])/ref_tt[n] );
//        cout << rcv.get_coord()[n] << '\t' << ref_tt[n] << '\t' << rcv.get_tt(0)[n] << '\n';
    }
    error /= (ref_tt.size() - 1);
    return error;
}


const char* models[] = {
    "./files/layers_fine2d.vtr",
    "./files/gradient_fine2d.vtr",
    "./files/layers_fine2d.vtu",
    "./files/gradient_fine2d.vtu"
};
const char* models_coarse[] = {
    "./files/layers_coarse2d.vtr",
    "./files/gradient_coarse2d.vtr",
    "./files/layers_coarse2d.vtu",
    "./files/gradient_gmsh2d.vtu"
};
const char* references[] = {
    "./files/sol_analytique_couches2d_tt.vtr",
    "./files/sol_analytique_gradient2d_tt.vtr",
    "./files/sol_analytique_couches2d_tt.vtr",
    "./files/sol_analytique_gradient2d_tt.vtr"
};
raytracing_method methods[] = { DYNAMIC_SHORTEST_PATH, FAST_SWEEPING, SHORTEST_PATH };


BOOST_DATA_TEST_CASE(
                     testGrid2D,
                     (bdata::make(models) ^ bdata::make(references)) * bdata::make(methods),
                     model, ref, method) {
    Src2D<double> src("./files/src2d.dat");
    src.init();
    Rcv2D<double> rcv("./files/rcv2d.dat");
    rcv.init(1);

    input_parameters par;
    par.method = method;
    switch(method) {
        case FAST_SWEEPING:
            par.weno3 = 1;
            break;
        case SHORTEST_PATH:
            par.nn[0] = 10;
            par.nn[1] = 10;
            par.nn[2] = 10;
            break;
        case DYNAMIC_SHORTEST_PATH:
            par.radius_tertiary_nodes = 3.0;
            par.nn[0] = 3;
            par.nn[1] = 3;
            par.nn[2] = 3;
            break;
        default:
            // do nothing
            break;
    }
    par.modelfile = model;
    Grid2D<double,uint32_t,sxz<double>> *g;
    if (string(model).find("vtr") != string::npos) {
        g = buildRectilinear2DfromVtr<double>(par, 1);
    } else {
        g = buildUnstructured2DfromVtu<double>(par, 1);
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

    BOOST_TEST(error < 0.02);
}

BOOST_DATA_TEST_CASE(
                     testGrid2D_ttrp,
                     (bdata::make(models_coarse) ^ bdata::make(references)) * bdata::make(methods),
                     model, ref, method) {
    Src2D<double> src("./files/src2d_in.dat");
    src.init();
    Rcv2D<double> rcv("./files/rcv2d_in.dat");
    rcv.init(1);

    input_parameters par;
    par.method = method;
    par.tt_from_rp = true;
    switch(method) {
        case FAST_SWEEPING:
            par.weno3 = 1;
            break;
        case SHORTEST_PATH:
            par.nn[0] = 10;
            par.nn[1] = 10;
            par.nn[2] = 10;
            break;
        case DYNAMIC_SHORTEST_PATH:
            par.radius_tertiary_nodes = 3.0;
            par.nn[0] = 3;
            par.nn[1] = 3;
            par.nn[2] = 3;
            break;
        default:
            // do nothing
            break;
    }
    par.modelfile = model;
    Grid2D<double,uint32_t,sxz<double>> *g;
    if (string(model).find("vtr") != string::npos) {
        g = buildRectilinear2DfromVtr<double>(par, 1);
    } else {
        g = buildUnstructured2DfromVtu<double>(par, 1);
    }
    vector<vector<sxz<double>>> r_data;
    try {
        g->raytrace(src.get_coord(), src.get_t0(), rcv.get_coord(), rcv.get_tt(0), r_data);
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        abort();
    }
    string filename = "./files/" + get_class_name(g) + "_rp.vtp";
    saveRayPaths(filename, r_data);
    
    double error = get_rel_error(ref, rcv);
    BOOST_TEST_MESSAGE( "\t\t" << get_class_name(g) << " ttrp - error = " << error );
    
    BOOST_TEST(error < 0.15);
}


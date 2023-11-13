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

#include <Eigen/Sparse>

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

vector<double> get_rel_error(const string& filename, const Rcv2D<double>& rcv,
                             const vector<vector<siv<double>>> &l_data,
                             const vector<double> &slowness) {
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

    Eigen::SparseMatrix<double> A(l_data.size(), slowness.size());
    vector<Eigen::Triplet<double>> coefficients;
    for ( size_t i=0; i<l_data.size(); ++i ) {
        for ( size_t j=0; j<l_data[i].size(); ++j ) {
            coefficients.push_back(Eigen::Triplet<double>(static_cast<int>(i),
                                                          static_cast<int>(l_data[i][j].i),
                                                          l_data[i][j].v));
        }
    }
    A.setFromTriplets(coefficients.begin(), coefficients.end());
    A.makeCompressed();
    Eigen::VectorXd s(slowness.size());
    for ( size_t n=0; n<slowness.size(); ++n ) {
        s[n] = slowness[n];
    }
    Eigen::VectorXd tt = A * s;

    // compute relative error
    vector<double> error(2, 0.0);
    for ( size_t n=1; n<ref_tt.size(); ++n ) {
        // start at 1 to avoid node at source location where tt = 0
        error[0] += abs( (ref_tt[n] - tt[n])/ref_tt[n] );
        error[1] += abs( (rcv.get_tt(0)[n] - tt[n])/rcv.get_tt(0)[n] );
    }
    error[0] /= (ref_tt.size() - 1);
    error[1] /= (ref_tt.size() - 1);
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
const char* models_L[] = {
    "./files/layers_coarse2d.vtr",
    "./files/layers_coarse2d.vtu"
};
const char* models_aniso[] = {
    "./files/elliptical_fine2d.vtr",
    "./files/elliptical_fine2d.vtu",
    "./files/weakly_an_fine2d.vtr",
    "./files/weakly_an_fine2d.vtu"
};
const char* references[] = {
    "./files/sol_analytique_couches2d_tt.vtr",
    "./files/sol_analytique_gradient2d_tt.vtr",
    "./files/sol_analytique_couches2d_tt.vtr",
    "./files/sol_analytique_gradient2d_tt.vtr"
};
const char* references_L[] = {
    "./files/sol_analytique_couches2d_tt.vtr",
    "./files/sol_analytique_couches2d_tt.vtr"
};
const char* references_aniso[] = {
    "./files/sol_analytique_elliptical_2d_tt.vtr",
    "./files/sol_analytique_elliptical_2d_tt.vtr",
    "./files/sol_analytique_weakly_an_2d_tt.vtr",
    "./files/sol_analytique_weakly_an_2d_tt.vtr"
};
raytracing_method methods[] = { DYNAMIC_SHORTEST_PATH, FAST_SWEEPING, SHORTEST_PATH };
raytracing_method methods_sp[] = { SHORTEST_PATH };


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

BOOST_DATA_TEST_CASE(
                     testGrid2D_L,
                     (bdata::make(models_L) ^ bdata::make(references_L)) * bdata::make(methods),
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
    vector<vector<siv<double>>> l_data;
    try {
        g->raytrace(src.get_coord(), src.get_t0(), rcv.get_coord(), rcv.get_tt(0), l_data);
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        abort();
    }

    vector<double> slowness;
    g->getSlowness(slowness);
    vector<double> error = get_rel_error(ref, rcv, l_data, slowness);
    BOOST_TEST_MESSAGE( "\t\t" << get_class_name(g) << " ttrp - error = " << error[0] << ", " << error[1] );

    BOOST_TEST(error[0] < 0.15);
    BOOST_TEST(error[1] < 0.001);
}

BOOST_DATA_TEST_CASE(
                     testGrid_aniso2D,
                     (bdata::make(models_aniso) ^ bdata::make(references_aniso)) * bdata::make(methods_sp),
                     model, ref, method) {
    Src2D<double> src("./files/src2d.dat");
    src.init();
    Rcv2D<double> rcv("./files/rcv2daniso.dat");
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
            par.radius_tertiary_nodes = 200.0;
            par.nn[0] = 5;
            par.nn[1] = 5;
            par.nn[2] = 5;
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
    
    auto start = string(model).find("files/") + 6;
    auto end = string(model).find("_");
    string type_aniso = string(model).substr(start, end-start);
    string filename = "./files/" + get_class_name(g) + "_tt_grid_" + type_aniso;
    g->saveTT(filename, 0, 0, 2);
    double error = get_rel_error(ref, rcv);
    BOOST_TEST_MESSAGE( "\t\t" << get_class_name(g) << " - error = " << error );

    BOOST_TEST(error < 0.012);
}


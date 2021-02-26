//
//  main.cpp
//  test_rgrid2d
//
//  Created by Bernard Giroux on 2021-02-24.
//  Copyright © 2021 Bernard Giroux. All rights reserved.
//

#define BOOST_TEST_MODULE Test Grid2D
#include <boost/test/included/unit_test.hpp>

#include <cmath>
#include <string>
#include <vector>

#include <vtkPointData.h>
#include <vtkRectilinearGrid.h>
#include <vtkSmartPointer.h>
#include <vtkXMLRectilinearGridReader.h>

#include "Grid2D.h"
#include "Rcv2D.h"
#include "Src2D.h"
#include "structs_ttcr.h"
#include "grids.h"

namespace ttcr {
    int verbose = 2;
}

using namespace std;
using namespace ttcr;

void get_ref_tt(const string& filename, const vector<sxz<double>>& x, vector<double>& tt) {
    vtkSmartPointer<vtkXMLRectilinearGridReader> reader =
    vtkSmartPointer<vtkXMLRectilinearGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    reader->GetOutput()->Register(reader);
    vtkRectilinearGrid *dataSet = reader->GetOutput();
    
    vtkPointData *pd = dataSet->GetPointData();
    tt.resize(x.size());
    for ( size_t n=0; n<x.size(); ++n ) {
        vtkIdType ii = dataSet->FindPoint(x[n].x, 0.0, x[n].z);
        tt[n] = static_cast<double>(pd->GetArray(0)->GetTuple1(ii));
    }
}

BOOST_AUTO_TEST_CASE(testGrid2Drcfs)
{
    Src2D<double> src("./files/src2d.dat");
    src.init();
    Rcv2D<double> rcv("./files/rcv2d.dat");
    rcv.init(1);

    input_parameters par;
    par.method = FAST_SWEEPING;
    par.weno3 = 1;
    par.modelfile = "./files/layers_fine2d.vtr";

    Grid2D<double,uint32_t,sxz<double>> *g = buildRectilinear2DfromVtr<double>(par, 1);
    try {
        g->raytrace(src.get_coord(), src.get_t0(), rcv.get_coord(), rcv.get_tt(0));
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        abort();
    }
    g->saveTT("./files/grid2drcfs_tt_grid", 0, 0, 2);
    vector<double> ref_tt;
    get_ref_tt("./files/sol_analytique_couches2d_tt.vtr", rcv.get_coord(), ref_tt);

    // compute relative error
    double error = 0.0;
    for ( size_t n=1; n<ref_tt.size(); ++n ) {
        // start at 1 to avoid node at source location where tt = 0
        error += abs( (ref_tt[n] - rcv.get_tt(0)[n])/ref_tt[n] );
//        cout << rcv.get_coord()[n] << '\t' << ref_tt[n] << '\t' << rcv.get_tt(0)[n] << '\n';
    }
    error /= (ref_tt.size() - 1);
    
    BOOST_TEST(error < 0.03);
}

BOOST_AUTO_TEST_CASE(testGrid2Drcsp)
{
    Src2D<double> src("./files/src2d.dat");
    src.init();
    Rcv2D<double> rcv("./files/rcv2d.dat");
    rcv.init(1);

    input_parameters par;
    par.method = SHORTEST_PATH;
    par.nn[0] = 10;
    par.nn[1] = 10;
    par.nn[2] = 10;
    par.modelfile = "./files/layers_fine2d.vtr";

    Grid2D<double,uint32_t,sxz<double>> *g = buildRectilinear2DfromVtr<double>(par, 1);
    try {
        g->raytrace(src.get_coord(), src.get_t0(), rcv.get_coord(), rcv.get_tt(0));
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        abort();
    }
    g->saveTT("./files/grid2drcsp_tt_grid", 0, 0, 2);
    vector<double> ref_tt;
    get_ref_tt("./files/sol_analytique_couches2d_tt.vtr", rcv.get_coord(), ref_tt);
    
    // compute relative error
    double error = 0.0;
    for ( size_t n=1; n<ref_tt.size(); ++n ) {
        // start at 1 to avoid node at source location where tt = 0
        error += abs( (ref_tt[n] - rcv.get_tt(0)[n])/ref_tt[n] );
//        cout << rcv.get_coord()[n] << '\t' << ref_tt[n] << '\t' << rcv.get_tt(0)[n] << '\n';
    }
    error /= (ref_tt.size() - 1);
    
    BOOST_TEST(error < 0.02);
}

BOOST_AUTO_TEST_CASE(testGrid2Drnsp)
{
    Src2D<double> src("./files/src2d.dat");
    src.init();
    Rcv2D<double> rcv("./files/rcv2d.dat");
    rcv.init(1);

    input_parameters par;
    par.method = SHORTEST_PATH;
    par.nn[0] = 10;
    par.nn[1] = 10;
    par.nn[2] = 10;
    par.modelfile = "./files/gradient_fine2d.vtr";

    Grid2D<double,uint32_t,sxz<double>> *g = buildRectilinear2DfromVtr<double>(par, 1);
    try {
        g->raytrace(src.get_coord(), src.get_t0(), rcv.get_coord(), rcv.get_tt(0));
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        abort();
    }
    g->saveTT("./files/grid2drnsp_tt_grid", 0, 0, 2);
    vector<double> ref_tt;
    get_ref_tt("./files/sol_analytique_gradient2d_tt.vtr", rcv.get_coord(), ref_tt);

    // compute relative error
    double error = 0.0;
    for ( size_t n=1; n<ref_tt.size(); ++n ) {
        // start at 1 to avoid node at source location where tt = 0
        error += abs( (ref_tt[n] - rcv.get_tt(0)[n])/ref_tt[n] );
    }
    error /= (ref_tt.size() - 1);

    BOOST_TEST(error < 0.001);
}

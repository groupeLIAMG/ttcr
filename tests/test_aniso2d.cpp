//
//  test_aniso2d.cpp
//  ttcr
//
//  Created by Bernard Giroux on 2026-08-24.
//  Copyright © 2026 Bernard Giroux. All rights reserved.
//

#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <vector>

#define BOOST_TEST_MODULE Test Anisotropy 2D
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

#include "Cell.h"
#include "Grid2D.h"
#include "Rcv2D.h"
#include "Src2D.h"
#include "structs_ttcr.h"
#include "grids.h"
#include "utils.h"

namespace ttcr {
    int verbose = 0;
    int gpu_profile = 0;
}

namespace bdata = boost::unit_test::data;

using namespace std;
using namespace ttcr;

template<typename T>
string get_class_name(const Grid2D<T,uint32_t,sxz<T>> *g) {
    string name = typeid(*g).name();
    size_t start = name.find("Grid2D");
    if (name.find("OpenCL") != string::npos) {
        return name.substr(start, 17);
    }
    name = name.substr(start, 11);
    if ((name[8] == 'f' && name[9] == 's') ||
        (name[8] == 's' && name[9] == 'p') ) {
        name.pop_back();
    }
    return name;
}

template<typename T>
double get_rel_error(const string& filename, const Rcv2D<T>& rcv) {
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
        double x[3] = {rcv.get_coord()[n].x, 0.0, rcv.get_coord()[n].z};
        vtkIdType ii = dataSet->FindPoint(x);
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

const char* models_aniso[] = {
    "./files/elliptical_fine2d.vtr",
    "./files/elliptical_fine2d.vtu",
    "./files/weakly_an_fine2d.vtr",
    "./files/weakly_an_fine2d.vtu"
};
const char* references_aniso[] = {
    "./files/sol_analytique_elliptical_2d_tt.vtr",
    "./files/sol_analytique_elliptical_2d_tt.vtr",
    "./files/sol_analytique_weakly_an_2d_tt.vtr",
    "./files/sol_analytique_weakly_an_2d_tt.vtr"
};
raytracing_method methods_sp[] = { SHORTEST_PATH };


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

    BOOST_TEST(error < 0.01);
}


// ---------------------------------------------------------------------------
//  The cells on their own, without a grid.  Minimal stand-ins for the node and
//  coordinate types the cell classes are templated on.
// ---------------------------------------------------------------------------

namespace {

    struct Pt2 { double x, z; };

    struct Nd2 {
        double x, z;
        double getX() const { return x; }
        double getZ() const { return z; }
        double getDistance(const Pt2& o) const { return std::hypot(x-o.x, z-o.z); }
        double getDistance(const Nd2& o) const { return std::hypot(x-o.x, z-o.z); }
    };

    typedef CellElliptical<double,Nd2,Pt2>          CellEll;
    typedef CellTiltedElliptical<double,Nd2,Pt2>    CellTEll;
    typedef CellVTI_SH<double,Nd2,Pt2>              CellSH;
    typedef CellTTI_SH<double,Nd2,Pt2>              CellTSH;
    typedef CellWeaklyAnelliptical<double,Nd2,Pt2>  CellWA;
    typedef CellVTI_PSV<double,Nd2,Pt2>             CellPSV;
    typedef CellTTI_PSV<double,Nd2,Pt2>             CellTPSV;

    vector<double> one(double x) { return vector<double>(1, x); }

    // group velocity of the first arrival, obtained by locating every branch of
    // the slowness surface by bisection and keeping the fastest.  Unlike a value
    // read from the tables it varies smoothly with the medium parameters, which
    // a finite difference requires.
    double ray_angle(double th, double vp, double vs, double e, double d, double sgn) {
        double v, dv;
        VTI_PSV_GroupVel<double>::phaseVelocity(th, vp, vs, e, d, sgn, v, dv);
        return th + std::atan2(dv, v);
    }
    double group_vel(double th, double vp, double vs, double e, double d, double sgn) {
        double v, dv;
        VTI_PSV_GroupVel<double>::phaseVelocity(th, vp, vs, e, d, sgn, v, dv);
        return std::sqrt(v*v + dv*dv);
    }
    double first_arrival_vel(double psi, double vp, double vs,
                             double e, double d, double sgn) {
        const int N = 20000;
        double best = 0.0;
        double prev = ray_angle(-M_PI_2, vp, vs, e, d, sgn) - psi;
        for ( int k=1; k<=N; ++k ) {
            const double th = -M_PI_2 + M_PI*k/N;
            const double cur = ray_angle(th, vp, vs, e, d, sgn) - psi;
            if ( (prev < 0.) != (cur < 0.) ) {
                double lo = -M_PI_2 + M_PI*(k-1)/N, hi = th;
                for ( int it=0; it<80; ++it ) {
                    const double mid = 0.5*(lo+hi);
                    if ( (ray_angle(lo,vp,vs,e,d,sgn)-psi < 0.) !=
                         (ray_angle(mid,vp,vs,e,d,sgn)-psi < 0.) ) hi = mid;
                    else lo = mid;
                }
                best = std::fmax(best, group_vel(0.5*(lo+hi), vp, vs, e, d, sgn));
            }
            prev = cur;
        }
        return best;
    }
}

// The traveltime across a segment is its length divided by the group velocity
// in the direction of the segment.  Check the tables against the group velocity
// of the first arrival, computed independently.
BOOST_AUTO_TEST_CASE(testGroupVelocityTables) {
    const double vp=3.094, vs=1.51, eps=0.256, dlt=-0.0505;
    for ( int ph=0; ph<2; ++ph ) {
        const double sgn = ph ? -1.0 : 1.0;
        VTI_PSV_GroupVel<double> gv;
        gv.build(one(vp), one(vs), one(eps), one(dlt), sgn);
        double worst = 0.0;
        for ( int pd=0; pd<=90; pd+=3 ) {
            const double psi = pd*M_PI/180.;
            const double got = gv.velocity(std::sin(psi), std::cos(psi), 0);
            const double want = first_arrival_vel(psi, vp, vs, eps, dlt, sgn);
            worst = std::fmax(worst, std::fabs(got/want - 1.));
        }
        BOOST_TEST_MESSAGE("\t\t" << (ph?"qSV":"qP") << " group velocity - error = " << worst);
        BOOST_TEST(worst < 1.e-5);
    }
    // an isotropic medium must give back its velocity, whatever the direction
    VTI_PSV_GroupVel<double> iso;
    iso.build(one(2.0), one(1.0), one(0.0), one(0.0), 1.0);
    for ( int pd=0; pd<=90; pd+=5 ) {
        const double a = pd*M_PI/180.;
        BOOST_TEST(std::fabs(iso.velocity(std::sin(a), std::cos(a), 0) - 2.0) < 1.e-9);
    }
    // and an elliptical one the ellipse, epsilon and delta being equal
    const double E = 0.2;
    VTI_PSV_GroupVel<double> ell;
    ell.build(one(vp), one(vs), one(E), one(E), 1.0);
    for ( int pd=0; pd<=90; pd+=5 ) {
        const double a = pd*M_PI/180., lx = std::sin(a), lz = std::cos(a);
        const double ah = vp*std::sqrt(1.+2.*E);
        const double want = 1.0/std::sqrt(lx*lx/(ah*ah) + lz*lz/(vp*vp));
        BOOST_TEST(std::fabs(ell.velocity(lx, lz, 0)/want - 1.) < 1.e-4);
    }
    // tables are shared by the cells describing the same medium
    VTI_PSV_GroupVel<double> shared;
    vector<double> vpv(1000, vp), vsv(1000, vs), ev(1000, eps), dv(1000, dlt);
    ev[7] = 0.3;
    shared.build(vpv, vsv, ev, dv, -1.0);
    BOOST_TEST(shared.size() == 2u);
}

// Every cell class must report the derivatives of the traveltime with respect to
// its medium parameters.  Swept over both signs of the ray angle: folding the
// angle into the first quadrant reverses the sense in which it varies, and a
// sweep over positive angles alone hides a sign error in the tilt derivative.
BOOST_AUTO_TEST_CASE(testSensitivityAgainstFiniteDifferences) {
    Nd2 src{0.0, 0.0};
    double worst_closed = 0.0, worst_tab = 0.0;

    for ( int pd=-85; pd<=85; pd+=10 ) {
        const double a = pd*M_PI/180.;
        Pt2 dst{std::sin(a), std::cos(a)};
        const double h = 1.e-5;

        {   // elliptical: slowness and anisotropy ratio
            const double s0=0.5, x0=1.2;
            CellEll c(1); c.setSlowness(one(s0)); c.setXi(one(x0));
            siv2<double> d; d.i = 0; c.computeDistance(src, dst, d);
            CellEll a1(1); a1.setSlowness(one(s0+h)); a1.setXi(one(x0));
            CellEll a2(1); a2.setSlowness(one(s0-h)); a2.setXi(one(x0));
            worst_closed = std::fmax(worst_closed, std::fabs(d.v -
                (a1.computeDt(src,dst,0)-a2.computeDt(src,dst,0))/(2*h)));
            CellEll b1(1); b1.setSlowness(one(s0)); b1.setXi(one(x0+h));
            CellEll b2(1); b2.setSlowness(one(s0)); b2.setXi(one(x0-h));
            worst_closed = std::fmax(worst_closed, std::fabs(d.v2 -
                (b1.computeDt(src,dst,0)-b2.computeDt(src,dst,0))/(2*h)));
        }
        {   // VTI SH: vertical S-wave velocity and gamma
            const double v0=1.8, g0=0.15;
            CellSH c(1); c.setVs0(one(v0)); c.setGamma(one(g0));
            siv2<double> d; d.i = 0; c.computeDistance(src, dst, d);
            CellSH a1(1); a1.setVs0(one(v0+h)); a1.setGamma(one(g0));
            CellSH a2(1); a2.setVs0(one(v0-h)); a2.setGamma(one(g0));
            worst_closed = std::fmax(worst_closed, std::fabs(d.v -
                (a1.computeDt(src,dst,0)-a2.computeDt(src,dst,0))/(2*h)));
            CellSH b1(1); b1.setVs0(one(v0)); b1.setGamma(one(g0+h));
            CellSH b2(1); b2.setVs0(one(v0)); b2.setGamma(one(g0-h));
            worst_closed = std::fmax(worst_closed, std::fabs(d.v2 -
                (b1.computeDt(src,dst,0)-b2.computeDt(src,dst,0))/(2*h)));
        }
        {   // tilted elliptical: the tilt angle as well
            const double s0=0.5, x0=1.2, t0=0.4;
            CellTEll c(1); c.setSlowness(one(s0)); c.setXi(one(x0)); c.setTiltAngle(one(t0));
            siv4<double> d; d.i = 0; c.computeDistance(src, dst, d);
            CellTEll a1(1); a1.setSlowness(one(s0)); a1.setXi(one(x0)); a1.setTiltAngle(one(t0+h));
            CellTEll a2(1); a2.setSlowness(one(s0)); a2.setXi(one(x0)); a2.setTiltAngle(one(t0-h));
            worst_closed = std::fmax(worst_closed, std::fabs(d.v3 -
                (a1.computeDt(src,dst,0)-a2.computeDt(src,dst,0))/(2*h)));
        }
        {   // tilted SH
            const double v0=1.8, g0=0.15, t0=0.4;
            CellTSH c(1); c.setVs0(one(v0)); c.setGamma(one(g0)); c.setTiltAngle(one(t0));
            siv4<double> d; d.i = 0; c.computeDistance(src, dst, d);
            CellTSH a1(1); a1.setVs0(one(v0)); a1.setGamma(one(g0)); a1.setTiltAngle(one(t0+h));
            CellTSH a2(1); a2.setVs0(one(v0)); a2.setGamma(one(g0)); a2.setTiltAngle(one(t0-h));
            worst_closed = std::fmax(worst_closed, std::fabs(d.v3 -
                (a1.computeDt(src,dst,0)-a2.computeDt(src,dst,0))/(2*h)));
        }
        {   // weakly anelliptical: slowness and the two coefficients
            const double s0=0.5, a2c=0.05, a4c=0.01;
            CellWA c(1); c.setSlowness(one(s0)); c.setS2(one(a2c)); c.setS4(one(a4c));
            siv4<double> d; d.i = 0; c.computeDistance(src, dst, d);
            CellWA b1(1); b1.setSlowness(one(s0)); b1.setS2(one(a2c+h)); b1.setS4(one(a4c));
            CellWA b2(1); b2.setSlowness(one(s0)); b2.setS2(one(a2c-h)); b2.setS4(one(a4c));
            worst_closed = std::fmax(worst_closed, std::fabs(d.v2 -
                (b1.computeDt(src,dst,0)-b2.computeDt(src,dst,0))/(2*h)));
        }
        {   // the tilt of a coupled qP cell.  The tables are sampled every tenth
            // of a degree, so the step has to be much larger than that
            const double vp=3.094, vs=1.51, e0=0.256, d0=-0.0505, t0=0.35, hh=0.02;
            CellTPSV c(1);
            c.setVp0(one(vp)); c.setVs0(one(vs)); c.setEpsilon(one(e0));
            c.setDelta(one(d0)); c.setTiltAngle(one(t0)); c.setPhase(1);
            siv5<double> d; d.i = 0; c.computeDistance(src, dst, d);
            CellTPSV a1(1);
            a1.setVp0(one(vp)); a1.setVs0(one(vs)); a1.setEpsilon(one(e0));
            a1.setDelta(one(d0)); a1.setTiltAngle(one(t0+hh)); a1.setPhase(1);
            CellTPSV a2(1);
            a2.setVp0(one(vp)); a2.setVs0(one(vs)); a2.setEpsilon(one(e0));
            a2.setDelta(one(d0)); a2.setTiltAngle(one(t0-hh)); a2.setPhase(1);
            const double fd = (a1.computeDt(src,dst,0)-a2.computeDt(src,dst,0))/(2*hh);
            worst_tab = std::fmax(worst_tab,
                std::fabs(d.v5-fd)/std::fmax(std::fabs(fd), 1.e-4));
        }
    }
    BOOST_TEST_MESSAGE("\t\tclosed-form derivatives - error = " << worst_closed);
    BOOST_TEST_MESSAGE("\t\ttilt of a tabulated cell - error = " << worst_tab);
    BOOST_TEST(worst_closed < 1.e-7);
    BOOST_TEST(worst_tab < 5.e-2);
}

// Identities that hold exactly, and so need no finite difference: the traveltime
// is homogeneous of degree one in a slowness and of degree minus one in a
// velocity, the group velocity of a coupled cell scaling with both of its
// vertical velocities together.
BOOST_AUTO_TEST_CASE(testEulerIdentities) {
    Nd2 src{0.0, 0.0};
    Pt2 dst{0.6, 0.8};

    CellEll e(1); e.setSlowness(one(0.5)); e.setXi(one(1.2));
    siv2<double> d1; d1.i = 0; e.computeDistance(src, dst, d1);
    BOOST_TEST(std::fabs(0.5*d1.v/e.computeDt(src,dst,0) - 1.) < 1.e-12);

    CellSH h(1); h.setVs0(one(1.8)); h.setGamma(one(0.15));
    siv2<double> d2; d2.i = 0; h.computeDistance(src, dst, d2);
    BOOST_TEST(std::fabs(1.8*d2.v/h.computeDt(src,dst,0) + 1.) < 1.e-12);

    CellWA w(1); w.setSlowness(one(0.5)); w.setS2(one(0.05)); w.setS4(one(0.01));
    siv4<double> d3; d3.i = 0; w.computeDistance(src, dst, d3);
    BOOST_TEST(std::fabs(0.5*d3.v/w.computeDt(src,dst,0) - 1.) < 1.e-12);

    CellPSV p(1);
    p.setVp0(one(3.094)); p.setVs0(one(1.51));
    p.setEpsilon(one(0.256)); p.setDelta(one(-0.0505)); p.setPhase(0);
    siv4<double> d4; d4.i = 0; p.computeDistance(src, dst, d4);
    BOOST_TEST(std::fabs((3.094*d4.v + 1.51*d4.v2)/p.computeDt(src,dst,0) + 1.) < 1.e-12);
}

// A tilted medium is the same medium seen from a rotated frame.
BOOST_AUTO_TEST_CASE(testTiltedAgainstVertical) {
    const double vp=3.094, vs=1.51, eps=0.256, dlt=-0.0505;
    const double vs0=1.8, gam=0.15;
    Nd2 src{0.0, 0.0};

    // a zero tilt angle must reproduce the cells with a vertical symmetry axis
    CellTPSV tp(1);
    tp.setVp0(one(vp)); tp.setVs0(one(vs)); tp.setEpsilon(one(eps));
    tp.setDelta(one(dlt)); tp.setTiltAngle(one(0.0)); tp.setPhase(0);
    CellPSV vp_(1);
    vp_.setVp0(one(vp)); vp_.setVs0(one(vs)); vp_.setEpsilon(one(eps));
    vp_.setDelta(one(dlt)); vp_.setPhase(0);
    CellTSH th(1); th.setVs0(one(vs0)); th.setGamma(one(gam)); th.setTiltAngle(one(0.0));
    CellSH  vh(1); vh.setVs0(one(vs0)); vh.setGamma(one(gam));
    for ( int pd=-180; pd<=180; pd+=3 ) {
        const double a = pd*M_PI/180.;
        Pt2 dst{std::sin(a), std::cos(a)};
        BOOST_TEST(tp.computeDt(src,dst,0) == vp_.computeDt(src,dst,0));
        BOOST_TEST(th.computeDt(src,dst,0) == vh.computeDt(src,dst,0));
    }

    // and tilting by an angle must shift the traveltime pattern by that angle
    for ( double tilt : {20.0, 40.0, -35.0} ) {
        const double b = tilt*M_PI/180.;
        CellTPSV t(1);
        t.setVp0(one(vp)); t.setVs0(one(vs)); t.setEpsilon(one(eps));
        t.setDelta(one(dlt)); t.setTiltAngle(one(b)); t.setPhase(0);
        double worst = 0.0;
        for ( int pd=-180; pd<=180; pd+=3 ) {
            const double a = pd*M_PI/180.;
            Pt2 dt_{std::sin(a), std::cos(a)};
            Pt2 dv_{std::sin(a+b), std::cos(a+b)};
            worst = std::fmax(worst, std::fabs(t.computeDt(src,dt_,0) -
                                               vp_.computeDt(src,dv_,0)));
        }
        BOOST_TEST(worst < 1.e-12);
    }

    // the SH phase velocity is elliptical, so a tilted SH medium is a tilted
    // ellipse described through another cell class
    const double b = 30.*M_PI/180., xi = std::sqrt(1.+2.*gam);
    CellTSH sh(1); sh.setVs0(one(vs0)); sh.setGamma(one(gam)); sh.setTiltAngle(one(b));
    CellTEll el(1);
    el.setSlowness(one(1./(vs0*xi)));
    el.setXi(one(xi));                 // setXi squares its argument
    el.setTiltAngle(one(b));
    for ( int pd=-180; pd<=180; pd+=3 ) {
        const double a = pd*M_PI/180.;
        Pt2 dst{std::sin(a), std::cos(a)};
        BOOST_TEST(std::fabs(sh.computeDt(src,dst,0) - el.computeDt(src,dst,0)) < 1.e-12);
    }
}

// The structs carrying the sensitivities.
BOOST_AUTO_TEST_CASE(testSivStructs) {
    siv4<double> a;
    BOOST_TEST((a.i==0 && a.v==0. && a.v2==0. && a.v3==0. && a.v4==0.));
    siv5<double> b;
    BOOST_TEST((b.i==0 && b.v==0. && b.v2==0. && b.v3==0. && b.v4==0. && b.v5==0.));

    // the parent walk accumulates the cells it crosses more than once
    siv4<double> x; x.i=7; x.v=1; x.v2=2; x.v3=3; x.v4=4;
    siv4<double> y; y.i=99; y.v=10; y.v2=20; y.v3=30; y.v4=40;
    x += y;
    BOOST_TEST((x.i==7u && x.v==11. && x.v2==22. && x.v3==33. && x.v4==44.));

    siv5<double> p; p.i=3; p.v=1; p.v2=2; p.v3=3; p.v4=4; p.v5=5;
    siv5<double> q; q.i=88; q.v=.5; q.v2=.5; q.v3=.5; q.v4=.5; q.v5=.5;
    p += q;
    BOOST_TEST((p.i==3u && p.v==1.5 && p.v5==5.5));

    // the assembly of L needs the rows ordered by cell
    vector<siv4<double>> v4(3);
    v4[0].i=5; v4[1].i=1; v4[2].i=3;
    sort(v4.begin(), v4.end(), CompareSiv4_i<double>());
    BOOST_TEST((v4[0].i==1u && v4[1].i==3u && v4[2].i==5u));
    vector<siv5<double>> v5(3);
    v5[0].i=9; v5[1].i=2; v5[2].i=4;
    sort(v5.begin(), v5.end(), CompareSiv5_i<double>());
    BOOST_TEST((v5[0].i==2u && v5[1].i==4u && v5[2].i==9u));

    // the number of parameters each cell class describes
    BOOST_TEST(CellEll::nParams  == 2u);
    BOOST_TEST(CellTEll::nParams == 3u);
    BOOST_TEST(CellSH::nParams   == 2u);
    BOOST_TEST(CellTSH::nParams  == 3u);
    BOOST_TEST(CellWA::nParams   == 3u);
    BOOST_TEST(CellPSV::nParams  == 4u);
    BOOST_TEST(CellTPSV::nParams == 5u);
}

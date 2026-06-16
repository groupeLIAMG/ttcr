// Standalone harness for the Grid3Dun getRaypath pmr-allocator prototype.
//
// It overrides global operator new/delete to count heap allocations, then runs
// one raypath-producing raytrace on a NODE-slowness tetrahedral mesh (which
// routes to Grid3Dunfs -> Grid3Dun::getRaypath, the code the prototype changes).
//
// The field solve allocates identically in the baseline and pmr builds, so the
// build-to-build *difference* in total allocations is exactly the per-ray-step
// neighbour-set allocations the stack-backed std::pmr::set removes.  Allocation
// counts are deterministic (no timing noise); wall time is also reported but is
// dominated by the one-time fast-sweeping solve.
//
// Build (from repo root), once per header variant:
//   c++ -std=c++17 -O3 -DVTK -I ttcr -I eigen-5.0.0 -I boost_1_91_0 \
//       -I/usr/local/VTK-9.6.2/include/vtk-9.6 -L/usr/local/VTK-9.6.2/lib \
//       -lvtksys-9.6 -lvtkIOXML-9.6 -lvtkCommonDataModel-9.6 -lvtkCommonCore-9.6 \
//       tests/bench_raypath.cpp ttcr/ttcr_io.cpp -o bench_raypath
//
// Run from the repo root (relative data paths):
//   ./bench_raypath [model.vtu] [rx_replication] [n_threads]
// Defaults: tests/files/gradient_medium.vtu  40  1
//
// n_threads > 1 builds n_threads identical sources, each with its own replicated
// receiver set, and calls the multi-source raytrace, which the library spreads
// one-source-per-thread across an internal pool.  All threads build neighbour
// sets concurrently, so this exercises global-allocator lock contention -- the
// regime where removing 94% of allocations should help wall time, not just
// memory.  Compare the pmr vs baseline wall-time speedup at 1 vs N threads.

#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <new>
#include <string>
#include <vector>

#include "Grid3D.h"
#include "Rcv.h"
#include "Src.h"
#include "structs_ttcr.h"
#include "grids.h"

// ---- global allocation counters ----------------------------------------------
namespace {
    std::atomic<uint64_t> g_nalloc{0};
    std::atomic<uint64_t> g_nbytes{0};
    bool g_count = false;   // only count inside the timed region
}

void* operator new(std::size_t n) {
    if (g_count) { g_nalloc.fetch_add(1, std::memory_order_relaxed);
                   g_nbytes.fetch_add(n, std::memory_order_relaxed); }
    void* p = std::malloc(n ? n : 1);
    if (!p) throw std::bad_alloc();
    return p;
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }

using namespace std;
using namespace ttcr;

int main(int argc, char* argv[]) {

    cout << std::unitbuf;

    const string model    = argc > 1 ? argv[1] : "tests/files/gradient_medium.vtu";
    const size_t rxRep     = argc > 2 ? stoul(argv[2]) : 40;
    const size_t nThreads  = argc > 3 ? stoul(argv[3]) : 1;

    Src<double> src("tests/files/src3d_in.dat");
    src.init();
    Rcv<double> rcv("tests/files/rcv3d_in.dat");
    rcv.init(1);

    const vector<sxyz<double>>& rcv0 = rcv.get_coord();
    vector<sxyz<double>> Rx1;
    Rx1.reserve(rcv0.size() * rxRep);
    for (size_t r = 0; r < rxRep; ++r)
        Rx1.insert(Rx1.end(), rcv0.begin(), rcv0.end());

    input_parameters par;
    par.method     = FAST_SWEEPING;     // node slowness -> Grid3Dunfs
    par.modelfile  = model;
    par.tt_from_rp = true;
    Grid3D<double, uint32_t>* g = buildUnstructured3DfromVtu<double>(par, nThreads);
    if (g == nullptr) { cerr << "build failed: " << model << '\n'; return 1; }

    // One source per thread; each source gets the same replicated receiver set.
    vector<vector<sxyz<double>>> Tx(nThreads, src.get_coord());
    vector<vector<double>>       t0(nThreads, src.get_t0());
    vector<vector<sxyz<double>>> Rx(nThreads, Rx1);
    // The multi-source raytrace indexes traveltimes[n]/r_data[n] without
    // resizing the outer vectors, so they must be pre-sized to the source count.
    vector<vector<double>>               tt(nThreads);
    vector<vector<vector<sxyz<double>>>> r_data(nThreads);

    cout << "model         : " << model << '\n';
    cout << "nodes         : " << g->getNumberOfNodes() << '\n';
    cout << "threads       : " << nThreads << '\n';
    cout << "sources       : " << nThreads << "  (1 per thread)\n";
    cout << "raypaths/src  : " << Rx1.size() << "  (" << rcv0.size()
         << " x " << rxRep << ")\n";
    cout << "raypaths total: " << Rx1.size() * nThreads << "\n\n";

    // ---- timed + counted region: nThreads concurrent solve+raypath streams --
    g_nalloc = 0; g_nbytes = 0; g_count = true;
    auto t0c = chrono::high_resolution_clock::now();
    try {
        g->raytrace(Tx, t0, Rx, tt, r_data);
    } catch (std::exception& e) {
        g_count = false; cerr << "raytrace threw: " << e.what() << '\n'; return 1;
    }
    auto t1c = chrono::high_resolution_clock::now();
    g_count = false;

    size_t npts = 0; long double s = 0.0L;
    for (const auto& src_rd : r_data)
        for (const auto& ray : src_rd) { npts += ray.size();
            for (const auto& p : ray) s += p.x + p.y + p.z; }

    cout << "wall time     : " << chrono::duration<double, milli>(t1c - t0c).count() << " ms\n";
    cout << "allocations   : " << g_nalloc.load() << '\n';
    cout << "alloc bytes   : " << g_nbytes.load() << '\n';
    cout << "raypath pts   : " << npts << '\n';
    cout << "coord cksum   : " << static_cast<double>(s) << '\n';

    delete g;
    return 0;
}

//
//  Grid2Drn_OpenCL.h
//  ttcr
//
//  OpenCL host wrapper for the 2D fast-sweeping kernels.
//  Mirrors Grid3Drn_OpenCL.h but for 2D grids:
//    - 4 sweep directions (±i, ±j) instead of 8
//    - node index i*(ncz+1)+j  (x outer, z inner)
//    - 4 kernel handles (basic, xz, weno3, weno3_xz)
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
 */

#ifndef ttcr_Grid2Drn_OpenCL_h
#define ttcr_Grid2Drn_OpenCL_h

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#include "Grid2Drn_kernels_src.h"

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "ttcr_t.h"   // extern int verbose, gpu_profile

namespace ttcr {

enum class SweepMode2D {
    BASIC,        ///< 1st-order, square cells (dx == dz)
    BASIC_XZ,     ///< 1st-order, non-square cells
    WENO3,        ///< WENO3, square cells
    WENO3_XZ      ///< WENO3, non-square cells
};

template<typename T1>
class OpenCLSweepSolver2D {
public:
    OpenCLSweepSolver2D() :
        platform(nullptr), device(nullptr), context(nullptr),
        queue(nullptr), program(nullptr),
        kernel_basic(nullptr), kernel_xz(nullptr),
        kernel_weno3(nullptr), kernel_weno3_xz(nullptr),
        d_tt_in(nullptr), d_tt_out(nullptr),
        d_slowness(nullptr), d_frozen(nullptr), d_plane_nodes(nullptr),
        maxlevel(0),
        ncx(0), ncz(0), dx(T1(0)), dz(T1(0)), num_nodes(0),
        initialized(false),
        current_mode(SweepMode2D::BASIC),
        optimal_local_size(64),
        profiling(false),
        prof_kernel_ns(0), prof_mirror_ns(0),
        prof_upload_ns(0), prof_download_ns(0),
        prof_launches(0), prof_wall_ms(0.0)
    {}

    ~OpenCLSweepSolver2D() { cleanup(); }

    // =========================================================================
    // Public interface
    // =========================================================================

    void initialize(size_t nx, size_t nz, T1 grid_dx, T1 grid_dz) {
        static_assert(std::is_floating_point<T1>::value &&
                      (sizeof(T1) == 4 || sizeof(T1) == 8),
                      "OpenCLSweepSolver2D requires T1 to be float or double");
        ncx       = nx;
        ncz       = nz;
        dx        = grid_dx;
        dz        = grid_dz;
        num_nodes = (ncx + 1) * (ncz + 1);

        setupSweepDirections();
        initializeOpenCL();
        buildKernels();
        allocateDeviceMemory();
        determineOptimalWorkGroupSize();
        buildPlaneNodeLists();
        setStaticKernelArgs();
        initialized = true;
    }

    void setSweepMode(SweepMode2D mode) { current_mode = mode; }
    void setProfiling(bool on)   { profiling = on; }
    bool isProfiling() const     { return profiling; }

    /// Run complete sweep cycles with full data upload and download.
    void runSweeps(std::vector<T1>& tt,
                   const std::vector<T1>& slowness,
                   const std::vector<bool>& frozen,
                   size_t num_sweeps = 1)
    {
        if (!initialized) throw std::runtime_error("2D OpenCL solver not initialized");
        if (tt.size() != num_nodes || slowness.size() != num_nodes ||
            frozen.size() != num_nodes)
            throw std::runtime_error("Array size mismatch in runSweeps2D");

        std::vector<unsigned char> frozen_u(num_nodes);
        for (size_t n = 0; n < num_nodes; ++n) frozen_u[n] = frozen[n] ? 1 : 0;

        uploadData(tt, slowness, frozen_u);
        for (size_t c = 0; c < num_sweeps; ++c) performSweepCycle();
        downloadData(tt);
        clFinish(queue);
    }

    /// Run sweep cycles without transferring data (iterative convergence loop).
    void runSweepsNoTransfer(size_t num_sweeps = 1) {
        if (!initialized) throw std::runtime_error("2D OpenCL solver not initialized");
        for (size_t c = 0; c < num_sweeps; ++c) performSweepCycle();
        clFinish(queue);
    }

    /// Upload travel times (for iterative loops after initial runSweeps).
    void uploadTravelTimes(const std::vector<T1>& tt) {
        if (!initialized) throw std::runtime_error("2D OpenCL solver not initialized");
        if (tt.size() != num_nodes) throw std::runtime_error("Array size mismatch");
        cl_event evt = nullptr;
        cl_event* ep = profiling ? &evt : nullptr;
        cl_int err = clEnqueueWriteBuffer(queue, d_tt_in, CL_TRUE, 0,
                                          num_nodes * sizeof(T1), tt.data(),
                                          0, nullptr, ep);
        checkError(err, "Uploading tt to d_tt_in");
        addEventTime(evt, prof_upload_ns);
        err = clEnqueueWriteBuffer(queue, d_tt_out, CL_TRUE, 0,
                                   num_nodes * sizeof(T1), tt.data(),
                                   0, nullptr, ep);
        checkError(err, "Uploading tt to d_tt_out");
        addEventTime(evt, prof_upload_ns);
    }

    /// Download the latest travel times from the GPU.
    void downloadTravelTimes(std::vector<T1>& tt) {
        if (!initialized) throw std::runtime_error("2D OpenCL solver not initialized");
        if (tt.size() != num_nodes) throw std::runtime_error("Array size mismatch");
        downloadData(tt);
    }

    std::string getDeviceInfo() const {
        if (!device) return "No device initialized";
        char   name[256], vendor[256];
        cl_ulong gmem, lmem;
        size_t   mwg;
        clGetDeviceInfo(device, CL_DEVICE_NAME,            sizeof(name),   name,   nullptr);
        clGetDeviceInfo(device, CL_DEVICE_VENDOR,          sizeof(vendor), vendor, nullptr);
        clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(gmem),  &gmem,  nullptr);
        clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE,  sizeof(lmem),  &lmem,  nullptr);
        clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(mwg), &mwg, nullptr);
        std::ostringstream oss;
        oss << "    Device: " << name   << "\n"
            << "    Vendor: " << vendor << "\n"
            << "    Global Memory: " << (gmem / (1024*1024)) << " MB\n"
            << "    Local Memory: "  << (lmem / 1024)        << " KB\n"
            << "    Max Work Group: " << mwg                 << "\n";
        return oss.str();
    }

    void reportProfile(std::ostream& os = std::cout) const {
        if (prof_launches == 0) {
            os << "  [OpenCL2D profile] no kernel launches recorded\n";
            return;
        }
        const double km  = prof_kernel_ns   * 1e-6;
        const double mm  = prof_mirror_ns   * 1e-6;
        const double upm = prof_upload_ns   * 1e-6;
        const double dnm = prof_download_ns * 1e-6;
        const double ovm = prof_wall_ms - km - mm;
        auto frac = [&](double ms){ return prof_wall_ms > 0.0 ? 100.0*ms/prof_wall_ms : 0.0; };
        os << "  [OpenCL2D profile]\n"
           << "      sweep wall   : " << prof_wall_ms << " ms\n"
           << "      kernel exec  : " << km  << " ms  (" << frac(km)  << "%)\n"
           << "      mirror copy  : " << mm  << " ms  (" << frac(mm)  << "%)\n"
           << "      other/overhd : " << ovm << " ms  (" << frac(ovm) << "%)\n"
           << "      upload       : " << upm << " ms\n"
           << "      download     : " << dnm << " ms\n"
           << "      launches     : " << prof_launches << "\n";
    }

private:
    // =========================================================================
    // OpenCL objects
    // =========================================================================
    cl_platform_id  platform;
    cl_device_id    device;
    cl_context      context;
    cl_command_queue queue;
    cl_program      program;
    cl_kernel       kernel_basic, kernel_xz, kernel_weno3, kernel_weno3_xz;
    cl_mem          d_tt_in, d_tt_out, d_slowness, d_frozen, d_plane_nodes;

    // Per-(direction, level) plane-node index lists.
    // 4 directions × num_nodes uints in d_plane_nodes.
    // level_offsets[d][L] is the start of level L in direction d's block.
    std::vector<std::vector<cl_uint>> level_offsets;
    int maxlevel;

    size_t ncx, ncz, num_nodes;
    T1     dx,  dz;
    bool   initialized;

    SweepMode2D current_mode;
    size_t      optimal_local_size;

    struct SweepDir2D { int i_dir, j_dir; };
    std::vector<SweepDir2D> sweep_directions;

    // Profiling
    bool     profiling;
    cl_ulong prof_kernel_ns, prof_mirror_ns, prof_upload_ns, prof_download_ns;
    size_t   prof_launches;
    double   prof_wall_ms;
    std::vector<cl_event> prof_events;

    // =========================================================================
    // Setup
    // =========================================================================

    void setupSweepDirections() {
        sweep_directions = {
            { +1, +1 },   // direction 0: i+, j+
            { -1, +1 },   // direction 1: i-, j+
            { -1, -1 },   // direction 2: i-, j-
            { +1, -1 },   // direction 3: i+, j-
        };
    }

    void initializeOpenCL() {
        cl_int err;
        err = clGetPlatformIDs(1, &platform, nullptr);
        checkError(err, "Getting platform");

        err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, nullptr);
        if (err != CL_SUCCESS) {
            if ( verbose ) std::cout << "No GPU found, trying CPU...\n";
            err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 1, &device, nullptr);
        }
        checkError(err, "Getting device");

        context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &err);
        checkError(err, "Creating context");

#ifdef CL_VERSION_2_0
        const cl_queue_properties qprops[] =
            { CL_QUEUE_PROPERTIES, CL_QUEUE_PROFILING_ENABLE, 0 };
        queue = clCreateCommandQueueWithProperties(context, device, qprops, &err);
#else
        queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &err);
#endif
        checkError(err, "Creating command queue");
    }

    void buildKernels() {
        // On-disk overrides for development; fall back to embedded source.
        const char* search[] = {
            "Grid2Drn_kernels.cl",
            "../Grid2Drn_kernels.cl",
            "../../Grid2Drn_kernels.cl",
            "./kernels/Grid2Drn_kernels.cl",
        };

        std::string src;
        bool found = false;
        for (auto& p : search) {
            std::ifstream f(p);
            if (f) {
                std::ostringstream ss;
                ss << f.rdbuf();
                src   = ss.str();
                found = true;
                break;
            }
        }
        if (!found) src = Grid2Drn_kernels_src;

        const char* cstr = src.c_str();
        size_t      len  = src.size();
        cl_int err;
        program = clCreateProgramWithSource(context, 1, &cstr, &len, &err);
        checkError(err, "Creating 2D program");

        std::string opts = "-cl-fast-relaxed-math";
        if (std::is_same<T1, double>::value) {
            cl_device_fp_config fp64 = 0;
            clGetDeviceInfo(device, CL_DEVICE_DOUBLE_FP_CONFIG,
                            sizeof(fp64), &fp64, nullptr);
            if (fp64 == 0)
                throw std::runtime_error(
                    "double requested but device has no cl_khr_fp64 support");
            opts += " -DUSE_DOUBLE=1";
        } else {
            opts += " -cl-single-precision-constant";
        }

        err = clBuildProgram(program, 1, &device, opts.c_str(), nullptr, nullptr);
        if (err != CL_SUCCESS) {
            size_t sz;
            clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &sz);
            std::vector<char> log(sz);
            clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sz, log.data(), nullptr);
            std::cerr << "2D kernel build failed:\n" << log.data() << "\n";
            throw std::runtime_error("2D kernel build failed");
        }

        kernel_basic     = clCreateKernel(program, "sweep_update_basic_2d",    &err);
        checkError(err, "Creating kernel_basic");
        kernel_xz        = clCreateKernel(program, "sweep_update_xz_2d",       &err);
        checkError(err, "Creating kernel_xz");
        kernel_weno3     = clCreateKernel(program, "sweep_update_weno3_2d",    &err);
        checkError(err, "Creating kernel_weno3");
        kernel_weno3_xz  = clCreateKernel(program, "sweep_update_weno3_xz_2d", &err);
        checkError(err, "Creating kernel_weno3_xz");
    }

    void allocateDeviceMemory() {
        cl_int  err;
        size_t  ms = num_nodes * sizeof(T1);
        size_t  fs = num_nodes * sizeof(unsigned char);

        d_tt_in    = clCreateBuffer(context, CL_MEM_READ_WRITE, ms, nullptr, &err);
        checkError(err, "Allocating d_tt_in");
        d_tt_out   = clCreateBuffer(context, CL_MEM_READ_WRITE, ms, nullptr, &err);
        checkError(err, "Allocating d_tt_out");
        d_slowness = clCreateBuffer(context, CL_MEM_READ_ONLY,  ms, nullptr, &err);
        checkError(err, "Allocating d_slowness");
        d_frozen   = clCreateBuffer(context, CL_MEM_READ_ONLY,  fs, nullptr, &err);
        checkError(err, "Allocating d_frozen");

        if ( verbose )
            std::cout << "\n  2D GPU: allocated "
                      << ((2*ms + ms + fs) / (1024*1024)) << " MB device memory\n";
    }

    void determineOptimalWorkGroupSize() {
        cl_int err = CL_SUCCESS;

        size_t device_max = 0;
        err = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE,
                              sizeof(device_max), &device_max, nullptr);
        checkError(err, "Querying max work group size");

        // Respect kernel-specific work-group limits (can be lower than the device maximum).
        size_t max_wg = device_max;
        for (cl_kernel k : { kernel_basic, kernel_xz, kernel_weno3, kernel_weno3_xz }) {
            size_t km = 0;
            if (clGetKernelWorkGroupInfo(k, device, CL_KERNEL_WORK_GROUP_SIZE,
                                         sizeof(km), &km, nullptr) == CL_SUCCESS && km > 0) {
                max_wg = std::min(max_wg, km);
            }
        }

        // Query preferred multiple from one of our kernels.
        size_t pref = 32;
#ifdef CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE
        size_t kp = 0;
        if (clGetKernelWorkGroupInfo(kernel_basic, device,
                                    CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
                                    sizeof(kp), &kp, nullptr) == CL_SUCCESS && kp > 0) {
            pref = kp;
        }
#endif

        // For 1D launches: 256 is a good default across vendors.
        optimal_local_size = std::min(size_t(256), max_wg);

        // Round down to preferred multiple (and keep within max_wg).
        if (pref > 1) {
            optimal_local_size = (optimal_local_size / pref) * pref;
        }
        if (optimal_local_size == 0) {
            optimal_local_size = std::min(pref, max_wg);
        }
        if (optimal_local_size == 0) {
            optimal_local_size = 1;
        }

        if ( verbose )
            std::cout << "  2D GPU work group size: " << optimal_local_size << "\n";
    }

    /// Precompute per-(direction, level) node-index lists.
    /// Direction d's block occupies [d*num_nodes, (d+1)*num_nodes) in the
    /// flat array; level_offsets[d][L] gives the start of level L within that.
    void buildPlaneNodeLists() {
        maxlevel = static_cast<int>(ncx + ncz);
        const int    nlev = maxlevel + 1;
        const size_t nx1  = ncx + 1;
        const size_t nz1  = ncz + 1;

        level_offsets.assign(4, std::vector<cl_uint>(nlev + 1, 0));
        std::vector<cl_uint> all_perm(4 * num_nodes);

        for (size_t d = 0; d < 4; ++d) {
            const auto&           dir = sweep_directions[d];
            std::vector<cl_uint>& off = level_offsets[d];

            std::vector<cl_uint> count(nlev, 0);
            for (size_t i = 0; i < nx1; ++i)
                for (size_t j = 0; j < nz1; ++j) {
                    const int ip = dir.i_dir > 0 ? (int)i : (int)(ncx - i);
                    const int jp = dir.j_dir > 0 ? (int)j : (int)(ncz - j);
                    ++count[ip + jp];
                }

            cl_uint run = 0;
            for (int L = 0; L < nlev; ++L) { off[L] = run; run += count[L]; }
            off[nlev] = run;

            std::vector<cl_uint> cursor(off.begin(), off.begin() + nlev);
            for (size_t i = 0; i < nx1; ++i)
                for (size_t j = 0; j < nz1; ++j) {
                    const int    ip  = dir.i_dir > 0 ? (int)i : (int)(ncx - i);
                    const int    jp  = dir.j_dir > 0 ? (int)j : (int)(ncz - j);
                    const size_t idx = i * nz1 + j;   // i*(ncz+1)+j
                    all_perm[d * num_nodes + cursor[ip + jp]++] = static_cast<cl_uint>(idx);
                }
        }

        cl_int err;
        d_plane_nodes = clCreateBuffer(context, CL_MEM_READ_ONLY,
                                       4 * num_nodes * sizeof(cl_uint), nullptr, &err);
        checkError(err, "Allocating 2D plane-node buffer");
        err = clEnqueueWriteBuffer(queue, d_plane_nodes, CL_TRUE, 0,
                                   4 * num_nodes * sizeof(cl_uint),
                                   all_perm.data(), 0, nullptr, nullptr);
        checkError(err, "Uploading 2D plane-node lists");
    }

    /// Bind the 9 arguments that never change; only plane_offset (arg 9) and
    /// plane_count (arg 10) vary per launch and are set in executeSweep().
    void setStaticKernelArgs() {
        const cl_uint ncx_u = static_cast<cl_uint>(ncx);
        const cl_uint ncz_u = static_cast<cl_uint>(ncz);
        for (cl_kernel k : { kernel_basic, kernel_xz, kernel_weno3, kernel_weno3_xz }) {
            cl_int  err = 0;
            cl_uint a   = 0;
            err |= clSetKernelArg(k, a++, sizeof(cl_mem), &d_tt_in);    // tt_in
            err |= clSetKernelArg(k, a++, sizeof(cl_mem), &d_tt_in);    // tt_out (same: Gauss-Seidel in-place)
            err |= clSetKernelArg(k, a++, sizeof(cl_mem), &d_slowness);
            err |= clSetKernelArg(k, a++, sizeof(cl_mem), &d_frozen);
            err |= clSetKernelArg(k, a++, sizeof(T1),     &dx);
            err |= clSetKernelArg(k, a++, sizeof(T1),     &dz);
            err |= clSetKernelArg(k, a++, sizeof(cl_uint), &ncx_u);
            err |= clSetKernelArg(k, a++, sizeof(cl_uint), &ncz_u);
            err |= clSetKernelArg(k, a++, sizeof(cl_mem), &d_plane_nodes);
            checkError(err, "Setting static 2D kernel args");
        }
    }

    // =========================================================================
    // Data transfer
    // =========================================================================

    void uploadData(const std::vector<T1>& tt,
                    const std::vector<T1>& slowness,
                    const std::vector<unsigned char>& frozen) {
        cl_int  err;
        cl_event evt = nullptr;
        cl_event* ep = profiling ? &evt : nullptr;

        err = clEnqueueWriteBuffer(queue, d_tt_in, CL_TRUE, 0,
                                   num_nodes*sizeof(T1), tt.data(), 0, nullptr, ep);
        checkError(err, "Uploading tt_in");
        addEventTime(evt, prof_upload_ns);

        err = clEnqueueWriteBuffer(queue, d_tt_out, CL_TRUE, 0,
                                   num_nodes*sizeof(T1), tt.data(), 0, nullptr, ep);
        checkError(err, "Uploading tt_out");
        addEventTime(evt, prof_upload_ns);

        err = clEnqueueWriteBuffer(queue, d_slowness, CL_TRUE, 0,
                                   num_nodes*sizeof(T1), slowness.data(), 0, nullptr, ep);
        checkError(err, "Uploading slowness");
        addEventTime(evt, prof_upload_ns);

        err = clEnqueueWriteBuffer(queue, d_frozen, CL_TRUE, 0,
                                   num_nodes*sizeof(unsigned char), frozen.data(),
                                   0, nullptr, ep);
        checkError(err, "Uploading frozen");
        addEventTime(evt, prof_upload_ns);
    }

    void downloadData(std::vector<T1>& tt) {
        cl_event  evt = nullptr;
        cl_event* ep  = profiling ? &evt : nullptr;
        cl_int err = clEnqueueReadBuffer(queue, d_tt_out, CL_TRUE, 0,
                                         num_nodes*sizeof(T1), tt.data(),
                                         0, nullptr, ep);
        checkError(err, "Downloading traveltimes");
        addEventTime(evt, prof_download_ns);
    }

    // =========================================================================
    // Sweep execution
    // =========================================================================

    cl_kernel selectKernel() const {
        switch (current_mode) {
            case SweepMode2D::BASIC:     return kernel_basic;
            case SweepMode2D::BASIC_XZ:  return kernel_xz;
            case SweepMode2D::WENO3:     return kernel_weno3;
            case SweepMode2D::WENO3_XZ:  return kernel_weno3_xz;
        }
        return kernel_basic;
    }

    void performSweepCycle() {
        cl_kernel cur = selectKernel();

        std::chrono::high_resolution_clock::time_point t0;
        if (profiling) t0 = std::chrono::high_resolution_clock::now();

        for (size_t dir = 0; dir < 4; ++dir)
            for (int level = 0; level <= maxlevel; ++level)
                executeSweep(cur, dir, level);

        // Mirror in-place result to d_tt_out so downloadData() is consistent.
        cl_event mevt = nullptr;
        cl_event* mep = profiling ? &mevt : nullptr;
        cl_int err = clEnqueueCopyBuffer(queue, d_tt_in, d_tt_out, 0, 0,
                                          num_nodes * sizeof(T1), 0, nullptr, mep);
        checkError(err, "Mirroring 2D GS result to d_tt_out");
        clFinish(queue);

        if (profiling) {
            auto t1 = std::chrono::high_resolution_clock::now();
            prof_wall_ms += std::chrono::duration<double,std::milli>(t1-t0).count();
            for (cl_event e : prof_events) {
                prof_kernel_ns += eventNs(e);
                clReleaseEvent(e);
            }
            prof_events.clear();
            if (mevt) { prof_mirror_ns += eventNs(mevt); clReleaseEvent(mevt); }
        }
    }

    void executeSweep(cl_kernel kernel, size_t dir, int level) {
        const cl_uint lstart = level_offsets[dir][level];
        const cl_uint lend   = level_offsets[dir][level + 1];
        const cl_uint count  = lend - lstart;
        if (count == 0) return;

        const cl_uint offset = static_cast<cl_uint>(dir * num_nodes) + lstart;

        static const cl_uint OFFSET_ARG = 9;
        static const cl_uint COUNT_ARG  = 10;
        cl_int err = 0;
        err |= clSetKernelArg(kernel, OFFSET_ARG, sizeof(cl_uint), &offset);
        err |= clSetKernelArg(kernel, COUNT_ARG,  sizeof(cl_uint), &count);
        checkError(err, "Setting per-plane 2D kernel args");

        size_t ls = optimal_local_size;
        if (ls < 1) ls = 1;
        size_t gs = ((count + ls - 1) / ls) * ls;

        cl_event  evt = nullptr;
        cl_event* ep  = profiling ? &evt : nullptr;
        err = clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &gs, &ls,
                                     0, nullptr, ep);
        checkError(err, "Launching 2D sweep kernel");
        if (profiling) { prof_events.push_back(evt); ++prof_launches; }
    }

    // =========================================================================
    // Utilities
    // =========================================================================

    void checkError(cl_int err, const std::string& op) {
        if (err != CL_SUCCESS) {
            std::ostringstream oss;
            oss << "OpenCL error during " << op << ": " << err;
            throw std::runtime_error(oss.str());
        }
    }

    static cl_ulong eventNs(cl_event evt) {
        cl_ulong s = 0, e = 0;
        clGetEventProfilingInfo(evt, CL_PROFILING_COMMAND_START, sizeof(s), &s, nullptr);
        clGetEventProfilingInfo(evt, CL_PROFILING_COMMAND_END,   sizeof(e), &e, nullptr);
        return (e > s) ? (e - s) : 0;
    }

    void addEventTime(cl_event evt, cl_ulong& acc) {
        if (!profiling || !evt) return;
        acc += eventNs(evt);
        clReleaseEvent(evt);
    }

    void cleanup() {
        if (profiling && prof_launches > 0) reportProfile(std::cout);
        for (cl_event e : prof_events) clReleaseEvent(e);
        prof_events.clear();

        if (d_tt_in)       clReleaseMemObject(d_tt_in);
        if (d_tt_out)      clReleaseMemObject(d_tt_out);
        if (d_slowness)    clReleaseMemObject(d_slowness);
        if (d_frozen)      clReleaseMemObject(d_frozen);
        if (d_plane_nodes) clReleaseMemObject(d_plane_nodes);

        if (kernel_basic)    clReleaseKernel(kernel_basic);
        if (kernel_xz)       clReleaseKernel(kernel_xz);
        if (kernel_weno3)    clReleaseKernel(kernel_weno3);
        if (kernel_weno3_xz) clReleaseKernel(kernel_weno3_xz);

        if (program) clReleaseProgram(program);
        if (queue)   clReleaseCommandQueue(queue);
        if (context) clReleaseContext(context);
    }
};

} // namespace ttcr

#endif // ttcr_Grid2Drn_OpenCL_h

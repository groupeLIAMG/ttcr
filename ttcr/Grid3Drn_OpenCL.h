/*
 * OpenCL Host Code for Grid3Drn Sweep Methods
 * 
 * This file provides C++ wrapper classes to execute the sweep and sweep_weno3
 * algorithms on GPU using OpenCL.
 * 
 * Usage:
 *   OpenCLSweepSolver solver;
 *   solver.initialize(ncx, ncy, ncz, dx, dy, dz);
 *   solver.setSweepType(SweepType::WENO3);
 *   solver.runSweeps(tt_data, slowness, frozen);
 */

#ifndef GRID3DRN_OPENCL_H
#define GRID3DRN_OPENCL_H

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

// Embedded OpenCL kernel source (auto-generated from Grid3Drn_kernels.cl by
// setup.py).  Provides ttcr::Grid3Drn_kernels_src so an installed wheel can
// build the kernels without locating the .cl file on disk.
#include "Grid3Drn_kernels_src.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <type_traits>
#include <chrono>

namespace ttcr {

// =============================================================================
// Sweep Type Enumeration
// =============================================================================

enum class SweepType {
    BASIC,      // First-order sweep
    WENO3       // Third-order WENO sweep
};

// =============================================================================
// OpenCL Sweep Solver Class
// =============================================================================

template<typename T1>
class OpenCLSweepSolver {
private:
    // OpenCL context and resources
    cl_platform_id platform;
    cl_device_id device;
    cl_context context;
    cl_command_queue queue;
    cl_program program;
    
    // Kernels
    cl_kernel kernel_basic;
    cl_kernel kernel_weno3;
    
    // Device buffers
    cl_mem d_tt_in;
    cl_mem d_tt_out;
    cl_mem d_slowness;
    cl_mem d_frozen;

    // Precomputed per-(direction,level) plane-node lists for the Gauss-Seidel
    // sweep.  d_plane_nodes holds the 8 direction permutations concatenated
    // (size 8*num_nodes uints): for direction d, the nodes of level L occupy
    // d_plane_nodes[d*num_nodes + level_offsets[d][L] .. d*num_nodes +
    // level_offsets[d][L+1]).  Built once in buildPlaneNodeLists() and reused
    // every sweep cycle, so each plane is dispatched as a 1-D range over
    // exactly its nodes instead of the whole grid box (see executeSweep).
    cl_mem d_plane_nodes;
    std::vector<std::vector<cl_uint>> level_offsets;  // [8][maxlevel+2]
    int maxlevel;

    // Grid parameters
    size_t ncx, ncy, ncz;
    T1 dx, dy, dz;
    size_t num_nodes;
    bool initialized;

    SweepType current_sweep_type;
    
    // Optimal work group size (determined during initialization)
    size_t optimal_local_size[3];
    
    // Sweep directions (8 combinations of ±1)
    struct SweepDirection {
        int i_dir, j_dir, k_dir;
        size_t i_start, j_start, k_start;
    };
    
    std::vector<SweepDirection> sweep_directions;

    // --- Profiling (opt-in, see setProfiling) -------------------------------
    // When enabled, kernel launches and host<->device transfers are timed using
    // OpenCL events (the queue is created with CL_QUEUE_PROFILING_ENABLE), and
    // the sweep cycles are wall-clock timed.  The breakdown is dumped at
    // cleanup() so a profiling run prints one summary at the end.  Disabled by
    // default, with zero added cost on the hot path (no events are attached).
    bool profiling;
    cl_ulong prof_kernel_ns;    // summed device execution time of sweep kernels
    cl_ulong prof_mirror_ns;    // summed device time of the per-cycle d_tt_in->d_tt_out copy
    cl_ulong prof_upload_ns;    // summed device time of host->device writes
    cl_ulong prof_download_ns;  // summed device time of device->host reads
    size_t   prof_launches;     // number of kernel launches (plane dispatches)
    double   prof_wall_ms;      // summed wall-clock time of performSweepCycle()
    std::vector<cl_event> prof_events;  // pending kernel events for the cycle

public:
    OpenCLSweepSolver() : 
        platform(nullptr), device(nullptr), context(nullptr), 
        queue(nullptr), program(nullptr),
        kernel_basic(nullptr), kernel_weno3(nullptr),
        d_tt_in(nullptr), d_tt_out(nullptr), d_slowness(nullptr), d_frozen(nullptr),
        d_plane_nodes(nullptr), maxlevel(0),
        ncx(0), ncy(0), ncz(0), num_nodes(0), initialized(false),
        current_sweep_type(SweepType::BASIC),
        profiling(false), prof_kernel_ns(0), prof_mirror_ns(0), prof_upload_ns(0),
        prof_download_ns(0), prof_launches(0), prof_wall_ms(0.0)
    {
    }
    
    ~OpenCLSweepSolver() {
        cleanup();
    }
    
    // =========================================================================
    // Initialization
    // =========================================================================
    
    void initialize(size_t nx, size_t ny, size_t nz,
                   T1 grid_dx, T1 grid_dy, T1 grid_dz) {
        // The kernel's real_t is selected to match T1 (see buildKernels); only
        // the 32-bit float / 64-bit double cases are supported.
        static_assert(std::is_floating_point<T1>::value &&
                      (sizeof(T1) == 4 || sizeof(T1) == 8),
                      "OpenCLSweepSolver requires T1 to be float or double");
        ncx = nx;
        ncy = ny;
        ncz = nz;
        dx = grid_dx;
        dy = grid_dy;
        dz = grid_dz;
        num_nodes = (ncx + 1) * (ncy + 1) * (ncz + 1);
        
        setupSweepDirections();

        // Initialize OpenCL
        initializeOpenCL();
        
        // Load and build kernels
        buildKernels();
        
        // Allocate device memory
        allocateDeviceMemory();

        // Determine optimal work group size
        determineOptimalWorkGroupSize();

        // Precompute the per-(direction,level) plane-node lists used by the
        // Gauss-Seidel sweep dispatch.
        buildPlaneNodeLists();

        // Bind the arguments that never change (buffers + grid scalars) once,
        // so the per-plane executeSweep() only sets plane_offset/plane_count.
        setStaticKernelArgs();

        initialized = true;
    }
    
    void setSweepType(SweepType type) {
        current_sweep_type = type;
    }

    // =========================================================================
    // Profiling (opt-in)
    // =========================================================================

    /// Enable/disable GPU profiling.  Safe to toggle at any time; the queue is
    /// always created profiling-capable, so this just controls whether events
    /// are attached and timed.
    void setProfiling(bool on) { profiling = on; }
    bool isProfiling() const { return profiling; }

    /// Zero the accumulated profiling counters.
    void resetProfile() {
        prof_kernel_ns = prof_mirror_ns = prof_upload_ns = prof_download_ns = 0;
        prof_launches = 0;
        prof_wall_ms = 0.0;
    }

    /// Print the accumulated profiling breakdown.  "kernel busy" is the
    /// fraction of sweep wall-clock time the device actually spent executing
    /// sweep kernels -- a low value means the run is launch/host-bound rather
    /// than compute-bound.
    void reportProfile(std::ostream& os = std::cout) const {
        if (prof_launches == 0) {
            os << "  [OpenCL profile] no kernel launches recorded\n";
            return;
        }
        const double kernel_ms   = prof_kernel_ns   * 1e-6;
        const double mirror_ms   = prof_mirror_ns   * 1e-6;
        const double upload_ms   = prof_upload_ns   * 1e-6;
        const double download_ms = prof_download_ns * 1e-6;
        // Everything in the sweep wall that is neither kernel execution nor the
        // per-cycle mirror copy: host/driver enqueue cost + inter-launch idle +
        // the per-cycle clFinish stalls.
        const double overhead_ms = prof_wall_ms - kernel_ms - mirror_ms;
        auto frac = [&](double ms){ return (prof_wall_ms > 0.0) ? 100.0*ms/prof_wall_ms : 0.0; };
        os << "  [OpenCL profile]\n"
           << "      sweep wall   : " << prof_wall_ms  << " ms\n"
           << "      kernel exec  : " << kernel_ms     << " ms  ("
                                      << frac(kernel_ms)   << "% of sweep wall)\n"
           << "      mirror copy  : " << mirror_ms     << " ms  ("
                                      << frac(mirror_ms)   << "% of sweep wall)\n"
           << "      other/overhd : " << overhead_ms   << " ms  ("
                                      << frac(overhead_ms) << "% of sweep wall)\n"
           << "      upload       : " << upload_ms     << " ms\n"
           << "      download     : " << download_ms   << " ms\n"
           << "      launches     : " << prof_launches << "\n";
    }

    // =========================================================================
    // Main Execution Method
    // =========================================================================
    
    /**
     * Execute 8-directional sweeps on GPU
     * 
     * @param tt        Travel time array (modified in-place)
     * @param slowness  Slowness values at each node
     * @param frozen    Boolean array indicating frozen nodes
     * @param num_sweeps Number of complete sweep cycles (default: 1)
     */
    void runSweeps(std::vector<T1>& tt,
                   const std::vector<T1>& slowness,
                   const std::vector<bool>& frozen,
                   size_t num_sweeps = 1) {
        
        if (!initialized) {
            throw std::runtime_error("OpenCL solver not initialized");
        }
        
        if (tt.size() != num_nodes || 
            slowness.size() != num_nodes || 
            frozen.size() != num_nodes) {
            throw std::runtime_error("Array size mismatch");
        }
        
        // Convert frozen from bool to uchar for OpenCL
        std::vector<unsigned char> frozen_uchar(num_nodes);
        for (size_t i = 0; i < num_nodes; ++i) {
            frozen_uchar[i] = frozen[i] ? 1 : 0;
        }
        
        // Upload data to device
        uploadData(tt, slowness, frozen_uchar);
        
        // Perform sweep cycles
        for (size_t cycle = 0; cycle < num_sweeps; ++cycle) {
            performSweepCycle();
        }

        // Download results
        downloadData(tt);
        
        // Ensure all OpenCL operations are complete before returning
        clFinish(queue);
    }
    
    /**
     * Run sweeps without uploading/downloading data (for iterative convergence loops)
     * Data must already be uploaded with uploadTravelTimes() or runSweeps()
     * Results must be downloaded with downloadTravelTimes()
     * 
     * @param num_sweeps Number of complete sweep cycles
     */
    void runSweepsNoTransfer(size_t num_sweeps = 1) {
        if (!initialized) {
            throw std::runtime_error("OpenCL solver not initialized");
        }
        
        // Perform sweep cycles
        for (size_t cycle = 0; cycle < num_sweeps; ++cycle) {
            performSweepCycle();
        }
        // Ensure all OpenCL operations are complete before returning
        clFinish(queue);
    }
    
    /**
     * Upload travel time data to GPU (for iterative convergence loops)
     */
    void uploadTravelTimes(const std::vector<T1>& tt) {
        if (!initialized) {
            throw std::runtime_error("OpenCL solver not initialized");
        }
        if (tt.size() != num_nodes) {
            throw std::runtime_error("Array size mismatch");
        }
        
        // Upload to BOTH buffers.  performSweepCycle() now syncs d_tt_in from
        // d_tt_out at its start, so d_tt_out must also hold the uploaded data;
        // otherwise the sync would clobber the freshly uploaded travel times
        // with whatever stale contents d_tt_out happened to hold.
        cl_event evt = nullptr;
        cl_event* ep = profiling ? &evt : nullptr;
        cl_int err = clEnqueueWriteBuffer(queue, d_tt_in, CL_TRUE, 0,
                                         num_nodes * sizeof(T1), tt.data(),
                                         0, nullptr, ep);
        checkError(err, "Uploading travel time data to d_tt_in");
        addEventTime(evt, prof_upload_ns);

        err = clEnqueueWriteBuffer(queue, d_tt_out, CL_TRUE, 0,
                                   num_nodes * sizeof(T1), tt.data(),
                                   0, nullptr, ep);
        checkError(err, "Uploading travel time data to d_tt_out");
        addEventTime(evt, prof_upload_ns);
    }
    
    /**
     * Download travel time results from GPU (for iterative convergence loops)
     */
    void downloadTravelTimes(std::vector<T1>& tt) {
        if (!initialized) {
            throw std::runtime_error("OpenCL solver not initialized");
        }
        if (tt.size() != num_nodes) {
            throw std::runtime_error("Array size mismatch");
        }
        
        downloadData(tt);
    }
    
    /**
     * Get device information string
     */
    std::string getDeviceInfo() const {
        if (!device) return "No device initialized";
        
        char device_name[256];
        char vendor[256];
        cl_ulong global_mem;
        cl_ulong local_mem;
        size_t max_work_group;
        
        clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(device_name), device_name, nullptr);
        clGetDeviceInfo(device, CL_DEVICE_VENDOR, sizeof(vendor), vendor, nullptr);
        clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(global_mem), &global_mem, nullptr);
        clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(local_mem), &local_mem, nullptr);
        clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(max_work_group), &max_work_group, nullptr);
        
        std::ostringstream oss;
        oss << "    Device: " << device_name << "\n"
            << "    Vendor: " << vendor << "\n"
            << "    Global Memory: " << (global_mem / (1024*1024)) << " MB\n"
            << "    Local Memory: " << (local_mem / 1024) << " KB\n"
            << "    Max Work Group Size: " << max_work_group << "\n";
        
        return oss.str();
    }
    
private:
    // =========================================================================
    // OpenCL Setup
    // =========================================================================
    
    void initializeOpenCL() {
        cl_int err;
        
        // Get platform
        err = clGetPlatformIDs(1, &platform, nullptr);
        checkError(err, "Getting platform");
        
        // Get device (prefer GPU)
        err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, nullptr);
        if (err != CL_SUCCESS) {
            std::cout << "No GPU found, trying CPU...\n";
            err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 1, &device, nullptr);
        }
        checkError(err, "Getting device");
        
        // Create context
        context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &err);
        checkError(err, "Creating context");
        
        // Create command queue.  Profiling is enabled so OpenCL event
        // timestamps are available when setProfiling(true) is used; this has
        // negligible cost when profiling is off (no events are attached).
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
        // Use the kernel source embedded at build time.  As a development
        // convenience, an on-disk Grid3Drn_kernels.cl in a few standard
        // locations takes precedence so the kernels can be tweaked without
        // recompiling the extension.
        std::string kernel_paths[] = {
            "Grid3Drn_kernels.cl",           // Current directory
            "../Grid3Drn_kernels.cl",        // Parent directory
            "../../Grid3Drn_kernels.cl",     // Two levels up
            "./kernels/Grid3Drn_kernels.cl", // Kernels subdirectory
        };

        std::string kernel_source;
        bool found = false;
        for (const auto& path : kernel_paths) {
            try {
                kernel_source = loadKernelSource(path);
                found = true;
                break;
            } catch (...) {
                // Try next path
            }
        }

        if (!found) {
            // Fall back to the source compiled into the library.
            kernel_source = Grid3Drn_kernels_src;
        }

        const char* source_cstr = kernel_source.c_str();
        size_t source_size = kernel_source.length();
        
        cl_int err;
        program = clCreateProgramWithSource(context, 1, &source_cstr, &source_size, &err);
        checkError(err, "Creating program");

        // Select the kernel's real_t to match the host scalar type T1, so the
        // device interprets the buffers and scalar args exactly as the host
        // packs them.  (Previously the kernel was always built as float, which
        // silently corrupted the double instantiation.)
        // OLD CODE:
        //     err = clBuildProgram(program, 1, &device, "-cl-fast-relaxed-math",
        //                          nullptr, nullptr);
        std::string build_opts = "-cl-fast-relaxed-math";
        if (std::is_same<T1, double>::value) {
            // Double build: make real_t = double (the kernel enables cl_khr_fp64
            // under USE_DOUBLE).  Fail loudly if the device has no fp64 support
            // rather than silently producing garbage.
            cl_device_fp_config fp64 = 0;
            clGetDeviceInfo(device, CL_DEVICE_DOUBLE_FP_CONFIG,
                            sizeof(fp64), &fp64, nullptr);
            if (fp64 == 0) {
                throw std::runtime_error("Double precision requested (T1 = double) "
                                         "but the OpenCL device reports no fp64 "
                                         "(cl_khr_fp64) support");
            }
            build_opts += " -DUSE_DOUBLE=1";
        } else {
            // Float build: unsuffixed FP literals in the kernel are 'double' by
            // language rule; force them to single precision so no float<->double
            // promotion occurs (the kernel deliberately avoids 'f' suffixes so
            // the same source serves the double build).
            build_opts += " -cl-single-precision-constant";
        }

        // Build program
        err = clBuildProgram(program, 1, &device, build_opts.c_str(), nullptr, nullptr);
        if (err != CL_SUCCESS) {
            // Get build log
            size_t log_size;
            clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &log_size);
            std::vector<char> build_log(log_size);
            clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, log_size, build_log.data(), nullptr);
            
            std::cerr << "Build failed:\n" << build_log.data() << "\n";
            throw std::runtime_error("Kernel build failed");
        }
        
        // Create kernel objects
        kernel_basic = clCreateKernel(program, "sweep_update_basic", &err);
        checkError(err, "Creating basic kernel");
        
        kernel_weno3 = clCreateKernel(program, "sweep_update_weno3", &err);
        checkError(err, "Creating weno3 kernel");
    }
    
    void allocateDeviceMemory() {
        cl_int err;
        size_t mem_size = num_nodes * sizeof(T1);
        size_t frozen_size = num_nodes * sizeof(unsigned char);
        
        // Allocate device buffers
        d_tt_in = clCreateBuffer(context, CL_MEM_READ_WRITE, mem_size, nullptr, &err);
        checkError(err, "Allocating tt_in buffer");
        
        d_tt_out = clCreateBuffer(context, CL_MEM_READ_WRITE, mem_size, nullptr, &err);
        checkError(err, "Allocating tt_out buffer");
        
        d_slowness = clCreateBuffer(context, CL_MEM_READ_ONLY, mem_size, nullptr, &err);
        checkError(err, "Allocating slowness buffer");
        
        d_frozen = clCreateBuffer(context, CL_MEM_READ_ONLY, frozen_size, nullptr, &err);
        checkError(err, "Allocating frozen buffer");
        
        if ( verbose )
            std::cout << "\n  Allocated " << ((2*mem_size + mem_size + frozen_size)/(1024*1024))
            << " MB device memory\n";
    }
    
    /**
     * Determine optimal work group size based on device capabilities
     * 
     * Strategy:
     * 1. Query device max work group size and preferred multiple
     * 2. Consider kernel-specific requirements
     * 3. Balance occupancy vs register pressure
     * 4. Ensure divisibility for 3D grid
     */
    void determineOptimalWorkGroupSize() {
        cl_int err;
        
        // Get device capabilities
        size_t max_work_group_size;
        err = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, 
                             sizeof(max_work_group_size), &max_work_group_size, nullptr);
        checkError(err, "Querying max work group size");
        
        size_t max_work_item_sizes[3];
        err = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_SIZES,
                             sizeof(max_work_item_sizes), max_work_item_sizes, nullptr);
        checkError(err, "Querying max work item sizes");
        
        // Get preferred work group size multiple (warp/wavefront size)
        size_t preferred_multiple = 32; // Default for NVIDIA
        
        // Try to query actual preferred size (OpenCL 1.1+)
        #ifdef CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE
        size_t kernel_preferred;
        err = clGetKernelWorkGroupInfo(kernel_weno3, device,
                                      CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
                                      sizeof(kernel_preferred), &kernel_preferred, nullptr);
        if (err == CL_SUCCESS) {
            preferred_multiple = kernel_preferred;
        }
        #endif
        
        // Get vendor to apply vendor-specific optimizations
        char vendor[256];
        clGetDeviceInfo(device, CL_DEVICE_VENDOR, sizeof(vendor), vendor, nullptr);
        std::string vendor_str(vendor);
        
        // =====================================================================
        // STRATEGY: Choose work group size based on vendor and capabilities
        // =====================================================================
        
        size_t target_threads;
        
        if (vendor_str.find("NVIDIA") != std::string::npos) {
            // NVIDIA GPUs: Maximize occupancy
            // - Warp size = 32
            // - Target 256-1024 threads per block
            // - Prefer 512 or 1024 for compute-heavy kernels
            target_threads = (max_work_group_size >= 1024) ? 1024 : 512;
            
        } else if (vendor_str.find("AMD") != std::string::npos || 
                   vendor_str.find("Advanced Micro Devices") != std::string::npos) {
            // AMD GPUs: Balance wavefront utilization
            // - Wavefront size = 64 (RDNA) or 32 (RDNA3)
            // - Prefer multiples of 64 for older, 32 for newer
            target_threads = (max_work_group_size >= 512) ? 512 : 256;
            
        } else if (vendor_str.find("Intel") != std::string::npos) {
            // Intel GPUs: Conservative approach
            // - SIMD width varies (8, 16, 32)
            // - Prefer 256 threads for good occupancy without register pressure
            target_threads = 256;
            
        } else if (vendor_str.find("Apple") != std::string::npos) {
            // Apple Silicon: Conservative
            // - Threadgroup size typically 256-512
            target_threads = 256;
            
        } else {
            // Unknown vendor: Safe default
            target_threads = std::min(max_work_group_size, size_t(256));
        }
        
        // Ensure we don't exceed device limits
        target_threads = std::min(target_threads, max_work_group_size);
        
        // =====================================================================
        // DECOMPOSE into 3D work group dimensions
        // =====================================================================
        
        // Strategy: Favor X and Y dimensions (better memory coalescing)
        // Common good configurations:
        //   1024 threads: (16, 16, 4) or (32, 8, 4)
        //    512 threads: (16, 8, 4) or (8, 8, 8)
        //    256 threads: (8, 8, 4) or (16, 4, 4)
        
        if (target_threads >= 1024) {
            // Large work groups: 16×16×4 = 1024
            optimal_local_size[0] = 16;
            optimal_local_size[1] = 16;
            optimal_local_size[2] = 4;
        } else if (target_threads >= 512) {
            // Medium-large: 16×8×4 = 512
            optimal_local_size[0] = 16;
            optimal_local_size[1] = 8;
            optimal_local_size[2] = 4;
        } else if (target_threads >= 256) {
            // Medium: 8×8×4 = 256
            optimal_local_size[0] = 8;
            optimal_local_size[1] = 8;
            optimal_local_size[2] = 4;
        } else {
            // Small: 8×4×4 = 128 (very conservative)
            optimal_local_size[0] = 8;
            optimal_local_size[1] = 4;
            optimal_local_size[2] = 4;
        }
        
        // Clamp to device maximum dimensions
        optimal_local_size[0] = std::min(optimal_local_size[0], max_work_item_sizes[0]);
        optimal_local_size[1] = std::min(optimal_local_size[1], max_work_item_sizes[1]);
        optimal_local_size[2] = std::min(optimal_local_size[2], max_work_item_sizes[2]);
        
        // Also clamp to the kernel-specific work group size limit, which can be
        // smaller than the device maximum when register pressure is high.
        for (cl_kernel k : { kernel_basic, kernel_weno3 }) {
            size_t kernel_wgs = 0;
            if (clGetKernelWorkGroupInfo(k, device, CL_KERNEL_WORK_GROUP_SIZE,
                                         sizeof(kernel_wgs), &kernel_wgs, nullptr) == CL_SUCCESS
                    && kernel_wgs > 0) {
                max_work_group_size = std::min(max_work_group_size, kernel_wgs);
            }
        }

        // Verify total size
        size_t total_size = optimal_local_size[0] * optimal_local_size[1] * optimal_local_size[2];
        if (total_size > max_work_group_size) {
            // Fallback: reduce Z dimension
            optimal_local_size[2] = max_work_group_size / (optimal_local_size[0] * optimal_local_size[1]);
            optimal_local_size[2] = std::max(size_t(1), optimal_local_size[2]);
        }
        
        if ( verbose ) {
            std::cout << "  Work group size: (" 
                      << optimal_local_size[0] << ", "
                      << optimal_local_size[1] << ", "
                      << optimal_local_size[2] << ") = "
                      << (optimal_local_size[0] * optimal_local_size[1] * optimal_local_size[2])
                      << " threads\n";
            std::cout << "  (Max possible: " << max_work_group_size 
                      << ", Vendor: " << vendor_str << ")\n";
        }
    }
    
    // =========================================================================
    // Data Transfer
    // =========================================================================
    
    void uploadData(const std::vector<T1>& tt,
                   const std::vector<T1>& slowness,
                   const std::vector<unsigned char>& frozen) {
        cl_int err;
        cl_event evt = nullptr;
        cl_event* ep = profiling ? &evt : nullptr;

        // Upload travel times to BOTH buffers (in case we swap)
        err = clEnqueueWriteBuffer(queue, d_tt_in, CL_TRUE, 0,
                                  num_nodes * sizeof(T1), tt.data(), 0, nullptr, ep);
        checkError(err, "Uploading tt data to d_tt_in");
        addEventTime(evt, prof_upload_ns);

        // Also initialize d_tt_out with the same data to avoid garbage values
        err = clEnqueueWriteBuffer(queue, d_tt_out, CL_TRUE, 0,
                                  num_nodes * sizeof(T1), tt.data(), 0, nullptr, ep);
        checkError(err, "Uploading tt data to d_tt_out");
        addEventTime(evt, prof_upload_ns);

        err = clEnqueueWriteBuffer(queue, d_slowness, CL_TRUE, 0,
                                  num_nodes * sizeof(T1), slowness.data(), 0, nullptr, ep);
        checkError(err, "Uploading slowness data");
        addEventTime(evt, prof_upload_ns);

        err = clEnqueueWriteBuffer(queue, d_frozen, CL_TRUE, 0,
                                  num_nodes * sizeof(unsigned char), frozen.data(), 0, nullptr, ep);
        checkError(err, "Uploading frozen data");
        addEventTime(evt, prof_upload_ns);
    }

    void downloadData(std::vector<T1>& tt) {
        cl_event evt = nullptr;
        cl_event* ep = profiling ? &evt : nullptr;
        cl_int err = clEnqueueReadBuffer(queue, d_tt_out, CL_TRUE, 0,
                                        num_nodes * sizeof(T1), tt.data(), 0, nullptr, ep);
        checkError(err, "Downloading results");
        addEventTime(evt, prof_download_ns);
    }
    
    // =========================================================================
    // Sweep Execution
    // =========================================================================
    
    void setupSweepDirections() {
        // Define all 8 sweep directions
        sweep_directions = {
            {+1, +1, +1,  0,     0,     0    },  // Direction 1
            {-1, +1, +1,  ncx,   0,     0    },  // Direction 2
            {+1, -1, +1,  0,     ncy,   0    },  // Direction 3
            {-1, -1, +1,  ncx,   ncy,   0    },  // Direction 4
            {+1, +1, -1,  0,     0,     ncz  },  // Direction 5
            {-1, +1, -1,  ncx,   0,     ncz  },  // Direction 6
            {+1, -1, -1,  0,     ncy,   ncz  },  // Direction 7
            {-1, -1, -1,  ncx,   ncy,   ncz  }   // Direction 8
        };
    }

    /**
     * Precompute, for each of the 8 sweep directions, the grid-node linear
     * indices grouped by diagonal level = ip + jp + kp (ip/jp/kp = node indices
     * oriented along the sweep direction).  The result is one permutation of all
     * num_nodes nodes per direction, ordered by level, stored back-to-back in a
     * single device buffer d_plane_nodes (direction d occupies
     * [d*num_nodes, (d+1)*num_nodes)).  level_offsets[d][L] gives the start of
     * level L within direction d's block; level_offsets[d][maxlevel+1] ==
     * num_nodes.  Built once and reused for every sweep cycle.
     *
     * Memory: 8 * num_nodes * sizeof(uint).  For very large grids this can be
     * reduced by exploiting the index-reflection symmetry between directions
     * (the -dir permutation is the +dir one with the axis index mirrored), but
     * it is left explicit here for clarity.
     */
    void buildPlaneNodeLists() {
        maxlevel = static_cast<int>(ncx + ncy + ncz);
        const int nlev = maxlevel + 1;             // levels 0 .. maxlevel
        const size_t nx1 = ncx + 1, ny1 = ncy + 1, nz1 = ncz + 1;

        level_offsets.assign(8, std::vector<cl_uint>(nlev + 1, 0));
        std::vector<cl_uint> all_perm(8 * num_nodes);

        for (size_t d = 0; d < 8; ++d) {
            const auto& dir = sweep_directions[d];
            std::vector<cl_uint>& off = level_offsets[d];

            // Counting sort of the nodes by level.
            std::vector<cl_uint> count(nlev, 0);
            for (size_t k = 0; k < nz1; ++k)
                for (size_t j = 0; j < ny1; ++j)
                    for (size_t i = 0; i < nx1; ++i) {
                        const int ip = dir.i_dir > 0 ? (int)i : (int)(ncx - i);
                        const int jp = dir.j_dir > 0 ? (int)j : (int)(ncy - j);
                        const int kp = dir.k_dir > 0 ? (int)k : (int)(ncz - k);
                        ++count[ip + jp + kp];
                    }

            // Prefix sum -> per-level start offsets (off[nlev] == num_nodes).
            cl_uint running = 0;
            for (int L = 0; L < nlev; ++L) { off[L] = running; running += count[L]; }
            off[nlev] = running;

            // Scatter the node indices into their level buckets.
            std::vector<cl_uint> cursor(off.begin(), off.begin() + nlev);
            for (size_t k = 0; k < nz1; ++k)
                for (size_t j = 0; j < ny1; ++j)
                    for (size_t i = 0; i < nx1; ++i) {
                        const int ip = dir.i_dir > 0 ? (int)i : (int)(ncx - i);
                        const int jp = dir.j_dir > 0 ? (int)j : (int)(ncy - j);
                        const int kp = dir.k_dir > 0 ? (int)k : (int)(ncz - k);
                        const size_t idx = (k * ny1 + j) * nx1 + i;
                        all_perm[d * num_nodes + cursor[ip + jp + kp]++] =
                            static_cast<cl_uint>(idx);
                    }
        }

        cl_int err;
        d_plane_nodes = clCreateBuffer(context, CL_MEM_READ_ONLY,
                                       8 * num_nodes * sizeof(cl_uint),
                                       nullptr, &err);
        checkError(err, "Allocating plane-node buffer");
        err = clEnqueueWriteBuffer(queue, d_plane_nodes, CL_TRUE, 0,
                                   8 * num_nodes * sizeof(cl_uint),
                                   all_perm.data(), 0, nullptr, nullptr);
        checkError(err, "Uploading plane-node lists");
    }

    // Set the kernel arguments that never change after initialize() once, for
    // both kernels.  Only plane_offset/plane_count (args 11 and 12) vary per
    // plane launch and are set in executeSweep().  On Apple's OpenCL->Metal
    // layer the per-enqueue host cost dominates the run (each tiny plane launch
    // costs ~27us regardless of its node count), so hoisting these 11 args out
    // of the inner loop removes 11*launches clSetKernelArg calls per cycle.
    //
    // The buffer handles (d_tt_in, d_slowness, d_frozen, d_plane_nodes) and the
    // grid scalars (dx/dy/dz, ncx/ncy/ncz) are fixed for the solver's lifetime;
    // uploadTravelTimes() only rewrites the *contents* of d_tt_in/d_tt_out, not
    // their handles, so these bindings stay valid across sweep cycles.
    void setStaticKernelArgs() {
        for (cl_kernel kernel : { kernel_basic, kernel_weno3 }) {
            cl_int err = 0;
            cl_uint a = 0;
            // In-place Gauss-Seidel: bind d_tt_in to BOTH the input and output
            // buffer arguments so a node sees its already-updated neighbours.
            err |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &d_tt_in);
            err |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &d_tt_in);
            err |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &d_slowness);
            err |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &d_frozen);
            err |= clSetKernelArg(kernel, a++, sizeof(T1), &dx);
            err |= clSetKernelArg(kernel, a++, sizeof(T1), &dy);
            err |= clSetKernelArg(kernel, a++, sizeof(T1), &dz);
            err |= clSetKernelArg(kernel, a++, sizeof(cl_uint), &ncx);
            err |= clSetKernelArg(kernel, a++, sizeof(cl_uint), &ncy);
            err |= clSetKernelArg(kernel, a++, sizeof(cl_uint), &ncz);
            err |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &d_plane_nodes);
            checkError(err, "Setting static kernel arguments");
        }
    }

    void performSweepCycle() {
        // Select kernel based on sweep type
        cl_kernel current_kernel;
        switch (current_sweep_type) {
            case SweepType::BASIC:
                current_kernel = kernel_basic;
                break;
            case SweepType::WENO3:
                current_kernel = kernel_weno3;
                break;
            default:
                current_kernel = kernel_basic;
        }

        // =====================================================================
        // OLD CODE (pure-Jacobi, full-grid passes with buffer swaps):
        //
        // Each of the 8 "sweeps" was an identical full-grid update reading
        // d_tt_in and writing d_tt_out, with the buffers swapped between
        // sweeps.  The sweep-direction arguments were ignored, so there was
        // no sweep ordering at all.  That breaks the causal upwind ordering
        // WENO3 needs (it converges only as a Gauss-Seidel sweep), which is
        // why the Jacobi version needed the monotone clamp band-aid and still
        // lost accuracy.  See the start-of-cycle sync that was needed to work
        // around the odd-swap stale-buffer issue.
        //
        // {
        //     cl_int err = clEnqueueCopyBuffer(queue, d_tt_out, d_tt_in,
        //                                      0, 0, num_nodes * sizeof(T1),
        //                                      0, nullptr, nullptr);
        //     checkError(err, "Syncing d_tt_in from d_tt_out at start of sweep cycle");
        //     clFinish(queue);
        // }
        // for (size_t dir = 0; dir < 8; ++dir) {
        //     executeSweep(current_kernel, dir);
        //     if (dir < 7) {
        //         std::swap(d_tt_in, d_tt_out);
        //     }
        // }
        // =====================================================================

        // NEW CODE: Gauss-Seidel plane-sweep ordering.
        //
        // Each of the 8 sweep directions is processed plane by plane in order
        // of increasing diagonal level = ip+jp+kp (ip/jp/kp oriented along the
        // sweep direction, see the kernel).  All nodes on a plane are mutually
        // independent and updated in parallel; planes are launched
        // sequentially so that a node's upwind neighbours (lower level) are
        // already updated within the same sweep.  Updates are done IN PLACE in
        // d_tt_in (d_tt_in is bound to both kernel buffer arguments), which is
        // what makes this a genuine Gauss-Seidel sweep rather than Jacobi.
        //
        // Each plane is now dispatched as a 1-D range over exactly its nodes,
        // using the per-(direction,level) lists precomputed in
        // buildPlaneNodeLists().  (Previously every plane launch dispatched the
        // full grid box and off-plane work-items returned at a gate -- O(n)
        // launches each spawning O(n^3) threads for O(n^2) useful updates.)
        std::chrono::high_resolution_clock::time_point t_begin;
        if (profiling) t_begin = std::chrono::high_resolution_clock::now();

        for (size_t dir = 0; dir < 8; ++dir) {
            for (int level = 0; level <= maxlevel; ++level) {
                executeSweep(current_kernel, dir, level);
            }
        }

        // Mirror the in-place result into d_tt_out so downloadData() (which
        // reads d_tt_out) and the "d_tt_out holds the latest result" invariant
        // continue to hold for the surrounding upload/download code.
        cl_event mirror_evt = nullptr;
        cl_event* mirror_ep = profiling ? &mirror_evt : nullptr;
        cl_int err = clEnqueueCopyBuffer(queue, d_tt_in, d_tt_out,
                                         0, 0, num_nodes * sizeof(T1),
                                         0, nullptr, mirror_ep);
        checkError(err, "Mirroring Gauss-Seidel result to d_tt_out");
        clFinish(queue);

        if (profiling) {
            auto t_end = std::chrono::high_resolution_clock::now();
            prof_wall_ms +=
                std::chrono::duration<double, std::milli>(t_end - t_begin).count();
            // All work of this cycle is complete (clFinish above), so the event
            // timestamps can be read without further synchronization.
            for (cl_event e : prof_events) {
                prof_kernel_ns += eventNs(e);
                clReleaseEvent(e);
            }
            prof_events.clear();
            if (mirror_evt) {
                prof_mirror_ns += eventNs(mirror_evt);
                clReleaseEvent(mirror_evt);
            }
        }
    }

    void executeSweep(cl_kernel kernel, size_t direction_idx, int level) {
        // Node range for this (direction, level) plane within the concatenated
        // d_plane_nodes buffer.  Skip the launch entirely if the plane is empty.
        const cl_uint level_start = level_offsets[direction_idx][level];
        const cl_uint level_end   = level_offsets[direction_idx][level + 1];
        const cl_uint plane_count = level_end - level_start;
        if (plane_count == 0) return;
        const cl_uint plane_offset =
            static_cast<cl_uint>(direction_idx * num_nodes) + level_start;

        // Only the two per-plane arguments are set here; the other 11 (buffers,
        // grid scalars, d_plane_nodes) are bound once in setStaticKernelArgs()
        // because they never change.  Their argument indices are 11 (offset)
        // and 12 (count) to match the kernel signature.
        //
        // OLD CODE (set all 13 args on every plane launch):
        // cl_uint arg_idx = 0;
        // err |= clSetKernelArg(kernel, arg_idx++, sizeof(cl_mem), &d_tt_in);   // in
        // err |= clSetKernelArg(kernel, arg_idx++, sizeof(cl_mem), &d_tt_in);   // out (in-place GS)
        // err |= clSetKernelArg(kernel, arg_idx++, sizeof(cl_mem), &d_slowness);
        // err |= clSetKernelArg(kernel, arg_idx++, sizeof(cl_mem), &d_frozen);
        // err |= clSetKernelArg(kernel, arg_idx++, sizeof(T1), &dx);
        // err |= clSetKernelArg(kernel, arg_idx++, sizeof(T1), &dy);
        // err |= clSetKernelArg(kernel, arg_idx++, sizeof(T1), &dz);
        // err |= clSetKernelArg(kernel, arg_idx++, sizeof(cl_uint), &ncx);
        // err |= clSetKernelArg(kernel, arg_idx++, sizeof(cl_uint), &ncy);
        // err |= clSetKernelArg(kernel, arg_idx++, sizeof(cl_uint), &ncz);
        // err |= clSetKernelArg(kernel, arg_idx++, sizeof(cl_mem), &d_plane_nodes);
        const cl_uint OFFSET_ARG = 11;
        const cl_uint COUNT_ARG  = 12;
        cl_int err = 0;
        err |= clSetKernelArg(kernel, OFFSET_ARG, sizeof(cl_uint), &plane_offset);
        err |= clSetKernelArg(kernel, COUNT_ARG,  sizeof(cl_uint), &plane_count);
        checkError(err, "Setting per-plane kernel arguments");

        // 1-D launch over exactly the plane's nodes.  Collapse the 3-D optimal
        // work-group size to a single dimension and round the global size up to
        // a multiple of it (the kernel guards the remainder with plane_count).
        size_t local_1d = optimal_local_size[0] * optimal_local_size[1] *
                          optimal_local_size[2];
        if (local_1d < 1) local_1d = 1;
        size_t global_1d = ((plane_count + local_1d - 1) / local_1d) * local_1d;

        // OLD CODE (3-D launch over the full grid box):
        // size_t local_work_size[3] = {
        //     optimal_local_size[0], optimal_local_size[1], optimal_local_size[2]
        // };
        // size_t global_work_size[3];
        // global_work_size[0] = ((ncx + 1 + local_work_size[0] - 1) / local_work_size[0]) * local_work_size[0];
        // global_work_size[1] = ((ncy + 1 + local_work_size[1] - 1) / local_work_size[1]) * local_work_size[1];
        // global_work_size[2] = ((ncz + 1 + local_work_size[2] - 1) / local_work_size[2]) * local_work_size[2];
        // err = clEnqueueNDRangeKernel(queue, kernel, 3, nullptr,
        //                             global_work_size, local_work_size,
        //                             0, nullptr, nullptr);

        // Launch kernel.  When profiling, attach an event and stash it; the
        // events are queried in performSweepCycle() after the cycle's clFinish
        // (so reading their timestamps adds no extra synchronization).
        cl_event evt = nullptr;
        cl_event* ep = profiling ? &evt : nullptr;
        err = clEnqueueNDRangeKernel(queue, kernel, 1, nullptr,
                                    &global_1d, &local_1d,
                                    0, nullptr, ep);
        checkError(err, "Launching kernel");
        if (profiling) {
            prof_events.push_back(evt);
            ++prof_launches;
        }

        // OLD CODE: drain the queue after every single plane launch.
        //
        //     clFinish(queue);
        //
        // This was a full host<->device round-trip per launch (8 dirs x
        // (ncx+ncy+ncz+1) levels per sweep cycle), and it dominated the run
        // time while the plane kernels themselves touch only a handful of
        // nodes.  It is NOT needed for correctness: the command queue is
        // in-order, so launch N+1 is guaranteed not to start until launch N
        // has completed, which is exactly the causal Gauss-Seidel ordering the
        // plane-sweep relies on.  A single clFinish() at the end of
        // performSweepCycle() (and the blocking read in downloadData()) is
        // enough to ensure all work is done before the host reads results.
        // clFinish(queue);
    }
    
    // =========================================================================
    // Utilities
    // =========================================================================
    
    std::string loadKernelSource(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open kernel file: " + filename);
        }
        
        std::ostringstream oss;
        oss << file.rdbuf();
        return oss.str();
    }
    
    void checkError(cl_int err, const std::string& operation) {
        if (err != CL_SUCCESS) {
            std::ostringstream oss;
            oss << "OpenCL error during " << operation << ": " << err;
            throw std::runtime_error(oss.str());
        }
    }

    // Device execution time of a completed event, in nanoseconds.
    static cl_ulong eventNs(cl_event evt) {
        cl_ulong s = 0, e = 0;
        clGetEventProfilingInfo(evt, CL_PROFILING_COMMAND_START, sizeof(s), &s, nullptr);
        clGetEventProfilingInfo(evt, CL_PROFILING_COMMAND_END,   sizeof(e), &e, nullptr);
        return (e > s) ? (e - s) : 0;
    }

    // Add a just-completed event's device time to an accumulator and release
    // it.  No-op when profiling is off (evt is then null).  Only valid for
    // events whose command has finished (e.g. blocking transfers).
    void addEventTime(cl_event evt, cl_ulong& accumulator) {
        if (!profiling || !evt) return;
        accumulator += eventNs(evt);
        clReleaseEvent(evt);
    }
    
    void cleanup() {
        // Emit the profiling summary (once, at teardown) for a profiling run.
        if (profiling && prof_launches > 0) {
            reportProfile(std::cout);
        }
        for (cl_event e : prof_events) clReleaseEvent(e);
        prof_events.clear();

        if (d_tt_in) clReleaseMemObject(d_tt_in);
        if (d_tt_out) clReleaseMemObject(d_tt_out);
        if (d_plane_nodes) clReleaseMemObject(d_plane_nodes);
        if (d_slowness) clReleaseMemObject(d_slowness);
        if (d_frozen) clReleaseMemObject(d_frozen);
        
        if (kernel_basic) clReleaseKernel(kernel_basic);
        if (kernel_weno3) clReleaseKernel(kernel_weno3);
        
        if (program) clReleaseProgram(program);
        if (queue) clReleaseCommandQueue(queue);
        if (context) clReleaseContext(context);
    }
};

} // namespace ttcr

#endif // GRID3DRN_OPENCL_H

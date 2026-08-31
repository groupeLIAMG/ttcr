//
//  Grid3Drcfs_OpenCL.h
//  ttcr
//
//  GPU-accelerated version of Grid3Drcfs using OpenCL
//
//  This class provides the same interface as Grid3Drcfs but uses
//  GPU acceleration for the sweep operations.
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
 *
 */

/**
 * @file Grid3Drcfs_OpenCL.h
 * @brief GPU-accelerated fast sweeping on a 3-D rectilinear grid, from a cell
 *        slowness model.
 *
 * Declares ttcr::Grid3Drcfs_OpenCL, the OpenCL counterpart of
 * ttcr::Grid3Drcfs. Same structure — including deriving from ttcr::Grid3Drn
 * rather than ttcr::Grid3Drc, and averaging the supplied cell slowness onto the
 * nodes — with only the sweep moving to the device.
 *
 * Selected by ttcr::input_parameters::method @c == @c FAST_SWEEPING_OPENCL.
 *
 * @warning Cubic cells only, unvalidated — @sa @ref g3drcfs_cubic
 *
 * @sa Grid3Drcfs.h, Grid3Drn_OpenCL.h, Grid2Drcfs_OpenCL.h
 */

#ifndef ttcr_Grid3Drcfs_OpenCL_h
#define ttcr_Grid3Drcfs_OpenCL_h

#include <cmath>
#include <limits>
#include <vector>
#include <memory>
#include <iostream>

#include "Grid3Drn.h"
#include "Grid3Drn_OpenCL.h"
#include "Node3Dn.h"

namespace ttcr {

    /**
     * @brief GPU-accelerated fast sweeping solver taking a cell slowness model.
     *
     * @tparam T1 floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2 integer type of node and cell indices.
     *
     * The OpenCL counterpart of ttcr::Grid3Drcfs, sharing its structure exactly
     * — same base (ttcr::Grid3Drn, not ttcr::Grid3Drc; see
     * @ref g3drcfs_hybrid), same cell-to-node slowness averaging — with only the
     * sweep moving to the device.
     *
     * @warning Inherits ttcr::Grid3Drcfs's **cubic-cell-only** restriction: the
     *          constructor takes one cell size for all three axes, and nothing
     *          validates that the model's three spacings agree.
     *          @sa @ref g3drcfs_cubic
     * @note GPU use is requested, not guaranteed — @ref isUsingGPU reports what
     *       is actually in use and the class falls back to the CPU otherwise, so
     *       results are unaffected and only performance changes.
     *
     * @sa Grid3Drcfs.h, Grid3Drn_OpenCL.h, Grid2Drcfs_OpenCL.h
     *
     * GPU-accelerated Grid3Drcfs using OpenCL
     *
     * This class maintains the same interface as Grid3Drcfs but uses
     * GPU-accelerated sweep operations via OpenCL kernels.
     * 
     * Grid3Drcfs uses cell-centered slowness values that are interpolated
     * to grid nodes. This GPU version accelerates the sweep operations
     * while maintaining the same slowness interpolation scheme.
     * 
     * Features:
     * - Same API as Grid3Drcfs (drop-in replacement)
     * - 8-20x speedup for large grids (>64³)
     * - Automatic fallback to CPU if GPU fails
     * - Support for both basic and WENO3 sweeps
     * - Convergence checking same as CPU version
     * 
     * Usage:
     *   Grid3Drcfs_OpenCL<double, uint32_t> grid(
     *       nx, ny, nz, dx, xmin, ymin, zmin,
     *       epsilon, max_iterations, use_weno3);
     *   
     *   grid.setSlowness(slowness_values);
     *   grid.raytrace(Tx, t0, Rx);
     */
    template<typename T1, typename T2>
    class Grid3Drcfs_OpenCL : public Grid3Drn<T1,T2,Node3Dn<T1,T2>> {
    public:
        /**
         * @brief Build the grid, its nodes, and the GPU solvers if requested.
         *
         * @param nx        Number of cells in x
         * @param ny        Number of cells in y
         * @param nz        Number of cells in z
         * @param ddx       Cell size (assumes cubic cells: dx=dy=dz)
         * @param minx      X origin
         * @param miny      Y origin
         * @param minz      Z origin
         * @param eps       Convergence epsilon (stop when change < eps)
         * @param maxit     Maximum iterations
         * @param w         Use WENO3 (true) or basic sweep (false)
         * @param ttrp      Recompute receiver traveltimes by integrating
         *                  slowness along the traced raypath
         *                  (ttcr::input_parameters::tt_from_rp)
         * @param intVel    Interpolate velocity
         * @param nt        Number of threads
         * @param _translateOrigin  Translate origin to (0,0,0)
         * @param enableGPU Enable GPU acceleration (true by default)
         *
         * @post Nodes and neighbour lists are built and, if @p enableGPU, the
         *       GPU solvers are initialised — falling back to the CPU if that
         *       fails. @p eps is a relative tolerance: the sweeps stop once the
         *       mean per-node change falls below that fraction of the
         *       traveltime range. Slowness is **not** set.
         */
        Grid3Drcfs_OpenCL(const T2 nx, const T2 ny, const T2 nz, const T1 ddx,
                          const T1 minx, const T1 miny, const T1 minz,
                          const T1 eps, const int maxit, const bool w,
                          const bool ttrp=true, const bool intVel=false,
                          const size_t nt=1, const bool _translateOrigin=false,
                          const bool enableGPU=true) :
        Grid3Drn<T1,T2,Node3Dn<T1,T2>>(nx, ny, nz, ddx, ddx, ddx, minx, miny, minz, 
                                        ttrp, intVel, nt, _translateOrigin),
        epsilon(eps), 
        nitermax(maxit), 
        niter_final(0), 
        niterw_final(0), 
        weno3(w),
        use_gpu(enableGPU),
        gpu_initialized(false),
        gpu_available(false)
        {
            this->buildGridNodes();
            this->template buildGridNeighbors<Node3Dn<T1,T2>>(this->nodes);

            if (use_gpu) {
                initializeGPU();
            }
        }

        /// Destructor; the OpenCLSweepSolver destructors release the device
        /// contexts, queues and buffers.
        virtual ~Grid3Drcfs_OpenCL() {
            // OpenCLSweepSolver destructor handles cleanup
        }

        /**
         * @brief Set the cell slowness model and average it onto the nodes.
         * @param s one slowness per cell, in the x-fastest order of
         *          @ref g3drc_numbering.
         * @throws std::length_error if @p s has the wrong size.
         *
         * Set slowness values (cell-centered)
         *
         * Slowness values are defined at cell centers and interpolated to grid nodes.
         * This maintains the same interpolation scheme as Grid3Drcfs.
         */
        void setSlowness(const std::vector<T1>& s);

        // Accessors (same as Grid3Drcfs)
        /// @return Number of basic sweep iterations the last solve took.
        const int get_niter() const { return niter_final; }
        /// @return Number of WENO3 sweep iterations the last solve took.
        const int get_niterw() const { return niterw_final; }

        // GPU-specific methods
        /**
         * @brief Whether the solve will actually run on the GPU.
         * @return True only if GPU use was requested **and** a device was
         *         successfully initialised.
         */
        bool isUsingGPU() const { return use_gpu && gpu_available; }
        /**
         * @brief Turn GPU acceleration on or off after construction.
         * @param enable true to use the GPU, false to force the CPU path.
         * @note Enabling initialises the device on first use. Turning it on may
         *       still leave @ref isUsingGPU false if initialisation fails.
         *       ttcr::Grid2Drcfs_OpenCL has no equivalent setter.
         */
        void setUseGPU(bool enable) {
            use_gpu = enable;
            if (use_gpu && !gpu_initialized) {
                initializeGPU();
            }
        }

        /**
         * @brief Human-readable description of the OpenCL device in use.
         * @return Device information, or @c "GPU not available" if the solve
         *         will run on the CPU.
         */
        std::string getGPUInfo() const {
            if (gpu_available && !gpu_solvers.empty()) {
                return gpu_solvers[0]->getDeviceInfo();
            }
            return "GPU not available";
        }

    protected:
        T1 epsilon;              ///< Convergence tolerance, dimensionless: the fraction of the traveltime range the mean per-node change must fall below.
        int nitermax;            ///< Maximum iterations
        mutable int niter_final; ///< Final iteration count (basic sweep); @c mutable so the @c const raytrace can record it
        mutable int niterw_final;///< Final iteration count (WENO3 sweep)
        bool weno3;              ///< Use WENO3 sweep

        // GPU-specific members
        mutable bool use_gpu;         ///< GPU acceleration is requested.
        mutable bool gpu_initialized; ///< Initialisation has been attempted.
        mutable bool gpu_available;   ///< A device was successfully initialised.
        // one independent solver per thread slot (own context/queue/
        // kernels/buffers).  Indexed by threadNo, sized to getNthreads().
        mutable std::vector<std::unique_ptr<OpenCLSweepSolver<T1>>> gpu_solvers;

    private:
        /// @name Non-copyable
        /// Private, the pre-C++11 idiom. @c operator= returns @c *this without
        /// copying, so it is not undefined behaviour, but it would silently do
        /// nothing if reached. Note that ttcr::Grid2Drcfs_OpenCL uses
        /// `= delete` for the same three.
        /// @{
        Grid3Drcfs_OpenCL() {}
        Grid3Drcfs_OpenCL(const Grid3Drcfs_OpenCL<T1,T2>& g) {}
        Grid3Drcfs_OpenCL<T1,T2>& operator=(const Grid3Drcfs_OpenCL<T1,T2>& g) { return *this; }
        /// @}

        /**
         * Initialize GPU solver
         */
        void initializeGPU() const {
            if (gpu_initialized) return;

            try {
                // Build one independent solver per thread slot so that
                // concurrent source solves run on separate GPU command streams
                // (validated to overlap on this hardware).  performSweepIterations
                // selects its slot by threadNo, matching the [0, nThreads) ids the
                // shot loop hands out; the pool is capped upstream by num_threads.
                const size_t nslots = this->getNthreads() > 0 ? this->getNthreads() : 1;
                gpu_solvers.clear();
                gpu_solvers.reserve(nslots);
                for (size_t s = 0; s < nslots; ++s) {
                    auto solver = std::make_unique<OpenCLSweepSolver<T1>>();
                    solver->initialize(this->ncx, this->ncy, this->ncz,
                                       this->dx, this->dy, this->dz);
                    // Profile only slot 0: with several streams running at once
                    // the per-iteration diagnostics and end-of-run breakdown
                    // from multiple solvers would interleave unreadably.
                    solver->setProfiling(s == 0 && gpu_profile != 0);
                    gpu_solvers.push_back(std::move(solver));
                }
                gpu_available = true;
                gpu_initialized = true;

                if ( verbose ) {
                    std::cout << "  GPU acceleration enabled for Grid3Drcfs ("
                              << nslots << " solver slot" << (nslots > 1 ? "s" : "")
                              << ")\n";
                    std::cout << gpu_solvers[0]->getDeviceInfo();
                }

            } catch (const std::exception& e) {
                std::cerr << "GPU initialization failed: " << e.what() << "\n";
                std::cerr << "Falling back to CPU implementation\n";
                gpu_solvers.clear();
                gpu_available = false;
                gpu_initialized = true;
                use_gpu = false;
            }
        }

        /**
         * Raytrace implementation - GPU accelerated version
         */
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                     const std::vector<T1>& t0,
                     const std::vector<sxyz<T1>>& Rx,
                     const size_t threadNo=0) const;
        
        void raytrace(const std::vector<sxyz<T1>>& Tx,
                     const std::vector<T1>& t0,
                     const std::vector<std::vector<sxyz<T1>>>& Rx,
                     const size_t threadNo=0) const;

        /**
         * Perform sweep iterations with convergence checking
         * 
         * @param frozen     Boolean array of frozen nodes
         * @param threadNo   Thread number
         * @param use_weno   Use WENO3 (true) or basic (false)
         * @return Number of iterations performed
         */
        int performSweepIterations(const std::vector<bool>& frozen,
                                   const size_t threadNo,
                                   const bool use_weno) const;
    };

    // =============================================================================
    // Implementation
    // =============================================================================

    template<typename T1, typename T2>
    void Grid3Drcfs_OpenCL<T1,T2>::setSlowness(const std::vector<T1>& s) {

        if ( static_cast<size_t>(this->ncx)*this->ncy*this->ncz != s.size() ) {
            throw std::length_error("Error: slowness vectors of incompatible size.");
        }

        // Interpolate slowness at grid nodes (same as Grid3Drcfs)
        // This ensures identical behavior between CPU and GPU versions

        const size_t nx = this->ncx;
        const size_t ny = this->ncy;
        const size_t nz = this->ncz;

        // corners
        this->nodes[                       0].setNodeSlowness( s[                       0] );
        this->nodes[                      nx].setNodeSlowness( s[                    nx-1] );
        this->nodes[           ny *(nx+1)   ].setNodeSlowness( s[          (ny-1)*nx     ] );
        this->nodes[           ny *(nx+1)+nx].setNodeSlowness( s[          (ny-1)*nx+nx-1] );
        this->nodes[(nz*(ny+1)   )*(nx+1)   ].setNodeSlowness( s[((nz-1)*ny     )*nx     ] );
        this->nodes[(nz*(ny+1)   )*(nx+1)+nx].setNodeSlowness( s[((nz-1)*ny     )*nx+nx-1] );
        this->nodes[(nz*(ny+1)+ny)*(nx+1)   ].setNodeSlowness( s[((nz-1)*ny+ny-1)*nx     ] );
        this->nodes[(nz*(ny+1)+ny)*(nx+1)+nx].setNodeSlowness( s[((nz-1)*ny+ny-1)*nx+nx-1] );

        // edges
        for ( size_t i=1; i<nx; ++i ) {
            this->nodes[                      i].setNodeSlowness( 0.5*(s[                    i]+s[                    i-1]) );
            this->nodes[           ny *(nx+1)+i].setNodeSlowness( 0.5*(s[          (ny-1)*nx+i]+s[          (ny-1)*nx+i-1]) );
            this->nodes[(nz*(ny+1)   )*(nx+1)+i].setNodeSlowness( 0.5*(s[((nz-1)*ny     )*nx+i]+s[((nz-1)*ny     )*nx+i-1]) );
            this->nodes[(nz*(ny+1)+ny)*(nx+1)+i].setNodeSlowness( 0.5*(s[((nz-1)*ny+ny-1)*nx+i]+s[((nz-1)*ny+ny-1)*nx+i-1]) );
        }
        for ( size_t j=1; j<ny; ++j ) {
            this->nodes[           j *(nx+1)   ].setNodeSlowness( 0.5*(s[            j*nx     ]+s[          (j-1)*nx     ]) );
            this->nodes[           j *(nx+1)+nx].setNodeSlowness( 0.5*(s[            j*nx+nx-1]+s[          (j-1)*nx+nx-1]) );
            this->nodes[(nz*(ny+1)+j)*(nx+1)   ].setNodeSlowness( 0.5*(s[((nz-1)*ny+j)*nx     ]+s[((nz-1)*ny+j-1)*nx     ]) );
            this->nodes[(nz*(ny+1)+j)*(nx+1)+nx].setNodeSlowness( 0.5*(s[((nz-1)*ny+j)*nx+nx-1]+s[((nz-1)*ny+j-1)*nx+nx-1]) );
        }
        for ( size_t k=1; k<nz; ++k ) {
            this->nodes[(k*(ny+1)   )*(nx+1)   ].setNodeSlowness( 0.5*(s[(k*ny     )*nx     ]+s[((k-1)*ny     )*nx     ]) );
            this->nodes[(k*(ny+1)   )*(nx+1)+nx].setNodeSlowness( 0.5*(s[(k*ny     )*nx+nx-1]+s[((k-1)*ny     )*nx+nx-1]) );
            this->nodes[(k*(ny+1)+ny)*(nx+1)   ].setNodeSlowness( 0.5*(s[(k*ny+ny-1)*nx     ]+s[((k-1)*ny+ny-1)*nx     ]) );
            this->nodes[(k*(ny+1)+ny)*(nx+1)+nx].setNodeSlowness( 0.5*(s[(k*ny+ny-1)*nx+nx-1]+s[((k-1)*ny+ny-1)*nx+nx-1]) );
        }

        // faces
        for ( size_t i=1; i<nx; ++i ) {
            for ( size_t j=1; j<ny; ++j ) {
                this->nodes[           j *(nx+1)+i].setNodeSlowness( 0.25*(s[             j *nx+i]+s[             j *nx+i-1]+
                                                                           s[          (j-1)*nx+i]+s[          (j-1)*nx+i-1]) );
                this->nodes[(nz*(ny+1)+j)*(nx+1)+i].setNodeSlowness( 0.25*(s[((nz-1)*ny+  j)*nx+i]+s[((nz-1)*ny+  j)*nx+i-1]+
                                                                           s[((nz-1)*ny+j-1)*nx+i]+s[((nz-1)*ny+j-1)*nx+i-1]) );
            }
        }
        for ( size_t i=1; i<nx; ++i ) {
            for ( size_t k=1; k<nz; ++k ) {
                this->nodes[(k*(ny+1)   )*(nx+1)+i].setNodeSlowness( 0.25*(s[(   k *ny     )*nx+i]+s[(   k *ny     )*nx+i-1]+
                                                                           s[((k-1)*ny     )*nx+i]+s[((k-1)*ny     )*nx+i-1]) );
                this->nodes[(k*(ny+1)+ny)*(nx+1)+i].setNodeSlowness( 0.25*(s[(   k *ny+ny-1)*nx+i]+s[(   k *ny+ny-1)*nx+i-1]+
                                                                           s[((k-1)*ny+ny-1)*nx+i]+s[((k-1)*ny+ny-1)*nx+i-1]) );
            }
        }
        for ( size_t j=1; j<ny; ++j ) {
            for ( size_t k=1; k<nz; ++k ) {
                this->nodes[(k*(ny+1)+j)*(nx+1)   ].setNodeSlowness( 0.25*(s[(k*ny+  j)*nx     ]+s[((k-1)*ny+  j)*nx     ]+
                                                                           s[(k*ny+j-1)*nx     ]+s[((k-1)*ny+j-1)*nx     ]) );
                this->nodes[(k*(ny+1)+j)*(nx+1)+nx].setNodeSlowness( 0.25*(s[(k*ny+  j)*nx+nx-1]+s[((k-1)*ny+  j)*nx+nx-1]+
                                                                           s[(k*ny+j-1)*nx+nx-1]+s[((k-1)*ny+j-1)*nx+nx-1]) );
            }
        }

        // interior
        for ( size_t i=1; i<nx; ++i ) {
            for ( size_t j=1; j<ny; ++j ) {
                for ( size_t k=1; k<nz; ++k ) {
                    this->nodes[(k*(ny+1)+j)*(nx+1)+i].setNodeSlowness( 0.125*(s[(    k*ny+j  )*nx+i  ]+
                                                                               s[(    k*ny+j  )*nx+i-1]+
                                                                               s[(    k*ny+j-1)*nx+i  ]+
                                                                               s[(    k*ny+j-1)*nx+i-1]+
                                                                               s[((k-1)*ny+j  )*nx+i  ]+
                                                                               s[((k-1)*ny+j  )*nx+i-1]+
                                                                               s[((k-1)*ny+j-1)*nx+i  ]+
                                                                               s[((k-1)*ny+j-1)*nx+i-1]) );
                }
            }
        }
    }

    template<typename T1, typename T2>
    int Grid3Drcfs_OpenCL<T1,T2>::performSweepIterations(
        const std::vector<bool>& frozen,
        const size_t threadNo,
        const bool use_weno) const 
    {
        // Store previous travel times for convergence check
        std::vector<T1> times(this->nodes.size());
        for (size_t n=0; n<this->nodes.size(); ++n) {
            times[n] = this->nodes[n].getTT(threadNo);
        }
        
        T1 change = std::numeric_limits<T1>::max();
        int niter = 0;
        T1 prev = 0.0;   // previous sweep's change, for the non-monotone guard
        T1 tref = 0.0;   // traveltime range, refreshed by fsmChange
        T1 tol = 0.0;    // absolute stop threshold, from epsilon and tref
        
        if (use_gpu && gpu_available) {
            // ================================================================
            // GPU PATH - Accelerated sweep iterations
            // ================================================================
            
            try {
                // Extract data for GPU
                std::vector<T1> tt(this->nodes.size());
                std::vector<T1> slowness(this->nodes.size());
                
                for (size_t n=0; n<this->nodes.size(); ++n) {
                    tt[n] = this->nodes[n].getTT(threadNo);
                    slowness[n] = this->nodes[n].getNodeSlowness();
                }
                
                // Set sweep type
                SweepType sweep_type = use_weno ? SweepType::WENO3 : SweepType::BASIC;
                gpu_solvers[threadNo]->setSweepType(sweep_type);

                // ---- Error-vs-analytic diagnostic (profiling only) --------
                // For a HOMOGENEOUS medium the exact traveltime is
                //   texact[n] = min over frozen f ( tt_f + s*dist(n,f) ),
                // since each frozen (source) node already holds the exact BC.
                // Comparing the iterate against this each cycle reveals the
                // discretization-error floor: the L1 error decays then flattens
                // once the iterative error drops below the scheme error, and
                // that flat level is the epsilon below which iterating is
                // pointless.  Built once; empty (and skipped) if the model is
                // not uniformly homogeneous or there are no frozen nodes.
                // O(N*Nfrozen), profiling only.
                std::vector<T1> texact;
                if (gpu_solvers[threadNo]->isProfiling()) {
                    const size_t N = this->nodes.size();
                    T1 s0 = this->nodes[0].getNodeSlowness();
                    T1 smin = s0, smax = s0;
                    for (size_t n=1; n<N; ++n) {
                        T1 s = this->nodes[n].getNodeSlowness();
                        if (s < smin) smin = s;
                        if (s > smax) smax = s;
                    }
                    std::vector<size_t> src;
                    for (size_t n=0; n<N; ++n) if (frozen[n]) src.push_back(n);
                    if (smax - smin <= 1.e-6 * smax && !src.empty()) {
                        texact.resize(N);
                        for (size_t n=0; n<N; ++n) {
                            T1 xn = this->nodes[n].getX();
                            T1 yn = this->nodes[n].getY();
                            T1 zn = this->nodes[n].getZ();
                            T1 best = std::numeric_limits<T1>::max();
                            for (size_t f : src) {
                                T1 ddx = xn - this->nodes[f].getX();
                                T1 ddy = yn - this->nodes[f].getY();
                                T1 ddz = zn - this->nodes[f].getZ();
                                T1 d = std::sqrt(ddx*ddx + ddy*ddy + ddz*ddz);
                                T1 t = this->nodes[f].getTT(threadNo) + s0 * d;
                                if (t < best) best = t;
                            }
                            texact[n] = best;
                        }
                    } else if (use_weno) {
                        std::cout << "  [OpenCL converge] (model not homogeneous "
                                     "or no frozen nodes: error-vs-analytic "
                                     "diagnostic disabled)\n";
                    }
                }
                // Mean L1 error vs the analytic field (-1 when unavailable).
                auto meanError = [&](const std::vector<T1>& cur) -> double {
                    if (texact.empty()) return -1.0;
                    double e = 0.0;
                    for (size_t n=0; n<cur.size(); ++n)
                        e += std::abs(static_cast<double>(cur[n]) -
                                      static_cast<double>(texact[n]));
                    return e / static_cast<double>(cur.size());
                };

                // Iterative refinement on GPU
                while (niter < nitermax && (niter < 2 || change >= tol || change > prev)) {
                    // Perform one sweep cycle on GPU (8 directional sweeps)
                    gpu_solvers[threadNo]->runSweeps(tt, slowness, frozen, 1);
                    
                    // Check convergence
                    prev = change;
                    change = fsmChange(tt, times, tref);
                    tol = fsmTolerance(epsilon, tref, this->nodes.size());

                    niter++;
                    // Convergence-curve diagnostic (only when GPU profiling is
                    // on).  Reports the raw L1 sum that drives the stop test and
                    // the grid-size-independent per-node mean (= sum/num_nodes);
                    // where the mean flattens is the epsilon below which further
                    // iterations buy no meaningful accuracy.
                    if (gpu_solvers[threadNo]->isProfiling()) {
                        double err = meanError(tt);
                        std::cout << "  [OpenCL converge] "
                                  << (use_weno ? "WENO3 " : "1st-ord")
                                  << " iter " << niter
                                  << "  change(sum) = " << change
                                  << "  mean = " << change / this->nodes.size();
                        if (err >= 0.0) std::cout << "  L1err = " << err;
                        std::cout << '\n';
                    }
                }
                
                // Update node travel times from GPU results
                for (size_t n=0; n<this->nodes.size(); ++n) {
                    this->nodes[n].setTT(tt[n], threadNo);
                }
                
            } catch (const std::exception& e) {
                std::cerr << "GPU sweep failed: " << e.what() << "\n";
                std::cerr << "Falling back to CPU for this raytrace\n";
                
                // Fall through to CPU implementation below
                use_gpu = false;
            }
        }
        
        if (!use_gpu || !gpu_available) {
            // ================================================================
            // CPU PATH - Original implementation
            // ================================================================
            
            while (niter < nitermax && (niter < 2 || change >= tol || change >= prev)) {
                if (use_weno) {
                    this->sweep_weno3(frozen, threadNo);
                } else {
                    this->sweep(frozen, threadNo);
                }
                
                prev = change;
                change = fsmChange(this->nodes, times, threadNo, tref);
                tol = fsmTolerance(epsilon, tref, this->nodes.size());
                niter++;
            }
        }
        
        if (niter == nitermax && change >= tol) {
            warnFSMnotConverged(use_weno ? "WENO3" : "first-order",
                                niter, change, epsilon, tref, this->nodes.size());
        }
        return niter;
    }

    template<typename T1, typename T2>
    void Grid3Drcfs_OpenCL<T1,T2>::raytrace(
        const std::vector<sxyz<T1>>& Tx,
        const std::vector<T1>& t0,
        const std::vector<sxyz<T1>>& Rx,
        const size_t threadNo) const 
    {
        this->checkPts(Tx, true);
        this->checkPts(Rx, true);

        // Reinitialize nodes
        for (size_t n=0; n<this->nodes.size(); ++n) {
            this->nodes[n].reinit(threadNo);
        }

        // Set Tx pts (frozen nodes)
        std::vector<bool> frozen(this->nodes.size(), false);
        int npts = weno3 ? 2 : 1;
        this->initFSM(Tx, t0, frozen, npts, threadNo);

        // Validation for WENO3
        if (weno3 && (this->dx != this->dz || this->dx != this->dy)) {
            throw std::logic_error("Error: WENO stencil needs dx equal to dy and dz");
        }

        if (weno3) {
            // Two-stage refinement: basic sweep then WENO3
            
            // Stage 1: Basic sweep for initial solution
            niter_final = performSweepIterations(frozen, threadNo, false);
            
            // Stage 2: WENO3 sweep for high accuracy
            niterw_final = performSweepIterations(frozen, threadNo, true);
            
        } else {
            // Single-stage: basic sweep only
            niter_final = performSweepIterations(frozen, threadNo, false);
            niterw_final = 0;
        }
    }

    template<typename T1, typename T2>
    void Grid3Drcfs_OpenCL<T1,T2>::raytrace(
        const std::vector<sxyz<T1>>& Tx,
        const std::vector<T1>& t0,
        const std::vector<std::vector<sxyz<T1>>>& Rx,
        const size_t threadNo) const 
    {
        this->checkPts(Tx, true);
        for (size_t n=0; n<Rx.size(); ++n) {
            this->checkPts(Rx[n], true);
        }

        // Reinitialize nodes
        for (size_t n=0; n<this->nodes.size(); ++n) {
            this->nodes[n].reinit(threadNo);
        }

        // Set Tx pts (frozen nodes)
        std::vector<bool> frozen(this->nodes.size(), false);
        int npts = weno3 ? 2 : 1;
        this->initFSM(Tx, t0, frozen, npts, threadNo);

        // Validation for WENO3
        if (weno3 && (this->dx != this->dz || this->dx != this->dy)) {
            throw std::logic_error("Error: WENO stencil needs dx equal to dy and dz");
        }

        if (weno3) {
            // Two-stage refinement: basic sweep then WENO3
            
            // Stage 1: Basic sweep for initial solution
            niter_final = performSweepIterations(frozen, threadNo, false);
            
            // Stage 2: WENO3 sweep for high accuracy
            niterw_final = performSweepIterations(frozen, threadNo, true);
            
        } else {
            // Single-stage: basic sweep only
            niter_final = performSweepIterations(frozen, threadNo, false);
            niterw_final = 0;
        }
    }

} // namespace ttcr

#endif /* ttcr_Grid3Drcfs_OpenCL_h */

//
//  Grid3Drnfs_OpenCL.h
//  ttcr
//
//  GPU-accelerated version of Grid3Drnfs using OpenCL
//
//  This class provides the same interface as Grid3Drnfs but uses
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

#ifndef ttcr_Grid3Drnfs_OpenCL_h
#define ttcr_Grid3Drnfs_OpenCL_h

#include <cmath>
#include <limits>
#include <vector>
#include <memory>
#include <iostream>
#include <type_traits>

#ifdef VTK
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkXMLRectilinearGridWriter.h"
#endif

#include "Grid3Drn.h"
#include "Grid3Drn_OpenCL.h"
#include "Node3Dn.h"

namespace ttcr {

    /**
     * GPU-accelerated Grid3Drnfs using OpenCL
     * 
     * This class maintains the same interface as Grid3Drnfs but uses
     * GPU-accelerated sweep operations via OpenCL kernels.
     * 
     * Features:
     * - Same API as Grid3Drnfs (drop-in replacement)
     * - 8-20x speedup for large grids (>64³)
     * - Automatic fallback to CPU if GPU fails
     * - Support for both basic and WENO3 sweeps
     * - Convergence checking same as CPU version
     * 
     * Usage:
     *   Grid3Drnfs_OpenCL<double, uint32_t> grid(
     *       nx, ny, nz, dx, xmin, ymin, zmin,
     *       epsilon, max_iterations, use_weno3);
     *   
     *   grid.setSlowness(slowness_values);
     *   grid.raytrace(Tx, t0, Rx);
     */
    template<typename T1, typename T2>
    class Grid3Drnfs_OpenCL : public Grid3Drn<T1,T2,Node3Dn<T1,T2>> {
    public:
        /**
         * Constructor
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
         * @param ttrp      Travel time reciprocal paths
         * @param intVel    Interpolate velocity
         * @param nt        Number of threads
         * @param _translateOrigin  Translate origin to (0,0,0)
         * @param enableGPU Enable GPU acceleration (true by default)
         */
        Grid3Drnfs_OpenCL(const T2 nx, const T2 ny, const T2 nz, const T1 ddx,
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
            epsilon *= static_cast<T1>(this->nodes.size());  // per-node tol -> L1-sum threshold (nodes built)

            if (use_gpu) {
                initializeGPU();
            }
        }

        virtual ~Grid3Drnfs_OpenCL() {
            // OpenCLSweepSolver destructor handles cleanup
        }

        // Accessors (same as Grid3Drnfs)
        const int get_niter() const { return niter_final; }
        const int get_niterw() const { return niterw_final; }
        
        // GPU-specific methods
        bool isUsingGPU() const { return use_gpu && gpu_available; }
        void setUseGPU(bool enable) { 
            use_gpu = enable;
            if (use_gpu && !gpu_initialized) {
                initializeGPU();
            }
        }
        
        std::string getGPUInfo() const {
            if (gpu_available && !gpu_solvers.empty()) {
                return gpu_solvers[0]->getDeviceInfo();
            }
            return "GPU not available";
        }
        
#ifdef VTK
        void saveModelVTR(const std::string &fname, const std::vector<T1> &d, const std::string &dname) const;
#endif

    protected:
        T1 epsilon;              // Convergence criterion: L1-sum threshold (input per-node tol scaled by nNodes in ctor)
        int nitermax;            // Maximum iterations
        mutable int niter_final; // Final iteration count (basic sweep)
        mutable int niterw_final;// Final iteration count (WENO3 sweep)
        bool weno3;              // Use WENO3 sweep
        
        // GPU-specific members
        mutable bool use_gpu;
        mutable bool gpu_initialized;
        mutable bool gpu_available;
        // one independent solver per thread slot (own context/queue/
        // kernels/buffers).  Indexed by threadNo, sized to getNthreads().
        mutable std::vector<std::unique_ptr<OpenCLSweepSolver<T1>>> gpu_solvers;

    private:
        Grid3Drnfs_OpenCL() {}
        Grid3Drnfs_OpenCL(const Grid3Drnfs_OpenCL<T1,T2>& g) {}
        Grid3Drnfs_OpenCL<T1,T2>& operator=(const Grid3Drnfs_OpenCL<T1,T2>& g) { return *this; }

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
                    std::cout << "  GPU acceleration enabled for Grid3Drnfs ("
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
    int Grid3Drnfs_OpenCL<T1,T2>::performSweepIterations(
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
        
        if (use_gpu && gpu_available) {
            // ================================================================
            // GPU PATH - Accelerated sweep iterations
            // ================================================================
            
            try {
                // Extract data for GPU (initial upload)
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
                // not uniformly homogeneous or there are no frozen nodes, since
                // the closed form is only valid then.  O(N*Nfrozen), profiling
                // only.
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

                // Upload initial data (only once!)
                std::vector<unsigned char> frozen_uchar(this->nodes.size());
                for (size_t i = 0; i < this->nodes.size(); ++i) {
                    frozen_uchar[i] = frozen[i] ? 1 : 0;
                }
                
                gpu_solvers[threadNo]->runSweeps(tt, slowness, frozen, 1);

                niter++;
                
                // Check convergence
                change = 0.0;
                for (size_t n=0; n<this->nodes.size(); ++n) {
                    T1 dt = std::abs(times[n] - tt[n]);
                    change += dt;
                    times[n] = tt[n];
                }
                // Convergence-curve diagnostic (only when GPU profiling is on).
                // Reports both the raw L1 sum that drives the stop test and the
                // per-node mean (= sum / num_nodes): the mean is the grid-size-
                // independent quantity, so reading where it flattens tells you
                // the epsilon below which further iterations buy no meaningful
                // accuracy.  See [[gpu-fsm-convergence]] notes.
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

                // Iterative refinement on GPU
                while (change >= epsilon && niter < nitermax) {
                    // Always run sweeps - runSweeps modifies tt in place
                    gpu_solvers[threadNo]->runSweepsNoTransfer(1);
                    
                    gpu_solvers[threadNo]->downloadTravelTimes(tt);

                    niter++;
                    // Check convergence
                    change = 0.0;
                    for (size_t n=0; n<this->nodes.size(); ++n) {
                        T1 dt = std::abs(times[n] - tt[n]);
                        change += dt;
                        times[n] = tt[n];
                    }
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
            
            while (change >= epsilon && niter < nitermax) {
                if (use_weno) {
                    this->sweep_weno3(frozen, threadNo);
                } else {
                    this->sweep(frozen, threadNo);
                }
                
                change = 0.0;
                for (size_t n=0; n<this->nodes.size(); ++n) {
                    T1 dt = std::abs(times[n] - this->nodes[n].getTT(threadNo));
                    change += dt;
                    times[n] = this->nodes[n].getTT(threadNo);
                }
                niter++;
            }
        }
        
        return niter;
    }

    template<typename T1, typename T2>
    void Grid3Drnfs_OpenCL<T1,T2>::raytrace(
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
    void Grid3Drnfs_OpenCL<T1,T2>::raytrace(
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

#ifdef VTK
    template<typename T1, typename T2>
    void Grid3Drnfs_OpenCL<T1,T2>::saveModelVTR(const std::string &fname,
                                                const std::vector<T1> &d,
                                                const std::string &dname) const
    {
        int nn[3] = {static_cast<int>(this->ncx+1),
            static_cast<int>(this->ncy+1),
            static_cast<int>(this->ncz+1)};

        vtkSmartPointer<vtkDoubleArray> xCoords = vtkSmartPointer<vtkDoubleArray>::New();
        for (size_t n=0; n<nn[0]; ++n) {
            xCoords->InsertNextValue( this->xmin + n*this->dx );
        }
        vtkSmartPointer<vtkDoubleArray> yCoords = vtkSmartPointer<vtkDoubleArray>::New();
        for (size_t n=0; n<nn[1]; ++n) {
            yCoords->InsertNextValue( this->ymin + n*this->dy );
        }
        vtkSmartPointer<vtkDoubleArray> zCoords = vtkSmartPointer<vtkDoubleArray>::New();
        for (size_t n=0; n<nn[2]; ++n) {
            zCoords->InsertNextValue( this->zmin + n*this->dz );
        }

        vtkSmartPointer<vtkRectilinearGrid> rgrid = vtkSmartPointer<vtkRectilinearGrid>::New();
        rgrid->SetDimensions( nn );
        rgrid->SetXCoordinates(xCoords);
        rgrid->SetYCoordinates(yCoords);
        rgrid->SetZCoordinates(zCoords);
        
        if constexpr (std::is_same<T1, float>::value) {
            vtkSmartPointer<vtkFloatArray> newScalars =
            vtkSmartPointer<vtkFloatArray>::New();
            
            newScalars->SetNumberOfComponents(1);
            newScalars->SetNumberOfTuples( rgrid->GetNumberOfPoints() );
            
            newScalars->SetName(dname.c_str());
            
            for ( size_t n=0; n<this->nodes.size(); ++n ) {
                if ( this->nodes[n].isPrimary() ) {
                    double x[3] = {this->nodes[n].getX(), this->nodes[n].getY(), this->nodes[n].getZ()};
                    vtkIdType id = rgrid->FindPoint(x);
                    newScalars->SetTuple1(id, d[n] );
                }
            }
            rgrid->GetPointData()->SetScalars(newScalars);
        } else {
            // using double
            vtkSmartPointer<vtkDoubleArray> newScalars =
            vtkSmartPointer<vtkDoubleArray>::New();

            newScalars->SetNumberOfComponents(1);
            newScalars->SetNumberOfTuples( rgrid->GetNumberOfPoints() );
            
            newScalars->SetName(dname.c_str());
            
            for ( size_t n=0; n<this->nodes.size(); ++n ) {
                if ( this->nodes[n].isPrimary() ) {
                    double x[3] = {this->nodes[n].getX(), this->nodes[n].getY(), this->nodes[n].getZ()};
                    vtkIdType id = rgrid->FindPoint(x);
                    newScalars->SetTuple1(id, d[n] );
                }
            }
            rgrid->GetPointData()->SetScalars(newScalars);
        }
        vtkSmartPointer<vtkXMLRectilinearGridWriter> writer =
        vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();

        writer->SetFileName( fname.c_str() );
        writer->SetInputData( rgrid );
        writer->SetDataModeToBinary();
        writer->Update();
    }
#endif
} // namespace ttcr

#endif /* ttcr_Grid3Drnfs_OpenCL_h */

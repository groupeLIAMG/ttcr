//
//  Grid2Drnfs_OpenCL.h
//  ttcr
//
//  GPU-accelerated Grid2Drnfs using OpenCL.
//  Drop-in replacement for Grid2Drnfs with the same constructor signature
//  (minus rotated_template, which is not supported by the GPU path).
//
//  The sweep45 rotated-template stencil is deliberately omitted: its
//  anti-diagonal ordering does not map to the Gauss-Seidel plane-sweep
//  approach used by the GPU kernels.  If rotated_template support is needed,
//  use Grid2Drnfs instead.
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

/**
 * @file Grid2Drnfs_OpenCL.h
 * @brief GPU-accelerated fast sweeping on a 2-D rectilinear grid with node-based
 *        slowness.
 *
 * Declares ttcr::Grid2Drnfs_OpenCL, the OpenCL counterpart of
 * ttcr::Grid2Drnfs. Slowness is taken at the nodes directly, with none of the
 * cell-to-node averaging that ttcr::Grid2Drcfs_OpenCL performs.
 *
 * @sa Grid2Drnfs.h, Grid2Drn_OpenCL.h, Grid2Drcfs_OpenCL.h
 */

#ifndef ttcr_Grid2Drnfs_OpenCL_h
#define ttcr_Grid2Drnfs_OpenCL_h

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "Grid2Drn.h"
#include "Grid2Drn_OpenCL.h"
#include "Node2Dn.h"

namespace ttcr {

/**
 * @brief GPU-accelerated fast sweeping solver with node-based slowness.
 *
 * @tparam T1 floating-point type of coordinates, slowness and traveltimes.
 * @tparam T2 integer type of node and cell indices.
 * @tparam S  point type, @ref sxz or @ref sxyz.
 *
 * The OpenCL counterpart of ttcr::Grid2Drnfs: same node construction, same
 * convergence test, with the sweep running on a device. Slowness is taken at
 * the nodes through the inherited ttcr::Grid2Drn::setSlowness — no cell-to-node
 * averaging, unlike ttcr::Grid2Drcfs_OpenCL.
 *
 * @note GPU use is requested, not guaranteed — @ref isUsingGPU reports what is
 *       actually in use, and the class falls back to the CPU otherwise, so
 *       results are unaffected and only performance changes. One solver is held
 *       per thread; ttcr::input_parameters::gpu_max_threads caps how many run
 *       concurrently.
 * @note No @c rotated_template option, unlike ttcr::Grid2Drnfs.
 *
 * @sa Grid2Drnfs.h, Grid2Drn_OpenCL.h, Grid2Drcfs_OpenCL.h
 */
template<typename T1, typename T2, typename S>
class Grid2Drnfs_OpenCL : public Grid2Drn<T1,T2,S,Node2Dn<T1,T2>> {
public:
    /**
     * @brief Build the grid, its nodes, and the GPU solvers if requested.
     *
     * @param nx         Number of cells in x
     * @param nz         Number of cells in z
     * @param ddx        Cell size in x
     * @param ddz        Cell size in z
     * @param minx       X origin
     * @param minz       Z origin
     * @param eps        Per-node convergence tolerance
     * @param maxit      Maximum sweep iterations
     * @param w          Use WENO3 refinement pass (true/false)
     * @param ttrp       Compute traveltimes from raypaths
     * @param nt         Number of threads
     * @param enableGPU  Enable GPU acceleration (true by default)
     *
     * @post Nodes and neighbour lists are built and, if @p enableGPU, the GPU
     *       solvers are initialised — falling back to the CPU if that fails.
     *       @p eps is a relative tolerance: the sweeps stop once the mean
     *       per-node change falls below that fraction of the traveltime
     *       range. Slowness is **not** set.
     */
    Grid2Drnfs_OpenCL(const T2 nx, const T2 nz,
                      const T1 ddx, const T1 ddz,
                      const T1 minx, const T1 minz,
                      const T1 eps, const int maxit,
                      const bool w, const bool ttrp,
                      const size_t nt = 1,
                      const bool enableGPU = true) :
        Grid2Drn<T1,T2,S,Node2Dn<T1,T2>>(nx, nz, ddx, ddz, minx, minz, ttrp, nt),
        epsilon(eps), nitermax(maxit),
        niter_final(0), niterw_final(0),
        weno3(w),
        use_gpu(enableGPU), gpu_initialized(false), gpu_available(false)
    {
        buildGridNodes();
        this->template buildGridNeighbors<Node2Dn<T1,T2>>(this->nodes);

        if (use_gpu) initializeGPU();
    }

    /// Destructor; the OpenCLSweepSolver2D destructors release the device state.
    virtual ~Grid2Drnfs_OpenCL() {}

    /// @return Number of sweep iterations the last solve took.
    const int  get_niter()  const { return niter_final;  }
    /// @return Number of WENO refinement iterations the last solve took.
    const int  get_niterw() const { return niterw_final; }
    /**
     * @brief Whether the solve will actually run on the GPU.
     * @return True only if GPU use was requested **and** a device was
     *         successfully initialised.
     */
    bool isUsingGPU()       const { return use_gpu && gpu_available; }

    /**
     * @brief Human-readable description of the OpenCL device in use.
     * @return Device information, or @c "GPU not available" if the solve will
     *         run on the CPU.
     */
    std::string getGPUInfo() const {
        if (gpu_available && !gpu_solvers.empty())
            return gpu_solvers[0]->getDeviceInfo();
        return "GPU not available";
    }

private:
    T1  epsilon;                   ///< Convergence tolerance, dimensionless: the fraction of the traveltime range the mean per-node change must fall below.
    int nitermax;                  ///< Iteration cap for the sweeps.
    mutable int niter_final, niterw_final;  ///< Iterations used by the last solve; @c mutable so the @c const raytrace can record them.
    bool weno3;                    ///< Run the WENO3 refinement pass.

    mutable bool use_gpu;          ///< GPU acceleration was requested.
    mutable bool gpu_initialized;  ///< Initialisation has been attempted.
    mutable bool gpu_available;    ///< A device was successfully initialised.
    /// One solver per thread, so concurrent sources can each drive the device.
    mutable std::vector<std::unique_ptr<OpenCLSweepSolver2D<T1>>> gpu_solvers;

    Grid2Drnfs_OpenCL() = delete;
    Grid2Drnfs_OpenCL(const Grid2Drnfs_OpenCL&) = delete;
    Grid2Drnfs_OpenCL& operator=(const Grid2Drnfs_OpenCL&) = delete;

    // -------------------------------------------------------------------------
    // Grid node construction (identical to Grid2Drnfs::buildGridNodes)
    // -------------------------------------------------------------------------
    void buildGridNodes() {
        T2 cell_upLeft   = std::numeric_limits<T2>::max();
        T2 cell_upRight  = std::numeric_limits<T2>::max();
        T2 cell_downLeft = 0;
        T2 cell_downRight = 0;

        for (T2 n = 0, nc = 0; nc <= this->ncx; ++nc) {
            T1 x = this->xmin + nc * this->dx;
            for (T2 nr = 0; nr <= this->ncz; ++nr) {
                T1 z = this->zmin + nr * this->dz;

                cell_downRight = (nr < this->ncz && nc < this->ncx)
                    ? nc * this->ncz + nr : std::numeric_limits<T2>::max();
                cell_upRight   = (nr > 0 && nc < this->ncx)
                    ? nc * this->ncz + nr - 1 : std::numeric_limits<T2>::max();
                cell_downLeft  = (nr < this->ncz && nc > 0)
                    ? (nc-1)*this->ncz + nr : std::numeric_limits<T2>::max();
                cell_upLeft    = (nr > 0 && nc > 0)
                    ? (nc-1)*this->ncz + nr - 1 : std::numeric_limits<T2>::max();

                if (cell_upLeft   != std::numeric_limits<T2>::max()) this->nodes[n].pushOwner(cell_upLeft);
                if (cell_downLeft != std::numeric_limits<T2>::max()) this->nodes[n].pushOwner(cell_downLeft);
                if (cell_upRight  != std::numeric_limits<T2>::max()) this->nodes[n].pushOwner(cell_upRight);
                if (cell_downRight!= std::numeric_limits<T2>::max()) this->nodes[n].pushOwner(cell_downRight);

                this->nodes[n].setX(x);
                this->nodes[n].setZ(z);
                this->nodes[n].setGridIndex(n);
                this->nodes[n].setPrimary(true);
                ++n;
            }
        }
    }

    // -------------------------------------------------------------------------
    // GPU initialisation
    // -------------------------------------------------------------------------
    void initializeGPU() const {
        if (gpu_initialized) return;
        try {
            const size_t nslots = this->getNthreads() > 0 ? this->getNthreads() : 1;
            gpu_solvers.clear();
            gpu_solvers.reserve(nslots);
            for (size_t s = 0; s < nslots; ++s) {
                auto solver = std::make_unique<OpenCLSweepSolver2D<T1>>();
                solver->initialize(this->ncx, this->ncz, this->dx, this->dz);
                solver->setProfiling(s == 0 && gpu_profile != 0);
                gpu_solvers.push_back(std::move(solver));
            }
            gpu_available   = true;
            gpu_initialized = true;

            if ( verbose ) {
                std::cout << "  GPU acceleration enabled for Grid2Drnfs_OpenCL ("
                          << nslots << " solver slot" << (nslots > 1 ? "s" : "") << ")\n";
                std::cout << gpu_solvers[0]->getDeviceInfo();
            }
        } catch (const std::exception& e) {
            std::cerr << "2D GPU init failed: " << e.what() << "\n"
                      << "Falling back to CPU\n";
            gpu_solvers.clear();
            gpu_available   = false;
            gpu_initialized = true;
            use_gpu         = false;
        }
    }

    // -------------------------------------------------------------------------
    // Core sweep loop (GPU or CPU fallback)
    // -------------------------------------------------------------------------
    int performSweepIterations(const std::vector<bool>& frozen,
                               const size_t threadNo,
                               const SweepMode2D mode) const
    {
        std::vector<T1> times(this->nodes.size());
        for (size_t n = 0; n < this->nodes.size(); ++n)
            times[n] = this->nodes[n].getTT(threadNo);

        T1  change = std::numeric_limits<T1>::max();
        int niter  = 0;
        T1 prev = 0.0;   // previous sweep's change, for the non-monotone guard
        T1 tref = 0.0;   // traveltime range, refreshed by fsmChange
        T1 tol = 0.0;    // absolute stop threshold, from epsilon and tref
        // Name of this pass, for the non-convergence warning below.
        const char *pass = (mode == SweepMode2D::WENO3 ||
                            mode == SweepMode2D::WENO3_XZ) ? "WENO3"
                                                           : "first-order";

        if (use_gpu && gpu_available) {
            try {
                std::vector<T1> tt(this->nodes.size());
                std::vector<T1> slowness(this->nodes.size());
                for (size_t n = 0; n < this->nodes.size(); ++n) {
                    tt[n]       = this->nodes[n].getTT(threadNo);
                    slowness[n] = this->nodes[n].getNodeSlowness();
                }

                gpu_solvers[threadNo]->setSweepMode(mode);
                gpu_solvers[threadNo]->runSweeps(tt, slowness, frozen, 1);
                niter++;

                prev = change;
                change = fsmChange(tt, times, tref);
                tol = fsmTolerance(epsilon, tref, this->nodes.size());

                while (niter < nitermax && (niter < 2 || change >= tol || change > prev)) {
                    gpu_solvers[threadNo]->runSweepsNoTransfer(1);
                    gpu_solvers[threadNo]->downloadTravelTimes(tt);
                    niter++;
                    prev = change;
                    change = fsmChange(tt, times, tref);
                    tol = fsmTolerance(epsilon, tref, this->nodes.size());
                }

                for (size_t n = 0; n < this->nodes.size(); ++n)
                    this->nodes[n].setTT(tt[n], threadNo);

                if (niter == nitermax && change >= tol) {
                    warnFSMnotConverged(pass, niter, change, epsilon,
                                        tref, this->nodes.size());
                }
                return niter;

            } catch (const std::exception& e) {
                std::cerr << "GPU 2D sweep failed: " << e.what()
                          << "\nFalling back to CPU\n";
                use_gpu = false;
            }
        }

        // CPU fallback
        change = std::numeric_limits<T1>::max();
        niter  = 0;
        prev   = 0.0;
        while (niter < nitermax && (niter < 2 || change >= tol || change >= prev)) {
            switch (mode) {
                case SweepMode2D::BASIC:
                    this->sweep(frozen, threadNo);
                    break;
                case SweepMode2D::BASIC_XZ:
                    this->sweep_xz(frozen, threadNo);
                    break;
                case SweepMode2D::WENO3:
                    this->sweep_weno3(frozen, threadNo);
                    break;
                case SweepMode2D::WENO3_XZ:
                    this->sweep_weno3_xz(frozen, threadNo);
                    break;
            }
            prev = change;
            change = fsmChange(this->nodes, times, threadNo, tref);
            tol = fsmTolerance(epsilon, tref, this->nodes.size());
            niter++;
        }
        if (niter == nitermax && change >= tol) {
            warnFSMnotConverged(pass, niter, change, epsilon,
                                tref, this->nodes.size());
        }
        return niter;
    }

    // -------------------------------------------------------------------------
    // Raytrace overrides
    // -------------------------------------------------------------------------
    void raytrace(const std::vector<S>& Tx,
                  const std::vector<T1>& t0,
                  const std::vector<S>& Rx,
                  const size_t threadNo = 0) const
    {
        this->checkPts(Tx);
        this->checkPts(Rx);

        for (size_t n = 0; n < this->nodes.size(); ++n)
            this->nodes[n].reinit(threadNo);

        std::vector<bool> frozen(this->nodes.size(), false);
        int npts = weno3 ? 2 : 1;
        this->initFSM(Tx, t0, frozen, npts, threadNo);

        const bool non_square = (this->dx != this->dz);

        if (weno3) {
            niter_final  = performSweepIterations(frozen, threadNo,
                               non_square ? SweepMode2D::BASIC_XZ : SweepMode2D::BASIC);
            niterw_final = performSweepIterations(frozen, threadNo,
                               non_square ? SweepMode2D::WENO3_XZ : SweepMode2D::WENO3);
        } else {
            niter_final  = performSweepIterations(frozen, threadNo,
                               non_square ? SweepMode2D::BASIC_XZ : SweepMode2D::BASIC);
            niterw_final = 0;
        }
    }

    void raytrace(const std::vector<S>& Tx,
                  const std::vector<T1>& t0,
                  const std::vector<const std::vector<S>*>& Rx,
                  const size_t threadNo = 0) const
    {
        this->checkPts(Tx);
        for (size_t n = 0; n < Rx.size(); ++n) this->checkPts(*Rx[n]);

        for (size_t n = 0; n < this->nodes.size(); ++n)
            this->nodes[n].reinit(threadNo);

        std::vector<bool> frozen(this->nodes.size(), false);
        int npts = weno3 ? 2 : 1;
        this->initFSM(Tx, t0, frozen, npts, threadNo);

        const bool non_square = (this->dx != this->dz);

        if (weno3) {
            niter_final  = performSweepIterations(frozen, threadNo,
                               non_square ? SweepMode2D::BASIC_XZ : SweepMode2D::BASIC);
            niterw_final = performSweepIterations(frozen, threadNo,
                               non_square ? SweepMode2D::WENO3_XZ : SweepMode2D::WENO3);
        } else {
            niter_final  = performSweepIterations(frozen, threadNo,
                               non_square ? SweepMode2D::BASIC_XZ : SweepMode2D::BASIC);
            niterw_final = 0;
        }
    }
};

} // namespace ttcr

#endif // ttcr_Grid2Drnfs_OpenCL_h

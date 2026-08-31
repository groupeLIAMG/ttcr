//
//  Grid2Drcfs_OpenCL.h
//  ttcr
//
//  GPU-accelerated Grid2Drcfs using OpenCL.
//  Cell-centered slowness variant: setSlowness() interpolates cell-averaged
//  values to the grid nodes (same scheme as Grid2Drcfs::setSlowness), then
//  the sweep operations run on the GPU exactly as in Grid2Drnfs_OpenCL.
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

#ifndef ttcr_Grid2Drcfs_OpenCL_h
#define ttcr_Grid2Drcfs_OpenCL_h

/**
 * @file Grid2Drcfs_OpenCL.h
 * @brief GPU-accelerated fast sweeping on a 2-D rectilinear grid, from a cell
 *        slowness model.
 *
 * Declares ttcr::Grid2Drcfs_OpenCL, the OpenCL counterpart of
 * ttcr::Grid2Drcfs. The class structure is identical — including deriving from
 * ttcr::Grid2Drn rather than ttcr::Grid2Drc, and averaging the supplied cell
 * slowness onto the nodes — and only the sweep itself moves to the device.
 *
 * Selected by ttcr::input_parameters::method @c == @c FAST_SWEEPING_OPENCL.
 * GPU use is a request, not a guarantee; see @ref g2drcfsocl_fallback.
 *
 * @sa Grid2Drcfs.h, Grid2Drn_OpenCL.h, Grid2Drc.h
 */

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#include "Grid2Drn.h"
#include "Grid2Drn_OpenCL.h"
#include "Node2Dn.h"

namespace ttcr {

/**
 * @brief GPU-accelerated fast sweeping solver taking a cell slowness model.
 *
 * @tparam T1 floating-point type of coordinates, slowness and traveltimes.
 * @tparam T2 integer type of node and cell indices.
 * @tparam S  point type, @ref sxz or @ref sxyz.
 *
 * The OpenCL counterpart of ttcr::Grid2Drcfs. It shares that class's structure
 * exactly — same base (ttcr::Grid2Drn, not ttcr::Grid2Drc; see
 * @ref g2drcfs_hybrid), same cell-to-node slowness averaging in
 * @ref setSlowness, same node construction — and differs only in running the
 * sweep on a GPU device.
 *
 * @section g2drcfsocl_fallback GPU availability
 * GPU use is requested at construction and may not be granted: the device may
 * be missing, or initialisation may fail. @ref isUsingGPU reports what is
 * actually in use and @ref getGPUInfo names the device. When the GPU is
 * unavailable the class falls back to the CPU path, so results are unaffected
 * and only performance changes.
 *
 * One solver object is held per thread, so several sources can be pushed to the
 * device concurrently; ttcr::input_parameters::gpu_max_threads caps how many.
 *
 * @note Unlike ttcr::Grid2Drcfs, this class has no @c rotated_template option,
 *       and it properly @c = @c delete s its default and copy operations rather
 *       than using the private-and-undefined idiom.
 *
 * @sa Grid2Drcfs.h, Grid2Drn_OpenCL.h
 */
template<typename T1, typename T2, typename S>
class Grid2Drcfs_OpenCL : public Grid2Drn<T1,T2,S,Node2Dn<T1,T2>> {
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
     * @param w          Use WENO3 refinement pass
     * @param ttrp       Compute traveltimes from raypaths
     * @param nt         Number of threads
     * @param enableGPU  Enable GPU acceleration (true by default)
     *
     * @post Nodes and neighbour lists are built and, if @p enableGPU, the GPU
     *       solvers are initialised — silently falling back to the CPU if that
     *       fails. @p eps is a relative tolerance: the sweeps stop once the
     *       mean per-node change falls below that fraction of the traveltime
     *       range. Slowness is **not** set.
     */
    Grid2Drcfs_OpenCL(const T2 nx, const T2 nz,
                      const T1 ddx, const T1 ddz,
                      const T1 minx, const T1 minz,
                      const T1 eps, const int maxit,
                      const bool w, const bool ttrp,
                      const size_t nt = 1,
                      const bool enableGPU = true) :
        Grid2Drn<T1,T2,S,Node2Dn<T1,T2>>(nx, nz, ddx, ddz, minx, minz, ttrp, nt),
        epsilon(eps), nitermax(maxit),
        niter_final(0), niterw_final(0),
        weno3(w), hasCellSlown(false),
        use_gpu(enableGPU), gpu_initialized(false), gpu_available(false)
    {
        buildGridNodes();
        this->template buildGridNeighbors<Node2Dn<T1,T2>>(this->nodes);

        if (use_gpu) initializeGPU();
    }

    virtual ~Grid2Drcfs_OpenCL() {}

    // -------------------------------------------------------------------------
    // Cell-slowness interface (same interpolation as Grid2Drcfs::setSlowness)
    // -------------------------------------------------------------------------
    /**
     * @brief Set the cell slowness model and average it onto the nodes.
     * @param s one slowness per cell, @f$n_{cx}n_{cz}@f$ values.
     * @throws std::length_error if @p s has the wrong size.
     * @note Uses the same averaging scheme as ttcr::Grid2Drcfs::setSlowness —
     *       corner nodes take one cell, edge nodes the mean of two, interior
     *       nodes the mean of four. @sa @ref g2drcfs_hybrid
     */
    void setSlowness(const std::vector<T1>& s) {
        if (static_cast<size_t>(this->ncx) * this->ncz != s.size())
            throw std::length_error("Error: slowness vectors of incompatible size.");

        cell_slowness = s;
        hasCellSlown  = true;

        const size_t nx = this->ncx;
        const size_t nz = this->ncz;

        // Four corners
        this->nodes[0].setNodeSlowness(s[0]);
        this->nodes[nz].setNodeSlowness(s[nz-1]);
        this->nodes[nx*(nz+1)].setNodeSlowness(s[nz*(nx-1)]);
        this->nodes[(nx+1)*(nz+1)-1].setNodeSlowness(s[nx*nz-1]);

        // Left and right edges
        for (size_t j = 1; j < nz; ++j) {
            this->nodes[j].setNodeSlowness(0.5*(s[j]+s[j-1]));
            this->nodes[nx*(nz+1)+j].setNodeSlowness(
                0.5*(s[nz*(nx-1)+j]+s[nz*(nx-1)+j-1]));
        }
        // Top and bottom edges
        for (size_t i = 1; i < nx; ++i) {
            this->nodes[i*(nz+1)].setNodeSlowness(0.5*(s[i*nz]+s[(i-1)*nz]));
            this->nodes[i*(nz+1)+nz].setNodeSlowness(
                0.5*(s[(i+1)*nz-1]+s[i*nz-1]));
        }
        // Interior nodes
        for (size_t i = 1; i < nx; ++i)
            for (size_t j = 1; j < nz; ++j)
                this->nodes[i*(nz+1)+j].setNodeSlowness(
                    0.25*(s[i*nz+j]+s[i*nz+j-1]+s[(i-1)*nz+j]+s[(i-1)*nz+j-1]));
    }

    /**
     * @brief Retrieve the cell slowness model.
     * @param[out] s the values passed to @ref setSlowness.
     * @note Returns the original cell values, not the averaged nodal ones.
     */
    void getSlowness(std::vector<T1>& s) const { s = cell_slowness; }

    /// @return Number of sweep iterations the last solve took.
    const int  get_niter()  const { return niter_final;  }
    /// @return Number of WENO refinement iterations the last solve took.
    const int  get_niterw() const { return niterw_final; }
    /// @return True once @ref setSlowness has been called.
    const bool hasCellSlowness() const { return hasCellSlown; }
    /**
     * @brief Whether the solve will actually run on the GPU.
     * @return True only if GPU use was requested **and** a device was
     *         successfully initialised. @sa @ref g2drcfsocl_fallback
     */
    bool isUsingGPU() const { return use_gpu && gpu_available; }

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
    bool hasCellSlown;             ///< Whether @ref cell_slowness has been populated.
    std::vector<T1> cell_slowness; ///< Copy of the cell slowness model as supplied.

    mutable bool use_gpu;          ///< GPU acceleration was requested.
    mutable bool gpu_initialized;  ///< Initialisation has been attempted.
    mutable bool gpu_available;    ///< A device was successfully initialised. @sa @ref g2drcfsocl_fallback
    /// One solver per thread, so concurrent sources can each drive the device.
    mutable std::vector<std::unique_ptr<OpenCLSweepSolver2D<T1>>> gpu_solvers;

    Grid2Drcfs_OpenCL() = delete;
    Grid2Drcfs_OpenCL(const Grid2Drcfs_OpenCL&) = delete;
    Grid2Drcfs_OpenCL& operator=(const Grid2Drcfs_OpenCL&) = delete;

    // -------------------------------------------------------------------------
    // Grid node construction (identical to Grid2Drcfs::buildGridNodes)
    // -------------------------------------------------------------------------
    void buildGridNodes() {
        T2 cell_upLeft    = std::numeric_limits<T2>::max();
        T2 cell_upRight   = std::numeric_limits<T2>::max();
        T2 cell_downLeft  = 0;
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

                if (cell_upLeft    != std::numeric_limits<T2>::max()) this->nodes[n].pushOwner(cell_upLeft);
                if (cell_downLeft  != std::numeric_limits<T2>::max()) this->nodes[n].pushOwner(cell_downLeft);
                if (cell_upRight   != std::numeric_limits<T2>::max()) this->nodes[n].pushOwner(cell_upRight);
                if (cell_downRight != std::numeric_limits<T2>::max()) this->nodes[n].pushOwner(cell_downRight);

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
                std::cout << "  GPU acceleration enabled for Grid2Drcfs_OpenCL ("
                          << nslots << " solver slot" << (nslots > 1 ? "s" : "") << ")\n";
                std::cout << gpu_solvers[0]->getDeviceInfo();
            }
        } catch (const std::exception& e) {
            std::cerr << "2D GPU init failed: " << e.what() << "\nFalling back to CPU\n";
            gpu_solvers.clear();
            gpu_available   = false;
            gpu_initialized = true;
            use_gpu         = false;
        }
    }

    // -------------------------------------------------------------------------
    // Core sweep loop
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
                case SweepMode2D::BASIC:     this->sweep(frozen, threadNo);         break;
                case SweepMode2D::BASIC_XZ:  this->sweep_xz(frozen, threadNo);      break;
                case SweepMode2D::WENO3:     this->sweep_weno3(frozen, threadNo);    break;
                case SweepMode2D::WENO3_XZ:  this->sweep_weno3_xz(frozen, threadNo); break;
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

#endif // ttcr_Grid2Drcfs_OpenCL_h

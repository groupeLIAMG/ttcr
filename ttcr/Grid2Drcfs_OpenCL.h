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

template<typename T1, typename T2, typename S>
class Grid2Drcfs_OpenCL : public Grid2Drn<T1,T2,S,Node2Dn<T1,T2>> {
public:
    /**
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
        epsilon *= static_cast<T1>(this->nodes.size());

        if (use_gpu) initializeGPU();
    }

    virtual ~Grid2Drcfs_OpenCL() {}

    // -------------------------------------------------------------------------
    // Cell-slowness interface (same interpolation as Grid2Drcfs::setSlowness)
    // -------------------------------------------------------------------------
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

    void getSlowness(std::vector<T1>& s) const { s = cell_slowness; }

    const int  get_niter()  const { return niter_final;  }
    const int  get_niterw() const { return niterw_final; }
    const bool hasCellSlowness() const { return hasCellSlown; }
    bool isUsingGPU() const { return use_gpu && gpu_available; }

    std::string getGPUInfo() const {
        if (gpu_available && !gpu_solvers.empty())
            return gpu_solvers[0]->getDeviceInfo();
        return "GPU not available";
    }

private:
    T1  epsilon;
    int nitermax;
    mutable int niter_final, niterw_final;
    bool weno3;
    bool hasCellSlown;
    std::vector<T1> cell_slowness;

    mutable bool use_gpu, gpu_initialized, gpu_available;
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

                change = 0;
                for (size_t n = 0; n < this->nodes.size(); ++n) {
                    change  += std::abs(times[n] - tt[n]);
                    times[n] = tt[n];
                }

                while (change >= epsilon && niter < nitermax) {
                    gpu_solvers[threadNo]->runSweepsNoTransfer(1);
                    gpu_solvers[threadNo]->downloadTravelTimes(tt);
                    niter++;
                    change = 0;
                    for (size_t n = 0; n < this->nodes.size(); ++n) {
                        change  += std::abs(times[n] - tt[n]);
                        times[n] = tt[n];
                    }
                }

                for (size_t n = 0; n < this->nodes.size(); ++n)
                    this->nodes[n].setTT(tt[n], threadNo);

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
        while (change >= epsilon && niter < nitermax) {
            switch (mode) {
                case SweepMode2D::BASIC:     this->sweep(frozen, threadNo);         break;
                case SweepMode2D::BASIC_XZ:  this->sweep_xz(frozen, threadNo);      break;
                case SweepMode2D::WENO3:     this->sweep_weno3(frozen, threadNo);    break;
                case SweepMode2D::WENO3_XZ:  this->sweep_weno3_xz(frozen, threadNo); break;
            }
            change = 0;
            for (size_t n = 0; n < this->nodes.size(); ++n) {
                T1 dt    = std::abs(times[n] - this->nodes[n].getTT(threadNo));
                change  += dt;
                times[n] = this->nodes[n].getTT(threadNo);
            }
            niter++;
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

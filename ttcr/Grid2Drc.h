//
//  Grid2Drc.h
//  ttcr
//
//  Created by Bernard Giroux on 15-12-23.
//  Copyright © 2015 Bernard Giroux. All rights reserved.
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
 * @file Grid2Drc.h
 * @brief Base class for 2-D rectilinear grids with cell-based slowness.
 *
 * Declares ttcr::Grid2Drc, the mid-level base shared by every 2-D rectilinear
 * solver whose slowness is constant within each cell: Grid2Drcsp
 * (shortest path), Grid2Drcfs (fast sweeping), Grid2Drcdsp (dynamic shortest
 * path) and Grid2Drcfs_OpenCL.
 *
 * It sits between the abstract ttcr::Grid2D interface and those solvers, and
 * owns everything that does not depend on which eikonal method is used: the
 * grid geometry, the node and cell containers, point location, slowness
 * accessors, traveltime interpolation and raypath tracing. Each derived class
 * adds only its @c raytrace implementation and whatever extra nodes that method
 * needs.
 *
 * @section g2drc_naming Naming
 * In the @c Grid2Drc family, @c r means **rectilinear** (a regular
 * grid with constant spacing @ref ttcr::Grid2Drc::dx and
 * @ref ttcr::Grid2Drc::dz, as opposed to an unstructured triangular mesh) and
 * @c c means slowness held **per cell** (as opposed to @c n, per node — see
 * Grid2Drn.h). The trailing @c sp / @c fs / @c dsp identifies the solver.
 *
 * @section g2drc_numbering Cell and node numbering
 * Cells are numbered **column-wise, z varying fastest**:
 * @f[ \mathrm{cellNo} = i\,n_{cz} + j @f]
 * with @f$i@f$ the x index in @f$[0, n_{cx})@f$ and @f$j@f$ the z index in
 * @f$[0, n_{cz})@f$. @ref ttcr::Grid2Drc::getCellNo maps a point to a cell and
 * @ref ttcr::Grid2Drc::getCellIJ inverts the numbering.
 *
 * Nodes are stored in one vector. The first @f$(n_{cx}+1)(n_{cz}+1)@f$ entries
 * are the **primary** nodes — the cell corners — and any **secondary** nodes a
 * solver adds along the edges follow them. This layout is why
 * @ref ttcr::Grid2Drc::getNumberOfNodes takes a flag, and why
 * @ref ttcr::Grid2Drc::dump_secondary can dump the tail of the vector.
 *
 * @section g2drc_aniso Anisotropy setters
 * The @c setXi, @c setTiltAngle, @c setVp0 … methods forward to the @c CELL
 * policy. Whether they do anything depends entirely on which policy class was
 * supplied: ttcr::Cell ignores them, while the anisotropic policies in Cell.h
 * store the values. Calling one on a grid whose cell type does not support it
 * throws — see Cell.h.
 *
 * @sa Grid2D.h, Grid2Drn.h, Cell.h, Grid2Drcsp.h, Grid2Drcfs.h, Grid2Drcdsp.h
 */

#ifndef ttcr_Grid2Drc_h
#define ttcr_Grid2Drc_h

#include <algorithm>
#include <cstring>
#include <exception>
#include <iostream>
#include <fstream>
#include <memory>
#include <mutex>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifdef VTK
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkXMLRectilinearGridWriter.h"
#endif

#include <boost/math/special_functions/sign.hpp>

#include "Grid2D.h"
#include "NodeKDTree2D.h"

namespace ttcr {

    /**
     * @brief 2-D rectilinear grid holding one slowness value per cell.
     *
     * @tparam T1   floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2   integer type of node and cell indices, normally @c uint32_t.
     * @tparam S    point type: @ref sxz for a planar grid, @ref sxyz for the
     *              draped 2-D grids embedded in 3-D (the @c 2Ds programs).
     * @tparam NODE node type, e.g. ttcr::Node2Dc or ttcr::Node2Dcsp — which one
     *              depends on the bookkeeping the derived solver needs.
     * @tparam CELL cell policy supplying the slowness model, from Cell.h. This
     *              is what makes the class isotropic or anisotropic: it provides
     *              @c getSlowness and @c computeDt, and the derived solvers call
     *              those rather than assuming a velocity law.
     *
     * @note Abstract in practice: it implements everything except @c raytrace,
     *       which each derived solver supplies.
     */
    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    class Grid2Drc : public Grid2D<T1,T2,S> {
    public:
        /**
         * @brief Build the grid geometry and allocate its nodes and cells.
         *
         * @param nx   number of cells along x.
         * @param nz   number of cells along z.
         * @param ddx  cell size along x.
         * @param ddz  cell size along z.
         * @param minx x coordinate of the grid origin.
         * @param minz z coordinate of the grid origin.
         * @param ttrp if true, receiver traveltimes are recomputed by
         *             integrating slowness along the traced raypath rather than
         *             interpolated from the nodal values
         *             (ttcr::input_parameters::tt_from_rp).
         * @param nt   number of threads; sizes the per-thread traveltime array
         *             of every node.
         *
         * @post The grid holds @f$(n_x+1)(n_z+1)@f$ primary nodes and
         *       @f$n_x n_z@f$ cells, and spans @f$[minx,\;minx+n_x\,ddx]@f$ by
         *       @f$[minz,\;minz+n_z\,ddz]@f$. Slowness is **not** set — call
         *       @ref setSlowness before raytracing.
         * @note Derived classes append their secondary nodes to @ref nodes
         *       afterwards; see @ref g2drc_numbering.
         */
        Grid2Drc(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
                 const T1 minx, const T1 minz, const bool ttrp,
                 const size_t nt=1) :
        Grid2D<T1,T2,S>(nx*nz, ttrp, nt),
        dx(ddx), dz(ddz), xmin(minx), zmin(minz),
        xmax(minx+nx*ddx), zmax(minz+nz*ddz),
        ncx(nx), ncz(nz),
        nodes(std::vector<NODE>( static_cast<size_t>(ncx+1) * (ncz+1), NODE(nt) )),
        cells(static_cast<size_t>(ncx)*ncz)
        {
        }

        /// Destructor.
        virtual ~Grid2Drc() {
        }

        /**
         * @brief Copy the cell slowness values out of the grid.
         * @param[out] slowness resized to @f$n_{cx}n_{cz}@f$ and filled, in the
         *                      cell numbering of @ref g2drc_numbering.
         */
        void getSlowness(std::vector<T1>& slowness) const {
            if (slowness.size() != static_cast<size_t>(ncx*ncz)) {
                slowness.resize(static_cast<size_t>(ncx*ncz));
            }
            for (size_t n=0; n<slowness.size(); ++n) {
                slowness[n] = cells.getSlowness(n);
            }
        }

        /**
         * @brief Set the slowness of every cell.
         * @param s one slowness per cell, in the numbering of
         *          @ref g2drc_numbering.
         * @throws std::length_error if @p s does not hold exactly
         *         @f$n_{cx}n_{cz}@f$ values, and whatever else the @c CELL
         *         policy raises.
         * @note The @c try / @c catch around the call rethrows unchanged and so
         *       has no effect; the exception propagates either way. The same
         *       idiom is repeated in every setter below.
         */
        void setSlowness(const std::vector<T1>& s) {
            try {
                cells.setSlowness( s );
            } catch (std::exception& e) {
                throw;
            }
        }
        /**
         * @brief Set the elliptical anisotropy ratio @f$\xi@f$ of every cell.
         * @param x one value per cell.
         * @throws std::exception if the @c CELL policy does not model
         *         @f$\xi@f$. @sa @ref g2drc_aniso
         */
        void setXi(const std::vector<T1>& x) {
            try {
                cells.setXi( x );
            } catch (std::exception& e) {
                throw;
            }
        }
        /**
         * @brief Set the symmetry-axis tilt angle of every cell.
         * @param t one angle per cell, in radians.
         * @throws std::exception if the @c CELL policy is not tilted.
         *         @sa @ref g2drc_aniso
         */
        void setTiltAngle(const std::vector<T1>& t) {
            try {
                cells.setTiltAngle( t );
            } catch (std::exception& e) {
                throw;
            }
        }
        /**
         * @brief Set the P-wave velocity along the symmetry axis, @f$V_{P0}@f$.
         * @param s one value per cell.
         * @throws std::exception if the @c CELL policy is not VTI.
         *         @sa @ref g2drc_aniso
         */
        void setVp0(const std::vector<T1>& s) {
            try {
                cells.setVp0(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        /**
         * @brief Set the S-wave velocity along the symmetry axis, @f$V_{S0}@f$.
         * @param s one value per cell.
         * @throws std::exception if the @c CELL policy is not VTI.
         *         @sa @ref g2drc_aniso
         */
        void setVs0(const std::vector<T1>& s) {
            try {
                cells.setVs0(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        /**
         * @brief Set Thomsen's @f$\delta@f$ for every cell.
         * @param s one value per cell.
         * @throws std::exception if the @c CELL policy is not VTI.
         *         @sa @ref g2drc_aniso
         */
        void setDelta(const std::vector<T1>& s) {
            try {
                cells.setDelta(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        /**
         * @brief Set Thomsen's @f$\epsilon@f$ for every cell.
         * @param s one value per cell.
         * @throws std::exception if the @c CELL policy is not VTI.
         *         @sa @ref g2drc_aniso
         */
        void setEpsilon(const std::vector<T1>& s) {
            try {
                cells.setEpsilon(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        /**
         * @brief Set Thomsen's @f$\gamma@f$ for every cell.
         * @param s one value per cell.
         * @throws std::exception if the @c CELL policy does not model SH
         *         anisotropy. @sa @ref g2drc_aniso
         */
        void setGamma(const std::vector<T1>& s) {
            try {
                cells.setGamma(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        /**
         * @brief Set the weakly-anelliptical coefficient @f$s_2@f$ of every cell.
         * @param s one value per cell.
         * @throws std::exception if the @c CELL policy is not weakly
         *         anelliptical. @sa @ref g2drc_aniso
         */
        void setS2(const std::vector<T1>& s) {
            try {
                cells.setS2(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        /**
         * @brief Set the weakly-anelliptical coefficient @f$s_4@f$ of every cell.
         * @param s one value per cell.
         * @throws std::exception if the @c CELL policy is not weakly
         *         anelliptical. @sa @ref g2drc_aniso
         */
        void setS4(const std::vector<T1>& s) {
            try {
                cells.setS4(s);
            } catch (std::exception& e) {
                throw;
            }
        }

        /**
         * @brief Number of nodes in the grid.
         * @param primary true to count only the primary (cell-corner) nodes,
         *                false (the default) to count every node including the
         *                secondary ones a solver has added.
         * @return The requested count. @sa @ref g2drc_numbering
         */
        size_t getNumberOfNodes(const bool primary=false) const {
            if ( primary ) {
                return static_cast<size_t>(ncx+1) * (ncz+1);
            } else {
                return nodes.size();
            }
        }
        /// @return Number of cells, @f$n_{cx}n_{cz}@f$.
        size_t getNumberOfCells() const { return static_cast<size_t>(ncx)*ncz; }

        /**
         * @brief Collect the traveltimes computed at the primary nodes.
         * @param[out] tt       resized to the primary-node count and filled in
         *                      node order.
         * @param[in]  threadNo thread whose solution to read; must match the
         *                      thread that ran the corresponding @c raytrace.
         * @note Secondary nodes are skipped, so the result is a plain
         *       @f$(n_{cx}+1)\times(n_{cz}+1)@f$ field regardless of how many
         *       extra nodes the solver used internally.
         */
        void getTT(std::vector<T1>& tt, const size_t threadNo=0) const {
            size_t nPrimary = static_cast<size_t>(ncx+1) * (ncz+1);
            tt.resize(nPrimary);
            size_t n = 0;
            for ( size_t nn=0; nn<nodes.size(); ++nn ) {
                if ( nodes[nn].isPrimary() ) {
                    tt[n++] = nodes[nn].getTT(threadNo);
                }
            }
        }

        /**
         * @brief Locate the cell containing a point.
         * @param pt point to locate.
         * @return Cell number, per @ref g2drc_numbering.
         * @note A point on the far boundary (within @ref small of @c xmax or
         *       @c zmax) is pulled back to the centre of the last cell, so the
         *       edge of the grid belongs to the cell inside it rather than
         *       indexing one past the end.
         * @warning A point outside the grid is **not** detected: the arithmetic
         *          simply yields an out-of-range cell number. Callers screen
         *          points with @ref checkPts first.
         */
        T2 getCellNo(const S& pt) const {
            T1 x = xmax-pt.x < small ? xmax-.5*dx : pt.x;
            T1 z = zmax-pt.z < small ? zmax-.5*dz : pt.z;
            T2 nx = static_cast<T2>( small + (x-xmin)/dx );
            T2 nz = static_cast<T2>( small + (z-zmin)/dz );
            return nx*ncz + nz;
        }

        /**
         * @brief Write the traveltime field to a file.
         * @param fname  output filename; an extension is appended per format.
         * @param all    if nonzero, include the secondary nodes as well as the
         *               primary ones.
         * @param nt     thread whose solution to write.
         * @param format 1 for plain text, 2 for VTK, 3 for a raw binary dump.
         * @note The VTK format requires the library to have been built with
         *       @c VTK defined.
         */
        void saveTT(const std::string &fname, const int all, const size_t nt=0,
                    const int format=1) const;

        /**
         * @brief Write the traveltime gradient field to a file.
         * @note **Not implemented** for this family — the body is empty, so the
         *       call silently does nothing. It exists only to satisfy the
         *       ttcr::Grid2D interface.
         */
        void saveTTgrad(const std::string &, const size_t nt=0,
                        const bool vtkFormat=0) const {}

        const T1 getXmin() const { return xmin; }  ///< @return x coordinate of the grid origin.
        const T1 getZmin() const { return zmin; }  ///< @return z coordinate of the grid origin.
        const T1 getDx() const { return dx; }      ///< @return Cell size along x.
        const T1 getDz() const { return dz; }      ///< @return Cell size along z.
        const T2 getNcx() const { return ncx; }    ///< @return Number of cells along x.
        const T2 getNcz() const { return ncz; }    ///< @return Number of cells along z.

        /**
         * @brief Write the coordinates of the secondary nodes to a stream.
         * @param os destination stream, one @c "x z" pair per line.
         * @note Dumps the tail of the node vector past the primary nodes, so it
         *       writes nothing on a grid whose solver adds none. Backs the
         *       @c -s command-line switch
         *       (ttcr::input_parameters::dump_secondary).
         */
        void dump_secondary(std::ofstream& os) const {
            size_t nPrimary = static_cast<size_t>(ncx+1) * (ncz+1);
            for ( size_t n=nPrimary; n<nodes.size(); ++n ) {
                os << nodes[n].getX() << ' ' << nodes[n].getZ() << '\n';
            }
        }
        
        /**
         * @brief Slowness at an arbitrary point.
         * @param pt point to sample.
         * @return Slowness of the cell containing @p pt.
         * @note Piecewise constant, not interpolated — every point inside a cell
         *       returns the same value. That is the defining property of the
         *       @c c family; the @c n grids interpolate between nodal values
         *       instead.
         */
        T1 computeSlowness(const S& pt) const {
            return cells.getSlowness(getCellNo(pt));
        }

    protected:
        T1 dx;           ///< cell size in x
        T1 dz;           ///< cell size in z
        T1 xmin;         ///< x origin of the grid
        T1 zmin;         ///< z origin of the grid
        T1 xmax;         ///< x end of the grid
        T1 zmax;         ///< z end of the grid
        T2 ncx;          ///< number of cells in x
        T2 ncz;          ///< number of cells in z

        /// Primary nodes first, then any secondary nodes the derived solver
        /// added. @c mutable because the traveltimes stored in them are updated
        /// by @c const raytracing methods. @sa @ref g2drc_numbering
        mutable std::vector<NODE> nodes;

        /// kd-tree over all nodes (primary + secondary), replacing the
        /// brute-force scan that checks whether a point coincides with a node.
        /// Built lazily on first use by @ref getNodeNo.
        // kd-tree over all nodes (primary + secondary), to replace the
        // brute-force scan checking whether a point coincides with a node.
        // Built lazily; node coordinates are fixed after construction and
        // queries are read-only, hence thread-safe.
        mutable std::unique_ptr<NodeKDTree2D<T1,T2>> kdtree;
        mutable std::once_flag kdtreeFlag;  ///< Guards the one-time @ref kdtree build.

        /**
         * @brief Find the node coinciding with a point, if any.
         * @param pt point to test.
         * @return Index of the node at @p pt, or @c T2's maximum if no node lies
         *         there.
         * @note Builds @ref kdtree on first call. @c std::call_once makes that
         *       safe under concurrency, and since node coordinates are fixed
         *       after construction and queries are read-only, later calls need
         *       no synchronisation.
         * @note "Coincides" means within @ref small, via the node's
         *       @c operator==; the kd-tree only narrows the search to the
         *       nearest candidate.
         */
        T2 getNodeNo(const S& pt) const {
            std::call_once(kdtreeFlag, [this]() {
                kdtree.reset(new NodeKDTree2D<T1,T2>(nodes, nodes.size()));
            });
            T2 nn = kdtree->findNearest(pt.x, pt.z);
            return nodes[nn] == pt ? nn : std::numeric_limits<T2>::max();
        }

        /// Slowness model, one value per cell, in the column-wise (z fastest)
        /// order of @ref g2drc_numbering. Its type is the @c CELL policy, so
        /// this member also carries any anisotropy parameters.
        CELL cells;   // column-wise (z axis) slowness vector of the cells

        /**
         * @brief Verify that every point lies inside the grid.
         * @param pts points to check, typically sources or receivers.
         * @throws std::runtime_error naming the offending point if one falls
         *         outside the grid bounds.
         * @note Called before raytracing, because @ref getCellNo would otherwise
         *       return an out-of-range cell for such a point.
         */
        void checkPts(const std::vector<S>& pts) const;

        /**
         * @brief Test whether a point lies inside a polygon.
         * @param p    point to test.
         * @param poly polygon vertices, in order.
         * @param N    number of vertices.
         * @return True if @p p is inside.
         */
        bool inPolygon(const S& p, const S poly[], const size_t N) const;

        /**
         * @brief Split a cell number into its x and z indices.
         * @param[in]  cellNo cell number.
         * @param[out] ind    @c i receives the x index, @c j the z index.
         * @sa @ref g2drc_numbering — this inverts @ref getCellNo.
         */
        void getCellIJ(const T2 cellNo, sij<T2>& ind) const {
            ind.i = cellNo / ncz;
            ind.j = cellNo - ncz * ind.i;
        }

        /**
         * @brief Compute the x and z cell indices of a point.
         * @param[in]  pt point to locate.
         * @param[out] i  x index.
         * @param[out] j  z index.
         * @warning Unchecked: a point outside the grid yields an out-of-range
         *          index, and on an unsigned @c T2 a point below the origin
         *          wraps to a huge value.
         */
        void getIJ(const S& pt, T2& i, T2& j) const {
            i = static_cast<T2>( small + (pt.x-xmin)/dx );
            j = static_cast<T2>( small + (pt.z-zmin)/dz );
        }

        /**
         * @brief Compute the x and z cell indices of a point, as signed values.
         * @param[in]  pt point to locate.
         * @param[out] i  x index.
         * @param[out] j  z index.
         * @note Overload returning @c ptrdiff_t, for callers that step off the
         *       grid and need to detect a negative index rather than wrap.
         */
        void getIJ(const S& pt, ptrdiff_t& i, ptrdiff_t& j) const {
            i = static_cast<ptrdiff_t>( small + (pt.x-xmin)/dx );
            j = static_cast<ptrdiff_t>( small + (pt.z-zmin)/dz );
        }

        /**
         * @brief Traveltime at a receiver.
         * @param Rx       receiver position.
         * @param threadNo thread whose solution to read.
         * @return Traveltime at @p Rx: the nodal value if it sits on a node,
         *         otherwise interpolated by @ref interpolateTraveltime.
         * @pre A @c raytrace call has populated the node traveltimes for
         *      @p threadNo.
         */
        T1 getTraveltime(const S& Rx, const size_t threadNo) const;

    private:
        /**
         * @brief Estimate the traveltime gradient at a point.
         * @param[out] g  gradient vector.
         * @param[in]  pt point at which to evaluate it.
         * @param[in]  nt thread whose traveltime field to differentiate.
         * @note The raypath tracers step down this gradient. @sa Grad.h
         */
        void grad(S &g, const S &pt, const size_t nt) const;

        /**
         * @brief Interpolate the traveltime field at an arbitrary point.
         * @param pt       point to evaluate.
         * @param threadNo thread whose solution to read.
         * @return Interpolated traveltime.
         */
        T1 interpolateTraveltime(const S& pt, const size_t threadNo) const;

        /**
         * @brief Traveltime obtained by integrating slowness along the raypath.
         *
         * Traces the ray from @p Rx back to the nearest source and accumulates
         * @f$\int s\,\mathrm{d}l@f$ along it, instead of reading the
         * interpolated nodal value. Selected by the @c ttrp constructor flag.
         *
         * @param Tx       source positions.
         * @param t0       origin time of each source.
         * @param Rx       receiver position.
         * @param threadNo thread whose solution to use.
         * @return Traveltime at @p Rx.
         * @note Usually more accurate than interpolation, since it follows the
         *       path rather than smoothing across cells, but it costs a full
         *       raypath trace per receiver.
         */
        T1 getTraveltimeFromRaypath(const std::vector<S>& Tx,
                                    const std::vector<T1>& t0,
                                    const S& Rx,
                                    const size_t threadNo) const final;

        /**
         * @brief Trace the raypath from a receiver back to the source.
         * @param[in]  Tx       source positions.
         * @param[in]  t0       origin time of each source.
         * @param[in]  Rx       receiver position.
         * @param[out] r_data   raypath points, receiver first and source last.
         * @param[out] tt       traveltime along the path.
         * @param[in]  threadNo thread whose solution to use.
         * @note Integrates down the traveltime gradient from @p Rx; the step is
         *       floored by ttcr::input_parameters::min_distance_rp.
         */
        void getRaypath(const std::vector<S>& Tx,
                        const std::vector<T1>& t0,
                        const S& Rx,
                        std::vector<S>& r_data,
                        T1 &tt,
                        const size_t threadNo) const final;

        /**
         * @brief Trace the raypath and record the length travelled in each cell.
         * @param[in]  Tx       source positions.
         * @param[in]  t0       origin time of each source.
         * @param[in]  Rx       receiver position.
         * @param[out] r_data   raypath points.
         * @param[out] l_data   per-cell segment lengths, as (cell index, length)
         *                      pairs — one row of the sensitivity matrix.
         * @param[out] tt       traveltime along the path.
         * @param[in]  threadNo thread whose solution to use.
         * @note This is the overload used when
         *       ttcr::input_parameters::saveM is set. @sa siv
         */
        void getRaypath(const std::vector<S>& Tx,
                        const std::vector<T1>& t0,
                        const S& Rx,
                        std::vector<S>& r_data,
                        std::vector<siv<T1>> &l_data,
                        T1 &tt,
                        const size_t threadNo) const final;

        /**
         * @brief Record the per-cell path lengths without returning the geometry.
         * @param[in]  Tx       source positions.
         * @param[in]  t0       origin time of each source.
         * @param[in]  Rx       receiver position.
         * @param[out] l_data   per-cell segment lengths.
         * @param[out] tt       traveltime along the path.
         * @param[in]  threadNo thread whose solution to use.
         * @note Same computation as the overload above, minus the raypath
         *       points — for callers that want only the sensitivity row.
         */
        void getRaypath(const std::vector<S>& Tx,
                        const std::vector<T1>& t0,
                        const S& Rx,
                        std::vector<siv<T1>> &l_data,
                        T1 &tt,
                        const size_t threadNo) const final;

    };

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Drc<T1,T2,S,NODE,CELL>::saveTT(const std::string& fname,
                                             const int all,
                                             const size_t nt,
                                             const int format) const {

        if ( format == 1 ) {
            std::string filename = fname+".dat";
            std::ofstream fout(filename.c_str());
            fout.precision(12);
            for ( T2 n=0; n<nodes.size(); ++n ) {
                fout << nodes[n].getX() << '\t'
                << nodes[n].getZ() << '\t'
                << nodes[n].getTT(nt) << '\n';
            }
            fout.close();
        } else  if ( format == 2 ) {
#ifdef VTK

            std::string filename = fname+".vtr";
            int nn[3] = {static_cast<int>(ncx+1), 1, static_cast<int>(ncz+1)};

            vtkSmartPointer<vtkDoubleArray> xCoords = vtkSmartPointer<vtkDoubleArray>::New();
            for (size_t n=0; n<nn[0]; ++n) {
                xCoords->InsertNextValue( xmin + n*dx );
            }
            vtkSmartPointer<vtkDoubleArray> yCoords = vtkSmartPointer<vtkDoubleArray>::New();
            yCoords->InsertNextValue( 0.0 );
            vtkSmartPointer<vtkDoubleArray> zCoords = vtkSmartPointer<vtkDoubleArray>::New();
            for (size_t n=0; n<nn[2]; ++n) {
                zCoords->InsertNextValue( zmin + n*dz );
            }

            vtkSmartPointer<vtkRectilinearGrid> rgrid = vtkSmartPointer<vtkRectilinearGrid>::New();
            rgrid->SetDimensions( nn );
            rgrid->SetXCoordinates(xCoords);
            rgrid->SetYCoordinates(yCoords);
            rgrid->SetZCoordinates(zCoords);


            vtkSmartPointer<vtkDoubleArray> newScalars =
            vtkSmartPointer<vtkDoubleArray>::New();

            newScalars->SetName("Travel time");
            newScalars->SetNumberOfComponents(1);
            newScalars->SetNumberOfTuples( rgrid->GetNumberOfPoints() );


            for ( size_t n=0; n<nodes.size(); ++n ) {
                if ( nodes[n].isPrimary() ) {
                    double x[3] = {nodes[n].getX(), 0.0, nodes[n].getZ()};
                    vtkIdType id = rgrid->FindPoint(x);
                    newScalars->SetTuple1(id, nodes[n].getTT(nt) );
                }
            }
            rgrid->GetPointData()->SetScalars(newScalars);

            vtkSmartPointer<vtkXMLRectilinearGridWriter> writer =
            vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();

            writer->SetFileName( filename.c_str() );
            writer->SetInputData( rgrid );
            writer->SetDataModeToBinary();
            writer->Update();
#else
            std::cerr << "VTK not included during compilation.\nNothing saved.\n";
#endif
        } else if ( format == 3 ) {
            std::string filename = fname+".bin";
            std::ofstream fout(filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
            for ( T2 n=0; n<nodes.size(); ++n ) {
                T1 tmp[] = { nodes[n].getX(), nodes[n].getZ(), nodes[n].getTT(nt) };
                fout.write( (char*)tmp, 3*sizeof(T1) );
            }
            fout.close();
        } else {
            throw std::runtime_error("Unsupported format for saving traveltimes");
        }
    }


    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Drc<T1,T2,S,NODE,CELL>::checkPts(const std::vector<S>& pts) const {
        for (size_t n=0; n<pts.size(); ++n) {
            if ( pts[n].x < xmin || pts[n].x > xmax ||
                pts[n].z < zmin || pts[n].z > zmax ) {
                std::ostringstream msg;
                msg << "Error: Point (" << pts[n].x << ", "<< pts[n] .z << ") outside grid.";
                throw std::runtime_error(msg.str());
            }
        }
    }


    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    bool Grid2Drc<T1,T2,S,NODE,CELL>::inPolygon(const S& p, const S poly[],
                                                const size_t N) const {
        bool c = false;
        for (size_t i = 0, j = N-1; i < N; j = i++) {
            if ((((poly[i].z <= p.z) && (p.z < poly[j].z)) ||
                 ((poly[j].z <= p.z) && (p.z < poly[i].z))) &&
                (p.x < (poly[j].x - poly[i].x) * (p.z - poly[i].z) / (poly[j].z - poly[i].z) + poly[i].x))
                c = !c;
        }
        return c;
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    T1 Grid2Drc<T1,T2,S,NODE,CELL>::interpolateTraveltime(const S& pt,
                                                          const size_t nt) const {

        const size_t nnz = ncz+1;

        // bilinear interpolation if not on node

        T1 tt;
        T2 i, k;

        getIJ(pt, i, k);
        if ( std::abs(pt.x - (xmin+i*dx))<small &&
            std::abs(pt.z - (zmin+k*dz))<small ) {
            // on node
            return nodes[i*nnz+k].getTT(nt);
        } else if ( std::abs(pt.x - (xmin+i*dx))<small ) {
            // on edge
            T1 t1 = nodes[i*nnz+k  ].getTT(nt);
            T1 t2 = nodes[i*nnz+k+1].getTT(nt);

            T1 w1 = (zmin+(k+1)*dz - pt.z)/dz;
            T1 w2 = (pt.z - (zmin+k*dz))/dz;

            tt = t1*w1 + t2*w2;
        } else if ( std::abs(pt.z - (zmin+k*dz))<small ) {
            // on edge
            T1 t1 = nodes[    i*nnz+k].getTT(nt);
            T1 t2 = nodes[(i+1)*nnz+k].getTT(nt);

            T1 w1 = (xmin+(i+1)*dx - pt.x)/dx;
            T1 w2 = (pt.x - (xmin+i*dx))/dx;

            tt = t1*w1 + t2*w2;
        } else {
            T1 t1 = nodes[    i*nnz+k  ].getTT(nt);
            T1 t2 = nodes[    i*nnz+k+1].getTT(nt);
            T1 t3 = nodes[(i+1)*nnz+k  ].getTT(nt);
            T1 t4 = nodes[(i+1)*nnz+k+1].getTT(nt);

            T1 w1 = (zmin+(k+1)*dz - pt.z)/dz;
            T1 w2 = (pt.z - (zmin+k*dz))/dz;

            t1 = t1*w1 + t2*w2;
            t2 = t3*w1 + t4*w2;

            w1 = (xmin+(i+1)*dx - pt.x)/dx;
            w2 = (pt.x - (xmin+i*dx))/dx;

            tt = t1*w1 + t2*w2;
        }
        return tt;
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Drc<T1,T2,S,NODE,CELL>::grad(S &g, const S &pt,
                                           const size_t nt) const {

        // compute travel time gradient at point pt

        T1 p1 = pt.x - dx/2.0;
        if (p1 < xmin) {
            p1 = xmin;
        }
        T1 p2 = p1 + dx;
        if (p2 > xmax) {
            p2 = xmax;
            p1 = xmax - dx;
        }
        g.x = (interpolateTraveltime({p2, pt.z}, nt) - interpolateTraveltime({p1, pt.z}, nt)) / dx;

        p1 = pt.z - dz/2.0;
        if (p1 < zmin) {
            p1 = zmin;
        }
        p2 = p1 + dz;
        if (p2 > zmax) {
            p2 = zmax;
            p1 = zmax - dz;
        }
        g.z = (interpolateTraveltime({pt.x, p2}, nt) - interpolateTraveltime({pt.x, p1}, nt)) / dz;
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    T1 Grid2Drc<T1,T2,S,NODE,CELL>::getTraveltime(const S& Rx, const size_t threadNo) const {

        T2 nn = getNodeNo( Rx );
        if ( nn != std::numeric_limits<T2>::max() ) {
            return nodes[nn].getTT(threadNo);
        }

        T2 cellNo = getCellNo( Rx );
        T2 neibNo = this->neighbors[cellNo][0];
        T1 dt = cells.computeDt(nodes[neibNo], Rx, cellNo);

        T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
        for ( size_t k=1; k< this->neighbors[cellNo].size(); ++k ) {
            neibNo = this->neighbors[cellNo][k];
            dt = cells.computeDt(nodes[neibNo], Rx, cellNo);
            if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
                traveltime =  nodes[neibNo].getTT(threadNo)+dt;
            }
        }
        return traveltime;
    }


    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    T1 Grid2Drc<T1,T2,S,NODE,CELL>::getTraveltimeFromRaypath(const std::vector<S>& Tx,
                                                             const std::vector<T1>& t0,
                                                             const S& Rx,
                                                             const size_t threadNo) const {

        T1 tt = 0.0;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return tt;
            }
        }

        S prev_pt( Rx );
        S curr_pt( Rx );
        // distance between opposite nodes of a voxel
        const T1 maxDist = sqrt( dx*dx + dz*dz );
        S g;

        bool reachedTx = false;
        while ( reachedTx == false ) {

            grad(g, curr_pt, threadNo);
            g *= -1.0;

            ptrdiff_t i, k;
            getIJ(curr_pt, i, k);

            // planes we will intersect
            T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
            T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

            if ( std::abs(xp-curr_pt.x)<small) {
                xp += dx*boost::math::sign(g.x);
            }
            if ( std::abs(zp-curr_pt.z)<small) {
                zp += dz*boost::math::sign(g.z);
            }

            // dist to planes
            T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
            T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

            if ( tx<tz ) { // closer to xp
                curr_pt += tx*g;
                curr_pt.x = xp;     // make sure we don't accumulate rounding errors
            } else {
                curr_pt += tz*g;
                curr_pt.z = zp;
            }

            if ( curr_pt.x < xmin || curr_pt.x > xmax ||
                curr_pt.z < zmin || curr_pt.z > zmax ) {
                //  we are going oustide the grid!
                // make gardient point along outside face
                if ( abs(g.x) > abs(g.z) ) {
                    g.x = boost::math::sign(g.x);
                    g.z = 0.0;
                } else {
                    g.x = 0.0;
                    g.z = boost::math::sign(g.z);
                }

                // put back previous coordinates
                curr_pt = prev_pt;

                // planes we will intersect
                T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
                T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

                if ( std::abs(xp-curr_pt.x)<small) {
                    xp += dx*boost::math::sign(g.x);
                }
                if ( std::abs(zp-curr_pt.z)<small) {
                    zp += dz*boost::math::sign(g.z);
                }

                // dist to planes
                T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
                T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

                if ( tx<tz ) { // closer to xp
                    curr_pt += tx*g;
                    curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                } else {
                    curr_pt += tz*g;
                    curr_pt.z = zp;
                }

                if ( curr_pt.x < xmin || curr_pt.x > xmax ||
                    curr_pt.z < zmin || curr_pt.z > zmax ) {
                    std::ostringstream msg;
                    msg << "Error while computing raypaths: going outside grid \n\
                    Rx: " << Rx << "\n\
                    Tx: " << Tx[0] << "\n";
                    for ( size_t ns=1; ns<Tx.size(); ++ns ) {
                        msg << "\
                        " << Tx[ns] << "\n";
                    }
                    throw std::runtime_error(msg.str());
                }
            }

            S mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
            T2 cellNo = getCellNo(mid_pt);
            tt += cells.computeDt(prev_pt, curr_pt, cellNo);
            prev_pt = curr_pt;

            // are we close enough to one the Tx nodes ?
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                T1 dist = curr_pt.getDistance( Tx[ns] );
                if ( dist < maxDist ) {

                    g = Tx[ns] - curr_pt;

                    T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
                    T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

                    if ( std::abs(xp-curr_pt.x)<small) {
                        xp += dx*boost::math::sign(g.x);
                    }
                    if ( std::abs(zp-curr_pt.z)<small) {
                        zp += dz*boost::math::sign(g.z);
                    }

                    // dist to planes
                    T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
                    T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

                    if ( tx<tz ) { // closer to xp
                        curr_pt += tx*g;
                        curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                    } else {
                        curr_pt += tz*g;
                        curr_pt.z = zp;
                    }

                    if ( curr_pt.getDistance(prev_pt) > dist ||  // we do not intersect
                        curr_pt == Tx[ns] ) {  // we have arrived
                        cellNo = getCellNo(Tx[ns]);
                        tt += cells.computeDt(Tx[ns], prev_pt, cellNo);
                    } else {
                        // to intersection
                        mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
                        cellNo = getCellNo(mid_pt);
                        tt += cells.computeDt(prev_pt, curr_pt, cellNo);
                        // to Tx
                        cellNo = getCellNo(Tx[ns]);
                        tt += cells.computeDt(Tx[ns], curr_pt, cellNo);
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
        return tt;
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Drc<T1,T2,S,NODE,CELL>::getRaypath(const std::vector<S>& Tx,
                                                 const std::vector<T1>& t0,
                                                 const S& Rx,
                                                 std::vector<S>& r_data,
                                                 T1 &tt,
                                                 const size_t threadNo) const {
        tt = 0.0;
        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        S curr_pt( Rx );
        // distance between opposite nodes of a voxel
        const T1 maxDist = sqrt( dx*dx + dz*dz );
        S g;

        bool reachedTx = false;
        while ( reachedTx == false ) {

            grad(g, curr_pt, threadNo);
            g *= -1.0;

            ptrdiff_t i, k;
            getIJ(curr_pt, i, k);

            // planes we will intersect
            T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
            T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

            if ( std::abs(xp-curr_pt.x)<small) {
                xp += dx*boost::math::sign(g.x);
            }
            if ( std::abs(zp-curr_pt.z)<small) {
                zp += dz*boost::math::sign(g.z);
            }

            // dist to planes
            T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
            T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

            if ( tx<tz ) { // closer to xp
                curr_pt += tx*g;
                curr_pt.x = xp;     // make sure we don't accumulate rounding errors
            } else {
                curr_pt += tz*g;
                curr_pt.z = zp;
            }

            if ( curr_pt.x < xmin || curr_pt.x > xmax ||
                curr_pt.z < zmin || curr_pt.z > zmax ) {
                //  we are going oustide the grid!
                // make gardient point along outside face
                if ( abs(g.x) > abs(g.z) ) {
                    g.x = boost::math::sign(g.x);
                    g.z = 0.0;
                } else {
                    g.x = 0.0;
                    g.z = boost::math::sign(g.z);
                }

                // put back previous coordinates
                curr_pt = r_data.back();

                // planes we will intersect
                T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
                T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

                if ( std::abs(xp-curr_pt.x)<small) {
                    xp += dx*boost::math::sign(g.x);
                }
                if ( std::abs(zp-curr_pt.z)<small) {
                    zp += dz*boost::math::sign(g.z);
                }

                // dist to planes
                T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
                T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

                if ( tx<tz ) { // closer to xp
                    curr_pt += tx*g;
                    curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                } else {
                    curr_pt += tz*g;
                    curr_pt.z = zp;
                }

                if ( curr_pt.x < xmin || curr_pt.x > xmax ||
                    curr_pt.z < zmin || curr_pt.z > zmax ) {
                    std::ostringstream msg;
                    msg << "Error while computing raypaths: going outside grid \n\
                    Rx: " << Rx << "\n\
                    Tx: " << Tx[0] << "\n";
                    for ( size_t ns=1; ns<Tx.size(); ++ns ) {
                        msg << "\
                        " << Tx[ns] << "\n";
                    }
                    throw std::runtime_error(msg.str());
                }
            }

            S mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
            T2 cellNo = getCellNo(mid_pt);
            tt += cells.computeDt(r_data.back(), curr_pt, cellNo);
            r_data.push_back( curr_pt );

            // are we close enough to one the Tx nodes ?
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                T1 dist = curr_pt.getDistance( Tx[ns] );
                if ( dist < maxDist ) {

                    g = Tx[ns] - curr_pt;

                    T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
                    T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

                    if ( std::abs(xp-curr_pt.x)<small) {
                        xp += dx*boost::math::sign(g.x);
                    }
                    if ( std::abs(zp-curr_pt.z)<small) {
                        zp += dz*boost::math::sign(g.z);
                    }

                    // dist to planes
                    T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
                    T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

                    if ( tx<tz ) { // closer to xp
                        curr_pt += tx*g;
                        curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                    } else {
                        curr_pt += tz*g;
                        curr_pt.z = zp;
                    }

                    if ( curr_pt.getDistance(r_data.back()) > dist ||  // we do not intersect
                        curr_pt == Tx[ns] ) {  // we have arrived
                        cellNo = getCellNo(Tx[ns]);
                        tt += cells.computeDt(Tx[ns], r_data.back(), cellNo);
                        r_data.push_back( Tx[ns] );
                    } else {
                        // to intersection
                        mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
                        cellNo = getCellNo(mid_pt);
                        tt += cells.computeDt(r_data.back(), curr_pt, cellNo);
                        r_data.push_back( curr_pt );
                        // to Tx
                        cellNo = getCellNo(Tx[ns]);
                        tt += cells.computeDt(Tx[ns], curr_pt, cellNo);
                        r_data.push_back( Tx[ns] );
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Drc<T1,T2,S,NODE,CELL>::getRaypath(const std::vector<S>& Tx,
                                                 const std::vector<T1>& t0,
                                                 const S& Rx,
                                                 std::vector<S>& r_data,
                                                 std::vector<siv<T1>> &l_data,
                                                 T1 &tt,
                                                 const size_t threadNo) const {
        tt = 0.0;
        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        S curr_pt( Rx );
        // distance between opposite nodes of a voxel
        const T1 maxDist = sqrt( dx*dx + dz*dz );
        S g;

        siv<T1> cell;
        bool reachedTx = false;
        while ( reachedTx == false ) {

            grad(g, curr_pt, threadNo);
            g *= -1.0;

            ptrdiff_t i, k;
            getIJ(curr_pt, i, k);

            // planes we will intersect
            T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
            T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

            if ( std::abs(xp-curr_pt.x)<small) {
                xp += dx*boost::math::sign(g.x);
            }
            if ( std::abs(zp-curr_pt.z)<small) {
                zp += dz*boost::math::sign(g.z);
            }

            // dist to planes
            T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
            T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

            if ( tx<tz ) { // closer to xp
                curr_pt += tx*g;
                curr_pt.x = xp;     // make sure we don't accumulate rounding errors
            } else {
                curr_pt += tz*g;
                curr_pt.z = zp;
            }

            if ( curr_pt.x < xmin || curr_pt.x > xmax ||
                curr_pt.z < zmin || curr_pt.z > zmax ) {
                //  we are going oustide the grid!
                std::ostringstream msg;
                msg << "Error while computing raypaths: going outside grid \n\
                Rx: " << Rx << "\n\
                Tx: " << Tx[0] << "\n";
                for ( size_t ns=1; ns<Tx.size(); ++ns ) {
                    msg << "\
                    " << Tx[ns] << "\n";
                }
                throw std::runtime_error(msg.str());
            }

            S mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
            cell.i = getCellNo(mid_pt);
            cell.v = curr_pt.getDistance(r_data.back());
            l_data.push_back(cell);
            tt += cells.computeDt(r_data.back(), curr_pt, cell.i);
            r_data.push_back( curr_pt );

            // are we close enough to one the Tx nodes ?
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                T1 dist = curr_pt.getDistance( Tx[ns] );
                if ( dist < maxDist ) {

                    g = Tx[ns] - curr_pt;

                    T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
                    T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

                    if ( std::abs(xp-curr_pt.x)<small) {
                        xp += dx*boost::math::sign(g.x);
                    }
                    if ( std::abs(zp-curr_pt.z)<small) {
                        zp += dz*boost::math::sign(g.z);
                    }

                    // dist to planes
                    T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
                    T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

                    if ( tx<tz ) { // closer to xp
                        curr_pt += tx*g;
                        curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                    } else {
                        curr_pt += tz*g;
                        curr_pt.z = zp;
                    }

                    if ( curr_pt.getDistance(r_data.back()) > dist ||  // we do not intersect
                        curr_pt == Tx[ns] ) {  // we have arrived
                        cell.i = getCellNo(Tx[ns]);
                        cell.v = Tx[ns].getDistance(r_data.back());
                        l_data.push_back(cell);
                        tt += cells.computeDt(Tx[ns], r_data.back(), cell.i);
                        r_data.push_back( Tx[ns] );
                    } else {
                        // to intersection
                        mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
                        cell.i = getCellNo(mid_pt);
                        cell.v = curr_pt.getDistance(r_data.back());
                        l_data.push_back(cell);
                        tt += cells.computeDt(r_data.back(), curr_pt, cell.i);
                        r_data.push_back( curr_pt );
                        // to Tx
                        size_t itmp = getCellNo(Tx[ns]);
                        if ( cell.i == itmp ) {
                            cell.v += Tx[ns].getDistance(r_data.back());
                        } else {
                            cell.i = itmp;
                            cell.v = Tx[ns].getDistance(r_data.back());
                        }
                        l_data.push_back(cell);
                        tt += cells.computeDt(Tx[ns], curr_pt, cell.i);
                        r_data.push_back( Tx[ns] );
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Drc<T1,T2,S,NODE,CELL>::getRaypath(const std::vector<S>& Tx,
                                                 const std::vector<T1>& t0,
                                                 const S& Rx,
                                                 std::vector<siv<T1>> &l_data,
                                                 T1 &tt,
                                                 const size_t threadNo) const {
        tt = 0.0;
        
        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }
        
        S curr_pt( Rx );
        S prev_pt( Rx );
        // distance between opposite nodes of a voxel
        const T1 maxDist = sqrt( dx*dx + dz*dz );
        S g;
        
        siv<T1> cell;
        bool reachedTx = false;
        while ( reachedTx == false ) {
            
            grad(g, curr_pt, threadNo);
            g *= -1.0;
            
            ptrdiff_t i, k;
            getIJ(curr_pt, i, k);
            
            // planes we will intersect
            T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
            T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));
            
            if ( std::abs(xp-curr_pt.x)<small) {
                xp += dx*boost::math::sign(g.x);
            }
            if ( std::abs(zp-curr_pt.z)<small) {
                zp += dz*boost::math::sign(g.z);
            }
            
            // dist to planes
            T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
            T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();
            
            if ( tx<tz ) { // closer to xp
                curr_pt += tx*g;
                curr_pt.x = xp;     // make sure we don't accumulate rounding errors
            } else {
                curr_pt += tz*g;
                curr_pt.z = zp;
            }
            
            if ( curr_pt.x < xmin || curr_pt.x > xmax ||
                curr_pt.z < zmin || curr_pt.z > zmax ) {
                //  we are going oustide the grid!
                std::ostringstream msg;
                msg << "Error while computing raypaths: going outside grid \n\
            Rx: " << Rx << "\n\
            Tx: " << Tx[0] << "\n";
                for ( size_t ns=1; ns<Tx.size(); ++ns ) {
                    msg << "\
                " << Tx[ns] << "\n";
                }
                throw std::runtime_error(msg.str());
            }
            
            S mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
            cell.i = getCellNo(mid_pt);
            cell.v = curr_pt.getDistance(prev_pt);
            l_data.push_back(cell);
            tt += cells.computeDt(prev_pt, curr_pt, cell.i);
            prev_pt = curr_pt;
            
            // are we close enough to one the Tx nodes ?
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                T1 dist = curr_pt.getDistance( Tx[ns] );
                if ( dist < maxDist ) {
                    
                    g = Tx[ns] - curr_pt;
                    
                    T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
                    T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));
                    
                    if ( std::abs(xp-curr_pt.x)<small) {
                        xp += dx*boost::math::sign(g.x);
                    }
                    if ( std::abs(zp-curr_pt.z)<small) {
                        zp += dz*boost::math::sign(g.z);
                    }
                    
                    // dist to planes
                    T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
                    T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();
                    
                    if ( tx<tz ) { // closer to xp
                        curr_pt += tx*g;
                        curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                    } else {
                        curr_pt += tz*g;
                        curr_pt.z = zp;
                    }
                    
                    if ( curr_pt.getDistance(prev_pt) > dist ||  // we do not intersect
                        curr_pt == Tx[ns] ) {  // we have arrived
                        cell.i = getCellNo(Tx[ns]);
                        cell.v = Tx[ns].getDistance(prev_pt);
                        l_data.push_back(cell);
                        tt += cells.computeDt(Tx[ns], prev_pt, cell.i);
                        
                    } else {
                        // to intersection
                        mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
                        cell.i = getCellNo(mid_pt);
                        cell.v = curr_pt.getDistance(prev_pt);
                        l_data.push_back(cell);
                        tt += cells.computeDt(prev_pt, curr_pt, cell.i);
                        prev_pt = curr_pt;
                        // to Tx
                        size_t itmp = getCellNo(Tx[ns]);
                        if ( cell.i == itmp ) {
                            cell.v += Tx[ns].getDistance(prev_pt);
                        } else {
                            cell.i = itmp;
                            cell.v = Tx[ns].getDistance(prev_pt);
                        }
                        l_data.push_back(cell);
                        tt += cells.computeDt(Tx[ns], curr_pt, cell.i);
                    }
                    
                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
    }

}

#endif /* Grid2Drc_h */

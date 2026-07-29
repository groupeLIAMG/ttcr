//
//  Grid2Drn.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-04-22.
//  Copyright (c) 2014 Bernard Giroux. All rights reserved.
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
 * @file Grid2Drn.h
 * @brief Base class for 2-D rectilinear grids with node-based slowness.
 *
 * Declares ttcr::Grid2Drn, the mid-level base for every 2-D rectilinear solver
 * whose slowness lives at the **nodes** rather than in the cells: Grid2Drnsp
 * (shortest path), Grid2Drnfs (fast sweeping), Grid2Drndsp (dynamic shortest
 * path), Grid2Drnfs_OpenCL — and also Grid2Drcfs and Grid2Drcfs_OpenCL, which
 * despite their @c rc names derive from here (see @ref g2drn_cellhook).
 *
 * It is the counterpart of ttcr::Grid2Drc and provides the same services —
 * geometry, node storage, point location, traveltime interpolation, raypath
 * tracing — plus the fast sweeping stencils, which live here because sweeping
 * needs nodal slowness.
 *
 * @section g2drn_vs_rc Node slowness versus cell slowness
 * The two families differ in where slowness is defined and, consequently, in
 * how a traveltime increment along a segment is computed:
 *
 * - **@c rn (here)** — one value per node. ttcr::Grid2Drn::computeDt averages
 *   the slowness at the two ends and multiplies by the distance,
 *   @f$\Delta t = \tfrac{1}{2}(s_1+s_2)\,\ell@f$ — trapezoidal integration of a
 *   field that varies linearly along the segment.
 * - **@c rc (ttcr::Grid2Drc)** — one value per cell, and the increment comes
 *   from the @c CELL policy, which may be anisotropic.
 *
 * A practical consequence: the @c rn grids have no @c CELL template parameter
 * and are isotropic, since a node carries a single scalar slowness.
 *
 * @section g2drn_cellhook The cell-slowness hook
 * Two virtual hooks, @ref ttcr::Grid2Drn::hasCellSlowness and
 * @ref ttcr::Grid2Drn::getCellSlowness, let a subclass declare that it really
 * does hold a piecewise-constant model. The base returns false and throws
 * respectively; ttcr::Grid2Drcfs overrides both because it accepts a cell model
 * and merely averages it onto the nodes for the solve.
 *
 * The hooks matter when a raypath is integrated: with them the integrator takes
 * the cell value at each segment midpoint, and without them it averages the
 * interpolated slowness at the segment ends. So a @c rc-named subclass sweeps on
 * smoothed nodal values but measures its raypaths against the original
 * piecewise-constant model.
 *
 * @section g2drn_sweeps Fast sweeping stencils
 * Several update stencils are provided and selected by the solver's options:
 * axis-aligned (@c sweep), 45-degree rotated (@c sweep45, and @c sweep_xz for
 * the mixed case) and their third-order WENO variants
 * (@c sweep_weno3, @c sweep_weno3_xz). Rotated stencils are enabled by
 * ttcr::input_parameters::rotated_template and WENO by
 * ttcr::input_parameters::weno3.
 *
 * @sa Grid2D.h, Grid2Drc.h, Grid3Drn.h, Grid2Drnsp.h, Grid2Drnfs.h
 */

#ifndef ttcr_Grid2Drn_h
#define ttcr_Grid2Drn_h

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifdef VTK
#include "vtkDoubleArray.h"
#include "vtkPolyData.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLRectilinearGridWriter.h"
#endif

#include <boost/math/special_functions/sign.hpp>

#include "Grid2D.h"
#include "Interpolator.h"

namespace ttcr {

    /**
     * @brief 2-D rectilinear grid holding one slowness value per node.
     *
     * @tparam T1   floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2   integer type of node and cell indices, normally @c uint32_t.
     * @tparam S    point type: @ref sxz for a planar grid, @ref sxyz for the
     *              draped 2-D grids embedded in 3-D.
     * @tparam NODE node type, e.g. ttcr::Node2Dn or ttcr::Node2Dnsp — whichever
     *              the derived solver needs. It must expose
     *              @c getNodeSlowness.
     *
     * @note No @c CELL parameter, unlike ttcr::Grid2Drc: a node carries a single
     *       scalar slowness, so this family is isotropic.
     * @note Abstract in practice — it implements everything except @c raytrace.
     */
    template<typename T1, typename T2, typename S, typename NODE>
    class Grid2Drn : public Grid2D<T1,T2,S> {
    public:
        /**
         * @brief Build the grid geometry and allocate its nodes.
         *
         * @param nx   number of cells along x.
         * @param nz   number of cells along z.
         * @param ddx  cell size along x.
         * @param ddz  cell size along z.
         * @param minx x coordinate of the grid origin.
         * @param minz z coordinate of the grid origin.
         * @param ttrp recompute receiver traveltimes by integrating slowness
         *             along the traced raypath.
         * @param nt   number of threads; sizes each node's traveltime array.
         *
         * @post The grid holds @f$(n_x+1)(n_z+1)@f$ primary nodes and
         *       spans @f$[minx,\;minx+n_x\,ddx]@f$ by
         *       @f$[minz,\;minz+n_z\,ddz]@f$. No cells are allocated — this
         *       family stores no per-cell data. Slowness is **not** set.
         * @note Derived classes append secondary nodes afterwards; the node
         *       vector is primary-first.
         */
        Grid2Drn(const T2 nx, const T2 nz, const T1 ddx, const T1 ddz,
                 const T1 minx, const T1 minz, const bool ttrp,
                 const size_t nt=1) :
        Grid2D<T1,T2,S>(nx*nz, ttrp, nt),
        dx(ddx), dz(ddz), xmin(minx), zmin(minz),
        xmax(minx+nx*ddx), zmax(minz+nz*ddz),
        ncx(nx), ncz(nz),
        nodes(std::vector<NODE>( static_cast<size_t>(ncx+1) * (ncz+1), NODE(nt) ))
        {
        }

        /// Destructor.
        virtual ~Grid2Drn() {
        }

        /**
         * @brief Set the slowness at every node.
         * @param s one slowness per node — **all** nodes, secondary ones
         *          included, so its size must equal
         *          @ref getNumberOfNodes() with the default argument.
         * @throws std::length_error if the sizes do not match.
         * @note Solvers that add secondary nodes override this with a version
         *       taking one value per *primary* node and interpolating onto the
         *       rest — see ttcr::Grid2Drnsp::setSlowness. This base version is
         *       used as-is only where no secondary nodes exist.
         * @warning Not symmetric with @ref getSlowness, which always returns
         *          primary-node values only.
         */
        void setSlowness(const std::vector<T1>& s) {
            if ( nodes.size() != s.size() ) {
                throw std::length_error("Error: slowness vectors of incompatible size.");
            }
            for ( size_t n=0; n<nodes.size(); ++n ) {
                nodes[n].setNodeSlowness( s[n] );
            }
        }

        /**
         * @brief Copy the primary-node slowness values out of the grid.
         * @param[out] slowness resized to @f$(n_{cx}+1)(n_{cz}+1)@f$ and filled
         *                      from the leading primary nodes.
         * @note Returns primary values only, whatever secondary nodes exist.
         */
        void getSlowness(std::vector<T1>& slowness) const {
            if (slowness.size() != static_cast<size_t>(ncx+1) * (ncz+1)) {
                slowness.resize(static_cast<size_t>(ncx+1) * (ncz+1));
            }
            for ( size_t n=0; n<slowness.size(); ++n ) {
                slowness[n] = nodes[n].getNodeSlowness();
            }
        }

        /**
         * @brief Number of nodes in the grid.
         * @param primary true to count only the primary (cell-corner) nodes,
         *                false (the default) to count every node.
         * @return The requested count.
         */
        size_t getNumberOfNodes(const bool primary=false) const {
            if ( primary ) {
                return static_cast<size_t>(ncx+1) * (ncz+1);
            } else {
                return nodes.size();
            }
        }

        /**
         * @brief Number of cells the grid geometry defines.
         * @return The product @f$n_{cx}n_{cz}@f$.
         * @note Geometric only — no per-cell data is stored by this family.
         */
        size_t getNumberOfCells() const { return static_cast<size_t>(ncx)*ncz; }

        /**
         * @brief Collect the traveltimes computed at the primary nodes.
         * @param[out] tt       resized to the primary-node count and filled.
         * @param[in]  threadNo thread whose solution to read.
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
         * @brief Write the traveltime field to a file.
         * @param fname  output filename.
         * @param all    if nonzero, include the secondary nodes.
         * @param nt     thread whose solution to write.
         * @param format 1 for plain text, 2 for VTK, 3 for a raw binary dump.
         */
        void saveTT(const std::string &fname, const int all, const size_t nt=0,
                    const int format=1) const;

        /**
         * @brief Write the traveltime gradient field to a file.
         * @param fname     output filename.
         * @param nt        thread whose solution to differentiate.
         * @param vtkFormat write VTK rather than plain text.
         * @note Unlike ttcr::Grid2Drc::saveTTgrad, which is an empty stub, this
         *       one is implemented — the gradient is well defined here because
         *       slowness is continuous across cell boundaries.
         */
        void saveTTgrad(const std::string &fname, const size_t nt=0,
                        const bool vtkFormat=0) const;

        const T1 getXmin() const { return xmin; }  ///< @return x coordinate of the grid origin.
        const T1 getZmin() const { return zmin; }  ///< @return z coordinate of the grid origin.
        const T1 getDx() const { return dx; }      ///< @return Cell size along x.
        const T1 getDz() const { return dz; }      ///< @return Cell size along z.
        const T2 getNcx() const { return ncx; }    ///< @return Number of cells along x.
        const T2 getNcz() const { return ncz; }    ///< @return Number of cells along z.

        /**
         * @brief Slowness at an arbitrary point.
         * @param pt point to sample.
         * @return Slowness interpolated from the surrounding nodes.
         * @note Continuous, unlike ttcr::Grid2Drc::computeSlowness which is
         *       piecewise constant per cell.
         */
        T1 computeSlowness(const S& pt) const;

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
        /// added. @c mutable because @c const raytracing methods update the
        /// traveltimes stored in them.
        mutable std::vector<NODE> nodes;

        /**
         * @brief Traveltime increment between two nodes.
         * @param source node the ray comes from.
         * @param node   node it reaches.
         * @return The increment @f$\tfrac{1}{2}(s_{source}+s_{node})\,\ell@f$ — the mean of
         *         the endpoint slownesses times the distance.
         * @note Trapezoidal integration, exact when slowness varies linearly
         *       along the segment. @sa @ref g2drn_vs_rc
         */
        T1 computeDt(const NODE& source, const NODE& node) const {
            return (node.getNodeSlowness()+source.getNodeSlowness())/2 * source.getDistance( node );
        }

        /**
         * @brief Traveltime increment from a node to an arbitrary point.
         * @param source node the ray comes from.
         * @param node   point it reaches.
         * @param slo    slowness at that point, which the caller must have
         *               obtained itself — a plain point carries none.
         * @return The increment @f$\tfrac{1}{2}(s_{source}+slo)\,\ell@f$.
         */
        T1 computeDt(const NODE& source, const S& node, T1 slo) const {
            return (slo+source.getNodeSlowness())/2 * source.getDistance( node );
        }

        /**
         * @brief Whether this grid also holds a piecewise-constant cell model.
         * @return False in the base.
         * @note Override point — see @ref g2drn_cellhook. Returning true commits
         *       the class to implementing @ref getCellSlowness.
         */
        virtual const bool hasCellSlowness() const { return false; }
        /**
         * @brief Slowness of one cell, for subclasses that hold a cell model.
         * @param cell_no cell number.
         * @return That cell's slowness.
         * @throws std::runtime_error in the base — reaching it means a subclass
         *         returned true from @ref hasCellSlowness without overriding
         *         this.
         */
        virtual const T1 getCellSlowness(const size_t cell_no) const {
            throw std::runtime_error("Method getCellSlowness should be implemented in subclass");
        }

        /**
         * @brief Verify that every point lies inside the grid.
         * @param pts points to check, typically sources or receivers.
         * @throws std::runtime_error naming the offending point.
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
         * @brief Traveltime gradient at a node, by grid index.
         * @param[out] g  gradient vector.
         * @param[in]  i  x index of the node.
         * @param[in]  j  z index of the node.
         * @param[in]  nt thread whose traveltime field to differentiate.
         */
        void grad(S &g, const size_t i, const size_t j, const size_t nt=0) const;

        /**
         * @brief Traveltime gradient at an arbitrary point.
         * @param[out] g  gradient vector.
         * @param[in]  pt point at which to evaluate it.
         * @param[in]  nt thread whose traveltime field to differentiate.
         * @note The raypath tracers step down this gradient.
         */
        void grad(S &g, const S &pt, const size_t nt=0) const;

        /**
         * @brief Trace the raypath from a receiver back to the source.
         * @param[in]  Tx       source positions.
         * @param[in]  Rx       receiver position.
         * @param[out] r_data   raypath points, receiver first and source last.
         * @param[in]  threadNo thread whose solution to use.
         * @note Integrates down the traveltime gradient. Slowness along the path
         *       comes from the cell model when @ref hasCellSlowness is true, and
         *       from the interpolated nodal field otherwise.
         *       @sa @ref g2drn_cellhook
         */
        void getRaypath(const std::vector<S>& Tx,
                        const S &Rx,
                        std::vector<S> &r_data,
                        const size_t threadNo=0) const;

        /**
         * @brief Locate the cell containing a point.
         * @param pt point to locate.
         * @return Cell number, z varying fastest.
         * @note Used to index a subclass's cell model via
         *       @ref getCellSlowness; this family stores nothing per cell
         *       itself.
         * @warning A point outside the grid is not detected.
         */
        T2 getCellNo(const S& pt) const {
            T1 x = xmax-pt.x < small ? xmax-.5*dx : pt.x;
            T1 z = zmax-pt.z < small ? zmax-.5*dz : pt.z;
            T2 nx = static_cast<T2>( small + (x-xmin)/dx );
            T2 nz = static_cast<T2>( small + (z-zmin)/dz );
            return nx*ncz + nz;
        }

        /**
         * @brief Split a cell number into its x and z indices.
         * @param[in]  cellNo cell number.
         * @param[out] ind    @c i receives the x index, @c j the z index.
         * @note Cells are numbered z-fastest here, as in ttcr::Grid2Drc — and
         *       unlike the 3-D families, which are x-fastest.
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
         * @note For callers that step off the grid and need a negative index
         *       rather than a wrap.
         */
        void getIJ(const S& pt, ptrdiff_t& i, ptrdiff_t& j) const {
            i = static_cast<ptrdiff_t>( small + (pt.x-xmin)/dx );
            j = static_cast<ptrdiff_t>( small + (pt.z-zmin)/dz );
        }

        /**
         * @name Fast sweeping passes
         *
         * One Gauss-Seidel pass over the grid in each of the alternating
         * directions, relaxing every node that is not frozen. The variants
         * differ in the stencil used; see @ref g2drn_sweeps.
         *
         * @param frozen   per-node flag; a frozen node holds a boundary or
         *                 source value and is never updated.
         * @param threadNo thread whose traveltime field to sweep.
         * @{
         */
        /// Axis-aligned first-order stencil.
        void sweep(const std::vector<bool>& frozen,
                   const size_t threadNo) const;
        /// 45-degree rotated stencil, adding the diagonal directions.
        void sweep45(const std::vector<bool>& frozen,
                     const size_t threadNo) const;
        /// Mixed x-z stencil, for the draped grids where the plane is not axis-aligned.
        void sweep_xz(const std::vector<bool>& frozen,
                      const size_t threadNo) const;
        /// Third-order WENO stencil.
        void sweep_weno3(const std::vector<bool>& frozen,
                         const size_t threadNo) const;
        /// Third-order WENO stencil, mixed x-z form.
        void sweep_weno3_xz(const std::vector<bool>& frozen,
                            const size_t threadNo) const;
        /// @}

        /**
         * @name Single-node updates
         *
         * Solve the local eikonal update at one node, given its already-relaxed
         * neighbours. The parameters are the node's x index, its z index and the
         * thread number. Each corresponds to the sweep of the same name.
         * @{
         */
        void update_node(const size_t, const size_t, const size_t=0) const;
        void update_node45(const size_t, const size_t, const size_t=0) const;
        void update_node_xz(const size_t, const size_t, const size_t=0) const;
        void update_node_weno3(const size_t, const size_t, const size_t=0) const;
        void update_node_weno3_xz(const size_t, const size_t, const size_t=0) const;
        /// @}

        /**
         * @brief Seed the fast sweeping solve around the sources.
         * @param[in]  Tx       source positions.
         * @param[in]  t0       origin time of each source.
         * @param[out] frozen   per-node flag; nodes near a source are given
         *                      analytic traveltimes and frozen, so the sweeps
         *                      cannot move them.
         * @param[in]  npts     half-width, in nodes, of the frozen region around
         *                      each source.
         * @param[in]  threadNo thread to initialise.
         * @note Sweeping cannot start from a point source directly — the local
         *       update needs neighbours that already hold valid times, which is
         *       what this frozen halo provides.
         */
        void initFSM(const std::vector<S>& Tx,
                     const std::vector<T1>& t0, std::vector<bool>& frozen,
                     const int npts, const size_t threadNo) const;

        /**
         * @brief Slowness at a point, interpolated from the surrounding nodes.
         * @param Rx point to sample.
         * @return Interpolated slowness.
         */
        T1 getSlowness(const S& Rx) const;

        /**
         * @brief Write the coordinates of the secondary nodes to a stream.
         * @param os destination stream, one @c "x z" pair per line.
         * @note Writes nothing on a grid whose solver adds no secondary nodes,
         *       which is the case for the sweeping solvers.
         */
        void dump_secondary(std::ofstream& os) const {
            size_t nPrimary = static_cast<size_t>(ncx+1) * (ncz+1);
            for ( size_t n=nPrimary; n<nodes.size(); ++n ) {
                os << nodes[n].getX() << ' ' << nodes[n].getZ() << '\n';
            }
        }

    private:
        Grid2Drn() {}
        Grid2Drn(const Grid2Drn<T1,T2,S,NODE>& g) {}
        Grid2Drn<T1,T2,S,NODE>& operator=(const Grid2Drn<T1,T2,S,NODE>& g) = delete;

        T1 getTraveltime(const S& Rx, const size_t threadNo) const;

        T1 getTraveltimeFromRaypath(const std::vector<S>& Tx,
                                    const std::vector<T1>& t0,
                                    const S& Rx,
                                    const size_t threadNo) const final;

        void getRaypath(const std::vector<S>& Tx,
                        const std::vector<T1>& t0,
                        const S& Rx,
                        std::vector<S>& r_data,
                        T1 &tt,
                        const size_t threadNo) const final;

        void getRaypath(const std::vector<S>& Tx,
                        const std::vector<T1>& t0,
                        const S& Rx,
                        std::vector<S>& r_data,
                        std::vector<siv<T1>> &l_data,
                        T1 &tt,
                        const size_t threadNo) const final;

        void getRaypath(const std::vector<S>& Tx,
                        const std::vector<T1>& t0,
                        const S& Rx,
                        std::vector<siv<T1>> &l_data,
                        T1 &tt,
                        const size_t threadNo) const final;
    };

    template<typename T1, typename T2, typename S, typename NODE>
    T1 Grid2Drn<T1,T2,S,NODE>::computeSlowness(const S& pt) const {
        const size_t nnx = ncx+1;
        const size_t nnz = ncz+1;

        // are we on an node or an edge?
        ptrdiff_t onX = -1;
        ptrdiff_t onZ = -1;
        for ( size_t n=0; n<nnx; ++n ) {
            if ( std::abs(pt.x - (xmin+n*dx)) < small2 ) {
                onX = n;
                break;
            }
        }
        for ( size_t n=0; n<nnz; ++n ) {
            if ( std::abs(pt.z - (zmin+n*dz)) < small2 ) {
                onZ = n;
                break;
            }
        }

        if ( onX!=-1 && onZ!=-1 ) {
            return nodes[onX*nnz + onZ].getNodeSlowness();

        } else if ( onX!=-1 ) {
            T2 k = static_cast<T2>( small + (pt.z-zmin)/dz );
            T1 s[2];
            T1 x[3];
            s[0] = nodes[onX*nnz + k  ].getNodeSlowness();
            s[1] = nodes[onX*nnz + k+1].getNodeSlowness();
            x[0] = pt.z;
            x[1] = zmin + k*dz;
            x[2] = zmin + (k+1)*dz;

            return Interpolator<T1>::linear(x, s);

        } else if ( onZ!=-1 ) {
            T2 i = static_cast<T2>( small + (pt.x-xmin)/dx );
            T1 s[2];
            T1 x[3];
            s[0] = nodes[i*nnz     + onZ].getNodeSlowness();
            s[1] = nodes[(i+1)*nnz + onZ].getNodeSlowness();
            x[0] = pt.x;
            x[1] = xmin + i*dx;
            x[2] = xmin + (i+1)*dx;

            return Interpolator<T1>::linear(x, s);

        } else {
            T2 i = static_cast<T2>( small + (pt.x-xmin)/dx );
            T2 k = static_cast<T2>( small + (pt.z-zmin)/dz );
            T1 s[4];
            T1 x[3];
            T1 z[3];

            s[0] = nodes[i*nnz     + k  ].getNodeSlowness();
            s[1] = nodes[i*nnz     + k+1].getNodeSlowness();
            s[2] = nodes[(i+1)*nnz + k  ].getNodeSlowness();
            s[3] = nodes[(i+1)*nnz + k+1].getNodeSlowness();
            x[0] = pt.x;
            z[0] = pt.z;
            x[1] = xmin + i*dx;
            z[1] = zmin + k*dz;
            x[2] = xmin + (i+1)*dx;
            z[2] = zmin + (k+1)*dz;

            return Interpolator<T1>::bilinear(x, z, s);

        }
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::checkPts(const std::vector<S>& pts) const {
        for (size_t n=0; n<pts.size(); ++n) {
            if ( pts[n].x < xmin || pts[n].x > xmax ||
                pts[n].z < zmin || pts[n].z > zmax ) {
                std::ostringstream msg;
                msg << "Error: Point (" << pts[n].x << ", "<< pts[n] .z << ") outside grid.";
                throw std::runtime_error(msg.str());
            }
        }
    }



    template<typename T1, typename T2, typename S, typename NODE>
    bool Grid2Drn<T1,T2,S,NODE>::inPolygon(const S& p, const S poly[], const size_t N) const {
        bool c = false;
        for (size_t i = 0, j = N-1; i < N; j = i++) {
            if ((((poly[i].z <= p.z) && (p.z < poly[j].z)) ||
                 ((poly[j].z <= p.z) && (p.z < poly[i].z))) &&
                (p.x < (poly[j].x - poly[i].x) * (p.z - poly[i].z) / (poly[j].z - poly[i].z) + poly[i].x))
                c = !c;
        }
        return c;
    }

    template<typename T1, typename T2, typename S, typename NODE>
    T1 Grid2Drn<T1,T2,S,NODE>::getTraveltime(const S &pt, const size_t nt) const {

        const size_t nnz = ncz+1;

        // bilinear interpolation if not on node

        T1 tt;
        T2 i, j;

        getIJ(pt, i, j);

        if ( std::abs(pt.x - (xmin+i*dx))<small && std::abs(pt.z - (zmin+j*dz))<small ) {
            // on node
            return nodes[i*nnz+j].getTT(nt);
        } else if ( std::abs(pt.x - (xmin+i*dx))<small ) {

            // on edge
            T1 t1 = nodes[i*nnz+j].getTT(nt);
            T1 t2 = nodes[i*nnz+j+1].getTT(nt);

            T1 w1 = (zmin+(j+1)*dz - pt.z)/dz;
            T1 w2 = (pt.z - (zmin+j*dz))/dz;

            tt = t1*w1 + t2*w2;

        } else if ( std::abs(pt.z - (zmin+j*dz))<small ) {

            // on edge
            T1 t1 = nodes[i*nnz+j].getTT(nt);
            T1 t2 = nodes[(i+1)*nnz+j].getTT(nt);

            T1 w1 = (xmin+(i+1)*dx - pt.x)/dx;
            T1 w2 = (pt.x - (xmin+i*dx))/dx;

            tt = t1*w1 + t2*w2;

        } else {

            T1 t1 = nodes[    i*nnz+j  ].getTT(nt);
            T1 t2 = nodes[(i+1)*nnz+j  ].getTT(nt);
            T1 t3 = nodes[    i*nnz+j+1].getTT(nt);
            T1 t4 = nodes[(i+1)*nnz+j+1].getTT(nt);

            T1 w1 = (xmin+(i+1)*dx - pt.x)/dx;
            T1 w2 = (pt.x - (xmin+i*dx))/dx;

            t1 = t1*w1 + t2*w2;
            t2 = t3*w1 + t4*w2;

            w1 = (zmin+(j+1)*dz - pt.z)/dz;
            w2 = (pt.z - (zmin+j*dz))/dz;

            tt = t1*w1 + t2*w2;
        }

        return tt;
    }



    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::saveTT(const std::string& fname, const int all,
                                        const size_t nt,
                                        const int format) const {
        if ( format == 1 ) {
            std::string filename = fname+".dat";
            std::ofstream fout(filename.c_str());
            fout.precision(12);
            for ( T2 n=0; n<nodes.size(); ++n ) {
                if ( nodes[n].isPrimary() || all==1 ) {

                    fout << nodes[n].getX() << '\t'
                    << nodes[n].getZ() << '\t'
                    << nodes[n].getTT(nt) << '\n';
                }
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
        } else if ( format == 3 ){
            std::string filename = fname+".bin";
            std::ofstream fout(filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
            for ( T2 n=0; n<nodes.size(); ++n ) {
                if ( nodes[n].isPrimary() || all==1 ) {
                    T1 tmp[] = { nodes[n].getX(), nodes[n].getZ(), nodes[n].getTT(nt) };
                    fout.write( (char*)tmp, 3*sizeof(T1) );
                }
            }
            fout.close();
        } else {
            throw std::runtime_error("Unsupported format for saving traveltimes");
        }
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::saveTTgrad(const std::string& fname,
                                            const size_t nt,
                                            const bool vtkFormat) const {

        if (vtkFormat) {
#ifdef VTK

            std::string filename = fname+".vtr";
            int nn[3] = {static_cast<int>(ncx), 1, static_cast<int>(ncz)};

            vtkSmartPointer<vtkDoubleArray> xCoords = vtkSmartPointer<vtkDoubleArray>::New();
            for (size_t n=0; n<nn[0]; ++n) {
                xCoords->InsertNextValue( xmin + (0.5+n)*dx );
            }
            vtkSmartPointer<vtkDoubleArray> yCoords = vtkSmartPointer<vtkDoubleArray>::New();
            yCoords->InsertNextValue( 0.0 );
            vtkSmartPointer<vtkDoubleArray> zCoords = vtkSmartPointer<vtkDoubleArray>::New();
            for (size_t n=0; n<nn[2]; ++n) {
                zCoords->InsertNextValue( zmin + (0.5+n)*dz );
            }

            vtkSmartPointer<vtkRectilinearGrid> rgrid = vtkSmartPointer<vtkRectilinearGrid>::New();
            rgrid->SetDimensions( nn );
            rgrid->SetXCoordinates(xCoords);
            rgrid->SetYCoordinates(yCoords);
            rgrid->SetZCoordinates(zCoords);


            vtkSmartPointer<vtkDoubleArray> grad_tt =
            vtkSmartPointer<vtkDoubleArray>::New();

            grad_tt->SetName("grad tt");
            grad_tt->SetNumberOfComponents(3);
            grad_tt->SetComponentName(0, "x");
            grad_tt->SetComponentName(1, "y");
            grad_tt->SetComponentName(2, "z");
            grad_tt->SetNumberOfTuples( rgrid->GetNumberOfPoints() );


            double x[3];
            x[1] = 0.0;
            for ( size_t i=0; i<ncx; ++i ) {
                for ( size_t j=0; j<ncz; ++j ) {
                    S g;
                    grad(g, i, j, nt);

                    x[0] = xmin + (i+0.5)*dx;
                    x[2] = zmin + (j+0.5)*dz;

                    vtkIdType id = rgrid->FindPoint(x);
                    grad_tt->SetTuple3(id, g.x, 0.0, g.z );
                }
            }
            rgrid->GetPointData()->SetVectors( grad_tt );


            vtkSmartPointer<vtkXMLRectilinearGridWriter> writer =
            vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();

            writer->SetFileName( filename.c_str() );
            writer->SetInputData( rgrid );
            writer->SetDataModeToBinary();
            writer->Update();
#else
            std::cerr << "VTK not included during compilation.\nNothing saved.\n";
#endif
        } else {
            std::string filename = fname+".dat";
            std::ofstream fout(filename.c_str());
            fout.precision(12);
            for ( size_t i=0; i<ncx; ++i ) {
                for ( size_t j=0; j<ncz; ++j ) {
                    S g;
                    grad(g, i, j, nt);

                    T1 x = xmin + (i+0.5)*dx;
                    T1 z = zmin + (j+0.5)*dz;

                    fout << x << ' ' << z << ' ' << g.x << ' ' << g.z << '\n';
                }
            }
            fout.close();
        }
    }




    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::grad(S& g, const size_t i, const size_t j,
                                      const size_t nt) const {

        // compute average gradient for cell (i,j)

        const size_t nnz = ncz+1;

        g.x = 0.5*(( nodes[(i+1)*nnz+j].getTT(nt)+nodes[(i+1)*nnz+j+1].getTT(nt) ) -
                   ( nodes[    i*nnz+j].getTT(nt)+nodes[    i*nnz+j+1].getTT(nt) ))/dx;
        g.z = 0.5*(( nodes[i*nnz+j+1].getTT(nt)+nodes[(i+1)*nnz+j+1].getTT(nt) ) -
                   ( nodes[i*nnz+j  ].getTT(nt)+nodes[(i+1)*nnz+j  ].getTT(nt) ))/dz;
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::grad(S &g, const S &pt,
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
        g.x = (getTraveltime({p2, pt.z}, nt) - getTraveltime({p1, pt.z}, nt)) / dx;

        p1 = pt.z - dz/2.0;
        if (p1 < zmin) {
            p1 = zmin;
        }
        p2 = p1 + dz;
        if (p2 > zmax) {
            p2 = zmax;
            p1 = zmax - dz;
        }
        g.z = (getTraveltime({pt.x, p2}, nt) - getTraveltime({pt.x, p1}, nt)) / dz;
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::getRaypath(const std::vector<S>& Tx,
                                            const S &Rx,
                                            std::vector<S> &r_data,
                                            const size_t threadNo) const {

        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
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

            r_data.push_back( curr_pt );

            // are we close enough to one the Tx nodes ?
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                if ( curr_pt.getDistance( Tx[ns] ) < maxDist ) {
                    r_data.push_back( Tx[ns] );
                    reachedTx = true;
                }
            }

        }
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::sweep(const std::vector<bool>& frozen,
                                       const size_t threadNo) const {


        // sweep first direction
        for ( size_t i=0; i<=ncx; ++i ) {
            for ( size_t j=0; j<=ncz; ++j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node(i, j, threadNo);
                }
            }
        }

        // sweep second direction
        for ( long int i=ncx; i>=0; --i ) {
            for ( size_t j=0; j<=ncz; ++j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node(i, j, threadNo);
                }
            }
        }

        // sweep third direction
        for ( long int i=ncx; i>=0; --i ) {
            for ( long int j=ncz; j>=0; --j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node(i, j, threadNo);
                }
            }
        }

        // sweep fourth direction
        for ( size_t i=0; i<=ncx; ++i ) {
            for ( long int j=ncz; j>=0; --j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node(i, j, threadNo);
                }
            }
        }
    }


    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::sweep45(const std::vector<bool>& frozen,
                                         const size_t threadNo) const {

        // sweep first direction
        for ( size_t i=0; i<=ncx; ++i ) {
            for ( size_t j=0; j<=ncz; ++j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node45(i, j, threadNo);
                }
            }
        }

        // sweep second direction
        for ( long int i=ncx; i>=0; --i ) {
            for ( size_t j=0; j<=ncz; ++j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node45(i, j, threadNo);
                }
            }
        }

        // sweep third direction
        for ( long int i=ncx; i>=0; --i ) {
            for ( long int j=ncz; j>=0; --j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node45(i, j, threadNo);
                }
            }
        }

        // sweep fourth direction
        for ( size_t i=0; i<=ncx; ++i ) {
            for ( long int j=ncz; j>=0; --j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node45(i, j, threadNo);
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::sweep_xz(const std::vector<bool>& frozen,
                                          const size_t threadNo) const {

        // sweep first direction
        for ( size_t i=0; i<=ncx; ++i ) {
            for ( size_t j=0; j<=ncz; ++j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_xz(i, j, threadNo);
                }
            }
        }

        // sweep second direction
        for ( long int i=ncx; i>=0; --i ) {
            for ( size_t j=0; j<=ncz; ++j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_xz(i, j, threadNo);
                }
            }
        }

        // sweep third direction
        for ( long int i=ncx; i>=0; --i ) {
            for ( long int j=ncz; j>=0; --j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_xz(i, j, threadNo);
                }
            }
        }

        // sweep fourth direction
        for ( size_t i=0; i<=ncx; ++i ) {
            for ( long int j=ncz; j>=0; --j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_xz(i, j, threadNo);
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::sweep_weno3(const std::vector<bool>& frozen,
                                             const size_t threadNo) const {

        // sweep first direction
        for ( size_t i=0; i<=ncx; ++i ) {
            for ( size_t j=0; j<=ncz; ++j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_weno3(i, j, threadNo);
                }
            }
        }

        // sweep second direction
        for ( long int i=ncx; i>=0; --i ) {
            for ( size_t j=0; j<=ncz; ++j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_weno3(i, j, threadNo);
                }
            }
        }

        // sweep third direction
        for ( long int i=ncx; i>=0; --i ) {
            for ( long int j=ncz; j>=0; --j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_weno3(i, j, threadNo);
                }
            }
        }

        // sweep fourth direction
        for ( size_t i=0; i<=ncx; ++i ) {
            for ( long int j=ncz; j>=0; --j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_weno3(i, j, threadNo);
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::sweep_weno3_xz(const std::vector<bool>& frozen,
                                                const size_t threadNo) const {

        // sweep first direction
        for ( size_t i=0; i<=ncx; ++i ) {
            for ( size_t j=0; j<=ncz; ++j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_weno3_xz(i, j, threadNo);
                }
            }
        }

        // sweep second direction
        for ( long int i=ncx; i>=0; --i ) {
            for ( size_t j=0; j<=ncz; ++j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_weno3_xz(i, j, threadNo);
                }
            }
        }

        // sweep third direction
        for ( long int i=ncx; i>=0; --i ) {
            for ( long int j=ncz; j>=0; --j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_weno3_xz(i, j, threadNo);
                }
            }
        }

        // sweep fourth direction
        for ( size_t i=0; i<=ncx; ++i ) {
            for ( long int j=ncz; j>=0; --j ) {
                if ( !frozen[ i*(ncz+1)+j ] ) {
                    update_node_weno3_xz(i, j, threadNo);
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::update_node(const size_t i, const size_t j,
                                             const size_t threadNo) const {

        T1 a, b, t;
        if (i==0)
            a = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo);
        else if (i==ncx)
            a = nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
        else {
            a = nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
            t = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo);
            a = a<t ? a : t;
        }

        if (j==0)
            b = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo);
        else if (j==ncz)
            b = nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
        else {
            b = nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
            t = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo);
            b = b<t ? b : t;
        }

        T1 fh = nodes[i*(ncz+1)+j].getNodeSlowness() * dx;

        if ( std::abs(a-b) >= fh )
            t = (a<b ? a : b) + fh;
        else
            t = 0.5*( a+b + sqrt(2.*fh*fh - (a-b)*(a-b) ) );

        if ( t<nodes[i*(ncz+1)+j].getTT(threadNo) )
            nodes[i*(ncz+1)+j].setTT(t,threadNo);

    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::update_node45(const size_t i, const size_t j,
                                               const size_t threadNo) const {
        // stencil rotated pi/4

        T1 a, b, t;
        if (i==0) {
            if (j!=ncz)
                a = nodes[ (i+1)*(ncz+1)+j+1 ].getTT(threadNo);
            else
                a = std::numeric_limits<T1>::max();
        } else if (i==ncx) {
            if (j!=0)
                a = nodes[ (i-1)*(ncz+1)+j-1 ].getTT(threadNo);
            else
                a = std::numeric_limits<T1>::max();
        } else {
            if (j!=ncz)
                a = nodes[ (i+1)*(ncz+1)+j+1 ].getTT(threadNo);
            else
                a = std::numeric_limits<T1>::max();
            if (j!=0)
                t = nodes[ (i-1)*(ncz+1)+j-1 ].getTT(threadNo);
            else
                t = std::numeric_limits<T1>::max();
            a = a<t ? a : t;
        }

        if (i==0) {
            if (j!=0)
                b = nodes[ (i+1)*(ncz+1)+j-1 ].getTT(threadNo);
            else
                b = std::numeric_limits<T1>::max();
        } else if (i==ncx) {
            if (j!=ncz)
                b = nodes[ (i-1)*(ncz+1)+j+1 ].getTT(threadNo);
            else
                b = std::numeric_limits<T1>::max();
        } else {
            if (j!=0)
                b = nodes[ (i+1)*(ncz+1)+j-1 ].getTT(threadNo);
            else
                b = std::numeric_limits<T1>::max();
            if (j!=ncz)
                t = nodes[ (i-1)*(ncz+1)+j+1 ].getTT(threadNo);
            else
                t = std::numeric_limits<T1>::max();
            b = b<t ? b : t;
        }

        T1 fh = 1.414213562373095 * nodes[i*(ncz+1)+j].getNodeSlowness() *
        dx;
        if ( std::abs(a-b) >= fh )
            t = (a<b ? a : b) + fh;
        else
            t = 0.5*( a+b + sqrt(2.*fh*fh - (a-b)*(a-b) ) );

        if ( t<nodes[i*(ncz+1)+j].getTT(threadNo) )
            nodes[i*(ncz+1)+j].setTT(t,threadNo);
    }


    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::update_node_xz(const size_t i, const size_t j,
                                                const size_t threadNo) const {

        T1 a, b, t;
        if (i==0)
            a = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo);
        else if (i==ncx)
            a = nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
        else {
            a = nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
            t = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo);
            a = a<t ? a : t;
        }

        if (j==0)
            b = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo);
        else if (j==ncz)
            b = nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
        else {
            b = nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
            t = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo);
            b = b<t ? b : t;
        }

        if ( a<b && ((b-a)/dx)>nodes[i*(ncz+1)+j].getNodeSlowness() ) {
            t = a + nodes[i*(ncz+1)+j].getNodeSlowness()*dx;
        } else if ( a>b && ((a-b)/dz)>nodes[i*(ncz+1)+j].getNodeSlowness() ) {
            t = b + nodes[i*(ncz+1)+j].getNodeSlowness()*dz;
        } else {
            T1 dx2 = dx*dx;
            T1 dz2 = dz*dz;
            T1 s2 = nodes[i*(ncz+1)+j].getNodeSlowness()*nodes[i*(ncz+1)+j].getNodeSlowness();
            t = (b*dx2 + a*dz2)/(dx2 + dz2) + sqrt((2.0*a*b*dx2*dz2 - a*a*dx2*dz2 -
                                                    b*b*dx2*dz2 + dx2*dx2*dz2*s2 +
                                                    dx2*dz2*dz2*s2)/((dx2 + dz2)*(dx2 + dz2)));
        }

        if ( t<nodes[i*(ncz+1)+j].getTT(threadNo) )
            nodes[i*(ncz+1)+j].setTT(t,threadNo);
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::update_node_weno3(const size_t i, const size_t j,
                                                   const size_t threadNo) const {

        // not valid if dx != dz

        //    @Article{zhang06,
        //        Title                    = {High Order Fast Sweeping Methods for Static {H}amilton–{J}acobi Equations},
        //        Author                   = {Yong-Tao Zhang and Hong-Kai Zhao and Jianliang Qian},
        //        Journal                  = {Journal of Scientific Computing},
        //        Year                     = {2006},
        //        Number                   = {1},
        //        Pages                    = {25--56},
        //        Volume                   = {29},
        //        DOI                      = {10.1007/s10915-005-9014-3},
        //        URL                      = {http://dx.doi.org/10.1007/s10915-005-9014-3}
        //        }


        T1 a, b, t;
        if (i==0) {
            a = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo);  // fist order
        } else if (i==1) {
            T1 num = nodes[ (i+2)*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) + nodes[ i*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 ap = (1.-w)*(nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo)-nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo))/(2.*dx) +
            w*(-nodes[ (i+2)*(ncz+1)+j ].getTT(threadNo) +4.*nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) -3.*nodes[ i*(ncz+1)+j ].getTT(threadNo))/(2.*dx);

            a = nodes[ i*(ncz+1)+j ].getTT(threadNo) + dx*ap;

            t = nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo); // first order for left
            a = a<t ? a : t;

        } else if (i==ncx) {
            a = nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
        } else if (i==ncx-1) {
            T1 num = nodes[ i*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-2)*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 am = (1.-w)*(nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo)-nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo))/(2.*dx) +
            w*(3.*nodes[ i*(ncz+1)+j ].getTT(threadNo) -4.*nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-2)*(ncz+1)+j ].getTT(threadNo))/(2.*dx);

            a = nodes[ i*(ncz+1)+j ].getTT(threadNo) - dx*am;

            t = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo); // first order for right
            a = a<t ? a : t;

        } else {
            T1 num = nodes[ (i+2)*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) + nodes[ i*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 ap = (1.-w)*(nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo)-nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo))/(2.*dx) +
            w*(-nodes[ (i+2)*(ncz+1)+j ].getTT(threadNo) +4.*nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) -3.*nodes[ i*(ncz+1)+j ].getTT(threadNo))/(2.*dx);

            num = nodes[ i*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-2)*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            w = 1./(1.+2.*r*r);

            T1 am = (1.-w)*(nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo)-nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo))/(2.*dx) +
            w*(3.*nodes[ i*(ncz+1)+j ].getTT(threadNo) -4.*nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-2)*(ncz+1)+j ].getTT(threadNo))/(2.*dx);

            a = nodes[ i*(ncz+1)+j ].getTT(threadNo) - dx*am < nodes[ i*(ncz+1)+j ].getTT(threadNo) + dx*ap ?
            nodes[ i*(ncz+1)+j ].getTT(threadNo) - dx*am : nodes[ i*(ncz+1)+j ].getTT(threadNo) + dx*ap;

        }

        if (j==0) {
            b = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo);
        } else if (j==1) {
            T1 num = nodes[ i*(ncz+1)+j+2 ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 bp = (1.-w)*(nodes[ i*(ncz+1)+j+1 ].getTT(threadNo)-nodes[ i*(ncz+1)+j-1 ].getTT(threadNo))/(2.*dx) +
            w*(-nodes[ i*(ncz+1)+j+2 ].getTT(threadNo) +4.*nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) -3.*nodes[ i*(ncz+1)+j ].getTT(threadNo))/(2.*dx);

            b = nodes[ i*(ncz+1)+j ].getTT(threadNo) + dx*bp;

            t = nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
            b = b<t ? b : t;

        } else if (j==ncz) {
            b = nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
        } else if (j==ncz-1) {
            T1 num = nodes[ i*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j-1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j-2 ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 bm = (1.-w)*(nodes[ i*(ncz+1)+j+1 ].getTT(threadNo)-nodes[ i*(ncz+1)+j-1 ].getTT(threadNo))/(2.*dx) +
            w*(3.*nodes[ i*(ncz+1)+j ].getTT(threadNo) -4.*nodes[ i*(ncz+1)+j-1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j-2 ].getTT(threadNo))/(2.*dx);

            b = nodes[ i*(ncz+1)+j ].getTT(threadNo) - dx*bm;

            t = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo); // first order for right
            b = b<t ? b : t;

        } else {
            T1 num = nodes[ i*(ncz+1)+j+2 ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 bp = (1.-w)*(nodes[ i*(ncz+1)+j+1 ].getTT(threadNo)-nodes[ i*(ncz+1)+j-1 ].getTT(threadNo))/(2.*dx) +
            w*(-nodes[ i*(ncz+1)+j+2 ].getTT(threadNo) +4.*nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) -3.*nodes[ i*(ncz+1)+j ].getTT(threadNo))/(2.*dx);

            num = nodes[ i*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j-1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j-2 ].getTT(threadNo);
            num *= num;
            r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            w = 1./(1.+2.*r*r);

            T1 bm = (1.-w)*(nodes[ i*(ncz+1)+j+1 ].getTT(threadNo)-nodes[ i*(ncz+1)+j-1 ].getTT(threadNo))/(2.*dx) +
            w*(3.*nodes[ i*(ncz+1)+j ].getTT(threadNo) -4.*nodes[ i*(ncz+1)+j-1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j-2 ].getTT(threadNo))/(2.*dx);

            b = nodes[ i*(ncz+1)+j ].getTT(threadNo) - dx*bm < nodes[ i*(ncz+1)+j ].getTT(threadNo) + dx*bp ?
            nodes[ i*(ncz+1)+j ].getTT(threadNo) - dx*bm : nodes[ i*(ncz+1)+j ].getTT(threadNo) + dx*bp;

        }

        T1 fh = nodes[i*(ncz+1)+j].getNodeSlowness() * dx;

        if ( std::abs(a-b) >= fh )
            t = (a<b ? a : b) + fh;
        else
            t = 0.5*( a+b + sqrt(2.*fh*fh - (a-b)*(a-b) ) );

        if ( t<nodes[i*(ncz+1)+j].getTT(threadNo) )
            nodes[i*(ncz+1)+j].setTT(t,threadNo);

        if ( t<nodes[i*(ncz+1)+j].getTT(threadNo) )
            nodes[i*(ncz+1)+j].setTT(t,threadNo);

    }


    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::update_node_weno3_xz(const size_t i, const size_t j,
                                                      const size_t threadNo) const {

        T1 a, b, t;
        if (i==0) {
            a = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo);  // fist order
        } else if (i==1) {
            T1 num = nodes[ (i+2)*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) + nodes[ i*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 ap = (1.-w)*(nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo)-nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo))/(2.*dx) +
            w*(-nodes[ (i+2)*(ncz+1)+j ].getTT(threadNo) +4.*nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) -3.*nodes[ i*(ncz+1)+j ].getTT(threadNo))/(2.*dx);

            a = nodes[ i*(ncz+1)+j ].getTT(threadNo) + dx*ap;

            t = nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo); // first order for left
            a = a<t ? a : t;

        } else if (i==ncx) {
            a = nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
        } else if (i==ncx-1) {
            T1 num = nodes[ i*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-2)*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 am = (1.-w)*(nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo)-nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo))/(2.*dx) +
            w*(3.*nodes[ i*(ncz+1)+j ].getTT(threadNo) -4.*nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-2)*(ncz+1)+j ].getTT(threadNo))/(2.*dx);

            a = nodes[ i*(ncz+1)+j ].getTT(threadNo) - dx*am;

            t = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo); // first order for right
            a = a<t ? a : t;

        } else {
            T1 num = nodes[ (i+2)*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) + nodes[ i*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 ap = (1.-w)*(nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo)-nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo))/(2.*dx) +
            w*(-nodes[ (i+2)*(ncz+1)+j ].getTT(threadNo) +4.*nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo) -3.*nodes[ i*(ncz+1)+j ].getTT(threadNo))/(2.*dx);

            num = nodes[ i*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-2)*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            w = 1./(1.+2.*r*r);

            T1 am = (1.-w)*(nodes[ (i+1)*(ncz+1)+j ].getTT(threadNo)-nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo))/(2.*dx) +
            w*(3.*nodes[ i*(ncz+1)+j ].getTT(threadNo) -4.*nodes[ (i-1)*(ncz+1)+j ].getTT(threadNo) + nodes[ (i-2)*(ncz+1)+j ].getTT(threadNo))/(2.*dx);

            a = nodes[ i*(ncz+1)+j ].getTT(threadNo) - dx*am < nodes[ i*(ncz+1)+j ].getTT(threadNo) + dx*ap ?
            nodes[ i*(ncz+1)+j ].getTT(threadNo) - dx*am : nodes[ i*(ncz+1)+j ].getTT(threadNo) + dx*ap;

        }

        if (j==0) {
            b = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo);
        } else if (j==1) {
            T1 num = nodes[ i*(ncz+1)+j+2 ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 bp = (1.-w)*(nodes[ i*(ncz+1)+j+1 ].getTT(threadNo)-nodes[ i*(ncz+1)+j-1 ].getTT(threadNo))/(2.*dz) +
            w*(-nodes[ i*(ncz+1)+j+2 ].getTT(threadNo) +4.*nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) -3.*nodes[ i*(ncz+1)+j ].getTT(threadNo))/(2.*dz);

            b = nodes[ i*(ncz+1)+j ].getTT(threadNo) + dz*bp;

            t = nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
            b = b<t ? b : t;

        } else if (j==ncz) {
            b = nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
        } else if (j==ncz-1) {
            T1 num = nodes[ i*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j-1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j-2 ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 bm = (1.-w)*(nodes[ i*(ncz+1)+j+1 ].getTT(threadNo)-nodes[ i*(ncz+1)+j-1 ].getTT(threadNo))/(2.*dz) +
            w*(3.*nodes[ i*(ncz+1)+j ].getTT(threadNo) -4.*nodes[ i*(ncz+1)+j-1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j-2 ].getTT(threadNo))/(2.*dz);

            b = nodes[ i*(ncz+1)+j ].getTT(threadNo) - dz*bm;

            t = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo); // first order for right
            b = b<t ? b : t;

        } else {
            T1 num = nodes[ i*(ncz+1)+j+2 ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j ].getTT(threadNo) + nodes[ i*(ncz+1)+j-1 ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 bp = (1.-w)*(nodes[ i*(ncz+1)+j+1 ].getTT(threadNo)-nodes[ i*(ncz+1)+j-1 ].getTT(threadNo))/(2.*dz) +
            w*(-nodes[ i*(ncz+1)+j+2 ].getTT(threadNo) +4.*nodes[ i*(ncz+1)+j+1 ].getTT(threadNo) -3.*nodes[ i*(ncz+1)+j ].getTT(threadNo))/(2.*dz);

            num = nodes[ i*(ncz+1)+j ].getTT(threadNo) -2.*nodes[ i*(ncz+1)+j-1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j-2 ].getTT(threadNo);
            num *= num;
            r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            w = 1./(1.+2.*r*r);

            T1 bm = (1.-w)*(nodes[ i*(ncz+1)+j+1 ].getTT(threadNo)-nodes[ i*(ncz+1)+j-1 ].getTT(threadNo))/(2.*dz) +
            w*(3.*nodes[ i*(ncz+1)+j ].getTT(threadNo) -4.*nodes[ i*(ncz+1)+j-1 ].getTT(threadNo) + nodes[ i*(ncz+1)+j-2 ].getTT(threadNo))/(2.*dz);

            b = nodes[ i*(ncz+1)+j ].getTT(threadNo) - dz*bm < nodes[ i*(ncz+1)+j ].getTT(threadNo) + dz*bp ?
            nodes[ i*(ncz+1)+j ].getTT(threadNo) - dz*bm : nodes[ i*(ncz+1)+j ].getTT(threadNo) + dz*bp;

        }

        if ( a<b && ((b-a)/dx)>nodes[i*(ncz+1)+j].getNodeSlowness() ) {
            t = a + nodes[i*(ncz+1)+j].getNodeSlowness()*dx;
        } else if ( a>b && ((a-b)/dz)>nodes[i*(ncz+1)+j].getNodeSlowness() ) {
            t = b + nodes[i*(ncz+1)+j].getNodeSlowness()*dz;
        } else {
            T1 dx2 = dx*dx;
            T1 dz2 = dz*dz;
            T1 s2 = nodes[i*(ncz+1)+j].getNodeSlowness()*nodes[i*(ncz+1)+j].getNodeSlowness();
            t = (b*dx2 + a*dz2)/(dx2 + dz2) + sqrt((2.0*a*b*dx2*dz2 - a*a*dx2*dz2 -
                                                    b*b*dx2*dz2 + dx2*dx2*dz2*s2 +
                                                    dx2*dz2*dz2*s2)/((dx2 + dz2)*(dx2 + dz2)));
        }

        if ( t<nodes[i*(ncz+1)+j].getTT(threadNo) )
            nodes[i*(ncz+1)+j].setTT(t,threadNo);
    }


    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::initFSM(const std::vector<S>& Tx,
                                         const std::vector<T1>& t0,
                                         std::vector<bool>& frozen,
                                         const int npts,
                                         const size_t threadNo) const {

        for (size_t n=0; n<Tx.size(); ++n) {
            bool found = false;
            for ( size_t nn=0; nn<nodes.size(); ++nn ) {
                if ( nodes[nn] == Tx[n] ) {
                    found = true;
                    nodes[nn].setTT( t0[n], threadNo );
                    frozen[nn] = true;

                    ptrdiff_t i = nn/(ncz+1);
                    ptrdiff_t j = nn - i*(ncz+1);

                    for ( ptrdiff_t ii=i-npts; ii<=i+npts; ++ii ) {
                        if ( ii>=0 && ii<=ncx ) {
                            for ( ptrdiff_t jj=j-npts; jj<=j+npts; ++jj ) {
                                if ( jj>=0 && jj<=ncz && !(ii==i && jj==j) ) {

                                    size_t nnn = ii*(ncz+1) + jj;

                                    T1 tt = t0[n] + nodes[nnn].getDistance(Tx[n]) * 0.5*(nodes[nnn].getNodeSlowness() + nodes[nn].getNodeSlowness());
                                    nodes[nnn].setTT( tt, threadNo );
                                    frozen[nnn] = true;
                                }
                            }
                        }
                    }

                    break;
                }
            }
            if ( found==false ) {

                // find cell where Tx resides
                ptrdiff_t cellNo = getCellNo(Tx[n]);

                ptrdiff_t i = cellNo/ncz;
                ptrdiff_t j = cellNo - i*ncz;

                for ( ptrdiff_t ii=i-(npts-1); ii<=i+npts; ++ii ) {
                    if ( ii>=0 && ii<=ncx ) {
                        for ( ptrdiff_t jj=j-(npts-1); jj<=j+npts; ++jj ) {
                            if ( jj>=0 && jj<=ncz ) {
                                size_t nnn = ii*(ncz+1) + jj;

                                T1 tt = t0[n] + nodes[nnn].getDistance(Tx[n]) * nodes[nnn].getNodeSlowness();
                                nodes[nnn].setTT( tt, threadNo );
                                frozen[nnn] = true;
                            }
                        }
                    }
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE>
    T1 Grid2Drn<T1,T2,S,NODE>::getSlowness(const S& pt) const {

        // bilinear interpolation if not on node

        T1 s;
        T2 i, j;

        getIJ(pt, i, j);

        if ( std::abs(pt.x - (xmin+i*dx))<small && std::abs(pt.z - (zmin+j*dz))<small ) {
            // on node
            return nodes[i*(ncz+1)+j].getNodeSlowness();
        } else if ( std::abs(pt.x - (xmin+i*dx))<small ) {

            // on edge
            T1 t1 = nodes[i*(ncz+1)+j].getNodeSlowness();
            T1 t2 = nodes[i*(ncz+1)+j+1].getNodeSlowness();

            T1 w1 = (zmin+(j+1)*dz - pt.z)/dz;
            T1 w2 = (pt.z - (zmin+j*dz))/dz;

            s = t1*w1 + t2*w2;

        } else if ( std::abs(pt.z - (zmin+j*dz))<small ) {

            // on edge
            T1 t1 = nodes[i*(ncz+1)+j].getNodeSlowness();
            T1 t2 = nodes[(i+1)*(ncz+1)+j].getNodeSlowness();

            T1 w1 = (xmin+(i+1)*dx - pt.x)/dx;
            T1 w2 = (pt.x - (xmin+i*dx))/dx;

            s = t1*w1 + t2*w2;

        } else {

            T1 t1 = nodes[    i*(ncz+1)+j  ].getNodeSlowness();
            T1 t2 = nodes[(i+1)*(ncz+1)+j  ].getNodeSlowness();
            T1 t3 = nodes[    i*(ncz+1)+j+1].getNodeSlowness();
            T1 t4 = nodes[(i+1)*(ncz+1)+j+1].getNodeSlowness();

            T1 w1 = (xmin+(i+1)*dx - pt.x)/dx;
            T1 w2 = (pt.x - (xmin+i*dx))/dx;

            t1 = t1*w1 + t2*w2;
            t2 = t3*w1 + t4*w2;

            w1 = (zmin+(j+1)*dz - pt.z)/dz;
            w2 = (pt.z - (zmin+j*dz))/dz;

            s = t1*w1 + t2*w2;
        }

        return s;
    }

    template<typename T1, typename T2, typename S, typename NODE>
    T1 Grid2Drn<T1,T2,S,NODE>::getTraveltimeFromRaypath(const std::vector<S>& Tx,
                                                        const std::vector<T1>& t0,
                                                        const S& Rx,
                                                        const size_t threadNo) const {

        T1 tt = 0.0;
        T1 s1=0.0, s2=0.0, slown=0.0;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return tt;
            }
        }

        S prev_pt( Rx );
        S curr_pt( Rx );
        if ( !this->hasCellSlowness() )
            s1 = getSlowness( curr_pt );
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

            if ( this->hasCellSlowness() ) {
                sxz<T1> mid_pt = static_cast<T1>(0.5) * (prev_pt + curr_pt);
                slown = this->getCellSlowness(getCellNo(mid_pt));
            } else {
                s2 = getSlowness( curr_pt );
                slown = 0.5*(s1 + s2);
                s1 = s2;
            }
            tt += slown * prev_pt.getDistance( curr_pt );
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
                        if ( this->hasCellSlowness() ) {
                            sxz<T1> mid_pt = static_cast<T1>(0.5) * (prev_pt + Tx[ns]);
                            slown = this->getCellSlowness(getCellNo(mid_pt));
                        } else {
                            s2 = getSlowness( Tx[ns] );
                            slown = 0.5*(s1 + s2);
                        }
                        tt += slown * prev_pt.getDistance( Tx[ns] );
                    } else {
                        if ( this->hasCellSlowness() ) {
                            sxz<T1> mid_pt = static_cast<T1>(0.5) * (prev_pt + curr_pt);
                            slown = this->getCellSlowness(getCellNo(mid_pt));
                            tt += slown * prev_pt.getDistance( curr_pt );
                            mid_pt = static_cast<T1>(0.5) * (Tx[ns] + curr_pt);
                            slown = this->getCellSlowness(getCellNo(mid_pt));
                            tt += slown * curr_pt.getDistance( Tx[ns] );
                        } else {
                            // to intersection
                            s2 = getSlowness( curr_pt );
                            tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
                            s1 = s2;
                            // to Tx
                            s2 = getSlowness( Tx[ns] );
                            tt += 0.5*(s1 + s2) * curr_pt.getDistance( Tx[ns] );
                        }
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
        return tt;
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::getRaypath(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            const S& Rx,
                                            std::vector<S>& r_data,
                                            T1 &tt,
                                            const size_t threadNo) const {

        tt = 0.0;
        T1 s1=0.0, s2=0.0, slown=0.0;

        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        S curr_pt( Rx );
        if ( !this->hasCellSlowness() )
            s1 = getSlowness( curr_pt );
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
            if ( this->hasCellSlowness() ) {
                sxz<T1> mid_pt = static_cast<T1>(0.5) * (r_data.back() + curr_pt);
                slown = this->getCellSlowness(getCellNo(mid_pt));
            } else {
                s2 = getSlowness( curr_pt );
                slown = 0.5*(s1 + s2);
                s1 = s2;
            }
            tt += slown * r_data.back().getDistance( curr_pt );
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
                        if ( this->hasCellSlowness() ) {
                            sxz<T1> mid_pt = static_cast<T1>(0.5) * (r_data.back() + Tx[ns]);
                            slown = this->getCellSlowness(getCellNo(mid_pt));
                        } else {
                            s2 = getSlowness( Tx[ns] );
                            slown = 0.5*(s1 + s2);
                        }
                        tt += slown * r_data.back().getDistance( Tx[ns] );
                        r_data.push_back( Tx[ns] );
                    } else {
                        if ( this->hasCellSlowness() ) {
                            sxz<T1> mid_pt = static_cast<T1>(0.5) * (r_data.back() + curr_pt);
                            slown = this->getCellSlowness(getCellNo(mid_pt));
                            tt += slown * r_data.back().getDistance( curr_pt );
                            mid_pt = static_cast<T1>(0.5) * (curr_pt + Tx[ns]);
                            slown = this->getCellSlowness(getCellNo(mid_pt));
                            tt += slown * curr_pt.getDistance( Tx[ns] );
                        } else {
                            // to intersection
                            s2 = getSlowness( curr_pt );
                            tt += 0.5*(s1 + s2) * r_data.back().getDistance( curr_pt );
                            r_data.push_back( curr_pt );
                            s1 = s2;
                            // to Tx
                            s2 = getSlowness( Tx[ns] );
                            tt += 0.5*(s1 + s2) * curr_pt.getDistance( Tx[ns] );
                            r_data.push_back( Tx[ns] );
                        }
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::getRaypath(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            const S& Rx,
                                            std::vector<S>& r_data,
                                            std::vector<siv<T1>> &l_data,
                                            T1 &tt,
                                            const size_t threadNo) const {

        tt = 0.0;
        T1 s1=0.0, s2=0.0, slown=0.0;

        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        S curr_pt( Rx );
        if ( !this->hasCellSlowness() )
            s1 = getSlowness( curr_pt );
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
            if ( this->hasCellSlowness() ) {
                slown = this->getCellSlowness(cell.i);
            } else {
                s2 = getSlowness( curr_pt );
                slown = 0.5*(s1 + s2);
                s1 = s2;
            }
            tt += slown * r_data.back().getDistance( curr_pt );
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
                        if ( this->hasCellSlowness() ) {
                            slown = this->getCellSlowness(cell.i);
                        } else {
                            s2 = getSlowness( Tx[ns] );
                            slown = 0.5*(s1 + s2);
                        }
                        tt += slown * r_data.back().getDistance( Tx[ns] );
                        r_data.push_back( Tx[ns] );
                    } else {
                        // to intersection
                        mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
                        cell.i = getCellNo(mid_pt);
                        cell.v = curr_pt.getDistance(r_data.back());
                        l_data.push_back(cell);
                        if ( this->hasCellSlowness() ) {
                            slown = this->getCellSlowness(cell.i);
                        } else {
                            s2 = getSlowness( curr_pt );
                            slown = 0.5*(s1 + s2);
                            s1 = s2;
                        }
                        tt += slown * r_data.back().getDistance( curr_pt );
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
                        if ( this->hasCellSlowness() ) {
                            slown = this->getCellSlowness(cell.i);
                        } else {
                            s2 = getSlowness( Tx[ns] );
                            slown = 0.5*(s1 + s2);
                        }
                        tt += slown * curr_pt.getDistance( Tx[ns] );
                        r_data.push_back( Tx[ns] );
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE>
    void Grid2Drn<T1,T2,S,NODE>::getRaypath(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            const S& Rx,
                                            std::vector<siv<T1>> &l_data,
                                            T1 &tt,
                                            const size_t threadNo) const {

        tt = 0.0;
        T1 s1=0.0, s2=0.0, slown=0.0;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        S curr_pt( Rx );
        S prev_pt( Rx );
        if ( !this->hasCellSlowness() )
            s1 = getSlowness( curr_pt );
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
            if ( this->hasCellSlowness() ) {
                slown = this->getCellSlowness(cell.i);
            } else {
                s2 = getSlowness( curr_pt );
                slown = 0.5*(s1 + s2);
                s1 = s2;
            }
            tt += slown * cell.v;
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
                        if ( this->hasCellSlowness() ) {
                            slown = this->getCellSlowness(cell.i);
                        } else {
                            s2 = getSlowness( Tx[ns] );
                            slown = 0.5*(s1 + s2);
                        }
                        tt += slown * cell.v;

                    } else {
                        // to intersection
                        mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
                        cell.i = getCellNo(mid_pt);
                        cell.v = curr_pt.getDistance(prev_pt);
                        l_data.push_back(cell);
                        if ( this->hasCellSlowness() ) {
                            slown = this->getCellSlowness(cell.i);
                        } else {
                            s2 = getSlowness( curr_pt );
                            slown = 0.5*(s1 + s2);
                            s1 = s2;
                        }
                        tt += slown * cell.v;

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
                        if ( this->hasCellSlowness() ) {
                            slown = this->getCellSlowness(cell.i);
                        } else {
                            s2 = getSlowness( Tx[ns] );
                            slown = 0.5*(s1 + s2);
                        }
                        tt += slown * cell.v;
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
    }
}

#endif

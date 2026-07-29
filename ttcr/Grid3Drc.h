//
//  Grid3Drc.h
//  ttcr
//
//  Created by Bernard Giroux on 08-04-24.
//
//  Modified by Benoit Larouche on 12-07-20
//  	: now support parallel raytracing from many source points
//  	  on the same 3D grid simultaneously, using OpenMP.
//  	  Secondary nodes are placed on every edge and face of the grid cells.
//
//  	  The velocity model is sampled for each cell and is constant inside the
//        cell.
//
//

//
// Copyright (C) 2012 Bernard Giroux, Benoît Larouche.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/**
 * @file Grid3Drc.h
 * @brief Base class for 3-D rectilinear grids with cell-based slowness.
 *
 * Declares ttcr::Grid3Drc, the 3-D counterpart of ttcr::Grid2Drc and the
 * mid-level base shared by Grid3Drcsp (shortest path), Grid3Drcfs (fast
 * sweeping), Grid3Drcdsp (dynamic shortest path) and Grid3Drcfs_OpenCL.
 *
 * It sits between the abstract ttcr::Grid3D interface and those solvers and
 * owns everything independent of the eikonal method: grid geometry, node and
 * cell containers, point location, slowness accessors and raypath tracing.
 *
 * @section g3drc_numbering Cell numbering — differs from 2-D
 * Cells are numbered with **x varying fastest**, then y, then z:
 * @f[ \mathrm{cellNo} = k\,n_{cx}n_{cy} + j\,n_{cx} + i @f]
 * so stepping one cell in x advances the index by 1, one in y by @f$n_{cx}@f$,
 * and one in z by @f$n_{cx}n_{cy}@f$.
 *
 * @warning This is the **opposite** convention to ttcr::Grid2Drc, where z
 *          varies fastest (@ref g2drc_numbering). Index arithmetic cannot be
 *          carried across from the 2-D classes unchanged. The
 *          "column-wise (z axis)" comment on @ref ttcr::Grid3Drc::cells is
 *          inherited from the 2-D file and does not describe this class.
 *
 * Nodes follow the same convention and are stored primary-first:
 * @f$(n_{cx}+1)(n_{cy}+1)(n_{cz}+1)@f$ cell corners, then any secondary nodes a
 * solver adds.
 *
 * @section g3drc_translate Origin translation
 * Unlike the 2-D family, these grids can be shifted so their origin sits at
 * (0,0,0) — enabled by ttcr::input_parameters::translateOrigin, and useful when
 * a survey's coordinates are large enough that differences lose precision.
 * Methods that take a point therefore also take a flag saying whether that
 * point has already been translated; passing it wrongly silently offsets the
 * result by the origin.
 *
 * @section g3drc_tol Tolerance
 * Point location here compares against @ref ttcr::small2 (@f$10^{-8}@f$),
 * whereas the 2-D classes use @ref ttcr::small (@f$10^{-4}@f$) — another
 * difference not to carry across.
 *
 * @sa Grid3D.h, Grid2Drc.h, Grid3Drn.h, Cell.h, Grid3Drcsp.h
 */

#ifndef ttcr_Grid3Drc_h
#define ttcr_Grid3Drc_h

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
#include <ctime>

#ifdef VTK
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkXMLRectilinearGridWriter.h"
#endif

#include <boost/math/special_functions/sign.hpp>

#include "Grid3D.h"
#include "NodeKDTree.h"

namespace ttcr {

    /**
     * @brief 3-D rectilinear grid holding one slowness value per cell.
     *
     * @tparam T1   floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2   integer type of node and cell indices, normally @c uint32_t.
     * @tparam NODE node type, e.g. ttcr::Node3Dc or ttcr::Node3Dcsp.
     * @tparam CELL cell policy supplying the slowness model, from Cell.h.
     *
     * @note There is no @c S template parameter: 3-D grids always use
     *       @ref sxyz, whereas the 2-D ones are parameterised on the point type
     *       so they can also be draped in 3-D.
     * @note Abstract in practice — it implements everything except @c raytrace.
     */
    template<typename T1, typename T2, typename NODE, typename CELL>
    class Grid3Drc : public Grid3D<T1,T2> {
    public:

        /**
         * @brief Build the grid geometry and allocate its nodes and cells.
         *
         * @param nx   number of cells along x.
         * @param ny   number of cells along y.
         * @param nz   number of cells along z.
         * @param ddx  cell size along x.
         * @param ddy  cell size along y.
         * @param ddz  cell size along z.
         * @param minx x coordinate of the grid origin.
         * @param miny y coordinate of the grid origin.
         * @param minz z coordinate of the grid origin.
         * @param ttrp recompute receiver traveltimes by integrating slowness
         *             along the traced raypath.
         * @param nt   number of threads; sizes each node's traveltime array.
         * @param _translateOrigin shift the grid so its origin is at (0,0,0).
         *             @sa @ref g3drc_translate
         *
         * @post The grid holds @f$(n_x+1)(n_y+1)(n_z+1)@f$ primary nodes and
         *       @f$n_x n_y n_z@f$ cells. Slowness is **not** set — call
         *       @ref setSlowness before raytracing.
         * @warning The cell count passed to the base and used to size @ref nodes
         *          is computed in @c T2 arithmetic (@c nx*ny*nz), not widened to
         *          @c size_t as it is elsewhere in the class, so a grid large
         *          enough to overflow a 32-bit index wraps silently.
         */
        Grid3Drc(const T2 nx, const T2 ny, const T2 nz,
                 const T1 ddx, const T1 ddy, const T1 ddz,
                 const T1 minx, const T1 miny, const T1 minz,
                 const bool ttrp, const size_t nt=1,
                 const bool _translateOrigin=false) :
        Grid3D<T1,T2>(ttrp, nx*ny*nz, nt, _translateOrigin),
        dx(ddx), dy(ddy), dz(ddz),
        xmin(minx), ymin(miny), zmin(minz),
        xmax(minx+nx*ddx), ymax(miny+ny*ddy), zmax(minz+nz*ddz),
        ncx(nx), ncy(ny), ncz(nz),
        nodes(std::vector<NODE>((nx+1)*(ny+1)*(nz+1), NODE(nt))),
        cells(CELL(nx*ny*nz))
        { }

        /// Destructor.
        virtual ~Grid3Drc() {}

        /**
         * @brief Copy the cell slowness values out of the grid.
         * @param[out] slowness resized to @f$n_{cx}n_{cy}n_{cz}@f$ and filled,
         *                      in the numbering of @ref g3drc_numbering.
         */
        void getSlowness(std::vector<T1>& slowness) const {
            if (slowness.size() != static_cast<size_t>(ncx)*ncy*ncz) {
                slowness.resize(static_cast<size_t>(ncx)*ncy*ncz);
            }
            for (size_t n=0; n<slowness.size(); ++n) {
                slowness[n] = cells.getSlowness(n);
            }
        }
        /**
         * @brief Set the slowness of every cell.
         * @param s one slowness per cell, in the numbering of
         *          @ref g3drc_numbering.
         * @throws std::length_error if @p s has the wrong size, plus whatever
         *         the @c CELL policy raises.
         * @note The @c try / @c catch rethrows unchanged and so has no effect.
         */
        void setSlowness(const std::vector<T1>& s) {
            try {
                cells.setSlowness( s );
            } catch (std::exception& e) {
                throw;
            }
        }
        /**
         * @brief Set the elliptical anisotropy ratio @f$\chi@f$ of every cell.
         * @param x one value per cell.
         * @throws std::exception if the @c CELL policy does not model
         *         @f$\chi@f$ — see Cell.h.
         * @note The 3-D base exposes only @f$\chi@f$ and @f$\psi@f$, whereas
         *       ttcr::Grid2Drc also forwards the Thomsen and
         *       weakly-anelliptical parameters.
         */
        void setChi(const std::vector<T1>& x) {
            cells.setChi( x );
        }
        /**
         * @brief Set the elliptical anisotropy ratio @f$\psi@f$ of every cell.
         * @param x one value per cell.
         * @throws std::exception if the @c CELL policy does not model
         *         @f$\psi@f$.
         */
        void setPsi(const std::vector<T1>& x) {
            cells.setPsi( x );
        }

        /**
         * @brief Total number of nodes, primary and secondary.
         * @return @ref nodes size.
         * @note Takes no "primary only" flag, unlike
         *       ttcr::Grid2Drc::getNumberOfNodes.
         */
        size_t getNumberOfNodes() const { return nodes.size(); }
        /// @return Number of cells, @f$n_{cx}n_{cy}n_{cz}@f$.
        size_t getNumberOfCells() const { return static_cast<size_t>(ncx)*ncy*ncz; }

        /**
         * @brief Collect the traveltimes computed at the primary nodes.
         * @param[out] tt       resized to the primary-node count and filled.
         * @param[in]  threadNo thread whose solution to read.
         * @note Relies on the primary nodes occupying the first
         *       @f$(n_{cx}+1)(n_{cy}+1)(n_{cz}+1)@f$ entries, so it simply reads
         *       that prefix rather than testing @c isPrimary as the 2-D version
         *       does.
         */
        void getTT(std::vector<T1>& tt, const size_t threadNo=0) const final {
            size_t nPrimary = static_cast<size_t>(ncx+1) * (ncy+1) * (ncz+1);
            tt.resize(nPrimary);
            for ( size_t n=0; n<nPrimary; ++n ) {
                tt[n] = nodes[n].getTT(threadNo);
            }
        }

        /**
         * @brief Trace the raypath from a receiver back to the source.
         * @param[in]  Tx       source positions.
         * @param[in]  t0       origin time of each source.
         * @param[in]  Rx       receiver position.
         * @param[out] r_data   raypath points, receiver first and source last.
         * @param[out] tt       traveltime along the path.
         * @param[in]  threadNo thread whose solution to use.
         * @pre A @c raytrace call has populated the node traveltimes for
         *      @p threadNo.
         */
        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const sxyz<T1> &Rx,
                        std::vector<sxyz<T1>> &r_data,
                        T1 &tt,
                        const size_t threadNo) const;

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
         * @brief Memory occupied by the neighbour lists, in bytes.
         * @return Total size of every per-cell neighbour vector.
         * @note Counts payload only — the per-vector bookkeeping is ignored, so
         *       this understates the real footprint.
         */
        size_t getNeighborsSize() const {
            size_t n_elem = 0;
            for ( size_t n=0; n<this->neighbors.size(); ++n ) {
                n_elem += this->neighbors[n].size();
            }
            return n_elem*sizeof(size_t);
        }

        /**
         * @brief Memory occupied by the nodes, in bytes.
         * @return Sum of every node's reported size.
         * @note Only as accurate as @c NODE::getSize, which overestimates for
         *       the non-shortest-path node types — see ttcr::Node3Dc::getSize.
         */
        size_t getNodesSize() const {
            size_t size = 0;
            for ( size_t n=0; n<nodes.size(); ++n ) {
                size += nodes[n].getSize();
            }
            return size;
        }

        const T1 getXmin() const { return xmin; }  ///< @return x coordinate of the grid origin.
        const T1 getYmin() const { return ymin; }  ///< @return y coordinate of the grid origin.
        const T1 getZmin() const { return zmin; }  ///< @return z coordinate of the grid origin.
        const T1 getDx() const { return dx; }      ///< @return Cell size along x.
        const T1 getDy() const { return dy; }      ///< @return Cell size along y.
        const T1 getDz() const { return dz; }      ///< @return Cell size along z.
        const T2 getNcx() const { return ncx; }    ///< @return Number of cells along x.
        const T2 getNcy() const { return ncy; }    ///< @return Number of cells along y.
        const T2 getNcz() const { return ncz; }    ///< @return Number of cells along z.

        /**
         * @brief Write the coordinates of the secondary nodes to a stream.
         * @param os destination stream, one @c "x y z" triple per line.
         * @note Dumps the tail of the node vector past the primary nodes, so it
         *       writes nothing when the solver adds none. Backs the @c -s
         *       command-line switch.
         */
        void dump_secondary(std::ofstream& os) const {
            size_t nPrimary = static_cast<size_t>(ncx+1) * (ncy+1) * (ncz+1);
            for ( size_t n=nPrimary; n<nodes.size(); ++n ) {
                os << nodes[n].getX() << ' ' << nodes[n].getY() << ' ' << nodes[n].getZ() << '\n';
            }
        }

        /**
         * @brief Slowness at an arbitrary point.
         * @param pt           point to sample. Taken by value, since it may be
         *                     translated internally.
         * @param isTranslated true if @p pt is already expressed in the
         *                     translated frame; false if it is in the caller's
         *                     original coordinates and must be shifted first.
         * @return Slowness of the cell containing @p pt.
         * @note Piecewise constant — every point inside a cell gives the same
         *       value. @sa @ref g3drc_translate for the flag's meaning.
         */
        T1 computeSlowness(sxyz<T1> pt, const bool isTranslated=false) const {
            if (this->translateOrigin == true && isTranslated == false) {
                pt -= this->origin;
            }
            return cells.getSlowness(getCellNo(pt));
        }

    protected:
        T1 dx;                   ///< cell size in x
        T1 dy;			         ///< cell size in y
        T1 dz;                   ///< cell size in z
        T1 xmin;                 ///< x origin of the grid
        T1 ymin;                 ///< y origin of the grid
        T1 zmin;                 ///< z origin of the grid
        T1 xmax;                 ///< x end of the grid
        T1 ymax;                 ///< y end of the grid
        T1 zmax;                 ///< z end of the grid
        T2 ncx;                  ///< number of cells in x
        T2 ncy;                  ///< number of cells in y
        T2 ncz;                  ///< number of cells in z

        /// Primary nodes first, then any secondary nodes the derived solver
        /// added. @c mutable because @c const raytracing methods update the
        /// traveltimes stored in them.
        mutable std::vector<NODE> nodes;

        /// Slowness model, one value per cell, in the **x-fastest** order of
        /// @ref g3drc_numbering. Its type is the @c CELL policy, so it also
        /// carries any anisotropy parameters.
        /// @warning The trailing source comment describing this as
        ///          "column-wise (z axis)" is inherited from Grid2Drc.h and is
        ///          wrong for 3-D; and @c Grid3Dcinterp no longer exists.
        CELL cells;   // column-wise (z axis) slowness vector of the cells, NOT used by Grid3Dcinterp

        /// kd-tree over all nodes (primary + secondary), replacing the
        /// brute-force scan that checks whether a point coincides with a node.
        /// Built lazily on first use by @ref getNodeNo.
        // kd-tree over all nodes (primary + secondary), to replace the
        // brute-force scan that checks whether a point coincides with a node.
        // Built lazily on first use; node coordinates are fixed after
        // construction and queries are read-only, hence thread-safe.
        mutable std::unique_ptr<NodeKDTree<T1,T2>> kdtree;
        mutable std::once_flag kdtreeFlag;  ///< Guards the one-time @ref kdtree build.

        /**
         * @brief Find the node coinciding with a point, if any.
         * @param pt point to test.
         * @return Index of the node at @p pt, or @c T2's maximum if none lies
         *         there.
         * @note Builds @ref kdtree on first call; @c std::call_once makes that
         *       safe under concurrency, and later queries need no
         *       synchronisation since node coordinates are then fixed.
         * @note "Coincides" means within @ref small, via the node's
         *       @c operator== — note this is @ref small, not the @ref small2
         *       used by @ref getCellNo in this class.
         */
        // Index of the node at pt, or npos if no node lies within "small" of pt
        // (same condition as Node::operator==).
        T2 getNodeNo(const sxyz<T1>& pt) const {
            std::call_once(kdtreeFlag, [this]() {
                kdtree.reset(new NodeKDTree<T1,T2>(nodes, nodes.size()));
            });
            T2 nn = kdtree->findNearest(pt.x, pt.y, pt.z);
            return nodes[nn] == pt ? nn : std::numeric_limits<T2>::max();
        }

        /**
         * @brief Create the primary nodes, and optionally secondary ones.
         * @param nsnx secondary nodes per cell edge along x.
         * @param nsny secondary nodes per cell edge along y.
         * @param nsnz secondary nodes per cell edge along z.
         * @note Lives in the base here, unlike the 2-D family where each solver
         *       supplies its own; the defaults of 0 give a primary-only grid,
         *       which is what the sweeping solvers want.
         */
        void buildGridNodes(const T2 nsnx=0, const T2 nsny=0, const T2 nsnz=0);

        /**
         * @brief Locate the cell containing a point.
         * @param pt point to locate.
         * @return Cell number, per @ref g3drc_numbering.
         * @note A point on the far boundary (within @ref small2 of @c xmax,
         *       @c ymax or @c zmax) is pulled back into the last cell.
         * @warning A point outside the grid is not detected — the arithmetic
         *          simply yields an out-of-range number. Screen with
         *          @ref checkPts first.
         * @pre @p pt is in the translated frame if origin translation is
         *      enabled. @sa @ref g3drc_translate
         */
        T2 getCellNo(const sxyz<T1>& pt) const {
            T1 x = xmax-pt.x < small2 ? xmax-.5*dx : pt.x;
            T1 y = ymax-pt.y < small2 ? ymax-.5*dy : pt.y;
            T1 z = zmax-pt.z < small2 ? zmax-.5*dz : pt.z;
            T2 nx = static_cast<T2>( small2 + (x-xmin)/dx );
            T2 ny = static_cast<T2>( small2 + (y-ymin)/dy );
            T2 nz = static_cast<T2>( small2 + (z-zmin)/dz );
            return ny*ncx + nz*(ncx*ncy) + nx;
        }


        /**
         * @brief Locate the cell containing a node.
         * @param node node to locate.
         * @return Cell number, per @ref g3drc_numbering.
         * @note Same computation as the point overload, reading the coordinates
         *       through the node's accessors. A node sitting exactly on a cell
         *       boundary belongs to the cell on the higher-index side, except at
         *       the far edge of the grid.
         */
        T2 getCellNo(const NODE& node) const {
            T1 x = xmax-node.getX() < small2 ? xmax-.5*dx : node.getX();
            T1 y = ymax-node.getY() < small2 ? ymax-.5*dy : node.getY();
            T1 z = zmax-node.getZ() < small2 ? zmax-.5*dz : node.getZ();
            T2 nx = static_cast<T2>( small2 + (x-xmin)/dx );
            T2 ny = static_cast<T2>( small2 + (y-ymin)/dy );
            T2 nz = static_cast<T2>( small2 + (z-zmin)/dz );
            return ny*ncx + nz*(ncx*ncy) + nx;
        }

        /**
         * @brief Split a cell number into its x, y and z indices.
         * @param[in]  cellNo cell number.
         * @param[out] ind    @c i receives the x index, @c j the y, @c k the z.
         * @sa @ref g3drc_numbering — this inverts @ref getCellNo.
         */
        void getCellIJK(const T2 cellNo, sijk<T2> &ind) const {
            ind.k = cellNo / (ncx*ncy);
            ind.j = (cellNo - ind.k*ncx*ncy) / ncx;
            ind.i = cellNo - ncx * ( ind.k*ncy + ind.j);
        }

        /**
         * @brief Compute the x, y and z cell indices of a point.
         * @param[in]  pt point to locate.
         * @param[out] i  x index.
         * @param[out] j  y index.
         * @param[out] k  z index.
         * @warning Unchecked: a point outside the grid yields an out-of-range
         *          index, and on an unsigned @c T2 a point below the origin
         *          wraps to a huge value.
         */
        void getIJK(const sxyz<T1>& pt, T2& i, T2& j, T2& k) const {
            i = static_cast<T2>( small2 + (pt.x-xmin)/dx );
            j = static_cast<T2>( small2 + (pt.y-ymin)/dy );
            k = static_cast<T2>( small2 + (pt.z-zmin)/dz );
        }

        /**
         * @brief Compute the x, y and z cell indices of a point, as signed values.
         * @param[in]  pt point to locate.
         * @param[out] i  x index.
         * @param[out] j  y index.
         * @param[out] k  z index.
         * @note Overload returning @c ptrdiff_t, for callers that step off the
         *       grid and need a negative index rather than a wrap.
         */
        void getIJK(const sxyz<T1>& pt, ptrdiff_t& i, ptrdiff_t& j, ptrdiff_t& k) const {
            i = static_cast<ptrdiff_t>( small2 + (pt.x-xmin)/dx );
            j = static_cast<ptrdiff_t>( small2 + (pt.y-ymin)/dy );
            k = static_cast<ptrdiff_t>( small2 + (pt.z-zmin)/dz );
        }

        /**
         * @brief Verify that every point lies inside the grid.
         * @param pts        points to check, typically sources or receivers.
         *                   Taken **by value**, since they may be translated.
         * @param translated true if @p pts are already in the translated frame.
         * @throws std::runtime_error naming the offending point if one falls
         *         outside the grid bounds.
         * @sa @ref g3drc_translate
         */
        void checkPts(std::vector<sxyz<T1>> pts, const bool translated=false) const;

        /**
         * @brief Traveltime at a receiver.
         * @param Rx       receiver position.
         * @param nodes    node vector to read from — passed explicitly so the
         *                 dynamic solvers can supply a set that includes their
         *                 temporary nodes.
         * @param threadNo thread whose solution to read.
         * @return Traveltime at @p Rx, interpolated if it does not sit on a node.
         */
        T1 getTraveltime(const sxyz<T1>& Rx,
                         const std::vector<NODE>& nodes,
                         const size_t threadNo) const;

        /**
         * @brief Traveltime at a receiver, also reporting where it came from.
         * @param[in]  Rx           receiver position.
         * @param[in]  nodes        node vector to read from.
         * @param[out] nodeParentRx index of the node the ray arrived from.
         * @param[out] cellParentRx index of the cell it crossed to get there.
         * @param[in]  threadNo     thread whose solution to read.
         * @return Traveltime at @p Rx.
         * @note The two output parameters are unnamed in this declaration; the
         *       names used here are those of the out-of-class definition. They
         *       are the starting point for walking a raypath back to the source.
         */
        T1 getTraveltime(const sxyz<T1>& Rx,
                         const std::vector<NODE>& nodes,
                         T2& nodeParentRx, T2& cellParentRx,
                         const size_t threadNo) const;

        /**
         * @brief Traveltime obtained by integrating slowness along the raypath.
         * @param Tx       source positions.
         * @param t0       origin time of each source.
         * @param Rx       receiver position.
         * @param threadNo thread whose solution to use.
         * @return Traveltime at @p Rx.
         * @note Selected by the @c ttrp constructor flag. More faithful than
         *       interpolation because it follows the path rather than smoothing
         *       across cells, at the cost of a full raypath trace per receiver.
         */
        T1 getTraveltimeFromRaypath(const std::vector<sxyz<T1>>& Tx,
                                    const std::vector<T1>& t0,
                                    const sxyz<T1>& Rx,
                                    const size_t threadNo) const;

        /**
         * @brief Record the per-cell path lengths without returning the geometry.
         * @param[in]  Tx       source positions.
         * @param[in]  t0       origin time of each source.
         * @param[in]  Rx       receiver position.
         * @param[out] l_data   per-cell segment lengths — one row of the
         *                      sensitivity matrix. @sa siv
         * @param[out] tt       traveltime along the path.
         * @param[in]  threadNo thread whose solution to use.
         */
        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const sxyz<T1> &Rx,
                        std::vector<siv<T1>> &l_data,
                        T1 &tt,
                        const size_t threadNo) const;

        /**
         * @brief Trace the raypath and record the length travelled in each cell.
         * @param[in]  Tx       source positions.
         * @param[in]  t0       origin time of each source.
         * @param[in]  Rx       receiver position.
         * @param[out] r_data   raypath points.
         * @param[out] l_data   per-cell segment lengths.
         * @param[out] tt       traveltime along the path.
         * @param[in]  threadNo thread whose solution to use.
         * @note Used when ttcr::input_parameters::saveM is set.
         */
        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const sxyz<T1> &Rx,
                        std::vector<sxyz<T1>> &r_data,
                        std::vector<siv<T1>> &l_data,
                        T1 &tt,
                        const size_t threadNo) const;

    private:
        /// @name Non-copyable
        /// Private, the pre-C++11 way of suppressing default construction and
        /// copying.
        /// @note @c operator= returns @c *this without copying anything, so it
        ///       is not undefined behaviour — but it would silently do nothing
        ///       if reached from inside the class. `= delete` would make that a
        ///       compile error.
        /// @{
        Grid3Drc() {}
        Grid3Drc(const Grid3Drc<T1,T2,NODE,CELL>& g) {}
        Grid3Drc<T1,T2,NODE,CELL>& operator=(const Grid3Drc<T1,T2,NODE,CELL>& g) { return *this; }
        /// @}

        T1 getTraveltime(const sxyz<T1>& pt,
                         const size_t threadNo) const final;

        void grad(sxyz<T1>& g, const sxyz<T1> &pt,
                  const size_t nt) const;

    };


    template<typename T1, typename T2, typename NODE, typename CELL>
    void Grid3Drc<T1,T2,NODE, CELL>::buildGridNodes(const T2 nsnx,
                                                    const T2 nsny,
                                                    const T2 nsnz) {

        if ( nsnx != 0 || nsny != 0 || nsnz != 0) {
            nodes.resize(// secondary nodes on the edges
                         ncx*nsnx*((ncy+1)*(ncz+1)) +
                         ncy*nsny*((ncx+1)*(ncz+1)) +
                         ncz*nsnz*((ncx+1)*(ncy+1)) +
                         // secondary nodes on the faces
                         (nsnx*nsny)*(ncx*ncy*(ncz+1))+
                         (nsnx*nsnz)*(ncx*ncz*(ncy+1))+
                         (nsny*nsnz)*(ncy*ncz*(ncx+1))+
                         // primary nodes
                         (ncx+1) * (ncy+1) * (ncz+1),
                         NODE(this->nThreads));
        }

        if ( this->translateOrigin ) {
            this->origin = {xmin, ymin, zmin};
            xmax -= xmin;
            ymax -= ymin;
            zmax -= zmin;
            xmin = 0.0;
            ymin = 0.0;
            zmin = 0.0;
        } else {
            this->origin = {0.0, 0.0, 0.0};
        }

        // Create the grid, assign a number for each node and find the owners
        // Nodes and cells are first indexed in z, then y, and x.
        // Secondary nodes are placed on the faces and edges of every cells.
        // Ex: the node in "node[A]=(i,j,k)" is followed by the node in
        // "node[A+1]=(i+dx,j,k)"

        T2 cXmYmZm;     // cell in the (x-,y-,z-) direction from the node
        T2 cXpYmZm;     // cell in the (x+,y-,z-) direction from the node
        T2 cXmYpZm;
        T2 cXpYpZm;
        T2 cXmYmZp;
        T2 cXpYmZp;
        T2 cXmYpZp;
        T2 cXpYpZp;

        T2 n = 0;
        for ( T2 nk=0; nk<=ncz; ++nk ) {

            T1 z = zmin + nk*dz;

            for ( T2 nj=0; nj<=ncy; ++nj ) {

                T1 y = ymin + nj*dy;

                for ( T2 ni=0; ni<=ncx; ++ni, ++n ){

                    T1 x = xmin + ni*dx;

                    // Find the adjacent cells for each primary node

                    if (ni < ncx && nj < ncy && nk < ncz){
                        cXpYpZp = nj*ncx + nk*(ncx*ncy) + ni;
                    }
                    else {
                        cXpYpZp = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj < ncy && nk < ncz){
                        cXmYpZp = nj*ncx + nk*(ncx*ncy) + ni - 1;
                    }
                    else {
                        cXmYpZp = std::numeric_limits<T2>::max();
                    }

                    if (ni < ncx && nj > 0 && nk < ncz){
                        cXpYmZp = (nj-1)*ncx + nk*(ncx*ncy) + ni;
                    }
                    else {
                        cXpYmZp = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj > 0 && nk < ncz){
                        cXmYmZp = (nj-1)*ncx + nk*(ncx * ncy) + ni - 1;
                    }
                    else {
                        cXmYmZp = std::numeric_limits<T2>::max();
                    }

                    if (ni < ncx && nj < ncy && nk > 0){
                        cXpYpZm = nj*ncx + (nk-1)*(ncx*ncy) + ni;
                    }
                    else {
                        cXpYpZm = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj < ncy && nk > 0){
                        cXmYpZm = nj*ncx + (nk-1)*(ncx*ncy) + ni - 1;
                    }
                    else {
                        cXmYpZm = std::numeric_limits<T2>::max();
                    }

                    if (ni < ncx && nj > 0 && nk > 0){
                        cXpYmZm = (nj-1)*ncx + (nk-1)*(ncx*ncy) + ni;
                    }
                    else {
                        cXpYmZm = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj > 0 && nk > 0){
                        cXmYmZm = (nj-1)*ncx + (nk-1)*(ncx*ncy) + ni-1;
                    }
                    else {
                        cXmYmZm = std::numeric_limits<T2>::max();
                    }


                    // Index the primary nodes owners

                    if ( cXmYmZm != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cXmYmZm );
                    }
                    if ( cXpYmZm != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cXpYmZm );
                    }
                    if ( cXmYpZm != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cXmYpZm );
                    }
                    if ( cXpYpZm != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cXpYpZm );
                    }
                    if ( cXmYmZp != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cXmYmZp );
                    }
                    if ( cXpYmZp != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cXpYmZp );
                    }
                    if ( cXmYpZp != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cXmYpZp );
                    }
                    if ( cXpYpZp != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cXpYpZp );
                    }

                    nodes[n].setXYZindex( x, y, z, n );
                    nodes[n].setPrimary(true);
                }
            }
        }

        if ( nsnx != 0 || nsny != 0 || nsnz != 0) {
            T1 dxs = dx/(nsnx+1);     // distance between secondary nodes in x
            T1 dys = dy/(nsny+1);
            T1 dzs = dz/(nsnz+1);

            for ( T2 nk=0; nk<=ncz; ++nk ) {

                T1 z = zmin + nk*dz;

                for ( T2 nj=0; nj<=ncy; ++nj ) {

                    T1 y = ymin + nj*dy;

                    for ( T2 ni=0; ni<=ncx; ++ni ){

                        T1 x = xmin + ni*dx;

                        // Find the adjacent cells for each primary node

                        if (ni < ncx && nj < ncy && nk < ncz){
                            cXpYpZp = nj*ncx + nk*(ncx*ncy) + ni;
                        }
                        else {
                            cXpYpZp = std::numeric_limits<T2>::max();
                        }

                        if (ni > 0 && nj < ncy && nk < ncz){
                            cXmYpZp = nj*ncx + nk*(ncx*ncy) + ni - 1;
                        }
                        else {
                            cXmYpZp = std::numeric_limits<T2>::max();
                        }

                        if (ni < ncx && nj > 0 && nk < ncz){
                            cXpYmZp = (nj-1)*ncx + nk*(ncx*ncy) + ni;
                        }
                        else {
                            cXpYmZp = std::numeric_limits<T2>::max();
                        }

                        if (ni > 0 && nj > 0 && nk < ncz){
                            cXmYmZp = (nj-1)*ncx + nk*(ncx * ncy) + ni - 1;
                        }
                        else {
                            cXmYmZp = std::numeric_limits<T2>::max();
                        }

                        if (ni < ncx && nj < ncy && nk > 0){
                            cXpYpZm = nj*ncx + (nk-1)*(ncx*ncy) + ni;
                        }
                        else {
                            cXpYpZm = std::numeric_limits<T2>::max();
                        }

                        if (ni > 0 && nj < ncy && nk > 0){
                            cXmYpZm = nj*ncx + (nk-1)*(ncx*ncy) + ni - 1;
                        }
                        else {
                            cXmYpZm = std::numeric_limits<T2>::max();
                        }

                        if (ni < ncx && nj > 0 && nk > 0){
                            cXpYmZm = (nj-1)*ncx + (nk-1)*(ncx*ncy) + ni;
                        }
                        else {
                            cXpYmZm = std::numeric_limits<T2>::max();
                        }

                        if (ni > 0 && nj > 0 && nk > 0){
                            cXmYmZm = (nj-1)*ncx + (nk-1)*(ncx*ncy) + ni-1;
                        }
                        else {
                            cXmYmZm = std::numeric_limits<T2>::max();
                        }

                        // Secondary nodes on x edge
                        if ( ni < ncx ) {
                            for (T2 ns=0; ns< nsnx; ++ns, ++n ) {

                                T1 xsv = xmin + ni* dx + (ns+1)*dxs;

                                if ( cXpYmZm != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXpYmZm );
                                }
                                if ( cXpYpZm != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXpYpZm );
                                }
                                if ( cXpYmZp != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXpYmZp );
                                }
                                if ( cXpYpZp != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXpYpZp );
                                }
                                nodes[n].setXYZindex( xsv, y, z, n );
                            }
                        }

                        // Secondary nodes on y edge
                        if ( nj < ncy ) {
                            for (T2 ns=0; ns< nsny; ++ns, ++n ) {

                                T1 ysv = ymin + nj* dy + (ns+1)*dys;

                                if ( cXmYpZm != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXmYpZm );
                                }
                                if ( cXpYpZm != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXpYpZm );
                                }
                                if ( cXmYpZp != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXmYpZp );
                                }
                                if ( cXpYpZp != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXpYpZp );
                                }
                                nodes[n].setXYZindex( x, ysv, z, n );
                            }
                        }

                        // Secondary nodes on z edge
                        if ( nk < ncz ) {
                            for (T2 ns=0; ns< nsnz; ++ns, ++n ) {

                                T1 zsv = zmin + nk* dz + (ns+1)*dzs;

                                if ( cXmYmZp != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXmYmZp );
                                }
                                if ( cXpYmZp != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXpYmZp );
                                }
                                if ( cXmYpZp != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXmYpZp );
                                }
                                if ( cXpYpZp != std::numeric_limits<T2>::max() )
                                {
                                    nodes[n].pushOwner( cXpYpZp );
                                }
                                nodes[n].setXYZindex( x, y, zsv, n );
                            }
                        }

                        // Secondary nodes on the xy0 planes
                        if ( ni < ncx && nj < ncy ) {
                            for ( T2 sy=0; sy < nsny; ++sy ) {
                                for ( T2 sx=0; sx < nsnx; ++sx, n++ ) {

                                    T1 ysv = ymin+ nj* dy+ (sy+1)*dys;
                                    T1 xsv = xmin+ ni* dx+ (sx+1)*dxs;

                                    if ( cXpYpZm != std::numeric_limits<T2>::max() )
                                    {
                                        nodes[n].pushOwner( cXpYpZm );
                                    }
                                    if ( cXpYpZp != std::numeric_limits<T2>::max() )
                                    {
                                        nodes[n].pushOwner( cXpYpZp );
                                    }
                                    nodes[n].setXYZindex( xsv, ysv, z, n );
                                }
                            }
                        }

                        // Secondary nodes on the x0z planes
                        if ( ni < ncx && nk < ncz ) {
                            for ( T2 sz=0; sz < nsnz; ++sz ) {
                                for ( T2 sx=0; sx < nsnx; ++sx, n++ ) {

                                    T1 zsv = zmin+ nk* dz+ (sz+1)*dzs;
                                    T1 xsv = xmin+ ni* dx+ (sx+1)*dxs;

                                    if ( cXpYmZp != std::numeric_limits<T2>::max() )
                                    {
                                        nodes[n].pushOwner( cXpYmZp );
                                    }
                                    if ( cXpYpZp != std::numeric_limits<T2>::max() )
                                    {
                                        nodes[n].pushOwner( cXpYpZp );
                                    }
                                    nodes[n].setXYZindex( xsv, y, zsv, n );
                                }
                            }
                        }

                        // Secondary nodes on the 0yz planes
                        if ( nj < ncy && nk < ncz ) {
                            for ( T2 sz=0; sz < nsnz; ++sz ) {
                                for ( T2 sy=0; sy < nsny; ++sy, n++ ) {

                                    T1 zsv = zmin+ nk* dz+ (sz+1)*dzs;
                                    T1 ysv = ymin+ nj* dy+ (sy+1)*dys;

                                    if ( cXmYpZp != std::numeric_limits<T2>::max() )
                                    {
                                        nodes[n].pushOwner( cXmYpZp );
                                    }
                                    if ( cXpYpZp != std::numeric_limits<T2>::max() )
                                    {
                                        nodes[n].pushOwner( cXpYpZp );
                                    }
                                    nodes[n].setXYZindex( x, ysv, zsv, n );
                                }
                            }
                        }
                    }
                }
            }
        }
        // sanity check
        if ( n != nodes.size() ) {
            std::cerr << "Error building grid, wrong number of nodes\n";
            abort();
        }
    }


    template<typename T1, typename T2, typename NODE, typename CELL>
    void Grid3Drc<T1,T2,NODE,CELL>::checkPts(std::vector<sxyz<T1>> pts, const bool translated) const {

        if (this->translateOrigin == true && translated == false) {
            for ( size_t n=0; n<pts.size(); ++n ) {
                pts[n] -= this->origin;
            }
        }

        // Check if the points from a vector are in the grid
        for ( size_t n=0; n<pts.size(); ++n ) {
            if ( pts[n].x < xmin || pts[n].x > xmax ||
                pts[n].y < ymin || pts[n].y > ymax ||
                pts[n].z < zmin || pts[n].z > zmax ) {
                std::ostringstream msg;
                msg << "Error: Point (" << pts[n] << ") outside grid.";
                throw std::runtime_error(msg.str());
            }
        }
    }


    template<typename T1, typename T2, typename NODE, typename CELL>
    T1 Grid3Drc<T1,T2,NODE,CELL>::getTraveltime(const sxyz<T1>& Rx,
                                                const std::vector<NODE>& nodes,
                                                const size_t threadNo) const {

        // Calculate and return the traveltime for a Rx point.
        T2 nn = getNodeNo( Rx );
        if ( nn != std::numeric_limits<T2>::max() ) {
            return nodes[nn].getTT(threadNo);
        }
        size_t cellNo = getCellNo( Rx );
        size_t neibNo = this->neighbors[cellNo][0];
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

    template<typename T1, typename T2, typename NODE, typename CELL>
    T1 Grid3Drc<T1,T2,NODE,CELL>::getTraveltime(const sxyz<T1>& Rx,
                                                const std::vector<NODE>& nodes,
                                                T2& nodeParentRx, T2& cellParentRx,
                                                const size_t threadNo) const {

        // Calculate and return the traveltime for a Rx point.
        T2 nn = getNodeNo( Rx );
        if ( nn != std::numeric_limits<T2>::max() ) {
            nodeParentRx = nodes[nn].getNodeParent(threadNo);
            cellParentRx = nodes[nn].getCellParent(threadNo);
            return nodes[nn].getTT(threadNo);
        }
        T2 cellNo = getCellNo( Rx );
        T2 neibNo = this->neighbors[cellNo][0];
        T1 dt = cells.computeDt(nodes[neibNo], Rx, cellNo);

        T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
        nodeParentRx = neibNo;
        cellParentRx = cellNo;
        for ( size_t k=1; k< this->neighbors[cellNo].size(); ++k ) {
            neibNo = this->neighbors[cellNo][k];
            dt = cells.computeDt(nodes[neibNo], Rx, cellNo);
            if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
                traveltime =  nodes[neibNo].getTT(threadNo)+dt;
                nodeParentRx = neibNo;
            }
        }
        return traveltime;
    }


    template<typename T1, typename T2, typename NODE, typename CELL>
    T1 Grid3Drc<T1,T2,NODE,CELL>::getTraveltime(const sxyz<T1> &pt,
                                                const size_t nt) const {

        const size_t nnx = ncx+1;
        const size_t nny = ncy+1;

        // trilinear interpolation if not on node

        T1 tt;
        T2 i, j, k;

        getIJK(pt, i, j, k);

        if ( std::abs(pt.x - (xmin+i*dx))<small2 &&
            std::abs(pt.y - (ymin+j*dy))<small2 &&
            std::abs(pt.z - (zmin+k*dz))<small2 ) {
            // on node
            return nodes[(k*nny+j)*nnx+i].getTT(nt);
        } else if ( std::abs(pt.x - (xmin+i*dx))<small2 &&
                   std::abs(pt.y - (ymin+j*dy))<small2 ) {
            // on edge
            T1 t1 = nodes[(    k*nny+j)*nnx+i].getTT(nt);
            T1 t2 = nodes[((k+1)*nny+j)*nnx+i].getTT(nt);

            T1 w1 = (zmin+(k+1)*dz - pt.z)/dz;
            T1 w2 = (pt.z - (zmin+k*dz))/dz;

            tt = t1*w1 + t2*w2;

        } else if ( std::abs(pt.x - (xmin+i*dx))<small2 &&
                   std::abs(pt.z - (zmin+k*dz))<small2 ) {
            // on edge
            T1 t1 = nodes[(k*nny+j  )*nnx+i].getTT(nt);
            T1 t2 = nodes[(k*nny+j+1)*nnx+i].getTT(nt);

            T1 w1 = (ymin+(j+1)*dy - pt.y)/dy;
            T1 w2 = (pt.y - (ymin+j*dy))/dy;

            tt = t1*w1 + t2*w2;

        } else if ( std::abs(pt.y - (ymin+j*dy))<small2 &&
                   std::abs(pt.z - (zmin+k*dz))<small2 ) {
            // on edge
            T1 t1 = nodes[(k*nny+j)*nnx+i  ].getTT(nt);
            T1 t2 = nodes[(k*nny+j)*nnx+i+1].getTT(nt);

            T1 w1 = (xmin+(i+1)*dx - pt.x)/dx;
            T1 w2 = (pt.x - (xmin+i*dx))/dx;

            tt = t1*w1 + t2*w2;

        } else if ( std::abs(pt.x - (xmin+i*dx))<small2 ) {
            // on YZ face
            T1 t1 = nodes[(    k*nny+j  )*nnx+i].getTT(nt);
            T1 t2 = nodes[((k+1)*nny+j  )*nnx+i].getTT(nt);
            T1 t3 = nodes[(    k*nny+j+1)*nnx+i].getTT(nt);
            T1 t4 = nodes[((k+1)*nny+j+1)*nnx+i].getTT(nt);

            T1 w1 = (zmin+(k+1)*dz - pt.z)/dz;
            T1 w2 = (pt.z - (zmin+k*dz))/dz;

            t1 = t1*w1 + t2*w2;
            t2 = t3*w1 + t4*w2;

            w1 = (ymin+(j+1)*dy - pt.y)/dy;
            w2 = (pt.y - (ymin+j*dy))/dy;

            tt = t1*w1 + t2*w2;

        } else if ( std::abs(pt.y - (ymin+j*dy))<small2 ) {
            // on XZ face
            T1 t1 = nodes[(    k*nny+j)*nnx+i  ].getTT(nt);
            T1 t2 = nodes[((k+1)*nny+j)*nnx+i  ].getTT(nt);
            T1 t3 = nodes[(    k*nny+j)*nnx+i+1].getTT(nt);
            T1 t4 = nodes[((k+1)*nny+j)*nnx+i+1].getTT(nt);

            T1 w1 = (zmin+(k+1)*dz - pt.z)/dz;
            T1 w2 = (pt.z - (zmin+k*dz))/dz;

            t1 = t1*w1 + t2*w2;
            t2 = t3*w1 + t4*w2;

            w1 = (xmin+(i+1)*dx - pt.x)/dx;
            w2 = (pt.x - (xmin+i*dx))/dx;

            tt = t1*w1 + t2*w2;

        } else if ( std::abs(pt.z - (zmin+k*dz))<small2 ) {
            // on XY face
            T1 t1 = nodes[(k*nny+j  )*nnx+i  ].getTT(nt);
            T1 t2 = nodes[(k*nny+j+1)*nnx+i  ].getTT(nt);
            T1 t3 = nodes[(k*nny+j  )*nnx+i+1].getTT(nt);
            T1 t4 = nodes[(k*nny+j+1)*nnx+i+1].getTT(nt);

            T1 w1 = (ymin+(j+1)*dy - pt.y)/dy;
            T1 w2 = (pt.y - (ymin+j*dy))/dy;

            t1 = t1*w1 + t2*w2;
            t2 = t3*w1 + t4*w2;

            w1 = (xmin+(i+1)*dx - pt.x)/dx;
            w2 = (pt.x - (xmin+i*dx))/dx;

            tt = t1*w1 + t2*w2;

        } else {
            T1 t1 = nodes[(    k*nny+j  )*nnx+i  ].getTT(nt);
            T1 t2 = nodes[((k+1)*nny+j  )*nnx+i  ].getTT(nt);
            T1 t3 = nodes[(    k*nny+j+1)*nnx+i  ].getTT(nt);
            T1 t4 = nodes[((k+1)*nny+j+1)*nnx+i  ].getTT(nt);
            T1 t5 = nodes[(    k*nny+j  )*nnx+i+1].getTT(nt);
            T1 t6 = nodes[((k+1)*nny+j  )*nnx+i+1].getTT(nt);
            T1 t7 = nodes[(    k*nny+j+1)*nnx+i+1].getTT(nt);
            T1 t8 = nodes[((k+1)*nny+j+1)*nnx+i+1].getTT(nt);

            T1 w1 = (zmin+(k+1)*dz - pt.z)/dz;
            T1 w2 = (pt.z - (zmin+k*dz))/dz;

            t1 = t1*w1 + t2*w2;
            t2 = t3*w1 + t4*w2;
            t3 = t5*w1 + t6*w2;
            t4 = t7*w1 + t8*w2;

            w1 = (ymin+(j+1)*dy - pt.y)/dy;
            w2 = (pt.y - (ymin+j*dy))/dy;

            t1 = t1*w1 + t2*w2;
            t2 = t3*w1 + t4*w2;

            w1 = (xmin+(i+1)*dx - pt.x)/dx;
            w2 = (pt.x - (xmin+i*dx))/dx;

            tt = t1*w1 + t2*w2;

        }

        return tt;
    }

    template<typename T1, typename T2, typename NODE, typename CELL>
    void Grid3Drc<T1,T2,NODE,CELL>::grad(sxyz<T1>& g, const sxyz<T1> &pt,
                                         const size_t nt) const {

        // compute travel time gradient at point pt

        T1 p1 = pt.x - dx/2.0;
        T1 p2 = p1 + dx;
        if ( p1 < xmin ) {  // check if on grid edge or out of grid
            p1 = xmin;  // shift pt to allow interpolating in getTraveltime
            p2 = p1 + dx;
        } else if ( p2 > xmax ) {
            p2 = xmax;
            p1 = p2 - dx;
        }
        g.x = (getTraveltime({p2, pt.y, pt.z}, nt) - getTraveltime({p1, pt.y, pt.z}, nt)) / dx;

        p1 = pt.y - dy/2.0;
        p2 = p1 + dy;
        if ( p1 < ymin ) {
            p1 = ymin;
            p2 = p1 + dy;
        } else if ( p2 > ymax ) {
            p2 = ymax;
            p1 = p2 - dy;
        }
        g.y = (getTraveltime({pt.x, p2, pt.z}, nt) - getTraveltime({pt.x, p1, pt.z}, nt)) / dy;

        p1 = pt.z - dz/2.0;
        p2 = p1 + dz;
        if ( p1 < zmin ) {
            p1 = zmin;
            p2 = p1 + dz;
        } else if ( p2 > zmax ) {
            p2 = zmax;
            p1 = p2 - dz;
        }
        g.z = (getTraveltime({pt.x, pt.y, p2}, nt) - getTraveltime({pt.x, pt.y, p1}, nt)) / dz;

    }

    template<typename T1, typename T2, typename NODE, typename CELL>
    T1 Grid3Drc<T1,T2,NODE,CELL>::getTraveltimeFromRaypath(const std::vector<sxyz<T1>>& Tx,
                                                           const std::vector<T1>& t0,
                                                           const sxyz<T1> &Rx,
                                                           const size_t threadNo) const {
        T1 tt = 0.0;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return tt;
            }
        }

        sxyz<T1> prev_pt( Rx );
        sxyz<T1> curr_pt( Rx );
        // distance between opposite nodes of a voxel
        const T1 maxDist = sqrt( dx*dx + dy*dy + dz*dz );
        sxyz<T1> g;

        bool reachedTx = false;
        while ( reachedTx == false ) {

            grad(g, curr_pt, threadNo);
            g *= -1.0;

            ptrdiff_t i, j, k;
            getIJK(curr_pt, i, j, k);

            // planes we will intersect
            T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
            T1 yp = ymin + dy*(j + (boost::math::sign(g.y)>0.0 ? 1.0 : 0.0));
            T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

            if ( std::abs(xp-curr_pt.x)<small2) {
                xp += dx*boost::math::sign(g.x);
            }
            if ( std::abs(yp-curr_pt.y)<small2) {
                yp += dy*boost::math::sign(g.y);
            }
            if ( std::abs(zp-curr_pt.z)<small2) {
                zp += dz*boost::math::sign(g.z);
            }

            // dist to planes
            T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
            T1 ty = g.y!=0.0 ? (yp - curr_pt.y)/g.y : std::numeric_limits<T1>::max();
            T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

            if ( tx<ty && tx<tz ) { // closer to xp
                curr_pt += tx*g;
                curr_pt.x = xp;     // make sure we don't accumulate rounding errors
            } else if ( ty<tz ) {
                curr_pt += ty*g;
                curr_pt.y = yp;
            } else {
                curr_pt += tz*g;
                curr_pt.z = zp;
            }

            if ( curr_pt.x < xmin || curr_pt.x > xmax ||
                curr_pt.y < ymin || curr_pt.y > ymax ||
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
            sxyz<T1> mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
            T2 cellNo = getCellNo(mid_pt);
            tt += cells.computeDt(prev_pt, curr_pt, cellNo);
            prev_pt = curr_pt;

            // are we close enough to one the Tx nodes ?
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                T1 dist = curr_pt.getDistance( Tx[ns] );
                if ( dist < maxDist ) {

                    g = Tx[ns] - curr_pt;
                    // check if we intersect a plane between curr_pt & Tx

                    getIJK(curr_pt, i, j, k);

                    xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
                    yp = ymin + dy*(j + (boost::math::sign(g.y)>0.0 ? 1.0 : 0.0));
                    zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

                    if ( std::abs(xp-curr_pt.x)<small2) {
                        xp += dx*boost::math::sign(g.x);
                    }
                    if ( std::abs(yp-curr_pt.y)<small2) {
                        yp += dy*boost::math::sign(g.y);
                    }
                    if ( std::abs(zp-curr_pt.z)<small2) {
                        zp += dz*boost::math::sign(g.z);
                    }

                    // dist to planes
                    tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
                    ty = g.y!=0.0 ? (yp - curr_pt.y)/g.y : std::numeric_limits<T1>::max();
                    tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

                    if ( tx<ty && tx<tz ) { // closer to xp
                        curr_pt += tx*g;
                        curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                    } else if ( ty<tz ) {
                        curr_pt += ty*g;
                        curr_pt.y = yp;
                    } else {
                        curr_pt += tz*g;
                        curr_pt.z = zp;
                    }

                    if ( curr_pt.getDistance(prev_pt) > dist ||  // we do not intersect a plane
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

    template<typename T1, typename T2, typename NODE, typename CELL>
    void Grid3Drc<T1,T2,NODE,CELL>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                               const std::vector<T1>& t0,
                                               const sxyz<T1> &Rx,
                                               std::vector<sxyz<T1>> &r_data,
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

        sxyz<T1> curr_pt( Rx );
        // distance between opposite nodes of a voxel
        const T1 maxDist = sqrt( dx*dx + dy*dy + dz*dz );
        sxyz<T1> g;
#ifdef DEBUG_RP
        std::cout << "Starting at " << curr_pt << '\n';
#endif

        bool reachedTx = false;
        while ( reachedTx == false ) {

            grad(g, curr_pt, threadNo);
            g *= -1.0;

            ptrdiff_t i, j, k;
            getIJK(curr_pt, i, j, k);

            // planes we will intersect
            T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
            T1 yp = ymin + dy*(j + (boost::math::sign(g.y)>0.0 ? 1.0 : 0.0));
            T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

            if ( std::abs(xp-curr_pt.x)<small2) {
                xp += dx*boost::math::sign(g.x);
            }
            if ( std::abs(yp-curr_pt.y)<small2) {
                yp += dy*boost::math::sign(g.y);
            }
            if ( std::abs(zp-curr_pt.z)<small2) {
                zp += dz*boost::math::sign(g.z);
            }

            // dist to planes
            T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
            T1 ty = g.y!=0.0 ? (yp - curr_pt.y)/g.y : std::numeric_limits<T1>::max();
            T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

            if ( tx<ty && tx<tz ) { // closer to xp
                curr_pt += tx*g;
                curr_pt.x = xp;     // make sure we don't accumulate rounding errors
            } else if ( ty<tz ) {
                curr_pt += ty*g;
                curr_pt.y = yp;
            } else {
                curr_pt += tz*g;
                curr_pt.z = zp;
            }
#ifdef DEBUG_RP
            std::cout << "Grad: " << g << "\t going to: " << curr_pt << '\n';
#endif

            if ( curr_pt.x < xmin || curr_pt.x > xmax ||
                curr_pt.y < ymin || curr_pt.y > ymax ||
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
            sxyz<T1> mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
            T2 cellNo = getCellNo(mid_pt);
            tt += cells.computeDt(r_data.back(), curr_pt, cellNo);
            r_data.push_back( curr_pt );

            // are we close enough to one the Tx nodes ?
            for ( size_t ns=0; ns<Tx.size(); ++ns ) {
                T1 dist = curr_pt.getDistance( Tx[ns] );
                if ( dist < maxDist ) {

                    g = Tx[ns] - curr_pt;
                    // check if we intersect a plane between curr_pt & Tx

                    getIJK(curr_pt, i, j, k);

                    xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
                    yp = ymin + dy*(j + (boost::math::sign(g.y)>0.0 ? 1.0 : 0.0));
                    zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

                    if ( std::abs(xp-curr_pt.x)<small2) {
                        xp += dx*boost::math::sign(g.x);
                    }
                    if ( std::abs(yp-curr_pt.y)<small2) {
                        yp += dy*boost::math::sign(g.y);
                    }
                    if ( std::abs(zp-curr_pt.z)<small2) {
                        zp += dz*boost::math::sign(g.z);
                    }

                    // dist to planes
                    tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
                    ty = g.y!=0.0 ? (yp - curr_pt.y)/g.y : std::numeric_limits<T1>::max();
                    tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

                    if ( tx<ty && tx<tz ) { // closer to xp
                        curr_pt += tx*g;
                        curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                    } else if ( ty<tz ) {
                        curr_pt += ty*g;
                        curr_pt.y = yp;
                    } else {
                        curr_pt += tz*g;
                        curr_pt.z = zp;
                    }

                    if ( curr_pt.getDistance(r_data.back()) > dist ||  // we do not intersect a plane
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

    template<typename T1, typename T2, typename NODE, typename CELL>
    void Grid3Drc<T1,T2,NODE,CELL>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                               const std::vector<T1>& t0,
                                               const sxyz<T1> &Rx,
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

        sxyz<T1> curr_pt( Rx ), prev_pt( Rx );
        // distance between opposite nodes of a voxel
        const T1 maxDist = sqrt( dx*dx + dy*dy + dz*dz );
        sxyz<T1> g;
#ifdef DEBUG_RP
        std::cout << "Starting at " << curr_pt << '\n';
#endif

        siv<T1> cell;
        bool reachedTx = false;
        while ( reachedTx == false ) {

            grad(g, curr_pt, threadNo);
            g *= -1.0;

            ptrdiff_t i, j, k;
            getIJK(curr_pt, i, j, k);

            // planes we will intersect
            T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
            T1 yp = ymin + dy*(j + (boost::math::sign(g.y)>0.0 ? 1.0 : 0.0));
            T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

            if ( std::abs(xp-curr_pt.x)<small2) {
                xp += dx*boost::math::sign(g.x);
            }
            if ( std::abs(yp-curr_pt.y)<small2) {
                yp += dy*boost::math::sign(g.y);
            }
            if ( std::abs(zp-curr_pt.z)<small2) {
                zp += dz*boost::math::sign(g.z);
            }

            // dist to planes
            T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
            T1 ty = g.y!=0.0 ? (yp - curr_pt.y)/g.y : std::numeric_limits<T1>::max();
            T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

            if ( tx<ty && tx<tz ) { // closer to xp
                curr_pt += tx*g;
                curr_pt.x = xp;     // make sure we don't accumulate rounding errors
            } else if ( ty<tz ) {
                curr_pt += ty*g;
                curr_pt.y = yp;
            } else {
                curr_pt += tz*g;
                curr_pt.z = zp;
            }
#ifdef DEBUG_RP
            std::cout << "Grad: " << g << "\t going to: " << curr_pt << '\n';
#endif

            if ( curr_pt.x < xmin || curr_pt.x > xmax ||
                curr_pt.y < ymin || curr_pt.y > ymax ||
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
            sxyz<T1> mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
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
                    // check if we intersect a plane between curr_pt & Tx

                    getIJK(curr_pt, i, j, k);

                    xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
                    yp = ymin + dy*(j + (boost::math::sign(g.y)>0.0 ? 1.0 : 0.0));
                    zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

                    if ( std::abs(xp-curr_pt.x)<small2) {
                        xp += dx*boost::math::sign(g.x);
                    }
                    if ( std::abs(yp-curr_pt.y)<small2) {
                        yp += dy*boost::math::sign(g.y);
                    }
                    if ( std::abs(zp-curr_pt.z)<small2) {
                        zp += dz*boost::math::sign(g.z);
                    }

                    // dist to planes
                    tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
                    ty = g.y!=0.0 ? (yp - curr_pt.y)/g.y : std::numeric_limits<T1>::max();
                    tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

                    if ( tx<ty && tx<tz ) { // closer to xp
                        curr_pt += tx*g;
                        curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                    } else if ( ty<tz ) {
                        curr_pt += ty*g;
                        curr_pt.y = yp;
                    } else {
                        curr_pt += tz*g;
                        curr_pt.z = zp;
                    }

                    if ( curr_pt.getDistance(prev_pt) > dist ||  // we do not intersect a plane
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

                        // to Tx
                        cell.i = getCellNo(Tx[ns]);
                        cell.v = Tx[ns].getDistance(curr_pt);
                        l_data.push_back(cell);
                        tt += cells.computeDt(Tx[ns], curr_pt, cell.i);
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
    }

    template<typename T1, typename T2, typename NODE, typename CELL>
    void Grid3Drc<T1,T2,NODE,CELL>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                               const std::vector<T1>& t0,
                                               const sxyz<T1> &Rx,
                                               std::vector<sxyz<T1>> &r_data,
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

        sxyz<T1> curr_pt( Rx );
        // distance between opposite nodes of a voxel
        const T1 maxDist = sqrt( dx*dx + dy*dy + dz*dz );
        sxyz<T1> g;
#ifdef DEBUG_RP
        std::cout << "Starting at " << curr_pt << '\n';
#endif

        siv<T1> cell;
        bool reachedTx = false;
        while ( reachedTx == false ) {

            grad(g, curr_pt, threadNo);
            g *= -1.0;

            ptrdiff_t i, j, k;
            getIJK(curr_pt, i, j, k);

            // planes we will intersect
            T1 xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
            T1 yp = ymin + dy*(j + (boost::math::sign(g.y)>0.0 ? 1.0 : 0.0));
            T1 zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

            if ( std::abs(xp-curr_pt.x)<small2) {
                xp += dx*boost::math::sign(g.x);
            }
            if ( std::abs(yp-curr_pt.y)<small2) {
                yp += dy*boost::math::sign(g.y);
            }
            if ( std::abs(zp-curr_pt.z)<small2) {
                zp += dz*boost::math::sign(g.z);
            }

            // dist to planes
            T1 tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
            T1 ty = g.y!=0.0 ? (yp - curr_pt.y)/g.y : std::numeric_limits<T1>::max();
            T1 tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

            if ( tx<ty && tx<tz ) { // closer to xp
                curr_pt += tx*g;
                curr_pt.x = xp;     // make sure we don't accumulate rounding errors
            } else if ( ty<tz ) {
                curr_pt += ty*g;
                curr_pt.y = yp;
            } else {
                curr_pt += tz*g;
                curr_pt.z = zp;
            }
#ifdef DEBUG_RP
            std::cout << "Grad: " << g << "\t going to: " << curr_pt << '\n';
#endif

            if ( curr_pt.x < xmin || curr_pt.x > xmax ||
                curr_pt.y < ymin || curr_pt.y > ymax ||
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
            sxyz<T1> mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
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
                    // check if we intersect a plane between curr_pt & Tx

                    getIJK(curr_pt, i, j, k);

                    xp = xmin + dx*(i + (boost::math::sign(g.x)>0.0 ? 1.0 : 0.0));
                    yp = ymin + dy*(j + (boost::math::sign(g.y)>0.0 ? 1.0 : 0.0));
                    zp = zmin + dz*(k + (boost::math::sign(g.z)>0.0 ? 1.0 : 0.0));

                    if ( std::abs(xp-curr_pt.x)<small2) {
                        xp += dx*boost::math::sign(g.x);
                    }
                    if ( std::abs(yp-curr_pt.y)<small2) {
                        yp += dy*boost::math::sign(g.y);
                    }
                    if ( std::abs(zp-curr_pt.z)<small2) {
                        zp += dz*boost::math::sign(g.z);
                    }

                    // dist to planes
                    tx = g.x!=0.0 ? (xp - curr_pt.x)/g.x : std::numeric_limits<T1>::max();
                    ty = g.y!=0.0 ? (yp - curr_pt.y)/g.y : std::numeric_limits<T1>::max();
                    tz = g.z!=0.0 ? (zp - curr_pt.z)/g.z : std::numeric_limits<T1>::max();

                    if ( tx<ty && tx<tz ) { // closer to xp
                        curr_pt += tx*g;
                        curr_pt.x = xp;     // make sure we don't accumulate rounding errors
                    } else if ( ty<tz ) {
                        curr_pt += ty*g;
                        curr_pt.y = yp;
                    } else {
                        curr_pt += tz*g;
                        curr_pt.z = zp;
                    }

                    if ( curr_pt.getDistance(r_data.back()) > dist ||  // we do not intersect a plane
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
                        cell.i = getCellNo(Tx[ns]);
                        cell.v = Tx[ns].getDistance(r_data.back());
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



    template<typename T1, typename T2, typename NODE, typename CELL>
    void Grid3Drc<T1,T2,NODE,CELL>::saveTT(const std::string &fname,
                                           const int all,
                                           const size_t nt,
                                           const int format) const {

        if ( format == 1 ) {
            std::string filename = fname+".dat";
            std::ofstream fout(filename.c_str());
            fout.precision(12);
            for ( T2 n=0; n<nodes.size(); ++n ) {
                if ( nodes[n].isPrimary() || all==1 ) {
                    fout << nodes[n].getX() << '\t'
                    << nodes[n].getY() << '\t'
                    << nodes[n].getZ() << '\t'
                    << nodes[n].getTT(nt) << '\n';
                }
            }
            fout.close();
        } else if ( format == 2 ) {
#ifdef VTK

            std::string filename = fname+".vtr";
            int nn[3] = {static_cast<int>(ncx+1), static_cast<int>(ncy+1), static_cast<int>(ncz+1)};

            vtkSmartPointer<vtkDoubleArray> xCoords = vtkSmartPointer<vtkDoubleArray>::New();
            for (size_t n=0; n<nn[0]; ++n) {
                xCoords->InsertNextValue( xmin + n*dx );
            }
            vtkSmartPointer<vtkDoubleArray> yCoords = vtkSmartPointer<vtkDoubleArray>::New();
            for (size_t n=0; n<nn[1]; ++n) {
                yCoords->InsertNextValue( ymin + n*dy );
            }
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
                    double x[3] = {nodes[n].getX(), nodes[n].getY(), nodes[n].getZ()};
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
                if ( nodes[n].isPrimary() || all==1 ) {
                    T1 tmp[] = { nodes[n].getX(), nodes[n].getY(), nodes[n].getZ(), nodes[n].getTT(nt) };
                    fout.write( (char*)tmp, 4*sizeof(T1) );
                }
            }
            fout.close();
        } else {
            throw std::runtime_error("Unsupported format for saving traveltimes");
        }
    }

}

#endif

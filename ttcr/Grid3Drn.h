//
//  Grid3Drn.h
//  ttcr
//
//  Created by Giroux Bernard on 12-08-15.
//  Copyright (c) 2012 INRS-ETE. All rights reserved.
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
 * @file Grid3Drn.h
 * @brief Base class for 3-D rectilinear grids with node-based slowness.
 *
 * Declares ttcr::Grid3Drn, the mid-level base for every 3-D rectilinear solver
 * whose slowness lives at the **nodes**: Grid3Drnsp (shortest path), Grid3Drnfs
 * (fast sweeping), Grid3Drndsp (dynamic shortest path), Grid3Drnfs_OpenCL —
 * and also Grid3Drcfs and Grid3Drcfs_OpenCL, which despite their @c rc names
 * derive from here and merely average a cell model onto the nodes first.
 *
 * It is the 3-D counterpart of ttcr::Grid2Drn and the node-slowness sibling of
 * ttcr::Grid3Drc. As in 2-D, a traveltime increment is the mean of the two
 * endpoint slownesses times the distance — trapezoidal integration of a field
 * that varies linearly along the segment — rather than a value drawn from a
 * @c CELL policy. The family therefore has no @c CELL template parameter and is
 * isotropic.
 *
 * @section g3drn_procvel Velocity versus slowness interpolation
 * Unlike the 2-D family, these grids carry a @c processVel flag
 * (ttcr::input_parameters::processVel). It changes what is interpolated between
 * nodes:
 *
 * - **false** — slowness is interpolated directly, @f$\sum_i w_i s_i@f$;
 * - **true** — velocity is interpolated and the result inverted,
 *   @f$1 / \sum_i w_i (1/s_i)@f$.
 *
 * The two agree only where slowness is constant; elsewhere they give different
 * fields, so the flag is a modelling choice rather than an implementation
 * detail. It is honoured throughout
 * ttcr::Grid3Drn::computeSlowness. @sa Interpolator.h, whose @c ...Vel variants
 * implement the same distinction.
 *
 * @section g3drn_numbering Node numbering
 * Nodes are numbered **x-fastest**, matching ttcr::Grid3Drc
 * (@ref g3drc_numbering) and opposite to the 2-D families. Primary nodes — the
 * @f$(n_{cx}+1)(n_{cy}+1)(n_{cz}+1)@f$ cell corners — come first, followed by
 * any secondary nodes a solver adds.
 *
 * @section g3drn_sweeps Fast sweeping stencils
 * Two update stencils are provided, @c sweep and @c sweep_weno3, selected by
 * ttcr::input_parameters::weno3. The 2-D base additionally offers rotated and
 * mixed x-z stencils; there are no 3-D equivalents, so
 * ttcr::input_parameters::rotated_template has no effect here.
 *
 * @sa Grid3D.h, Grid3Drc.h, Grid2Drn.h, Grid3Drnsp.h, Grid3Drnfs.h
 */

#ifndef ttcr_Grid3Drn_h
#define ttcr_Grid3Drn_h

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <memory>
#include <mutex>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <ctime>

#include <boost/math/special_functions/sign.hpp>

#ifdef VTK
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"
#include "vtkXMLRectilinearGridWriter.h"
#endif

#include "Grid3D.h"
#include "Interpolator.h"
#include "NodeKDTree.h"

namespace ttcr {

    template<typename T1, typename T2, typename NODE>
    /**
     * @brief 3-D rectilinear grid holding one slowness value per node.
     *
     * @tparam T1   floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2   integer type of node and cell indices, normally @c uint32_t.
     * @tparam NODE node type, e.g. ttcr::Node3Dn or ttcr::Node3Dnsp. It must
     *              expose @c getNodeSlowness.
     *
     * @note No @c CELL parameter, unlike ttcr::Grid3Drc: a node carries a single
     *       scalar slowness, so this family is isotropic.
     * @note Abstract in practice — it implements everything except @c raytrace.
     */
    class Grid3Drn : public Grid3D<T1,T2> {
    public:

        /**
         * @brief Build the grid geometry and allocate its nodes.
         *
         * @param nx      number of cells along x.
         * @param ny      number of cells along y.
         * @param nz      number of cells along z.
         * @param ddx     cell size along x.
         * @param ddy     cell size along y.
         * @param ddz     cell size along z.
         * @param minx    x coordinate of the grid origin.
         * @param miny    y coordinate of the grid origin.
         * @param minz    z coordinate of the grid origin.
         * @param ttrp    recompute receiver traveltimes by integrating slowness
         *                along the traced raypath.
         * @param procVel interpolate velocity rather than slowness.
         *                @sa @ref g3drn_procvel
         * @param nt      number of threads; sizes each node's traveltime array.
         * @param _translateOrigin shift the grid so its origin is at (0,0,0).
         *
         * @post The grid holds @f$(n_x+1)(n_y+1)(n_z+1)@f$ primary nodes. No
         *       cells are allocated — this family stores no per-cell data.
         *       Slowness is **not** set.
         * @warning The cell count passed to the base is computed in @c T2
         *          arithmetic (@c nx*ny*nz), so a grid large enough to overflow
         *          a 32-bit index wraps silently — the same caveat as
         *          ttcr::Grid3Drc.
         */
        /* Constructor Format:
         Grid3Drn<T1,T2>::Grid3Drn(nb cells in x, nb cells in y, nb cells in z,
         x cells size, y cells size, z cells size,
         x origin, y origin, z origin,
         index of the thread)
         */
        Grid3Drn(const T2 nx, const T2 ny, const T2 nz,
                 const T1 ddx, const T1 ddy, const T1 ddz,
                 const T1 minx, const T1 miny, const T1 minz,
                 const bool ttrp, const bool procVel, const size_t nt=1,
                 const bool _translateOrigin=false) :
        Grid3D<T1,T2>(ttrp, nx*ny*nz, nt, _translateOrigin),
        dx(ddx), dy(ddy), dz(ddz),
        xmin(minx), ymin(miny), zmin(minz),
        xmax(minx+nx*ddx), ymax(miny+ny*ddy), zmax(minz+nz*ddz),
        ncx(nx), ncy(ny), ncz(nz), processVel(procVel),
        nodes(std::vector<NODE>((nx+1)*(ny+1)*(nz+1), NODE(nt)))
        { }

        /// Destructor.
        virtual ~Grid3Drn() {}

        /**
         * @brief Set the slowness at every node.
         * @param s one slowness per node — **all** nodes, secondary ones
         *          included, so its size must equal @ref getNumberOfNodes.
         * @throws std::length_error if the sizes do not match.
         * @note Virtual, and overridden by the solvers that add secondary nodes
         *       (ttcr::Grid3Drnsp, ttcr::Grid3Drndsp) with a version taking one
         *       value per *primary* node and interpolating onto the rest, and by
         *       ttcr::Grid3Drcfs with one taking cell values.
         * @warning Not symmetric with @ref getSlowness, which returns
         *          primary-node values only.
         */
        virtual void setSlowness(const std::vector<T1>& s) {
            if ( nodes.size() != s.size() ) {
                throw std::length_error("Error: slowness vectors of incompatible size.");
            }
            for ( size_t n=0; n<nodes.size(); ++n ) {
                nodes[n].setNodeSlowness( s[n] );
            }
        }
        /**
         * @brief Copy the primary-node slowness values out of the grid.
         * @param[out] slowness resized to @f$(n_{cx}+1)(n_{cy}+1)(n_{cz}+1)@f$
         *                      and filled from the leading primary nodes.
         */
        void getSlowness(std::vector<T1>& slowness) const {
            if (slowness.size() != static_cast<size_t>(ncx+1) * (ncy+1) * (ncz+1)) {
                slowness.resize(static_cast<size_t>(ncx+1) * (ncy+1) * (ncz+1));
            }
            for ( size_t n=0; n<slowness.size(); ++n ) {
                slowness[n] = nodes[n].getNodeSlowness();
            }
        }

        /// @return Total number of nodes, primary and secondary.
        size_t getNumberOfNodes() const { return nodes.size(); }
        /**
         * @brief Number of cells the grid geometry defines.
         * @return The product @f$n_{cx}n_{cy}n_{cz}@f$.
         * @note Geometric only — no per-cell data is stored by this family.
         */
        size_t getNumberOfCells() const { return static_cast<size_t>(ncx)*ncy*ncz; }

        /**
         * @brief Collect the traveltimes computed at the primary nodes.
         * @param[out] tt       resized to the primary-node count and filled.
         * @param[in]  threadNo thread whose solution to read.
         * @note Reads the leading prefix of the node vector, relying on primary
         *       nodes coming first.
         */
        void getTT(std::vector<T1>& tt, const size_t threadNo=0) const final {
            size_t nPrimary = static_cast<size_t>(ncx+1) * (ncy+1) * (ncz+1);
            tt.resize(nPrimary);
            for ( size_t n=0; n<nPrimary; ++n ) {
                tt[n] = nodes[n].getTT(threadNo);
            }
        }

        /**
         * @brief Memory occupied by the neighbour lists, in bytes.
         * @return Total payload of the per-cell neighbour vectors; the
         *         per-vector bookkeeping is not counted, so this understates the
         *         real footprint.
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
         *       the non-shortest-path node types.
         */
        size_t getNodesSize() const {
            size_t size = 0;
            for ( size_t n=0; n<nodes.size(); ++n ) {
                size += nodes[n].getSize();
            }
            return size;
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
         * @brief Read a traveltime field back from a file.
         * @param fname  input filename.
         * @param all    whether the file includes the secondary nodes.
         * @param nt     thread whose slot to populate.
         * @param format format the file was written in, matching @ref saveTT.
         * @note The counterpart of @ref saveTT, letting a solved field be reused
         *       without recomputing it. No 2-D equivalent exists.
         */
        void loadTT(const std::string &fname, const int all, const size_t nt=0,
                    const int format=1) const;

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
         * @note Writes nothing when the solver adds no secondary nodes.
         */
        void dump_secondary(std::ofstream& os) const {
            size_t nPrimary = static_cast<size_t>(ncx+1) * (ncy+1) * (ncz+1);
            for ( size_t n=nPrimary; n<nodes.size(); ++n ) {
                os << nodes[n].getX() << ' ' << nodes[n].getY() << ' ' << nodes[n].getZ() << '\n';
            }
        }

        /**
         * @brief Slowness at an arbitrary point, interpolated from the nodes.
         * @param pt           point to sample. Taken by value, since it may be
         *                     translated internally.
         * @param isTranslated true if @p pt is already in the translated frame.
         * @return Interpolated slowness.
         * @note Continuous, unlike ttcr::Grid3Drc::computeSlowness which is
         *       piecewise constant. Whether slowness or velocity is the
         *       interpolated quantity depends on @ref processVel —
         *       @sa @ref g3drn_procvel
         */
        T1 computeSlowness(sxyz<T1> pt, const bool isTranslated=false) const;

#ifdef VTK
        /**
         * @brief Write the model as a VTK rectilinear grid.
         * @param fname        output filename.
         * @param saveSlowness true to write slowness, false to write velocity.
         * @note Available only when the library is built with @c VTK defined.
         */
        void saveModelVTR(const std::string &fname,
                          const bool saveSlowness=true) const;
#endif

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

        /// Interpolate velocity rather than slowness. @sa @ref g3drn_procvel
        bool processVel;

        /// Primary nodes first, then any secondary nodes the derived solver
        /// added. @c mutable because @c const raytracing methods update the
        /// traveltimes stored in them.
        mutable std::vector<NODE> nodes;

        // kd-tree over all nodes (primary + secondary), to replace the
        // brute-force scan that checks whether a point coincides with a node.
        // Built lazily on first use; node coordinates are fixed after
        // construction and queries are read-only, hence thread-safe.
        mutable std::unique_ptr<NodeKDTree<T1,T2>> kdtree;
        mutable std::once_flag kdtreeFlag;

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
         * @note The defaults of 0 give a primary-only grid, which is what the
         *       sweeping solvers want.
         */
        void buildGridNodes(const T2 nsnx=0, const T2 nsny=0, const T2 nsnz=0);

        /**
         * @brief Fill in the slowness of the secondary nodes.
         * @pre The primary nodes already carry their slowness values.
         * @note Secondary nodes lie along cell edges, so their slowness is
         *       interpolated from the primary nodes at the ends of that edge —
         *       honouring @ref processVel. @sa @ref g3drn_procvel
         */
        void interpSecondary();

        /**
         * @brief Locate the cell containing a point.
         * @param pt point to locate.
         * @return Cell number, x varying fastest (@ref g3drn_numbering).
         * @note Used to index a subclass's cell model; this family stores
         *       nothing per cell itself.
         * @warning A point outside the grid is not detected.
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
         * @return Cell number, x varying fastest.
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
         * @sa @ref g3drn_numbering — this inverts @ref getCellNo.
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
         * @note For callers that step off the grid and need a negative index
         *       rather than a wrap.
         */
        void getIJK(const sxyz<T1>& pt, ptrdiff_t& i, ptrdiff_t& j, ptrdiff_t& k) const {
            i = static_cast<ptrdiff_t>( small2 + (pt.x-xmin)/dx );
            j = static_cast<ptrdiff_t>( small2 + (pt.y-ymin)/dy );
            k = static_cast<ptrdiff_t>( small2 + (pt.z-zmin)/dz );
        }

        /**
         * @brief Verify that every point lies inside the grid.
         * @param pts        points to check. Taken **by value**, since they may
         *                   be translated.
         * @param translated true if @p pts are already in the translated frame.
         * @throws std::runtime_error naming the offending point.
         */
        void checkPts(std::vector<sxyz<T1>> pts, const bool translated=false) const;

        /**
         * @brief Traveltime increment between two nodes.
         * @param source node the ray comes from.
         * @param node   node it reaches.
         * @return The increment @f$\tfrac{1}{2}(s_{source}+s_{node})\,\ell@f$ — trapezoidal
         *         integration along the segment.
         * @note Independent of @ref processVel — the increment always averages
         *       slowness, whatever quantity @ref computeSlowness interpolates.
         */
        T1 computeDt(const NODE& source, const NODE& node) const {
            return (node.getNodeSlowness()+source.getNodeSlowness())/2. * source.getDistance( node );
        }

        /**
         * @brief Traveltime increment from a node to an arbitrary point.
         * @param source node the ray comes from.
         * @param node   point it reaches.
         * @param slo    slowness at that point, which the caller must supply —
         *               a plain point carries none.
         * @return The increment @f$\tfrac{1}{2}(s_{source}+slo)\,\ell@f$.
         */
        T1 computeDt(const NODE& source, const sxyz<T1>& node, T1 slo) const {
            return (slo+source.getNodeSlowness())/2. * source.getDistance( node );
        }

        /**
         * @brief Test whether a value is close to a whole number.
         * @param value value to test.
         * @return True if it lies within @ref small of an integer.
         * @note Used to decide whether a point falls exactly on a grid line, so
         *       interpolation can be reduced by one dimension.
         * @warning Compares the **signed** remainder against @ref small, so a
         *          value just below an integer gives a negative remainder and
         *          passes trivially. Only values just above are really tested.
         */
        bool isNearInt( double value ) const {
            return ( remainder(value, 1.)  <= small );
        }

        /**
         * @brief Traveltime at a point.
         * @param pt point to evaluate.
         * @param nt thread whose solution to read.
         * @return Traveltime, interpolated if @p pt does not sit on a node.
         */
        T1 getTraveltime(const sxyz<T1> &pt, const size_t nt) const;


        /**
         * @brief Traveltime at a receiver, also reporting where it came from.
         * @param[in]  Rx           receiver position.
         * @param[out] nodeParentRx index of the node the ray arrived from.
         * @param[out] cellParentRx index of the cell it crossed to get there.
         * @param[in]  threadNo     thread whose solution to read.
         * @return Traveltime at @p Rx.
         * @note The two outputs are the starting point for walking a raypath
         *       back to the source. They are unnamed in this declaration; the
         *       names used here are those of the out-of-class definition.
         */
        T1 getTraveltime(const sxyz<T1>& Rx,
                         T2& nodeParentRx, T2& cellParentRx,
                         const size_t threadNo) const;


        /**
         * @brief Traveltime gradient at a node, by grid index.
         * @param[out] g       gradient vector.
         * @param[in]  i       x index of the node.
         * @param[in]  j       y index of the node.
         * @param[in]  k       z index of the node.
         * @param[in]  nt      thread whose traveltime field to differentiate.
         */
        void grad(sxyz<T1>& g, const size_t i, const size_t j, const size_t k,
                  const size_t nt) const;

        /**
         * @brief Second-order traveltime gradient at an arbitrary point.
         * @param[out] g  gradient vector.
         * @param[in]  pt point at which to evaluate it.
         * @param[in]  nt thread whose traveltime field to differentiate.
         * @note More accurate but wider-stencilled than @ref grad; which one is
         *       used depends on ttcr::input_parameters::raypath_method.
         *       @sa Grad.h
         */
        void gradO2(sxyz<T1>& g, const sxyz<T1> &pt, const size_t nt) const;
        /**
         * @brief Traveltime gradient at an arbitrary point.
         * @param[out] g  gradient vector.
         * @param[in]  pt point at which to evaluate it.
         * @param[in]  nt thread whose traveltime field to differentiate.
         * @note The raypath tracers step down this gradient.
         */
        void grad(sxyz<T1>& g, const sxyz<T1> &pt, const size_t nt) const;

        /**
         * @brief Traveltime obtained by integrating slowness along the raypath.
         * @param Tx       source positions.
         * @param t0       origin time of each source.
         * @param Rx       receiver position.
         * @param threadNo thread whose solution to use.
         * @return Traveltime at @p Rx.
         * @note Selected by the @c ttrp constructor flag; follows the path
         *       rather than smoothing across cells, at the cost of a raypath
         *       trace per receiver.
         */
        T1 getTraveltimeFromRaypath(const std::vector<sxyz<T1>>& Tx,
                                    const std::vector<T1>& t0,
                                    const sxyz<T1> &Rx,
                                    const size_t threadNo) const;

        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const sxyz<T1> &Rx,
                        std::vector<sxyz<T1>> &r_data,
                        const size_t threadNo) const;

        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const sxyz<T1> &Rx,
                        std::vector<sxyz<T1>> &r_data,
                        T1 &tt,
                        const size_t threadNo) const;

        /**
         * @brief Trace the raypath, recording node-level sensitivities.
         * @param[in]  Tx       source positions.
         * @param[in]  t0       origin time of each source.
         * @param[in]  Rx       receiver position.
         * @param[out] m_data   sensitivity entries as (row, node, value)
         *                      triples. @sa sijv
         * @param[in]  RxNo     receiver number, used as the matrix row index.
         * @param[out] tt       traveltime along the path.
         * @param[in]  threadNo thread whose solution to use.
         * @note @ref sijv rather than @ref siv because slowness lives at the
         *       nodes here: a segment's sensitivity spreads over the nodes it
         *       interpolates between, so an explicit row index is needed.
         */
        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const sxyz<T1>& Rx,
                        std::vector<sijv<T1>>& m_data,
                        const size_t RxNo,
                        T1 &tt,
                        const size_t threadNo) const;

        /**
         * @brief Trace the raypath, returning geometry and node sensitivities.
         * @param[in]  Tx       source positions.
         * @param[in]  t0       origin time of each source.
         * @param[in]  Rx       receiver position.
         * @param[out] r_data   raypath points.
         * @param[out] m_data   sensitivity entries. @sa sijv
         * @param[in]  RxNo     receiver number, used as the matrix row index.
         * @param[out] tt       traveltime along the path.
         * @param[in]  threadNo thread whose solution to use.
         */
        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const sxyz<T1>& Rx,
                        std::vector<sxyz<T1>>& r_data,
                        std::vector<sijv<T1>>& m_data,
                        const size_t RxNo,
                        T1 &tt,
                        const size_t threadNo) const;

        /**
         * @brief Record the per-cell path lengths without returning the geometry.
         * @param[in]  Tx       source positions.
         * @param[in]  t0       origin time of each source.
         * @param[in]  Rx       receiver position.
         * @param[out] l_data   per-cell segment lengths. @sa siv
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
         * @param[out] l_data   per-cell segment lengths. @sa siv
         * @param[out] tt       traveltime along the path.
         * @param[in]  threadNo thread whose solution to use.
         */
        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const sxyz<T1> &Rx,
                        std::vector<sxyz<T1>> &r_data,
                        std::vector<siv<T1>> &l_data,
                        T1 &tt,
                        const size_t threadNo) const;

        /**
         * @name Fast sweeping passes
         *
         * One Gauss-Seidel pass over the grid in each alternating direction,
         * relaxing every node that is not frozen. Only two stencils exist in
         * 3-D; see @ref g3drn_sweeps.
         *
         * @param frozen   per-node flag; a frozen node holds a source value and
         *                 is never updated.
         * @param threadNo thread whose traveltime field to sweep.
         * @{
         */
        /// Axis-aligned first-order stencil.
        void sweep(const std::vector<bool>& frozen,
                   const size_t threadNo) const;
        /// Third-order WENO stencil.
        void sweep_weno3(const std::vector<bool>& frozen,
                         const size_t threadNo) const;
        /// @}

        /**
         * @name Single-node updates
         *
         * Solve the local eikonal update at one node given its already-relaxed
         * neighbours. The parameters are the node's x, y and z indices and the
         * thread number.
         * @{
         */
        void update_node(const size_t, const size_t, const size_t, const size_t=0) const;
        void update_node_weno3(const size_t, const size_t, const size_t, const size_t=0) const;
        /// @}

        /**
         * @brief Seed the fast sweeping solve around the sources.
         * @param[in]  Tx       source positions.
         * @param[in]  t0       origin time of each source.
         * @param[out] frozen   per-node flag; nodes near a source are given
         *                      analytic traveltimes and frozen.
         * @param[in]  npts     half-width, in nodes, of the frozen region around
         *                      each source.
         * @param[in]  threadNo thread to initialise.
         * @note Sweeping cannot start from a point source directly — the local
         *       update needs neighbours already holding valid times, which this
         *       frozen halo provides.
         */
        void initFSM(const std::vector<sxyz<T1>>& Tx,
                     const std::vector<T1>& t0,
                     std::vector<bool>& frozen,
                     const int npts,
                     const size_t threadNo) const;

    private:
        /// @name Non-copyable
        /// Private, the pre-C++11 idiom. @c operator= returns @c *this without
        /// copying, so it is not undefined behaviour, but it would silently do
        /// nothing if reached.
        /// @{
        Grid3Drn() {}
        Grid3Drn(const Grid3Drn<T1,T2,NODE>& g) {}
        Grid3Drn<T1,T2,NODE>& operator=(const Grid3Drn<T1,T2,NODE>& g) { return *this; }
        /// @}

        /**
         * @brief One-dimensional third-order WENO upwind derivative.
         * @param v0 traveltime at stencil point 0.
         * @param v1 traveltime at stencil point 1.
         * @param v2 traveltime at stencil point 2.
         * @param v3 traveltime at stencil point 3.
         * @param v4 traveltime at stencil point 4.
         * @param dx      node spacing along the direction being differenced.
         * @param forward true for the forward-biased stencil, false for backward.
         * @return The WENO-weighted derivative estimate.
         * @note The nonlinear weights fall away from a stencil that straddles a
         *       kink in the traveltime field, which is what keeps the scheme
         *       third-order in smooth regions without oscillating at wavefront
         *       crossings.
         */
        T1 weno3_upwind(const T1 v0, const T1 v1, const T1 v2, const T1 v3, const T1 v4, const T1 dx, bool forward) const;

    public:
        /**
         * @copydoc Grid3D::getTraveltimeGradient
         * @note Uses the second-order @ref gradO2 rather than the fourth-order
         *       @ref grad.  The wider stencil of the latter reaches across the
         *       strongly curved part of the field near the source, and measured
         *       against an analytic constant-gradient solution it puts the
         *       take-off direction out by 5--9 degrees where @ref gradO2 stays
         *       under 2.
         */
        void getTraveltimeGradient(sxyz<T1>& g, const sxyz<T1>& pt,
                                   const size_t threadNo) const override {
            gradO2(g, pt, threadNo);
        }

        /// Mean node spacing, the length scale used by @ref Grid3D::computeH.
        const T1 getAverageEdgeLength() const override {
            return (dx + dy + dz) / static_cast<T1>(3.0);
        }

    };


    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::buildGridNodes(const T2 nsnx, const T2 nsny, const T2 nsnz) {

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

        // Create the grid, assign a number for each node, determine the type of the node and find the owners
        // Nodes and cells are first indexed in z, then y, and x.
        // Secondary nodes are placed on the faces and edges of every cells.
        // Ex: the node in "node[A]=(i,j,k)" is followed by the node in "node[A+1]=(i+dx,j,k)"

        T2 cell_XmYmZm;     // cell in the (x-,y-,z-) direction from the node
        T2 cell_XpYmZm;     // cell in the (x+,y-,z-) direction from the node
        T2 cell_XmYpZm;
        T2 cell_XpYpZm;
        T2 cell_XmYmZp;
        T2 cell_XpYmZp;
        T2 cell_XmYpZp;
        T2 cell_XpYpZp;

        T2 n=0;
        for ( T2 nk=0; nk<=ncz; ++nk ) {

            T1 z = zmin + nk*dz;

            for ( T2 nj=0; nj<=ncy; ++nj ) {

                T1 y = ymin + nj*dy;

                for ( T2 ni=0; ni<=ncx; ++ni, ++n ){

                    T1 x = xmin + ni*dx;

                    // Find the adjacent cells for each primary node

                    if (ni < ncx && nj < ncy && nk < ncz){
                        cell_XpYpZp = nj*ncx + nk*(ncx*ncy) + ni;
                    }
                    else {
                        cell_XpYpZp = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj < ncy && nk < ncz){
                        cell_XmYpZp = nj*ncx + nk*(ncx*ncy) + ni - 1;
                    }
                    else {
                        cell_XmYpZp = std::numeric_limits<T2>::max();
                    }

                    if (ni < ncx && nj > 0 && nk < ncz){
                        cell_XpYmZp = (nj-1)*ncx + nk*(ncx*ncy) + ni;
                    }
                    else {
                        cell_XpYmZp = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj > 0 && nk < ncz){
                        cell_XmYmZp = (nj-1)*ncx + nk*(ncx * ncy) + ni - 1;
                    }
                    else {
                        cell_XmYmZp = std::numeric_limits<T2>::max();
                    }

                    if (ni < ncx && nj < ncy && nk > 0){
                        cell_XpYpZm = nj*ncx + (nk-1)*(ncx*ncy) + ni;
                    }
                    else {
                        cell_XpYpZm = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj < ncy && nk > 0){
                        cell_XmYpZm = nj*ncx + (nk-1)*(ncx*ncy) + ni - 1;
                    }
                    else {
                        cell_XmYpZm = std::numeric_limits<T2>::max();
                    }

                    if (ni < ncx && nj > 0 && nk > 0){
                        cell_XpYmZm = (nj-1)*ncx + (nk-1)*(ncx*ncy) + ni;
                    }
                    else {
                        cell_XpYmZm = std::numeric_limits<T2>::max();
                    }

                    if (ni > 0 && nj > 0 && nk > 0){
                        cell_XmYmZm = (nj-1)*ncx + (nk-1)*(ncx*ncy) + ni - 1;
                    }
                    else {
                        cell_XmYmZm = std::numeric_limits<T2>::max();
                    }


                    // Index the primary nodes owners

                    if (cell_XmYmZm != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cell_XmYmZm );
                    }
                    if (cell_XpYmZm != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cell_XpYmZm );
                    }
                    if (cell_XmYpZm != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cell_XmYpZm );
                    }
                    if (cell_XpYpZm != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cell_XpYpZm );
                    }
                    if (cell_XmYmZp != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cell_XmYmZp );
                    }
                    if (cell_XpYmZp != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cell_XpYmZp );
                    }
                    if (cell_XmYpZp != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cell_XmYpZp );
                    }
                    if (cell_XpYpZp != std::numeric_limits<T2>::max() ) {
                        nodes[n].pushOwner( cell_XpYpZp );
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
                            cell_XpYpZp = nj*ncx + nk*(ncx*ncy) + ni;
                        }
                        else {
                            cell_XpYpZp = std::numeric_limits<T2>::max();
                        }

                        if (ni > 0 && nj < ncy && nk < ncz){
                            cell_XmYpZp = nj*ncx + nk*(ncx*ncy) + ni - 1;
                        }
                        else {
                            cell_XmYpZp = std::numeric_limits<T2>::max();
                        }

                        if (ni < ncx && nj > 0 && nk < ncz){
                            cell_XpYmZp = (nj-1)*ncx + nk*(ncx*ncy) + ni;
                        }
                        else {
                            cell_XpYmZp = std::numeric_limits<T2>::max();
                        }

                        if (ni > 0 && nj > 0 && nk < ncz){
                            cell_XmYmZp = (nj-1)*ncx + nk*(ncx * ncy) + ni - 1;
                        }
                        else {
                            cell_XmYmZp = std::numeric_limits<T2>::max();
                        }

                        if (ni < ncx && nj < ncy && nk > 0){
                            cell_XpYpZm = nj*ncx + (nk-1)*(ncx*ncy) + ni;
                        }
                        else {
                            cell_XpYpZm = std::numeric_limits<T2>::max();
                        }

                        if (ni > 0 && nj < ncy && nk > 0){
                            cell_XmYpZm = nj*ncx + (nk-1)*(ncx*ncy) + ni - 1;
                        }
                        else {
                            cell_XmYpZm = std::numeric_limits<T2>::max();
                        }

                        if (ni < ncx && nj > 0 && nk > 0){
                            cell_XpYmZm = (nj-1)*ncx + (nk-1)*(ncx*ncy) + ni;
                        }
                        else {
                            cell_XpYmZm = std::numeric_limits<T2>::max();
                        }

                        if (ni > 0 && nj > 0 && nk > 0){
                            cell_XmYmZm = (nj-1)*ncx + (nk-1)*(ncx*ncy) + ni - 1;
                        }
                        else {
                            cell_XmYmZm = std::numeric_limits<T2>::max();
                        }

                        // Secondary nodes on x edge
                        if ( ni < ncx ) {
                            for (T2 ns=0; ns< nsnx; ++ns, ++n ) {

                                T1 xsv = xmin + ni*dx + (ns+1)*dxs;

                                if ( cell_XpYmZm != std::numeric_limits<T2>::max() ) {
                                    nodes[n].pushOwner( cell_XpYmZm );
                                }
                                if ( cell_XpYpZm != std::numeric_limits<T2>::max() ) {
                                    nodes[n].pushOwner( cell_XpYpZm );
                                }
                                if ( cell_XpYmZp != std::numeric_limits<T2>::max() ) {
                                    nodes[n].pushOwner( cell_XpYmZp );
                                }
                                if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
                                    nodes[n].pushOwner( cell_XpYpZp );
                                }
                                nodes[n].setXYZindex( xsv, y, z, n );

                                if (nj >0 && nk>0){
                                    nodes[n].setPrimary(false);
                                }
                                else if (nj==0 && nk>0){
                                    nodes[n].setPrimary(false);
                                }
                                else if (nj>0 && nk==0){
                                    nodes[n].setPrimary(false);
                                }
                                else if (nj==0 && nk==0){
                                    nodes[n].setPrimary(false);
                                }
                            }
                        }

                        // Secondary nodes on y edge
                        if ( nj < ncy ) {
                            for (T2 ns=0; ns< nsny; ++ns, ++n ) {

                                T1 ysv = ymin + nj* dy + (ns+1)*dys;

                                if ( cell_XmYpZm != std::numeric_limits<T2>::max() ) {
                                    nodes[n].pushOwner( cell_XmYpZm );
                                }
                                if ( cell_XpYpZm != std::numeric_limits<T2>::max() ) {
                                    nodes[n].pushOwner( cell_XpYpZm );
                                }
                                if ( cell_XmYpZp != std::numeric_limits<T2>::max() ) {
                                    nodes[n].pushOwner( cell_XmYpZp );
                                }
                                if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
                                    nodes[n].pushOwner( cell_XpYpZp );
                                }
                                nodes[n].setXYZindex( x, ysv, z, n );

                                if (ni >0 && nk>0){
                                    nodes[n].setPrimary(false);
                                }
                                else if (ni>0 && nk==0){
                                    nodes[n].setPrimary(false);
                                }
                                else if (ni==0 && nk>0){
                                    nodes[n].setPrimary(false);
                                }
                                else if (ni==0 && nk==0){
                                    nodes[n].setPrimary(false);
                                }
                            }
                        }

                        // Secondary nodes on z edge
                        if ( nk < ncz ) {
                            for (T2 ns=0; ns< nsnz; ++ns, ++n ) {

                                T1 zsv = zmin + nk* dz + (ns+1)*dzs;

                                if ( cell_XmYmZp != std::numeric_limits<T2>::max() ) {
                                    nodes[n].pushOwner( cell_XmYmZp );
                                }
                                if ( cell_XpYmZp != std::numeric_limits<T2>::max() ) {
                                    nodes[n].pushOwner( cell_XpYmZp );
                                }
                                if ( cell_XmYpZp != std::numeric_limits<T2>::max() ) {
                                    nodes[n].pushOwner( cell_XmYpZp );
                                }
                                if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
                                    nodes[n].pushOwner( cell_XpYpZp );
                                }
                                nodes[n].setXYZindex( x, y, zsv, n );

                                if (ni >0 && nj>0){
                                    nodes[n].setPrimary(false);
                                }
                                else if (ni>0 && nj==0){
                                    nodes[n].setPrimary(false);
                                }
                                else if (ni==0 && nj>0){
                                    nodes[n].setPrimary(false);
                                }
                                else if (ni==0 && nj==0){
                                    nodes[n].setPrimary(false);
                                }
                            }
                        }

                        // Secondary nodes on the xy0 planes
                        if ( ni < ncx && nj < ncy ) {
                            for (T2 sy=0; sy < nsny; ++sy){
                                for (T2 sx=0; sx < nsnx; ++sx, n++) {

                                    T1 ysv = ymin + nj*dy + (sy+1)*dys;
                                    T1 xsv = xmin + ni*dx + (sx+1)*dxs;

                                    if ( cell_XpYpZm != std::numeric_limits<T2>::max() ) {
                                        nodes[n].pushOwner( cell_XpYpZm );
                                    }
                                    if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
                                        nodes[n].pushOwner( cell_XpYpZp );
                                    }
                                    nodes[n].setXYZindex( xsv, ysv, z, n );

                                    if (nk>0){
                                        nodes[n].setPrimary(false);
                                    }
                                    else if (nk==0){
                                        nodes[n].setPrimary(false);
                                    }
                                }
                            }
                        }

                        // Secondary nodes on the x0z planes
                        if ( ni < ncx && nk < ncz ) {
                            for(T2 sz=0; sz < nsnz; ++sz){
                                for(T2 sx=0; sx < nsnx; ++sx, n++){

                                    T1 zsv = zmin + nk*dz + (sz+1)*dzs;
                                    T1 xsv = xmin + ni*dx + (sx+1)*dxs;

                                    if ( cell_XpYmZp != std::numeric_limits<T2>::max() ) {
                                        nodes[n].pushOwner( cell_XpYmZp );
                                    }
                                    if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
                                        nodes[n].pushOwner( cell_XpYpZp );
                                    }
                                    nodes[n].setXYZindex( xsv, y, zsv, n );

                                    if (nj>0){
                                        nodes[n].setPrimary(false);
                                    }
                                    else if (nj==0){
                                        nodes[n].setPrimary(false);
                                    }
                                }
                            }
                        }

                        // Secondary nodes on the 0yz planes
                        if ( nj < ncy && nk < ncz ) {
                            for(T2 sz=0; sz < nsnz; ++sz){
                                for(T2 sy=0; sy < nsny; ++sy, n++){

                                    T1 zsv = zmin + nk*dz + (sz+1)*dzs;
                                    T1 ysv = ymin + nj*dy + (sy+1)*dys;

                                    if ( cell_XmYpZp != std::numeric_limits<T2>::max() ) {
                                        nodes[n].pushOwner( cell_XmYpZp );
                                    }
                                    if ( cell_XpYpZp != std::numeric_limits<T2>::max() ) {
                                        nodes[n].pushOwner( cell_XpYpZp );
                                    }
                                    nodes[n].setXYZindex( x, ysv, zsv, n );

                                    if (ni>0){
                                        nodes[n].setPrimary(false);
                                    }
                                    else if (ni==0){
                                        nodes[n].setPrimary(false);
                                    }
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

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::interpSecondary() {
        T2 nPrimary = (ncx+1)*(ncy+1)*(ncz+1);
        for ( size_t n=nPrimary; n<nodes.size(); ++n ) {
            nodes[n].setNodeSlowness( computeSlowness({nodes[n].getX(),
                nodes[n].getY(), nodes[n].getZ()}, true) );
        }
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::checkPts(std::vector<sxyz<T1>> pts, const bool translated) const {

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


    template<typename T1, typename T2, typename NODE>
    T1 Grid3Drn<T1,T2,NODE>::getTraveltime(const sxyz<T1> &pt, const size_t nt) const {

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


    template<typename T1, typename T2, typename NODE>
    T1 Grid3Drn<T1,T2,NODE>::getTraveltime(const sxyz<T1>& Rx,
                                           T2& nodeParentRx, T2& cellParentRx,
                                           const size_t threadNo) const {

        // Calculate and return the traveltime for a Rx point.
        T2 nn = getNodeNo( Rx );
        if ( nn != std::numeric_limits<T2>::max() ) {
            nodeParentRx = nodes[nn].getNodeParent(threadNo);
            cellParentRx = nodes[nn].getCellParent(threadNo);
            return nodes[nn].getTT(threadNo);
        }
        //If Rx is not on a node:
        T1 slo = computeSlowness( Rx, true );

        T2 cellNo = getCellNo( Rx );
        T2 neibNo = this->neighbors[cellNo][0];
        T1 dt = computeDt(nodes[neibNo], Rx, slo);

        T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
        nodeParentRx = neibNo;
        cellParentRx = cellNo;
        for ( size_t k=1; k< this->neighbors[cellNo].size(); ++k ) {
            neibNo = this->neighbors[cellNo][k];
            dt = computeDt(nodes[neibNo], Rx, slo);
            if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
                traveltime =  nodes[neibNo].getTT(threadNo)+dt;
                nodeParentRx = neibNo;
            }
        }
        return traveltime;
    }

    
    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::grad(sxyz<T1>& g, const size_t i, const size_t j, const size_t k,
                                    const size_t nt) const {

        // compute average gradient for voxel (i,j,k)

        const size_t nnx = ncx+1;
        const size_t nny = ncy+1;

        g.x = 0.25*(nodes[(    k*nny+j  )*nnx+i+1].getTT(nt) - nodes[(    k*nny+j  )*nnx+i  ].getTT(nt) +
                    nodes[(    k*nny+j+1)*nnx+i+1].getTT(nt) - nodes[(    k*nny+j+1)*nnx+i  ].getTT(nt) +
                    nodes[((k+1)*nny+j  )*nnx+i+1].getTT(nt) - nodes[((k+1)*nny+j  )*nnx+i  ].getTT(nt) +
                    nodes[((k+1)*nny+j+1)*nnx+i+1].getTT(nt) - nodes[((k+1)*nny+j+1)*nnx+i  ].getTT(nt))/dx;
        g.y = 0.25*(nodes[(    k*nny+j+1)*nnx+i  ].getTT(nt) - nodes[(    k*nny+j  )*nnx+i  ].getTT(nt) +
                    nodes[(    k*nny+j+1)*nnx+i+1].getTT(nt) - nodes[(    k*nny+j  )*nnx+i+1].getTT(nt) +
                    nodes[((k+1)*nny+j+1)*nnx+i  ].getTT(nt) - nodes[((k+1)*nny+j  )*nnx+i  ].getTT(nt) +
                    nodes[((k+1)*nny+j+1)*nnx+i+1].getTT(nt) - nodes[((k+1)*nny+j  )*nnx+i+1].getTT(nt))/dy;
        g.z = 0.25*(nodes[((k+1)*nny+j  )*nnx+i  ].getTT(nt) - nodes[(    k*nny+j  )*nnx+i  ].getTT(nt) +
                    nodes[((k+1)*nny+j  )*nnx+i+1].getTT(nt) - nodes[(    k*nny+j  )*nnx+i+1].getTT(nt) +
                    nodes[((k+1)*nny+j+1)*nnx+i  ].getTT(nt) - nodes[(    k*nny+j+1)*nnx+i  ].getTT(nt) +
                    nodes[((k+1)*nny+j+1)*nnx+i+1].getTT(nt) - nodes[(    k*nny+j+1)*nnx+i+1].getTT(nt))/dz;
    }


    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::gradO2(sxyz<T1>& g, const sxyz<T1> &pt,
                                      const size_t nt) const {

        // compute travel time gradient (2nd order centered operator) at point pt

        T1 p1 = pt.x - dx/2.0;
        T1 p2 = p1 + dx;
        if ( p1 <= xmin ) {  // check if on grid edge or out of grid
            p1 = pt.x + dx/2.0;  // shift pt to allow interpolating in getTraveltime
            p2 = p1 + dx;
        } else if ( p2 >= xmax ) {
            p2 = pt.x - dx/2.0;
            p1 = p2 - dx;
        }
        g.x = (getTraveltime({p2, pt.y, pt.z}, nt) - getTraveltime({p1, pt.y, pt.z}, nt)) / dx;

        p1 = pt.y - dy/2.0;
        p2 = p1 + dy;
        if ( p1 <= ymin ) {
            p1 = pt.y + dy/2.0;
            p2 = p1 + dy;
        } else if ( p2 >= ymax ) {
            p2 = pt.y - dy/2.0;
            p1 = p2 - dy;
        }
        g.y = (getTraveltime({pt.x, p2, pt.z}, nt) - getTraveltime({pt.x, p1, pt.z}, nt)) / dy;

        p1 = pt.z - dz/2.0;
        p2 = p1 + dz;
        if ( p1 <= zmin ) {
            p1 = pt.z + dz/2.0;
            p2 = p1 + dz;
        } else if ( p2 >= zmax ) {
            p2 = pt.z - dz/2.0;
            p1 = p2 - dz;
        }
        g.z = (getTraveltime({pt.x, pt.y, p2}, nt) - getTraveltime({pt.x, pt.y, p1}, nt)) / dz;

    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::grad(sxyz<T1>& g, const sxyz<T1> &pt,
                                    const size_t nt) const {

        // compute travel time gradient (4th order centered operator) at point pt

        static const T1 k1 = 1./24.;
        static const T1 k2 = 9./8.;

        T1 p1 = pt.x - dx;
        T1 p2 = p1 + 0.5*dx;
        T1 p3 = p1 + 1.5*dx;
        T1 p4 = p1 + 2.0*dx;
        if ( p1 <= xmin ) {  // check if on grid edge or out of grid
            p1 = xmin;  // shift pt to allow interpolating in getTraveltime
            p2 = p1 + 0.5*dx;
            p3 = p1 + 1.5*dx;
            p4 = p1 + 2.0*dx;
        } else if ( p4 >= xmax ) {
            p4 = xmax;
            p3 = p4 - 0.5*dx;
            p2 = p4 - 1.5*dx;
            p1 = p4 - 2.0*dx;
        }
        g.x = (k1 * getTraveltime({p1, pt.y, pt.z}, nt) -
               k2 * getTraveltime({p2, pt.y, pt.z}, nt) +
               k2 * getTraveltime({p3, pt.y, pt.z}, nt) -
               k1 * getTraveltime({p4, pt.y, pt.z}, nt)) / dx;

        p1 = pt.y - dy/2.0;
        p2 = p1 + 0.5*dy;
        p3 = p1 + 1.5*dy;
        p4 = p1 + 2.0*dy;
        if ( p1 <= ymin ) {
            p1 = ymin;
            p2 = p1 + 0.5*dy;
            p3 = p1 + 1.5*dy;
            p4 = p1 + 2.0*dy;
        } else if ( p4 >= ymax ) {
            p4 = ymax;
            p3 = p4 - 0.5*dy;
            p2 = p4 - 1.5*dy;
            p1 = p4 - 2.0*dy;
        }
        g.y = (k1 * getTraveltime({pt.x, p1, pt.z}, nt) -
               k2 * getTraveltime({pt.x, p2, pt.z}, nt) +
               k2 * getTraveltime({pt.x, p3, pt.z}, nt) -
               k1 * getTraveltime({pt.x, p4, pt.z}, nt)) / dy;

        p1 = pt.z - dz/2.0;
        p2 = p1 + 0.5*dz;
        p3 = p1 + 1.5*dz;
        p4 = p1 + 2.0*dz;
        if ( p1 <= zmin ) {
            p1 = zmin;
            p2 = p1 + 0.5*dz;
            p3 = p1 + 1.5*dz;
            p4 = p1 + 2.0*dz;
        } else if ( p4 >= zmax ) {
            p4 = zmax;
            p3 = p4 - 0.5*dz;
            p2 = p4 - 1.5*dz;
            p1 = p4 - 2.0*dz;
        }
        g.z = (k1 * getTraveltime({pt.x, pt.y, p1}, nt) -
               k2 * getTraveltime({pt.x, pt.y, p2}, nt) +
               k2 * getTraveltime({pt.x, pt.y, p3}, nt) -
               k1 * getTraveltime({pt.x, pt.y, p4}, nt)) / dz;
    }

    template<typename T1, typename T2, typename NODE>
    T1 Grid3Drn<T1,T2,NODE>::getTraveltimeFromRaypath(const std::vector<sxyz<T1>>& Tx,
                                                      const std::vector<T1>& t0,
                                                      const sxyz<T1> &Rx,
                                                      const size_t threadNo) const {
        T1 tt = 0.0;
        T1 s1, s2;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return tt;
            }
        }

        sxyz<T1> prev_pt( Rx );
        sxyz<T1> curr_pt( Rx );
        s1 = computeSlowness( curr_pt, true );
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

            s2 = computeSlowness( curr_pt, true );
            tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
            s1 = s2;
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
                        s2 = computeSlowness( Tx[ns], true );
                        tt += t0[ns] + 0.5*(s1 + s2) * prev_pt.getDistance( Tx[ns] );
                    } else {
                        // to intersection
                        s2 = computeSlowness( curr_pt, true );
                        tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
                        s1 = s2;
                        // to Tx
                        s2 = computeSlowness( Tx[ns], true );
                        tt += t0[ns] + 0.5*(s1 + s2) * curr_pt.getDistance( Tx[ns] );
                    }

                    reachedTx = true;
                }
            }
        }

        return tt;
    }


    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                          const sxyz<T1> &Rx,
                                          std::vector<sxyz<T1>> &r_data,
                                          const size_t threadNo) const {

        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                return;
            }
        }

        sxyz<T1> curr_pt( Rx );
        // distance between opposite nodes of a voxel
        const T1 maxDist = sqrt( dx*dx + dy*dy + dz*dz );
        sxyz<T1> g;
#ifdef DEBUG_RP
        std::cout << "\n\nStarting raypath computation\n  curr_pt = " << curr_pt << '\n';
#endif
        bool reachedTx = false;
        while ( reachedTx == false ) {

            grad(g, curr_pt, threadNo);
            g *= -1.0;
#ifdef DEBUG_RP
            std::cout << "  g = " << g << '\n';
#endif
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
            std::cout << "  curr_pt = " << curr_pt << '\n';
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

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                          const std::vector<T1>& t0,
                                          const sxyz<T1> &Rx,
                                          std::vector<sxyz<T1>> &r_data,
                                          T1 &tt,
                                          const size_t threadNo) const {
        tt = 0.0;
        T1 s1, s2;

        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        sxyz<T1> curr_pt( Rx );
        s1 = computeSlowness( curr_pt, true );
        // distance between opposite nodes of a voxel
        const T1 maxDist = sqrt( dx*dx + dy*dy + dz*dz );
        sxyz<T1> g;
#ifdef DEBUG_RP
        std::cout << '\n' << curr_pt << '\n';
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
#ifdef DEBUG_RP
            std::cout << curr_pt << '\n';
#endif
            s2 = computeSlowness( curr_pt, true );
            tt += 0.5*(s1 + s2) * r_data.back().getDistance( curr_pt );
            s1 = s2;
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

                    if ( curr_pt.getDistance(r_data.back()) > dist  ||  // we do not intersect a plane
                        curr_pt == Tx[ns] ) {  // we have arrived

                        s2 = computeSlowness( Tx[ns], true );
                        tt += t0[ns] + 0.5*(s1 + s2) * r_data.back().getDistance( Tx[ns] );
                        r_data.push_back( Tx[ns] );
#ifdef DEBUG_RP
                        std::cout << Tx[ns] << "\n\n";
#endif
                    } else {
                        // to intersection
                        s2 = computeSlowness( curr_pt, true );
                        tt += 0.5*(s1 + s2) * r_data.back().getDistance( curr_pt );
                        r_data.push_back( curr_pt );
#ifdef DEBUG_RP
                        std::cout << curr_pt << '\n';
#endif
                        s1 = s2;
                        // to Tx
                        s2 = computeSlowness( Tx[ns], true );
                        tt += t0[ns] + 0.5*(s1 + s2) * curr_pt.getDistance( Tx[ns] );
                        r_data.push_back( Tx[ns] );
#ifdef DEBUG_RP
                        std::cout << Tx[ns] << "\n\n";
#endif
                    }

                    reachedTx = true;
                }
            }
        }
    }


    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                          const std::vector<T1>& t0,
                                          const sxyz<T1> &Rx,
                                          std::vector<sijv<T1>>& m_data,
                                          const size_t RxNo,
                                          T1 &tt,
                                          const size_t threadNo) const {
        tt = 0.0;
        T1 s1, s2;

        const size_t nnx = ncx+1;
        const size_t nny = ncy+1;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                return;
            }
        }

        sxyz<T1> curr_pt( Rx ), prev_pt( Rx ), mid_pt;
        s1 = computeSlowness( curr_pt, true );
        sijv<T1> m;
        m.i = RxNo;

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

            s2 = computeSlowness( curr_pt, true );
            tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
            s1 = s2;
            prev_pt = curr_pt;

            // compute terms of matrix M
            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
            T1 s = computeSlowness(mid_pt, true);
            s *= s;
            T1 ds = curr_pt.getDistance( prev_pt );

            size_t ix = (mid_pt.x-xmin)/dx;
            size_t iy = (mid_pt.y-ymin)/dy;
            size_t iz = (mid_pt.z-zmin)/dz;
            for ( size_t ii=0; ii<2; ++ii ) {
                for ( size_t jj=0; jj<2; ++jj ) {
                    for ( size_t kk=0; kk<2; ++kk ) {
                        size_t iv = ix+ii;
                        size_t jv = iy+jj;
                        size_t kv = iz+kk;
                        T1 dvdv = (1. - std::abs(mid_pt.x - iv*dx)/dx) *
                        (1. - std::abs(mid_pt.y - jv*dy)/dy) *
                        (1. - std::abs(mid_pt.z - kv*dz)/dz);

                        m.j = (kv*nny+jv)*nnx+iv;
                        m.v = -s * ds * dvdv;

                        bool found = false;
                        for ( size_t nm=0; nm<m_data.size(); ++nm ) {
                            if ( m_data[nm].j == m.j ) {
                                m_data[nm].v += m.v;
                                found = true;
                                break;
                            }
                        }
                        if ( found == false ) {
                            m_data.push_back(m);
                        }
                    }
                }
            }


            // are we close enough to one of the Tx nodes ?
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

                    if ( curr_pt.getDistance(prev_pt) > dist  ||  // we do not intersect a plane
                        curr_pt == Tx[ns] ) {  // we have arrived

                        s2 = computeSlowness( Tx[ns], true );
                        tt += t0[ns] + 0.5*(s1 + s2) * prev_pt.getDistance( Tx[ns] );

                        // compute terms of matrix M
                        mid_pt = static_cast<T1>(0.5)*(Tx[ns] + prev_pt);
                        s = computeSlowness(mid_pt, true);
                        s *= s;
                        ds = Tx[ns].getDistance( prev_pt );

                        ix = (mid_pt.x-xmin)/dx;
                        iy = (mid_pt.y-ymin)/dy;
                        iz = (mid_pt.z-zmin)/dz;
                        for ( size_t ii=0; ii<2; ++ii ) {
                            for ( size_t jj=0; jj<2; ++jj ) {
                                for ( size_t kk=0; kk<2; ++kk ) {
                                    size_t iv = ix+ii;
                                    size_t jv = iy+jj;
                                    size_t kv = iz+kk;
                                    T1 dvdv = (1. - std::abs(mid_pt.x - iv*dx)/dx) *
                                    (1. - std::abs(mid_pt.y - jv*dy)/dy) *
                                    (1. - std::abs(mid_pt.z - kv*dz)/dz);

                                    m.j = (kv*nny+jv)*nnx+iv;
                                    m.v = -s * ds * dvdv;

                                    bool found = false;
                                    for ( size_t nm=0; nm<m_data.size(); ++nm ) {
                                        if ( m_data[nm].j == m.j ) {
                                            m_data[nm].v += m.v;
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found == false ) {
                                        m_data.push_back(m);
                                    }
                                }
                            }
                        }
                    } else {
                        // to intersection
                        s2 = computeSlowness( curr_pt, true );
                        tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
                        s1 = s2;

                        // compute terms of matrix M
                        mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                        s = computeSlowness(mid_pt, true);
                        s *= s;
                        ds = curr_pt.getDistance( prev_pt );

                        ix = (mid_pt.x-xmin)/dx;
                        iy = (mid_pt.y-ymin)/dy;
                        iz = (mid_pt.z-zmin)/dz;
                        for ( size_t ii=0; ii<2; ++ii ) {
                            for ( size_t jj=0; jj<2; ++jj ) {
                                for ( size_t kk=0; kk<2; ++kk ) {
                                    size_t iv = ix+ii;
                                    size_t jv = iy+jj;
                                    size_t kv = iz+kk;
                                    T1 dvdv = (1. - std::abs(mid_pt.x - iv*dx)/dx) *
                                    (1. - std::abs(mid_pt.y - jv*dy)/dy) *
                                    (1. - std::abs(mid_pt.z - kv*dz)/dz);

                                    m.j = (kv*nny+jv)*nnx+iv;
                                    m.v = -s * ds * dvdv;

                                    bool found = false;
                                    for ( size_t nm=0; nm<m_data.size(); ++nm ) {
                                        if ( m_data[nm].j == m.j ) {
                                            m_data[nm].v += m.v;
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found == false ) {
                                        m_data.push_back(m);
                                    }

                                }
                            }
                        }

                        // to Tx
                        s2 = computeSlowness( Tx[ns], true );
                        tt += t0[ns] + 0.5*(s1 + s2) * curr_pt.getDistance( Tx[ns] );

                        // compute terms of matrix M
                        mid_pt = static_cast<T1>(0.5)*(Tx[ns] + curr_pt);
                        s = computeSlowness(mid_pt, true);
                        s *= s;
                        ds = Tx[ns].getDistance( curr_pt );

                        ix = (mid_pt.x-xmin)/dx;
                        iy = (mid_pt.y-ymin)/dy;
                        iz = (mid_pt.z-zmin)/dz;
                        for ( size_t ii=0; ii<2; ++ii ) {
                            for ( size_t jj=0; jj<2; ++jj ) {
                                for ( size_t kk=0; kk<2; ++kk ) {
                                    size_t iv = ix+ii;
                                    size_t jv = iy+jj;
                                    size_t kv = iz+kk;
                                    T1 dvdv = (1. - std::abs(mid_pt.x - iv*dx)/dx) *
                                    (1. - std::abs(mid_pt.y - jv*dy)/dy) *
                                    (1. - std::abs(mid_pt.z - kv*dz)/dz);

                                    m.j = (kv*nny+jv)*nnx+iv;
                                    m.v = -s * ds * dvdv;

                                    bool found = false;
                                    for ( size_t nm=0; nm<m_data.size(); ++nm ) {
                                        if ( m_data[nm].j == m.j ) {
                                            m_data[nm].v += m.v;
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found == false ) {
                                        m_data.push_back(m);
                                    }
                                }
                            }
                        }
                    }

                    reachedTx = true;
                }
            }
        }
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                          const std::vector<T1>& t0,
                                          const sxyz<T1> &Rx,
                                          std::vector<siv<T1>> &l_data,
                                          T1 &tt,
                                          const size_t threadNo) const {
        tt = 0.0;
        T1 s1, s2;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        sxyz<T1> curr_pt( Rx ), prev_pt( Rx );
        s1 = computeSlowness( curr_pt, true );
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
            s2 = computeSlowness( curr_pt, true );
            tt += 0.5*(s1 + s2) * cell.v;
            s1 = s2;
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
                        s2 = computeSlowness( Tx[ns], true );
                        tt += 0.5*(s1 + s2) * cell.v;

                    } else {
                        // to intersection
                        mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
                        cell.i = getCellNo(mid_pt);
                        cell.v = curr_pt.getDistance(prev_pt);
                        l_data.push_back(cell);
                        s2 = computeSlowness( curr_pt, true );
                        tt += 0.5*(s1 + s2) * cell.v;

                        s1 = s2;
                        // to Tx
                        cell.i = getCellNo(Tx[ns]);
                        cell.v = Tx[ns].getDistance(curr_pt);
                        l_data.push_back(cell);
                        s2 = computeSlowness( Tx[ns], true );
                        tt += 0.5*(s1 + s2) * cell.v;
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                          const std::vector<T1>& t0,
                                          const sxyz<T1> &Rx,
                                          std::vector<sxyz<T1>> &r_data,
                                          std::vector<siv<T1>> &l_data,
                                          T1 &tt,
                                          const size_t threadNo) const {
        tt = 0.0;
        T1 s1, s2;
        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        sxyz<T1> curr_pt( Rx );
        s1 = computeSlowness( curr_pt, true );
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
            s2 = computeSlowness( curr_pt, true );
            tt += 0.5*(s1 + s2) * cell.v;
            s1 = s2;
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
                        s2 = computeSlowness( Tx[ns], true );
                        tt += 0.5*(s1 + s2) * cell.v;
                        r_data.push_back( Tx[ns] );
                    } else {
                        // to intersection
                        mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
                        cell.i = getCellNo(mid_pt);
                        cell.v = curr_pt.getDistance(r_data.back());
                        l_data.push_back(cell);
                        s2 = computeSlowness( curr_pt, true );
                        tt += 0.5*(s1 + s2) * cell.v;
                        r_data.push_back( curr_pt );

                        s1 = s2;
                        // to Tx
                        cell.i = getCellNo(Tx[ns]);
                        cell.v = Tx[ns].getDistance(r_data.back());
                        l_data.push_back(cell);
                        s2 = computeSlowness( Tx[ns], true );
                        tt += 0.5*(s1 + s2) * cell.v;
                        r_data.push_back( Tx[ns] );
                    }

                    tt += t0[ns];
                    reachedTx = true;
                }
            }
        }
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                          const std::vector<T1>& t0,
                                          const sxyz<T1> &Rx,
                                          std::vector<sxyz<T1>>& r_data,
                                          std::vector<sijv<T1>>& m_data,
                                          const size_t RxNo,
                                          T1 &tt,
                                          const size_t threadNo) const {
        tt = 0.0;
        T1 s1, s2;

        const size_t nnx = ncx+1;
        const size_t nny = ncy+1;

        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                return;
            }
        }

        sxyz<T1> curr_pt( Rx ), prev_pt, mid_pt;
        s1 = computeSlowness( curr_pt, true );
        sijv<T1> m;
        m.i = RxNo;

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

            prev_pt = r_data.back();
            s2 = computeSlowness( curr_pt, true );
            tt += 0.5*(s1 + s2) * r_data.back().getDistance( curr_pt );
            s1 = s2;
            r_data.push_back( curr_pt );

            // compute terms of matrix M
            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
            T1 s = computeSlowness(mid_pt, true);
            s *= s;
            T1 ds = curr_pt.getDistance( prev_pt );

            size_t ix = (mid_pt.x-xmin)/dx;
            size_t iy = (mid_pt.y-ymin)/dy;
            size_t iz = (mid_pt.z-zmin)/dz;
            for ( size_t ii=0; ii<2; ++ii ) {
                for ( size_t jj=0; jj<2; ++jj ) {
                    for ( size_t kk=0; kk<2; ++kk ) {
                        size_t iv = ix+ii;
                        size_t jv = iy+jj;
                        size_t kv = iz+kk;
                        T1 dvdv = (1. - std::abs(mid_pt.x - iv*dx)/dx) *
                        (1. - std::abs(mid_pt.y - jv*dy)/dy) *
                        (1. - std::abs(mid_pt.z - kv*dz)/dz);

                        m.j = (kv*nny+jv)*nnx+iv;
                        m.v = -s * ds * dvdv;

                        bool found = false;
                        for ( size_t nm=0; nm<m_data.size(); ++nm ) {
                            if ( m_data[nm].j == m.j ) {
                                m_data[nm].v += m.v;
                                found = true;
                                break;
                            }
                        }
                        if ( found == false ) {
                            m_data.push_back(m);
                        }
                    }
                }
            }


            // are we close enough to one of the Tx nodes ?
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

                    if ( curr_pt.getDistance(r_data.back()) > dist  ||  // we do not intersect a plane
                        curr_pt == Tx[ns] ) {  // we have arrived

                        prev_pt = r_data.back();
                        s2 = computeSlowness( Tx[ns], true );
                        tt += t0[ns] + 0.5*(s1 + s2) * r_data.back().getDistance( Tx[ns] );
                        r_data.push_back( Tx[ns] );

                        // compute terms of matrix M
                        mid_pt = static_cast<T1>(0.5)*(Tx[ns] + prev_pt);
                        s = computeSlowness(mid_pt, true);
                        s *= s;
                        ds = Tx[ns].getDistance( prev_pt );

                        ix = (mid_pt.x-xmin)/dx;
                        iy = (mid_pt.y-ymin)/dy;
                        iz = (mid_pt.z-zmin)/dz;
                        for ( size_t ii=0; ii<2; ++ii ) {
                            for ( size_t jj=0; jj<2; ++jj ) {
                                for ( size_t kk=0; kk<2; ++kk ) {
                                    size_t iv = ix+ii;
                                    size_t jv = iy+jj;
                                    size_t kv = iz+kk;
                                    T1 dvdv = (1. - std::abs(mid_pt.x - iv*dx)/dx) *
                                    (1. - std::abs(mid_pt.y - jv*dy)/dy) *
                                    (1. - std::abs(mid_pt.z - kv*dz)/dz);

                                    m.j = (kv*nny+jv)*nnx+iv;
                                    m.v = -s * ds * dvdv;

                                    bool found = false;
                                    for ( size_t nm=0; nm<m_data.size(); ++nm ) {
                                        if ( m_data[nm].j == m.j ) {
                                            m_data[nm].v += m.v;
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found == false ) {
                                        m_data.push_back(m);
                                    }
                                }
                            }
                        }
                    } else {
                        // to intersection
                        s2 = computeSlowness( curr_pt, true );
                        tt += 0.5*(s1 + s2) * r_data.back().getDistance( curr_pt );
                        r_data.push_back( curr_pt );
                        s1 = s2;

                        prev_pt = r_data.back();
                        // compute terms of matrix M
                        mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                        s = computeSlowness(mid_pt, true);
                        s *= s;
                        ds = curr_pt.getDistance( prev_pt );

                        ix = (mid_pt.x-xmin)/dx;
                        iy = (mid_pt.y-ymin)/dy;
                        iz = (mid_pt.z-zmin)/dz;
                        for ( size_t ii=0; ii<2; ++ii ) {
                            for ( size_t jj=0; jj<2; ++jj ) {
                                for ( size_t kk=0; kk<2; ++kk ) {
                                    size_t iv = ix+ii;
                                    size_t jv = iy+jj;
                                    size_t kv = iz+kk;
                                    T1 dvdv = (1. - std::abs(mid_pt.x - iv*dx)/dx) *
                                    (1. - std::abs(mid_pt.y - jv*dy)/dy) *
                                    (1. - std::abs(mid_pt.z - kv*dz)/dz);

                                    m.j = (kv*nny+jv)*nnx+iv;
                                    m.v = -s * ds * dvdv;

                                    bool found = false;
                                    for ( size_t nm=0; nm<m_data.size(); ++nm ) {
                                        if ( m_data[nm].j == m.j ) {
                                            m_data[nm].v += m.v;
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found == false ) {
                                        m_data.push_back(m);
                                    }

                                }
                            }
                        }

                        // to Tx
                        prev_pt = r_data.back();
                        s2 = computeSlowness( Tx[ns], true );
                        tt += t0[ns] + 0.5*(s1 + s2) * curr_pt.getDistance( Tx[ns] );
                        r_data.push_back( Tx[ns] );

                        // compute terms of matrix M
                        mid_pt = static_cast<T1>(0.5)*(Tx[ns] + prev_pt);
                        s = computeSlowness(mid_pt, true);
                        s *= s;
                        ds = Tx[ns].getDistance( prev_pt );

                        ix = (mid_pt.x-xmin)/dx;
                        iy = (mid_pt.y-ymin)/dy;
                        iz = (mid_pt.z-zmin)/dz;
                        for ( size_t ii=0; ii<2; ++ii ) {
                            for ( size_t jj=0; jj<2; ++jj ) {
                                for ( size_t kk=0; kk<2; ++kk ) {
                                    size_t iv = ix+ii;
                                    size_t jv = iy+jj;
                                    size_t kv = iz+kk;
                                    T1 dvdv = (1. - std::abs(mid_pt.x - iv*dx)/dx) *
                                    (1. - std::abs(mid_pt.y - jv*dy)/dy) *
                                    (1. - std::abs(mid_pt.z - kv*dz)/dz);

                                    m.j = (kv*nny+jv)*nnx+iv;
                                    m.v = -s * ds * dvdv;

                                    bool found = false;
                                    for ( size_t nm=0; nm<m_data.size(); ++nm ) {
                                        if ( m_data[nm].j == m.j ) {
                                            m_data[nm].v += m.v;
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found == false ) {
                                        m_data.push_back(m);
                                    }
                                }
                            }
                        }
                    }

                    reachedTx = true;
                }
            }
        }
    }

    template<typename T1, typename T2, typename NODE>
    T1 Grid3Drn<T1,T2,NODE>::computeSlowness(sxyz<T1> pt, const bool isTranslated) const {

        if (this->translateOrigin == true && isTranslated == false) {
            pt -= this->origin;
        }
        const size_t nnx = ncx+1;
        const size_t nny = ncy+1;
        const size_t nnz = ncz+1;

        // are we on an node, an edge or a face?
        ptrdiff_t onX = -1;
        ptrdiff_t onY = -1;
        ptrdiff_t onZ = -1;
        for ( size_t n=0; n<nnx; ++n ) {
            if ( std::abs(pt.x - (xmin+n*dx)) < small2 ) {
                onX = n;
                break;
            }
        }
        for ( size_t n=0; n<nny; ++n ) {
            if ( std::abs(pt.y - (ymin+n*dy)) < small2 ) {
                onY = n;
                break;
            }
        }
        for ( size_t n=0; n<nnz; ++n ) {
            if ( std::abs(pt.z - (zmin+n*dz)) < small2 ) {
                onZ = n;
                break;
            }
        }

        if ( onX!=-1 && onY!=-1 && onZ!=-1 ) {
            return nodes[(onZ*nny+onY)*nnx+onX].getNodeSlowness();
        } else if ( onX!=-1 && onY!=-1 ) {
            T2 k = static_cast<T2>( small + (pt.z-zmin)/dz );
            T1 s[2];
            T1 x[3];
            if ( processVel ) {
                s[0] = 1.0 / nodes[(k*nny+onY)*nnx+onX].getNodeSlowness();
                s[1] = 1.0 / nodes[((k+1)*nny+onY)*nnx+onX].getNodeSlowness();
            } else {
                s[0] = nodes[(k*nny+onY)*nnx+onX].getNodeSlowness();
                s[1] = nodes[((k+1)*nny+onY)*nnx+onX].getNodeSlowness();
            }
            x[0] = pt.z;
            x[1] = zmin + k*dz;
            x[2] = zmin + (k+1)*dz;

            if ( processVel )
                return 1.0 / Interpolator<T1>::linear(x, s);
            else
                return Interpolator<T1>::linear(x, s);

        } else if ( onX!=-1 && onZ!=-1 ) {
            T2 j = static_cast<T2>( small + (pt.y-ymin)/dy );
            T1 s[2];
            T1 x[3];
            if ( processVel ) {
                s[0] = 1.0 / nodes[(onZ*nny+j)*nnx+onX].getNodeSlowness();
                s[1] = 1.0 / nodes[(onZ*nny+j+1)*nnx+onX].getNodeSlowness();
            } else {
                s[0] = nodes[(onZ*nny+j)*nnx+onX].getNodeSlowness();
                s[1] = nodes[(onZ*nny+j+1)*nnx+onX].getNodeSlowness();
            }
            x[0] = pt.y;
            x[1] = ymin + j*dy;
            x[2] = ymin + (j+1)*dy;

            if ( processVel )
                return 1.0 / Interpolator<T1>::linear(x, s);
            else
                return Interpolator<T1>::linear(x, s);

        } else if ( onY!=-1 && onZ!=-1 ) {
            T2 i = static_cast<T2>( small + (pt.x-xmin)/dx );
            T1 s[2];
            T1 x[3];
            if ( processVel ) {
                s[0] = 1.0 / nodes[(onZ*nny+onY)*nnx+i].getNodeSlowness();
                s[1] = 1.0 / nodes[(onZ*nny+onY)*nnx+i+1].getNodeSlowness();
            } else {
                s[0] = nodes[(onZ*nny+onY)*nnx+i].getNodeSlowness();
                s[1] = nodes[(onZ*nny+onY)*nnx+i+1].getNodeSlowness();
            }
            x[0] = pt.x;
            x[1] = xmin + i*dx;
            x[2] = xmin + (i+1)*dx;

            if ( processVel )
                return 1.0 / Interpolator<T1>::linear(x, s);
            else
                return Interpolator<T1>::linear(x, s);
        } else if ( onX!=-1 ) {
            T2 j = static_cast<T2>( small + (pt.y-ymin)/dy );
            T2 k = static_cast<T2>( small + (pt.z-zmin)/dz );
            T1 s[4];
            T1 x[3];
            T1 y[3];
            if ( processVel ) {
                s[0] = 1.0 / nodes[((k  )*nny+j  )*nnx+onX].getNodeSlowness();
                s[1] = 1.0 / nodes[((k+1)*nny+j  )*nnx+onX].getNodeSlowness();
                s[2] = 1.0 / nodes[((k  )*nny+j+1)*nnx+onX].getNodeSlowness();
                s[3] = 1.0 / nodes[((k+1)*nny+j+1)*nnx+onX].getNodeSlowness();
            } else {
                s[0] = nodes[((k  )*nny+j  )*nnx+onX].getNodeSlowness();
                s[1] = nodes[((k+1)*nny+j  )*nnx+onX].getNodeSlowness();
                s[2] = nodes[((k  )*nny+j+1)*nnx+onX].getNodeSlowness();
                s[3] = nodes[((k+1)*nny+j+1)*nnx+onX].getNodeSlowness();
            }
            x[0] = pt.y;
            y[0] = pt.z;
            x[1] = ymin + j*dy;
            y[1] = zmin + k*dz;
            x[2] = ymin + (j+1)*dy;
            y[2] = zmin + (k+1)*dz;

            if ( processVel )
                return 1.0 / Interpolator<T1>::bilinear(x, y, s);
            else
                return Interpolator<T1>::bilinear(x, y, s);

        } else if ( onY!=-1 ) {
            T2 i = static_cast<T2>( small + (pt.x-xmin)/dx );
            T2 k = static_cast<T2>( small + (pt.z-zmin)/dz );
            T1 s[4];
            T1 x[3];
            T1 y[3];
            if ( processVel ) {
                s[0] = 1.0 / nodes[((k  )*nny+onY)*nnx+i  ].getNodeSlowness();
                s[1] = 1.0 / nodes[((k+1)*nny+onY)*nnx+i  ].getNodeSlowness();
                s[2] = 1.0 / nodes[((k  )*nny+onY)*nnx+i+1].getNodeSlowness();
                s[3] = 1.0 / nodes[((k+1)*nny+onY)*nnx+i+1].getNodeSlowness();
            } else {
                s[0] = nodes[((k  )*nny+onY)*nnx+i  ].getNodeSlowness();
                s[1] = nodes[((k+1)*nny+onY)*nnx+i  ].getNodeSlowness();
                s[2] = nodes[((k  )*nny+onY)*nnx+i+1].getNodeSlowness();
                s[3] = nodes[((k+1)*nny+onY)*nnx+i+1].getNodeSlowness();
            }
            x[0] = pt.x;
            y[0] = pt.z;
            x[1] = xmin + i*dx;
            y[1] = zmin + k*dz;
            x[2] = xmin + (i+1)*dx;
            y[2] = zmin + (k+1)*dz;

            if ( processVel )
                return 1.0 / Interpolator<T1>::bilinear(x, y, s);
            else
                return Interpolator<T1>::bilinear(x, y, s);

        } else if ( onZ!=-1 ) {
            T2 i = static_cast<T2>( small + (pt.x-xmin)/dx );
            T2 j = static_cast<T2>( small + (pt.y-ymin)/dy );
            T1 s[4];
            T1 x[3];
            T1 y[3];
            if ( processVel ) {
                s[0] = 1.0 / nodes[(onZ*nny+j  )*nnx+i  ].getNodeSlowness();
                s[1] = 1.0 / nodes[(onZ*nny+j+1)*nnx+i  ].getNodeSlowness();
                s[2] = 1.0 / nodes[(onZ*nny+j  )*nnx+i+1].getNodeSlowness();
                s[3] = 1.0 / nodes[(onZ*nny+j+1)*nnx+i+1].getNodeSlowness();
            } else {
                s[0] = nodes[(onZ*nny+j  )*nnx+i  ].getNodeSlowness();
                s[1] = nodes[(onZ*nny+j+1)*nnx+i  ].getNodeSlowness();
                s[2] = nodes[(onZ*nny+j  )*nnx+i+1].getNodeSlowness();
                s[3] = nodes[(onZ*nny+j+1)*nnx+i+1].getNodeSlowness();
            }
            x[0] = pt.x;
            y[0] = pt.y;
            x[1] = xmin + i*dx;
            y[1] = ymin + j*dy;
            x[2] = xmin + (i+1)*dx;
            y[2] = ymin + (j+1)*dy;

            if ( processVel )
                return 1.0 / Interpolator<T1>::bilinear(x, y, s);
            else
                return Interpolator<T1>::bilinear(x, y, s);

        } else {
            T2 i = static_cast<T2>( small + (pt.x-xmin)/dx );
            T2 j = static_cast<T2>( small + (pt.y-ymin)/dy );
            T2 k = static_cast<T2>( small + (pt.z-zmin)/dz );
            T1 s[8];
            T1 x[3];
            T1 y[3];
            T1 z[3];

            if ( processVel ) {
                s[0] = 1.0 / nodes[((k  )*nny+j  )*nnx+i  ].getNodeSlowness();
                s[1] = 1.0 / nodes[((k+1)*nny+j  )*nnx+i  ].getNodeSlowness();
                s[2] = 1.0 / nodes[((k  )*nny+j+1)*nnx+i  ].getNodeSlowness();
                s[3] = 1.0 / nodes[((k+1)*nny+j+1)*nnx+i  ].getNodeSlowness();
                s[4] = 1.0 / nodes[((k  )*nny+j  )*nnx+i+1].getNodeSlowness();
                s[5] = 1.0 / nodes[((k+1)*nny+j  )*nnx+i+1].getNodeSlowness();
                s[6] = 1.0 / nodes[((k  )*nny+j+1)*nnx+i+1].getNodeSlowness();
                s[7] = 1.0 / nodes[((k+1)*nny+j+1)*nnx+i+1].getNodeSlowness();
            } else {
                s[0] = nodes[((k  )*nny+j  )*nnx+i  ].getNodeSlowness();
                s[1] = nodes[((k+1)*nny+j  )*nnx+i  ].getNodeSlowness();
                s[2] = nodes[((k  )*nny+j+1)*nnx+i  ].getNodeSlowness();
                s[3] = nodes[((k+1)*nny+j+1)*nnx+i  ].getNodeSlowness();
                s[4] = nodes[((k  )*nny+j  )*nnx+i+1].getNodeSlowness();
                s[5] = nodes[((k+1)*nny+j  )*nnx+i+1].getNodeSlowness();
                s[6] = nodes[((k  )*nny+j+1)*nnx+i+1].getNodeSlowness();
                s[7] = nodes[((k+1)*nny+j+1)*nnx+i+1].getNodeSlowness();
            }
            x[0] = pt.x;
            y[0] = pt.y;
            z[0] = pt.z;
            x[1] = xmin + i*dx;
            y[1] = ymin + j*dy;
            z[1] = zmin + k*dz;
            x[2] = xmin + (i+1)*dx;
            y[2] = ymin + (j+1)*dy;
            z[2] = zmin + (k+1)*dz;

            if ( processVel )
                return 1.0 / Interpolator<T1>::trilinear(x, y, z, s);
            else
                return Interpolator<T1>::trilinear(x, y, z, s);

        }

    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::saveTT(const std::string & fname, const int all,
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

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::loadTT(const std::string & fname, const int all,
                                      const size_t nt,
                                      const int format) const {

        if ( format == 1 ) {
            std::string filename = fname+".dat";
            std::ifstream fin(filename.c_str());
            T1 x, y, z, tt;
            for ( T2 n=0; n<nodes.size(); ++n ) {
                if ( nodes[n].isPrimary() || all==1 ) {
                    fin >> x >> y >> z >> tt;
                    nodes[n].setTT(tt, nt);
                }
            }
            fin.close();
        } else if ( format == 2 ) {
#ifdef VTK
            std::string filename = fname+".vtr";

            vtkSmartPointer<vtkXMLRectilinearGridReader> reader =
            vtkSmartPointer<vtkXMLRectilinearGridReader>::New();

            reader->SetFileName( filename.c_str() );
            reader->Update();

            for ( size_t n=0; n<nodes.size(); ++n ) {
                if ( nodes[n].isPrimary() ) {
                    double x[3] = {nodes[n].getX(), nodes[n].getY(), nodes[n].getZ()};
                    vtkIdType id = reader->GetOutput()->FindPoint(x);
                    nodes[n].setTT(reader->GetOutput()->GetPointData()->GetArray("Travel time")->GetTuple1(id), nt);
                }
            }
#else
            std::cerr << "VTK not included during compilation.\n";
#endif
        } else if ( format == 3 ) {
            std::string filename = fname+".bin";
            std::ifstream fin(filename.c_str(),std::ios::binary);
            for ( T2 n=0; n<nodes.size(); ++n ) {
                if ( nodes[n].isPrimary() || all==1 ) {
                    T1 tmp[4];
                    fin.read( (char*)tmp, 4*sizeof(T1) );
                    nodes[n].setTT(tmp[3], nt);
                }
            }
            fin.close();
        } else {
            throw std::runtime_error("Unsupported format for traveltimes");
        }
    }


    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::sweep(const std::vector<bool>& frozen,
                                     const size_t threadNo) const {

        // sweep first direction
        for ( size_t k=0; k<=ncz; ++k ) {
            for ( size_t j=0; j<=ncy; ++j ) {
                for ( size_t i=0; i<=ncx; ++i ) {
                    if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                        update_node(i, j, k, threadNo);
                    }
                }
            }
        }
        // sweep second direction
        for ( size_t k=0; k<=ncz; ++k ) {
            for ( size_t j=0; j<=ncy; ++j ) {
                for ( long int i=ncx; i>=0; --i ) {
                    if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                        update_node(i, j, k, threadNo);
                    }
                }
            }
        }
        // sweep third direction
        for ( size_t k=0; k<=ncz; ++k ) {
            for ( long int j=ncy; j>=0; --j ) {
                for ( size_t i=0; i<=ncx; ++i ) {
                    if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                        update_node(i, j, k, threadNo);
                    }
                }
            }
        }
        // sweep fourth direction
        for ( size_t k=0; k<=ncz; ++k ) {
            for ( long int j=ncy; j>=0; --j ) {
                for ( long int i=ncx; i>=0; --i ) {
                    if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                        update_node(i, j, k, threadNo);
                    }
                }
            }
        }
        // sweep fifth direction
        for ( long int k=ncz; k>=0; --k ) {
            for ( size_t j=0; j<=ncy; ++j ) {
                for ( size_t i=0; i<=ncx; ++i ) {
                    if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                        update_node(i, j, k, threadNo);
                    }
                }
            }
        }
        // sweep sixth direction
        for ( long int k=ncz; k>=0; --k ) {
            for ( size_t j=0; j<=ncy; ++j ) {
                for ( long int i=ncx; i>=0; --i ) {
                    if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                        update_node(i, j, k, threadNo);
                    }
                }
            }
        }
        // sweep seventh direction
        for ( long int k=ncz; k>=0; --k ) {
            for ( long int j=ncy; j>=0; --j ) {
                for ( size_t i=0; i<=ncx; ++i ) {
                    if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                        update_node(i, j, k, threadNo);
                    }
                }
            }
        }
        // sweep eighth direction
        for ( long int k=ncz; k>=0; --k ) {
            for ( long int j=ncy; j>=0; --j ) {
                for ( long int i=ncx; i>=0; --i ) {
                    if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                        update_node(i, j, k, threadNo);
                    }
                }
            }
        }
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::update_node(const size_t i, const size_t j, const size_t k,
                                           const size_t threadNo) const {
        T1 a1, a2, a3, t;

        if (k==0)
            a1 = nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
        else if (k==ncz)
            a1 = nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
        else {
            a1 = nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
            t  = nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
            a1 = a1<t ? a1 : t;
        }

        if (j==0)
            a2 = nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo);
        else if (j==ncy)
            a2 = nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo);
        else {
            a2 = nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo);
            t  = nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo);
            a2 = a2<t ? a2 : t;
        }

        if (i==0)
            a3 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo);
        else if (i==ncx)
            a3 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo);
        else {
            a3 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo);
            t  = nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo);
            a3 = a3<t ? a3 : t;
        }

        if ( a1>a2 ) std::swap(a1, a2);
        if ( a1>a3 ) std::swap(a1, a3);
        if ( a2>a3 ) std::swap(a2, a3);

        T1 fh = nodes[(k*(ncy+1)+j)*(ncx+1)+i].getNodeSlowness() * dx;

        t = a1 + fh;
        if ( t > a2 ) {

            t = 0.5*(a1+a2+sqrt(2.*fh*fh - (a1-a2)*(a1-a2)));

            if ( t > a3 ) {

                t = 1./3. * ((a1 + a2 + a3) + sqrt(-2.*a1*a1 + 2.*a1*a2 - 2.*a2*a2 +
                                                   2.*a1*a3 + 2.*a2*a3 -
                                                   2.*a3*a3 + 3.*fh*fh));

            }
        }

        if ( t<nodes[(k*(ncy+1)+j)*(ncx+1)+i].getTT(threadNo) )
            nodes[(k*(ncy+1)+j)*(ncx+1)+i].setTT(t,threadNo);

    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::sweep_weno3(const std::vector<bool>& frozen,
                                           const size_t threadNo) const {
        // sweep first direction
        for ( size_t k=0; k<=ncz; ++k ) {
            for ( size_t j=0; j<=ncy; ++j ) {
                for ( size_t i=0; i<=ncx; ++i ) {
                    if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                        update_node_weno3(i, j, k, threadNo);
                    }
                }
            }
        }
        // sweep second direction
        for ( size_t k=0; k<=ncz; ++k ) {
            for ( size_t j=0; j<=ncy; ++j ) {
                for ( long int i=ncx; i>=0; --i ) {
                    if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                        update_node_weno3(i, j, k, threadNo);
                    }
                }
            }
        }
        // sweep third direction
        for ( size_t k=0; k<=ncz; ++k ) {
            for ( long int j=ncy; j>=0; --j ) {
                for ( size_t i=0; i<=ncx; ++i ) {
                    if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                        update_node_weno3(i, j, k, threadNo);
                    }
                }
            }
        }
        // sweep fourth direction
        for ( size_t k=0; k<=ncz; ++k ) {
            for ( long int j=ncy; j>=0; --j ) {
                for ( long int i=ncx; i>=0; --i ) {
                    if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                        update_node_weno3(i, j, k, threadNo);
                    }
                }
            }
        }
        // sweep fifth direction
        for ( long int k=ncz; k>=0; --k ) {
            for ( size_t j=0; j<=ncy; ++j ) {
                for ( size_t i=0; i<=ncx; ++i ) {
                    if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                        update_node_weno3(i, j, k, threadNo);
                    }
                }
            }
        }
        // sweep sixth direction
        for ( long int k=ncz; k>=0; --k ) {
            for ( size_t j=0; j<=ncy; ++j ) {
                for ( long int i=ncx; i>=0; --i ) {
                    if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                        update_node_weno3(i, j, k, threadNo);
                    }
                }
            }
        }
        // sweep seventh direction
        for ( long int k=ncz; k>=0; --k ) {
            for ( long int j=ncy; j>=0; --j ) {
                for ( size_t i=0; i<=ncx; ++i ) {
                    if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                        update_node_weno3(i, j, k, threadNo);
                    }
                }
            }
        }
        // sweep eighth direction
        for ( long int k=ncz; k>=0; --k ) {
            for ( long int j=ncy; j>=0; --j ) {
                for ( long int i=ncx; i>=0; --i ) {
                    if ( !frozen[ (k*(ncy+1)+j)*(ncx+1)+i ] ) {
                        update_node_weno3(i, j, k, threadNo);
                    }
                }
            }
        }
    }

    template<typename T1, typename T2, typename NODE>
    T1 Grid3Drn<T1,T2,NODE>::weno3_upwind(const T1 v0, const T1 v1, const T1 v2, const T1 v3,
                                          const T1 v4, const T1 dx, bool forward) const {
    
        const T1 eps = std::numeric_limits<T1>::epsilon();
        
        if (forward) {
            // Forward differencing: ap = d/dx approximation
            const T1 num = (v4 - 2.0 * v3 + v2);
            const T1 den = (v3 - 2.0 * v2 + v1);
            const T1 r = (eps + num * num) / (eps + den * den);
            const T1 w = 1.0 / (1.0 + 2.0 * r * r);
            
            const T1 ap = (1.0 - w) * (v3 - v1) / (2.0 * dx) +
                             w * (-v4 + 4.0 * v3 - 3.0 * v2) / (2.0 * dx);
            
            return v2 + dx * ap;
        } else {
            // Backward differencing: am = -d/dx approximation
            const T1 num = (v2 - 2.0 * v1 + v0);
            const T1 den = (v3 - 2.0 * v2 + v1);
            const T1 r = (eps + num * num) / (eps + den * den);
            const T1 w = 1.0 / (1.0 + 2.0 * r * r);
            
            const T1 am = (1.0 - w) * (v3 - v1) / (2.0 * dx) +
                             w * (3.0 * v2 - 4.0 * v1 + v0) / (2.0 * dx);
            
            return v2 - dx * am;
        }
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::update_node_weno3(const size_t i,
                                                 const size_t j,
                                                 const size_t k,
                                                 const size_t threadNo) const {
        T1 a1, a2, a3, t;

        // ========== K direction ==========
        if (k==0) {
            a1 = nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);  // first order
        } else if (k==1) {
            a1 = weno3_upwind(0.0,  // v0 not used forward
                              nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo),
                              nodes[ (    k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo),
                              nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo),
                              nodes[ ((k+2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo),
                              dx, true);
            
            /* OLD CODE:
            T1 num = nodes[ ((k+2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
                  2.*nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
                     nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
                  2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
                     nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 ap = (1.-w)*(nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo)-
                            nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx) +
            w*(-nodes[ ((k+2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
               4.*nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
               3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx);

            a1 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) + dx*ap;
             */
            
            t = nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo); // first order for left
            a1 = a1<t ? a1 : t;

        } else if (k==ncz) {
            a1 = nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
        } else if (k==ncz-1) {
            a1 = weno3_upwind(nodes[ ((k-2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo),
                              nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo),
                              nodes[ (    k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo),
                              nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo),
                              0.0, // v4 not used backward
                              dx, false);
            
            /* OLD CODE:
            T1 num = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
            2.*nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
            nodes[ ((k-2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
            2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
            nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 am = (1.-w)*(nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo)-
                            nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx) +
            w*(3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
               4.*nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
               nodes[ ((k-2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx);

            a1 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) - dx*am;
             */

            t = nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo); // first order for right
            a1 = a1<t ? a1 : t;

        } else {
            // Forward direction
            a1 = weno3_upwind(nodes[ ((k-2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo),
                              nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo),
                              nodes[ (    k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo),
                              nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo),
                              nodes[ ((k+2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo),
                              dx, true);
            
            // Backward direction
            t = weno3_upwind(nodes[ ((k-2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo),
                             nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo),
                             nodes[ (    k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo),
                             nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo),
                             nodes[ ((k+2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo),
                             dx, false);
            
            /* OLD CODE:
            T1 num = nodes[ ((k+2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
                  2.*nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
                     nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
            2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
            nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 ap = (1.-w)*(nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo)-
                            nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx) +
            w*(-nodes[ ((k+2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
               4.*nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
               3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx);

            a1 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) + dx*ap;

            num = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
            2.*nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
            nodes[ ((k-2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
            num *= num;
            r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            w = 1./(1.+2.*r*r);

            T1 am = (1.-w)*(nodes[ ((k+1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo)-
                            nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx) +
            w*(3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
               4.*nodes[ ((k-1)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
               nodes[ ((k-2)*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx);

            t = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) - dx*am;
             */

            a1 = a1<t ? a1 : t;

        }

        // ========== J direction ==========
        if (j==0) {
            a2 = nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo);
        } else if (j==1) {
            a2 = weno3_upwind(0.0, // v0 not used forward
                              nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo),
                              nodes[ (k*(ncy+1)+j  )*(ncx+1)+i ].getTT(threadNo),
                              nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo),
                              nodes[ (k*(ncy+1)+j+2)*(ncx+1)+i ].getTT(threadNo),
                              dx, true);
            
            /* OLD CODE:
            T1 num = nodes[ (k*(ncy+1)+j+2)*(ncx+1)+i ].getTT(threadNo) -
            2.*nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo) +
            nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo) -
            2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
            nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 ap = (1.-w)*(nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo)-
                            nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo))/(2.*dx) +
            w*(-nodes[ (k*(ncy+1)+j+2)*(ncx+1)+i ].getTT(threadNo) +
               4.*nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo) -
               3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx);

            a2 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) + dx*ap;
             */

            t = nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo); // first order for left
            a2 = a2<t ? a2 : t;

        } else if (j==ncy) {
            a2 = nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo);
        } else if (j==ncy-1) {
            a2 = weno3_upwind(nodes[ (k*(ncy+1)+j-2)*(ncx+1)+i ].getTT(threadNo),
                              nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo),
                              nodes[ (k*(ncy+1)+j  )*(ncx+1)+i ].getTT(threadNo),
                              nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo),
                              0.0, // v4 not used backward
                              dx, false);
                        
            /* OLD CODE:
            T1 num = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
            2.*nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo) +
            nodes[ (k*(ncy+1)+j-2)*(ncx+1)+i ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo) -
            2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
            nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 am = (1.-w)*(nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo)-
                            nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo))/(2.*dx) +
            w*(3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
               4.*nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo) +
               nodes[ (k*(ncy+1)+j-2)*(ncx+1)+i ].getTT(threadNo))/(2.*dx);

            a2 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) - dx*am;
             */

            t = nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo); // first order for right
            a2 = a2<t ? a2 : t;

        } else {
            // Forward direction
            a2 = weno3_upwind(nodes[ (k*(ncy+1)+j-2)*(ncx+1)+i ].getTT(threadNo),
                              nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo),
                              nodes[ (k*(ncy+1)+j  )*(ncx+1)+i ].getTT(threadNo),
                              nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo),
                              nodes[ (k*(ncy+1)+j+2)*(ncx+1)+i ].getTT(threadNo),
                              dx, true);
            
            // Backward direction
            t = weno3_upwind(nodes[ (k*(ncy+1)+j-2)*(ncx+1)+i ].getTT(threadNo),
                             nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo),
                             nodes[ (k*(ncy+1)+j  )*(ncx+1)+i ].getTT(threadNo),
                             nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo),
                             nodes[ (k*(ncy+1)+j+2)*(ncx+1)+i ].getTT(threadNo),
                             dx, false);
            
            /* OLD CODE:
            T1 num = nodes[ (k*(ncy+1)+j+2)*(ncx+1)+i ].getTT(threadNo) -
            2.*nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo) +
            nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo) -
            2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
            nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 ap = (1.-w)*(nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo)-
                            nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo))/(2.*dx) +
            w*(-nodes[ (k*(ncy+1)+j+2)*(ncx+1)+i ].getTT(threadNo) +
               4.*nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo) -
               3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx);

            a2 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) + dx*ap;

            num = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
            2.*nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo) +
            nodes[ (k*(ncy+1)+j-2)*(ncx+1)+i ].getTT(threadNo);
            num *= num;
            r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            w = 1./(1.+2.*r*r);

            T1 am = (1.-w)*(nodes[ (k*(ncy+1)+j+1)*(ncx+1)+i ].getTT(threadNo)-
                            nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo))/(2.*dx) +
            w*(3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
               4.*nodes[ (k*(ncy+1)+j-1)*(ncx+1)+i ].getTT(threadNo) +
               nodes[ (k*(ncy+1)+j-2)*(ncx+1)+i ].getTT(threadNo))/(2.*dx);

            t = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) - dx*am;
             */

            a2 = a2<t ? a2 : t;

        }

        // ========== I direction ==========
        if (i==0) {
            a3 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo);
        } else if (i==1) {
            a3 = weno3_upwind(0.0, // v0 not used forward
                              nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo),
                              nodes[ (k*(ncy+1)+j)*(ncx+1)+i   ].getTT(threadNo),
                              nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo),
                              nodes[ (k*(ncy+1)+j)*(ncx+1)+i+2 ].getTT(threadNo),
                              dx, true);
                        
            /* OLD CODE:
            T1 num = nodes[ (k*(ncy+1)+j)*(ncx+1)+i+2 ].getTT(threadNo) -
            2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo) +
            nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo) -
            2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
            nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 ap = (1.-w)*(nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo)-
                            nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo))/(2.*dx) +
            w*(-nodes[ (k*(ncy+1)+j)*(ncx+1)+i+2 ].getTT(threadNo) +
               4.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo) -
               3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx);

            a3 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) + dx*ap;
             */

            t = nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo); // first order for left
            a3 = a3<t ? a3 : t;

        } else if (i==ncx) {
            a3 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo);
        } else if (i==ncx-1) {
            a3 = weno3_upwind(nodes[ (k*(ncy+1)+j)*(ncx+1)+i-2 ].getTT(threadNo),
                              nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo),
                              nodes[ (k*(ncy+1)+j)*(ncx+1)+i   ].getTT(threadNo),
                              nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo),
                              0.0, // v4 not used backward
                              dx, false);
                        
            /* OLD CODE:
            T1 num = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
            2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo) +
            nodes[ (k*(ncy+1)+j)*(ncx+1)+i-2 ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo) -
            2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
            nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 am = (1.-w)*(nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo)-
                            nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo))/(2.*dx) +
            w*(3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
               4.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo) +
               nodes[ (k*(ncy+1)+j)*(ncx+1)+i-2 ].getTT(threadNo))/(2.*dx);

            a3 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) - dx*am;
             */

            t = nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo); // first order for right
            a3 = a3<t ? a3 : t;

        } else {
            // Forward direction
            a3 = weno3_upwind(nodes[ (k*(ncy+1)+j)*(ncx+1)+i-2 ].getTT(threadNo),
                              nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo),
                              nodes[ (k*(ncy+1)+j)*(ncx+1)+i ]  .getTT(threadNo),
                              nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo),
                              nodes[ (k*(ncy+1)+j)*(ncx+1)+i+2 ].getTT(threadNo),
                              dx, true);
            
            // Backward direction
            t = weno3_upwind(nodes[ (k*(ncy+1)+j)*(ncx+1)+i-2 ].getTT(threadNo),
                             nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo),
                             nodes[ (k*(ncy+1)+j)*(ncx+1)+i   ].getTT(threadNo),
                             nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo),
                             nodes[ (k*(ncy+1)+j)*(ncx+1)+i+2 ].getTT(threadNo),
                             dx, false);

            /* OLD CODE:
            T1 num = nodes[ (k*(ncy+1)+j)*(ncx+1)+i+2 ].getTT(threadNo) -
            2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo) +
            nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo);
            num *= num;
            T1 den = nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo) -
            2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) +
            nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo);
            den *= den;
            T1 r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            T1 w = 1./(1.+2.*r*r);

            T1 ap = (1.-w)*(nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo)-
                            nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo))/(2.*dx) +
            w*(-nodes[ (k*(ncy+1)+j)*(ncx+1)+i+2 ].getTT(threadNo) +
               4.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo) -
               3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo))/(2.*dx);

            a3 = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) + dx*ap;

            num = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
            2.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo) +
            nodes[ (k*(ncy+1)+j)*(ncx+1)+i-2 ].getTT(threadNo);
            num *= num;
            r = (std::numeric_limits<T1>::epsilon()+num)/(std::numeric_limits<T1>::epsilon()+den);
            w = 1./(1.+2.*r*r);

            T1 am = (1.-w)*(nodes[ (k*(ncy+1)+j)*(ncx+1)+i+1 ].getTT(threadNo)-
                            nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo))/(2.*dx) +
            w*(3.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) -
               4.*nodes[ (k*(ncy+1)+j)*(ncx+1)+i-1 ].getTT(threadNo) +
               nodes[ (k*(ncy+1)+j)*(ncx+1)+i-2 ].getTT(threadNo))/(2.*dx);

            t = nodes[ (k*(ncy+1)+j)*(ncx+1)+i ].getTT(threadNo) - dx*am;
             */

            a3 = a3<t ? a3 : t;
        }

        if ( a1>a2 ) std::swap(a1, a2);
        if ( a1>a3 ) std::swap(a1, a3);
        if ( a2>a3 ) std::swap(a2, a3);

        T1 fh = nodes[(k*(ncy+1)+j)*(ncx+1)+i].getNodeSlowness() * dx;

        t = a1 + fh;
        if ( t > a2 ) {

            t = 0.5*(a1+a2+sqrt(2.*fh*fh - (a1-a2)*(a1-a2)));

            if ( t > a3 ) {

                t = 1./3. * ((a1 + a2 + a3) + sqrt(-2.*a1*a1 + 2.*a1*a2 -
                                                   2.*a2*a2 + 2.*a1*a3 + 2.*a2*a3 -
                                                   2.*a3*a3 + 3.*fh*fh));

            }
        }

        if ( t<nodes[(k*(ncy+1)+j)*(ncx+1)+i].getTT(threadNo) )
            nodes[(k*(ncy+1)+j)*(ncx+1)+i].setTT(t,threadNo);

    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::initFSM(const std::vector<sxyz<T1>>& Tx,
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

                    ptrdiff_t k = nn/((ncy+1)*(ncx+1));
                    ptrdiff_t j = (nn-k*(ncy+1)*(ncx+1))/(ncx+1);
                    ptrdiff_t i = nn - (k*(ncy+1)+j)*(ncx+1);

                    for ( ptrdiff_t kk=k-npts; kk<=k+npts; ++kk ) {
                        if ( kk>=0 && kk<=ncz ) {
                            for ( ptrdiff_t jj=j-npts; jj<=j+npts; ++jj ) {
                                if ( jj>=0 && jj<=ncy ) {
                                    for ( ptrdiff_t ii=i-npts; ii<=i+npts; ++ii ) {
                                        if ( ii>=0 && ii<=ncx && !(ii==i && jj==j && kk==k) ) {

                                            size_t nnn = (kk*(ncy+1)+jj)*(ncx+1)+ii;
                                            T1 tt = t0[n] + nodes[nnn].getDistance(Tx[n]) * nodes[nnn].getNodeSlowness();
                                            nodes[nnn].setTT( tt, threadNo );
                                            frozen[nnn] = true;

                                        }
                                    }
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

                ptrdiff_t k = cellNo/(ncy*ncx);
                ptrdiff_t j = (cellNo-k*ncy*ncx)/ncx;
                ptrdiff_t i = cellNo - (k*ncy+j)*ncx;

                for ( ptrdiff_t kk=k-(npts-1); kk<=k+npts; ++kk ) {
                    if ( kk>=0 && kk<=ncz ) {
                        for ( ptrdiff_t jj=j-(npts-1); jj<=j+npts; ++jj ) {
                            if ( jj>=0 && jj<=ncy ) {
                                for ( ptrdiff_t ii=i-(npts-1); ii<=i+npts; ++ii ) {
                                    if ( ii>=0 && ii<=ncx && !(ii==i && jj==j && kk==k) ) {

                                        size_t nnn = (kk*(ncy+1)+jj)*(ncx+1)+ii;
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
        }
    }

#ifdef VTK
    template<typename T1, typename T2, typename NODE>
    void Grid3Drn<T1,T2,NODE>::saveModelVTR(const std::string &fname,
                                            const bool saveSlowness) const {

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

        newScalars->SetNumberOfComponents(1);
        newScalars->SetNumberOfTuples( rgrid->GetNumberOfPoints() );

        if ( saveSlowness ) {
            newScalars->SetName("Slowness");
            for ( size_t n=0; n<nodes.size(); ++n ) {
                if ( nodes[n].isPrimary() ) {
                    double x[3] = {nodes[n].getX(), nodes[n].getY(), nodes[n].getZ()};
                    vtkIdType id = rgrid->FindPoint(x);
                    newScalars->SetTuple1(id, nodes[n].getNodeSlowness() );
                }
            }
        } else {
            newScalars->SetName("Velocity");
            for ( size_t n=0; n<nodes.size(); ++n ) {
                if ( nodes[n].isPrimary() ) {
                    double x[3] = {nodes[n].getX(), nodes[n].getY(), nodes[n].getZ()};
                    vtkIdType id = rgrid->FindPoint(x);
                    newScalars->SetTuple1(id, 1./nodes[n].getNodeSlowness() );
                }
            }
        }
        rgrid->GetPointData()->SetScalars(newScalars);

        vtkSmartPointer<vtkXMLRectilinearGridWriter> writer =
        vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();

        writer->SetFileName( fname.c_str() );
        writer->SetInputData( rgrid );
        writer->SetDataModeToBinary();
        writer->Update();
    }
#endif


}

#endif

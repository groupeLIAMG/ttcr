//
//  Grid2Duc.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-02-24.
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
 * @file Grid2Duc.h
 * @brief Base class for 2-D unstructured triangular meshes with cell-based
 *        slowness.
 *
 * Declares ttcr::Grid2Duc, the mid-level base shared by every 2-D unstructured
 * solver whose slowness is constant within each triangle: Grid2Ducsp
 * (shortest path), Grid2Ducfs (fast sweeping), Grid2Ducfm (fast marching) and
 * Grid2Ducdsp (dynamic shortest path).
 *
 * It is the unstructured counterpart of ttcr::Grid2Drc. The @c u in the name
 * means the mesh is **unstructured** — an arbitrary triangulation rather than a
 * regular grid — and @c c that slowness is held **per cell**, as opposed to
 * @c n (per node, see Grid2Dun.h).
 *
 * @section g2duc_vs_rect What changes on an unstructured mesh
 * Losing the regular grid costs the conveniences the rectilinear classes rely
 * on, and each has a replacement here:
 *
 * - **Point location** has no index arithmetic. ttcr::Grid2Duc::getCellNo first
 *   asks a kd-tree for the nearest node and tests only the triangles incident to
 *   it, falling back to an exhaustive scan if that misses.
 * - **Cell geometry** is not uniform, so each triangle caches its interior
 *   angles and edge lengths in a ttcr::triangleElemAngle rather than recomputing
 *   them from @c dx and @c dz.
 * - **Obtuse triangles** break the local update; see @ref g2duc_obtuse.
 *
 * @section g2duc_obtuse Obtuse triangles and virtual nodes
 * The local eikonal update assumes the wavefront reaches a vertex from within
 * the angular sector spanned by its two already-known neighbours. When the angle
 * at that vertex exceeds 90 degrees the assumption fails: the true first arrival
 * can come from outside the triangle, and a naive update violates causality and
 * overestimates the traveltime.
 *
 * ttcr::Grid2Duc::processObtuse repairs this once, at construction. For each
 * obtuse angle it finds the triangle across the opposite edge and builds a
 * ttcr::virtualNode — a support point obtained by unfolding into that neighbour,
 * so the update can be posed over a wider, non-obtuse sector. The results are
 * cached in ttcr::Grid2Duc::virtualNodes, keyed by triangle, and
 * ttcr::Grid2Duc::localSolver consults them whenever it meets an obtuse angle.
 *
 * Triangles on the boundary of the domain have no opposite neighbour to unfold
 * into, so **no correction is applied there** and the traveltime at such a
 * vertex keeps the obtuse-angle error.
 *
 * @sa Grid2D.h, Grid2Drc.h, Grid2Dun.h, Cell.h, Grid2Ducsp.h
 */

#ifndef ttcr_Grid2Duc_h
#define ttcr_Grid2Duc_h

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <memory>
#include <memory_resource>
#include <mutex>
#include <set>
#include <stdexcept>
#include <vector>

#ifdef VTK
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdocumentation"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkProbeFilter.h"
#include "vtkRectilinearGrid.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSmartPointer.h"
#include "vtkTriangle.h"
#include "vtkXMLRectilinearGridWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"
#pragma clang diagnostic pop
#endif

#include <boost/math/special_functions/sign.hpp>

#include "Grid2D.h"
#include "Grad.h"
#include "NodeKDTree2D.h"
#include "utils.h"

namespace ttcr {

    /**
     * @brief 2-D triangular mesh holding one slowness value per triangle.
     *
     * @tparam T1   floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2   integer type of node and cell indices, normally @c uint32_t.
     * @tparam S    point type: @ref sxz for a planar mesh, @ref sxyz for a
     *              draped mesh embedded in 3-D (the @c 2Ds programs).
     * @tparam NODE node type, e.g. ttcr::Node2Dc or ttcr::Node2Dcsp.
     * @tparam CELL cell policy supplying the slowness model, from Cell.h —
     *              what makes the class isotropic or anisotropic.
     *
     * @note Abstract in practice: it implements everything except @c raytrace.
     */
    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    class Grid2Duc : public Grid2D<T1,T2,S> {
    public:
        /**
         * @brief Build the mesh from a node list and a triangle list.
         *
         * @param no   node coordinates; these become the primary nodes, in the
         *             order given.
         * @param tri  triangles, each naming three node indices into @p no.
         * @param ttrp recompute receiver traveltimes by integrating slowness
         *             along the traced raypath.
         * @param nt   number of threads; sizes each node's traveltime array.
         *
         * @post Nodes and cells are allocated and the triangles copied into
         *       @ref triangles as ttcr::triangleElemAngle. Their angles and edge
         *       lengths are **not** yet computed, nor are the node ownership
         *       lists or the obtuse-angle corrections — the derived class calls
         *       @ref buildGridNodes and @ref processObtuse. Slowness is not set.
         * @note Unlike the rectilinear grids there is no geometry to derive:
         *       the mesh is supplied wholesale, so there are no cell-size or
         *       origin parameters.
         */
        Grid2Duc(const std::vector<S>& no,
                 const std::vector<triangleElem<T2>>& tri, const bool ttrp,
                 const size_t nt=1) :
        Grid2D<T1,T2,S>(tri.size(), ttrp, nt),
        nPrimary(static_cast<T2>(no.size())),
        nodes(std::vector<NODE>(no.size(), NODE(nt))),
        triangles(), virtualNodes(), cells(tri.size())
        {
            for (auto it=tri.begin(); it!=tri.end(); ++it) {
                triangles.push_back( *it );
            }
        }

        /// Destructor.
        virtual ~Grid2Duc() {}

        /**
         * @brief Set the slowness of every triangle.
         * @param s one slowness per triangle, in the order the triangles were
         *          supplied to the constructor.
         * @throws std::length_error if @p s does not match the triangle count,
         *         plus whatever the @c CELL policy raises.
         * @note The @c try / @c catch rethrows unchanged and so has no effect;
         *       the same idiom is repeated in the setters below.
         */
        void setSlowness(const std::vector<T1>& s) {
            try {
                cells.setSlowness( s );
            } catch (std::exception& e) {
                throw;
            }
        }
        /**
         * @brief Set the elliptical anisotropy ratio @f$\xi@f$ of every triangle.
         * @param x one value per triangle.
         * @throws std::exception if the @c CELL policy does not model @f$\xi@f$.
         *         Whether these anisotropy setters do anything depends entirely
         *         on which policy was supplied — see Cell.h.
         */
        void setXi(const std::vector<T1>& x) {
            try {
                cells.setXi( x );
            } catch (std::exception& e) {
                throw;
            }
        }
        /**
         * @brief Set the symmetry-axis tilt angle of every triangle.
         * @param t one angle per triangle, in radians.
         * @throws std::exception if the @c CELL policy is not tilted.
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
         * @param s one value per triangle.
         * @throws std::exception if the @c CELL policy is not VTI.
         */
        /**
         * @brief Select the phase to model in a transversely isotropic medium.
         * @param p 1 for the qP wave, any other value for the qSV wave.
         * @throws std::exception if the @c CELL policy describes no such choice.
         */
        void setPhase(const int p) {
            try {
                cells.setPhase(p);
            } catch (std::exception& e) {
                throw;
            }
        }
        void setVp0(const std::vector<T1>& s) {
            try {
                cells.setVp0(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        /**
         * @brief Set the S-wave velocity along the symmetry axis, @f$V_{S0}@f$.
         * @param s one value per triangle.
         * @throws std::exception if the @c CELL policy is not VTI.
         */
        void setVs0(const std::vector<T1>& s) {
            try {
                cells.setVs0(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        /**
         * @brief Set Thomsen's @f$\delta@f$ for every triangle.
         * @param s one value per triangle.
         * @throws std::exception if the @c CELL policy is not VTI.
         */
        void setDelta(const std::vector<T1>& s) {
            try {
                cells.setDelta(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        /**
         * @brief Set Thomsen's @f$\epsilon@f$ for every triangle.
         * @param s one value per triangle.
         * @throws std::exception if the @c CELL policy is not VTI.
         */
        void setEpsilon(const std::vector<T1>& s) {
            try {
                cells.setEpsilon(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        /**
         * @brief Set Thomsen's @f$\gamma@f$ for every triangle.
         * @param s one value per triangle.
         * @throws std::exception if the @c CELL policy does not model SH anisotropy.
         */
        void setGamma(const std::vector<T1>& s) {
            try {
                cells.setGamma(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        /**
         * @brief Set the weakly-anelliptical coefficient @f$s_2@f$ of every triangle.
         * @param s one value per triangle.
         * @throws std::exception if the @c CELL policy is not weakly anelliptical.
         */
        void setS2(const std::vector<T1>& s) {
            try {
                cells.setS2(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        /**
         * @brief Set the weakly-anelliptical coefficient @f$s_4@f$ of every triangle.
         * @param s one value per triangle.
         * @throws std::exception if the @c CELL policy is not weakly anelliptical.
         */
        void setS4(const std::vector<T1>& s) {
            try {
                cells.setS4(s);
            } catch (std::exception& e) {
                throw;
            }
        }
        /**
         * @brief Copy the triangle slowness values out of the mesh.
         * @param[out] s resized to the triangle count and filled in triangle
         *               order.
         */
        void getSlowness(std::vector<T1>& s) const {
            if (s.size() != triangles.size()) {
                s.resize(triangles.size());
            }
            for (size_t n=0; n<s.size(); ++n) {
                s[n] = cells.getSlowness(n);
            }
        }

        /**
         * @brief Write a traveltime directly into one node.
         * @param tt traveltime to store.
         * @param nn node index.
         * @param nt thread whose slot to write.
         * @warning Unchecked and non-@c const: it bypasses the solver entirely.
         *          Intended for seeding a boundary or source value; using it
         *          after a solve corrupts the field.
         */
        void setTT(const T1 tt, const size_t nn, const size_t nt=0) {
            nodes[nn].setTT(tt, nt);
        }

        /**
         * @brief Number of nodes in the mesh.
         * @param primary true to count only the primary nodes — those supplied
         *                to the constructor — false (the default) to include the
         *                secondary nodes a solver has added.
         * @return The requested count.
         */
        size_t getNumberOfNodes(const bool primary=false) const {
            if ( primary ) {
                return nPrimary;
            } else {
                return nodes.size();
            }
        }
        /// @return Number of triangles.
        size_t getNumberOfCells() const { return triangles.size(); }

        /**
         * @brief Collect the traveltimes computed at the primary nodes.
         * @param[out] tt       resized to @ref nPrimary and filled.
         * @param[in]  threadNo thread whose solution to read.
         */
        void getTT(std::vector<T1>& tt, const size_t threadNo=0) const final {
            tt.resize(nPrimary);
            for ( size_t n=0; n<nPrimary; ++n ) {
                tt[n] = nodes[n].getTT(threadNo);
            }
        }

        /**
         * @name Mesh bounding box
         *
         * Computed by scanning every node on each call — an unstructured mesh
         * has no stored extent, unlike the rectilinear grids where these are
         * plain members. Cache the result if you need it in a loop.
         * @{
         */
        const T1 getXmin() const {
            T1 xmin = nodes[0].getX();
            for ( auto it=nodes.begin(); it!=nodes.end(); ++it ) {
                xmin = xmin<it->getX() ? xmin : it->getX();
            }
            return xmin;
        }
        const T1 getXmax() const {
            T1 xmax = nodes[0].getX();
            for ( auto it=nodes.begin(); it!=nodes.end(); ++it ) {
                xmax = xmax>it->getX() ? xmax : it->getX();
            }
            return xmax;
        }
        const T1 getZmin() const {
            T1 zmin = nodes[0].getZ();
            for ( auto it=nodes.begin(); it!=nodes.end(); ++it ) {
                zmin = zmin<it->getZ() ? zmin : it->getZ();
            }
            return zmin;
        }
        const T1 getZmax() const {
            T1 zmax = nodes[0].getZ();
            for ( auto it=nodes.begin(); it!=nodes.end(); ++it ) {
                zmax = zmax>it->getZ() ? zmax : it->getZ();
            }
            return zmax;
        }
        /// @}

        /**
         * @brief Write the traveltime field to a file.
         * @param fname  output filename.
         * @param all    if nonzero, include the secondary nodes.
         * @param nt     thread whose solution to write.
         * @param format 1 for plain text, 2 for VTK, 3 for a raw binary dump.
         */
        void saveTT(const std::string &fname, const int all, const size_t nt=0,
                    const int format=1) const;

#ifdef VTK
        /**
         * @brief Write the mesh and its model as a VTK unstructured grid.
         * @param fname              output filename.
         * @param saveSlowness       write slowness rather than velocity.
         * @param savePhysicalEntity also write each triangle's physical-entity
         *                           tag, the reflector/region label carried over
         *                           from the Gmsh model.
         * @note Available only when built with @c VTK defined. The natural
         *       format here, since the mesh has no rectilinear structure.
         */
        void saveModelVTU(const std::string &fname, const bool saveSlowness=true,
                          const bool savePhysicalEntity=false) const;
        /**
         * @brief Resample the model onto a rectilinear grid and write it as VTK.
         * @param fname        output filename.
         * @param d            cell sizes of the target grid, as a 3-element array.
         * @param saveSlowness write slowness rather than velocity.
         * @note Lossy: the triangulation is sampled onto a regular grid, so
         *       sharp interfaces are stair-stepped. For viewing only.
         */
        void saveModelVTR(const std::string &fname, const double* d,
                          const bool saveSlowness=true) const;
#endif

        /**
         * @brief Compute the area of every triangle.
         * @param[out] area resized to the triangle count and filled.
         * @note Used to area-weight the cell-to-node interpolation, so a large
         *       triangle counts for more than a small one at a shared vertex.
         */
        void calculateArea(std::vector<T1> &area) const;
        /**
         * @brief Interpolate the cell slowness model onto the nodes, in place.
         * @param[out] s resized to the node count and filled with the
         *               area-weighted mean of the incident triangles' slowness.
         * @note The mesh keeps its per-cell model; this only produces a nodal
         *       view of it, for output or for a solver that needs one.
         */
        void interpolateAtNodes(std::vector<T1> &s) const;
        /**
         * @brief Interpolate an arbitrary per-cell field onto the nodes.
         * @param[in]  cellVal one value per triangle.
         * @param[out] nodeVal resized to the node count and filled with the
         *                     area-weighted mean at each node.
         */
        void interpolateAtNodes(const std::vector<T1> &cellVal,
                                std::vector<T1> &nodeVal) const;

        /**
         * @brief Write the coordinates of the secondary nodes to a stream.
         * @param os destination stream, one @c "x z" pair per line.
         * @note Writes nothing when the solver adds none.
         */
        void dump_secondary(std::ofstream& os) const {
            for ( size_t n=nPrimary; n<nodes.size(); ++n ) {
                os << nodes[n].getX() << ' ' << nodes[n].getZ() << '\n';
            }
        }

        /**
         * @brief Copy out the primary node coordinates.
         * @param[out] _nodes resized to @ref nPrimary and filled.
         * @note Secondary nodes are omitted, so the result matches the node list
         *       originally supplied to the constructor.
         */
        void getNodes(std::vector<S>& _nodes) const final {
            _nodes.resize(nPrimary);
            for ( size_t n=0; n<nPrimary; ++n ) {
                _nodes[n] = nodes[n];
            }
        }
        /**
         * @brief Copy out the triangle connectivity, as fixed-size arrays.
         * @param[out] tri resized to the triangle count; each entry holds the
         *                 three node indices.
         */
        void getTriangles(std::vector<std::array<T2, 3>>& tri) const final {
            tri.resize(triangles.size());
            for ( size_t n=0; n<triangles.size(); ++n ) {
                tri[n] = {triangles[n].i[0], triangles[n].i[1], triangles[n].i[2]};
            }
        }
        /**
         * @brief Copy out the triangle connectivity, as vectors.
         * @param[out] tri resized to the triangle count; each entry holds the
         *                 three node indices.
         * @note Same data as the @c std::array overload, in the shape the Cython
         *       bindings want.
         */
        void getTriangles(std::vector<std::vector<T2>>& tri) const final {
            tri.resize(triangles.size());
            for ( size_t n=0; n<triangles.size(); ++n ) {
                tri[n].resize(3);
                tri[n] = {triangles[n].i[0], triangles[n].i[1], triangles[n].i[2]};
            }
        }

        /**
         * @brief Mean edge length over the whole mesh.
         * @return The average, a natural length scale for a mesh with no
         *         uniform cell size — used to turn relative tolerances and
         *         radii into absolute distances.
         */
        const T1 getAverageEdgeLength() const;

    protected:
        /// Number of primary nodes, i.e. the mesh vertices supplied to the
        /// constructor. Nodes at or beyond this index are secondary.
        T2 nPrimary;
        /// Primary nodes first, then any secondary nodes the derived solver
        /// added. @c mutable because @c const raytracing methods update the
        /// traveltimes stored in them.
        mutable std::vector<NODE> nodes;
        /// The triangles, each caching its interior angles and edge lengths —
        /// the unstructured stand-in for the rectilinear grids' @c dx / @c dz.
        std::vector<triangleElemAngle<T1,T2>> triangles;
        /// Unfolded support points for obtuse triangles, keyed by triangle
        /// index. Empty for a mesh with no obtuse angles.
        /// @sa @ref g2duc_obtuse
        std::map<T2, virtualNode<T1,NODE>> virtualNodes;

        /// Slowness model, one value per triangle. Its type is the @c CELL
        /// policy, so this member also carries any anisotropy parameters.
        CELL cells;

        /**
         * @brief Position the primary nodes and build their ownership lists.
         * @param no node coordinates.
         * @param nt number of threads.
         * @post Every node knows which triangles are incident to it, which is
         *       what makes @ref getCellNo and @ref processObtuse possible.
         */
        void buildGridNodes(const std::vector<S>& no,
                            const size_t nt);

        /**
         * @brief Build the primary nodes plus secondary nodes along each edge.
         * @param no    node coordinates.
         * @param nsecondary number of secondary nodes per triangle edge.
         * @param nt    number of threads.
         * @note The overload used by the shortest-path solvers, which need the
         *       extra nodes for ray-direction resolution.
         */
        void buildGridNodes(const std::vector<S>& no,
                            const T2 nsecondary,
                            const size_t nt);

        // kd-tree over all nodes, for fast point location and on-node lookup.
        // Built lazily; node coordinates are fixed after construction and
        // queries are read-only, hence thread-safe.
        mutable std::unique_ptr<NodeKDTree2D<T1,T2>> kdtree;
        mutable std::once_flag kdtreeFlag;

        /**
         * @brief Index of the node closest to a point.
         * @param pt query point.
         * @return Nearest node index — always a valid one, whether or not the
         *         node actually coincides with @p pt.
         * @note Builds @ref kdtree on first call; @c std::call_once makes that
         *       safe under concurrency, and later queries need no
         *       synchronisation since node coordinates are then fixed.
         */
        T2 getNearestNode(const S& pt) const {
            std::call_once(kdtreeFlag, [this]() {
                kdtree.reset(new NodeKDTree2D<T1,T2>(nodes, nodes.size()));
            });
            return kdtree->findNearest(pt.x, pt.z);
        }

        /**
         * @brief Locate the triangle containing a point.
         * @param pt point to locate.
         * @return Index of the containing triangle.
         * @throws std::runtime_error naming the point if it lies outside the
         *         mesh — unlike the rectilinear grids, which return an
         *         out-of-range index silently.
         * @note Two-stage: the containing triangle is almost always incident to
         *       the nearest node, so those few are tested first via the kd-tree,
         *       and only a miss falls back to scanning every triangle. Point
         *       location has no closed form on an unstructured mesh.
         *       @sa @ref g2duc_vs_rect
         * @note Where several incident triangles match — a point exactly on a
         *       shared edge or vertex — the lowest index wins, making the result
         *       deterministic.
         */
        T2 getCellNo(const S& pt) const {
            // The containing triangle is, in the common case, incident to the
            // nearest node, so check that node's triangles first.
            T2 nn = getNearestNode(pt);
            T2 cell = std::numeric_limits<T2>::max();
            for ( T2 t : nodes[nn].getOwners() ) {
                if ( insideTriangle(pt, t) && t < cell ) {
                    cell = t;
                }
            }
            if ( cell != std::numeric_limits<T2>::max() ) {
                return cell;
            }
            // fall back to the exhaustive search
            for ( T2 n=0; n<triangles.size(); ++n ) {
                if ( insideTriangle(pt, n) ) {
                    return n;
                }
            }
            std::ostringstream msg;
            msg << "Point " << pt << " cannot be found in mesh.";
            throw std::runtime_error(msg.str());
        }

        /**
         * @brief Cached classification of the source points.
         *
         * Locating a source — deciding whether it sits on a node and which
         * triangle holds it — costs a node scan and a @ref getCellNo call.
         * Without caching that would repeat for every (source, receiver) pair;
         * with it, once per source.
         */
        // Cached classification of the Tx points (initTxVars), so the per-Tx
        // node scan + getCellNo runs once per source rather than once per
        // (source, receiver) pair.  Indexed by threadNo (each thread owns a
        // slot in parallel runs); recomputed only when its Tx changes.
        struct txInfo_t {
            std::vector<sxz<T1>> tx;            ///< The source positions this entry describes.
            std::vector<bool> txOnNode;         ///< Whether each source coincides with a node.
            std::vector<T2> txNode;             ///< Node index for the sources that do.
            std::vector<T2> txCell;             ///< Containing triangle for the sources that do not.
            std::vector<std::vector<T2>> txCells; ///< Triangles adjacent to each source.
            bool valid = false;                 ///< False until first populated.
        };
        /// One ttcr::Grid2Duc::txInfo_t per thread, so parallel runs with different sources
        /// do not invalidate each other's cache. @c mutable because the @c const
        /// raytracing methods fill it.
        mutable std::vector<txInfo_t> txInfoCache;

        /**
         * @brief Source classification for a thread, computing it if stale.
         * @param Tx       source positions.
         * @param threadNo thread whose cache slot to use.
         * @return The cached ttcr::Grid2Duc::txInfo_t, recomputed via @ref initTxVars only
         *         if this thread's entry is unset or describes different sources.
         * @note The cache is keyed on the @p Tx vector itself, compared
         *       element-wise, so moving a source invalidates it correctly.
         */
        const txInfo_t& getTxInfo(const std::vector<sxz<T1>>& Tx,
                                  const size_t threadNo) const {
            if ( txInfoCache.size() <= threadNo ) {
                txInfoCache.resize( threadNo + 1 );
            }
            txInfo_t& ti = txInfoCache[threadNo];
            if ( ti.valid && ti.tx == Tx ) {
                return ti;
            }
            ti.tx = Tx;
            ti.txOnNode.assign( Tx.size(), false );
            ti.txNode.assign( Tx.size(), 0 );
            ti.txCell.assign( Tx.size(), 0 );
            ti.txCells.assign( Tx.size(), std::vector<T2>() );
            initTxVars(Tx, ti.txOnNode, ti.txNode, ti.txCell, ti.txCells);
            ti.valid = true;
            return ti;
        }

        /**
         * @brief Traveltime at a receiver.
         * @param Rx       receiver position.
         * @param threadNo thread whose solution to read.
         * @return Traveltime, interpolated if @p Rx does not sit on a node.
         */
        T1 getTraveltime(const S& Rx,
                         const size_t threadNo) const;

        /**
         * @brief Traveltime at a receiver, also reporting where it came from.
         * @param[in]  Rx           receiver position.
         * @param[out] nodeParentRx node the ray arrived from.
         * @param[out] cellParentRx triangle it crossed to get there.
         * @param[in]  threadNo     thread whose solution to read.
         * @return Traveltime at @p Rx.
         * @note The two outputs start the walk back along a raypath.
         */
        T1 getTraveltime(const S& Rx,
                         T2& nodeParentRx,
                         T2& cellParentRx,
                         const size_t threadNo) const;


        /**
         * @brief Verify that every point lies inside the mesh.
         * @param pts points to check, typically sources or receivers.
         * @throws std::runtime_error naming the offending point.
         * @note Two overloads so a draped mesh can be checked with either point
         *       type.
         */
        void checkPts(const std::vector<sxz<T1>>& pts) const;
        /// @copydoc checkPts(const std::vector<sxz<T1>>&) const
        void checkPts(const std::vector<sxyz<T1>>& pts) const;

        /**
         * @brief Test whether a point lies inside a given triangle.
         * @param pt  point to test.
         * @param tri triangle index.
         * @return True if inside.
         * @note The primitive underlying @ref getCellNo.
         */
        bool insideTriangle(const sxz<T1>& pt, const T2 tri) const;
        /// @copydoc insideTriangle(const sxz<T1>&, const T2) const
        bool insideTriangle(const sxyz<T1>& pt, const T2 tri) const;

        /**
         * @brief Build the virtual nodes that correct obtuse triangles.
         *
         * Run once after the mesh is built. For each interior angle above 90
         * degrees it unfolds into the triangle across the opposite edge and
         * records a ttcr::virtualNode in ttcr::virtualNodes, keyed by triangle.
         *
         * @note Triangles on the domain boundary have no opposite neighbour, so
         *       no correction is stored for them and the obtuse-angle error
         *       remains there. @sa @ref g2duc_obtuse
         */
        void processObtuse();

        /**
         * @brief Solve the local eikonal update at one vertex.
         * @param vertexC  node to update.
         * @param threadNo thread whose traveltimes to read and write.
         * @note Consults ttcr::virtualNodes when the angle at @p vertexC is
         *       obtuse, and uses the triangle's own cached angles otherwise.
         */
        void localSolver(NODE *vertexC, const size_t threadNo) const;


        /**
         * @brief Classify the source points against the mesh.
         * @param[in]  Tx       source positions.
         * @param[out] txOnNode whether each source coincides with a node.
         * @param[out] txNode   node index for those that do.
         * @param[out] txCell   containing triangle for those that do not.
         * @param[out] txCells  triangles adjacent to each source.
         * @note Called through @ref getTxInfo, which caches the result per
         *       thread — do not call it directly in a loop.
         */
        void initTxVars(const std::vector<sxz<T1>>& Tx,
                        std::vector<bool>& txOnNode,
                        std::vector<T2>& txNode,
                        std::vector<T2>& txCell,
                        std::vector<std::vector<T2>>& txCells) const;

        void getRaypath(const std::vector<sxz<T1>>& Tx,
                        const sxz<T1> &Rx,
                        std::vector<sxz<T1>> &r_data,
                        const size_t threadNo) const;

        void getRaypath(const std::vector<sxz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const sxz<T1>& Rx,
                        std::vector<sxz<T1>>& r_data,
                        T1 &tt,
                        const size_t threadNo) const;

        void getRaypath(const std::vector<sxz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const sxz<T1>& Rx,
                        std::vector<sxz<T1>>& r_data,
                        std::vector<siv<T1>> &l_data,
                        T1 &tt,
                        const size_t threadNo) const;

        void getRaypath(const std::vector<sxz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const sxz<T1>& Rx,
                        std::vector<siv<T1>> &l_data,
                        T1 &tt,
                        const size_t threadNo) const;

        T1 getTraveltimeFromRaypath(const std::vector<sxz<T1>>& Tx,
                                    const std::vector<T1>& t0,
                                    const sxz<T1>& Rx,
                                    const size_t threadNo) const;

        bool findIntersection(const T2 i0, const T2 i1,
                              const sxz<T1> &g,
                              sxz<T1> &curr_pt) const;

        T2 findNextCell1(const T2 i0, const T2 i1, const T2 nodeNo) const;
        T2 findNextCell2(const T2 i0, const T2 i1, const T2 cellNo) const;

        template<typename SetT> void getNeighborNodes(const T2 cellNo, SetT &nnodes) const;

    };

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Duc<T1,T2,S,NODE,CELL>::buildGridNodes(const std::vector<S>& no,
                                                const size_t nt) {
        // primary nodes
        for ( T2 n=0; n<no.size(); ++n ) {
            nodes[n].setXYZindex( no[n], n );
            nodes[n].setPrimary(true);
        }
        for ( T2 ntri=0; ntri<triangles.size(); ++ntri ) {
            for ( size_t nl=0; nl<3; ++nl ) {
                // push owner for primary nodes
                nodes[ triangles[ntri].i[nl] ].pushOwner( ntri );
            }

            // distance between node 1 & 2 (opposite of node 0)
            T1 a = nodes[ triangles[ntri].i[1] ].getDistance( nodes[ triangles[ntri].i[2] ] );

            // distance between node 0 & 2 (opposite of node 1)
            T1 b = nodes[ triangles[ntri].i[0] ].getDistance( nodes[ triangles[ntri].i[2] ] );

            // distance between node 0 & 1 (opposite of node 2]
            T1 c = nodes[ triangles[ntri].i[0] ].getDistance( nodes[ triangles[ntri].i[1] ] );

            triangles[ntri].l[0] = a;
            triangles[ntri].l[1] = b;
            triangles[ntri].l[2] = c;

            // angle at node 0
            triangles[ntri].a[0] = acos((b*b + c*c - a*a)/(2.*b*c));

            // angle at node 1
            triangles[ntri].a[1] = acos((c*c + a*a - b*b)/(2.*a*c));

            // angle at node 2
            triangles[ntri].a[2] = acos((a*a + b*b - c*c)/(2.*a*b));

        }
    }


    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Duc<T1,T2,S,NODE,CELL>::buildGridNodes(const std::vector<S>& no,
                                                const T2 nsecondary,
                                                const size_t nt) {

        // primary nodes
        for ( T2 n=0; n<no.size(); ++n ) {
            nodes[n].setXYZindex( no[n], n );
            nodes[n].setPrimary(true);
        }
        T2 nNodes = static_cast<T2>(nodes.size());

        std::map<std::array<T2,2>,std::vector<T2>> lineMap;
        std::array<T2,2> lineKey;
        typename std::map<std::array<T2,2>,std::vector<T2>>::iterator lineIt;

        size_t estLineNo = (triangles.size()+triangles.size()/10) * 3/2;
        nodes.reserve( nNodes + estLineNo*nsecondary );

        // edge nodes
        NODE tmpNode(nt);
        for ( T2 ntri=0; ntri<triangles.size(); ++ntri ) {

            for ( size_t nl=0; nl<3; ++nl ) {

                // push owner for primary nodes
                nodes[ triangles[ntri].i[nl] ].pushOwner( ntri );

                if ( nsecondary>0 ) {

                    lineKey = { triangles[ntri].i[nl],
                        triangles[ntri].i[(nl+1)%3] };
                    std::sort(lineKey.begin(), lineKey.end());

                    lineIt = lineMap.find( lineKey );
                    if ( lineIt == lineMap.end() ) {
                        // not found, insert new pair
                        lineMap[ lineKey ] = std::vector<T2>(nsecondary);
                    } else {
                        for ( size_t n=0; n<lineIt->second.size(); ++n ) {
                            // setting owners
                            nodes[ lineIt->second[n] ].pushOwner( ntri );
                        }
                        continue;
                    }

                    S d = (no[lineKey[1]]-no[lineKey[0]])/static_cast<T1>(nsecondary+1);

                    for ( size_t n2=0; n2<nsecondary; ++n2 ) {
                        tmpNode.setXYZindex(no[lineKey[0]]+static_cast<T1>(1+n2)*d,
                                            nNodes );
                        lineMap[lineKey][n2] = nNodes++;
                        nodes.push_back( tmpNode );
                        nodes.back().pushOwner( ntri );
                    }
                }
            }
        }

        nodes.shrink_to_fit();
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    T1 Grid2Duc<T1,T2,S,NODE,CELL>::getTraveltime(const S& Rx,
                                             const size_t threadNo) const {

        T2 nn = getNearestNode( Rx );
        if ( nodes[nn] == Rx ) {
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


    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    T1 Grid2Duc<T1,T2,S,NODE,CELL>::getTraveltime(const S& Rx,
                                             T2& nodeParentRx,
                                             T2& cellParentRx,
                                             const size_t threadNo) const {

        T2 nn = getNearestNode( Rx );
        if ( nodes[nn] == Rx ) {
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



    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Duc<T1,T2,S,NODE,CELL>::checkPts(const std::vector<sxz<T1>>& pts) const {

        for (size_t n=0; n<pts.size(); ++n) {
            bool found = false;
            // check first if point is on a node
            for ( T2 nt=0; nt<nodes.size(); ++nt ) {
                if ( nodes[nt] == pts[n]) {
                    found = true;
                    break;
                }
            }
            if ( found == false ) {
                for ( T2 nt=0; nt<triangles.size(); ++nt ) {
                    if ( insideTriangle(pts[n], nt) ) {
                        found = true;
                    }
                }
            }
            if ( found == false ) {
                std::ostringstream msg;
                msg << "Error: Point no " << n << " (" << pts[n].x << ", "<< pts[n] .z << ") outside mesh.";
                throw std::runtime_error(msg.str());
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Duc<T1,T2,S,NODE,CELL>::checkPts(const std::vector<sxyz<T1>>& pts) const {

        for (size_t n=0; n<pts.size(); ++n) {
            bool found = false;
            // check first if point is on a node
            for ( T2 nt=0; nt<nodes.size(); ++nt ) {
                if ( nodes[nt] == pts[n]) {
                    found = true;
                    break;
                }
            }
            if ( found == false ) {
                for ( T2 nt=0; nt<triangles.size(); ++nt ) {
                    if ( insideTriangle(pts[n], nt) ) {
                        found = true;
                    }
                }
            }
            if ( found == false ) {
                std::ostringstream msg;
                msg << "Error: Point no " << n << " (" << pts[n].x << ", "<< pts[n] .y << ", "<< pts[n] .z << ") outside mesh.";
                throw std::runtime_error(msg.str());
            }
        }
    }


    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    bool Grid2Duc<T1,T2,S,NODE,CELL>::insideTriangle(const sxz<T1>& v, const T2 nt) const {

        if ( !testInTriangleBoundingBox(&(nodes[ triangles[nt].i[0] ]),
                                        &(nodes[ triangles[nt].i[1] ]),
                                        &(nodes[ triangles[nt].i[2] ]), v) ) {
            return false;
        }

        // from http://mathworld.wolfram.com/TriangleInterior.html

        sxz<T1> v0 = { nodes[ triangles[nt].i[0] ].getX(),
            nodes[ triangles[nt].i[0] ].getZ() };
        sxz<T1> v1 = { nodes[ triangles[nt].i[1] ].getX()-v0.x,
            nodes[ triangles[nt].i[1] ].getZ()-v0.z };
        sxz<T1> v2 = { nodes[ triangles[nt].i[2] ].getX()-v0.x,
            nodes[ triangles[nt].i[2] ].getZ()-v0.z };

        T1 invDenom = 1. / det(v1, v2);
        T1 a = (det(v, v2) - det(v0, v2)) * invDenom;
        T1 b = -(det(v, v1) - det(v0, v1)) * invDenom;
        return (a >= -small3) && (b >= -small3) && (a + b <= 1.+small3);
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    bool Grid2Duc<T1,T2,S,NODE,CELL>::insideTriangle(const sxyz<T1>& p, const T2 nt) const {

        if ( !testInTriangleBoundingBox(&(nodes[ triangles[nt].i[0] ]),
                                        &(nodes[ triangles[nt].i[1] ]),
                                        &(nodes[ triangles[nt].i[2] ]), p) ) {
            return false;
        }

        sxyz<T1> a = { nodes[ triangles[nt].i[0] ].getX(),
            nodes[ triangles[nt].i[0] ].getY(),
            nodes[ triangles[nt].i[0] ].getZ() };
        sxyz<T1> b = { nodes[ triangles[nt].i[1] ].getX(),
            nodes[ triangles[nt].i[1] ].getY(),
            nodes[ triangles[nt].i[1] ].getZ() };
        sxyz<T1> c = { nodes[ triangles[nt].i[2] ].getX(),
            nodes[ triangles[nt].i[2] ].getY(),
            nodes[ triangles[nt].i[2] ].getZ() };

        // Translate point and triangle so that point lies at origin
        a -= p; b -= p; c -= p;
        // Compute normal vectors for triangles pab and pbc
        sxyz<T1> u = cross(b, c);
        sxyz<T1> v = cross(c, a);
        // Make sure they are both pointing in the same direction
        if (dot(u, v) < static_cast<T1>(0.0)) return false;
        // Compute normal vector for triangle pca
        sxyz<T1> w = cross(a, b);
        // Make sure it points in the same direction as the first two
        if (dot(u, w) < static_cast<T1>(0.0)) return false;
        // Otherwise P must be in (or on) the triangle
        return true;
    }


    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Duc<T1,T2,S,NODE,CELL>::saveTT(const std::string &fname, const int all,
                                        const size_t nt, const int format) const {

        if ( format == 1 ) {
            std::string filename = fname+".dat";
            std::ofstream fout(filename.c_str());
            fout.precision(12);
            T2 nMax = nPrimary;
            if ( all == 1 ) {
                nMax = static_cast<T2>(nodes.size());
            }
            for ( T2 n=0; n<nMax; ++n ) {
                fout << nodes[n].getX() << '\t'
                << nodes[n].getZ() << '\t'
                << nodes[n].getTT(nt) << '\n';
            }
            fout.close();
        } else if ( format == 2) {
#ifdef VTK
            std::string filename = fname+".vtu";

            vtkSmartPointer<vtkUnstructuredGrid> ugrid =
            vtkSmartPointer<vtkUnstructuredGrid>::New();

            vtkSmartPointer<vtkPoints> newPts =
            vtkSmartPointer<vtkPoints>::New();
            vtkSmartPointer<vtkDoubleArray> newScalars =
            vtkSmartPointer<vtkDoubleArray>::New();

            newScalars->SetName("Travel time");

            double xyz[3];
            T2 nMax = nPrimary;  // only primary are saved
            for (size_t n=0; n<nMax; ++n) {
                xyz[0] = nodes[n].getX();
                xyz[1] = nodes[n].getY();
                xyz[2] = nodes[n].getZ();
                newPts->InsertPoint(n, xyz);
                newScalars->InsertValue(n, nodes[n].getTT(nt) );
            }

            ugrid->SetPoints(newPts);
            ugrid->GetPointData()->SetScalars(newScalars);

            vtkSmartPointer<vtkTriangle> tri =
            vtkSmartPointer<vtkTriangle>::New();
            for (size_t n=0; n<triangles.size(); ++n) {
                tri->GetPointIds()->SetId(0, triangles[n].i[0] );
                tri->GetPointIds()->SetId(1, triangles[n].i[1] );
                tri->GetPointIds()->SetId(2, triangles[n].i[2] );

                ugrid->InsertNextCell( tri->GetCellType(), tri->GetPointIds() );
            }
            vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
            vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

            writer->SetFileName( filename.c_str() );
            writer->SetInputData( ugrid );
            writer->SetDataModeToBinary();
            writer->Update();
#else
            std::cerr << "VTK not included during compilation.\nNothing saved.\n";
#endif
        } else if ( format == 3 ) {
            std::string filename = fname+".bin";
            std::ofstream fout(filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
            T2 nMax = nPrimary;
            if ( all == 1 ) {
                nMax = static_cast<T2>(nodes.size());
            }
            for ( T2 n=0; n<nMax; ++n ) {
                T1 tmp[] = { nodes[n].getX(), nodes[n].getZ(), nodes[n].getTT(nt) };
                fout.write( (char*)tmp, 3*sizeof(T1) );
            }
            fout.close();
        } else {
            throw std::runtime_error("Unsupported format for saving traveltimes");
        }
    }

#ifdef VTK

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Duc<T1,T2,S,NODE,CELL>::saveModelVTU(const std::string &fname,
                                              const bool saveSlowness,
                                              const bool savePhysicalEntity) const {

        vtkSmartPointer<vtkUnstructuredGrid> ugrid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();

        vtkSmartPointer<vtkPoints> newPts =
        vtkSmartPointer<vtkPoints>::New();

        double xyz[3];
        T2 nMax = nPrimary;  // only primary are saved
        for (size_t n=0; n<nMax; ++n) {
            xyz[0] = nodes[n].getX();
            xyz[1] = nodes[n].getY();
            xyz[2] = nodes[n].getZ();
            newPts->InsertPoint(n, xyz);
        }

        ugrid->SetPoints(newPts);

        vtkSmartPointer<vtkTriangle> tri =
        vtkSmartPointer<vtkTriangle>::New();
        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();
        if ( saveSlowness ) {
            data->SetName("Slowness");

            for (size_t n=0; n<triangles.size(); ++n) {
                tri->GetPointIds()->SetId(0, triangles[n].i[0] );
                tri->GetPointIds()->SetId(1, triangles[n].i[1] );
                tri->GetPointIds()->SetId(2, triangles[n].i[2] );

                ugrid->InsertNextCell( tri->GetCellType(), tri->GetPointIds() );
                data->InsertNextValue( cells.getSlowness(n) );
            }
        } else {
            data->SetName("Velocity");

            for (size_t n=0; n<triangles.size(); ++n) {
                tri->GetPointIds()->SetId(0, triangles[n].i[0] );
                tri->GetPointIds()->SetId(1, triangles[n].i[1] );
                tri->GetPointIds()->SetId(2, triangles[n].i[2] );

                ugrid->InsertNextCell( tri->GetCellType(), tri->GetPointIds() );
                data->InsertNextValue( 1./cells.getSlowness(n) );
            }
        }

        ugrid->GetCellData()->SetScalars(data);

        vtkSmartPointer<vtkIntArray> data_pe = vtkSmartPointer<vtkIntArray>::New();
        if ( savePhysicalEntity ) {
            data_pe->SetName("Physical entity");
            for (size_t n=0; n<triangles.size(); ++n) {
                data_pe->InsertNextValue(triangles[n].physical_entity );
            }
            ugrid->GetCellData()->AddArray(data_pe);
        }

        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

        writer->SetFileName( fname.c_str() );
        writer->SetInputData( ugrid );
        writer->SetDataModeToBinary();
        writer->Update();
    }


    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Duc<T1,T2,S,NODE,CELL>::saveModelVTR(const std::string &fname,
                                              const double *d,
                                              const bool saveSlowness) const {

        double x[] = { nodes[0].getX(), nodes[0].getX(),
            0.0, 0.0,
            nodes[0].getZ(), nodes[0].getZ() };

        for (size_t n=1; n<nPrimary; ++n) {
            x[0] = x[0] < nodes[n].getX() ? x[0] : nodes[n].getX();
            x[1] = x[1] > nodes[n].getX() ? x[1] : nodes[n].getX();

            x[4] = x[4] < nodes[n].getZ() ? x[4] : nodes[n].getZ();
            x[5] = x[5] > nodes[n].getZ() ? x[5] : nodes[n].getZ();
        }

        int nn[3];
        nn[0] = 1 + (x[1]-x[0])/d[0];
        nn[1] = 1;
        nn[2] = 1 + (x[5]-x[4])/d[2];

        vtkSmartPointer<vtkDoubleArray> xCoords = vtkSmartPointer<vtkDoubleArray>::New();
        for (size_t n=0; n<nn[0]; ++n) {
            xCoords->InsertNextValue( x[0] + n*d[0] );
        }
        vtkSmartPointer<vtkDoubleArray> yCoords = vtkSmartPointer<vtkDoubleArray>::New();
        for (size_t n=0; n<nn[1]; ++n) {
            yCoords->InsertNextValue( x[2] + n*d[1] );
        }
        vtkSmartPointer<vtkDoubleArray> zCoords = vtkSmartPointer<vtkDoubleArray>::New();
        for (size_t n=0; n<nn[2]; ++n) {
            zCoords->InsertNextValue( x[4] + n*d[2] );
        }

        vtkSmartPointer<vtkRectilinearGrid> rgrid = vtkSmartPointer<vtkRectilinearGrid>::New();
        rgrid->SetDimensions( nn );
        rgrid->SetXCoordinates(xCoords);
        rgrid->SetYCoordinates(yCoords);
        rgrid->SetZCoordinates(zCoords);

        vtkSmartPointer<vtkDoubleArray> data = vtkSmartPointer<vtkDoubleArray>::New();

        sxz<T1> pt;
        if ( saveSlowness ) {
            data->SetName("Slowness");
            for ( size_t n=0; n<rgrid->GetNumberOfCells(); ++n ) {

                rgrid->GetCell(n)->GetBounds(x);
                pt.x = 0.5*(x[0]+x[1]);
                pt.z = 0.5*(x[4]+x[5]);

                for ( T2 nt=0; nt<triangles.size(); ++nt ) {
                    if ( insideTriangle(pt, nt) ) {
                        data->InsertNextValue( cells.getSlowness(n) );
                        break;
                    }
                }
            }
        } else {
            data->SetName("Velocity");
            for ( size_t n=0; n<rgrid->GetNumberOfCells(); ++n ) {

                rgrid->GetCell(n)->GetBounds(x);
                pt.x = 0.5*(x[0]+x[1]);
                pt.z = 0.5*(x[4]+x[5]);

                for ( T2 nt=0; nt<triangles.size(); ++nt ) {
                    if ( insideTriangle(pt, nt) ) {
                        data->InsertNextValue( 1./cells.getSlowness(n) );
                        break;
                    }
                }
            }
        }

        rgrid->GetCellData()->SetScalars( data );

        vtkSmartPointer<vtkXMLRectilinearGridWriter> writer =
        vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();

        writer->SetFileName( fname.c_str() );
        writer->SetInputData( rgrid );
        writer->SetDataModeToBinary();
        writer->Update();

    }
#endif


    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Duc<T1,T2,S,NODE,CELL>::processObtuse() {

        const double pi2 = pi / 2.;

        for ( T2 ntri=0; ntri<triangles.size(); ++ntri ) {

            for ( T2 n=0; n<3; ++n ) {
                if ( triangles[ntri].a[n] > pi2 ) {

                    // look for opposite triangle

                    T2 i0 = triangles[ntri].i[n];
                    T2 i1 = triangles[ntri].i[(n+1)%3];
                    T2 i2 = triangles[ntri].i[(n+2)%3];

                    T2 oppositeTriangle = 0;
                    bool found = false;
                    for ( size_t n1=0; n1<nodes[i1].getOwners().size(); ++n1) {
                        for ( size_t n2=0; n2<nodes[i2].getOwners().size(); ++n2) {
                            if ( nodes[i2].getOwners()[n2] == nodes[i1].getOwners()[n1]) {
                                oppositeTriangle = nodes[i2].getOwners()[n2];
                                found = true;
                                break;
                            }

                        }
                        if ( found ) break;
                    }

                    if ( !found ) continue; // no opposite triangle, must be on edge of domain.  No correction applied.


                    // find opposite node
                    T2 i3 = triangles[oppositeTriangle].i[0];
                    if ( i3 == i1 || i3 == i2 )
                        i3 = triangles[oppositeTriangle].i[1];
                    else if ( i3 == i1 || i3 == i2 )
                        i3 = triangles[oppositeTriangle].i[2];

                    virtualNode<T1,NODE> vn;

                    // keep i1 & try replacing i2 with i3
                    vn.node1 = &(nodes[i1]);
                    vn.node2 = &(nodes[i3]);

                    // distance between node 1 & 3 (opposite of node 0)
                    T1 a = nodes[i1].getDistance( nodes[i3] );

                    // distance between node 0 & 3 (opposite of node 1)
                    T1 b = nodes[i0].getDistance( nodes[i3] );

                    // distance between node 0 & 1 (opposite of node 3)
                    T1 c = nodes[i0].getDistance( nodes[i1] );

                    // angle at node 0
                    T1 a0 = acos((b*b + c*c - a*a)/(2.*b*c));



                    if ( a0 > pi2 ) { // still obtuse -> replace i1 instead of i2 with i3

                        vn.node1 = &(nodes[i3]);
                        vn.node2 = &(nodes[i2]);

                        // distance between node 2 & 3 (opposite of node 0)
                        a = nodes[i2].getDistance( nodes[i3]);

                        // distance between node 0 & 2 (opposite of node 1)
                        b = nodes[i0].getDistance( nodes[i2]);

                        // distance between node 0 & 3 (opposite of node 2)
                        c = nodes[i0].getDistance( nodes[i3]);

                        a0 = acos((b*b + c*c - a*a)/(2.*b*c));


                    }

                    vn.a[0] = a0;
                    vn.a[1] = acos((c*c + a*a - b*b)/(2.*a*c));
                    vn.a[2] = acos((a*a + b*b - c*c)/(2.*a*b));

                    vn.e[0] = a;
                    vn.e[1] = b;
                    vn.e[2] = c;

                    virtualNodes.insert( std::pair<T2, virtualNode<T1,NODE>>(ntri, vn) );

                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Duc<T1,T2,S,NODE,CELL>::localSolver(NODE *vertexC,
                                             const size_t threadNo) const {

        static const double pi2 = pi / 2.;
        T2 i0, i1, i2;
        NODE *vertexA, *vertexB;
        T1 a, b, c, alpha, beta;

        for ( size_t no=0; no<vertexC->getOwners().size(); ++no ) {

            T2 triangleNo = vertexC->getOwners()[no];

            for ( i0=0; i0<3; ++i0 ) {
                if ( vertexC->getGridIndex() == triangles[triangleNo].i[i0] ) break;
            }

            if ( triangles[triangleNo].a[i0] > pi/2 && !virtualNodes.empty() ) {

                virtualNode<T1,NODE> vn = virtualNodes.at(triangleNo);

                vertexA = vn.node1;
                vertexB = vn.node2;

                c = vn.e[0];
                a = vn.e[1];
                b = vn.e[2];

                alpha = vn.a[2];
                beta = vn.a[1];
            } else {

                i1 = (i0+1)%3;
                i2 = (i0+2)%3;

                vertexA = &(nodes[triangles[triangleNo].i[i1]]);
                vertexB = &(nodes[triangles[triangleNo].i[i2]]);

                c = triangles[triangleNo].l[i0];
                a = triangles[triangleNo].l[i1];
                b = triangles[triangleNo].l[i2];

                alpha = triangles[triangleNo].a[i2];
                beta = triangles[triangleNo].a[i1];
            }

            if ( std::abs(vertexB->getTT(threadNo)-vertexA->getTT(threadNo)) <= c*cells.getSlowness(triangleNo)) {

                T1 theta = asin( std::abs(vertexB->getTT(threadNo)-vertexA->getTT(threadNo))/
                                (c*cells.getSlowness(triangleNo)) );

                if ( ((0.>alpha-pi2?0.:alpha-pi2)<=theta && theta<=(pi2-beta) ) ||
                    ((alpha-pi2)<=theta && theta<=(0.<pi2-beta?0.:pi2-beta)) ) {
                    T1 h = a*sin(alpha-theta);
                    T1 H = b*sin(beta+theta);

                    T1 t = 0.5*(h*cells.getSlowness(triangleNo)+vertexB->getTT(threadNo)) +
                    0.5*(H*cells.getSlowness(triangleNo)+vertexA->getTT(threadNo));

                    if ( t<vertexC->getTT(threadNo) )
                        vertexC->setTT(t, threadNo);
                } else {
                    T1 t = vertexA->getTT(threadNo) + b*cells.getSlowness(triangleNo);
                    t = t<vertexB->getTT(threadNo) + a*cells.getSlowness(triangleNo) ? t :
                    vertexB->getTT(threadNo) + a*cells.getSlowness(triangleNo);
                    if ( t<vertexC->getTT(threadNo) )
                        vertexC->setTT(t, threadNo);
                }
            } else {
                T1 t = vertexA->getTT(threadNo) + b*cells.getSlowness(triangleNo);
                t = t<vertexB->getTT(threadNo) + a*cells.getSlowness(triangleNo) ? t :
                vertexB->getTT(threadNo) + a*cells.getSlowness(triangleNo);
                if ( t<vertexC->getTT(threadNo) )
                    vertexC->setTT(t, threadNo);
            }
        }
    }


    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Duc<T1,T2,S,NODE,CELL>::initTxVars(const std::vector<sxz<T1>>& Tx,
                                                 std::vector<bool>& txOnNode,
                                                 std::vector<T2>& txNode,
                                                 std::vector<T2>& txCell,
                                                 std::vector<std::vector<T2>>& txCells) const {
        for ( size_t nt=0; nt<Tx.size(); ++nt ) {
            T2 nn = getNearestNode( Tx[nt] );
            if ( nodes[nn] == Tx[nt] ) {
                txOnNode[nt] = true;
                txNode[nt] = nn;
            }
        }
        for ( size_t nt=0; nt<Tx.size(); ++nt ) {
            if ( !txOnNode[nt] ) {
                T2 txc = getCellNo( Tx[nt] );
                txCell[nt] = txc;
                // find cells surrounding txCell
                // this is because gradient is inaccurate when we're getting
                // close to Tx, so we want to check if we are reaching one of
                // the surrounding cells
                std::set<T2> indices;
                indices.insert(triangles[txc].i[0]);
                indices.insert(triangles[txc].i[1]);
                indices.insert(triangles[txc].i[2]);
                for ( T2 ntr=0; ntr<triangles.size(); ++ntr ) {
                    if ( ntr == txc ) {
                        continue;
                    }
                    if (indices.find(triangles[ntr].i[0])!=indices.end() ||
                        indices.find(triangles[ntr].i[1])!=indices.end() ||
                        indices.find(triangles[ntr].i[2])!=indices.end()) {
                        txCells[nt].push_back(ntr);
                    }
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Duc<T1,T2,S,NODE,CELL>::getRaypath(const std::vector<sxz<T1>>& Tx,
                                            const sxz<T1> &Rx,
                                            std::vector<sxz<T1>> &r_data,
                                            const size_t threadNo) const {

        T1 minDist = small;
        r_data.push_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                return;
            }
        }

        const txInfo_t& txi = getTxInfo(Tx, threadNo);
        const std::vector<bool>& txOnNode = txi.txOnNode;
        const std::vector<T2>& txNode = txi.txNode;
        const std::vector<T2>& txCell = txi.txCell;
        const std::vector<std::vector<T2>>& txCells = txi.txCells;

        T2 cellNo = 0;
        T2 nodeNo = 0;
        sxz<T1> curr_pt( Rx );

        bool onNode = false;
        bool reachedTx = false;
        bool onEdge = false;
        std::array<T2,2> edgeNodes;
        Grad2D_ls_so<T1,NODE> grad2d;

        {
            T2 nn = getNearestNode( curr_pt );
            if ( nodes[nn] == curr_pt ) {
                if ( nodes[nn].isPrimary() ) {
                    nodeNo = nn;
                    onNode = true;
                } else {
                    onEdge = true;
                    cellNo = getCellNo( curr_pt );
                    if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[0]], nodes[triangles[cellNo].i[1]]) ) {
                        edgeNodes[0] = triangles[cellNo].i[0];
                        edgeNodes[1] = triangles[cellNo].i[1];
                    } else if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[0]], nodes[triangles[cellNo].i[2]]) ) {
                        edgeNodes[0] = triangles[cellNo].i[0];
                        edgeNodes[1] = triangles[cellNo].i[2];
                    } else {
                        edgeNodes[0] = triangles[cellNo].i[1];
                        edgeNodes[1] = triangles[cellNo].i[2];
                    }
                }
            }
        }
        if ( !onNode ) {
            cellNo = getCellNo( curr_pt );
            for ( auto nt=0; nt<txCell.size(); ++nt ) {
                if ( cellNo == txCell[nt] ) {
                    r_data.push_back( Tx[nt] );
                    reachedTx = true;
                    break;
                }
            }
            // check if on edge
            if ( !onEdge ) {
                if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[0]], nodes[triangles[cellNo].i[1]]) ) {
                    onEdge = true;
                    edgeNodes[0] = triangles[cellNo].i[0];
                    edgeNodes[1] = triangles[cellNo].i[1];
                } else if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[0]], nodes[triangles[cellNo].i[2]]) ) {
                    onEdge = true;
                    edgeNodes[0] = triangles[cellNo].i[0];
                    edgeNodes[1] = triangles[cellNo].i[2];
                } else if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[1]], nodes[triangles[cellNo].i[2]]) ) {
                    onEdge = true;
                    edgeNodes[0] = triangles[cellNo].i[1];
                    edgeNodes[1] = triangles[cellNo].i[2];
                }
            }
        }

        if ( onNode ) {
            for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                //  check if cell is (one of) TxCell(s)
                for (size_t nt=0; nt<Tx.size(); ++nt) {
                    if ( *nc == txCell[nt] ) {
                        r_data.push_back( Tx[nt] );
                        reachedTx = true;
                        break;
                    }
                }
            }
        } else {
            for (size_t nt=0; nt<Tx.size(); ++nt) {
                if ( cellNo == txCell[nt] ) {
                    r_data.push_back( Tx[nt] );
                    reachedTx = true;
                    break;
                }
            }
        }

        sxz<T1> g;
        while ( reachedTx == false ) {

            if ( onNode ) {

                // find cell for which gradient intersect opposing segment
                bool foundIntersection = false;
                std::vector<sxz<T1>> grads;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                    T2 nb[2];
                    size_t n=0;
                    for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                        if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                            nb[n++] = *nn;
                        }
                    }
                    if ( nb[0]>nb[1] ) std::swap(nb[0], nb[1]);

                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    getNeighborNodes(*nc, nnodes);

                    g = grad2d.compute(nnodes, threadNo);

                    sxz<T1> v1 = { nodes[ nb[0] ].getX() - nodes[ nodeNo ].getX(),
                        nodes[ nb[0] ].getZ() - nodes[ nodeNo ].getZ() };
                    sxz<T1> v2 = { nodes[ nb[1] ].getX() - nodes[ nodeNo ].getX(),
                        nodes[ nb[1] ].getZ() - nodes[ nodeNo ].getZ() };

                    g.normalize();
                    v1.normalize();
                    v2.normalize();

                    T1 theta1 = acos( dot(v1, g) );
                    T1 theta2 = acos( dot(v1, v2) );

                    if ( theta1 > theta2 ) {
                        grads.push_back( g );
                        continue;
                    }

                    if ( boost::math::sign( cross(v1, g) ) != boost::math::sign( cross(v1, v2) ) ) {
                        grads.push_back( g );
                        continue;
                    }

                    foundIntersection = true;

                    bool break_flag = findIntersection(nb[0], nb[1], g, curr_pt);

                    r_data.push_back( curr_pt );
                    if ( break_flag ) break;

                    onEdge = true;
                    edgeNodes[0] = nb[0];
                    edgeNodes[1] = nb[1];

                    // find next cell
                    cellNo = findNextCell1(nb[0], nb[1], nodeNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge for Rx "
                        << Rx.x << ' ' << Rx.z << std::endl;
                        r_data.resize(1);
                        r_data[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }

                if ( foundIntersection == false ) {

                    g = { 0., 0. };
                    // check if we are on a node close to Tx
                    bool closeToTx = false;
                    for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                        if ( txOnNode[nt] ) {
                            for ( auto nc=nodes[txNode[nt]].getOwners().begin(); nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                                if (find(this->neighbors[*nc].begin(), this->neighbors[*nc].end(), nodeNo) != this->neighbors[*nc].end()) {
                                    closeToTx = true;
                                    g = Tx[nt] - curr_pt;
                                    break;
                                }
                            }
                        } else {
                            // check if surrounding triangles include nodeNo
                            for ( size_t no=0; no<3; ++no ) {
                                T2 node = triangles[txCell[nt]].i[no];
                                for ( auto nc=nodes[node].getOwners().begin(); nc!=nodes[node].getOwners().end(); ++nc ) {
                                    if (find(this->neighbors[*nc].begin(), this->neighbors[*nc].end(), nodeNo) != this->neighbors[*nc].end()) {
                                        closeToTx = true;
                                        g = Tx[nt] - curr_pt;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if ( !closeToTx ) {
                        // compute average gradient
                        for ( size_t n=0; n<grads.size(); ++n ) {
                            g.x += grads[n].x;
                            g.z += grads[n].z;
                        }
                        g.x /= grads.size();
                        g.z /= grads.size();
                    }

                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                        T2 nb[2];
                        size_t n=0;
                        for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                            if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                                nb[n++] = *nn;
                            }
                        }
                        if ( nb[0]>nb[1] ) std::swap(nb[0], nb[1]);

                        sxz<T1> v1 = { nodes[ nb[0] ].getX() - nodes[ nodeNo ].getX(),
                            nodes[ nb[0] ].getZ() - nodes[ nodeNo ].getZ() };
                        sxz<T1> v2 = { nodes[ nb[1] ].getX() - nodes[ nodeNo ].getX(),
                            nodes[ nb[1] ].getZ() - nodes[ nodeNo ].getZ() };

                        g.normalize();
                        v1.normalize();
                        v2.normalize();

                        T1 theta1 = acos( dot(v1, g) );
                        T1 theta2 = acos( dot(v1, v2) );

                        if ( theta1 > theta2 ) {
                            continue;
                        }

                        if ( boost::math::sign( cross(v1, g) ) != boost::math::sign( cross(v1, v2) ) ) {
                            continue;
                        }

                        foundIntersection = true;

                        bool break_flag = findIntersection(nb[0], nb[1], g, curr_pt);

                        r_data.push_back( curr_pt );
                        if ( break_flag ) break;

                        onEdge = true;
                        edgeNodes[0] = nb[0];
                        edgeNodes[1] = nb[1];

                        // find next cell
                        cellNo = findNextCell1(nb[0], nb[1], nodeNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            std::cout << "\n\nWarning: finding raypath failed to converge for Rx "
                            << Rx.x << ' ' << Rx.z << std::endl;
                            r_data.resize(1);
                            r_data[0] = Rx;
                            reachedTx = true;
                        }
                        break;
                    }
                }
                if ( foundIntersection == false ) {

                    // we might be on a node on the outer limit of the mesh, with
                    // a gradient pointing slightly outside the mesh

                    // find node closest to gradient vector
                    sxz<T1> tentativeNode;
                    T1 distance = std::numeric_limits<T1>::max();
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                        for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                            if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                                // compute distance
                                sxz<T1> tmp_node = sxz<T1>(nodes[*nn]);
                                T1 tmp = distPointToLine(curr_pt, curr_pt+g, tmp_node );
                                if ( tmp < distance ) {
                                    // make sure we point in the same direction
                                    sxz<T1> tmp_vec = tmp_node - curr_pt;
                                    tmp_vec.normalize();
                                    if ( acos( dot(tmp_vec, g) ) < 0.5235 ) {
                                        // within 30°
                                        distance = tmp;
                                        tentativeNode = nodes[*nn];
                                    }
                                }
                            }
                        }
                    }

                    // check if distance is "small", i.e. less than 1/3 of edge length
                    if ( distance < 0.33 * curr_pt.getDistance(tentativeNode) ) {
                        curr_pt = tentativeNode;
                        r_data.push_back( curr_pt );
                        onNode = true;
                    } else {
                        std::cout << "\n\nWarning: finding raypath failed to converge for Rx "
                        << Rx.x << ' ' << Rx.z << std::endl;
                        r_data.resize(1);
                        r_data[0] = Rx;
                        reachedTx = true;
                    }
                }

            } else {
                // on edge

                std::byte nnodes_buf[8192];
                std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                std::pmr::set<NODE*> nnodes{&nnodes_pool};
                getNeighborNodes(cellNo, nnodes);

                g = grad2d.compute(nnodes, threadNo);

                for (size_t n=0; n<txCells.size(); ++n) {
                    for (auto txn=txCells[n].begin(); txn!=txCells[n].end(); ++txn) {
                        if (cellNo == *txn) {
                            g = Tx[n] - curr_pt;
                            break;
                        }
                    }
                }

                g.normalize();

                // we have 3 segments that we might intersect
                T2 ind[3][2] = { {this->neighbors[cellNo][0], this->neighbors[cellNo][1]},
                    {this->neighbors[cellNo][0], this->neighbors[cellNo][2]},
                    {this->neighbors[cellNo][1], this->neighbors[cellNo][2]} };

                for ( size_t ns=0; ns<3; ++ns ) {
                    if ( ind[ns][0]>ind[ns][1] )
                        std::swap( ind[ns][0], ind[ns][1] );
                }



                sxz<T1> pt_i;
                T1 m1, b1, m2, b2;
                bool foundIntersection = false;
                for ( size_t ns=0; ns<3; ++ns ) {


                    // equation of the edge segment
                    T1 den = nodes[ ind[ns][1] ].getX() - nodes[ ind[ns][0] ].getX();

                    if ( den == 0.0 ) {
                        m1 = INFINITY;
                        b1 = nodes[ ind[ns][1] ].getX();
                    } else {
                        m1 = ( nodes[ ind[ns][1] ].getZ() - nodes[ ind[ns][0] ].getZ() ) / den;
                        b1 = nodes[ ind[ns][1] ].getZ() - m1*nodes[ ind[ns][1] ].getX();
                    }

                    // equation of the vector starting at curr_pt & pointing along gradient
                    if ( g.x == 0.0 ) {
                        m2 = INFINITY;
                        b2 = curr_pt.x;
                    } else {
                        m2 = g.z/g.x;
                        b2 = curr_pt.z - m2*curr_pt.x;
                    }

                    if ( onEdge && ind[ns][0]==edgeNodes[0] && ind[ns][1]==edgeNodes[1] ) {

                        if ( std::abs(m1-m2)<small ) {
                            // curr_pt is on an edge and gradient is along the edge
                            // den is the direction of vector P0->P1 along x
                            if ( boost::math::sign(den) == boost::math::sign(g.x) ) {
                                curr_pt.x = nodes[ ind[ns][1] ].getX();
                                curr_pt.z = nodes[ ind[ns][1] ].getZ();
                                r_data.push_back( curr_pt );
                                foundIntersection = true;
                                break;
                            } else {
                                curr_pt.x = nodes[ ind[ns][0] ].getX();
                                curr_pt.z = nodes[ ind[ns][0] ].getZ();
                                r_data.push_back( curr_pt );
                                foundIntersection = true;
                                break;
                            }

                        }
                        continue;
                    }
                    // intersection of edge segment & gradient vector
                    if ( m1 == INFINITY ) {
                        pt_i.x = b1;
                        pt_i.z = m2*pt_i.x + b2;
                    } else if ( m2 == INFINITY ) {
                        pt_i.x = b2;
                        pt_i.z = m1*pt_i.x + b1;
                    } else {
                        pt_i.x = (b2-b1)/(m1-m2);
                        pt_i.z = m2*pt_i.x + b2;
                    }

                    sxz<T1> vec(pt_i.x-curr_pt.x, pt_i.z-curr_pt.z);
                    if ( dot(vec, g) <= 0.0 ) {
                        // we are not pointing in the same direction
                        continue;
                    }

                    if (((pt_i.x<=nodes[ ind[ns][1] ].getX() && pt_i.x>=nodes[ ind[ns][0] ].getX()) ||
                         (pt_i.x>=nodes[ ind[ns][1] ].getX() && pt_i.x<=nodes[ ind[ns][0] ].getX()) ||
                         (abs(pt_i.x-nodes[ ind[ns][1] ].getX()) < small2) ||
                         (abs(pt_i.x-nodes[ ind[ns][0] ].getX()) < small2) ) &&
                        ((pt_i.z<=nodes[ ind[ns][0] ].getZ() && pt_i.z>=nodes[ ind[ns][1] ].getZ()) ||
                         (pt_i.z>=nodes[ ind[ns][0] ].getZ() && pt_i.z<=nodes[ ind[ns][1] ].getZ()) ||
                         (abs(pt_i.z-nodes[ ind[ns][1] ].getZ()) < small2) ||
                         (abs(pt_i.z-nodes[ ind[ns][0] ].getZ()) < small2)))
                    {
                        foundIntersection = true;
                        r_data.push_back( pt_i );
                        curr_pt = pt_i;

                        onEdge = true;
                        edgeNodes[0] = ind[ns][0];
                        edgeNodes[1] = ind[ns][1];

                        // find next cell
                        cellNo = findNextCell2(ind[ns][0], ind[ns][1], cellNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            std::cout << "\n\nWarning: finding raypath failed to converge for Rx "
                            << Rx.x << ' ' << Rx.z << std::endl;
                            r_data.resize(1);
                            r_data[0] = Rx;
                            reachedTx = true;
                        }
                        break;
                    }

                }
                if ( foundIntersection == false ) {

                    // we must be on an edge with gradient pointing slightly outside triangle
                    sxz<T1> vec(nodes[ edgeNodes[1] ].getX() - nodes[ edgeNodes[0] ].getX(),
                                nodes[ edgeNodes[1] ].getZ() - nodes[ edgeNodes[0] ].getZ());

                    if ( dot(vec, g) > 0.0 ) {
                        curr_pt.x = nodes[ edgeNodes[1] ].getX();
                        curr_pt.z = nodes[ edgeNodes[1] ].getZ();
                        r_data.push_back( curr_pt );
                        foundIntersection = true;
                    } else {
                        curr_pt.x = nodes[ edgeNodes[0] ].getX();
                        curr_pt.z = nodes[ edgeNodes[0] ].getZ();
                        r_data.push_back( curr_pt );
                        foundIntersection = true;
                    }
                }
            }

            onNode = false;
            {
                T2 nn = getNearestNode( curr_pt );
                if ( nodes[nn] == curr_pt && nodes[nn].isPrimary() ) {
                    nodeNo = nn;
                    onNode = true;
                    onEdge = false;
                }
            }

            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( curr_pt.getDistance( Tx[nt] ) < minDist ) {
                        reachedTx = true;
                        break;
                    }
                }
            } else {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        for ( auto nc=nodes[txNode[nt]].getOwners().begin(); nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                            if ( cellNo == *nc ) {
                                r_data.push_back( Tx[nt] );
                                reachedTx = true;
                                break;
                            }
                        }
                    } else {
                        if ( cellNo == txCell[nt] ) {
                            r_data.push_back( Tx[nt] );
                            reachedTx = true;
                        }
                    }
                    if ( reachedTx ) break;
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Duc<T1,T2,S,NODE,CELL>::getRaypath(const std::vector<sxz<T1>>& Tx,
                                            const std::vector<T1>& t0,
                                            const sxz<T1> &Rx,
                                            std::vector<sxz<T1>> &r_data,
                                            T1 &tt,
                                            const size_t threadNo) const {
        T1 minDist = small;
        r_data.push_back( Rx );
        tt = 0.0;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        const txInfo_t& txi = getTxInfo(Tx, threadNo);
        const std::vector<bool>& txOnNode = txi.txOnNode;
        const std::vector<T2>& txNode = txi.txNode;
        const std::vector<T2>& txCell = txi.txCell;
        const std::vector<std::vector<T2>>& txCells = txi.txCells;

        T2 cellNo = 0;
        T2 nodeNo = 0;
        sxz<T1> curr_pt( Rx );

        bool onNode = false;
        bool reachedTx = false;
        bool onEdge = false;
        std::array<T2,2> edgeNodes;
        Grad2D_ls_so<T1,NODE> grad2d;

        {
            T2 nn = getNearestNode( curr_pt );
            if ( nodes[nn] == curr_pt ) {
                if ( nodes[nn].isPrimary() ) {
                    nodeNo = nn;
                    onNode = true;
                } else {
                    onEdge = true;
                    cellNo = getCellNo( curr_pt );
                    if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[0]], nodes[triangles[cellNo].i[1]]) ) {
                        edgeNodes[0] = triangles[cellNo].i[0];
                        edgeNodes[1] = triangles[cellNo].i[1];
                    } else if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[0]], nodes[triangles[cellNo].i[2]]) ) {
                        edgeNodes[0] = triangles[cellNo].i[0];
                        edgeNodes[1] = triangles[cellNo].i[2];
                    } else {
                        edgeNodes[0] = triangles[cellNo].i[1];
                        edgeNodes[1] = triangles[cellNo].i[2];
                    }
                }
            }
        }
        if ( !onNode ) {
            cellNo = getCellNo( curr_pt );
            for ( size_t nt=0; nt<txCell.size(); ++nt ) {
                if ( cellNo == txCell[nt] ) {
                    tt += t0[nt] + cells.computeDt(Tx[nt], r_data.back(), cellNo); // slowness[cellNo] * r_data.back().getDistance( Tx[nt] );
                    r_data.push_back( Tx[nt] );
                    reachedTx = true;
                    break;
                }
            }
            // check if on edge
            if ( !onEdge ) {
                if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[0]], nodes[triangles[cellNo].i[1]]) ) {
                    onEdge = true;
                    edgeNodes[0] = triangles[cellNo].i[0];
                    edgeNodes[1] = triangles[cellNo].i[1];
                } else if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[0]], nodes[triangles[cellNo].i[2]]) ) {
                    onEdge = true;
                    edgeNodes[0] = triangles[cellNo].i[0];
                    edgeNodes[1] = triangles[cellNo].i[2];
                } else if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[1]], nodes[triangles[cellNo].i[2]]) ) {
                    onEdge = true;
                    edgeNodes[0] = triangles[cellNo].i[1];
                    edgeNodes[1] = triangles[cellNo].i[2];
                }
            }
        }

        if ( onNode ) {
            for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                //  check if cell is (one of) TxCell(s)
                for (size_t nt=0; nt<Tx.size(); ++nt) {
                    if ( *nc == txCell[nt] ) {
                        tt += t0[nt] + cells.computeDt(Tx[nt], r_data.back(), *nc); //slowness[*nc] * r_data.back().getDistance( Tx[nt] );
                        r_data.push_back( Tx[nt] );
                        reachedTx = true;
                        break;
                    }
                }
            }
        } else {
            for (size_t nt=0; nt<Tx.size(); ++nt) {
                if ( cellNo == txCell[nt] ) {
                    tt += t0[nt] + cells.computeDt(Tx[nt], r_data.back(), cellNo); //slowness[cellNo] * r_data.back().getDistance( Tx[nt] );
                    r_data.push_back( Tx[nt] );
                    reachedTx = true;
                    break;
                }
            }
        }

        sxz<T1> g;
        while ( reachedTx == false ) {

            if ( onNode ) {

                // find cell for which gradient intersect opposing segment
                bool foundIntersection = false;
                std::vector<sxz<T1>> grads;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                    T2 nb[2];
                    size_t n=0;
                    for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                        if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                            nb[n++] = *nn;
                        }
                    }
                    if ( nb[0]>nb[1] ) std::swap(nb[0], nb[1]);

                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    getNeighborNodes(*nc, nnodes);

                    g = grad2d.compute(nnodes, threadNo);

                    sxz<T1> v1 = { nodes[ nb[0] ].getX() - nodes[ nodeNo ].getX(),
                        nodes[ nb[0] ].getZ() - nodes[ nodeNo ].getZ() };
                    sxz<T1> v2 = { nodes[ nb[1] ].getX() - nodes[ nodeNo ].getX(),
                        nodes[ nb[1] ].getZ() - nodes[ nodeNo ].getZ() };

                    g.normalize();
                    v1.normalize();
                    v2.normalize();

                    T1 theta1 = acos( dot(v1, g) );
                    T1 theta2 = acos( dot(v1, v2) );

                    if ( theta1 > theta2 ) {
                        grads.push_back( g );
                        continue;
                    }

                    if ( boost::math::sign( cross(v1, g) ) != boost::math::sign( cross(v1, v2) ) ) {
                        grads.push_back( g );
                        continue;
                    }

                    foundIntersection = true;

                    //  check if cell is (one of) TxCell(s)
                    for (size_t nt=0; nt<Tx.size(); ++nt) {
                        if ( *nc == txCell[nt] ) {
                            tt += t0[nt] + cells.computeDt(Tx[nt], r_data.back(), *nc); //slowness[*nc] * r_data.back().getDistance( Tx[nt] );
                            r_data.push_back( Tx[nt] );
                            reachedTx = true;
                            break;
                        }
                    }
                    if ( reachedTx ) {
                        break;
                    }

                    bool break_flag = findIntersection(nb[0], nb[1], g, curr_pt);

                    tt += cells.computeDt(curr_pt, r_data.back(), *nc); //slowness[*nc] * r_data.back().getDistance( curr_pt );
                    r_data.push_back( curr_pt );
                    if ( break_flag ) break;

                    onEdge = true;
                    edgeNodes[0] = nb[0];
                    edgeNodes[1] = nb[1];

                    // find next cell
                    cellNo = findNextCell1(nb[0], nb[1], nodeNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge for Rx "
                        << Rx.x << ' ' << Rx.z << std::endl;
                        r_data.resize(1);
                        r_data[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }

                if ( foundIntersection == false ) {

                    g = { 0., 0. };
                    // check if we are on a node close to Tx
                    bool closeToTx = false;
                    for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                        if ( txOnNode[nt] ) {
                            for ( auto nc=nodes[txNode[nt]].getOwners().begin(); nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                                if (find(this->neighbors[*nc].begin(), this->neighbors[*nc].end(), nodeNo) != this->neighbors[*nc].end()) {
                                    closeToTx = true;
                                    g = Tx[nt] - curr_pt;
                                    break;
                                }
                            }
                        } else {
                            // check if surrounding triangles include nodeNo
                            for ( size_t no=0; no<3; ++no ) {
                                T2 node = triangles[txCell[nt]].i[no];
                                for ( auto nc=nodes[node].getOwners().begin(); nc!=nodes[node].getOwners().end(); ++nc ) {
                                    if (find(this->neighbors[*nc].begin(), this->neighbors[*nc].end(), nodeNo) != this->neighbors[*nc].end()) {
                                        closeToTx = true;
                                        g = Tx[nt] - curr_pt;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if ( !closeToTx ) {
                        // compute average gradient
                        for ( size_t n=0; n<grads.size(); ++n ) {
                            g.x += grads[n].x;
                            g.z += grads[n].z;
                        }
                        g.x /= grads.size();
                        g.z /= grads.size();
                    }

                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                        T2 nb[2];
                        size_t n=0;
                        for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                            if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                                nb[n++] = *nn;
                            }
                        }
                        if ( nb[0]>nb[1] ) std::swap(nb[0], nb[1]);

                        sxz<T1> v1 = { nodes[ nb[0] ].getX() - nodes[ nodeNo ].getX(),
                            nodes[ nb[0] ].getZ() - nodes[ nodeNo ].getZ() };
                        sxz<T1> v2 = { nodes[ nb[1] ].getX() - nodes[ nodeNo ].getX(),
                            nodes[ nb[1] ].getZ() - nodes[ nodeNo ].getZ() };

                        g.normalize();
                        v1.normalize();
                        v2.normalize();

                        T1 theta1 = acos( dot(v1, g) );
                        T1 theta2 = acos( dot(v1, v2) );

                        if ( theta1 > theta2 ) {
                            continue;
                        }

                        if ( boost::math::sign( cross(v1, g) ) != boost::math::sign( cross(v1, v2) ) ) {
                            continue;
                        }

                        foundIntersection = true;

                        bool break_flag = findIntersection(nb[0], nb[1], g, curr_pt);

                        tt += cells.computeDt(curr_pt, r_data.back(), *nc); //slowness[*nc] * r_data.back().getDistance( curr_pt );
                        r_data.push_back( curr_pt );
                        if ( break_flag ) break;

                        onEdge = true;
                        edgeNodes[0] = nb[0];
                        edgeNodes[1] = nb[1];

                        // find next cell
                        cellNo = findNextCell1(nb[0], nb[1], nodeNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            std::cout << "\n\nWarning: finding raypath failed to converge for Rx "
                            << Rx.x << ' ' << Rx.z << std::endl;
                            r_data.resize(1);
                            r_data[0] = Rx;
                            reachedTx = true;
                        }
                        break;
                    }
                }
                if ( foundIntersection == false ) {

                    // we might be on a node on the outer limit of the mesh, with
                    // a gradient pointing slightly outside the mesh

                    // find node closest to gradient vector
                    sxz<T1> tentativeNode;
                    T1 distance = std::numeric_limits<T1>::max();
                    T2 common_cell = 0;
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                        for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                            if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                                // compute distance
                                sxz<T1> tmp_node = sxz<T1>(nodes[*nn]);
                                T1 tmp = distPointToLine(curr_pt, curr_pt+g, tmp_node );
                                if ( tmp < distance ) {
                                    // make sure we point in the same direction
                                    sxz<T1> tmp_vec = tmp_node - curr_pt;
                                    tmp_vec.normalize();
                                    if ( acos( dot(tmp_vec, g) ) < 0.5235 ) {
                                        // within 30°
                                        distance = tmp;
                                        tentativeNode = nodes[*nn];
                                        common_cell = *nc;  // keep track of cell
                                    }
                                }
                            }
                        }
                    }

                    // check if distance is "small", i.e. less than 1/3 of edge length
                    if ( distance < 0.33 * curr_pt.getDistance(tentativeNode) ) {
                        curr_pt = tentativeNode;
                        tt += cells.computeDt(curr_pt, r_data.back(), common_cell); //slowness[common_cell] * r_data.back().getDistance( curr_pt );
                        r_data.push_back( curr_pt );
                        onNode = true;
                    } else {
                        std::cout << "\n\nWarning: finding raypath failed to converge for Rx "
                        << Rx.x << ' ' << Rx.z << std::endl;
                        r_data.resize(1);
                        r_data[0] = Rx;
                        reachedTx = true;
                    }
                }

            } else {
                // on edge

                std::byte nnodes_buf[8192];
                std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                std::pmr::set<NODE*> nnodes{&nnodes_pool};
                getNeighborNodes(cellNo, nnodes);

                g = grad2d.compute(nnodes, threadNo);

                for (size_t n=0; n<txCells.size(); ++n) {
                    for (auto txn=txCells[n].begin(); txn!=txCells[n].end(); ++txn) {
                        if (cellNo == *txn) {
                            g = Tx[n] - curr_pt;
                            break;
                        }
                    }
                }

                g.normalize();

                // we have 3 segments that we might intersect
                T2 ind[3][2] = { {this->neighbors[cellNo][0], this->neighbors[cellNo][1]},
                    {this->neighbors[cellNo][0], this->neighbors[cellNo][2]},
                    {this->neighbors[cellNo][1], this->neighbors[cellNo][2]} };

                for ( size_t ns=0; ns<3; ++ns ) {
                    if ( ind[ns][0]>ind[ns][1] )
                        std::swap( ind[ns][0], ind[ns][1] );
                }

                sxz<T1> pt_i;
                T1 m1, b1, m2, b2;
                bool foundIntersection = false;
                for ( size_t ns=0; ns<3; ++ns ) {

                    // equation of the edge segment
                    T1 den = nodes[ ind[ns][1] ].getX() - nodes[ ind[ns][0] ].getX();

                    if ( den == 0.0 ) {
                        m1 = INFINITY;
                        b1 = nodes[ ind[ns][1] ].getX();
                    } else {
                        m1 = ( nodes[ ind[ns][1] ].getZ() - nodes[ ind[ns][0] ].getZ() ) / den;
                        b1 = nodes[ ind[ns][1] ].getZ() - m1*nodes[ ind[ns][1] ].getX();
                    }

                    // equation of the vector starting at curr_pt & pointing along gradient
                    if ( g.x == 0.0 ) {
                        m2 = INFINITY;
                        b2 = curr_pt.x;
                    } else {
                        m2 = g.z/g.x;
                        b2 = curr_pt.z - m2*curr_pt.x;
                    }

                    if ( onEdge && ind[ns][0]==edgeNodes[0] && ind[ns][1]==edgeNodes[1] ) {

                        if ( std::abs(m1-m2)<small ) {
                            // curr_pt is on an edge and gradient is along the edge
                            // den is the direction of vector P0->P1 along x
                            if ( boost::math::sign(den) == boost::math::sign(g.x) ) {
                                curr_pt.x = nodes[ ind[ns][1] ].getX();
                                curr_pt.z = nodes[ ind[ns][1] ].getZ();
                                sxz<T1> mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
                                cellNo = getCellNo(mid_pt);
                                tt += cells.computeDt(curr_pt, r_data.back(), cellNo); //slowness[cellNo] * r_data.back().getDistance( curr_pt );
                                r_data.push_back( curr_pt );
                                foundIntersection = true;
                                break;
                            } else {
                                curr_pt.x = nodes[ ind[ns][0] ].getX();
                                curr_pt.z = nodes[ ind[ns][0] ].getZ();
                                sxz<T1> mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
                                cellNo = getCellNo(mid_pt);
                                tt += cells.computeDt(curr_pt, r_data.back(), cellNo); //slowness[cellNo] * r_data.back().getDistance( curr_pt );
                                r_data.push_back( curr_pt );
                                foundIntersection = true;
                                break;
                            }

                        }
                        continue;
                    }
                    // intersection of edge segment & gradient vector
                    if ( m1 == INFINITY ) {
                        pt_i.x = b1;
                        pt_i.z = m2*pt_i.x + b2;
                    } else if ( m2 == INFINITY ) {
                        pt_i.x = b2;
                        pt_i.z = m1*pt_i.x + b1;
                    } else {
                        pt_i.x = (b2-b1)/(m1-m2);
                        pt_i.z = m2*pt_i.x + b2;
                    }

                    sxz<T1> vec(pt_i.x-curr_pt.x, pt_i.z-curr_pt.z);
                    if ( dot(vec, g) <= 0.0 ) {
                        // we are not pointing in the same direction
                        continue;
                    }

                    if (((pt_i.x<=nodes[ ind[ns][1] ].getX() && pt_i.x>=nodes[ ind[ns][0] ].getX()) ||
                         (pt_i.x>=nodes[ ind[ns][1] ].getX() && pt_i.x<=nodes[ ind[ns][0] ].getX()) ||
                         (abs(pt_i.x-nodes[ ind[ns][1] ].getX()) < small2) ||
                         (abs(pt_i.x-nodes[ ind[ns][0] ].getX()) < small2) ) &&
                        ((pt_i.z<=nodes[ ind[ns][0] ].getZ() && pt_i.z>=nodes[ ind[ns][1] ].getZ()) ||
                         (pt_i.z>=nodes[ ind[ns][0] ].getZ() && pt_i.z<=nodes[ ind[ns][1] ].getZ()) ||
                         (abs(pt_i.z-nodes[ ind[ns][1] ].getZ()) < small2) ||
                         (abs(pt_i.z-nodes[ ind[ns][0] ].getZ()) < small2)))
                    {
                        foundIntersection = true;
                        sxz<T1> mid_pt = static_cast<T1>(0.5)*(r_data.back() + pt_i);
                        cellNo = getCellNo(mid_pt);
                        tt += cells.computeDt(pt_i, r_data.back(), cellNo); //slowness[cellNo] * r_data.back().getDistance( pt_i );
                        r_data.push_back( pt_i );
                        curr_pt = pt_i;

                        onEdge = true;
                        edgeNodes[0] = ind[ns][0];
                        edgeNodes[1] = ind[ns][1];

                        // find next cell
                        cellNo = findNextCell2(ind[ns][0], ind[ns][1], cellNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            std::cout << "\n\nWarning: finding raypath failed to converge for Rx "
                            << Rx.x << ' ' << Rx.z << std::endl;
                            r_data.resize(1);
                            r_data[0] = Rx;
                            reachedTx = true;
                        }
                        break;
                    }

                }
                if ( foundIntersection == false ) {

                    // we must be on an edge with gradient pointing slightly outside triangle
                    sxz<T1> vec(nodes[ edgeNodes[1] ].getX() - nodes[ edgeNodes[0] ].getX(),
                                nodes[ edgeNodes[1] ].getZ() - nodes[ edgeNodes[0] ].getZ());

                    if ( dot(vec, g) > 0.0 ) {
                        curr_pt.x = nodes[ edgeNodes[1] ].getX();
                        curr_pt.z = nodes[ edgeNodes[1] ].getZ();
                        sxz<T1> mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
                        cellNo = getCellNo(mid_pt);
                        tt += cells.computeDt(curr_pt, r_data.back(), cellNo); //slowness[cellNo] * r_data.back().getDistance( curr_pt );
                        r_data.push_back( curr_pt );
                        foundIntersection = true;
                    } else {
                        curr_pt.x = nodes[ edgeNodes[0] ].getX();
                        curr_pt.z = nodes[ edgeNodes[0] ].getZ();
                        sxz<T1> mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
                        cellNo = getCellNo(mid_pt);
                        tt += cells.computeDt(curr_pt, r_data.back(), cellNo); //slowness[cellNo] * r_data.back().getDistance( curr_pt );
                        r_data.push_back( curr_pt );
                        foundIntersection = true;
                    }
                }

            }

            onNode = false;
            {
                T2 nn = getNearestNode( curr_pt );
                if ( nodes[nn] == curr_pt && nodes[nn].isPrimary() ) {
                    nodeNo = nn;
                    onNode = true;
                    onEdge = false;
                }
            }

            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( curr_pt.getDistance( Tx[nt] ) < minDist ) {
                        reachedTx = true;
                        break;
                    }
                }
            } else {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        for ( auto nc=nodes[txNode[nt]].getOwners().begin(); nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                            if ( cellNo == *nc ) {
                                tt += t0[nt] + cells.computeDt(Tx[nt], r_data.back(), cellNo); //slowness[cellNo] * r_data.back().getDistance( Tx[nt] );
                                r_data.push_back( Tx[nt] );
                                reachedTx = true;
                                break;
                            }
                        }
                    } else {
                        if ( cellNo == txCell[nt] ) {
                            tt += t0[nt] + cells.computeDt(Tx[nt], r_data.back(), cellNo); //slowness[cellNo] * r_data.back().getDistance( Tx[nt] );
                            r_data.push_back( Tx[nt] );
                            reachedTx = true;
                        }
                    }
                    if ( reachedTx ) break;
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Duc<T1,T2,S,NODE,CELL>::getRaypath(const std::vector<sxz<T1>>& Tx,
                                            const std::vector<T1>& t0,
                                            const sxz<T1> &Rx,
                                            std::vector<sxz<T1>> &r_data,
                                            std::vector<siv<T1>> &l_data,
                                            T1 &tt,
                                            const size_t threadNo) const {
        T1 minDist = small;
        r_data.push_back( Rx );
        tt = 0.0;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        const txInfo_t& txi = getTxInfo(Tx, threadNo);
        const std::vector<bool>& txOnNode = txi.txOnNode;
        const std::vector<T2>& txNode = txi.txNode;
        const std::vector<T2>& txCell = txi.txCell;
        const std::vector<std::vector<T2>>& txCells = txi.txCells;

        T2 cellNo = 0;
        T2 nodeNo = 0;
        sxz<T1> curr_pt( Rx );

        bool onNode = false;
        bool reachedTx = false;
        bool onEdge = false;
        std::array<T2,2> edgeNodes;
        Grad2D_ls_so<T1,NODE> grad2d;
        siv<T1> cell;

        {
            T2 nn = getNearestNode( curr_pt );
            if ( nodes[nn] == curr_pt ) {
                if ( nodes[nn].isPrimary() ) {
                    nodeNo = nn;
                    onNode = true;
                } else {
                    onEdge = true;
                    cellNo = getCellNo( curr_pt );
                    if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[0]], nodes[triangles[cellNo].i[1]]) ) {
                        edgeNodes[0] = triangles[cellNo].i[0];
                        edgeNodes[1] = triangles[cellNo].i[1];
                    } else if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[0]], nodes[triangles[cellNo].i[2]]) ) {
                        edgeNodes[0] = triangles[cellNo].i[0];
                        edgeNodes[1] = triangles[cellNo].i[2];
                    } else {
                        edgeNodes[0] = triangles[cellNo].i[1];
                        edgeNodes[1] = triangles[cellNo].i[2];
                    }
                }
            }
        }
        if ( !onNode ) {
            cellNo = getCellNo( curr_pt );
            for ( size_t nt=0; nt<txCell.size(); ++nt ) {
                if ( cellNo == txCell[nt] ) {
                    cell.i = cellNo;
                    cell.v = Tx[nt].getDistance(r_data.back());
                    l_data.push_back(cell);
                    tt += t0[nt] + cells.computeDt(Tx[nt], r_data.back(), cellNo); //slowness[cellNo] * cell.v;
                    r_data.push_back( Tx[nt] );
                    reachedTx = true;
                    break;
                }
            }
            // check if on edge
            if ( !onEdge ) {
                if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[0]], nodes[triangles[cellNo].i[1]]) ) {
                    onEdge = true;
                    edgeNodes[0] = triangles[cellNo].i[0];
                    edgeNodes[1] = triangles[cellNo].i[1];
                } else if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[0]], nodes[triangles[cellNo].i[2]]) ) {
                    onEdge = true;
                    edgeNodes[0] = triangles[cellNo].i[0];
                    edgeNodes[1] = triangles[cellNo].i[2];
                } else if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[1]], nodes[triangles[cellNo].i[2]]) ) {
                    onEdge = true;
                    edgeNodes[0] = triangles[cellNo].i[1];
                    edgeNodes[1] = triangles[cellNo].i[2];
                }
            }
        }

        if ( onNode ) {
            for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                //  check if cell is (one of) TxCell(s)
                for (size_t nt=0; nt<Tx.size(); ++nt) {
                    if ( *nc == txCell[nt] ) {
                        cell.i = *nc;
                        cell.v = Tx[nt].getDistance(r_data.back());
                        l_data.push_back(cell);
                        tt += t0[nt] + cells.computeDt(Tx[nt], r_data.back(), *nc); //slowness[*nc] * cell.v;
                        r_data.push_back( Tx[nt] );
                        reachedTx = true;
                        break;
                    }
                }
            }
        } else {
            for (size_t nt=0; nt<Tx.size(); ++nt) {
                if ( cellNo == txCell[nt] ) {
                    cell.i = cellNo;
                    cell.v = Tx[nt].getDistance(r_data.back());
                    l_data.push_back(cell);
                    tt += t0[nt] + cells.computeDt(Tx[nt], r_data.back(), cellNo); //slowness[cellNo] * cell.v;
                    r_data.push_back( Tx[nt] );
                    reachedTx = true;
                    break;
                }
            }
        }

        sxz<T1> g;
        while ( reachedTx == false ) {

            if ( onNode ) {

                // find cell for which gradient intersect opposing segment
                bool foundIntersection = false;
                std::vector<sxz<T1>> grads;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                    T2 nb[2];
                    size_t n=0;
                    for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                        if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                            nb[n++] = *nn;
                        }
                    }
                    if ( nb[0]>nb[1] ) std::swap(nb[0], nb[1]);

                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    getNeighborNodes(*nc, nnodes);

                    g = grad2d.compute(nnodes, threadNo);

                    sxz<T1> v1 = { nodes[ nb[0] ].getX() - nodes[ nodeNo ].getX(),
                        nodes[ nb[0] ].getZ() - nodes[ nodeNo ].getZ() };
                    sxz<T1> v2 = { nodes[ nb[1] ].getX() - nodes[ nodeNo ].getX(),
                        nodes[ nb[1] ].getZ() - nodes[ nodeNo ].getZ() };

                    g.normalize();
                    v1.normalize();
                    v2.normalize();

                    T1 theta1 = acos( dot(v1, g) );
                    T1 theta2 = acos( dot(v1, v2) );

                    if ( theta1 > theta2 ) {
                        grads.push_back( g );
                        continue;
                    }

                    if ( boost::math::sign( cross(v1, g) ) != boost::math::sign( cross(v1, v2) ) ) {
                        grads.push_back( g );
                        continue;
                    }

                    foundIntersection = true;

                    //  check if cell is (one of) TxCell(s)
                    for (size_t nt=0; nt<Tx.size(); ++nt) {
                        if ( *nc == txCell[nt] ) {
                            cell.i = *nc;
                            cell.v = Tx[nt].getDistance(r_data.back());
                            l_data.push_back(cell);
                            tt += t0[nt] + cells.computeDt(Tx[nt], r_data.back(), *nc); //slowness[*nc] * cell.v;
                            r_data.push_back( Tx[nt] );
                            reachedTx = true;
                            break;
                        }
                    }
                    if ( reachedTx ) {
                        break;
                    }

                    bool break_flag = findIntersection(nb[0], nb[1], g, curr_pt);

                    cell.i = *nc;
                    cell.v = r_data.back().getDistance( curr_pt );
                    l_data.push_back(cell);
                    tt += cells.computeDt(curr_pt, r_data.back(), *nc); //slowness[*nc] * cell.v;
                    r_data.push_back( curr_pt );
                    if ( break_flag ) break;

                    onEdge = true;
                    edgeNodes[0] = nb[0];
                    edgeNodes[1] = nb[1];

                    // find next cell
                    cellNo = findNextCell1(nb[0], nb[1], nodeNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge for Rx "
                        << Rx.x << ' ' << Rx.z << std::endl;
                        r_data.resize(1);
                        r_data[0] = Rx;
                        l_data.resize(0);
                        reachedTx = true;
                    }
                    break;
                }

                if ( foundIntersection == false ) {

                    g = { 0., 0. };
                    // check if we are on a node close to Tx
                    bool closeToTx = false;
                    for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                        if ( txOnNode[nt] ) {
                            for ( auto nc=nodes[txNode[nt]].getOwners().begin(); nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                                if (find(this->neighbors[*nc].begin(), this->neighbors[*nc].end(), nodeNo) != this->neighbors[*nc].end()) {
                                    closeToTx = true;
                                    g = Tx[nt] - curr_pt;
                                    break;
                                }
                            }
                        } else {
                            // check if surrounding triangles include nodeNo
                            for ( size_t no=0; no<3; ++no ) {
                                T2 node = triangles[txCell[nt]].i[no];
                                for ( auto nc=nodes[node].getOwners().begin(); nc!=nodes[node].getOwners().end(); ++nc ) {
                                    if (find(this->neighbors[*nc].begin(), this->neighbors[*nc].end(), nodeNo) != this->neighbors[*nc].end()) {
                                        closeToTx = true;
                                        g = Tx[nt] - curr_pt;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if ( !closeToTx ) {
                        // compute average gradient
                        for ( size_t n=0; n<grads.size(); ++n ) {
                            g.x += grads[n].x;
                            g.z += grads[n].z;
                        }
                        g.x /= grads.size();
                        g.z /= grads.size();
                    }

                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                        T2 nb[2];
                        size_t n=0;
                        for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                            if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                                nb[n++] = *nn;
                            }
                        }
                        if ( nb[0]>nb[1] ) std::swap(nb[0], nb[1]);

                        sxz<T1> v1 = { nodes[ nb[0] ].getX() - nodes[ nodeNo ].getX(),
                            nodes[ nb[0] ].getZ() - nodes[ nodeNo ].getZ() };
                        sxz<T1> v2 = { nodes[ nb[1] ].getX() - nodes[ nodeNo ].getX(),
                            nodes[ nb[1] ].getZ() - nodes[ nodeNo ].getZ() };

                        g.normalize();
                        v1.normalize();
                        v2.normalize();

                        T1 theta1 = acos( dot(v1, g) );
                        T1 theta2 = acos( dot(v1, v2) );

                        if ( theta1 > theta2 ) {
                            continue;
                        }

                        if ( boost::math::sign( cross(v1, g) ) != boost::math::sign( cross(v1, v2) ) ) {
                            continue;
                        }

                        foundIntersection = true;

                        bool break_flag = findIntersection(nb[0], nb[1], g, curr_pt);

                        cell.i = *nc;
                        cell.v = r_data.back().getDistance( curr_pt );
                        l_data.push_back(cell);
                        tt += cells.computeDt(curr_pt, r_data.back(), *nc); //slowness[*nc] * cell.v;
                        r_data.push_back( curr_pt );
                        if ( break_flag ) break;

                        onEdge = true;
                        edgeNodes[0] = nb[0];
                        edgeNodes[1] = nb[1];

                        // find next cell
                        cellNo = findNextCell1(nb[0], nb[1], nodeNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            std::cout << "\n\nWarning: finding raypath failed to converge for Rx "
                            << Rx.x << ' ' << Rx.z << std::endl;
                            r_data.resize(1);
                            r_data[0] = Rx;
                            l_data.resize(0);
                            reachedTx = true;
                        }
                        break;
                    }
                }
                if ( foundIntersection == false ) {

                    // we might be on a node on the outer limit of the mesh, with
                    // a gradient pointing slightly outside the mesh

                    // find node closest to gradient vector
                    sxz<T1> tentativeNode;
                    T1 distance = std::numeric_limits<T1>::max();
                    T2 common_cell = 0;
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                        for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                            if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                                // compute distance
                                sxz<T1> tmp_node = sxz<T1>(nodes[*nn]);
                                T1 tmp = distPointToLine(curr_pt, curr_pt+g, tmp_node );
                                if ( tmp < distance ) {
                                    // make sure we point in the same direction
                                    sxz<T1> tmp_vec = tmp_node - curr_pt;
                                    tmp_vec.normalize();
                                    if ( acos( dot(tmp_vec, g) ) < 0.5235 ) {
                                        // within 30°
                                        distance = tmp;
                                        tentativeNode = nodes[*nn];
                                        common_cell = *nc;  // keep track of cell
                                    }
                                }
                            }
                        }
                    }

                    // check if distance is "small", i.e. less than 1/3 of edge length
                    if ( distance < 0.33 * curr_pt.getDistance(tentativeNode) ) {
                        curr_pt = tentativeNode;
                        cell.i = common_cell;
                        cell.v = r_data.back().getDistance( curr_pt );
                        l_data.push_back(cell);
                        tt += cells.computeDt(curr_pt, r_data.back(), common_cell); //slowness[common_cell] * cell.v;
                        r_data.push_back( curr_pt );
                        onNode = true;
                    } else {
                        std::cout << "\n\nWarning: finding raypath failed to converge for Rx "
                        << Rx.x << ' ' << Rx.z << std::endl;
                        r_data.resize(1);
                        r_data[0] = Rx;
                        l_data.resize(0);
                        reachedTx = true;
                    }
                }

            } else {
                // on edge

                std::byte nnodes_buf[8192];
                std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                std::pmr::set<NODE*> nnodes{&nnodes_pool};
                getNeighborNodes(cellNo, nnodes);

                g = grad2d.compute(nnodes, threadNo);

                for (size_t n=0; n<txCells.size(); ++n) {
                    for (auto txn=txCells[n].begin(); txn!=txCells[n].end(); ++txn) {
                        if (cellNo == *txn) {
                            g = Tx[n] - curr_pt;
                            break;
                        }
                    }
                }

                g.normalize();

                // we have 3 segments that we might intersect
                T2 ind[3][2] = { {this->neighbors[cellNo][0], this->neighbors[cellNo][1]},
                    {this->neighbors[cellNo][0], this->neighbors[cellNo][2]},
                    {this->neighbors[cellNo][1], this->neighbors[cellNo][2]} };

                for ( size_t ns=0; ns<3; ++ns ) {
                    if ( ind[ns][0]>ind[ns][1] )
                        std::swap( ind[ns][0], ind[ns][1] );
                }

                sxz<T1> pt_i;
                T1 m1, b1, m2, b2;
                bool foundIntersection = false;
                for ( size_t ns=0; ns<3; ++ns ) {

                    // equation of the edge segment
                    T1 den = nodes[ ind[ns][1] ].getX() - nodes[ ind[ns][0] ].getX();

                    if ( den == 0.0 ) {
                        m1 = INFINITY;
                        b1 = nodes[ ind[ns][1] ].getX();
                    } else {
                        m1 = ( nodes[ ind[ns][1] ].getZ() - nodes[ ind[ns][0] ].getZ() ) / den;
                        b1 = nodes[ ind[ns][1] ].getZ() - m1*nodes[ ind[ns][1] ].getX();
                    }

                    // equation of the vector starting at curr_pt & pointing along gradient
                    if ( g.x == 0.0 ) {
                        m2 = INFINITY;
                        b2 = curr_pt.x;
                    } else {
                        m2 = g.z/g.x;
                        b2 = curr_pt.z - m2*curr_pt.x;
                    }

                    if ( onEdge && ind[ns][0]==edgeNodes[0] && ind[ns][1]==edgeNodes[1] ) {

                        if ( std::abs(m1-m2)<small ) {
                            // curr_pt is on an edge and gradient is along the edge
                            // den is the direction of vector P0->P1 along x
                            if ( boost::math::sign(den) == boost::math::sign(g.x) ) {
                                curr_pt.x = nodes[ ind[ns][1] ].getX();
                                curr_pt.z = nodes[ ind[ns][1] ].getZ();
                                sxz<T1> mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
                                cell.i = getCellNo(mid_pt);
                                cell.v = r_data.back().getDistance( curr_pt );
                                l_data.push_back(cell);
                                tt += cells.computeDt(curr_pt, r_data.back(), cell.i); //slowness[cell.i] * cell.v;
                                r_data.push_back( curr_pt );
                                foundIntersection = true;
                                break;
                            } else {
                                curr_pt.x = nodes[ ind[ns][0] ].getX();
                                curr_pt.z = nodes[ ind[ns][0] ].getZ();
                                sxz<T1> mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
                                cell.i = getCellNo(mid_pt);
                                cell.v = r_data.back().getDistance( curr_pt );
                                l_data.push_back(cell);
                                tt += cells.computeDt(curr_pt, r_data.back(), cell.i); //slowness[cell.i] * cell.v;
                                r_data.push_back( curr_pt );
                                foundIntersection = true;
                                break;
                            }

                        }
                        continue;
                    }
                    // intersection of edge segment & gradient vector
                    if ( m1 == INFINITY ) {
                        pt_i.x = b1;
                        pt_i.z = m2*pt_i.x + b2;
                    } else if ( m2 == INFINITY ) {
                        pt_i.x = b2;
                        pt_i.z = m1*pt_i.x + b1;
                    } else {
                        pt_i.x = (b2-b1)/(m1-m2);
                        pt_i.z = m2*pt_i.x + b2;
                    }

                    sxz<T1> vec(pt_i.x-curr_pt.x, pt_i.z-curr_pt.z);
                    if ( dot(vec, g) <= 0.0 ) {
                        // we are not pointing in the same direction
                        continue;
                    }

                    if (((pt_i.x<=nodes[ ind[ns][1] ].getX() && pt_i.x>=nodes[ ind[ns][0] ].getX()) ||
                         (pt_i.x>=nodes[ ind[ns][1] ].getX() && pt_i.x<=nodes[ ind[ns][0] ].getX()) ||
                         (abs(pt_i.x-nodes[ ind[ns][1] ].getX()) < small2) ||
                         (abs(pt_i.x-nodes[ ind[ns][0] ].getX()) < small2) ) &&
                        ((pt_i.z<=nodes[ ind[ns][0] ].getZ() && pt_i.z>=nodes[ ind[ns][1] ].getZ()) ||
                         (pt_i.z>=nodes[ ind[ns][0] ].getZ() && pt_i.z<=nodes[ ind[ns][1] ].getZ()) ||
                         (abs(pt_i.z-nodes[ ind[ns][1] ].getZ()) < small2) ||
                         (abs(pt_i.z-nodes[ ind[ns][0] ].getZ()) < small2)))
                    {
                        foundIntersection = true;
                        sxz<T1> mid_pt = static_cast<T1>(0.5)*(r_data.back() + pt_i);
                        cell.i = getCellNo(mid_pt);
                        cell.v = r_data.back().getDistance( pt_i );
                        l_data.push_back(cell);
                        tt += cells.computeDt(pt_i, r_data.back(), cell.i); //slowness[cell.i] * cell.v;
                        r_data.push_back( pt_i );
                        curr_pt = pt_i;

                        onEdge = true;
                        edgeNodes[0] = ind[ns][0];
                        edgeNodes[1] = ind[ns][1];

                        // find next cell
                        cellNo = findNextCell2(ind[ns][0], ind[ns][1], cellNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            std::cout << "\n\nWarning: finding raypath failed to converge for Rx "
                            << Rx.x << ' ' << Rx.z << std::endl;
                            r_data.resize(1);
                            r_data[0] = Rx;
                            l_data.resize(0);
                            reachedTx = true;
                        }
                        break;
                    }

                }
                if ( foundIntersection == false ) {

                    // we must be on an edge with gradient pointing slightly outside triangle
                    sxz<T1> vec(nodes[ edgeNodes[1] ].getX() - nodes[ edgeNodes[0] ].getX(),
                                nodes[ edgeNodes[1] ].getZ() - nodes[ edgeNodes[0] ].getZ());

                    if ( dot(vec, g) > 0.0 ) {
                        curr_pt.x = nodes[ edgeNodes[1] ].getX();
                        curr_pt.z = nodes[ edgeNodes[1] ].getZ();
                        sxz<T1> mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
                        cell.i = getCellNo(mid_pt);
                        cell.v = r_data.back().getDistance( curr_pt );
                        l_data.push_back(cell);
                        tt += cells.computeDt(curr_pt, r_data.back(), cell.i); //slowness[cell.i] * cell.v;
                        r_data.push_back( curr_pt );
                        foundIntersection = true;
                    } else {
                        curr_pt.x = nodes[ edgeNodes[0] ].getX();
                        curr_pt.z = nodes[ edgeNodes[0] ].getZ();
                        sxz<T1> mid_pt = static_cast<T1>(0.5)*(r_data.back() + curr_pt);
                        cell.i = getCellNo(mid_pt);
                        cell.v = r_data.back().getDistance( curr_pt );
                        l_data.push_back(cell);
                        tt += cells.computeDt(curr_pt, r_data.back(), cell.i); //slowness[cell.i] * cell.v;
                        r_data.push_back( curr_pt );
                        foundIntersection = true;
                    }
                }
            }

            onNode = false;
            {
                T2 nn = getNearestNode( curr_pt );
                if ( nodes[nn] == curr_pt && nodes[nn].isPrimary() ) {
                    nodeNo = nn;
                    onNode = true;
                    onEdge = false;
                }
            }

            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( curr_pt.getDistance( Tx[nt] ) < minDist ) {
                        reachedTx = true;
                        break;
                    }
                }
            } else {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        for ( auto nc=nodes[txNode[nt]].getOwners().begin(); nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                            if ( cellNo == *nc ) {
                                cell.i = cellNo;
                                cell.v = Tx[nt].getDistance(r_data.back());
                                l_data.push_back(cell);
                                tt += t0[nt] + cells.computeDt(Tx[nt], r_data.back(), cellNo); //slowness[cellNo] * cell.v;
                                r_data.push_back( Tx[nt] );
                                reachedTx = true;
                                break;
                            }
                        }
                    } else {
                        if ( cellNo == txCell[nt] ) {
                            cell.i = cellNo;
                            cell.v = Tx[nt].getDistance(r_data.back());
                            l_data.push_back(cell);
                            tt += t0[nt] + cells.computeDt(Tx[nt], r_data.back(), cellNo); //slowness[cellNo] * cell.v;
                            r_data.push_back( Tx[nt] );
                            reachedTx = true;
                        }
                    }
                    if ( reachedTx ) break;
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Duc<T1,T2,S,NODE,CELL>::getRaypath(const std::vector<sxz<T1>>& Tx,
                                            const std::vector<T1>& t0,
                                            const sxz<T1> &Rx,
                                            std::vector<siv<T1>> &l_data,
                                            T1 &tt,
                                            const size_t threadNo) const {
        T1 minDist = small;
        tt = 0.0;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        const txInfo_t& txi = getTxInfo(Tx, threadNo);
        const std::vector<bool>& txOnNode = txi.txOnNode;
        const std::vector<T2>& txNode = txi.txNode;
        const std::vector<T2>& txCell = txi.txCell;
        const std::vector<std::vector<T2>>& txCells = txi.txCells;

        T2 cellNo = 0;
        T2 nodeNo = 0;
        sxz<T1> prev_pt( Rx );
        sxz<T1> curr_pt( Rx );

        bool onNode = false;
        bool reachedTx = false;
        bool onEdge = false;
        std::array<T2,2> edgeNodes;
        Grad2D_ls_so<T1,NODE> grad2d;
        siv<T1> cell;

        {
            T2 nn = getNearestNode( curr_pt );
            if ( nodes[nn] == curr_pt ) {
                if ( nodes[nn].isPrimary() ) {
                    nodeNo = nn;
                    onNode = true;
                } else {
                    onEdge = true;
                    cellNo = getCellNo( curr_pt );
                    if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[0]], nodes[triangles[cellNo].i[1]]) ) {
                        edgeNodes[0] = triangles[cellNo].i[0];
                        edgeNodes[1] = triangles[cellNo].i[1];
                    } else if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[0]], nodes[triangles[cellNo].i[2]]) ) {
                        edgeNodes[0] = triangles[cellNo].i[0];
                        edgeNodes[1] = triangles[cellNo].i[2];
                    } else {
                        edgeNodes[0] = triangles[cellNo].i[1];
                        edgeNodes[1] = triangles[cellNo].i[2];
                    }
                }
            }
        }
        if ( !onNode ) {
            cellNo = getCellNo( curr_pt );
            for ( size_t nt=0; nt<txCell.size(); ++nt ) {
                if ( cellNo == txCell[nt] ) {
                    cell.i = cellNo;
                    cell.v = Tx[nt].getDistance(prev_pt);
                    l_data.push_back(cell);
                    tt += t0[nt] + cells.computeDt(Tx[nt], prev_pt, cellNo); //slowness[cellNo] * cell.v;
                    reachedTx = true;
                    break;
                }
            }
            // check if on edge
            if ( !onEdge ) {
                if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[0]], nodes[triangles[cellNo].i[1]]) ) {
                    onEdge = true;
                    edgeNodes[0] = triangles[cellNo].i[0];
                    edgeNodes[1] = triangles[cellNo].i[1];
                } else if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[0]], nodes[triangles[cellNo].i[2]]) ) {
                    onEdge = true;
                    edgeNodes[0] = triangles[cellNo].i[0];
                    edgeNodes[1] = triangles[cellNo].i[2];
                } else if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[1]], nodes[triangles[cellNo].i[2]]) ) {
                    onEdge = true;
                    edgeNodes[0] = triangles[cellNo].i[1];
                    edgeNodes[1] = triangles[cellNo].i[2];
                }
            }
        }

        if ( onNode ) {
            for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                //  check if cell is (one of) TxCell(s)
                for (size_t nt=0; nt<Tx.size(); ++nt) {
                    if ( *nc == txCell[nt] ) {
                        cell.i = *nc;
                        cell.v = Tx[nt].getDistance(prev_pt);
                        l_data.push_back(cell);
                        tt += t0[nt] + cells.computeDt(Tx[nt], prev_pt, *nc); //slowness[*nc] * cell.v;
                        reachedTx = true;
                        break;
                    }
                }
            }
        } else {
            for (size_t nt=0; nt<Tx.size(); ++nt) {
                if ( cellNo == txCell[nt] ) {
                    cell.i = cellNo;
                    cell.v = Tx[nt].getDistance(prev_pt);
                    l_data.push_back(cell);
                    tt += t0[nt] + cells.computeDt(Tx[nt], prev_pt, cellNo); //slowness[cellNo] * cell.v;
                    reachedTx = true;
                    break;
                }
            }
        }

        sxz<T1> g;
        while ( reachedTx == false ) {

            if ( onNode ) {

                // find cell for which gradient intersect opposing segment
                bool foundIntersection = false;
                std::vector<sxz<T1>> grads;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                    T2 nb[2];
                    size_t n=0;
                    for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                        if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                            nb[n++] = *nn;
                        }
                    }
                    if ( nb[0]>nb[1] ) std::swap(nb[0], nb[1]);

                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    getNeighborNodes(*nc, nnodes);

                    g = grad2d.compute(nnodes, threadNo);

                    sxz<T1> v1 = { nodes[ nb[0] ].getX() - nodes[ nodeNo ].getX(),
                        nodes[ nb[0] ].getZ() - nodes[ nodeNo ].getZ() };
                    sxz<T1> v2 = { nodes[ nb[1] ].getX() - nodes[ nodeNo ].getX(),
                        nodes[ nb[1] ].getZ() - nodes[ nodeNo ].getZ() };

                    g.normalize();
                    v1.normalize();
                    v2.normalize();

                    T1 theta1 = acos( dot(v1, g) );
                    T1 theta2 = acos( dot(v1, v2) );

                    if ( theta1 > theta2 ) {
                        grads.push_back( g );
                        continue;
                    }

                    if ( boost::math::sign( cross(v1, g) ) != boost::math::sign( cross(v1, v2) ) ) {
                        grads.push_back( g );
                        continue;
                    }

                    foundIntersection = true;

                    //  check if cell is (one of) TxCell(s)
                    for (size_t nt=0; nt<Tx.size(); ++nt) {
                        if ( *nc == txCell[nt] ) {
                            cell.i = *nc;
                            cell.v = Tx[nt].getDistance(prev_pt);
                            l_data.push_back(cell);
                            tt += t0[nt] + cells.computeDt(Tx[nt], prev_pt, *nc); //slowness[*nc] * cell.v;
                            reachedTx = true;
                            break;
                        }
                    }
                    if ( reachedTx ) {
                        break;
                    }

                    bool break_flag = findIntersection(nb[0], nb[1], g, curr_pt);

                    cell.i = *nc;
                    cell.v = prev_pt.getDistance( curr_pt );
                    l_data.push_back(cell);
                    tt += cells.computeDt(curr_pt, prev_pt, *nc); //slowness[*nc] * cell.v;
                    prev_pt = curr_pt;
                    if ( break_flag ) break;

                    onEdge = true;
                    edgeNodes[0] = nb[0];
                    edgeNodes[1] = nb[1];

                    // find next cell
                    cellNo = findNextCell1(nb[0], nb[1], nodeNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge for Rx "
                        << Rx.x << ' ' << Rx.z << std::endl;
                        l_data.resize(0);
                        reachedTx = true;
                    }
                    break;
                }

                if ( foundIntersection == false ) {

                    g = { 0., 0. };
                    // check if we are on a node close to Tx
                    bool closeToTx = false;
                    for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                        if ( txOnNode[nt] ) {
                            for ( auto nc=nodes[txNode[nt]].getOwners().begin(); nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                                if (find(this->neighbors[*nc].begin(), this->neighbors[*nc].end(), nodeNo) != this->neighbors[*nc].end()) {
                                    closeToTx = true;
                                    g = Tx[nt] - curr_pt;
                                    break;
                                }
                            }
                        } else {
                            // check if surrounding triangles include nodeNo
                            for ( size_t no=0; no<3; ++no ) {
                                T2 node = triangles[txCell[nt]].i[no];
                                for ( auto nc=nodes[node].getOwners().begin(); nc!=nodes[node].getOwners().end(); ++nc ) {
                                    if (find(this->neighbors[*nc].begin(), this->neighbors[*nc].end(), nodeNo) != this->neighbors[*nc].end()) {
                                        closeToTx = true;
                                        g = Tx[nt] - curr_pt;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if ( !closeToTx ) {
                        // compute average gradient
                        for ( size_t n=0; n<grads.size(); ++n ) {
                            g.x += grads[n].x;
                            g.z += grads[n].z;
                        }
                        g.x /= grads.size();
                        g.z /= grads.size();
                    }

                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                        T2 nb[2];
                        size_t n=0;
                        for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                            if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                                nb[n++] = *nn;
                            }
                        }
                        if ( nb[0]>nb[1] ) std::swap(nb[0], nb[1]);

                        sxz<T1> v1 = { nodes[ nb[0] ].getX() - nodes[ nodeNo ].getX(),
                            nodes[ nb[0] ].getZ() - nodes[ nodeNo ].getZ() };
                        sxz<T1> v2 = { nodes[ nb[1] ].getX() - nodes[ nodeNo ].getX(),
                            nodes[ nb[1] ].getZ() - nodes[ nodeNo ].getZ() };

                        g.normalize();
                        v1.normalize();
                        v2.normalize();

                        T1 theta1 = acos( dot(v1, g) );
                        T1 theta2 = acos( dot(v1, v2) );

                        if ( theta1 > theta2 ) {
                            continue;
                        }

                        if ( boost::math::sign( cross(v1, g) ) != boost::math::sign( cross(v1, v2) ) ) {
                            continue;
                        }

                        foundIntersection = true;

                        bool break_flag = findIntersection(nb[0], nb[1], g, curr_pt);

                        cell.i = *nc;
                        cell.v = prev_pt.getDistance( curr_pt );
                        l_data.push_back(cell);
                        tt += cells.computeDt(prev_pt, curr_pt, *nc); //slowness[*nc] * cell.v;
                        prev_pt = curr_pt;
                        if ( break_flag ) break;

                        onEdge = true;
                        edgeNodes[0] = nb[0];
                        edgeNodes[1] = nb[1];

                        // find next cell
                        cellNo = findNextCell1(nb[0], nb[1], nodeNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            std::cout << "\n\nWarning: finding raypath failed to converge for Rx "
                            << Rx.x << ' ' << Rx.z << std::endl;
                            l_data.resize(0);
                            reachedTx = true;
                        }
                        break;
                    }
                }
                if ( foundIntersection == false ) {

                    // we might be on a node on the outer limit of the mesh, with
                    // a gradient pointing slightly outside the mesh

                    // find node closest to gradient vector
                    sxz<T1> tentativeNode;
                    T1 distance = std::numeric_limits<T1>::max();
                    T2 common_cell = 0;
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                        for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                            if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                                // compute distance
                                sxz<T1> tmp_node = sxz<T1>(nodes[*nn]);
                                T1 tmp = distPointToLine(curr_pt, curr_pt+g, tmp_node );
                                if ( tmp < distance ) {
                                    // make sure we point in the same direction
                                    sxz<T1> tmp_vec = tmp_node - curr_pt;
                                    tmp_vec.normalize();
                                    if ( acos( dot(tmp_vec, g) ) < 0.5235 ) {
                                        // within 30°
                                        distance = tmp;
                                        tentativeNode = nodes[*nn];
                                        common_cell = *nc;  // keep track of cell
                                    }
                                }
                            }
                        }
                    }

                    // check if distance is "small", i.e. less than 1/3 of edge length
                    if ( distance < 0.33 * curr_pt.getDistance(tentativeNode) ) {
                        curr_pt = tentativeNode;
                        cell.i = common_cell;
                        cell.v = prev_pt.getDistance( curr_pt );
                        l_data.push_back(cell);
                        tt += cells.computeDt(prev_pt, curr_pt, common_cell); //slowness[common_cell] * cell.v;
                        prev_pt = curr_pt;
                        onNode = true;
                    } else {
                        std::cout << "\n\nWarning: finding raypath failed to converge for Rx "
                        << Rx.x << ' ' << Rx.z << std::endl;
                        l_data.resize(0);
                        reachedTx = true;
                    }
                }

            } else {
                // on edge

                std::byte nnodes_buf[8192];
                std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                std::pmr::set<NODE*> nnodes{&nnodes_pool};
                getNeighborNodes(cellNo, nnodes);

                g = grad2d.compute(nnodes, threadNo);

                for (size_t n=0; n<txCells.size(); ++n) {
                    for (auto txn=txCells[n].begin(); txn!=txCells[n].end(); ++txn) {
                        if (cellNo == *txn) {
                            g = Tx[n] - curr_pt;
                            break;
                        }
                    }
                }

                g.normalize();

                // we have 3 segments that we might intersect
                T2 ind[3][2] = { {this->neighbors[cellNo][0], this->neighbors[cellNo][1]},
                    {this->neighbors[cellNo][0], this->neighbors[cellNo][2]},
                    {this->neighbors[cellNo][1], this->neighbors[cellNo][2]} };

                for ( size_t ns=0; ns<3; ++ns ) {
                    if ( ind[ns][0]>ind[ns][1] )
                        std::swap( ind[ns][0], ind[ns][1] );
                }

                sxz<T1> pt_i;
                T1 m1, b1, m2, b2;
                bool foundIntersection = false;
                for ( size_t ns=0; ns<3; ++ns ) {

                    // equation of the edge segment
                    T1 den = nodes[ ind[ns][1] ].getX() - nodes[ ind[ns][0] ].getX();

                    if ( den == 0.0 ) {
                        m1 = INFINITY;
                        b1 = nodes[ ind[ns][1] ].getX();
                    } else {
                        m1 = ( nodes[ ind[ns][1] ].getZ() - nodes[ ind[ns][0] ].getZ() ) / den;
                        b1 = nodes[ ind[ns][1] ].getZ() - m1*nodes[ ind[ns][1] ].getX();
                    }

                    // equation of the vector starting at curr_pt & pointing along gradient
                    if ( g.x == 0.0 ) {
                        m2 = INFINITY;
                        b2 = curr_pt.x;
                    } else {
                        m2 = g.z/g.x;
                        b2 = curr_pt.z - m2*curr_pt.x;
                    }

                    if ( onEdge && ind[ns][0]==edgeNodes[0] && ind[ns][1]==edgeNodes[1] ) {

                        if ( std::abs(m1-m2)<small ) {
                            // curr_pt is on an edge and gradient is along the edge
                            // den is the direction of vector P0->P1 along x
                            if ( boost::math::sign(den) == boost::math::sign(g.x) ) {
                                curr_pt.x = nodes[ ind[ns][1] ].getX();
                                curr_pt.z = nodes[ ind[ns][1] ].getZ();
                                sxz<T1> mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
                                cell.i = getCellNo(mid_pt);
                                cell.v = prev_pt.getDistance( curr_pt );
                                l_data.push_back(cell);
                                tt += cells.computeDt(prev_pt, curr_pt, cell.i); //slowness[cell.i] * cell.v;
                                prev_pt = curr_pt;
                                foundIntersection = true;
                                break;
                            } else {
                                curr_pt.x = nodes[ ind[ns][0] ].getX();
                                curr_pt.z = nodes[ ind[ns][0] ].getZ();
                                sxz<T1> mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
                                cell.i = getCellNo(mid_pt);
                                cell.v = prev_pt.getDistance( curr_pt );
                                l_data.push_back(cell);
                                tt += cells.computeDt(prev_pt, curr_pt, cell.i); //slowness[cell.i] * cell.v;
                                prev_pt = curr_pt;
                                foundIntersection = true;
                                break;
                            }

                        }
                        continue;
                    }
                    // intersection of edge segment & gradient vector
                    if ( m1 == INFINITY ) {
                        pt_i.x = b1;
                        pt_i.z = m2*pt_i.x + b2;
                    } else if ( m2 == INFINITY ) {
                        pt_i.x = b2;
                        pt_i.z = m1*pt_i.x + b1;
                    } else {
                        pt_i.x = (b2-b1)/(m1-m2);
                        pt_i.z = m2*pt_i.x + b2;
                    }

                    sxz<T1> vec(pt_i.x-curr_pt.x, pt_i.z-curr_pt.z);
                    if ( dot(vec, g) <= 0.0 ) {
                        // we are not pointing in the same direction
                        continue;
                    }

                    if (((pt_i.x<=nodes[ ind[ns][1] ].getX() && pt_i.x>=nodes[ ind[ns][0] ].getX()) ||
                         (pt_i.x>=nodes[ ind[ns][1] ].getX() && pt_i.x<=nodes[ ind[ns][0] ].getX()) ||
                         (abs(pt_i.x-nodes[ ind[ns][1] ].getX()) < small2) ||
                         (abs(pt_i.x-nodes[ ind[ns][0] ].getX()) < small2) ) &&
                        ((pt_i.z<=nodes[ ind[ns][0] ].getZ() && pt_i.z>=nodes[ ind[ns][1] ].getZ()) ||
                         (pt_i.z>=nodes[ ind[ns][0] ].getZ() && pt_i.z<=nodes[ ind[ns][1] ].getZ()) ||
                         (abs(pt_i.z-nodes[ ind[ns][1] ].getZ()) < small2) ||
                         (abs(pt_i.z-nodes[ ind[ns][0] ].getZ()) < small2)))
                    {
                        foundIntersection = true;
                        sxz<T1> mid_pt = static_cast<T1>(0.5)*(prev_pt + pt_i);
                        cell.i = getCellNo(mid_pt);
                        cell.v = prev_pt.getDistance( pt_i );
                        l_data.push_back(cell);
                        tt += cells.computeDt(prev_pt, pt_i, cell.i); //slowness[cell.i] * cell.v;
                        prev_pt = pt_i;
                        curr_pt = pt_i;

                        onEdge = true;
                        edgeNodes[0] = ind[ns][0];
                        edgeNodes[1] = ind[ns][1];

                        // find next cell
                        cellNo = findNextCell2(ind[ns][0], ind[ns][1], cellNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            std::cout << "\n\nWarning: finding raypath failed to converge for Rx "
                            << Rx.x << ' ' << Rx.z << std::endl;
                            l_data.resize(0);
                            reachedTx = true;
                        }
                        break;
                    }

                }
                if ( foundIntersection == false ) {

                    // we must be on an edge with gradient pointing slightly outside triangle
                    sxz<T1> vec(nodes[ edgeNodes[1] ].getX() - nodes[ edgeNodes[0] ].getX(),
                                nodes[ edgeNodes[1] ].getZ() - nodes[ edgeNodes[0] ].getZ());

                    if ( dot(vec, g) > 0.0 ) {
                        curr_pt.x = nodes[ edgeNodes[1] ].getX();
                        curr_pt.z = nodes[ edgeNodes[1] ].getZ();
                        sxz<T1> mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
                        cell.i = getCellNo(mid_pt);
                        cell.v = prev_pt.getDistance( curr_pt );
                        l_data.push_back(cell);
                        tt += cells.computeDt(prev_pt, curr_pt, cell.i); //slowness[cell.i] * cell.v;
                        prev_pt = curr_pt;
                        foundIntersection = true;
                    } else {
                        curr_pt.x = nodes[ edgeNodes[0] ].getX();
                        curr_pt.z = nodes[ edgeNodes[0] ].getZ();
                        sxz<T1> mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
                        cell.i = getCellNo(mid_pt);
                        cell.v = prev_pt.getDistance( curr_pt );
                        l_data.push_back(cell);
                        tt += cells.computeDt(prev_pt, curr_pt, cell.i); //slowness[cell.i] * cell.v;
                        prev_pt = curr_pt;
                        foundIntersection = true;
                    }
                }
            }

            onNode = false;
            {
                T2 nn = getNearestNode( curr_pt );
                if ( nodes[nn] == curr_pt && nodes[nn].isPrimary() ) {
                    nodeNo = nn;
                    onNode = true;
                    onEdge = false;
                }
            }

            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( curr_pt.getDistance( Tx[nt] ) < minDist ) {
                        reachedTx = true;
                        break;
                    }
                }
            } else {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        for ( auto nc=nodes[txNode[nt]].getOwners().begin(); nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                            if ( cellNo == *nc ) {
                                cell.i = cellNo;
                                cell.v = Tx[nt].getDistance(prev_pt);
                                l_data.push_back(cell);
                                tt += t0[nt] + cells.computeDt(Tx[nt], prev_pt, cellNo); //slowness[cellNo] * cell.v;
                                reachedTx = true;
                                break;
                            }
                        }
                    } else {
                        if ( cellNo == txCell[nt] ) {
                            cell.i = cellNo;
                            cell.v = Tx[nt].getDistance(prev_pt);
                            l_data.push_back(cell);
                            tt += t0[nt] + cells.computeDt(Tx[nt], prev_pt, cellNo); //slowness[cellNo] * cell.v;
                            reachedTx = true;
                        }
                    }
                    if ( reachedTx ) break;
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    T1 Grid2Duc<T1,T2,S,NODE,CELL>::getTraveltimeFromRaypath(const std::vector<sxz<T1>>& Tx,
                                                        const std::vector<T1>& t0,
                                                        const sxz<T1> &Rx,
                                                        const size_t threadNo) const {

        T1 minDist = small;
        T1 tt = 0.0;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return tt;
            }
        }

        const txInfo_t& txi = getTxInfo(Tx, threadNo);
        const std::vector<bool>& txOnNode = txi.txOnNode;
        const std::vector<T2>& txNode = txi.txNode;
        const std::vector<T2>& txCell = txi.txCell;
        const std::vector<std::vector<T2>>& txCells = txi.txCells;

        T2 cellNo = 0;
        T2 nodeNo = 0;
        sxz<T1> prev_pt( Rx );
        sxz<T1> curr_pt( Rx );

        bool onNode = false;
        bool reachedTx = false;
        bool onEdge = false;
        std::array<T2,2> edgeNodes;
        Grad2D_ls_so<T1,NODE> grad2d;

        {
            T2 nn = getNearestNode( curr_pt );
            if ( nodes[nn] == curr_pt ) {
                if ( nodes[nn].isPrimary() ) {
                    nodeNo = nn;
                    onNode = true;
                } else {
                    onEdge = true;
                    cellNo = getCellNo( curr_pt );
                    if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[0]], nodes[triangles[cellNo].i[1]]) ) {
                        edgeNodes[0] = triangles[cellNo].i[0];
                        edgeNodes[1] = triangles[cellNo].i[1];
                    } else if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[0]], nodes[triangles[cellNo].i[2]]) ) {
                        edgeNodes[0] = triangles[cellNo].i[0];
                        edgeNodes[1] = triangles[cellNo].i[2];
                    } else {
                        edgeNodes[0] = triangles[cellNo].i[1];
                        edgeNodes[1] = triangles[cellNo].i[2];
                    }
                }
            }
        }
        if ( !onNode ) {
            cellNo = getCellNo( curr_pt );
            for ( size_t nt=0; nt<txCell.size(); ++nt ) {
                if ( cellNo == txCell[nt] ) {
                    tt += t0[nt] + cells.computeDt(Tx[nt], prev_pt, cellNo); //slowness[cellNo] * prev_pt.getDistance( Tx[nt] );
                    reachedTx = true;
                    break;
                }
            }
            // check if on edge
            if ( !onEdge ) {
                if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[0]], nodes[triangles[cellNo].i[1]]) ) {
                    onEdge = true;
                    edgeNodes[0] = triangles[cellNo].i[0];
                    edgeNodes[1] = triangles[cellNo].i[1];
                } else if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[0]], nodes[triangles[cellNo].i[2]]) ) {
                    onEdge = true;
                    edgeNodes[0] = triangles[cellNo].i[0];
                    edgeNodes[1] = triangles[cellNo].i[2];
                } else if ( areCollinear(curr_pt, nodes[triangles[cellNo].i[1]], nodes[triangles[cellNo].i[2]]) ) {
                    onEdge = true;
                    edgeNodes[0] = triangles[cellNo].i[1];
                    edgeNodes[1] = triangles[cellNo].i[2];
                }
            }
        }

        if ( onNode ) {
            for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                //  check if cell is (one of) TxCell(s)
                for (size_t nt=0; nt<Tx.size(); ++nt) {
                    if ( *nc == txCell[nt] ) {
                        tt += t0[nt] + cells.computeDt(Tx[nt], prev_pt, *nc); //slowness[*nc] * prev_pt.getDistance( Tx[nt] );
                        reachedTx = true;
                        break;
                    }
                }
            }
        } else {
            for (size_t nt=0; nt<Tx.size(); ++nt) {
                if ( cellNo == txCell[nt] ) {
                    tt += t0[nt] + cells.computeDt(Tx[nt], prev_pt, cellNo); //slowness[cellNo] * prev_pt.getDistance( Tx[nt] );
                    reachedTx = true;
                    break;
                }
            }
        }

        sxz<T1> g;
        while ( reachedTx == false ) {

            if ( onNode ) {

                // find cell for which gradient intersect opposing segment
                bool foundIntersection = false;
                std::vector<sxz<T1>> grads;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                    T2 nb[2];
                    size_t n=0;
                    for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                        if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                            nb[n++] = *nn;
                        }
                    }
                    if ( nb[0]>nb[1] ) std::swap(nb[0], nb[1]);

                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    getNeighborNodes(*nc, nnodes);

                    g = grad2d.compute(nnodes, threadNo);

                    sxz<T1> v1 = { nodes[ nb[0] ].getX() - nodes[ nodeNo ].getX(),
                        nodes[ nb[0] ].getZ() - nodes[ nodeNo ].getZ() };
                    sxz<T1> v2 = { nodes[ nb[1] ].getX() - nodes[ nodeNo ].getX(),
                        nodes[ nb[1] ].getZ() - nodes[ nodeNo ].getZ() };

                    g.normalize();
                    v1.normalize();
                    v2.normalize();

                    T1 theta1 = acos( dot(v1, g) );
                    T1 theta2 = acos( dot(v1, v2) );

                    if ( theta1 > theta2 ) {
                        grads.push_back( g );
                        continue;
                    }

                    if ( boost::math::sign( cross(v1, g) ) != boost::math::sign( cross(v1, v2) ) ) {
                        grads.push_back( g );
                        continue;
                    }

                    foundIntersection = true;

                    //  check if cell is (one of) TxCell(s)
                    for (size_t nt=0; nt<Tx.size(); ++nt) {
                        if ( *nc == txCell[nt] ) {
                            tt += t0[nt] + cells.computeDt(Tx[nt], prev_pt, *nc); //slowness[*nc] * prev_pt.getDistance( Tx[nt] );
                            reachedTx = true;
                            break;
                        }
                    }
                    if ( reachedTx ) {
                        break;
                    }

                    bool break_flag = findIntersection(nb[0], nb[1], g, curr_pt);

                    tt += cells.computeDt(prev_pt, curr_pt, *nc); //slowness[*nc] * prev_pt.getDistance( curr_pt );
                    prev_pt = curr_pt;
                    if ( break_flag ) break;

                    onEdge = true;
                    edgeNodes[0] = nb[0];
                    edgeNodes[1] = nb[1];

                    // find next cell
                    cellNo = findNextCell1(nb[0], nb[1], nodeNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge for Rx "
                        << Rx.x << ' ' << Rx.z << std::endl;
                        reachedTx = true;
                    }
                    break;
                }

                if ( foundIntersection == false ) {

                    g = { 0., 0. };
                    // check if we are on a node close to Tx
                    bool closeToTx = false;
                    for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                        if ( txOnNode[nt] ) {
                            for ( auto nc=nodes[txNode[nt]].getOwners().begin(); nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                                if (find(this->neighbors[*nc].begin(), this->neighbors[*nc].end(), nodeNo) != this->neighbors[*nc].end()) {
                                    closeToTx = true;
                                    g = Tx[nt] - curr_pt;
                                    break;
                                }
                            }
                        } else {
                            // check if surrounding triangles include nodeNo
                            for ( size_t no=0; no<3; ++no ) {
                                T2 node = triangles[txCell[nt]].i[no];
                                for ( auto nc=nodes[node].getOwners().begin(); nc!=nodes[node].getOwners().end(); ++nc ) {
                                    if (find(this->neighbors[*nc].begin(), this->neighbors[*nc].end(), nodeNo) != this->neighbors[*nc].end()) {
                                        closeToTx = true;
                                        g = Tx[nt] - curr_pt;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if ( !closeToTx ) {
                        // compute average gradient
                        for ( size_t n=0; n<grads.size(); ++n ) {
                            g.x += grads[n].x;
                            g.z += grads[n].z;
                        }
                        g.x /= grads.size();
                        g.z /= grads.size();
                    }

                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                        T2 nb[2];
                        size_t n=0;
                        for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                            if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                                nb[n++] = *nn;
                            }
                        }
                        if ( nb[0]>nb[1] ) std::swap(nb[0], nb[1]);

                        sxz<T1> v1 = { nodes[ nb[0] ].getX() - nodes[ nodeNo ].getX(),
                            nodes[ nb[0] ].getZ() - nodes[ nodeNo ].getZ() };
                        sxz<T1> v2 = { nodes[ nb[1] ].getX() - nodes[ nodeNo ].getX(),
                            nodes[ nb[1] ].getZ() - nodes[ nodeNo ].getZ() };

                        g.normalize();
                        v1.normalize();
                        v2.normalize();

                        T1 theta1 = acos( dot(v1, g) );
                        T1 theta2 = acos( dot(v1, v2) );

                        if ( theta1 > theta2 ) {
                            continue;
                        }

                        if ( boost::math::sign( cross(v1, g) ) != boost::math::sign( cross(v1, v2) ) ) {
                            continue;
                        }

                        foundIntersection = true;

                        bool break_flag = findIntersection(nb[0], nb[1], g, curr_pt);

                        tt += cells.computeDt(prev_pt, curr_pt, *nc); //slowness[*nc] * prev_pt.getDistance( curr_pt );
                        prev_pt = curr_pt;
                        if ( break_flag ) break;

                        onEdge = true;
                        edgeNodes[0] = nb[0];
                        edgeNodes[1] = nb[1];

                        // find next cell
                        cellNo = findNextCell1(nb[0], nb[1], nodeNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            std::cout << "\n\nWarning: finding raypath failed to converge for Rx "
                            << Rx.x << ' ' << Rx.z << std::endl;
                            reachedTx = true;
                        }
                        break;
                    }
                }
                if ( foundIntersection == false ) {
                    // we might be on a node on the outer limit of the mesh, with
                    // a gradient pointing slightly outside the mesh

                    // find node closest to gradient vector
                    sxz<T1> tentativeNode;
                    T1 distance = std::numeric_limits<T1>::max();
                    T2 common_cell = 0;
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                        for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                            if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                                // compute distance
                                sxz<T1> tmp_node = sxz<T1>(nodes[*nn]);
                                T1 tmp = distPointToLine(curr_pt, curr_pt+g, tmp_node );
                                if ( tmp < distance ) {
                                    // make sure we point in the same direction
                                    sxz<T1> tmp_vec = tmp_node - curr_pt;
                                    tmp_vec.normalize();
                                    if ( acos( dot(tmp_vec, g) ) < 0.5235 ) {
                                        // within 30°
                                        distance = tmp;
                                        tentativeNode = nodes[*nn];
                                        common_cell = *nc;  // keep track of cell
                                    }
                                }
                            }
                        }
                    }

                    // check if distance is "small", i.e. less than 1/3 of edge length
                    if ( distance < 0.33 * curr_pt.getDistance(tentativeNode) ) {
                        curr_pt = tentativeNode;
                        tt += cells.computeDt(prev_pt, curr_pt, common_cell); //slowness[common_cell] * prev_pt.getDistance( curr_pt );
                        prev_pt = curr_pt;
                        onNode = true;
                    } else {
                        std::cout << "\n\nWarning: finding raypath failed to converge for Rx "
                        << Rx.x << ' ' << Rx.z << std::endl;
                        reachedTx = true;
                    }
                }

            } else {
                // on edge

                std::byte nnodes_buf[8192];
                std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                std::pmr::set<NODE*> nnodes{&nnodes_pool};
                getNeighborNodes(cellNo, nnodes);

                g = grad2d.compute(nnodes, threadNo);

                for (size_t n=0; n<txCells.size(); ++n) {
                    for (auto txn=txCells[n].begin(); txn!=txCells[n].end(); ++txn) {
                        if (cellNo == *txn) {
                            g = Tx[n] - curr_pt;
                            break;
                        }
                    }
                }

                g.normalize();

                // we have 3 segments that we might intersect
                T2 ind[3][2] = { {this->neighbors[cellNo][0], this->neighbors[cellNo][1]},
                    {this->neighbors[cellNo][0], this->neighbors[cellNo][2]},
                    {this->neighbors[cellNo][1], this->neighbors[cellNo][2]} };

                for ( size_t ns=0; ns<3; ++ns ) {
                    if ( ind[ns][0]>ind[ns][1] )
                        std::swap( ind[ns][0], ind[ns][1] );
                }

                sxz<T1> pt_i;
                T1 m1, b1, m2, b2;
                bool foundIntersection = false;
                for ( size_t ns=0; ns<3; ++ns ) {

                    // equation of the edge segment
                    T1 den = nodes[ ind[ns][1] ].getX() - nodes[ ind[ns][0] ].getX();

                    if ( den == 0.0 ) {
                        m1 = INFINITY;
                        b1 = nodes[ ind[ns][1] ].getX();
                    } else {
                        m1 = ( nodes[ ind[ns][1] ].getZ() - nodes[ ind[ns][0] ].getZ() ) / den;
                        b1 = nodes[ ind[ns][1] ].getZ() - m1*nodes[ ind[ns][1] ].getX();
                    }

                    // equation of the vector starting at curr_pt & pointing along gradient
                    if ( g.x == 0.0 ) {
                        m2 = INFINITY;
                        b2 = curr_pt.x;
                    } else {
                        m2 = g.z/g.x;
                        b2 = curr_pt.z - m2*curr_pt.x;
                    }

                    if ( onEdge && ind[ns][0]==edgeNodes[0] && ind[ns][1]==edgeNodes[1] ) {

                        if ( std::abs(m1-m2)<small ) {
                            // curr_pt is on an edge and gradient is along the edge
                            // den is the direction of vector P0->P1 along x
                            if ( boost::math::sign(den) == boost::math::sign(g.x) ) {
                                curr_pt.x = nodes[ ind[ns][1] ].getX();
                                curr_pt.z = nodes[ ind[ns][1] ].getZ();
                                sxz<T1> mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
                                cellNo = getCellNo(mid_pt);
                                tt += cells.computeDt(prev_pt, curr_pt, cellNo); //slowness[cellNo] * prev_pt.getDistance( curr_pt );
                                prev_pt = curr_pt;
                                foundIntersection = true;
                                break;
                            } else {
                                curr_pt.x = nodes[ ind[ns][0] ].getX();
                                curr_pt.z = nodes[ ind[ns][0] ].getZ();
                                sxz<T1> mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
                                cellNo = getCellNo(mid_pt);
                                tt += cells.computeDt(prev_pt, curr_pt, cellNo); //slowness[cellNo] * prev_pt.getDistance( curr_pt );
                                prev_pt = curr_pt;
                                foundIntersection = true;
                                break;
                            }

                        }
                        continue;
                    }
                    // intersection of edge segment & gradient vector
                    if ( m1 == INFINITY ) {
                        pt_i.x = b1;
                        pt_i.z = m2*pt_i.x + b2;
                    } else if ( m2 == INFINITY ) {
                        pt_i.x = b2;
                        pt_i.z = m1*pt_i.x + b1;
                    } else {
                        pt_i.x = (b2-b1)/(m1-m2);
                        pt_i.z = m2*pt_i.x + b2;
                    }

                    sxz<T1> vec(pt_i.x-curr_pt.x, pt_i.z-curr_pt.z);
                    if ( dot(vec, g) <= 0.0 ) {
                        // we are not pointing in the same direction
                        continue;
                    }

                    if (((pt_i.x<=nodes[ ind[ns][1] ].getX() && pt_i.x>=nodes[ ind[ns][0] ].getX()) ||
                         (pt_i.x>=nodes[ ind[ns][1] ].getX() && pt_i.x<=nodes[ ind[ns][0] ].getX()) ||
                         (abs(pt_i.x-nodes[ ind[ns][1] ].getX()) < small2) ||
                         (abs(pt_i.x-nodes[ ind[ns][0] ].getX()) < small2) ) &&
                        ((pt_i.z<=nodes[ ind[ns][0] ].getZ() && pt_i.z>=nodes[ ind[ns][1] ].getZ()) ||
                         (pt_i.z>=nodes[ ind[ns][0] ].getZ() && pt_i.z<=nodes[ ind[ns][1] ].getZ()) ||
                         (abs(pt_i.z-nodes[ ind[ns][1] ].getZ()) < small2) ||
                         (abs(pt_i.z-nodes[ ind[ns][0] ].getZ()) < small2)))
                    {
                        foundIntersection = true;
                        sxz<T1> mid_pt = static_cast<T1>(0.5)*(prev_pt + pt_i);
                        cellNo = getCellNo(mid_pt);
                        tt += cells.computeDt(prev_pt, pt_i, cellNo); //slowness[cellNo] * prev_pt.getDistance( pt_i );
                        prev_pt = pt_i;
                        curr_pt = pt_i;

                        onEdge = true;
                        edgeNodes[0] = ind[ns][0];
                        edgeNodes[1] = ind[ns][1];

                        // find next cell
                        cellNo = findNextCell2(ind[ns][0], ind[ns][1], cellNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            std::cout << "\n\nWarning: finding raypath failed to converge for Rx "
                            << Rx.x << ' ' << Rx.z << std::endl;
                            reachedTx = true;
                        }
                        break;
                    }

                }
                if ( foundIntersection == false ) {

                    // we must be on an edge with gradient pointing slightly outside triangle
                    sxz<T1> vec(nodes[ edgeNodes[1] ].getX() - nodes[ edgeNodes[0] ].getX(),
                                nodes[ edgeNodes[1] ].getZ() - nodes[ edgeNodes[0] ].getZ());

                    if ( dot(vec, g) > 0.0 ) {
                        curr_pt.x = nodes[ edgeNodes[1] ].getX();
                        curr_pt.z = nodes[ edgeNodes[1] ].getZ();
                        sxz<T1> mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
                        cellNo = getCellNo(mid_pt);
                        tt += cells.computeDt(prev_pt, curr_pt, cellNo); //slowness[cellNo] * prev_pt.getDistance( curr_pt );
                        prev_pt = curr_pt;
                        foundIntersection = true;
                    } else {
                        curr_pt.x = nodes[ edgeNodes[0] ].getX();
                        curr_pt.z = nodes[ edgeNodes[0] ].getZ();
                        sxz<T1> mid_pt = static_cast<T1>(0.5)*(prev_pt + curr_pt);
                        cellNo = getCellNo(mid_pt);
                        tt += cells.computeDt(prev_pt, curr_pt, cellNo); //slowness[cellNo] * prev_pt.getDistance( curr_pt );
                        prev_pt = curr_pt;
                        foundIntersection = true;
                    }
                }

            }

            onNode = false;
            {
                T2 nn = getNearestNode( curr_pt );
                if ( nodes[nn] == curr_pt && nodes[nn].isPrimary() ) {
                    nodeNo = nn;
                    onNode = true;
                    onEdge = false;
                }
            }

            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( curr_pt.getDistance( Tx[nt] ) < minDist ) {
                        reachedTx = true;
                        break;
                    }
                }
            } else {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        for ( auto nc=nodes[txNode[nt]].getOwners().begin(); nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                            if ( cellNo == *nc ) {
                                tt += t0[nt] + cells.computeDt(Tx[nt], prev_pt, cellNo); //slowness[cellNo] * prev_pt.getDistance( Tx[nt] );
                                reachedTx = true;
                                break;
                            }
                        }
                    } else {
                        if ( cellNo == txCell[nt] ) {
                            tt += t0[nt] + cells.computeDt(Tx[nt], prev_pt, cellNo); //slowness[cellNo] * prev_pt.getDistance( Tx[nt] );
                            reachedTx = true;
                        }
                    }
                    if ( reachedTx ) break;
                }
            }
        }
        return tt;
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    bool Grid2Duc<T1,T2,S,NODE,CELL>::findIntersection(const T2 i0, const T2 i1,
                                                       const sxz<T1> &g,
                                                       sxz<T1> &curr_pt) const {

        // equation of the vector starting at curr_pt & pointing along gradient
        T1 m2, b2;
        if ( g.x == 0.0 ) {
            m2 = INFINITY;
            b2 = curr_pt.x;
        } else {
            m2 = g.z/g.x;
            b2 = curr_pt.z - m2*curr_pt.x;
        }

        // is gradient direction the same as one of the two edges

        // slope of 1st edge segment
        T1 den = nodes[ i0 ].getX() - curr_pt.x;

        T1 m1;
        if ( den == 0.0 ) m1 = INFINITY;
        else m1 = ( nodes[ i0 ].getZ() - curr_pt.z ) / den;

        if ( m1 == m2 ) {
            curr_pt.x = nodes[ i0 ].getX();
            curr_pt.z = nodes[ i0 ].getZ();
            return true;
        }

        // slope of 2nd edge segment
        den = nodes[ i1 ].getX() - curr_pt.x;
        if ( den == 0.0 ) m1 = INFINITY;
        else m1 = ( nodes[ i1 ].getZ() - curr_pt.z ) / den;

        if ( m1 == m2 ) {
            curr_pt.x = nodes[ i1 ].getX();
            curr_pt.z = nodes[ i1 ].getZ();
            return true;
        }

        // slope of opposing edge segment
        den = nodes[ i1 ].getX() - nodes[ i0 ].getX();
        T1 b1;
        if ( den == 0.0 ) {
            m1 = INFINITY;
            b1 = nodes[ i1 ].getX();
        } else {
            m1 = ( nodes[ i1 ].getZ() - nodes[ i0 ].getZ() ) / den;
            b1 = nodes[ i1 ].getZ() - m1*nodes[ i1 ].getX();
        }

        sxz<T1> pt_i;
        // intersection of edge segment & gradient vector
        if ( m1 == INFINITY ) {
            pt_i.x = b1;
            pt_i.z = m2*pt_i.x + b2;
        } else if ( m2 == INFINITY ) {
            pt_i.x = b2;
            pt_i.z = m1*pt_i.x + b1;
        } else {
            pt_i.x = (b2-b1)/(m1-m2);
            pt_i.z = m2*pt_i.x + b2;
        }

        curr_pt = pt_i;

        return false;
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    T2 Grid2Duc<T1,T2,S,NODE,CELL>::findNextCell1(const T2 i0, const T2 i1, const T2 nodeNo) const {
        std::vector<T2> cells;
        for ( auto nc0=nodes[i0].getOwners().begin(); nc0!=nodes[i0].getOwners().end(); ++nc0 ) {
            if ( std::find(nodes[i1].getOwners().begin(),
                           nodes[i1].getOwners().end(),
                           *nc0) != nodes[i1].getOwners().end()) {
                cells.push_back( *nc0 );
            }
        }
        if ( cells.size() == 1 ) {
            // we are on external edge
            return cells[0];
        }
        for ( auto nc0=nodes[nodeNo].getOwners().begin(); nc0!=nodes[nodeNo].getOwners().end(); ++nc0 ) {
            if ( *nc0 == cells[0] ) {
                return cells[1];
            } else if ( *nc0 == cells[1] ) {
                return cells[0];
            }
        }
        return std::numeric_limits<T2>::max();
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    T2 Grid2Duc<T1,T2,S,NODE,CELL>::findNextCell2(const T2 i0, const T2 i1, const T2 cellNo) const {
        std::vector<T2> cells;
        for ( auto nc0=nodes[i0].getOwners().begin(); nc0!=nodes[i0].getOwners().end(); ++nc0 ) {
            if ( std::find(nodes[i1].getOwners().begin(),
                           nodes[i1].getOwners().end(),
                           *nc0) != nodes[i1].getOwners().end()) {
                cells.push_back( *nc0 );
            }
        }
        if ( cells.size() == 1 ) {
            // we are on external edge
            return cells[0];
        }
        if ( cellNo == cells[0] ) {
            return cells[1];
        } else if ( cellNo == cells[1] ) {
            return cells[0];
        }
        return std::numeric_limits<T2>::max();
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    template<typename SetT>
    void Grid2Duc<T1,T2,S,NODE,CELL>::getNeighborNodes(const T2 cellNo,
                                                       SetT &nnodes) const {

        for ( size_t n=0; n<3; ++n ) {
            T2 nodeNo = this->neighbors[cellNo][n];
            nnodes.insert( &(nodes[nodeNo]) );

            for ( auto nc=nodes[nodeNo].getOwners().cbegin(); nc!=nodes[nodeNo].getOwners().cend(); ++nc ) {
                for ( size_t nn=0; nn<3; ++nn ) {
                    nnodes.insert( &(nodes[ this->neighbors[*nc][nn] ]) );
                }
            }
        }
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Duc<T1,T2,S,NODE,CELL>::calculateArea(std::vector<T1> &area) const {
        if ( area.size() != triangles.size() ) {
            area.resize( triangles.size() );
        }
        for ( size_t n=0; n<area.size(); ++n ) {
            area[n] =  abs( 0.5 * (nodes[triangles[n].i[0]].getX()*(nodes[triangles[n].i[1]].getY()-nodes[triangles[n].i[2]].getY()) +
                                   nodes[triangles[n].i[1]].getX()*(nodes[triangles[n].i[2]].getY()-nodes[triangles[n].i[0]].getY()) +
                                   nodes[triangles[n].i[2]].getX()*(nodes[triangles[n].i[0]].getY()-nodes[triangles[n].i[1]].getY())));
        }
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Duc<T1,T2,S,NODE,CELL>::interpolateAtNodes(std::vector<T1> &interpolated) const {
        // interpolated: values interpolated at nodes

        if ( interpolated.size() != nodes.size() ) {
            interpolated.resize( nodes.size() );
        }

        // calculate area of triangles
        static std::vector<T1> area;
        if ( area.size() == 0 ) {
            this->calculateArea(area);
        }

        for ( size_t n=0; n<nodes.size(); ++n ) {
            T1 totalArea = area[nodes[n].getOwners()[0]];
            interpolated[n] = cells.getSlowness(nodes[n].getOwners()[0]) * totalArea;
            for ( size_t nn=1; nn<nodes[n].getOwners().size(); ++nn ) {
                interpolated[n] += cells.getSlowness(nodes[n].getOwners()[nn]) * area[nodes[n].getOwners()[nn]];
                totalArea += area[nodes[n].getOwners()[nn]];
            }
            interpolated[n] /= totalArea;
        }
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    void Grid2Duc<T1,T2,S,NODE,CELL>::interpolateAtNodes(const std::vector<T1> &field,
                                                         std::vector<T1> &interpolated) const {
        // field: values defined at cells
        // interpolated: values interpolated at nodes

        if ( field.size() != triangles.size() ) {
            throw std::length_error("Error: field vector of incompatible size.");
        }
        if ( interpolated.size() != nodes.size() ) {
            interpolated.resize( nodes.size() );
        }

        // calculate area of triangles
        static std::vector<T1> area;
        if ( area.size() == 0 ) {
            this->calculateArea(area);
        }

        for ( size_t n=0; n<nodes.size(); ++n ) {
            T1 totalArea = area[nodes[n].getOwners()[0]];
            interpolated[n] = field[nodes[n].getOwners()[0]] * totalArea;
            for ( size_t nn=1; nn<nodes[n].getOwners().size(); ++nn ) {
                interpolated[n] += field[nodes[n].getOwners()[nn]] * area[nodes[n].getOwners()[nn]];
                totalArea += area[nodes[n].getOwners()[nn]];
            }
            interpolated[n] /= totalArea;
        }
    }

    template<typename T1, typename T2, typename S, typename NODE, typename CELL>
    const T1 Grid2Duc<T1,T2,S,NODE,CELL>::getAverageEdgeLength() const {
        std::set<std::array<T2,2>> edges;
        typename std::set<std::array<T2,2>>::iterator edgIt;
        T2 iNodes[3][2] = {
            {0,1},
            {0,2},
            {1,2}
        };
        T1 sum = 0.0;
        for (size_t ntri=0; ntri<triangles.size(); ++ntri) {
            for (size_t n=0; n<3; ++n) {
                std::array<T2, 2> edgei = {triangles[ntri].i[iNodes[n][0]],
                    triangles[ntri].i[iNodes[n][1]]};
                std::sort(edgei.begin(), edgei.end());
                edgIt = edges.find(edgei);
                if ( edgIt  == edges.end() ) {
                    T1 d = nodes[edgei[0]].getDistance(nodes[edgei[1]]);
                    sum += d;
                    edges.insert(edgei);
                }
            }
        }
        return (sum/edges.size());
    }

}

#endif

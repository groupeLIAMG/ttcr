//
//  Grid3Dun.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-04-21.
//  Copyright (c) 2014 Bernard Giroux. All rights reserved.
//  Copyright (c) 2018 Bernard Giroux, Maher Nasr. All rights reserved.
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
 * @file Grid3Dun.h
 * @brief Base class for 3-D unstructured tetrahedral meshes with node-based
 *        slowness.
 *
 * Declares ttcr::Grid3Dun, the mid-level base shared by every 3-D unstructured
 * solver whose slowness lives at the **nodes**: Grid3Dunsp (shortest path),
 * Grid3Dunfs (fast sweeping), Grid3Dunfm (fast marching) and Grid3Dundsp
 * (dynamic shortest path).
 *
 * It is the last corner of the 2x2 of 3-D bases: unstructured like
 * ttcr::Grid3Duc, node-slowness like ttcr::Grid3Drn. At roughly 12 000 lines it
 * is the largest file in the project, almost all of it raypath tracing.
 *
 * @section g3dun_vs_uc Node slowness versus cell slowness
 * As in every other family, the difference is where slowness lives and how a
 * traveltime increment is formed:
 *
 * - **@c un (here)** — one value per node. ttcr::Grid3Dun::computeDt averages
 *   the slowness at the two ends and multiplies by the distance,
 *   @f$\Delta t = \tfrac{1}{2}(s_1+s_2)\,\ell@f$.
 * - **@c uc (ttcr::Grid3Duc)** — one value per tetrahedron, drawn from a
 *   @c CELL policy that may be anisotropic.
 *
 * So this family has no @c CELL template parameter and is isotropic, and its
 * slowness field is continuous across element boundaries rather than piecewise
 * constant. Sampling that continuous field is what
 * ttcr::Grid3Dun::computeSlowness does, in two forms: at an arbitrary point
 * (locate the tetrahedron, then interpolate over its four primary nodes), or at
 * a point already known to lie on a node, an edge or a face — the form the
 * raypath tracer uses, since it always knows which.
 *
 * @section g3dun_procvel Velocity versus slowness interpolation
 * Interpolating slowness and interpolating velocity give different answers,
 * because the mean of the reciprocals is not the reciprocal of the mean. The
 * @c processVel flag, set from ttcr::input_parameters::processVel, selects
 * which: when true the node velocities are interpolated and the result
 * inverted; when false slowness is interpolated directly. The distinction runs
 * right through the file — it is tested at over a hundred sites, wherever a
 * value is needed between nodes.
 *
 * This mirrors @ref g3drn_procvel in the rectilinear node-slowness family. Only
 * the node-slowness families offer the choice; with cell slowness there is
 * nothing to interpolate.
 *
 * @section g3dun_update The local update, and an unused alternative
 * ttcr::Grid3Dun::localUpdate3D computes the traveltime at a tetrahedron's
 * fourth vertex from the other three (Lelievre et al., 2011), falling back to
 * ttcr::Grid3Dun::localUpdate2D on each face when the three-dimensional stencil
 * is not causal. This is what the sweeping and marching solvers call.
 *
 * There is no obtuse-triangle machinery here, unlike the 2-D unstructured
 * families (@ref g2duc_obtuse): the face-wise fallback plays that role.
 *
 * A second, older stencil — ttcr::Grid3Dun::local3Dsolver,
 * ttcr::Grid3Dun::local2Dsolver and ttcr::Grid3Dun::solveEq23, after Qian et
 * al. (2007) — is present but **dead**: nothing outside the three of them calls
 * any of them. The same trio, likewise dead, appears in ttcr::Grid3Duc; see
 * @ref g3duc_update.
 *
 * @section g3dun_raypath Raypath tracing
 * Raypaths are traced backwards from receiver to source, stepping cell by cell
 * along the negated traveltime gradient. The gradient estimator is chosen by
 * the @c rp_method constructor argument:
 *
 * | @c rp_method | estimator |
 * | :----------: | :-------- |
 * | 0 | ttcr::Grad3D_ls_fo — least squares, first order |
 * | 1 | ttcr::Grad3D_ls_so — least squares, second order |
 * | 2 | ttcr::Grad3D_ab — averaged over the cells sharing the point |
 *
 * Methods 0 and 1 estimate a gradient at a point; method 2 differs enough that
 * several code paths branch on @c rp_method @c < @c 2.
 *
 * The tracer tracks whether the current point sits on a node, on an edge or on
 * a face, since each case has its own way of finding the next cell.
 * @c min_dist, from ttcr::input_parameters::min_distance_rp, is the tolerance
 * for snapping an intersection point onto a nearby node rather than leaving it
 * a hair away on the edge.
 *
 * @section g3dun_blti An unused second raypath implementation
 * The file also carries a complete alternative raypath implementation, the
 * @c blti group (Nasr, 2018): ttcr::Grid3Dun::getTraveltime_blti,
 * ttcr::Grid3Dun::getRaypath_blti and the helpers @c blti_raytrace,
 * @c blti2D_raytrace and @c blti_solver_around_source, together with their own
 * @c txInfoCacheBlti / @c getTxInfoBlti source cache. Rather than following a
 * gradient it solves for the ray's exit point on each face directly.
 *
 * It is **dead code**: neither entry point is called anywhere in the C++ sources
 * or the Cython bindings. The group spans roughly 1 900 lines, about 15 % of the
 * file.
 *
 * @section g3dun_txinfo Source classification cache
 * Like ttcr::Grid3Duc, this class caches the classification of the Tx points —
 * whether each lies on a node, an edge or a face, and which cells own it — in
 * @c txInfoCache, indexed by thread number so no locking is needed. The setup
 * scans every node, so doing it once per source rather than once per
 * (source, receiver) pair matters. See @ref g3duc_txinfo.
 *
 * @section g3dun_tomo Tomography operators
 * Two methods build matrices for inversion rather than forward modelling:
 *
 * - ttcr::Grid3Dun::computeD — interpolation weights mapping node slowness onto
 *   arbitrary points, one row per point.
 * - ttcr::Grid3Dun::computeK — first- or second-order spatial derivative
 *   operators at every primary node, from a weighted least-squares Taylor fit
 *   over the surrounding nodes. Used to regularise an inversion.
 *
 * @sa Grid3D.h, Grid3Duc.h, Grid3Drn.h, Grid2Dun.h, Grid3Dunsp.h
 */

#ifndef ttcr_Grid3Dun_h
#define ttcr_Grid3Dun_h

#include <cassert>
#include <cmath>
#include <cstddef>

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <memory_resource>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#ifdef VTK
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkProbeFilter.h"
#include "vtkRectilinearGrid.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSmartPointer.h"
#include "vtkTetra.h"
#include "vtkXMLRectilinearGridWriter.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkXMLUnstructuredGridWriter.h"
#endif

#include <Eigen/Dense>

#include "Grad.h"
#include "Grid3D.h"
#include "Interpolator.h"
#include "NodeKDTree.h"
#include "utils.h"

namespace ttcr {

    /**
     * @brief 3-D tetrahedral mesh holding one slowness value per node.
     *
     * Holds the geometry (nodes and tetrahedra), the slowness field, and
     * everything the derived solvers share: point location, the local
     * traveltime update, raypath tracing and the tomography operators. It does
     * not itself propagate a wavefront — that is what the four derived classes
     * differ in.
     *
     * @tparam T1 floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2 unsigned integer type used for node and cell indices.
     * @tparam NODE node type; ttcr::Node3Dnsp for the shortest-path solver
     *              (it carries the graph edges), ttcr::Node3Dn for the others.
     *
     * @sa Grid3Dunsp.h, Grid3Dunfs.h, Grid3Dunfm.h, Grid3Dundsp.h
     */
    template<typename T1, typename T2, typename NODE>
    class Grid3Dun : public Grid3D<T1,T2> {
    public:
        /**
         * @brief Build the mesh from its nodes and tetrahedra.
         * @param no  primary node coordinates.
         * @param tet tetrahedra, given as quadruples of indices into @p no.
         * @param rp  raypath gradient estimator, 0, 1 or 2; see
         *            @ref g3dun_raypath.
         * @param procVel interpolate velocity rather than slowness; see
         *            @ref g3dun_procvel.
         * @param ttrp compute traveltimes by integrating along the raypath
         *            instead of reading them off the grid.
         * @param md  minimum distance for snapping a raypath point onto a node.
         * @param nt  number of threads; per-thread storage is sized from it.
         * @param _translateOrigin shift the mesh so its lower corner sits at the
         *            origin, which helps conditioning when coordinates are large.
         *
         * @note Only the primary nodes are created here. Secondary nodes are
         *       added by the derived class, which alone knows how many it wants;
         *       see ttcr::Grid3Dun::buildGridNodes.
         */
        Grid3Dun(const std::vector<sxyz<T1>>& no,
                 const std::vector<tetrahedronElem<T2>>& tet,
                 const int rp,
                 const bool procVel,
                 const bool ttrp,
                 const T1 md,
                 const size_t nt=1,
                 const bool _translateOrigin=false) :
        Grid3D<T1,T2>(ttrp, tet.size(), nt, _translateOrigin),
        rp_method(rp), processVel(procVel),
        nPrimary(static_cast<T2>(no.size())),
        source_radius(0.0), min_dist(md),
        nodes(std::vector<NODE>(no.size(), NODE(nt))),
        tetrahedra(tet),
        txInfoCache(nt),
        txInfoCacheBlti(nt)
        {}

        virtual ~Grid3Dun() {}

        /**
         * @name Slowness accessors
         * @{
         */

        /**
         * @brief Give every node the same slowness.
         * @param s the uniform value.
         * @note Applies to **all** nodes, secondary ones included, so it needs
         *       no companion interpolation step — every node ends up with the
         *       same value whichever kind it is.
         */
        void setSlowness(const T1 s) {
            for ( size_t n=0; n<nodes.size(); ++n ) {
                nodes[n].setNodeSlowness( s );
            }
        }

        /**
         * @brief Set the slowness at every node, from a raw array.
         * @param s  array of slowness values, one per node.
         * @param ns length of @p s.
         * @throws std::length_error if @p ns differs from the node count.
         *
         * @note The count must cover **all** nodes, secondary ones included, so
         *       this is not the overload a caller working in terms of the
         *       primary nodes wants; the derived classes provide that.
         */
        void setSlowness(const T1 *s, const size_t ns) {
            if ( nodes.size() != ns ) {
                throw std::length_error("Error: slowness vectors of incompatible size.");
            }
            for ( size_t n=0; n<nodes.size(); ++n ) {
                nodes[n].setNodeSlowness( s[n] );
            }
        }

        /**
         * @brief Set the slowness at every node.
         * @param s one value per node, secondary nodes included.
         * @throws std::length_error if @p s does not match the node count.
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
         * @brief Retrieve the slowness at the primary nodes.
         * @param[out] slowness resized to the primary node count if needed.
         *
         * @note Asymmetric with the setters: this returns only the primary
         *       values, being the model parameters, while they take all nodes.
         */
        void getSlowness(std::vector<T1>& slowness) const {
            if (slowness.size() != nPrimary) {
                slowness.resize(nPrimary);
            }
            for ( size_t n=0; n<slowness.size(); ++n ) {
                slowness[n] = nodes[n].getNodeSlowness();
            }
        }
        /** @} */

        /**
         * @brief Set the radius around a source used to seed the initial
         *        traveltimes.
         * @param r radius; 0 (the default) seeds only the immediate neighbours.
         *
         * Controls how widely the wavefront is initialised before sweeping or
         * marching begins. Every node within @p r of the source is given the
         * straight-ray estimate @f$t_0 + \tfrac{1}{2}(s_{tx}+s_n)\,d@f$, kept
         * only if it improves on the node's current value. With @p r left at
         * zero the seeding reaches just the nodes sharing a tetrahedron with the
         * source.
         *
         * Widening it helps where the wavefront curvature near a point source is
         * too strong for the stencils to resolve. Only the source node itself is
         * frozen; the seeded nodes are ordinary nodes that the solver may still
         * lower.
         *
         * @note Used by the sweeping and marching solvers. The shortest-path
         *       solvers ignore it, and ttcr::Grid3Dundsp instead sets it to its
         *       own tertiary-node radius.
         */
        void setSourceRadius(const double r) { source_radius = r; }

        /**
         * @brief Set the traveltime at one node.
         * @param tt the value.
         * @param nn node index.
         * @param nt thread number.
         */
        void setTT(const T1 tt, const size_t nn, const size_t nt=0) {
            nodes[nn].setTT(tt, nt);
        }

        /**
         * @brief Total node count, secondary nodes included.
         * @return the number of nodes.
         * @note Compare @c nPrimary, the count of mesh vertices alone.
         */
        size_t getNumberOfNodes() const { return nodes.size(); }

        /**
         * @brief Retrieve the traveltimes at the primary nodes.
         * @param[out] tt resized to the primary node count.
         * @param threadNo thread number.
         * @note Secondary nodes are omitted: they are an artefact of the
         *       discretisation, not part of the model.
         */
        void getTT(std::vector<T1>& tt, const size_t threadNo=0) const final {
            tt.resize(nPrimary);
            for ( size_t n=0; n<nPrimary; ++n ) {
                tt[n] = nodes[n].getTT(threadNo);
            }
        }

        /**
         * @name Bounding box
         * Extent of the mesh, computed by scanning every node on each call.
         * @{
         */
        /** @brief Smallest node x coordinate. */
        const T1 getXmin() const {
            T1 xmin = nodes[0].getX();
            for ( auto it=nodes.begin(); it!=nodes.end(); ++it ) {
                xmin = xmin<it->getX() ? xmin : it->getX();
            }
            return xmin;
        }
        /** @brief Largest node x coordinate. */
        const T1 getXmax() const {
            T1 xmax = nodes[0].getX();
            for ( auto it=nodes.begin(); it!=nodes.end(); ++it ) {
                xmax = xmax>it->getX() ? xmax : it->getX();
            }
            return xmax;
        }
        /** @brief Smallest node y coordinate. */
        const T1 getYmin() const {
            T1 ymin = nodes[0].getY();
            for ( auto it=nodes.begin(); it!=nodes.end(); ++it ) {
                ymin = ymin<it->getY() ? ymin : it->getY();
            }
            return ymin;
        }
        /** @brief Largest node y coordinate. */
        const T1 getYmax() const {
            T1 ymax = nodes[0].getY();
            for ( auto it=nodes.begin(); it!=nodes.end(); ++it ) {
                ymax = ymax>it->getY() ? ymax : it->getY();
            }
            return ymax;
        }
        /** @brief Smallest node z coordinate. */
        const T1 getZmin() const {
            T1 zmin = nodes[0].getZ();
            for ( auto it=nodes.begin(); it!=nodes.end(); ++it ) {
                zmin = zmin<it->getZ() ? zmin : it->getZ();
            }
            return zmin;
        }
        /** @brief Largest node z coordinate. */
        const T1 getZmax() const {
            T1 zmax = nodes[0].getZ();
            for ( auto it=nodes.begin(); it!=nodes.end(); ++it ) {
                zmax = zmax>it->getZ() ? zmax : it->getZ();

            }
            return zmax;
        }
        /** @} */

        /**
         * @brief Write the traveltime field to disk.
         * @param fname  base filename; the extension follows from @p format.
         * @param all    if 1, include the secondary nodes.
         * @param nt     thread number.
         * @param format output format:
         *               | value | file | contents |
         *               | :---: | :--- | :------- |
         *               | 1 | @c .dat | ASCII, one @c "x y z t" line per node |
         *               | 2 | @c .vtu | VTK unstructured grid |
         *               | 3 | @c .bin | raw binary |
         *
         * @note @p all is honoured by formats 1 and 3 only; the VTK output
         *       always writes just the primary nodes, since the secondary ones
         *       are not vertices of any cell.
         */
        void saveTT(const std::string &fname, const int all, const size_t nt=0,
                    const int format=1) const;
        /**
         * @brief Read back a traveltime field written by @ref saveTT.
         * @param fname  base filename.
         * @param all    if 1, expect the secondary nodes too.
         * @param nt     thread number.
         * @param format must match the one used to write.
         *
         * @note Declared @c const yet it writes into the nodes; it can do so
         *       because @c nodes is @c mutable.
         */
        void loadTT(const std::string &fname, const int all, const size_t nt=0,
                    const int format=1) const;

#ifdef VTK
        /**
         * @brief Write the mesh and its slowness model to a VTK @c .vtu file.
         * @param fname filename.
         * @param saveSlowness write slowness; when false, velocity is written.
         * @param savePhysicalEntity also write the tetrahedra's physical-entity
         *        tags, as carried over from the Gmsh mesh.
         */
        void saveModelVTU(const std::string &fname, const bool saveSlowness=true,
                          const bool savePhysicalEntity=false) const;
#endif

        /**
         * @brief Write the coordinates of the secondary nodes, for debugging.
         * @param os open output stream; one @c "x y z" line per node.
         */
        void dump_secondary(std::ofstream& os) const {
            for ( size_t n=nPrimary; n<nodes.size(); ++n ) {
                os << nodes[n].getX() << ' ' << nodes[n].getY() << ' ' << nodes[n].getZ() << '\n';
            }
        }

        /**
         * @brief Build the interpolation weights mapping node slowness onto a
         *        set of points.
         * @param pts        the points; taken by value since they may be
         *                   translated in place.
         * @param[out] d_data one row per point, each a list of
         *                   (point, node, weight) triplets.
         * @param translated set when @p pts are already in translated
         *                   coordinates, so the origin shift is not applied
         *                   twice.
         *
         * A point coinciding with a node gets a single unit weight; otherwise
         * the weights are the barycentric coordinates within its tetrahedron.
         * Used to relate a measurement to the model parameters in tomography.
         *
         * @sa computeK
         */
        void computeD(std::vector<sxyz<T1>> pts,
                      std::vector<std::vector<sijv<T1>>> &d_data,
                      const bool translated=false) const;

        /**
         * @brief Build spatial derivative operators at every primary node.
         * @param[out] k_data three matrices, for @f$\partial/\partial x@f$,
         *             @f$\partial/\partial y@f$ and @f$\partial/\partial z@f$,
         *             each with one row of (node, weight) pairs per primary node.
         * @param order derivative order, 1 or 2.
         * @param taylorSeriesOrder order of the Taylor expansion fitted, 1 or 2.
         * @param weighting weight the neighbours by inverse distance.
         * @param s0inside treat the centre node's own value as an unknown of the
         *        fit rather than as fixed; selects the @c buildA2 /
         *        @c fill_k_data2 variants.
         * @param additionnalPoints widen the neighbourhood beyond the minimum
         *        needed for the fit, which helps where the mesh is irregular.
         *
         * At each node the surrounding nodes are collected and a weighted
         * least-squares Taylor fit is solved for the derivative coefficients.
         * The resulting operators regularise an inversion by penalising
         * roughness.
         *
         * @throws std::runtime_error if @p order or @p taylorSeriesOrder is
         *         outside {1, 2}, or if a second-order derivative is asked of a
         *         first-order expansion, which cannot supply it.
         * @sa computeD, getSurroundingNodes, buildA
         */
        void computeK(std::vector<std::vector<std::vector<siv<T1>>>>& k_data,
                      const int order=2, const int taylorSeriesOrder=2,
                      const bool weighting=1, const bool s0inside=0,
                      const int additionnalPoints=0) const;

        /**
         * @brief Mean length of the tetrahedron edges.
         * @return the average.
         *
         * A characteristic length for the mesh, used by ttcr::Grid3Dundsp to
         * express its tertiary-node radius as a multiple of the local element
         * size rather than as an absolute distance.
         */
        const T1 getAverageEdgeLength() const;

    protected:
        int rp_method;   ///< Raypath gradient estimator; see @ref g3dun_raypath.
        bool processVel; ///< Interpolate velocity rather than slowness; see @ref g3dun_procvel.
        T2 nPrimary;     ///< Number of primary nodes, i.e. mesh vertices. They occupy the first @c nPrimary slots of @c nodes.
        T1 source_radius;///< Radius for seeding traveltimes around a source; see @ref setSourceRadius.
        T1 min_dist;     ///< Tolerance for snapping a raypath point onto a node.
        mutable std::vector<NODE> nodes;             ///< Primary nodes first, then the secondary ones. @c mutable so that @c const methods can write traveltimes.
        std::vector<tetrahedronElem<T2>> tetrahedra; ///< The elements, each holding four primary node indices.

        // kd-tree over the primary nodes, for fast point location in getCellNo.
        // Built lazily on first use (std::call_once) since node coordinates are
        // fixed after construction; queries are read-only and thread-safe.
        mutable std::unique_ptr<NodeKDTree<T1,T2>> kdtree;
        mutable std::once_flag kdtreeFlag;

        T2 getNearestNode(const sxyz<T1>& pt) const {
            std::call_once(kdtreeFlag, [this]() {
                kdtree.reset(new NodeKDTree<T1,T2>(nodes, nPrimary));
            });
            return kdtree->findNearest(pt.x, pt.y, pt.z);
        }

        /**
         * @brief Cached classification of the Tx points.
         *
         * For a given source the Tx are constant while looping over all
         * receivers, so this setup -- which scans every node and calls
         * @ref getCellNo per Tx -- need only be done once per source rather than
         * once per (source, receiver) pair. See @ref g3dun_txinfo.
         */
        struct txInfo_t {
            std::vector<sxyz<T1>> tx;        ///< The Tx the rest of the struct describes; compared against to detect staleness.
            std::vector<bool> txOnNode;      ///< Whether each Tx coincides with a node.
            std::vector<bool> txOnEdge;      ///< Whether each Tx lies on an edge.
            std::vector<bool> txOnFace;      ///< Whether each Tx lies on a face.
            std::vector<T2> txNode;          ///< Node index, where @c txOnNode.
            std::vector<T2> txCell;          ///< Index of the tetrahedron containing each Tx.
            std::vector<std::array<T2,2>> txEdges; ///< The two nodes of the edge, where @c txOnEdge.
            std::vector<std::array<T2,3>> txFaces; ///< The three nodes of the face, where @c txOnFace.
            std::vector<std::vector<T2>> txNeighborCells; ///< Cells sharing a node with each Tx.
            bool valid = false;              ///< Set once the slot has been filled.
        };
        /**
         * @brief Per-thread Tx classification cache.
         * Indexed by thread number, so parallel runs need no locking: each
         * thread owns a distinct slot. A slot is recomputed only when its Tx
         * changes.
         */
        mutable std::vector<txInfo_t> txInfoCache;
        /**
         * @brief Separate cache for the @c blti raypath methods, whose Tx setup
         *        differs: it matches any (primary or secondary) node and
         *        computes only the owning cell and neighbour cells, with no
         *        edge/face classification.
         * @note Dead along with the rest of the @c blti group; see
         *       @ref g3dun_blti.
         */
        mutable std::vector<txInfo_t> txInfoCacheBlti;

        /**
         * @brief Classify the Tx points, reusing the cached result when possible.
         * @param Tx       source coordinates.
         * @param threadNo thread number, selecting the cache slot.
         * @return the cache entry, recomputed only if @p Tx differs from the one
         *         it holds.
         */
        const txInfo_t& getTxInfo(const std::vector<sxyz<T1>>& Tx,
                                  const size_t threadNo) const;

        /**
         * @brief Tx classification for the @c blti methods.
         * @param Tx       source coordinates.
         * @param threadNo thread number.
         * @return the cache entry.
         * @note Dead code; see @ref g3dun_blti.
         */
        const txInfo_t& getTxInfoBlti(const std::vector<sxyz<T1>>& Tx,
                                      const size_t threadNo) const;

        /**
         * @brief Traveltime increment between two nodes.
         * @param source the node departed from.
         * @param node   the node arrived at.
         * @return the trapezoidal rule along the straight segment,
         *         @f$\tfrac{1}{2}(s_1+s_2)\,\ell@f$.
         *
         * This averaging of the endpoint slownesses is what distinguishes the
         * node-slowness families from the cell-slowness ones, where the single
         * cell value applies; see @ref g3dun_vs_uc.
         */
        T1 computeDt(const NODE& source, const NODE& node) const {
            return (node.getNodeSlowness()+source.getNodeSlowness())/2 * source.getDistance( node );
        }

        /**
         * @brief Traveltime increment from a node to an arbitrary point.
         * @param source the node departed from.
         * @param node   the point arrived at.
         * @param slo    slowness at @p node, which the caller must supply since
         *               a bare point carries none.
         * @return the increment @f$\tfrac{1}{2}(s_1+s_2)\,\ell@f$.
         */
        T1 computeDt(const NODE& source, const sxyz<T1>& node, T1 slo) const {
            return (slo+source.getNodeSlowness())/2 * source.getDistance( node );
        }

        /**
         * @brief Traveltime at an arbitrary point, interpolated from the nodes.
         * @param Rx       the point.
         * @param threadNo thread number.
         * @return the traveltime; read directly if @p Rx sits on a node,
         *         otherwise interpolated over its tetrahedron.
         */
        T1 getTraveltime(const sxyz<T1>& Rx,
                         const size_t threadNo) const;

        /**
         * @brief Traveltime at a point, using a caller-supplied node set.
         * @param Rx       the point.
         * @param nodes    the nodes to interpolate from.
         * @param threadNo thread number.
         * @return the traveltime.
         *
         * @note The @p nodes parameter shadows the member of the same name. It
         *       lets ttcr::Grid3Dundsp pass its augmented set, which includes
         *       the tertiary nodes inserted around the source.
         */
        T1 getTraveltime(const sxyz<T1>& Rx,
                         const std::vector<NODE>& nodes,
                         const size_t threadNo) const;

        /**
         * @brief Verify that points lie within the mesh.
         * @param pts        the points; taken by value as they may be
         *                   translated in place.
         * @param translated set when @p pts are already in translated
         *                   coordinates.
         * @throws std::runtime_error naming the offending point if one falls
         *         outside every tetrahedron.
         */
        void checkPts(std::vector<sxyz<T1>> pts, const bool translated=false) const;

        /**
         * @brief Whether a point lies inside a tetrahedron, by determinants.
         * @param pt      the point.
         * @param tetraNo the tetrahedron.
         * @return true if inside.
         *
         * Forms the five signed volumes @f$D_0 \ldots D_4@f$ obtained by
         * substituting @p pt for each vertex in turn; the point is inside when
         * they all share the sign of @f$D_0@f$. A zero among them puts the point
         * on a face, an edge or a vertex.
         *
         * @sa insideTetrahedron2, which decides the same question differently.
         */
        bool insideTetrahedron(const sxyz<T1>& pt, const T2 tetraNo) const;
        /**
         * @brief Whether a point lies inside a tetrahedron, by face tests.
         * @param pt      the point.
         * @param tetraNo the tetrahedron.
         * @return true if inside.
         *
         * Asks of each of the four faces whether @p pt falls on the same side as
         * the opposite vertex.
         */
        bool insideTetrahedron2(const sxyz<T1>& pt, const T2 tetraNo) const;

        /**
         * @brief Find the tetrahedron containing a point.
         * @param pt the point.
         * @return the cell index.
         *
         * Starts from the nearest primary node — a kd-tree query — and searches
         * only the tetrahedra owning it, keeping the one whose split volumes
         * best match its own. Restricting the search this way is what keeps
         * point location off the critical path.
         */
        T2 getCellNo(const sxyz<T1>& pt) const;

        /**
         * @name Node construction
         * Called from a derived class's constructor, which alone knows how many
         * secondary nodes it wants.
         * @{
         */
        /**
         * @brief Build the primary nodes only.
         * @param no coordinates of the mesh vertices.
         * @param nt number of threads.
         */
        void buildGridNodes(const std::vector<sxyz<T1>>& no, const size_t nt);
        /**
         * @brief Build the primary nodes and insert secondary ones.
         * @param no          coordinates of the mesh vertices.
         * @param nsecondary  number of secondary nodes per edge.
         * @param nt          number of threads.
         *
         * Secondary nodes are placed along the edges and across the faces,
         * refining the discretisation without changing the mesh. They follow the
         * primary nodes in @c nodes, so the first @c nPrimary slots keep their
         * meaning.
         */
        void buildGridNodes(const std::vector<sxyz<T1>>& no,
                            const T2 nsecondary, const size_t nt);
        /** @} */

        /**
         * @brief Update the traveltime at one vertex from a neighbouring
         *        tetrahedron.
         * @param vertexC  the node to update; lowered only if the new value is
         *                 smaller.
         * @param threadNo thread number.
         *
         * The three-dimensional local solver of Lelievre et al. (2011): for each
         * tetrahedron owning @p vertexC, the traveltime is extrapolated from the
         * other three vertices. Where that stencil is not causal — the incoming
         * ray would enter through the wrong face — it falls back to
         * @ref localUpdate2D on each face in turn.
         *
         * This is the workhorse of the sweeping and marching solvers; see
         * @ref g3dun_update.
         */
        void localUpdate3D(NODE *vertexC, const size_t threadNo) const;

        /**
         * @brief Update a vertex from one triangular face.
         * @param vertexA  first known vertex.
         * @param vertexB  second known vertex.
         * @param vertexC  the vertex being updated.
         * @param tetraNo  the tetrahedron the face belongs to.
         * @param threadNo thread number.
         * @return the candidate traveltime, or infinity if the face yields no
         *         causal solution.
         *
         * The two-dimensional fallback used by @ref localUpdate3D.
         */
        T1 localUpdate2D(const NODE *vertexA,
                         const NODE *vertexB,
                         const NODE *vertexC,
                         const T2 tetraNo,
                         const size_t threadNo) const;

        /**
         * @name Unused local solver (Qian et al., 2007)
         * An older alternative to @ref localUpdate3D. **Dead code**: nothing
         * outside these three functions calls any of them. The same trio,
         * likewise dead, appears in ttcr::Grid3Duc. See @ref g3dun_update.
         * @{
         */
        /**
         * @brief Update a vertex by the Qian stencil.
         * @param vertexD  the node to update.
         * @param threadNo thread number.
         */
        void local3Dsolver(NODE *vertexD, const size_t threadNo) const;

        /**
         * @brief Face-wise fallback for @ref local3Dsolver.
         * @param vertexA  first known vertex.
         * @param vertexB  second known vertex.
         * @param vertexC  the vertex being updated.
         * @param tetraNo  the tetrahedron.
         * @param threadNo thread number.
         * @return the candidate traveltime.
         */
        T1 local2Dsolver(const NODE *vertexA,
                         const NODE *vertexB,
                         const NODE *vertexC,
                         const T2 tetraNo,
                         const size_t threadNo) const;

        /**
         * @brief Solve the quadratic system arising in @ref local3Dsolver.
         * @param a      first coefficient triplet.
         * @param b      second coefficient triplet.
         * @param[out] n up to two solution vectors.
         * @return the number of real solutions found, 0, 1 or 2.
         */
        int solveEq23(const T1 a[], const T1 b[], T1 n[][3]) const;
        /** @} */

        /**
         * @brief Traveltime obtained by integrating slowness along the raypath.
         * @param Tx       source coordinates.
         * @param t0       source excitation times.
         * @param Rx       receiver coordinate.
         * @param threadNo thread number.
         * @return the integrated traveltime.
         *
         * An alternative to reading the traveltime off the grid, selected by
         * ttcr::input_parameters::tt_from_rp. It traces the raypath and
         * accumulates @f$\int s\,d\ell@f$ along it, which avoids the
         * interpolation error of sampling the traveltime field at a point that
         * is not a node, at the cost of tracing a ray per receiver.
         */
        T1 getTraveltimeFromRaypath(const std::vector<sxyz<T1>>& Tx,
                                    const std::vector<T1>& t0,
                                    const sxyz<T1>& Rx,
                                    const size_t threadNo) const;

        // keep the base-class getRaypath overloads visible; the helper
        // overloads below would otherwise hide them (-Woverloaded-virtual)
        using Grid3D<T1,T2>::getRaypath;

        /**
         * @name Raypath tracing
         * Five overloads over the same traversal, differing only in what they
         * accumulate: the ray geometry (@p r_data), the tomography matrix row
         * (@p m_data), the traveltime (@p tt), or some combination. All trace
         * backwards from @p Rx to the nearest @p Tx along the negated
         * traveltime gradient; see @ref g3dun_raypath.
         * @{
         */
        /**
         * @brief Trace a raypath.
         * @param Tx            source coordinates.
         * @param Rx            receiver coordinate.
         * @param[out] r_data   the raypath, from receiver to source.
         * @param threadNo      thread number.
         */
        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const sxyz<T1>& Rx,
                        std::vector<sxyz<T1>>& r_data,
                        const size_t threadNo) const;

        /**
         * @brief Trace a raypath and integrate the traveltime along it.
         * @param Tx            source coordinates.
         * @param t0            source excitation times.
         * @param Rx            receiver coordinate.
         * @param[out] r_data   the raypath.
         * @param[out] tt       the integrated traveltime.
         * @param threadNo      thread number.
         */
        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const sxyz<T1>& Rx,
                        std::vector<sxyz<T1>>& r_data,
                        T1 &tt,
                        const size_t threadNo) const;

        /**
         * @brief Accumulate a tomography matrix row and the traveltime,
         *        discarding the ray geometry.
         * @param Tx            source coordinates.
         * @param t0            source excitation times.
         * @param Rx            receiver coordinate.
         * @param[out] m_data   (row, node, weight) triplets giving the ray's
         *                      sensitivity to each node slowness.
         * @param RxNo          receiver index, used as the row number.
         * @param[out] tt       the integrated traveltime.
         * @param threadNo      thread number.
         */
        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const sxyz<T1>& Rx,
                        std::vector<sijv<T1>>& m_data,
                        const size_t RxNo,
                        T1 &tt,
                        const size_t threadNo) const;

        /**
         * @brief Trace a raypath and accumulate a tomography matrix row.
         * @param Tx            source coordinates.
         * @param Rx            receiver coordinate.
         * @param[out] r_data   the raypath.
         * @param[out] m_data   the matrix row.
         * @param RxNo          receiver index.
         * @param threadNo      thread number.
         */
        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const sxyz<T1>& Rx,
                        std::vector<sxyz<T1>>& r_data,
                        std::vector<sijv<T1>>& m_data,
                        const size_t RxNo,
                        const size_t threadNo) const;

        /**
         * @brief Trace a raypath, accumulating both the matrix row and the
         *        traveltime.
         * @param Tx            source coordinates.
         * @param t0            source excitation times.
         * @param Rx            receiver coordinate.
         * @param[out] r_data   the raypath.
         * @param[out] m_data   the matrix row.
         * @param RxNo          receiver index.
         * @param[out] tt       the integrated traveltime.
         * @param threadNo      thread number.
         */
        void getRaypath(const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const sxyz<T1>& Rx,
                        std::vector<sxyz<T1>>& r_data,
                        std::vector<sijv<T1>>& m_data,
                        const size_t RxNo,
                        T1 &tt,
                        const size_t threadNo) const;
        /** @} */

        /**
         * @name Unused raypath implementation (blti)
         * A second, self-contained raypath implementation that solves for the
         * ray's exit point on each face rather than following a gradient.
         * **Dead code**: neither entry point is called anywhere in the C++
         * sources or the Cython bindings. See @ref g3dun_blti.
         * @{
         */
        /**
         * @brief Traveltime along a @c blti raypath.
         * @param Tx       source coordinates.
         * @param t0       source excitation times.
         * @param Rx       receiver coordinate.
         * @param threadNo thread number.
         * @return the integrated traveltime.
         */
        T1 getTraveltime_blti(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const sxyz<T1>& Rx,
                              const size_t threadNo) const;

        /**
         * @brief Trace a @c blti raypath.
         * @param Tx            source coordinates.
         * @param t0            source excitation times.
         * @param Rx            receiver coordinate.
         * @param[out] r_data   the raypath.
         * @param[out] tt       the integrated traveltime.
         * @param threadNo      thread number.
         */
        void getRaypath_blti(const std::vector<sxyz<T1>>& Tx,
                             const std::vector<T1>& t0,
                             const sxyz<T1>& Rx,
                             std::vector<sxyz<T1>>& r_data,
                             T1 &tt,
                             const size_t threadNo) const;

        /**
         * @brief Locate the ray's crossing of a face near the source.
         * @param Source            the source point.
         * @param curr_pt           current point on the ray.
         * @param face              the three nodes of the face.
         * @param[out] barycenters  barycentric coordinates of the crossing.
         * @return true if the straight segment from @p curr_pt to @p Source
         *         crosses @p face.
         */
        bool blti_solver_around_source(const sxyz<T1>& Source,
                                       const sxyz<T1>& curr_pt,
                                       const std::array<T2,3>& face,
                                       std::array<T1,3>& barycenters) const;
        /**
         * @brief Advance the ray across one face.
         * @param curr_pt       current point.
         * @param faces         the three nodes of the face.
         * @param[out] next_pt  the exit point.
         * @param threadNo      thread number.
         * @param s             slowness at @p curr_pt.
         * @return true if a valid exit point was found within the face.
         */
        bool blti_raytrace(const sxyz<T1>& curr_pt,
                           const std::array<T2, 3>& faces,
                           sxyz<T1>& next_pt,
                           const size_t threadNo,
                           const T1& s) const;
        /**
         * @brief Advance the ray along one edge.
         * @param curr_pt       current point.
         * @param node1         first node of the edge.
         * @param node2         second node of the edge.
         * @param[out] next_pt  the exit point.
         * @param threadNo      thread number.
         * @param s             slowness at @p curr_pt.
         * @return true if a valid exit point was found.
         */
        bool blti2D_raytrace(const sxyz<T1>& curr_pt,
                             const T2& node1,
                             const T2& node2,
                             sxyz<T1>& next_pt,
                             const size_t threadNo,
                             const T1& s) const;
        /** @} */

        /**
         * @name Locating a raypath point
         * The tracer must know at every step whether the current point sits on a
         * node, on an edge or on a face, since each case offers a different set
         * of cells to step into.
         * @{
         */
        /**
         * @brief Classify a point against the nodes of one face.
         * @param[in,out] curr_pt the point; snapped onto a node or edge when
         *                within @c min_dist of one.
         * @param ind            the three nodes to test against.
         * @param[out] onNode    set when the point coincides with a node.
         * @param[out] nodeNo    that node, when @p onNode.
         * @param[out] onEdge    set when the point lies on an edge.
         * @param[out] edgeNodes that edge's two nodes, when @p onEdge.
         * @param[out] onFace    set when the point lies on a face.
         * @param[out] faceNodes that face's three nodes, when @p onFace.
         * @return true if the point was located.
         */
        bool check_pt_location(sxyz<T1>& curr_pt,
                               const std::array<T2,3>& ind,
                               bool& onNode,
                               T2& nodeNo,
                               bool& onEdge,
                               std::array<T2,2>& edgeNodes,
                               bool& onFace,
                               std::array<T2,3>& faceNodes) const;

        /**
         * @brief Classify a point against a wider set of candidates.
         * @param[in,out] curr_pt the point.
         * @param ind1           candidate nodes to test first.
         * @param ind2           the three nodes of a face, tested after.
         * @param[out] onNode    set when the point coincides with a node.
         * @param[out] nodeNo    that node.
         * @param[out] onEdge    set when the point lies on an edge.
         * @param[out] edgeNodes that edge's two nodes.
         * @param[out] onFace    set when the point lies on a face.
         * @param[out] faceNodes that face's three nodes.
         * @return true if the point was located.
         */
        bool check_pt_location(sxyz<T1>& curr_pt,
                               const std::vector<T2>& ind1,
                               const std::array<T2,3>& ind2,
                               bool& onNode,
                               T2& nodeNo,
                               bool& onEdge,
                               std::array<T2,2>& edgeNodes,
                               bool& onFace,
                               std::array<T2,3>& faceNodes) const;
        /** @} */

        /**
         * @name Accumulating a tomography matrix row
         * Called once per raypath segment to spread the segment's length over
         * the nodes that influence the slowness along it.
         * @{
         */
        /**
         * @brief Add one segment's contribution to a matrix row.
         * @param[in,out] m_data the row being built.
         * @param m              scratch triplet, reused across calls.
         * @param allNodes       nodes whose slowness affects this segment.
         * @param mid_pt         segment midpoint, where the weights are
         *                       evaluated.
         * @param ds             segment length.
         */
        void update_m_data(std::vector<sijv<T1>>& m_data,
                           sijv<T1>& m,
                           const std::set<T2>& allNodes,
                           const sxyz<T1>& mid_pt,
                           const T1 ds) const ;

        /**
         * @brief Add one segment's contribution, scaled by squared slowness.
         * @param[in,out] m_data the row being built.
         * @param m              scratch triplet.
         * @param allNodes       nodes whose slowness affects this segment.
         * @param mid_pt         segment midpoint.
         * @param ds             segment length.
         * @param s_sq           squared slowness at the midpoint.
         *
         * The extra factor is what @ref g3dun_procvel demands: when velocity is
         * the interpolated quantity, the derivative with respect to a node value
         * picks up an @f$s^2@f$ from differentiating the reciprocal.
         */
        void update_m_data(std::vector<sijv<T1>>& m_data,
                           sijv<T1>& m,
                           const std::set<T2>& allNodes,
                           const sxyz<T1>& mid_pt,
                           const T1 ds,
                           const T1 s_sq) const ;
        /** @} */

        /**
         * @name Ray/element intersection
         * @{
         */
        /**
         * @brief Intersect a ray leaving a node with a triangle.
         * @param iO         node the ray leaves from.
         * @param vec        ray direction.
         * @param iA         first triangle node.
         * @param iB         second triangle node.
         * @param iC         third triangle node.
         * @param[out] pt_i  the intersection point.
         * @return true if the ray meets the triangle within its bounds.
         */
        bool intersectVecTriangle(const T2 iO, const sxyz<T1> &vec,
                                  const T2 iA, T2 iB, T2 iC,
                                  sxyz<T1> &pt_i) const;
        /**
         * @brief Intersect a ray leaving an arbitrary point with a triangle.
         * @param O          origin of the ray.
         * @param vec        ray direction.
         * @param iA         first triangle node.
         * @param iB         second triangle node.
         * @param iC         third triangle node.
         * @param[out] pt_i  the intersection point.
         * @return true if the ray meets the triangle within its bounds.
         */
        bool intersectVecTriangle(const sxyz<T1> &O, const sxyz<T1> &vec,
                                  const T2 iA, T2 iB, T2 iC,
                                  sxyz<T1> &pt_i) const;

        /**
         * @brief Intersect a ray with the edges of a face.
         * @param curr_pt        origin of the ray.
         * @param g              ray direction.
         * @param faceNodes      the face's three nodes.
         * @param[out] pt_i      the intersection point.
         * @param[out] edgeNodes the two nodes of the edge crossed.
         * @return true if an edge was crossed.
         *
         * Used when the ray leaves a face through one of its sides rather than
         * through its interior.
         */
        bool intersectVecEdge(const sxyz<T1>& curr_pt,
                              const sxyz<T1>& g,
                              std::array<T2,3>& faceNodes,
                              sxyz<T1>&  pt_i,
                              std::array<T2,2>& edgeNodes) const;
        /** @} */

        /**
         * @name Stepping to the next cell
         * @{
         */
        /**
         * @brief Find a cell sharing a face, given a node it must contain.
         * @param faceNodes the face's three nodes.
         * @param nodeNo    a node the cell must own.
         * @return the cell index, or the cell count if none qualifies.
         */
        T2 findAdjacentCell1(const std::array<T2,3> &faceNodes, const T2 nodeNo) const;
        /**
         * @brief Find the cell on the far side of a face.
         * @param faceNodes the face's three nodes.
         * @param cellNo    the cell being left.
         * @return the neighbouring cell, or the cell count if @p faceNodes lies
         *         on the mesh boundary.
         */
        T2 findAdjacentCell2(const std::array<T2,3> &faceNodes, const T2 cellNo) const;
        /**
         * @brief Find the cell on the far side of a face, disambiguating by
         *        point.
         * @param faceNodes the face's three nodes.
         * @param cellNo    the cell being left.
         * @param curr_pt   the current point, used to pick among candidates when
         *                  the face is degenerate or the point sits on an edge.
         * @return the neighbouring cell.
         */
        T2 findAdjacentCell2(const std::array<T2,3> &faceNodes,
                             const T2 & cellNo,
                             const sxyz<T1>& curr_pt ) const;
        /** @} */

        /**
         * @brief Collect the nodes sharing a cell with a given node.
         * @tparam SetT     any set-like container; the caller chooses, which
         *                  lets a polymorphic-memory-resource set be used to
         *                  keep the allocation off the heap in hot loops.
         * @param nodeNo    the node.
         * @param[out] nodes the neighbours found.
         */
        template<typename SetT> void getNeighborNodes(const T2 nodeNo, SetT& nodes) const;

        /**
         * @brief Group the nodes surrounding each of a set of nodes into
         *        triangles.
         * @param nodes     the nodes of interest.
         * @param[out] nb   for each input node, the triangles of neighbours
         *                  around it.
         *
         * Supplies ttcr::Grad3D_ab with the triangle fan it averages a gradient
         * over.
         */
        void getNeighborNodesAB(const std::vector<NODE*>& nodes,
                                std::vector<std::vector<std::array<NODE*,3>>>& nb) const;

        /**
         * @brief Collect enough nodes around one node to fit a Taylor expansion.
         * @param nodeNumber     the centre node.
         * @param minNbrPoints   minimum number of neighbours required; the
         *                       search widens through successive rings until it
         *                       is met.
         * @param[out] surroundingNodes the neighbourhood.
         *
         * @sa computeK, which needs 4 points for a first-order expansion and 10
         *     for a second-order one.
         */
        void getSurroundingNodes(const T2 nodeNumber,
                                 const T2 minNbrPoints,
                                 std::set<T2>& surroundingNodes) const;

        /**
         * @name Least-squares fit for the derivative operators
         * Build the design matrix of a Taylor expansion about one node, then
         * turn its pseudo-inverse into operator rows. The plain and @c 2
         * variants differ in whether the centre node's own value is an unknown
         * of the fit — see @c s0inside in @ref computeK.
         * @{
         */
        /**
         * @brief Build the design matrix, centre node value taken as known.
         * @param nodeNumber       the centre node.
         * @param surroundingNodes its neighbourhood.
         * @param weighting        weight rows by inverse distance.
         * @param order            order of the Taylor expansion, 1 or 2.
         * @param[out] A           the design matrix.
         * @param[out] W           the diagonal weight matrix.
         */
        void buildA(const T2 nodeNumber,
                    const std::set<T2>& surroundingNodes,
                    const bool weighting,
                    const int order,
                    Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>& A,
                    Eigen::Matrix<T1,Eigen::Dynamic, Eigen::Dynamic>& W) const;

        /**
         * @brief Build the design matrix, centre node value treated as unknown.
         * @param nodeNumber       the centre node.
         * @param surroundingNodes its neighbourhood.
         * @param weighting        weight rows by inverse distance.
         * @param order            order of the Taylor expansion, 1 or 2.
         * @param[out] A           the design matrix, one column wider than
         *                         @ref buildA gives.
         * @param[out] W           the diagonal weight matrix.
         */
        void buildA2(const T2 nodeNumber,
                     const std::set<T2>& surroundingNodes,
                     const bool weighting,
                     const int order,
                     Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>& A,
                     Eigen::Matrix<T1,Eigen::Dynamic, Eigen::Dynamic>& W) const;

        /**
         * @brief Copy fitted coefficients into the operator rows.
         * @param nodeNo           the centre node, indexing the row.
         * @param surroundingNodes its neighbourhood, indexing the columns.
         * @param i                row of @p Acoefs holding the x derivative.
         * @param j                row holding the y derivative.
         * @param k                row holding the z derivative.
         * @param Acoefs           the pseudo-inverse of the design matrix.
         * @param[out] k_data      the three operators being built.
         */
        void fill_k_data(const T2 nodeNo, const std::set<T2>& surroundingNodes,
                         const int i, const int j, const int k,
                         const Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>& Acoefs,
                         std::vector<std::vector<std::vector<siv<T1>>>>& k_data) const;

        /**
         * @brief Copy fitted coefficients into the operator rows, including the
         *        centre node's own weight.
         * @param nodeNo           the centre node.
         * @param surroundingNodes its neighbourhood.
         * @param i                row of @p Acoefs holding the x derivative.
         * @param j                row holding the y derivative.
         * @param k                row holding the z derivative.
         * @param Acoefs           the pseudo-inverse of the design matrix.
         * @param[out] k_data      the three operators being built.
         */
        void fill_k_data2(const T2 nodeNo, const std::set<T2>& surroundingNodes,
                          const int i, const int j, const int k,
                          const Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>& Acoefs,
                          std::vector<std::vector<std::vector<siv<T1>>>>& k_data) const;
        /** @} */

        /**
         * @brief Print a cell, a point within it and a gradient, for debugging.
         * @param cellNo the cell.
         * @param pt     the point.
         * @param g      the gradient.
         */
        void plotCell(const T2 cellNo, const sxyz<T1> &pt, const sxyz<T1> &g) const;

        /**
         * @name Sampling the slowness field
         * @{
         */
        /**
         * @brief Slowness at an arbitrary point.
         * @param pt           the point; taken by value as it may be translated
         *                     in place.
         * @param isTranslated set when @p pt is already in translated
         *                     coordinates.
         * @return the interpolated slowness.
         *
         * Locates the containing tetrahedron and interpolates over its four
         * primary nodes, honouring @c processVel.
         */
        T1 computeSlowness(sxyz<T1> pt, const bool isTranslated=false) const;
        /**
         * @brief Slowness at a point whose position is already classified.
         * @param curr_pt   the point.
         * @param onNode    set when it coincides with a node.
         * @param nodeNo    that node, when @p onNode.
         * @param onEdge    set when it lies on an edge.
         * @param edgeNodes that edge's two nodes.
         * @param faceNodes the face's three nodes, used when neither of the
         *                  above holds.
         * @return the slowness: read directly at a node, interpolated linearly
         *         along an edge, bilinearly across a face.
         *
         * The form the raypath tracer uses, since it always knows where the
         * point sits and so can skip the point location the other overload does.
         */
        T1 computeSlowness(const sxyz<T1>& curr_pt,
                           const bool onNode,
                           const T2 nodeNo,
                           const bool onEdge,
                           const std::array<T2,2>& edgeNodes,
                           const std::array<T2,3>& faceNodes) const;
        /** @} */

        /**
         * @name Filling in the secondary nodes
         * The model is defined at the primary nodes only, so after it is set the
         * secondary nodes must be given values consistent with it. Which of the
         * two a derived class calls follows from @c processVel.
         * @{
         */
        /**
         * @brief Interpolate slowness onto the secondary nodes.
         * @param nSecondary number of secondary nodes per edge.
         */
        void interpSlownessSecondary(const T2 nSecondary);
        /**
         * @brief Interpolate velocity onto the secondary nodes, storing the
         *        reciprocal.
         * @param nSecondary number of secondary nodes per edge.
         */
        void interpVelocitySecondary(const T2 nSecondary);
        /** @} */

        /**
         * @brief Print the tracer's state at one step, for debugging.
         * @param curr_pt   current point.
         * @param g         current gradient.
         * @param onNode    whether the point is on a node.
         * @param onEdge    whether it is on an edge.
         * @param onFace    whether it is on a face.
         * @param cellNo    the current cell.
         * @param nodeNo    the current node, where applicable.
         * @param edgeNodes the current edge, where applicable.
         * @param faceNodes the current face, where applicable.
         */
        void printRaypathData(const sxyz<T1>& curr_pt,
                              const sxyz<T1>& g,
                              const bool onNode,
                              const bool onEdge,
                              const bool onFace,
                              const T2 cellNo,
                              const T2 nodeNo,
                              const std::array<T2,2> &edgeNodes,
                              const std::array<T2,3> &faceNodes) const;

        /**
         * @brief The four primary nodes of a tetrahedron.
         * @param cellNo the cell.
         * @return its four vertices.
         *
         * @c neighbors[cellNo] lists secondary nodes alongside primary ones, so
         * the vertices have to be filtered out of it.
         */
        std::array<T2,4> getPrimary(const T2 cellNo) const {
            size_t i = 0;
            std::array<T2,4> tmp;
            for (size_t n=0; n<this->neighbors[cellNo].size(); ++n) {
                if ( nodes[this->neighbors[cellNo][n]].isPrimary() )
                    tmp[i++] = this->neighbors[cellNo][n];
                if ( i == 4 ) break;
            }
            return tmp;
        }

        /**
         * @name Steering into the source
         * Near the source the traveltime gradient is unreliable — the field is
         * steepest and most curved exactly there — so once the ray is within one
         * cell of a Tx the gradient is replaced outright by the direction
         * straight to it. Without this the ray tends to circle the source
         * without reaching it.
         * @{
         */
        /**
         * @brief Aim the gradient at a source, if the current cell touches one.
         * @param curr_pt   current point.
         * @param[in,out] g the gradient; overwritten with @c Tx-curr_pt on a hit.
         * @param cellNo    the cell the ray is in.
         * @param Tx        source coordinates.
         * @param txCell    the cell containing each Tx.
         *
         * A hit means @p cellNo is owned by one of the primary nodes of a Tx
         * cell, i.e. the two cells share at least a vertex.
         */
        void checkCloseToTx(const sxyz<T1>& curr_pt,
                            sxyz<T1>& g,
                            const T2 cellNo,
                            const std::vector<sxyz<T1>>& Tx,
                            const std::vector<T2>& txCell) const {

            for (size_t nt=0; nt<Tx.size(); ++nt) {
                std::array<T2,4> itmp = getPrimary(txCell[nt]);
                for (size_t n=0; n<4; ++n){
                    for (auto nc=nodes[itmp[n]].getOwners().begin(); nc!=nodes[itmp[n]].getOwners().end(); ++nc){
                        if ( *nc==cellNo ) {
                            g = Tx[nt]-curr_pt;
                            return;
                        }
                    }
                }
            }
        }

        /**
         * @brief Aim the gradient at a source, if the current edge touches one.
         * @param curr_pt   current point, known to lie on an edge.
         * @param[in,out] g the gradient; overwritten on a hit.
         * @param edgeNodes the edge's two nodes.
         * @param Tx        source coordinates.
         * @param txCell    the cell containing each Tx.
         *
         * A hit means an endpoint of the edge is a vertex of a Tx cell.
         */
        void checkCloseToTx(const sxyz<T1>& curr_pt,
                            sxyz<T1>& g,
                            const std::array<T2,2> &edgeNodes,
                            const std::vector<sxyz<T1>>& Tx,
                            const std::vector<T2>& txCell) const {

            for (size_t nt=0; nt<Tx.size(); ++nt) {
                std::array<T2,4> itmp = getPrimary(txCell[nt]);
                for (size_t n=0; n<4; ++n){
                    if ( itmp[n]==edgeNodes[0] || itmp[n]==edgeNodes[1] ) {
                        g = Tx[nt]-curr_pt;
                        return;
                    }
                }
            }
        }

        /** @} */

        /**
         * @name Projecting onto a face
         * When the ray sits on a face or an edge, the raw gradient generally
         * points out of the mesh or back into the cell just left. Projecting it
         * into the plane of a face keeps the ray on the surface it is
         * constrained to follow.
         * @{
         */
        /**
         * @brief Project a gradient into the plane of a face.
         * @param g         the gradient.
         * @param faceNodes the face's three nodes.
         * @return the component of @p g lying in the face.
         *
         * Computed as @f$\mathbf{n}\times(\mathbf{g}\times\mathbf{n})@f$ with
         * @f$\mathbf{n}@f$ the unit face normal.
         */
        sxyz<T1> projectOnFace(const sxyz<T1>& g, const std::array<T2,3>& faceNodes) const {

            // calculate normal to the face
            sxyz<T1> v1 = {nodes[faceNodes[0]].getX()-nodes[faceNodes[1]].getX(),
                nodes[faceNodes[0]].getY()-nodes[faceNodes[1]].getY(),
                nodes[faceNodes[0]].getZ()-nodes[faceNodes[1]].getZ()};
            sxyz<T1> v2 = {nodes[faceNodes[0]].getX()-nodes[faceNodes[2]].getX(),
                nodes[faceNodes[0]].getY()-nodes[faceNodes[2]].getY(),
                nodes[faceNodes[0]].getZ()-nodes[faceNodes[2]].getZ()};
            sxyz<T1> n = cross(v1, v2);
            n.normalize();
            return cross(n, cross(g, n));
        }

        /**
         * @brief Project a gradient onto the best face around an edge, and find
         *        where it leaves.
         * @param curr_pt            the point, on an edge.
         * @param g                  the gradient.
         * @param[in,out] edgeNodes  the edge; updated if the ray moves onto
         *                           another one.
         * @param cells              the cells sharing the edge, the candidates.
         * @param[out] pt_i          where the projected direction leaves the
         *                           chosen face.
         * @return the projected direction.
         *
         * An edge is shared by several faces, so unlike the face overload this
         * one must choose: it keeps the face whose projection stays closest to
         * @p g.
         */
        sxyz<T1> projectOnFace(const sxyz<T1>& curr_pt,
                               const sxyz<T1>& g,
                               std::array<T2,2>& edgeNodes,
                               const std::vector<T2> cells,
                               sxyz<T1>&  pt_i) const;

        /**
         * @brief Project a gradient onto a face around a node.
         * @param curr_pt           the point, on a node.
         * @param nodeNo            that node.
         * @param[in,out] g         the gradient; replaced by its projection.
         * @param[out] edgeNodes    the edge the ray moves onto.
         * @param[out] pt_i         where it leaves the face.
         * @return true if a suitable face was found.
         */
        bool projectOnFace(const sxyz<T1>& curr_pt,
                           const T2 nodeNo,
                           sxyz<T1>& g,
                           std::array<T2,2>& edgeNodes,
                           sxyz<T1>& pt_i) const;
        /**
         * @brief Project a gradient onto a face around a node, choosing by
         *        traveltime.
         * @param curr_pt           the point, on a node.
         * @param nodeNo            that node.
         * @param[in,out] g         the gradient; replaced by its projection.
         * @param[out] edgeNodes    the edge the ray moves onto.
         * @param[out] pt_i         where it leaves the face.
         * @param threadNo          thread number, needed to read traveltimes.
         * @return true if a suitable face was found.
         *
         * Having the traveltimes to hand lets this overload break ties by
         * preferring the face the wavefront actually arrived through.
         */
        bool projectOnFace(const sxyz<T1>& curr_pt,
                           const T2 nodeNo,
                           sxyz<T1>& g,
                           std::array<T2,2>& edgeNodes,
                           sxyz<T1>& pt_i,
                           const size_t & threadNo ) const;
        /**
         * @brief Project a gradient onto a face around an edge, choosing by
         *        traveltime.
         * @param curr_pt            the point, on an edge.
         * @param g                  the gradient.
         * @param[in,out] edgeNodes  the edge; updated if the ray moves on.
         * @param cells              the cells sharing the edge.
         * @param[out] pt_i          where the projected direction leaves the
         *                           face.
         * @param threadNo           thread number.
         * @return the projected direction.
         */
        sxyz<T1> projectOnFace(const sxyz<T1>& curr_pt,
                               const sxyz<T1>& g,
                               std::array<T2,2>& edgeNodes,
                               const std::vector<T2> cells,
                               sxyz<T1>&  pt_i,
                               const size_t & threadNo) const;
        /** @} */
    };

    template<typename T1, typename T2, typename NODE>
    T2 Grid3Dun<T1,T2,NODE>::getCellNo(const sxyz<T1>& pt) const {
        T2 closestNode = getNearestNode(pt);
        T1 minVolumeDiff = std::numeric_limits<T1>::max();
        T2 cell=0;

        for (auto tet=nodes[closestNode].getOwners().begin(); tet!=nodes[closestNode].getOwners().end(); ++tet) {
            T2 celli = *tet;
            for (size_t n=0;n<4;++n){
                T2 neighborNode = this->neighbors[celli][n];
                for (auto tet2=nodes[neighborNode].getOwners().begin(); tet2!=nodes[neighborNode].getOwners().end(); ++tet2) {
                    sxyz<T1> v1 = { nodes[ this->neighbors[*tet2][0] ]};
                    sxyz<T1> v2 = { nodes[ this->neighbors[*tet2][1] ]};
                    sxyz<T1> v3 = { nodes[ this->neighbors[*tet2][2] ]};
                    sxyz<T1> v4 = { nodes[ this->neighbors[*tet2][3] ]};

                    T1 D0 = 1.e6*det4(v1, v2, v3, v4);
                    T1 D1 = 1.e6*det4(pt, v2, v3, v4);
                    T1 D2 = 1.e6*det4(v1, pt, v3, v4);
                    T1 D3 = 1.e6*det4(v1, v2, pt, v4);
                    T1 D4 = 1.e6*det4(v1, v2, v3, pt);

                    T1 VolumeDiff = std::abs(std::abs(D0)-std::abs(D1)-std::abs(D2)-std::abs(D3)-std::abs(D4));
                    if (VolumeDiff < minVolumeDiff) {
                        minVolumeDiff = VolumeDiff;
                        cell = *tet2;
                    }
                }
            }
        }
        return cell;
    }

    template<typename T1, typename T2, typename NODE>
    const typename Grid3Dun<T1,T2,NODE>::txInfo_t&
    Grid3Dun<T1,T2,NODE>::getTxInfo(const std::vector<sxyz<T1>>& Tx,
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
        ti.txOnEdge.assign( Tx.size(), false );
        ti.txOnFace.assign( Tx.size(), false );
        ti.txNode.assign( Tx.size(), 0 );
        ti.txCell.assign( Tx.size(), 0 );
        ti.txEdges.assign( Tx.size(), std::array<T2,2>{{0, 0}} );
        ti.txFaces.assign( Tx.size(), std::array<T2,3>{{0, 0, 0}} );
        ti.txNeighborCells.assign( Tx.size(), std::vector<T2>() );

        for ( size_t nt=0; nt<Tx.size(); ++nt ) {
            for ( T2 nn=0; nn<nodes.size(); ++nn ) {
                if ( nodes[nn].isPrimary() ) {
                    if ( nodes[nn] == Tx[nt] ) {
                        ti.txOnNode[nt] = true;
                        ti.txNode[nt] = nn;
                        break;
                    }
                }
            }
        }
        for ( size_t nt=0; nt<Tx.size(); ++nt ) {
            if ( !ti.txOnNode[nt] ) {
                ti.txCell[nt] = getCellNo( Tx[nt] );

                std::array<T2,4> itmp = getPrimary(ti.txCell[nt]);
                // find adjacent cells
                const T2 ind[6][2] = {
                    {itmp[0], itmp[1]},
                    {itmp[0], itmp[2]},
                    {itmp[0], itmp[3]},
                    {itmp[1], itmp[2]},
                    {itmp[1], itmp[3]},
                    {itmp[2], itmp[3]} };

                for ( size_t nedge=0; nedge<6; ++nedge ) {
                    for ( auto nc0=nodes[ind[nedge][0]].getOwners().begin(); nc0!=nodes[ind[nedge][0]].getOwners().end(); ++nc0 ) {
                        if ( std::find(nodes[ind[nedge][1]].getOwners().begin(), nodes[ind[nedge][1]].getOwners().end(), *nc0)!=nodes[ind[nedge][1]].getOwners().end() )
                            ti.txNeighborCells[nt].push_back( *nc0 );
                    }
                }
                // check if on edge
                for ( size_t nedge=0; nedge<6; ++nedge ) {
                    if ( distSqPointToSegment( &nodes[ind[nedge][0]], &nodes[ind[nedge][1]], Tx[nt]) < small2 ) {
                        ti.txOnEdge[nt] = true;
                        ti.txEdges[nt][0] = ind[nedge][0];
                        ti.txEdges[nt][1] = ind[nedge][1];
                        break;
                    }
                }
                if ( !ti.txOnEdge[nt] ) {
                    // check if on face
                    const T2 indf[4][3] = {
                        {itmp[0], itmp[1], itmp[2]},
                        {itmp[0], itmp[1], itmp[3]},
                        {itmp[0], itmp[2], itmp[3]},
                        {itmp[1], itmp[2], itmp[3]} };
                    for ( size_t nface=0; nface<4; ++nface ) {
                        if ( testInTriangle(&nodes[indf[nface][0]], &nodes[indf[nface][1]], &nodes[indf[nface][2]], Tx[nt]) ) {
                            ti.txOnFace[nt] = true;
                            ti.txFaces[nt][0] = indf[nface][0];
                            ti.txFaces[nt][1] = indf[nface][1];
                            ti.txFaces[nt][2] = indf[nface][2];
                            break;
                        }
                    }
                }
            }
        }
        ti.valid = true;
        return ti;
    }

    template<typename T1, typename T2, typename NODE>
    const typename Grid3Dun<T1,T2,NODE>::txInfo_t&
    Grid3Dun<T1,T2,NODE>::getTxInfoBlti(const std::vector<sxyz<T1>>& Tx,
                                        const size_t threadNo) const {

        if ( txInfoCacheBlti.size() <= threadNo ) {
            txInfoCacheBlti.resize( threadNo + 1 );
        }
        txInfo_t& ti = txInfoCacheBlti[threadNo];
        if ( ti.valid && ti.tx == Tx ) {
            return ti;
        }

        ti.tx = Tx;
        ti.txOnNode.assign( Tx.size(), false );
        ti.txNode.assign( Tx.size(), 0 );
        ti.txCell.assign( Tx.size(), 0 );
        ti.txNeighborCells.assign( Tx.size(), std::vector<T2>() );

        for ( size_t nt=0; nt<Tx.size(); ++nt ) {
            for ( T2 nn=0; nn<nodes.size(); ++nn ) {
                if ( nodes[nn] == Tx[nt] ) {
                    ti.txOnNode[nt] = true;
                    ti.txNode[nt] = nn;
                    break;
                }
            }
        }
        for ( size_t nt=0; nt<Tx.size(); ++nt ) {
            if ( !ti.txOnNode[nt] ) {
                ti.txCell[nt] = getCellNo( Tx[nt] );

                std::array<T2,4> itmp = getPrimary(ti.txCell[nt]);
                // find adjacent cells
                const T2 ind[6][2] = {
                    {itmp[0], itmp[1]},
                    {itmp[0], itmp[2]},
                    {itmp[0], itmp[3]},
                    {itmp[1], itmp[2]},
                    {itmp[1], itmp[3]},
                    {itmp[2], itmp[3]} };

                for ( size_t nedge=0; nedge<6; ++nedge ) {
                    for ( auto nc0=nodes[ind[nedge][0]].getOwners().begin(); nc0!=nodes[ind[nedge][0]].getOwners().end(); ++nc0 ) {
                        if ( std::find(nodes[ind[nedge][1]].getOwners().begin(), nodes[ind[nedge][1]].getOwners().end(), *nc0)!=nodes[ind[nedge][1]].getOwners().end() )
                            ti.txNeighborCells[nt].push_back( *nc0 );
                    }
                }
            }
        }
        ti.valid = true;
        return ti;
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::buildGridNodes(const std::vector<sxyz<T1>>& no,
                                              const size_t nt) {

        if ( this->translateOrigin ) {
            T1 xmin = no[0].x;
            T1 ymin = no[0].y;
            T1 zmin = no[0].z;
            for ( T2 n=1; n<no.size(); ++n ) {
                xmin = xmin < no[n].x ? xmin : no[n].x;
                ymin = ymin < no[n].y ? ymin : no[n].y;
                zmin = zmin < no[n].z ? zmin : no[n].z;
            }
            this->origin = {xmin, ymin, zmin};
        } else {
            this->origin = {0.0, 0.0, 0.0};
        }

        // primary nodes
        for ( T2 n=0; n<no.size(); ++n ) {
            nodes[n].setXYZindex(no[n].x - this->origin.x,
                                 no[n].y - this->origin.y,
                                 no[n].z - this->origin.z, n );
            nodes[n].setPrimary(true);
        }

        //
        //              1
        //            ,/|`\
        //          ,/  |  `\
        //        ,0    '.   `4
        //      ,/       1     `\
        //    ,/         |       `\
        //   0-----5-----'.--------3
        //    `\.         |      ,/
        //       `\.      |     3
        //          `2.   '. ,/
        //             `\. |/
        //                `2
        //
        //
        //  triangle 0:  0-1  1-2  2-0     (first occurence of segment underlined)
        //               ---  ---  ---
        //  triangle 1:  1-2  2-3  3-1
        //                    ---  ---
        //  triangle 2:  0-2  2-3  3-0
        //                         ---
        //  triangle 3:  0-1  1-3  3-0


        for ( T2 ntet=0; ntet<tetrahedra.size(); ++ntet ) {

            // for each triangle
            for ( T2 ntri=0; ntri<4; ++ntri ) {

                // push owner for primary nodes
                nodes[ tetrahedra[ntet].i[ntri] ].pushOwner( ntet );
            }
        }
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::buildGridNodes(const std::vector<sxyz<T1>>& no,
                                              const T2 nsecondary,
                                              const size_t nt) {

        if ( this->translateOrigin ) {
            T1 xmin = no[0].x;
            T1 ymin = no[0].y;
            T1 zmin = no[0].z;
            for ( T2 n=1; n<no.size(); ++n ) {
                xmin = xmin < no[n].x ? xmin : no[n].x;
                ymin = ymin < no[n].y ? ymin : no[n].y;
                zmin = zmin < no[n].z ? zmin : no[n].z;
            }
            this->origin = {xmin, ymin, zmin};
        } else {
            this->origin = {0.0, 0.0, 0.0};
        }

        // primary nodes
        for ( size_t n=0; n<no.size(); ++n ) {
            nodes[n].setXYZindex( no[n].x, no[n].y, no[n].z, static_cast<T2>(n) );
            nodes[n].setPrimary(true);
        }
        T2 nNodes = static_cast<T2>(nodes.size());

        size_t nFaceNodes = 0;
        for ( size_t n=1; n<=(nsecondary-1); ++n ) {
            nFaceNodes += n;
        }

        std::map<std::array<T2,2>,std::vector<T2>> lineMap;
        std::array<T2,2> lineKey;
        typename std::map<std::array<T2,2>,std::vector<T2>>::iterator lineIt;

        size_t estLineNo = (tetrahedra.size()+tetrahedra.size()/10) * 6/2;
        size_t estFaceNo = (tetrahedra.size()+tetrahedra.size()/10) * 4/2;
        nodes.reserve( nNodes + estLineNo*nsecondary + estFaceNo*nFaceNodes );

        T2 iNodes[4][3] = {
            {0,1,2},  // (relative) indices of nodes of 1st triangle
            {1,2,3},  // (relative) indices of nodes of 2nd triangle
            {0,2,3},  // (relative) indices of nodes of 3rd triangle
            {0,1,3}   // (relative) indices of nodes of 4th triangle
        };

        //
        //              1
        //            ,/|`\
        //          ,/  |  `\
        //        ,0    '.   `4
        //      ,/       1     `\
        //    ,/         |       `\
        //   0-----5-----'.--------3
        //    `\.         |      ,/
        //       `\.      |     3
        //          `2.   '. ,/
        //             `\. |/
        //                `2
        //
        //
        //  triangle 0:  0-1  1-2  2-0     (first occurence of segment underlined)
        //               ---  ---  ---
        //  triangle 1:  1-2  2-3  3-1
        //                    ---  ---
        //  triangle 2:  0-2  2-3  3-0
        //                         ---
        //  triangle 3:  0-1  1-3  3-0

        if ( verbose>1 && nsecondary > 0 ) {
            std::cout << '\n';
        }

        // edge nodes
        NODE tmpNode(nt);
        for ( T2 ntet=0; ntet<tetrahedra.size(); ++ntet ) {

            if ( verbose>1 && nsecondary > 0 ) {
                std::cout << "\r  Building edge nodes: " << (100*ntet)/tetrahedra.size() << "%";
                std::cout.flush();
            }

            // for each triangle
            for ( T2 ntri=0; ntri<4; ++ntri ) {

                // push owner for primary nodes
                nodes[ tetrahedra[ntet].i[ntri] ].pushOwner( ntet );

                if ( nsecondary > 0 ) {
                    // start from ntri to avoid redundancy
                    for ( size_t nl=ntri; nl<3; ++nl ) {

                        lineKey = {tetrahedra[ntet].i[ iNodes[ntri][nl] ],
                            tetrahedra[ntet].i[ iNodes[ntri][(nl+1)%3] ]};
                        std::sort(lineKey.begin(), lineKey.end());

                        lineIt = lineMap.find( lineKey );
                        if ( lineIt == lineMap.end() ) {
                            // not found, insert new pair
                            lineMap[ lineKey ] = std::vector<T2>(nsecondary);
                        } else {
                            for ( size_t n=0; n<lineIt->second.size(); ++n ) {
                                // setting owners
                                nodes[ lineIt->second[n] ].pushOwner( ntet );
                            }
                            continue;
                        }

                        sxyz<T1> d = (no[lineKey[1]]-no[lineKey[0]])/static_cast<T1>(nsecondary+1);

                        for ( size_t n2=0; n2<nsecondary; ++n2 ) {
                            tmpNode.setXYZindex(no[lineKey[0]].x+(1+n2)*d.x,
                                                no[lineKey[0]].y+(1+n2)*d.y,
                                                no[lineKey[0]].z+(1+n2)*d.z,
                                                nNodes );
                            lineMap[lineKey][n2] = nNodes++;
                            nodes.push_back( tmpNode );
                            nodes.back().pushOwner( ntet );
                        }
                    }
                }
            }
        }

        if ( verbose>1 && nsecondary > 0 ) {
            std::cout << "\r  Building edge nodes: 100%\n";
        }

        if ( nsecondary > 1 ) {

            std::map<std::array<T2,3>,std::vector<T2>> faceMap;
            std::array<T2,3> faceKey;
            typename std::map<std::array<T2,3>,std::vector<T2>>::iterator faceIt;

            ptrdiff_t ncut = nsecondary - 1;

            for ( T2 ntet=0; ntet<tetrahedra.size(); ++ntet ) {

                if ( verbose>1 ) {
                    std::cout << "\r  Building face nodes: " << (100*ntet)/tetrahedra.size() << "%";
                    std::cout.flush();
                }

                // for each triangle
                for ( T2 ntri=0; ntri<4; ++ntri ) {

                    faceKey = {tetrahedra[ntet].i[ iNodes[ntri][0] ],
                        tetrahedra[ntet].i[ iNodes[ntri][1] ],
                        tetrahedra[ntet].i[ iNodes[ntri][2] ]};
                    std::sort(faceKey.begin(), faceKey.end());

                    faceIt = faceMap.find( faceKey );
                    if ( faceIt == faceMap.end() ) {
                        // not found, insert new pair
                        faceMap[ faceKey ] = std::vector<T2>(nFaceNodes);
                    } else {
                        for ( size_t n=0; n<faceIt->second.size(); ++n ) {
                            // setting owners
                            nodes[ faceIt->second[n] ].pushOwner( ntet );
                        }
                        continue;
                    }

                    sxyz<T1> d1 = (no[faceKey[1]]-no[faceKey[0]])/static_cast<T1>(nsecondary+1);
                    sxyz<T1> d2 = (no[faceKey[1]]-no[faceKey[2]])/static_cast<T1>(nsecondary+1);

                    size_t ifn = 0;
                    for ( ptrdiff_t n=0; n<ncut; ++n ) {

                        sxyz<T1> pt1 = no[faceKey[0]]+static_cast<T1>(1+n)*d1;
                        sxyz<T1> pt2 = no[faceKey[2]]+static_cast<T1>(1+n)*d2;

                        size_t nseg = ncut+1-n;

                        sxyz<T1> d = (pt2-pt1)/static_cast<T1>(nseg);

                        for ( size_t n2=0; n2<nseg-1; ++n2 ) {
                            tmpNode.setXYZindex(pt1.x+(1+n2)*d.x,
                                                pt1.y+(1+n2)*d.y,
                                                pt1.z+(1+n2)*d.z,
                                                nNodes );
                            faceMap[faceKey][ifn++] = nNodes++;
                            nodes.push_back( tmpNode );
                            nodes.back().pushOwner( ntet );
                        }
                    }
                }
            }
        }
        if ( verbose>1 && nsecondary > 0 ) {
            std::cout << "\r  Building face nodes: 100%\n";
        }

        nodes.shrink_to_fit();

        if ( this->translateOrigin ) {
            for (auto node=nodes.begin(); node!=nodes.end(); ++node) {
                *node -= this->origin;
            }
        }
    }

    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dun<T1,T2,NODE>::getTraveltime(const sxyz<T1>& Rx,
                                           const size_t threadNo) const {

        for ( size_t nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == Rx ) {
                return nodes[nn].getTT(threadNo);
            }
        }
        // TODO : proceed by interpolation
        // If Rx is not on a node
        T1 slo = computeSlowness( Rx, true );

        T2 cellNo = getCellNo( Rx );

        T2 neibNo = this->neighbors[cellNo][0];
        T1 dt = computeDt(nodes[neibNo], Rx, slo);

        T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
        for ( size_t k=1; k< this->neighbors[cellNo].size(); ++k ) {
            neibNo = this->neighbors[cellNo][k];
            dt = computeDt(nodes[neibNo], Rx, slo);
            if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
                traveltime =  nodes[neibNo].getTT(threadNo)+dt;
            }
        }
        return traveltime;
    }


    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dun<T1,T2,NODE>::getTraveltime(const sxyz<T1>& Rx,
                                           const std::vector<NODE>& nodes,
                                           const size_t threadNo) const {

        for ( size_t nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == Rx ) {
                return nodes[nn].getTT(threadNo);
            }
        }
        //If Rx is not on a node:
        T1 slo = computeSlowness( Rx, true );

        T2 cellNo = getCellNo( Rx );

        T2 neibNo = this->neighbors[cellNo][0];
        T1 dt = computeDt(nodes[neibNo], Rx, slo);

        T1 traveltime = nodes[neibNo].getTT(threadNo)+dt;
        for ( size_t k=1; k< this->neighbors[cellNo].size(); ++k ) {
            neibNo = this->neighbors[cellNo][k];
            dt = computeDt(nodes[neibNo], Rx, slo);
            if ( traveltime > nodes[neibNo].getTT(threadNo)+dt ) {
                traveltime =  nodes[neibNo].getTT(threadNo)+dt;
            }
        }
        return traveltime;
    }


    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::checkPts(std::vector<sxyz<T1>> pts, const bool translated) const {
        
        if (this->translateOrigin == true && translated == false) {
            for ( size_t n=0; n<pts.size(); ++n ) {
                pts[n] -= this->origin;
            }
        }

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
                // check if inside tetrahedra
                for ( T2 nt=0; nt<tetrahedra.size(); ++nt ) {
                    if ( insideTetrahedron(pts[n], nt) ) {
                        found = true;
                        break;
                    }
                }
            }
            if ( found == false ) {
                std::ostringstream msg;
                msg << "Error: Point (" << pts[n] << ") outside grid.";
                throw std::runtime_error(msg.str());
            }
        }
    }

    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::insideTetrahedron2(const sxyz<T1>& p, const T2 nt) const {

        sxyz<T1> v1 = { nodes[ tetrahedra[nt].i[0] ].getX(),
            nodes[ tetrahedra[nt].i[0] ].getY(),
            nodes[ tetrahedra[nt].i[0] ].getZ() };

        sxyz<T1> v2 = { nodes[ tetrahedra[nt].i[1] ].getX(),
            nodes[ tetrahedra[nt].i[1] ].getY(),
            nodes[ tetrahedra[nt].i[1] ].getZ() };

        sxyz<T1> v3 = { nodes[ tetrahedra[nt].i[2] ].getX(),
            nodes[ tetrahedra[nt].i[2] ].getY(),
            nodes[ tetrahedra[nt].i[2] ].getZ() };

        sxyz<T1> v4 = { nodes[ tetrahedra[nt].i[3] ].getX(),
            nodes[ tetrahedra[nt].i[3] ].getY(),
            nodes[ tetrahedra[nt].i[3] ].getZ() };

        return sameSide(v1, v2, v3, v4, p) && sameSide(v2, v3, v4, v1, p) &&
        sameSide(v3, v4, v1, v2, p) && sameSide(v4, v1, v2, v3, p);
    }

    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::insideTetrahedron(const sxyz<T1>& v, const T2 nt) const {


        // from http://steve.hollasch.net/cgindex/geometry/ptintet.html

        sxyz<T1> v1 = { nodes[ tetrahedra[nt].i[0] ].getX(),
            nodes[ tetrahedra[nt].i[0] ].getY(),
            nodes[ tetrahedra[nt].i[0] ].getZ() };

        sxyz<T1> v2 = { nodes[ tetrahedra[nt].i[1] ].getX(),
            nodes[ tetrahedra[nt].i[1] ].getY(),
            nodes[ tetrahedra[nt].i[1] ].getZ() };

        sxyz<T1> v3 = { nodes[ tetrahedra[nt].i[2] ].getX(),
            nodes[ tetrahedra[nt].i[2] ].getY(),
            nodes[ tetrahedra[nt].i[2] ].getZ() };

        sxyz<T1> v4 = { nodes[ tetrahedra[nt].i[3] ].getX(),
            nodes[ tetrahedra[nt].i[3] ].getY(),
            nodes[ tetrahedra[nt].i[3] ].getZ() };

        T1 D0 = det4(v1, v2, v3, v4);
        T1 D1 = det4( v, v2, v3, v4);
        T1 D2 = det4(v1,  v, v3, v4);
        T1 D3 = det4(v1, v2,  v, v4);
        T1 D4 = det4(v1, v2, v3,  v);

        int t1 = (signum(D0)==signum(D1));
        int t2 = (signum(D0)==signum(D2));
        int t3 = (signum(D0)==signum(D3));
        int t4 = (signum(D0)==signum(D4));

        bool it1, it2, it3, it4;
        it1 = it2 = it3 = it4 = 0;

        if ( std::abs(D1)<small ) {
            // points are coplanar, check if pt is inside triangle
            it1 = testInTriangle(v2, v3, v4, v);
        }
        if ( std::abs(D2)<small2 ) {
            it2 = testInTriangle(v1, v3, v4, v);
        }

        if ( std::abs(D3)<small2 ) {
            it3 = testInTriangle(v1, v2, v4, v);
        }

        if ( std::abs(D4)<small2 ) {
            it4 = testInTriangle(v1, v2, v3, v);
        }

        return ( t1 && t2 && t3 && t4 ) || it1 || it2 || it3 || it4;
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::saveTT(const std::string &fname, const int all,
                                      const size_t nt, const int format) const {

        if ( format == 1 ) {
            std::string filename = fname+".dat";
            std::ofstream fout(filename.c_str());
            T2 nMax = nPrimary;
            if ( all == 1 ) {
                nMax = static_cast<T2>(nodes.size());
            }
            for ( T2 n=0; n<nMax; ++n ) {
                fout << nodes[n].getX() << '\t'
                << nodes[n].getY() << '\t'
                << nodes[n].getZ() << '\t'
                << nodes[n].getTT(nt) << '\n';
            }
            fout.close();
        } else if ( format == 2 ) {
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

            vtkSmartPointer<vtkTetra> tet =
            vtkSmartPointer<vtkTetra>::New();
            for (size_t n=0; n<tetrahedra.size(); ++n) {
                tet->GetPointIds()->SetId(0, tetrahedra[n].i[0] );
                tet->GetPointIds()->SetId(1, tetrahedra[n].i[1] );
                tet->GetPointIds()->SetId(2, tetrahedra[n].i[2] );
                tet->GetPointIds()->SetId(3, tetrahedra[n].i[3] );

                ugrid->InsertNextCell( tet->GetCellType(), tet->GetPointIds() );
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
                T1 tmp[] = { nodes[n].getX(), nodes[n].getY(), nodes[n].getZ(), nodes[n].getTT(nt) };
                fout.write( (char*)tmp, 4*sizeof(T1) );
            }
            fout.close();
        } else {
            throw std::runtime_error("Unsupported format for saving traveltimes");
        }
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::loadTT(const std::string &fname, const int all,
                                      const size_t nt, const int format) const {

        if ( format == 1 ) {
            std::string filename = fname+".dat";
            std::ifstream fin(filename.c_str());
            T2 nMax = nPrimary;
            if ( all == 1 ) {
                nMax = static_cast<T2>(nodes.size());
            }
            T1 x, y, z, tt;
            for ( T2 n=0; n<nMax; ++n ) {
                fin >> x >> y >> z >> tt;
                nodes[n].setTT(tt, nt);
            }
            fin.close();
        } else if ( format == 2 ) {
#ifdef VTK
            std::string filename = fname+".vtu";

            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

            reader->SetFileName( filename.c_str() );
            reader->Update();

            T2 nMax = nPrimary;  // only primary were saved
            for (size_t n=0; n<nMax; ++n) {
                nodes[n].setTT(reader->GetOutput()->GetPointData()->GetArray("Travel time")->GetTuple1(n), nt);
            }
#else
            std::cerr << "VTK not included during compilation.\n";
#endif
        } else if ( format == 3 ) {
            std::string filename = fname+".bin";
            std::ifstream fin(filename.c_str(), std::ios::binary);
            T2 nMax = nPrimary;
            if ( all == 1 ) {
                nMax = static_cast<T2>(nodes.size());
            }
            for ( T2 n=0; n<nMax; ++n ) {
                T1 tmp[4];
                fin.read( (char*)tmp, 4*sizeof(T1) );
                nodes[n].setTT(tmp[3], nt);
            }
            fin.close();
        } else {
            throw std::runtime_error("Unsupported format for traveltimes");
        }
    }

#ifdef VTK
    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::saveModelVTU(const std::string &fname,
                                            const bool saveSlowness,
                                            const bool savePhysicalEntity) const {

        vtkSmartPointer<vtkUnstructuredGrid> ugrid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();

        vtkSmartPointer<vtkPoints> newPts =
        vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkDoubleArray> newScalars =
        vtkSmartPointer<vtkDoubleArray>::New();

        double xyz[3];
        T2 nMax = nPrimary;  // only primary are saved
        if ( saveSlowness ) {
            newScalars->SetName("Slowness");

            for (size_t n=0; n<nMax; ++n) {
                xyz[0] = nodes[n].getX();
                xyz[1] = nodes[n].getY();
                xyz[2] = nodes[n].getZ();
                newPts->InsertPoint(n, xyz);
                newScalars->InsertValue(n, nodes[n].getNodeSlowness() );
            }
        } else {
            newScalars->SetName("Velocity");

            for (size_t n=0; n<nMax; ++n) {
                xyz[0] = nodes[n].getX();
                xyz[1] = nodes[n].getY();
                xyz[2] = nodes[n].getZ();
                newPts->InsertPoint(n, xyz);
                newScalars->InsertValue(n, static_cast<T1>(1.0)/nodes[n].getNodeSlowness() );
            }
        }

        ugrid->SetPoints(newPts);
        ugrid->GetPointData()->SetScalars(newScalars);

        vtkSmartPointer<vtkTetra> tet =
        vtkSmartPointer<vtkTetra>::New();

        for (size_t n=0; n<tetrahedra.size(); ++n) {
            tet->GetPointIds()->SetId(0, tetrahedra[n].i[0] );
            tet->GetPointIds()->SetId(1, tetrahedra[n].i[1] );
            tet->GetPointIds()->SetId(2, tetrahedra[n].i[2] );
            tet->GetPointIds()->SetId(3, tetrahedra[n].i[3] );

            ugrid->InsertNextCell( tet->GetCellType(), tet->GetPointIds() );
        }

        vtkSmartPointer<vtkIntArray> data_pe = vtkSmartPointer<vtkIntArray>::New();
        if ( savePhysicalEntity ) {
            data_pe->SetName("Physical entity");
            for (size_t n=0; n<tetrahedra.size(); ++n) {
                data_pe->InsertNextValue(tetrahedra[n].physical_entity );
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

#endif



    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::localUpdate3D(NODE *vertexD,
                                             const size_t threadNo) const {

        // method of Lelievre et al. 2011

        T2 iA, iB, iC, iD;
        NODE *vertexA, *vertexB, *vertexC;

        for ( size_t no=0; no<vertexD->getOwners().size(); ++no ) {

            T2 tetNo = vertexD->getOwners()[no];

            for ( iD=0; iD<4; ++iD ) {
                if ( vertexD->getGridIndex() == tetrahedra[tetNo].i[iD] ) break;
            }

            iA = (iD+1)%4;
            iB = (iD+2)%4;
            iC = (iD+3)%4;
            vertexA = &(nodes[tetrahedra[tetNo].i[iA]]);
            vertexB = &(nodes[tetrahedra[tetNo].i[iB]]);
            vertexC = &(nodes[tetrahedra[tetNo].i[iC]]);

            if ( vertexA->getTT(threadNo) > vertexB->getTT(threadNo) ) {
                std::swap(iA, iB);
                std::swap(vertexA, vertexB);
            }
            if ( vertexA->getTT(threadNo) > vertexC->getTT(threadNo) ) {
                std::swap(iA, iC);
                std::swap(vertexA, vertexC);
            }
            if ( vertexB->getTT(threadNo) > vertexC->getTT(threadNo) ) {
                std::swap(iB, iC);
                std::swap(vertexB, vertexC);
            }

            if ( vertexA->getTT(threadNo) == std::numeric_limits<T1>::max() ) {
                continue;
            }

            T1 tABC = std::numeric_limits<T1>::max();

            if ( vertexB->getTT(threadNo) != std::numeric_limits<T1>::max() &&
                vertexC->getTT(threadNo) != std::numeric_limits<T1>::max() ) {

                T1 u = vertexB->getTT(threadNo) - vertexA->getTT(threadNo);
                T1 v = vertexC->getTT(threadNo) - vertexA->getTT(threadNo);

                sxyz<T1> v_b = { vertexC->getX() - vertexA->getX(),
                    vertexC->getY() - vertexA->getY(),
                    vertexC->getZ() - vertexA->getZ() };
                sxyz<T1> v_c = { vertexB->getX() - vertexA->getX(),
                    vertexB->getY() - vertexA->getY(),
                    vertexB->getZ() - vertexA->getZ() };

                // vector normal to plane
                sxyz<T1> v_n = cross(v_b, v_c);

                T1 b = norm( v_b );
                T1 c = norm( v_c );
                T1 d2 = dot(v_b, v_c);

                T1 alpha = acos( d2 / (b*c) );

                T1 phi = c*b*sin(alpha);

                // check for negative value
                T1 w_tilde = vertexD->getNodeSlowness()*vertexD->getNodeSlowness()*phi*phi -
                u*u*b*b - v*v*c*c + 2.*u*v*d2;
                if ( w_tilde > 0.0 ) {

                    w_tilde = sqrt( w_tilde );

                    // Point (ξ_0 , ζ_0 ) is the normalized projection of node D onto face ABC
                    // project D on plane

                    T1 d_tmp = -vertexA->getX()*v_n.x - vertexA->getY()*v_n.y - vertexA->getZ()*v_n.z;

                    T1 k = -(d_tmp + v_n.x*vertexD->getX() + v_n.y*vertexD->getY() + v_n.z*vertexD->getZ())/
                    norm2(v_n);

                    sxyz<T1> pt;   // -> Point (ξ_0 , ζ_0 )
                    pt.x = vertexD->getX() + k*v_n.x;
                    pt.y = vertexD->getY() + k*v_n.y;
                    pt.z = vertexD->getZ() + k*v_n.z;

                    T1 rho0 = vertexD->getDistance( pt );

                    sxyz<T1> v_pt = {pt.x-vertexA->getX(), pt.y-vertexA->getY(), pt.z-vertexA->getZ()};

                    T1 xi0;
                    T1 zeta0;
                    projNorm(v_b/b, v_c/c, v_pt, xi0, zeta0);
                    if ( xi0 < 0.0 || zeta0 < 0.0 ) {
                        // this should not happen unless we have incorrect triangle
                        continue;
                    }

                    T1 beta = u*b*b - v*d2;
                    T1 gamma = v*c*c - u*d2;

                    T1 xi_tilde = -std::abs(beta)*rho0/(phi*w_tilde);
                    T1 zeta_tilde = -std::abs(gamma)*rho0/(phi*w_tilde);

                    T1 xi = xi_tilde + xi0;
                    T1 zeta = zeta_tilde + zeta0;

                    if ( 0.<xi && xi<1. && 0.<zeta && zeta<1. && 0.<(xi+zeta) && (xi+zeta)<1. ) {
                        tABC = vertexA->getTT(threadNo) + u*xi0 + v*zeta0 + w_tilde*rho0/phi;
                    }
                }
            }

            T1 t = vertexA->getTT(threadNo) + vertexD->getNodeSlowness() * vertexD->getDistance( *vertexA );
            if ( t < tABC ) tABC = t;
            t = vertexB->getTT(threadNo) + vertexD->getNodeSlowness() * vertexD->getDistance( *vertexB );
            if ( t < tABC ) tABC = t;
            t = vertexC->getTT(threadNo) + vertexD->getNodeSlowness() * vertexD->getDistance( *vertexC );
            if ( t < tABC ) tABC = t;

            t = localUpdate2D(vertexA, vertexB, vertexD, tetNo, threadNo);
            if ( t < tABC ) tABC = t;
            t = localUpdate2D(vertexA, vertexC, vertexD, tetNo, threadNo);
            if ( t < tABC ) tABC = t;
            t = localUpdate2D(vertexB, vertexC, vertexD, tetNo, threadNo);
            if ( t < tABC ) tABC = t;

            if ( tABC<vertexD->getTT(threadNo) )
                vertexD->setTT(tABC, threadNo);

        }
    }


    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dun<T1,T2,NODE>::localUpdate2D(const NODE *vertexA,
                                           const NODE *vertexB,
                                           const NODE *vertexC,
                                           const T2 tetNo,
                                           const size_t threadNo) const {

        if ( vertexB->getTT(threadNo)==std::numeric_limits<T1>::max() &&
            vertexA->getTT(threadNo)==std::numeric_limits<T1>::max() ) {
            return std::numeric_limits<T1>::max();
        }
        T1 t;

        T1 u = vertexB->getTT(threadNo) - vertexA->getTT(threadNo);

        sxyz<T1> v_b = { vertexC->getX() - vertexA->getX(),
            vertexC->getY() - vertexA->getY(),
            vertexC->getZ() - vertexA->getZ() };
        sxyz<T1> v_c = { vertexB->getX() - vertexA->getX(),
            vertexB->getY() - vertexA->getY(),
            vertexB->getZ() - vertexA->getZ() };

        T1 c = norm( v_c );

        T1 w2 = vertexC->getNodeSlowness()*vertexC->getNodeSlowness()*c*c - u*u;
        if ( w2 < 0.0 ) return std::numeric_limits<T1>::max();

        T1 w = sqrt( w2 );

        T1 k = dot(v_b,v_c)/dot(v_c,v_c);
        sxyz<T1> pt;
        pt.x = vertexA->getX() + k*v_c.x;
        pt.y = vertexA->getY() + k*v_c.y;
        pt.z = vertexA->getZ() + k*v_c.z;

        T1 rho0 = vertexC->getDistance( pt );
        T1 xi0 = k;

        T1 xi = xi0 - u*rho0/(w*c);

        if ( 0.<xi && xi<1. ) {
            t = vertexA->getTT(threadNo) + u*xi0 + w*rho0/c;
        } else {
            t = std::numeric_limits<T1>::max();
        }

        return t;
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::local3Dsolver(NODE *vertexD,
                                             const size_t threadNo) const {

        // Méthode de Qian et al. 2007

        T2 iA, iB, iC, iD;
        NODE *vertexA, *vertexB, *vertexC;
        T1 AB, AC;

        for ( size_t no=0; no<vertexD->getOwners().size(); ++no ) {

            T2 tetNo = vertexD->getOwners()[no];

            for ( iD=0; iD<4; ++iD ) {
                if ( vertexD->getGridIndex() == tetrahedra[tetNo].i[iD] ) break;
            }

            iA = (iD+1)%4;
            iB = (iD+2)%4;
            iC = (iD+3)%4;
            vertexA = &(nodes[tetrahedra[tetNo].i[iA]]);
            vertexB = &(nodes[tetrahedra[tetNo].i[iB]]);
            vertexC = &(nodes[tetrahedra[tetNo].i[iC]]);

            if ( vertexA->getTT(threadNo) > vertexB->getTT(threadNo) ) {
                std::swap(iA, iB);
                std::swap(vertexA, vertexB);
            }
            if ( vertexA->getTT(threadNo) > vertexC->getTT(threadNo) ) {
                std::swap(iA, iC);
                std::swap(vertexA, vertexC);
            }

            if ( vertexA->getTT(threadNo) == std::numeric_limits<T1>::max() ) {
                continue;
            }

            AB = vertexA->getDistance( *vertexB );
            AC = vertexA->getDistance( *vertexC );

            bool apply2Dsolvers = true;

            if (std::abs(vertexB->getTT(threadNo)-vertexA->getTT(threadNo))<=AB*vertexD->getNodeSlowness() &&
                std::abs(vertexC->getTT(threadNo)-vertexA->getTT(threadNo))<=AC*vertexD->getNodeSlowness()) {

                // Qian et al, 2007, eq 2.3

                T1 ab[4], ac[4], n[2][3];

                // vec(AB)
                ab[0] = vertexB->getX()-vertexA->getX();
                ab[1] = vertexB->getY()-vertexA->getY();
                ab[2] = vertexB->getZ()-vertexA->getZ();

                ab[3] = (vertexB->getTT(threadNo)-vertexA->getTT(threadNo)) / vertexD->getNodeSlowness();

                // vec(AC)
                ac[0] = vertexC->getX()-vertexA->getX();
                ac[1] = vertexC->getY()-vertexA->getY();
                ac[2] = vertexC->getZ()-vertexA->getZ();

                ac[3] = (vertexC->getTT(threadNo)-vertexA->getTT(threadNo)) / vertexD->getNodeSlowness();

                int rv = solveEq23(ab, ac, n);

                if ( rv == 1 ) {

                    for ( size_t ns=0; ns<2; ++ns ) {
                        //
                        // find pt E
                        //

                        // plane vec(AB) cross vec(AC) passing by A: ax + by + cz + d = 0

                        T1 a = ab[1]*ac[2] - ac[1]*ab[2];
                        T1 b = ac[0]*ab[2] - ab[0]*ac[2];
                        T1 c = ab[0]*ac[1] - ac[0]*ab[1];

                        T1 d = -vertexA->getX()*a - vertexA->getY()*b - vertexA->getZ()*c;

                        T1 k = -(d + a*vertexD->getX() + b*vertexD->getY() + c*vertexD->getZ())/  // TODO check here if vertex D
                        (a*n[ns][0] + b*n[ns][1] + c*n[ns][2]);

                        sxyz<T1> E;
                        E.x = vertexD->getX() + k*n[ns][0];
                        E.y = vertexD->getY() + k*n[ns][1];
                        E.z = vertexD->getZ() + k*n[ns][2];

                        if ( testInTriangle(vertexA, vertexB, vertexC, E) ) {

                            // find point on wavefront plane

                            a = n[ns][0];
                            b = n[ns][1];
                            c = n[ns][2];
                            d = -vertexA->getX()*a - vertexA->getY()*b - vertexA->getZ()*c;

                            k = -(d + a*vertexD->getX() + b*vertexD->getY() + c*vertexD->getZ())/
                            (a*n[ns][0] + b*n[ns][1] + c*n[ns][2]);

                            sxyz<T1> pt;
                            pt.x = vertexD->getX() + k*n[ns][0];
                            pt.y = vertexD->getY() + k*n[ns][1];
                            pt.z = vertexD->getZ() + k*n[ns][2];

                            sxyz<T1> AD;
                            AD.x = vertexD->getX() - vertexA->getX();
                            AD.y = vertexD->getY() - vertexA->getY();
                            AD.z = vertexD->getZ() - vertexA->getZ();

                            T1 d2 = vertexD->getDistance( E );
                            T1 d3 = vertexD->getDistance( pt );
                            T1 d4 = std::abs( AD.x*n[ns][0] + AD.y*n[ns][1] + AD.z*n[ns][2] );

                            if ( std::abs(d3-d4)>small ) {
                                std::cout << " d3 ne d4: " << d3 << '\t' << d4 << '\t' << d2 << '\n';
                            }

                            T1 t = vertexA->getTT(threadNo) +
                            d3*vertexD->getNodeSlowness();

                            if ( t<vertexD->getTT(threadNo) )
                                vertexD->setTT(t, threadNo);

                            apply2Dsolvers = false;
                            break;
                        }
                    }
                }
            }

            if ( apply2Dsolvers ) {
                T1 tABD = local2Dsolver(vertexA, vertexB, vertexD, tetNo, threadNo);
                T1 tACD = local2Dsolver(vertexA, vertexC, vertexD, tetNo, threadNo);
                T1 tBCD = local2Dsolver(vertexB, vertexC, vertexD, tetNo, threadNo);

                T1 t = tABD < tACD ? tABD : tACD;
                t = t < tBCD ? t : tBCD;

                if ( t<vertexD->getTT(threadNo) )
                    vertexD->setTT(t, threadNo);
            }

        }
    }

    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dun<T1,T2,NODE>::local2Dsolver(const NODE *vertexA,
                                           const NODE *vertexB,
                                           const NODE *vertexC,
                                           const T2 tetraNo,
                                           const size_t threadNo) const {
        static const double pi2 = pi / 2.;

        if ( vertexB->getTT(threadNo)==std::numeric_limits<T1>::max() &&
            vertexA->getTT(threadNo)==std::numeric_limits<T1>::max() ) {
            return std::numeric_limits<T1>::max();
        }

        T1 t;

        T1 a = vertexB->getDistance( *vertexC );
        T1 b = vertexA->getDistance( *vertexC );
        T1 c = vertexA->getDistance( *vertexB );
        if ( std::abs(vertexB->getTT(threadNo)-vertexA->getTT(threadNo))<= c*vertexC->getNodeSlowness() ) {

            T1 theta = asin( std::abs(vertexB->getTT(threadNo)-vertexA->getTT(threadNo))/
                            (c*vertexC->getNodeSlowness()) );

            T1 gamma = acos((a*a + b*b - c*c)/(2.*a*b));

            if ( gamma > pi2 ) {
                std::cout << "*** Obtuse angle: " << gamma*57.2957795 << " ***\n";
            } else {
                std::cout << "Accute angle: " << gamma*57.2957795 << " \n";
            }

            T1 beta  = acos((b*b + c*c - a*a)/(2.*b*c));
            T1 alpha = acos((a*a + c*c - b*b)/(2.*a*c));

            if ( ((0.>alpha-pi2?0.:alpha-pi2)<=theta && theta<=(pi2-beta) ) ||
                ((alpha-pi2)<=theta && theta<=(0.<pi2-beta?0.:pi2-beta)) ) {
                T1 h = a*sin(alpha-theta);
                T1 H = b*sin(beta+theta);
                t = 0.5*(h*vertexC->getNodeSlowness() + vertexB->getTT(threadNo)) +
                0.5 *(H*vertexC->getNodeSlowness() + vertexA->getTT(threadNo));

            } else {
                t = vertexA->getTT(threadNo) + b*vertexC->getNodeSlowness();
                t = t < vertexB->getTT(threadNo) + a*vertexC->getNodeSlowness() ? t : vertexB->getTT(threadNo) + a*vertexC->getNodeSlowness();
            }
        } else {
            t = vertexA->getTT(threadNo) + b*vertexC->getNodeSlowness();
            t = t < vertexB->getTT(threadNo) + a*vertexC->getNodeSlowness() ? t : vertexB->getTT(threadNo) + a*vertexC->getNodeSlowness();
        }
        t = t<vertexC->getTT(threadNo) ? t : vertexC->getTT(threadNo);

        return t;
    }


    template<typename T1, typename T2, typename NODE>
    int Grid3Dun<T1,T2,NODE>::solveEq23(const T1 a[], const T1 b[], T1 n[][3]) const {
        // returns 0 if no solution
        //         1 if solutions exist

        T1 a02 = a[0]*a[0];
        T1 a12 = a[1]*a[1];
        T1 a22 = a[2]*a[2];
        T1 a32 = a[3]*a[3];
        T1 b02 = b[0]*b[0];
        T1 b12 = b[1]*b[1];
        T1 b22 = b[2]*b[2];
        T1 b32 = b[3]*b[3];
        T1 a23 = a[2]*a[2]*a[2];

        T1 s1 = (a[2]*b[1] - a[1]*b[2])*(a[2]*b[1] - a[1]*b[2])*
        (a02*(b12 + b22) - a32*(b02 + b12 + b22) + 2*a[0]*a[3]*b[0]*b[3] - a02*b32 -
         2*a[1]*b[1]*(a[0]*b[0] + a[2]*b[2] - a[3]*b[3]) + 2*a[2]*b[2]*
         (-(a[0]*b[0]) + a[3]*b[3]) + a22*(b02 + b12 - b32) + a12*(b02 + b22 - b32));

        if ( s1 < 0.0 ) {
            return 0;
        } else {

            T1 d1 = (a22*(b02 + b12) - 2*a[0]*a[2]*b[0]*b[2] -
                     2*a[1]*b[1]*(a[0]*b[0] + a[2]*b[2])  +
                     a12*(b02 + b22) + a02*(b12 + b22));
            T1 d2 = ((a[2]*b[1] - a[1]*b[2]) *(a22*(b02 + b12) -
                                               2*a[0]*a[2]*b[0]*b[2] -
                                               2*a[1]*b[1]*(a[0]*b[0] + a[2]*b[2]) +
                                               a12*(b02 + b22) + a02*(b12 + b22)));

            if ( d1==0.0 || d2==0.0 ) return 0;

            s1 = sqrt(s1);

            n[0][0] = (a[0]*a[3]*(b12 + b22) + a12*b[0]*b[3] - a[0]*a[2]*b[2]*b[3] -
                       a[1]*b[1]*(a[3]*b[0] + a[0]*b[3]) +
                       a[2]*b[0]*(-(a[3]*b[2]) + a[2]*b[3])  - s1) / d1;

            n[0][1] = (a[2]*a[3]*b[1]*(-(a[0]*b[0]*b[1])  + a[1]*(b02 + 2*b22)) +
                       a23*b12*b[3] + a[2]*(a[0]*b[1]*(-(a[1]*b[0])  + a[0]*b[1]) +
                                            a12*b22)*b[3] -
                       a22*b[1]*b[2]*(a[3]*b[1] + 2*a[1]*b[3]) -
                       a[1]*b[2]*(a[1]*a[3]*(b02 + b22) - a[0]*a[1]*b[0]*b[3] +
                                  a[0]*b[1]*(-(a[3]*b[0])  + a[0]*b[3]) ) +
                       a[2]*b[0]*s1 - a[0]*b[2]*s1) / d2;

            n[0][2] = (a[1]*b22*(a[3]*(a[0]*b[0] + a[1]*b[1])  - (a02 + a12)*b[3]) +
                       a22*b[1]*(a[3]*(b02 + b12) - (a[0]*b[0] + a[1]*b[1]) *b[3]) +
                       a[2]*b[2]*(-(a[1]*a[3]*(b02 + 2*b12)) + a[0]*a[1]*b[0]*b[3] +
                                  2*a12*b[1]*b[3] + a[0]*b[1]*(-(a[3]*b[0]) + a[0]*b[3]) ) -
                       a[1]*b[0]*s1 + a[0]*b[1]*s1)/ d2;

            n[1][0] = (a[0]*a[3]*(b12 + b22) + a12*b[0]*b[3] - a[0]*a[2]*b[2]*b[3] -
                       a[1]*b[1]*(a[3]*b[0] + a[0]*b[3]) +
                       a[2]*b[0]*(-(a[3]*b[2]) + a[2]*b[3])  + s1) / d1;

            n[1][1] = (a[2]*a[3]*b[1]*(-(a[0]*b[0]*b[1])  + a[1]*(b02 + 2*b22)) +
                       a23*b12*b[3] + a[2]*(a[0]*b[1]*(-(a[1]*b[0])  + a[0]*b[1]) +
                                            a12*b22)*b[3] - a22*b[1]*b[2]*(a[3]*b[1] +
                                                                           2*a[1]*b[3]) -
                       a[1]*b[2]*(a[1]*a[3]*(b02 + b22) - a[0]*a[1]*b[0]*b[3] +
                                  a[0]*b[1]*(-(a[3]*b[0])  + a[0]*b[3]) ) -
                       a[2]*b[0]*s1 + a[0]*b[2]*s1) / d2;

            n[1][2] = (a[1]*b22*(a[3]*(a[0]*b[0] + a[1]*b[1])  - (a02 + a12)*b[3]) +
                       a22*b[1]*(a[3]*(b02 + b12) - (a[0]*b[0] + a[1]*b[1]) *b[3]) +
                       a[2]*b[2]*(-(a[1]*a[3]*(b02 + 2*b12)) + a[0]*a[1]*b[0]*b[3] +
                                  2*a12*b[1]*b[3] + a[0]*b[1]*(-(a[3]*b[0])  + a[0]*b[3]) ) +
                       a[1]*b[0]*s1 - a[0]*b[1]*s1) / d2;
        }
        return 1;
    }


    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dun<T1,T2,NODE>::getTraveltimeFromRaypath(const std::vector<sxyz<T1>>& Tx,
                                                      const std::vector<T1>& t0,
                                                      const sxyz<T1> &Rx,
                                                      const size_t threadNo) const {
        T1 tt = 0.0;

        T1 minDist = small;
        T1 s1=0.0, s2=0.0;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                return t0[ns];
            }
        }

        const txInfo_t& txi = getTxInfo(Tx, threadNo);
        const std::vector<bool>& txOnNode = txi.txOnNode;
        const std::vector<bool>& txOnEdge = txi.txOnEdge;
        const std::vector<bool>& txOnFace = txi.txOnFace;
        const std::vector<T2>& txNode = txi.txNode;
        const std::vector<T2>& txCell = txi.txCell;
        const std::vector<std::array<T2,2>>& txEdges = txi.txEdges;
        const std::vector<std::array<T2,3>>& txFaces = txi.txFaces;
        const std::vector<std::vector<T2>>& txNeighborCells = txi.txNeighborCells;
#ifdef DEBUG_RP
        std::cout << "\n\n\n*** RP debug data - Source\n";
        std::vector<std::vector<sxyz<T1>>> r_data(1);
        r_data[0].push_back(Rx);
        for ( size_t nt=0; nt<Tx.size(); ++nt ) {
            std::cout << "   src no: " << nt << '\n';
            if ( txOnNode[nt] ) {
                std::cout << "     onNode\n"
                << "\t    txNode: " << txNode[nt] << '\n'
                << "\t    coords: " << nodes[txNode[nt]] << '\n';
            } else if (txOnEdge[nt] ) {
                std::cout << "     onEdge\n"
                << "\t    edge vertices no: " << txEdges[nt][0] << ' ' << txEdges[nt][1] << '\n'
                << "\t  vertices: " << nodes[txEdges[nt][0]] << '\n'
                << "\t          : " << nodes[txEdges[nt][1]] << '\n';
            } else if (txOnEdge[nt] ) {
                std::cout << "     onFace\n"
                << "\t    face vertices no: " << txFaces[nt][0] << ' ' << txFaces[nt][1] << ' ' << txFaces[nt][2] << '\n'
                << "\t  vertices: " << nodes[txFaces[nt][0]] << '\n'
                << "\t          : " << nodes[txFaces[nt][1]] << '\n'
                << "\t          : " << nodes[txFaces[nt][2]] << '\n';
            } else {
                std::array<T2,4> itmp = getPrimary(txCell[nt]);
                std::cout << "     inCell\n"
                << "\t    cellNo: " << txCell[nt] << '\n'
                << "\t  vertices: " << nodes[itmp[0]] << '\n'
                << "\t          : " << nodes[itmp[1]] << '\n'
                << "\t          : " << nodes[itmp[2]] << '\n'
                << "\t          : " << nodes[itmp[3]] << '\n';
            }
            std::cout << '\n';
        }
#endif

        T2 cellNo=0, nodeNo=0;
        sxyz<T1> curr_pt( Rx ), prev_pt( Rx );
        bool atRx = true;

        bool onNode = false;
        bool onEdge = false;
        bool onFace = false;
        std::array<T2,2> edgeNodes{ {0, 0} };
        std::array<T2,3> faceNodes{ {0, 0, 0} };
        Grad3D<T1,NODE>* grad3d = nullptr;
        if ( rp_method == 0 ) {
            grad3d = new Grad3D_ls_fo<T1,NODE>();
        } else if ( rp_method == 1 ) {
            grad3d = new Grad3D_ls_so<T1,NODE>();
        } else if ( rp_method == 2 ) {
            grad3d = new Grad3D_ab<T1,NODE>();
        }
        bool reachedTx = false;

        for ( T2 nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == curr_pt ) {
                nodeNo = nn;
                onNode = true;
                s1 = nodes[nodeNo].getNodeSlowness();
                break;
            }
        }
        if ( !onNode ) {
            cellNo = getCellNo( curr_pt );

            std::array<T2,4> itmp = getPrimary(cellNo);
            // find adjacent cells
            const T2 ind[6][2] = {
                {itmp[0], itmp[1]},
                {itmp[0], itmp[2]},
                {itmp[0], itmp[3]},
                {itmp[1], itmp[2]},
                {itmp[1], itmp[3]},
                {itmp[2], itmp[3]} };

            for ( size_t n=0; n<6; ++n ) {
                if ( areCollinear(curr_pt, nodes[ind[n][0]], nodes[ind[n][1]]) ) {
                    onEdge = true;
                    edgeNodes[0] = ind[n][0];
                    edgeNodes[1] = ind[n][1];
                    if ( processVel )
                        s1 = Interpolator<T1>::linearVel(curr_pt,
                                                         nodes[edgeNodes[0]],
                                                         nodes[edgeNodes[1]]);
                    else
                        s1 = Interpolator<T1>::linear(curr_pt,
                                                      nodes[edgeNodes[0]],
                                                      nodes[edgeNodes[1]]);
                    break;
                }
            }
            if ( !onEdge ) {
                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }

                for ( size_t n=0; n<4; ++n ) {
                    if ( areCoplanar(curr_pt, nodes[ind[n][0]], nodes[ind[n][1]], nodes[ind[n][2]]) ) {
                        onFace = true;
                        faceNodes[0] = ind[n][0];
                        faceNodes[1] = ind[n][1];
                        faceNodes[2] = ind[n][2];
                        if ( processVel )
                            s1 = Interpolator<T1>::bilinearTriangleVel(curr_pt,
                                                                       nodes[faceNodes[0]],
                                                                       nodes[faceNodes[1]],
                                                                       nodes[faceNodes[2]]);
                        else
                            s1 = Interpolator<T1>::bilinearTriangle(curr_pt,
                                                                    nodes[faceNodes[0]],
                                                                    nodes[faceNodes[1]],
                                                                    nodes[faceNodes[2]]);
                        //  faceNodes shoud not be assigned, face was not intersected
                        faceNodes = {0, 0, 0};
                        break;
                    }
                }
                if ( !onFace ) {
                    if ( processVel )
                        s1 = Interpolator<T1>::trilinearTriangleVel(curr_pt,
                                                                    nodes[itmp[0]],
                                                                    nodes[itmp[1]],
                                                                    nodes[itmp[2]],
                                                                    nodes[itmp[3]]);
                    else
                        s1 = Interpolator<T1>::trilinearTriangle(curr_pt,
                                                                 nodes[itmp[0]],
                                                                 nodes[itmp[1]],
                                                                 nodes[itmp[2]],
                                                                 nodes[itmp[3]]);

                }
            }
        }

        for ( size_t nt=0; nt<txCell.size(); ++nt ) {
            if ( cellNo == txCell[nt] ) {
                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( processVel )
                    s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                nodes[itmp[0]],
                                                                nodes[itmp[1]],
                                                                nodes[itmp[2]],
                                                                nodes[itmp[3]]);
                else
                    s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                             nodes[itmp[0]],
                                                             nodes[itmp[1]],
                                                             nodes[itmp[2]],
                                                             nodes[itmp[3]]);

                tt += t0[nt] + 0.5*(s1 + s2) * curr_pt.getDistance( Tx[nt] );
                reachedTx = true;
                break;
            }
        }

        sxyz<T1> g;
        while ( reachedTx == false ) {
            if ( onNode ) {
#ifdef DEBUG_RP
                printRaypathData(curr_pt, g, onNode, onEdge, onFace, cellNo,
                                 nodeNo, edgeNodes, faceNodes);
#endif
                if ( rp_method < 2 ) {
                    // find cells common to edge
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                        getNeighborNodes(*nc, nnodes);
                    }
                    // compute gradient with nodes from all common cells
                    g = grad3d->compute(curr_pt, nodes[nodeNo].getTT(threadNo), nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(1);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[nodeNo]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }

                // find cell for which gradient intersect opposing face
                bool foundIntersection = false;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                    std::array<T2,3> nb;
                    size_t n=0;
                    for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                        if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                            nb[n++] = *nn;
                        }
                    }
                    std::sort(nb.begin(), nb.end());

                    foundIntersection = intersectVecTriangle( nodeNo, g, nb[0], nb[1], nb[2], curr_pt);
                    if ( !foundIntersection ) {
                        continue;
                    }

                    //  check if cell is (one of) TxCell(s)
                    for (size_t nt=0; nt<Tx.size(); ++nt) {
                        if ( *nc == txCell[nt] ) {
                            std::array<T2,4> itmp = getPrimary(*nc);

                            if ( processVel )
                                s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                            nodes[itmp[0]],
                                                                            nodes[itmp[1]],
                                                                            nodes[itmp[2]],
                                                                            nodes[itmp[3]]);
                            else
                                s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                         nodes[itmp[0]],
                                                                         nodes[itmp[1]],
                                                                         nodes[itmp[2]],
                                                                         nodes[itmp[3]]);

                            tt += t0[nt] + 0.5*(s1 + s2) * curr_pt.getDistance( Tx[nt] );
#ifdef DEBUG_RP
                            r_data[0].push_back(Tx[nt]);
#endif
                            reachedTx = true;
                            break;
                        }
                    }
                    if ( reachedTx ) {
                        break;
                    }

                    bool break_flag = check_pt_location(curr_pt, nb, onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);
#ifdef DEBUG_RP
                    r_data[0].push_back(curr_pt);
#endif
                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
                    s1 = s2;
                    prev_pt = curr_pt;

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell1(faceNodes, nodeNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath (onNode) failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        reachedTx = true;
                    }
                    break;
                }

                if ( foundIntersection == false ) {
#ifdef DEBUG_RP
                    std::cout << "\n\nWarning: raypath (onNode) likely going outside mesh for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << '\n';
                    std::cout << "         Projecting gradient on external face and resuming" << std::endl;
#endif
                    // projet gradient on face and find intersection on next edge
                    sxyz<T1> pt_i;
                    foundIntersection = projectOnFace(curr_pt, nodeNo, g, edgeNodes, pt_i, threadNo);

                    if ( foundIntersection == false ) {
                        std::cout << "\n\nWarning: finding raypath (onNode) failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        reachedTx = true;
                        break;
                    }

                    // find which cell we are in
                    for ( size_t nc=0; nc<tetrahedra.size(); ++nc ) {
                        std::array<T2,4> itet = {tetrahedra[nc].i[0],
                            tetrahedra[nc].i[1],
                            tetrahedra[nc].i[2],
                            tetrahedra[nc].i[3]};
                        // because we are at the surface, there is only one cell with nodeNo & edgeNodes
                        if (std::find(itet.begin(), itet.end(), nodeNo) != itet.end() &&
                            std::find(itet.begin(), itet.end(), edgeNodes[0]) != itet.end() &&
                            std::find(itet.begin(), itet.end(), edgeNodes[1]) != itet.end()) {
                            cellNo = static_cast<T2>(nc);
                            break;
                        }
                    }

                    //  check if cell is (one of) TxCell(s)
                    for (size_t nt=0; nt<Tx.size(); ++nt) {
                        if ( cellNo == txCell[nt] ) {
                            std::array<T2,4> itmp = getPrimary(cellNo);

                            if ( processVel )
                                s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                            nodes[itmp[0]],
                                                                            nodes[itmp[1]],
                                                                            nodes[itmp[2]],
                                                                            nodes[itmp[3]]);
                            else
                                s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                         nodes[itmp[0]],
                                                                         nodes[itmp[1]],
                                                                         nodes[itmp[2]],
                                                                         nodes[itmp[3]]);

                            tt += t0[nt] + 0.5*(s1 + s2) * prev_pt.getDistance( Tx[nt] );
#ifdef DEBUG_RP
                            r_data[0].push_back(Tx[nt]);
#endif
                            reachedTx = true;
                            break;
                        }
                    }
                    if ( reachedTx ) {
                        break;
                    }
                    // we might be on one of the nodes
                    if ( pt_i.getDistance(nodes[edgeNodes[0]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[0];
                        pt_i = nodes[edgeNodes[0]];
                        onEdge = false;
                    } else if ( pt_i.getDistance(nodes[edgeNodes[1]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[1];
                        pt_i = nodes[edgeNodes[1]];
                        onEdge = false;
                    } else {
                        onEdge = true;
                        onNode = false;
                    }

                    curr_pt = pt_i;
#ifdef DEBUG_RP
                    r_data[0].push_back(curr_pt);
#endif

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
                    s1 = s2;
                    prev_pt = curr_pt;
                }

            } else if ( onEdge ) {
#ifdef DEBUG_RP
                printRaypathData(curr_pt, g, onNode, onEdge, onFace, cellNo,
                                 nodeNo, edgeNodes, faceNodes);
#endif
                // find cells common to edge
                std::vector<T2> cells;
                for ( auto nc0=nodes[edgeNodes[0]].getOwners().begin(); nc0!=nodes[edgeNodes[0]].getOwners().end(); ++nc0 ) {
                    if ( std::find(nodes[edgeNodes[1]].getOwners().begin(), nodes[edgeNodes[1]].getOwners().end(), *nc0)!=nodes[edgeNodes[1]].getOwners().end() ) {
                        cells.push_back( *nc0 );
                    }
                }
                if ( rp_method < 2 ) {
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    for (size_t n=0; n<cells.size(); ++n ) {
                        getNeighborNodes(cells[n], nnodes);
                    }
                    T1 d01 = nodes[edgeNodes[0]].getDistance(nodes[edgeNodes[1]]);
                    T1 w0 = curr_pt.getDistance(nodes[edgeNodes[1]]) / d01;
                    T1 w1 = curr_pt.getDistance(nodes[edgeNodes[0]]) / d01;
                    T1 curr_t = nodes[edgeNodes[0]].getTT(threadNo)*w0 + nodes[edgeNodes[1]].getTT(threadNo)*w1;
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(2);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[edgeNodes[0]]);
                    ref_pt[1] = &(nodes[edgeNodes[1]]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, edgeNodes, Tx, txCell);

                bool foundIntersection = false;
                for (size_t n=0; n<cells.size(); ++n ) {

                    T2 testCellNo = cells[n];

                    // there are 2 faces that might be intersected
                    std::array<T2,2> edgeNodes2;
                    size_t n2=0;
                    for ( auto nn=this->neighbors[testCellNo].begin(); nn!= this->neighbors[testCellNo].end(); ++nn ) {
                        if ( *nn!=edgeNodes[0] && *nn!=edgeNodes[1] && nodes[*nn].isPrimary() ) {
                            edgeNodes2[n2++] = *nn;
                        }
                    }

                    sxyz<T1> pt_i;
                    T2 itmpNode;
                    foundIntersection = intersectVecTriangle(curr_pt, g,
                                                             edgeNodes[0],
                                                             edgeNodes2[0],
                                                             edgeNodes2[1], pt_i);
                    itmpNode = edgeNodes[0];
                    if ( !foundIntersection ) {
                        foundIntersection = intersectVecTriangle(curr_pt, g,
                                                                 edgeNodes[1],
                                                                 edgeNodes2[0],
                                                                 edgeNodes2[1], pt_i);
                        itmpNode = edgeNodes[1];
                    }
                    if ( !foundIntersection ) {
                        continue;
                    }

                    cellNo = testCellNo;
                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, this->neighbors[cellNo],
                                                        {itmpNode, edgeNodes2[0], edgeNodes2[1]},
                                                        onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);
#ifdef DEBUG_RP
                    r_data[0].push_back(curr_pt);
#endif

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
                    s1 = s2;
                    prev_pt = curr_pt;

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath (onEdge) failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        reachedTx = true;
                    }
                    break;
                }
                if ( foundIntersection == false ) {
                    // we must be going outside the mesh
                    // hack: project gradient on face and find intersection
#ifdef DEBUG_RP
                    std::cout << "\n\nWarning: raypath (onEdge) likely going outside mesh for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << '\n';
                    std::cout << "         Projecting gradient on external face and resuming" << std::endl;
#endif
                    sxyz<T1> pt_i;
                    g = projectOnFace(curr_pt, g, edgeNodes, cells, pt_i, threadNo);
                    if ( g.x==0.0 && g.y==0.0 && g.z==0.0 ) {
                        foundIntersection = false;
                    } else {
                        foundIntersection = true;
                    }

                    if ( foundIntersection == false || curr_pt == pt_i ) {
                        std::cout << "\n\nWarning: finding raypath (onEdge) failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        reachedTx = true;
                        break;
                    }

                    // we might be on one of the nodes
                    if ( pt_i.getDistance(nodes[edgeNodes[0]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[0];
                        pt_i = nodes[edgeNodes[0]];
                        onEdge = false;
                    } else if ( pt_i.getDistance(nodes[edgeNodes[1]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[1];
                        pt_i = nodes[edgeNodes[1]];
                        onEdge = false;
                    }

                    curr_pt = pt_i;
#ifdef DEBUG_RP
                    r_data[0].push_back(curr_pt);
#endif

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
                    s1 = s2;
                    prev_pt = curr_pt;
                }

            } else if ( onFace ) { // on Face
#ifdef DEBUG_RP
                printRaypathData(curr_pt, g, onNode, onEdge, onFace, cellNo,
                                 nodeNo, edgeNodes, faceNodes);
#endif
                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( rp_method < 2 ) {
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    getNeighborNodes(cellNo, nnodes);
                    T1 curr_t;
                    if ( atRx ) {
                        curr_t = Interpolator<T1>::trilinearTime(curr_pt,
                                                                 nodes[itmp[0]],
                                                                 nodes[itmp[1]],
                                                                 nodes[itmp[2]],
                                                                 nodes[itmp[3]],
                                                                 threadNo);
                        atRx = false;
                    } else {
                        curr_t = Interpolator<T1>::bilinearTime(curr_pt,
                                                                nodes[faceNodes[0]],
                                                                nodes[faceNodes[1]],
                                                                nodes[faceNodes[2]],
                                                                threadNo);
                    }
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(3);
                    if ( atRx ) {
                        ref_pt[0] = &(nodes[itmp[0]]);
                        ref_pt[1] = &(nodes[itmp[1]]);
                        ref_pt[2] = &(nodes[itmp[2]]);
                        ref_pt.push_back( &(nodes[itmp[3]]) );
                        atRx = false;
                    } else {
                        ref_pt[0] = &(nodes[faceNodes[0]]);
                        ref_pt[1] = &(nodes[faceNodes[1]]);
                        ref_pt[2] = &(nodes[faceNodes[2]]);
                    }
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }
                // there are 3 faces that might be intersected

                bool foundIntersection = false;
                for ( size_t n=0; n<4; ++n ) {
                    if ( ind[n] == faceNodes ) continue;

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);

                    if ( !foundIntersection )
                        continue;

                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);
#ifdef DEBUG_RP
                    r_data[0].push_back(curr_pt);
#endif

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
                    s1 = s2;
                    prev_pt = curr_pt;

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath (onFace) failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        reachedTx = true;
                    }
                    break;
                }

                if ( foundIntersection == false ) {

                    // we must be on an face with gradient pointing slightly outward tetrahedron
                    // return in other cell but keep gradient
                    cellNo = findAdjacentCell2(faceNodes, cellNo);

                    std::array<T2,4> itmp = getPrimary(cellNo);
                    ind[0] = { { itmp[0], itmp[1], itmp[2] } };
                    ind[1] = { { itmp[0], itmp[1], itmp[3] } };
                    ind[2] = { { itmp[0], itmp[2], itmp[3] } };
                    ind[3] = { { itmp[1], itmp[2], itmp[3] } };

                    for ( size_t n=0; n<4; ++n ) {
                        std::sort( ind[n].begin(), ind[n].end() );
                    }

                    for ( size_t n=0; n<4; ++n ) {
                        if ( ind[n] == faceNodes ) continue;

                        sxyz<T1> pt_i;
                        foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                                 ind[n][1], ind[n][2],
                                                                 pt_i);

                        if ( !foundIntersection ) {
                            continue;
                        }
                        curr_pt = pt_i;

                        bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                            nodeNo, onEdge, edgeNodes,
                                                            onFace, faceNodes);
#ifdef DEBUG_RP
                        r_data[0].push_back(curr_pt);
#endif

                        s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                             faceNodes);

                        tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
                        s1 = s2;
                        prev_pt = curr_pt;

                        if ( break_flag ) break;

                        // find next cell
                        cellNo = findAdjacentCell2(faceNodes, cellNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            std::cout << "\n\nWarning: finding raypath (onFace) failed to converge (cell not found) for Rx "
                            << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                            tt = 0.0;
                            reachedTx = true;
                        }
                        break;
                    }
                }
                if ( foundIntersection == false ) {
                    // we must be going outside the mesh
                    // hack: project gradient on face and continue
#ifdef DEBUG_RP
                    std::cout << "\n\nWarning: raypath (onFace) likely going outside mesh for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << '\n';
                    std::cout << "         Projecting gradient on external face and resuming" << std::endl;
#endif
                    g = projectOnFace(g, faceNodes);

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecEdge(curr_pt, g, faceNodes, pt_i, edgeNodes);

                    if ( foundIntersection == false ) {
                        std::cout << "\n\nWarning: finding raypath (onFace) failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        reachedTx = true;
                        break;
                    }

                    // we might be on one of the nodes
                    if ( pt_i.getDistance(nodes[edgeNodes[0]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[0];
                        pt_i = nodes[edgeNodes[0]];
                        onEdge = false;
                        onFace = false;
                    } else if ( pt_i.getDistance(nodes[edgeNodes[1]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[1];
                        pt_i = nodes[edgeNodes[1]];
                        onEdge = false;
                        onFace = false;
                    } else {
                        onEdge = true;
                        onFace = false;
                    }
                    curr_pt = pt_i;
#ifdef DEBUG_RP
                    r_data[0].push_back(curr_pt);
#endif

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
                    s1 = s2;
                    prev_pt = curr_pt;
                }
            } else { // at Rx, somewhere in a tetrahedron
#ifdef DEBUG_RP
                printRaypathData(curr_pt, g, onNode, onEdge, onFace, cellNo,
                                 nodeNo, edgeNodes, faceNodes);
#endif
                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( rp_method < 2 ) {
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    getNeighborNodes(cellNo, nnodes);
                    T1 curr_t = Interpolator<T1>::trilinearTime(curr_pt,
                                                                nodes[itmp[0]],
                                                                nodes[itmp[1]],
                                                                nodes[itmp[2]],
                                                                nodes[itmp[3]],
                                                                threadNo);
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(4);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[itmp[0]]);
                    ref_pt[1] = &(nodes[itmp[1]]);
                    ref_pt[2] = &(nodes[itmp[2]]);
                    ref_pt[3] = &(nodes[itmp[3]]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }
                bool foundIntersection = false;
                for ( size_t n=0; n<4; ++n ) {

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);

                    if ( !foundIntersection )
                        continue;

                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);
#ifdef DEBUG_RP
                    r_data[0].push_back(curr_pt);
#endif

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    tt += 0.5*(s1 + s2) * prev_pt.getDistance( curr_pt );
                    s1 = s2;
                    prev_pt = curr_pt;

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath (inCell) failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        reachedTx = true;
                    }
                    break;
                }

                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath (inCell) failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    tt = 0.0;
                    reachedTx = true;
                }
            }

            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        if ( curr_pt.getDistance( Tx[nt] ) < minDist ) {
                            tt += t0[nt];
                            reachedTx = true;
                            break;
                        }
                    } else if ( txOnEdge[nt] ) {
                        if ( curr_pt.getDistance(nodes[txEdges[nt][0]]) < minDist ||
                            curr_pt.getDistance(nodes[txEdges[nt][1]]) < minDist ) {

                            if ( processVel )
                                s2 = Interpolator<T1>::linearVel(Tx[nt],
                                                                 nodes[txEdges[nt][0]],
                                                                 nodes[txEdges[nt][1]]);
                            else
                                s2 = Interpolator<T1>::linear(Tx[nt],
                                                              nodes[txEdges[nt][0]],
                                                              nodes[txEdges[nt][1]]);

                            tt += t0[nt] + 0.5*(s1 + s2) * prev_pt.getDistance( Tx[nt] );
                            reachedTx = true;
#ifdef DEBUG_RP
                            r_data[0].push_back(Tx[nt]);
#endif
                            break;
                        }
                    } else if ( txOnFace[nt] ) {
                        if ( curr_pt.getDistance(nodes[txFaces[nt][0]]) < minDist ||
                            curr_pt.getDistance(nodes[txFaces[nt][1]]) < minDist ||
                            curr_pt.getDistance(nodes[txFaces[nt][2]]) < minDist ) {

                            if ( processVel )
                                s2 = Interpolator<T1>::bilinearTriangleVel(Tx[nt],
                                                                           nodes[txFaces[nt][0]],
                                                                           nodes[txFaces[nt][1]],
                                                                           nodes[txFaces[nt][2]]);
                            else
                                s2 = Interpolator<T1>::bilinearTriangle(Tx[nt],
                                                                        nodes[txFaces[nt][0]],
                                                                        nodes[txFaces[nt][1]],
                                                                        nodes[txFaces[nt][2]]);

                            tt += t0[nt] + 0.5*(s1 + s2) * prev_pt.getDistance( Tx[nt] );
                            reachedTx = true;
#ifdef DEBUG_RP
                            r_data[0].push_back(Tx[nt]);
#endif
                            break;
                        }
                    }
                }
            } if ( onEdge ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        if ( txNode[nt] == edgeNodes[0] || txNode[nt] == edgeNodes[1] ) {
                            s2 = nodes[txNode[nt]].getNodeSlowness();
                            tt += t0[nt] + 0.5*(s1 + s2) * prev_pt.getDistance( Tx[nt] );
                            reachedTx = true;
#ifdef DEBUG_RP
                            r_data[0].push_back(Tx[nt]);
#endif
                            break;
                        }
                    } else {
                        std::array<T2,4> itmp = getPrimary(txCell[nt]);
                        // find adjacent cells
                        const T2 ind[6][2] = {
                            {itmp[0], itmp[1]},
                            {itmp[0], itmp[2]},
                            {itmp[0], itmp[3]},
                            {itmp[1], itmp[2]},
                            {itmp[1], itmp[3]},
                            {itmp[2], itmp[3]} };
                        for ( size_t ne=0; ne<6; ++ne ) {
                            if ( (ind[ne][0] == edgeNodes[0] && ind[ne][1] == edgeNodes[1]) ||
                                (ind[ne][0] == edgeNodes[1] && ind[ne][1] == edgeNodes[0]) ) {
                                if ( processVel )
                                    s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                nodes[itmp[0]],
                                                                                nodes[itmp[1]],
                                                                                nodes[itmp[2]],
                                                                                nodes[itmp[3]]);
                                else
                                    s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                             nodes[itmp[0]],
                                                                             nodes[itmp[1]],
                                                                             nodes[itmp[2]],
                                                                             nodes[itmp[3]]);

                                tt += t0[nt] + 0.5*(s1 + s2) * prev_pt.getDistance( Tx[nt] );
                                reachedTx = true;
#ifdef DEBUG_RP
                                r_data[0].push_back(Tx[nt]);
#endif
                                break;
                            }
                        }
                    }
                }
            } else {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        for ( auto nc=nodes[txNode[nt]].getOwners().begin(); nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                            if ( cellNo == *nc ) {
                                std::array<T2,4> itmp = getPrimary(cellNo);

                                if ( processVel )
                                    s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                nodes[itmp[0]],
                                                                                nodes[itmp[1]],
                                                                                nodes[itmp[2]],
                                                                                nodes[itmp[3]]);
                                else
                                    s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                             nodes[itmp[0]],
                                                                             nodes[itmp[1]],
                                                                             nodes[itmp[2]],
                                                                             nodes[itmp[3]]);

                                tt += t0[nt] + 0.5*(s1 + s2) * prev_pt.getDistance( Tx[nt] );
#ifdef DEBUG_RP
                                r_data[0].push_back(Tx[nt]);
#endif
                                reachedTx = true;
                                break;
                            }
                        }
                    } else {
                        if ( cellNo == txCell[nt] ) {
                            std::array<T2,4> itmp = getPrimary(cellNo);

                            if ( processVel )
                                s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                            nodes[itmp[0]],
                                                                            nodes[itmp[1]],
                                                                            nodes[itmp[2]],
                                                                            nodes[itmp[3]]);
                            else
                                s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                         nodes[itmp[0]],
                                                                         nodes[itmp[1]],
                                                                         nodes[itmp[2]],
                                                                         nodes[itmp[3]]);

                            tt += t0[nt] + 0.5*(s1 + s2) * prev_pt.getDistance( Tx[nt] );
#ifdef DEBUG_RP
                            r_data[0].push_back(Tx[nt]);
#endif
                            reachedTx = true;
                        } else {
                            for ( size_t nn=0; nn<txNeighborCells[nt].size(); ++nn ) {
                                if ( cellNo == txNeighborCells[nt][nn] ) {
                                    std::array<T2,4> itmp = getPrimary(txCell[nt]);
                                    std::array<T2,3> ind[4] = {
                                        { { itmp[0], itmp[1], itmp[2] } },
                                        { { itmp[0], itmp[1], itmp[3] } },
                                        { { itmp[0], itmp[2], itmp[3] } },
                                        { { itmp[1], itmp[2], itmp[3] } }
                                    };
                                    bool found = false;
                                    for ( size_t n=0; n<4; ++n ) {
                                        std::sort( ind[n].begin(), ind[n].end() );
                                        if ( faceNodes == ind[n] ) {
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found ) {
                                        std::array<T2,4> itmp = getPrimary(cellNo);

                                        if ( processVel )
                                            s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                        nodes[itmp[0]],
                                                                                        nodes[itmp[1]],
                                                                                        nodes[itmp[2]],
                                                                                        nodes[itmp[3]]);
                                        else
                                            s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                                     nodes[itmp[0]],
                                                                                     nodes[itmp[1]],
                                                                                     nodes[itmp[2]],
                                                                                     nodes[itmp[3]]);

                                        tt += t0[nt] + 0.5*(s1 + s2) * prev_pt.getDistance( Tx[nt] );
                                        reachedTx = true;
#ifdef DEBUG_RP
                                        r_data[0].push_back(Tx[nt]);
#endif
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if ( reachedTx ) break;
                }
            }
        }
#ifdef DEBUG_RP
        std::ostringstream fname;
        fname << "raypath_" << Rx.x << '_' << Rx.y << '_' << Rx.z << ".vtp";
        saveRayPaths(fname.str(), r_data);
#endif
        delete grad3d;
        return tt;
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::printRaypathData(const sxyz<T1>& curr_pt,
                                                const sxyz<T1>& g,
                                                const bool onNode,
                                                const bool onEdge,
                                                const bool onFace,
                                                const T2 cellNo,
                                                const T2 nodeNo,
                                                const std::array<T2,2> &edgeNodes,
                                                const std::array<T2,3> &faceNodes) const {
        std::array<T2,4> itmp = getPrimary(cellNo);
        std::cout << "\n*** RP debug data\n   curr_pt: " << curr_pt << '\n'
        << "         g: " << g << '\n'
        << "    cellNo: " << cellNo << '\n'
        << "  vertices: " << nodes[itmp[0]] << '\n'
        << "          : " << nodes[itmp[1]] << '\n'
        << "          : " << nodes[itmp[2]] << '\n'
        << "          : " << nodes[itmp[3]] << '\n';
        if ( onNode ) {
            std::cout << "\tonNode\n"
            << "\t    nodeNo: " << nodeNo << '\n'
            << "\t    coords: " << nodes[nodeNo] << '\n';
        }
        if ( onEdge ) {
            std::cout << "\tonEdge\n"
            << "\t    edgeNo: " << edgeNodes[0] << ' ' << edgeNodes[1] << '\n'
            << "\t    coords: " << nodes[edgeNodes[0]] << '\n'
            << "\t    coords: " << nodes[edgeNodes[1]] << '\n';
        }
        if ( onFace ) {
            std::cout << "\tonFace\n"
            << "\t    faceNo: " << faceNodes[0] << ' ' << faceNodes[1] << ' ' << faceNodes[2] << '\n'
            << "\t    coords: " << nodes[faceNodes[0]] << '\n'
            << "\t    coords: " << nodes[faceNodes[1]] << '\n'
            << "\t    coords: " << nodes[faceNodes[2]] << '\n';
        }
        std::cout << std::endl;
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                          const sxyz<T1> &Rx,
                                          std::vector<sxyz<T1>> &r_data,
                                          const size_t threadNo) const {

        T1 minDist = small;
        std::vector<sxyz<T1>> r_tmp;
        r_tmp.emplace_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                return;
            }
        }

        const txInfo_t& txi = getTxInfo(Tx, threadNo);
        const std::vector<bool>& txOnNode = txi.txOnNode;
        const std::vector<bool>& txOnEdge = txi.txOnEdge;
        const std::vector<bool>& txOnFace = txi.txOnFace;
        const std::vector<T2>& txNode = txi.txNode;
        const std::vector<T2>& txCell = txi.txCell;
        const std::vector<std::array<T2,2>>& txEdges = txi.txEdges;
        const std::vector<std::array<T2,3>>& txFaces = txi.txFaces;
        const std::vector<std::vector<T2>>& txNeighborCells = txi.txNeighborCells;

        T2 cellNo, nodeNo;
        sxyz<T1> curr_pt( Rx );

        bool onNode = false;
        bool onEdge = false;
        bool onFace = false;
        std::array<T2,2> edgeNodes{ {0, 0} };
        std::array<T2,3> faceNodes{ {0, 0, 0} };
        Grad3D<T1,NODE>* grad3d = nullptr;
        if ( rp_method == 0 ) {
            grad3d = new Grad3D_ls_fo<T1,NODE>();
        } else if ( rp_method == 1 ) {
            grad3d = new Grad3D_ls_so<T1,NODE>();
        } else if ( rp_method == 2 ) {
            grad3d = new Grad3D_ab<T1,NODE>();
        }
        bool reachedTx = false;

        for ( T2 nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == curr_pt ) {
                nodeNo = nn;
                onNode = true;
                break;
            }
        }
        if ( !onNode ) {
            cellNo = getCellNo( curr_pt );

            std::array<T2,4> itmp = getPrimary(cellNo);
            T2 ind[6][2] = {
                {itmp[0], itmp[1]},
                {itmp[0], itmp[2]},
                {itmp[0], itmp[3]},
                {itmp[1], itmp[2]},
                {itmp[1], itmp[3]},
                {itmp[2], itmp[3]} };

            for ( size_t n=0; n<6; ++n ) {
                if ( areCollinear(curr_pt, nodes[ind[n][0]], nodes[ind[n][1]]) ) {
                    onEdge = true;
                    edgeNodes[0] = ind[n][0];
                    edgeNodes[1] = ind[n][1];
                    break;
                }
            }
            if ( !onEdge ) {
                std::array<T2,4> itmp = getPrimary(cellNo);
                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }

                for ( size_t n=0; n<4; ++n ) {
                    if ( areCoplanar(curr_pt, nodes[ind[n][0]], nodes[ind[n][1]], nodes[ind[n][2]]) ) {
                        onFace = true;
                        //  faceNodes shoud not be assigned, face was not intersected
                        break;
                    }
                }
            }
        }

        for ( auto nt=0; nt<txCell.size(); ++nt ) {
            if ( cellNo == txCell[nt] ) {
                r_tmp.push_back( Tx[nt] );
                reachedTx = true;
                break;
            }
        }

        sxyz<T1> g;
        while ( reachedTx == false ) {

            if ( onNode ) {

                if ( rp_method < 2 ) {
                    // find cells common to edge
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                        getNeighborNodes(*nc, nnodes);
                    }
                    // compute gradient with nodes from all common cells
                    g = grad3d->compute(curr_pt, nodes[nodeNo].getTT(threadNo), nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(1);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[nodeNo]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }

                // find cell for which gradient intersect opposing face
                bool foundIntersection = false;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                    std::array<T2,3> nb;
                    size_t n=0;
                    for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                        if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                            nb[n++] = *nn;
                        }
                    }
                    std::sort(nb.begin(), nb.end());

                    foundIntersection = intersectVecTriangle( nodeNo, g, nb[0], nb[1], nb[2], curr_pt);
                    if ( !foundIntersection ) {
                        continue;
                    }

                    //  check if cell is (one of) TxCell(s)
                    for (size_t nt=0; nt<Tx.size(); ++nt) {
                        if ( *nc == txCell[nt] ) {
                            std::array<T2,4> itmp = getPrimary(*nc);

                            r_tmp.push_back( Tx[nt] );
                            reachedTx = true;
                            break;
                        }
                    }
                    if ( reachedTx ) {
                        break;
                    }

                    bool break_flag = check_pt_location(curr_pt, nb, onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    r_tmp.push_back( curr_pt );

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell1(faceNodes, nodeNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }

                if ( foundIntersection == false ) {
                    // projet gradient on face and find intersection on next edge
                    sxyz<T1> pt_i;
                    foundIntersection = projectOnFace(curr_pt, nodeNo, g, edgeNodes, pt_i, threadNo);

                    if ( foundIntersection == false ) {
                        std::cout << "\n\nWarning: finding raypath (onNode) failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                        break;
                    }

                    // find which cell we are in
                    for ( size_t nc=0; nc<tetrahedra.size(); ++nc ) {
                        std::array<T2,4> itet = {tetrahedra[nc].i[0],
                            tetrahedra[nc].i[1],
                            tetrahedra[nc].i[2],
                            tetrahedra[nc].i[3]};
                        // because we are at the surface, there is only one cell with nodeNo & edgeNodes
                        if (std::find(itet.begin(), itet.end(), nodeNo) != itet.end() &&
                            std::find(itet.begin(), itet.end(), edgeNodes[0]) != itet.end() &&
                            std::find(itet.begin(), itet.end(), edgeNodes[1]) != itet.end()) {
                            cellNo = static_cast<T2>(nc);
                            break;
                        }
                    }

                    //  check if cell is (one of) TxCell(s)
                    for (size_t nt=0; nt<Tx.size(); ++nt) {
                        if ( cellNo == txCell[nt] ) {
                            std::array<T2,4> itmp = getPrimary(cellNo);

                            r_tmp.push_back( Tx[nt] );
                            reachedTx = true;
                            break;
                        }
                    }
                    if ( reachedTx ) {
                        break;
                    }
                    // we might be on one of the nodes
                    if ( pt_i.getDistance(nodes[edgeNodes[0]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[0];
                        pt_i = nodes[edgeNodes[0]];
                        onEdge = false;
                    } else if ( pt_i.getDistance(nodes[edgeNodes[1]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[1];
                        pt_i = nodes[edgeNodes[1]];
                        onEdge = false;
                    } else {
                        onEdge = true;
                        onNode = false;
                    }

                    curr_pt = pt_i;
                    r_tmp.push_back( curr_pt );
                }

            } else if ( onEdge ) {

                // find cells common to edge
                std::vector<T2> cells;
                for ( auto nc0=nodes[edgeNodes[0]].getOwners().begin(); nc0!=nodes[edgeNodes[0]].getOwners().end(); ++nc0 ) {
                    if ( std::find(nodes[edgeNodes[1]].getOwners().begin(), nodes[edgeNodes[1]].getOwners().end(), *nc0)!=nodes[edgeNodes[1]].getOwners().end() ) {
                        cells.push_back( *nc0 );
                    }
                }
                if ( rp_method < 2 ) {
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    for (size_t n=0; n<cells.size(); ++n ) {
                        getNeighborNodes(cells[n], nnodes);
                    }
                    T1 d01 = nodes[edgeNodes[0]].getDistance(nodes[edgeNodes[1]]);
                    T1 w0 = curr_pt.getDistance(nodes[edgeNodes[1]]) / d01;
                    T1 w1 = curr_pt.getDistance(nodes[edgeNodes[0]]) / d01;
                    T1 curr_t = nodes[edgeNodes[0]].getTT(threadNo)*w0 + nodes[edgeNodes[1]].getTT(threadNo)*w1;
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(2);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[edgeNodes[0]]);
                    ref_pt[1] = &(nodes[edgeNodes[1]]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, edgeNodes, Tx, txCell);

                bool foundIntersection = false;
                for (size_t n=0; n<cells.size(); ++n ) {

                    T2 testCellNo = cells[n];

                    // there are 2 faces that might be intersected
                    std::array<T2,2> edgeNodes2;
                    size_t n2=0;
                    for ( auto nn=this->neighbors[testCellNo].begin(); nn!= this->neighbors[testCellNo].end(); ++nn ) {
                        if ( *nn!=edgeNodes[0] && *nn!=edgeNodes[1] && nodes[*nn].isPrimary() ) {
                            edgeNodes2[n2++] = *nn;
                        }
                    }

                    sxyz<T1> pt_i;
                    T2 itmpNode;
                    foundIntersection = intersectVecTriangle(curr_pt, g,
                                                             edgeNodes[0],
                                                             edgeNodes2[0],
                                                             edgeNodes2[1], pt_i);
                    itmpNode = edgeNodes[0];
                    if ( !foundIntersection ) {
                        foundIntersection = intersectVecTriangle(curr_pt, g,
                                                                 edgeNodes[1],
                                                                 edgeNodes2[0],
                                                                 edgeNodes2[1], pt_i);
                        itmpNode = edgeNodes[1];
                    }
                    if ( !foundIntersection ) {
                        continue;
                    }

                    cellNo = testCellNo;
                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, this->neighbors[cellNo],
                                                        {itmpNode, edgeNodes2[0], edgeNodes2[1]},
                                                        onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    r_tmp.push_back( curr_pt );

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }
                if ( foundIntersection == false ) {
                    sxyz<T1> pt_i;
                    g = projectOnFace(curr_pt, g, edgeNodes, cells, pt_i,threadNo);
                    if ( g.x==0.0 && g.y==0.0 && g.z==0.0 ) {
                        foundIntersection = false;
                    } else {
                        foundIntersection = true;
                    }

                    if ( foundIntersection == false || curr_pt == pt_i ) {
                        std::cout << "\n\nWarning: finding raypath (onEdge) failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                        break;
                    }

                    // we might be on one of the nodes
                    if ( pt_i.getDistance(nodes[edgeNodes[0]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[0];
                        pt_i = nodes[edgeNodes[0]];
                        onEdge = false;
                    } else if ( pt_i.getDistance(nodes[edgeNodes[1]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[1];
                        pt_i = nodes[edgeNodes[1]];
                        onEdge = false;
                    }

                    curr_pt = pt_i;
                    r_tmp.push_back(curr_pt);
                }

            } else if ( onFace ) { // on Face

                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( rp_method < 2 ) {
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    getNeighborNodes(cellNo, nnodes);
                    T1 curr_t;
                    if ( r_tmp.size() <= 1 ) {
                        curr_t = Interpolator<T1>::trilinearTime(curr_pt,
                                                                 nodes[itmp[0]],
                                                                 nodes[itmp[1]],
                                                                 nodes[itmp[2]],
                                                                 nodes[itmp[3]],
                                                                 threadNo);
                    } else {
                        curr_t = Interpolator<T1>::bilinearTime(curr_pt,
                                                                nodes[faceNodes[0]],
                                                                nodes[faceNodes[1]],
                                                                nodes[faceNodes[2]],
                                                                threadNo);
                    }
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(3);
                    if ( r_tmp.size() <= 1 ) {
                        ref_pt[0] = &(nodes[itmp[0]]);
                        ref_pt[1] = &(nodes[itmp[1]]);
                        ref_pt[2] = &(nodes[itmp[2]]);
                        ref_pt.push_back( &(nodes[itmp[3]]) );
                    } else {
                        ref_pt[0] = &(nodes[faceNodes[0]]);
                        ref_pt[1] = &(nodes[faceNodes[1]]);
                        ref_pt[2] = &(nodes[faceNodes[2]]);
                    }
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }
                // there are 3 faces that might be intersected

                bool foundIntersection = false;
                for ( size_t n=0; n<4; ++n ) {
                    if ( ind[n] == faceNodes ) continue;

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);

                    if ( !foundIntersection )
                        continue;

                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    r_tmp.push_back( curr_pt );

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }

                if ( foundIntersection == false ) {

                    // we must be on an face with gradient pointing slightly outward tetrahedron
                    // return in other cell but keep gradient
                    cellNo = findAdjacentCell2(faceNodes, cellNo);

                    std::array<T2,4> itmp = getPrimary(cellNo);
                    ind[0] = { { itmp[0], itmp[1], itmp[2] } };
                    ind[1] = { { itmp[0], itmp[1], itmp[3] } };
                    ind[2] = { { itmp[0], itmp[2], itmp[3] } };
                    ind[3] = { { itmp[1], itmp[2], itmp[3] } };

                    for ( size_t n=0; n<4; ++n ) {
                        std::sort( ind[n].begin(), ind[n].end() );
                    }

                    for ( size_t n=0; n<4; ++n ) {
                        if ( ind[n] == faceNodes ) continue;

                        sxyz<T1> pt_i;
                        foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                                 ind[n][1], ind[n][2],
                                                                 pt_i);

                        if ( !foundIntersection ) {
                            continue;
                        }
                        curr_pt = pt_i;

                        bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                            nodeNo, onEdge, edgeNodes,
                                                            onFace, faceNodes);

                        r_tmp.push_back( curr_pt );

                        if ( break_flag ) break;

                        // find next cell
                        cellNo = findAdjacentCell2(faceNodes, cellNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                            << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                            r_tmp.resize(1);
                            r_tmp[0] = Rx;
                            reachedTx = true;
                        }
                        break;
                    }
                }
                if ( foundIntersection == false ) {
                    g = projectOnFace(g, faceNodes);

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecEdge(curr_pt, g, faceNodes, pt_i, edgeNodes);

                    if ( foundIntersection == false ) {
                        std::cout << "\n\nWarning: finding raypath (onFace) failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                        break;
                    }

                    // we might be on one of the nodes
                    if ( pt_i.getDistance(nodes[edgeNodes[0]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[0];
                        pt_i = nodes[edgeNodes[0]];
                        onEdge = false;
                        onFace = false;
                    } else if ( pt_i.getDistance(nodes[edgeNodes[1]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[1];
                        pt_i = nodes[edgeNodes[1]];
                        onEdge = false;
                        onFace = false;
                    } else {
                        onEdge = true;
                        onFace = false;
                    }
                    curr_pt = pt_i;
                    r_tmp.push_back( curr_pt );
                }
            } else { // at Rx, somewhere in a tetrahedron

                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( rp_method < 2 ) {
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    getNeighborNodes(cellNo, nnodes);
                    T1 curr_t = Interpolator<T1>::trilinearTime(curr_pt,
                                                                nodes[itmp[0]],
                                                                nodes[itmp[1]],
                                                                nodes[itmp[2]],
                                                                nodes[itmp[3]],
                                                                threadNo);
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(4);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[itmp[0]]);
                    ref_pt[1] = &(nodes[itmp[1]]);
                    ref_pt[2] = &(nodes[itmp[2]]);
                    ref_pt[3] = &(nodes[itmp[3]]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }
                // there are 3 faces that might be intersected

                bool foundIntersection = false;
                for ( size_t n=0; n<4; ++n ) {

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);

                    if ( !foundIntersection )
                        continue;

                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    r_tmp.push_back( curr_pt );

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }
                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath within cell failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_data.resize(1);
                    r_data[0] = Rx;
                    reachedTx = true;
                }
            }

            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        if ( curr_pt.getDistance( Tx[nt] ) < minDist ) {
                            reachedTx = true;
                            break;
                        }
                    } else if ( txOnEdge[nt] ) {
                        if ( curr_pt.getDistance(nodes[txEdges[nt][0]]) < minDist ||
                            curr_pt.getDistance(nodes[txEdges[nt][1]]) < minDist ) {

                            r_tmp.push_back(Tx[nt]);
                            reachedTx = true;
                            break;
                        }
                    } else if ( txOnFace[nt] ) {
                        if ( curr_pt.getDistance(nodes[txFaces[nt][0]]) < minDist ||
                            curr_pt.getDistance(nodes[txFaces[nt][1]]) < minDist ||
                            curr_pt.getDistance(nodes[txFaces[nt][2]]) < minDist ) {

                            r_tmp.push_back(Tx[nt]);
                            reachedTx = true;
                            break;
                        }
                    }
                }
            } if ( onEdge ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        if ( txNode[nt] == edgeNodes[0] || txNode[nt] == edgeNodes[1] ) {
                            r_tmp.push_back( Tx[nt] );
                            reachedTx = true;
                            break;
                        }
                    } else {
                        std::array<T2,4> itmp = getPrimary(txCell[nt]);
                        // find adjacent cells
                        const T2 ind[6][2] = {
                            {itmp[0], itmp[1]},
                            {itmp[0], itmp[2]},
                            {itmp[0], itmp[3]},
                            {itmp[1], itmp[2]},
                            {itmp[1], itmp[3]},
                            {itmp[2], itmp[3]} };
                        for ( size_t ne=0; ne<6; ++ne ) {
                            if ( (ind[ne][0] == edgeNodes[0] && ind[ne][1] == edgeNodes[1]) ||
                                (ind[ne][0] == edgeNodes[1] && ind[ne][1] == edgeNodes[0]) ) {
                                r_tmp.push_back( Tx[nt] );
                                reachedTx = true;
                                break;
                            }
                        }
                    }
                }
            } else {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        for ( auto nc=nodes[txNode[nt]].getOwners().begin(); nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                            if ( cellNo == *nc ) {
                                r_tmp.push_back( Tx[nt] );
                                reachedTx = true;
                                break;
                            }
                        }
                    } else {
                        if ( cellNo == txCell[nt] ) {
                            r_tmp.push_back( Tx[nt] );
                            reachedTx = true;
                        } else {
                            for ( size_t nn=0; nn<txNeighborCells[nt].size(); ++nn ) {
                                if ( cellNo == txNeighborCells[nt][nn] ) {
                                    std::array<T2,4> itmp = getPrimary(txCell[nt]);
                                    std::array<T2,3> ind[4] = {
                                        { { itmp[0], itmp[1], itmp[2] } },
                                        { { itmp[0], itmp[1], itmp[3] } },
                                        { { itmp[0], itmp[2], itmp[3] } },
                                        { { itmp[1], itmp[2], itmp[3] } }
                                    };
                                    bool found = false;
                                    for ( size_t n=0; n<4; ++n ) {
                                        std::sort( ind[n].begin(), ind[n].end() );
                                        if ( faceNodes == ind[n] ) {
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found ) {
                                        r_tmp.push_back( Tx[nt] );
                                        reachedTx = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if ( reachedTx ) break;
                }
            }
        }
        // for inversion, the order should be from Tx to Rx, so we reorder...
        size_t npts = r_tmp.size();
        r_data.resize( npts );
        for ( size_t nn=0; nn<npts; ++nn ) {
            r_data[nn] = r_tmp[ npts-1-nn ];
        }

        delete grad3d;
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                          const std::vector<T1>& t0,
                                          const sxyz<T1> &Rx,
                                          std::vector<sxyz<T1>> &r_data,
                                          T1 &tt,
                                          const size_t threadNo) const {

        T1 minDist = small;
        std::vector<sxyz<T1>> r_tmp;
        r_tmp.emplace_back( Rx );
        tt = 0.0;
        T1 s1=0.0, s2=0.0;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        const txInfo_t& txi = getTxInfo(Tx, threadNo);
        const std::vector<bool>& txOnNode = txi.txOnNode;
        const std::vector<bool>& txOnEdge = txi.txOnEdge;
        const std::vector<bool>& txOnFace = txi.txOnFace;
        const std::vector<T2>& txNode = txi.txNode;
        const std::vector<T2>& txCell = txi.txCell;
        const std::vector<std::array<T2,2>>& txEdges = txi.txEdges;
        const std::vector<std::array<T2,3>>& txFaces = txi.txFaces;
        const std::vector<std::vector<T2>>& txNeighborCells = txi.txNeighborCells;
#ifdef DEBUG_RP
        std::cout << "\n*** RP debug data - Source\n";
        for ( size_t nt=0; nt<Tx.size(); ++nt ) {
            std::cout << "   src no: " << nt << '\n';
            if ( txOnNode[nt] ) {
                std::cout << "     onNode\n"
                << "\t    txNode: " << txNode[nt] << '\n'
                << "\t    coords: " << nodes[txNode[nt]] << '\n';
            } else if (txOnEdge[nt] ) {
                std::cout << "     onEdge\n"
                << "\t    edge vertices no: " << txEdges[nt][0] << ' ' << txEdges[nt][1] << '\n'
                << "\t  vertices: " << nodes[txEdges[nt][0]] << '\n'
                << "\t          : " << nodes[txEdges[nt][1]] << '\n';
            } else if (txOnEdge[nt] ) {
                std::cout << "     onFace\n"
                << "\t    face vertices no: " << txFaces[nt][0] << ' ' << txFaces[nt][1] << ' ' << txFaces[nt][2] << '\n'
                << "\t  vertices: " << nodes[txFaces[nt][0]] << '\n'
                << "\t          : " << nodes[txFaces[nt][1]] << '\n'
                << "\t          : " << nodes[txFaces[nt][2]] << '\n';
            } else {
                std::array<T2,4> itmp = getPrimary(txCell[nt]);
                std::cout << "     inCell\n"
                << "\t    cellNo: " << txCell[nt] << '\n'
                << "\t  vertices: " << nodes[itmp[0]] << '\n'
                << "\t          : " << nodes[itmp[1]] << '\n'
                << "\t          : " << nodes[itmp[2]] << '\n'
                << "\t          : " << nodes[itmp[3]] << '\n';
            }
            std::cout << '\n';
        }
#endif

        T2 cellNo=0, nodeNo=0;
        sxyz<T1> curr_pt( Rx );

        bool onNode = false;
        bool onEdge = false;
        bool onFace = false;
        std::array<T2,2> edgeNodes{ {0, 0} };
        std::array<T2,3> faceNodes{ {0, 0, 0} };
        Grad3D<T1,NODE>* grad3d = nullptr;
        if ( rp_method == 0 ) {
            grad3d = new Grad3D_ls_fo<T1,NODE>();
        } else if ( rp_method == 1 ) {
            grad3d = new Grad3D_ls_so<T1,NODE>();
        } else if ( rp_method == 2 ) {
            grad3d = new Grad3D_ab<T1,NODE>();
        }
        bool reachedTx = false;

        for ( T2 nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == curr_pt ) {
                nodeNo = nn;
                onNode = true;
                s1 = nodes[nodeNo].getNodeSlowness();
                break;
            }
        }
        if ( !onNode ) {
            cellNo = getCellNo( curr_pt );

            std::array<T2,4> itmp = getPrimary(cellNo);
            T2 ind[6][2] = {
                {itmp[0], itmp[1]},
                {itmp[0], itmp[2]},
                {itmp[0], itmp[3]},
                {itmp[1], itmp[2]},
                {itmp[1], itmp[3]},
                {itmp[2], itmp[3]} };

            for ( size_t n=0; n<6; ++n ) {
                if ( areCollinear(curr_pt, nodes[ind[n][0]], nodes[ind[n][1]]) ) {
                    onEdge = true;
                    edgeNodes[0] = ind[n][0];
                    edgeNodes[1] = ind[n][1];

                    if ( processVel )
                        s1 = Interpolator<T1>::linearVel(r_tmp.back(),
                                                         nodes[edgeNodes[0]],
                                                         nodes[edgeNodes[1]]);
                    else
                        s1 = Interpolator<T1>::linear(r_tmp.back(),
                                                      nodes[edgeNodes[0]],
                                                      nodes[edgeNodes[1]]);

                    break;
                }
            }
            if ( !onEdge ) {
                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }

                for ( size_t n=0; n<4; ++n ) {
                    if ( areCoplanar(curr_pt, nodes[ind[n][0]], nodes[ind[n][1]], nodes[ind[n][2]]) ) {
                        onFace = true;
                        faceNodes[0] = ind[n][0];
                        faceNodes[1] = ind[n][1];
                        faceNodes[2] = ind[n][2];

                        if ( processVel )
                            s1 = Interpolator<T1>::bilinearTriangleVel(curr_pt,
                                                                       nodes[faceNodes[0]],
                                                                       nodes[faceNodes[1]],
                                                                       nodes[faceNodes[2]]);
                        else
                            s1 = Interpolator<T1>::bilinearTriangle(curr_pt,
                                                                    nodes[faceNodes[0]],
                                                                    nodes[faceNodes[1]],
                                                                    nodes[faceNodes[2]]);
                        //  faceNodes shoud not be assigned, face was not intersected
                        faceNodes = {0, 0, 0};
                        break;
                    }
                }
                if ( !onFace ) {
                    if ( processVel )
                        s1 = Interpolator<T1>::trilinearTriangleVel(curr_pt,
                                                                    nodes[itmp[0]],
                                                                    nodes[itmp[1]],
                                                                    nodes[itmp[2]],
                                                                    nodes[itmp[3]]);
                    else
                        s1 = Interpolator<T1>::trilinearTriangle(curr_pt,
                                                                 nodes[itmp[0]],
                                                                 nodes[itmp[1]],
                                                                 nodes[itmp[2]],
                                                                 nodes[itmp[3]]);

                }
            }
        }

        for ( size_t nt=0; nt<txCell.size(); ++nt ) {
            if ( cellNo == txCell[nt] ) {
                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( processVel )
                    s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                nodes[itmp[0]],
                                                                nodes[itmp[1]],
                                                                nodes[itmp[2]],
                                                                nodes[itmp[3]]);
                else
                    s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                             nodes[itmp[0]],
                                                             nodes[itmp[1]],
                                                             nodes[itmp[2]],
                                                             nodes[itmp[3]]);

                tt += t0[nt] + 0.5*(s1 + s2) * r_tmp.back().getDistance( Tx[nt] );
                r_tmp.push_back( Tx[nt] );
                reachedTx = true;
                break;
            }
        }

        sxyz<T1> g;
        size_t max_itertions = tetrahedra.size();
        size_t itertions = 0;
        while ( reachedTx == false ) {

            if ( onNode ) {
#ifdef DEBUG_RP
                printRaypathData(curr_pt, g, onNode, onEdge, onFace, cellNo,
                                 nodeNo, edgeNodes, faceNodes);
#endif
                if ( rp_method < 2 ) {
                    // find cells common to edge
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                        getNeighborNodes(*nc, nnodes);
                    }
                    // compute gradient with nodes from all common cells
                    g = grad3d->compute(curr_pt, nodes[nodeNo].getTT(threadNo), nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(1);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[nodeNo]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }

                // find cell for which gradient intersect opposing face
                bool foundIntersection = false;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                    std::array<T2,3> nb;
                    size_t n=0;
                    for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                        if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                            nb[n++] = *nn;
                        }
                    }
                    std::sort(nb.begin(), nb.end());
#ifdef DEBUG_RP
                    if ( verbose > 1 ) {
                        std::cout << "import numpy as np\nimport matplotlib.pyplot as plt\nfrom mpl_toolkits.mplot3d import Axes3D\n\n";
                        std::cout << "fig = plt.figure()\nax = fig.add_subplot(111, projection='3d')\n";
                        std::cout << "pt = np.array([" << nodes[nodeNo].getX() << ", " << nodes[nodeNo].getY() << ", " << nodes[nodeNo].getZ() << "])\n";
                        std::cout << "g = np.array([" << g.x << ", " << g.y << ", " << g.z << "])\n";
                        std::cout << "pt2 = pt + 0.1*g\n";
                        std::cout << "c1 = np.array([" << nodes[nb[0]].getX() << ", " << nodes[nb[0]].getY() << ", " << nodes[nb[0]].getZ() << "])\n";
                        std::cout << "c2 = np.array([" << nodes[nb[1]].getX() << ", " << nodes[nb[1]].getY() << ", " << nodes[nb[1]].getZ() << "])\n";
                        std::cout << "c3 = np.array([" << nodes[nb[2]].getX() << ", " << nodes[nb[2]].getY() << ", " << nodes[nb[2]].getZ() << "])\n";
                        std::cout << "ax.plot([c1[0], c2[0], c3[0], c1[0]], [c1[1], c2[1], c3[1], c1[1]], [c1[2], c2[2], c3[2], c1[2]], 'k')\n";
                        std::cout << "ax.plot([pt[0], pt2[0]], [pt[1], pt2[1]], [pt[2], pt2[2]], 'r')\n";
                        std::cout << "ax.scatter(pt[0], pt[1], pt[2], c='r')\n\n\n";
                    }
#endif
                    foundIntersection = intersectVecTriangle( nodeNo, g, nb[0], nb[1], nb[2], curr_pt);
                    if ( !foundIntersection ) {
                        continue;
                    }

                    //  check if cell is (one of) TxCell(s)
                    for (size_t nt=0; nt<Tx.size(); ++nt) {
                        if ( *nc == txCell[nt] ) {
                            std::array<T2,4> itmp = getPrimary(*nc);

                            if ( processVel )
                                s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                            nodes[itmp[0]],
                                                                            nodes[itmp[1]],
                                                                            nodes[itmp[2]],
                                                                            nodes[itmp[3]]);
                            else
                                s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                         nodes[itmp[0]],
                                                                         nodes[itmp[1]],
                                                                         nodes[itmp[2]],
                                                                         nodes[itmp[3]]);

                            tt += t0[nt] + 0.5*(s1 + s2) * r_tmp.back().getDistance( Tx[nt] );
                            r_tmp.push_back( Tx[nt] );
                            reachedTx = true;
                            break;
                        }
                    }
                    if ( reachedTx ) {
                        break;
                    }

                    bool break_flag = check_pt_location(curr_pt, nb, onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    tt += 0.5*(s1 + s2) * r_tmp.back().getDistance( curr_pt );
                    s1 = s2;
                    r_tmp.push_back( curr_pt );

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell1(faceNodes, nodeNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }

                if ( foundIntersection == false ) {
                    // projet gradient on face and find intersection on next edge
                    sxyz<T1> pt_i;
                    foundIntersection = projectOnFace(curr_pt, nodeNo, g, edgeNodes, pt_i, threadNo);

                    if ( foundIntersection == false ) {
                        std::cout << "\n\nWarning: finding raypath (onNode) failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                        break;
                    }
                    // find which cell we are in
                    for ( size_t nc=0; nc<tetrahedra.size(); ++nc ) {
                        std::array<T2,4> itet = {tetrahedra[nc].i[0],
                            tetrahedra[nc].i[1],
                            tetrahedra[nc].i[2],
                            tetrahedra[nc].i[3]};
                        // because we are at the surface, there is only one cell with nodeNo & edgeNodes
                        if (std::find(itet.begin(), itet.end(), nodeNo) != itet.end() &&
                            std::find(itet.begin(), itet.end(), edgeNodes[0]) != itet.end() &&
                            std::find(itet.begin(), itet.end(), edgeNodes[1]) != itet.end()) {
                            cellNo = static_cast<T2>(nc);
                            break;
                        }
                    }

                    //  check if cell is (one of) TxCell(s)
                    for (size_t nt=0; nt<Tx.size(); ++nt) {
                        if ( cellNo == txCell[nt] ) {
                            std::array<T2,4> itmp = getPrimary(cellNo);

                            if ( processVel )
                                s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                            nodes[itmp[0]],
                                                                            nodes[itmp[1]],
                                                                            nodes[itmp[2]],
                                                                            nodes[itmp[3]]);
                            else
                                s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                         nodes[itmp[0]],
                                                                         nodes[itmp[1]],
                                                                         nodes[itmp[2]],
                                                                         nodes[itmp[3]]);

                            tt += t0[nt] + 0.5*(s1 + s2) * r_tmp.back().getDistance( Tx[nt] );
                            r_tmp.push_back( Tx[nt] );
                            reachedTx = true;
                            break;
                        }
                    }
                    if ( reachedTx ) {
                        break;
                    }

                    // we might be on one of the nodes
                    if ( pt_i.getDistance(nodes[edgeNodes[0]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[0];
                        pt_i = nodes[edgeNodes[0]];
                        onEdge = false;
                    } else if ( pt_i.getDistance(nodes[edgeNodes[1]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[1];
                        pt_i = nodes[edgeNodes[1]];
                        onEdge = false;
                    } else {
                        onEdge = true;
                        onNode = false;
                    }

                    curr_pt = pt_i;

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    tt += 0.5*(s1 + s2) * r_tmp.back().getDistance( curr_pt );
                    s1 = s2;
                    r_tmp.push_back( curr_pt );
                }

            } else if ( onEdge ) {
#ifdef DEBUG_RP
                printRaypathData(curr_pt, g, onNode, onEdge, onFace, cellNo,
                                 nodeNo, edgeNodes, faceNodes);
#endif
                // find cells common to edge
                std::vector<T2> cells;
                for ( auto nc0=nodes[edgeNodes[0]].getOwners().begin(); nc0!=nodes[edgeNodes[0]].getOwners().end(); ++nc0 ) {
                    if ( std::find(nodes[edgeNodes[1]].getOwners().begin(), nodes[edgeNodes[1]].getOwners().end(), *nc0)!=nodes[edgeNodes[1]].getOwners().end() ) {
                        cells.push_back( *nc0 );
                    }
                }
                if ( rp_method < 2 ) {
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    for (size_t n=0; n<cells.size(); ++n ) {
                        getNeighborNodes(cells[n], nnodes);
                    }
                    T1 d01 = nodes[edgeNodes[0]].getDistance(nodes[edgeNodes[1]]);
                    T1 w0 = curr_pt.getDistance(nodes[edgeNodes[1]]) / d01;
                    T1 w1 = curr_pt.getDistance(nodes[edgeNodes[0]]) / d01;
                    T1 curr_t = nodes[edgeNodes[0]].getTT(threadNo)*w0 + nodes[edgeNodes[1]].getTT(threadNo)*w1;
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(2);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[edgeNodes[0]]);
                    ref_pt[1] = &(nodes[edgeNodes[1]]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, edgeNodes, Tx, txCell);

                bool foundIntersection = false;
                for (size_t n=0; n<cells.size(); ++n ) {

                    cellNo = cells[n];

                    // there are 2 faces that might be intersected
                    std::array<T2,2> edgeNodes2;
                    size_t n2=0;
                    for ( auto nn=this->neighbors[cellNo].begin(); nn!= this->neighbors[cellNo].end(); ++nn ) {
                        if ( *nn!=edgeNodes[0] && *nn!=edgeNodes[1] && nodes[*nn].isPrimary() ) {
                            edgeNodes2[n2++] = *nn;
                        }
                    }

                    sxyz<T1> pt_i;
                    T2 itmpNode;
                    foundIntersection = intersectVecTriangle(curr_pt, g,
                                                             edgeNodes[0],
                                                             edgeNodes2[0],
                                                             edgeNodes2[1], pt_i);
                    itmpNode = edgeNodes[0];
                    if ( !foundIntersection ) {
                        foundIntersection = intersectVecTriangle(curr_pt, g,
                                                                 edgeNodes[1],
                                                                 edgeNodes2[0],
                                                                 edgeNodes2[1], pt_i);
                        itmpNode = edgeNodes[1];
                    }
                    if ( !foundIntersection ) {
                        continue;
                    }

                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, this->neighbors[cellNo],
                                                        {itmpNode, edgeNodes2[0], edgeNodes2[1]},
                                                        onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);
                    tt += 0.5*(s1 + s2) * r_tmp.back().getDistance( curr_pt );
                    s1 = s2;
                    r_tmp.push_back( curr_pt );

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }
                if ( foundIntersection == false ) {
                    sxyz<T1> pt_i;
                    g = projectOnFace(curr_pt, g, edgeNodes, cells, pt_i, threadNo);
                    if ( g.x==0.0 && g.y==0.0 && g.z==0.0 ) {
                        foundIntersection = false;
                    } else {
                        foundIntersection = true;
                    }

                    if ( foundIntersection == false || curr_pt == pt_i ) {
                        std::cout << "\n\nWarning: finding raypath (onEdge) failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                        break;
                    }

                    // we might be on one of the nodes
                    if ( pt_i.getDistance(nodes[edgeNodes[0]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[0];
                        pt_i = nodes[edgeNodes[0]];
                        onEdge = false;
                    } else if ( pt_i.getDistance(nodes[edgeNodes[1]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[1];
                        pt_i = nodes[edgeNodes[1]];
                        onEdge = false;
                    }

                    curr_pt = pt_i;

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    tt += 0.5*(s1 + s2) * r_tmp.back().getDistance( curr_pt );
                    s1 = s2;
                    r_tmp.push_back( curr_pt );
                }

            } else if ( onFace ) { // on Face
#ifdef DEBUG_RP
                printRaypathData(curr_pt, g, onNode, onEdge, onFace, cellNo,
                                 nodeNo, edgeNodes, faceNodes);
#endif
                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( rp_method < 2 ) {
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    getNeighborNodes(cellNo, nnodes);
                    T1 curr_t;
                    if ( r_tmp.size() <= 1 ) {
                        curr_t = Interpolator<T1>::trilinearTime(curr_pt,
                                                                 nodes[itmp[0]],
                                                                 nodes[itmp[1]],
                                                                 nodes[itmp[2]],
                                                                 nodes[itmp[3]],
                                                                 threadNo);
                    } else {
                        curr_t = Interpolator<T1>::bilinearTime(curr_pt,
                                                                nodes[faceNodes[0]],
                                                                nodes[faceNodes[1]],
                                                                nodes[faceNodes[2]],
                                                                threadNo);
                    }
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(3);
                    if ( r_tmp.size() <= 1 ) {
                        ref_pt[0] = &(nodes[itmp[0]]);
                        ref_pt[1] = &(nodes[itmp[1]]);
                        ref_pt[2] = &(nodes[itmp[2]]);
                        ref_pt.push_back( &(nodes[itmp[3]]) );
                    } else {
                        ref_pt[0] = &(nodes[faceNodes[0]]);
                        ref_pt[1] = &(nodes[faceNodes[1]]);
                        ref_pt[2] = &(nodes[faceNodes[2]]);
                    }
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }
                // there are 3 faces that might be intersected

                bool foundIntersection = false;
                for ( size_t n=0; n<4; ++n ) {
                    if ( ind[n] == faceNodes ) continue;

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);

                    if ( !foundIntersection )
                        continue;

                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    tt += 0.5*(s1 + s2) * r_tmp.back().getDistance( curr_pt );
                    s1 = s2;
                    r_tmp.push_back( curr_pt );

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }

                if ( foundIntersection == false ) {

                    // we must be on an face with gradient pointing slightly outward tetrahedron
                    // return in other cell but keep gradient
                    cellNo = findAdjacentCell2(faceNodes, cellNo);

                    std::array<T2,4> itmp = getPrimary(cellNo);
                    ind[0] = { { itmp[0], itmp[1], itmp[2] } };
                    ind[1] = { { itmp[0], itmp[1], itmp[3] } };
                    ind[2] = { { itmp[0], itmp[2], itmp[3] } };
                    ind[3] = { { itmp[1], itmp[2], itmp[3] } };

                    for ( size_t n=0; n<4; ++n ) {
                        std::sort( ind[n].begin(), ind[n].end() );
                    }

                    for ( size_t n=0; n<4; ++n ) {
                        if ( ind[n] == faceNodes ) continue;

                        sxyz<T1> pt_i;
                        foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                                 ind[n][1], ind[n][2],
                                                                 pt_i);

                        if ( !foundIntersection ) {
                            continue;
                        }
                        curr_pt = pt_i;

                        bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                            nodeNo, onEdge, edgeNodes,
                                                            onFace, faceNodes);

                        s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                             faceNodes);

                        tt += 0.5*(s1 + s2) * r_tmp.back().getDistance( curr_pt );
                        s1 = s2;
                        r_tmp.push_back( curr_pt );

                        if ( break_flag ) break;

                        // find next cell
                        cellNo = findAdjacentCell2(faceNodes, cellNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                            << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                            tt = 0.0;
                            r_tmp.resize(1);
                            r_tmp[0] = Rx;
                            reachedTx = true;
                        }
                        break;
                    }
                }
                if ( foundIntersection == false ) {
                    g = projectOnFace(g, faceNodes);

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecEdge(curr_pt, g, faceNodes, pt_i, edgeNodes);

                    if ( foundIntersection == false ) {
                        std::cout << "\n\nWarning: finding raypath (onFace) failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                        break;
                    }

                    // we might be on one of the nodes
                    if ( pt_i.getDistance(nodes[edgeNodes[0]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[0];
                        pt_i = nodes[edgeNodes[0]];
                        onEdge = false;
                        onFace = false;
                    } else if ( pt_i.getDistance(nodes[edgeNodes[1]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[1];
                        pt_i = nodes[edgeNodes[1]];
                        onEdge = false;
                        onFace = false;
                    } else {
                        onEdge = true;
                        onFace = false;
                    }
                    curr_pt = pt_i;

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    tt += 0.5*(s1 + s2) * r_tmp.back().getDistance( curr_pt );
                    s1 = s2;
                    r_tmp.push_back( curr_pt );
                }
            } else { // at Rx, somewhere in a tetrahedron
#ifdef DEBUG_RP
                printRaypathData(curr_pt, g, onNode, onEdge, onFace, cellNo,
                                 nodeNo, edgeNodes, faceNodes);
#endif
                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( rp_method < 2 ) {
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    getNeighborNodes(cellNo, nnodes);
                    T1 curr_t = Interpolator<T1>::trilinearTime(curr_pt,
                                                                nodes[itmp[0]],
                                                                nodes[itmp[1]],
                                                                nodes[itmp[2]],
                                                                nodes[itmp[3]],
                                                                threadNo);
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(4);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[itmp[0]]);
                    ref_pt[1] = &(nodes[itmp[1]]);
                    ref_pt[2] = &(nodes[itmp[2]]);
                    ref_pt[3] = &(nodes[itmp[3]]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }
                // there are 4 faces that might be intersected

                bool foundIntersection = false;
                for ( size_t n=0; n<4; ++n ) {

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);

                    if ( !foundIntersection )
                        continue;

                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    tt += 0.5*(s1 + s2) * r_tmp.back().getDistance( curr_pt );
                    s1 = s2;
                    r_tmp.push_back( curr_pt );

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        reachedTx = true;
                    }
                    break;
                }
                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath within cell failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    tt = 0.0;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    reachedTx = true;
                }
            }

            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        if ( curr_pt.getDistance( Tx[nt] ) < minDist ) {
                            tt += t0[nt];
                            reachedTx = true;
                            break;
                        }
                    } else if ( txOnEdge[nt] ) {
                        if ( curr_pt.getDistance(nodes[txEdges[nt][0]]) < minDist ||
                            curr_pt.getDistance(nodes[txEdges[nt][1]]) < minDist ) {

                            if ( processVel )
                                s2 = Interpolator<T1>::linearVel(Tx[nt],
                                                                 nodes[txEdges[nt][0]],
                                                                 nodes[txEdges[nt][1]]);
                            else
                                s2 = Interpolator<T1>::linear(Tx[nt],
                                                              nodes[txEdges[nt][0]],
                                                              nodes[txEdges[nt][1]]);

                            tt += t0[nt] + 0.5*(s1 + s2) * r_tmp.back().getDistance( Tx[nt] );
                            reachedTx = true;
                            r_tmp.push_back(Tx[nt]);
                            break;
                        }
                    } else if ( txOnFace[nt] ) {
                        if ( curr_pt.getDistance(nodes[txFaces[nt][0]]) < minDist ||
                            curr_pt.getDistance(nodes[txFaces[nt][1]]) < minDist ||
                            curr_pt.getDistance(nodes[txFaces[nt][2]]) < minDist ) {

                            if ( processVel )
                                s2 = Interpolator<T1>::bilinearTriangleVel(Tx[nt],
                                                                           nodes[txFaces[nt][0]],
                                                                           nodes[txFaces[nt][1]],
                                                                           nodes[txFaces[nt][2]]);
                            else
                                s2 = Interpolator<T1>::bilinearTriangle(Tx[nt],
                                                                        nodes[txFaces[nt][0]],
                                                                        nodes[txFaces[nt][1]],
                                                                        nodes[txFaces[nt][2]]);

                            tt += t0[nt] + 0.5*(s1 + s2) * r_tmp.back().getDistance( Tx[nt] );
                            reachedTx = true;
                            r_tmp.push_back(Tx[nt]);
                            break;
                        }
                    }
                }
            } if ( onEdge ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        if ( txNode[nt] == edgeNodes[0] || txNode[nt] == edgeNodes[1] ) {
                            s2 = nodes[txNode[nt]].getNodeSlowness();
                            tt += t0[nt] + 0.5*(s1 + s2) * r_tmp.back().getDistance( Tx[nt] );
                            r_tmp.push_back( Tx[nt] );
                            reachedTx = true;
                            break;
                        }
                    } else {
                        std::array<T2,4> itmp = getPrimary(txCell[nt]);
                        // find adjacent cells
                        const T2 ind[6][2] = {
                            {itmp[0], itmp[1]},
                            {itmp[0], itmp[2]},
                            {itmp[0], itmp[3]},
                            {itmp[1], itmp[2]},
                            {itmp[1], itmp[3]},
                            {itmp[2], itmp[3]} };
                        for ( size_t ne=0; ne<6; ++ne ) {
                            if ( (ind[ne][0] == edgeNodes[0] && ind[ne][1] == edgeNodes[1]) ||
                                (ind[ne][0] == edgeNodes[1] && ind[ne][1] == edgeNodes[0]) ) {
                                if ( processVel )
                                    s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                nodes[itmp[0]],
                                                                                nodes[itmp[1]],
                                                                                nodes[itmp[2]],
                                                                                nodes[itmp[3]]);
                                else
                                    s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                             nodes[itmp[0]],
                                                                             nodes[itmp[1]],
                                                                             nodes[itmp[2]],
                                                                             nodes[itmp[3]]);

                                tt += t0[nt] + 0.5*(s1 + s2) * r_tmp.back().getDistance( Tx[nt] );
                                r_tmp.push_back( Tx[nt] );
                                reachedTx = true;
                                break;
                            }
                        }
                    }
                }
            } else {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        for ( auto nc=nodes[txNode[nt]].getOwners().begin(); nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                            if ( cellNo == *nc ) {
                                std::array<T2,4> itmp = getPrimary(cellNo);

                                if ( processVel )
                                    s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                nodes[itmp[0]],
                                                                                nodes[itmp[1]],
                                                                                nodes[itmp[2]],
                                                                                nodes[itmp[3]]);
                                else
                                    s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                             nodes[itmp[0]],
                                                                             nodes[itmp[1]],
                                                                             nodes[itmp[2]],
                                                                             nodes[itmp[3]]);

                                tt += t0[nt] + 0.5*(s1 + s2) * r_tmp.back().getDistance( Tx[nt] );
                                r_tmp.push_back( Tx[nt] );
                                reachedTx = true;
                                break;
                            }
                        }
                    } else {
                        if ( cellNo == txCell[nt] ) {
                            std::array<T2,4> itmp = getPrimary(cellNo);
                            if ( processVel )
                                s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                            nodes[itmp[0]],
                                                                            nodes[itmp[1]],
                                                                            nodes[itmp[2]],
                                                                            nodes[itmp[3]]);
                            else
                                s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                         nodes[itmp[0]],
                                                                         nodes[itmp[1]],
                                                                         nodes[itmp[2]],
                                                                         nodes[itmp[3]]);

                            tt += t0[nt] + 0.5*(s1 + s2) * r_tmp.back().getDistance( Tx[nt] );
                            r_tmp.push_back( Tx[nt] );
                            reachedTx = true;
                            break;
                        } else {
                            for ( size_t nn=0; nn<txNeighborCells[nt].size(); ++nn ) {
                                if ( cellNo == txNeighborCells[nt][nn] ) {
                                    std::array<T2,4> itmp = getPrimary(txCell[nt]);
                                    std::array<T2,3> ind[4] = {
                                        { { itmp[0], itmp[1], itmp[2] } },
                                        { { itmp[0], itmp[1], itmp[3] } },
                                        { { itmp[0], itmp[2], itmp[3] } },
                                        { { itmp[1], itmp[2], itmp[3] } }
                                    };
                                    bool found = false;
                                    for ( size_t n=0; n<4; ++n ) {
                                        std::sort( ind[n].begin(), ind[n].end() );
                                        if ( faceNodes == ind[n] ) {
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found ) {
                                        itmp = getPrimary(cellNo);
                                        if ( processVel )
                                            s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                        nodes[itmp[0]],
                                                                                        nodes[itmp[1]],
                                                                                        nodes[itmp[2]],
                                                                                        nodes[itmp[3]]);
                                        else
                                            s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                                     nodes[itmp[0]],
                                                                                     nodes[itmp[1]],
                                                                                     nodes[itmp[2]],
                                                                                     nodes[itmp[3]]);

                                        tt += t0[nt] + 0.5*(s1 + s2) * r_tmp.back().getDistance( Tx[nt] );
                                        r_tmp.push_back( Tx[nt] );
                                        reachedTx = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if ( reachedTx ) break;
                }
            }
            itertions ++;
            if (itertions > max_itertions){
                std::cout << "\n\nWarning: raypath failed to converge for Rx : "
                << Rx.x << ' ' << Rx.y << ' ' << Rx.z << " (infinite loop risk)"<<std::endl;
                tt = 0.0;
                r_tmp.resize(1);
                r_tmp[0] = Rx;
                break; // break while loop
            }
        }
        // for inversion, the order should be from Tx to Rx, so we reorder...
        size_t npts = r_tmp.size();
        r_data.resize( npts );
        for ( size_t nn=0; nn<npts; ++nn ) {
            r_data[nn] = r_tmp[ npts-1-nn ];
        }

        delete grad3d;
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                          const std::vector<T1>& t0,
                                          const sxyz<T1> &Rx,
                                          std::vector<sijv<T1>>& m_data,
                                          const size_t RxNo,
                                          T1 &tt,
                                          const size_t threadNo) const {
        // Important note: m_data holds terms of traveltime derivative w/r to slowness (not velocity)
        T1 minDist = small;
        tt = 0.0;
        T1 s1=0.0, s2=0.0;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        const txInfo_t& txi = getTxInfo(Tx, threadNo);
        const std::vector<bool>& txOnNode = txi.txOnNode;
        const std::vector<bool>& txOnEdge = txi.txOnEdge;
        const std::vector<bool>& txOnFace = txi.txOnFace;
        const std::vector<T2>& txNode = txi.txNode;
        const std::vector<T2>& txCell = txi.txCell;
        const std::vector<std::array<T2,2>>& txEdges = txi.txEdges;
        const std::vector<std::array<T2,3>>& txFaces = txi.txFaces;
        const std::vector<std::vector<T2>>& txNeighborCells = txi.txNeighborCells;

        T2 cellNo=0, nodeNo=0, nodeNoPrev=0;
        sxyz<T1> curr_pt( Rx ), mid_pt, prev_pt( Rx );
        bool atRx = true;
        sijv<T1> m;
        m.i = RxNo;

        bool onNode = false;
        bool onEdge = false;
        bool onFace = false;
        std::array<T2,2> edgeNodes{ {0, 0} }, edgeNodesPrev;
        std::array<T2,3> faceNodes{ {0, 0, 0} }, faceNodesPrev;
        Grad3D<T1,NODE>* grad3d = nullptr;
        if ( rp_method == 0 ) {
            grad3d = new Grad3D_ls_fo<T1,NODE>();
        } else if ( rp_method == 1 ) {
            grad3d = new Grad3D_ls_so<T1,NODE>();
        } else if ( rp_method == 2 ) {
            grad3d = new Grad3D_ab<T1,NODE>();
        }
        bool reachedTx = false;

        for ( T2 nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == curr_pt ) {
                nodeNo = nn;
                onNode = true;
                s1 = nodes[nodeNo].getNodeSlowness();
                break;
            }
        }
        if ( !onNode ) {
            cellNo = getCellNo( curr_pt );

            std::array<T2,4> itmp = getPrimary(cellNo);
            T2 ind[6][2] = {
                {itmp[0], itmp[1]},
                {itmp[0], itmp[2]},
                {itmp[0], itmp[3]},
                {itmp[1], itmp[2]},
                {itmp[1], itmp[3]},
                {itmp[2], itmp[3]} };

            for ( size_t n=0; n<6; ++n ) {
                if ( areCollinear(curr_pt, nodes[ind[n][0]], nodes[ind[n][1]]) ) {
                    onEdge = true;
                    edgeNodes[0] = ind[n][0];
                    edgeNodes[1] = ind[n][1];

                    if ( processVel )
                        s1 = Interpolator<T1>::linearVel(curr_pt,
                                                         nodes[edgeNodes[0]],
                                                         nodes[edgeNodes[1]]);
                    else
                        s1 = Interpolator<T1>::linear(curr_pt,
                                                      nodes[edgeNodes[0]],
                                                      nodes[edgeNodes[1]]);

                    break;
                }
            }
            if ( !onEdge ) {
                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }

                for ( size_t n=0; n<4; ++n ) {
                    if ( areCoplanar(curr_pt, nodes[ind[n][0]], nodes[ind[n][1]], nodes[ind[n][2]]) ) {
                        onFace = true;
                        faceNodes[0] = ind[n][0];
                        faceNodes[1] = ind[n][1];
                        faceNodes[2] = ind[n][2];

                        if ( processVel )
                            s1 = Interpolator<T1>::bilinearTriangleVel(curr_pt,
                                                                       nodes[faceNodes[0]],
                                                                       nodes[faceNodes[1]],
                                                                       nodes[faceNodes[2]]);
                        else {
                            s1 = Interpolator<T1>::bilinearTriangle(curr_pt,
                                                                    nodes[faceNodes[0]],
                                                                    nodes[faceNodes[1]],
                                                                    nodes[faceNodes[2]]);
                        }
                        //  faceNodes shoud not be assigned, face was not intersected
                        faceNodes = {0, 0, 0};
                        break;
                    }
                }
                if ( !onFace ) {
                    if ( processVel )
                        s1 = Interpolator<T1>::trilinearTriangleVel(curr_pt,
                                                                    nodes[itmp[0]],
                                                                    nodes[itmp[1]],
                                                                    nodes[itmp[2]],
                                                                    nodes[itmp[3]]);
                    else
                        s1 = Interpolator<T1>::trilinearTriangle(curr_pt,
                                                                 nodes[itmp[0]],
                                                                 nodes[itmp[1]],
                                                                 nodes[itmp[2]],
                                                                 nodes[itmp[3]]);
                }
            }
        }

        for ( size_t nt=0; nt<txCell.size(); ++nt ) {
            if ( cellNo == txCell[nt] ) {
                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( processVel )
                    s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                nodes[itmp[0]],
                                                                nodes[itmp[1]],
                                                                nodes[itmp[2]],
                                                                nodes[itmp[3]]);
                else
                    s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                             nodes[itmp[0]],
                                                             nodes[itmp[1]],
                                                             nodes[itmp[2]],
                                                             nodes[itmp[3]]);

                mid_pt = static_cast<T1>(0.5)*(curr_pt + Tx[nt]);
                T1 ds = curr_pt.getDistance( Tx[nt] );

                std::set<T2> allNodes;
                allNodes.insert( itmp[0] );
                allNodes.insert( itmp[1] );
                allNodes.insert( itmp[2] );
                allNodes.insert( itmp[3] );
                if ( processVel ) {
                    T1 s_sq = computeSlowness(mid_pt, true);
                    s_sq *= s_sq;  // slowness squared
                    update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                } else {
                    update_m_data(m_data, m, allNodes, mid_pt, ds);
                }

                tt += t0[nt] + 0.5*(s1 + s2) * ds;
                reachedTx = true;
                break;
            }
        }
        T1 ds;
        sxyz<T1> g;
        size_t max_itertions = tetrahedra.size();
        size_t itertions = 0;
        while ( reachedTx == false ) {
            if ( onNode ) {

                // find cells common to edge
                std::byte nnodes_buf[8192];
                std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                std::pmr::set<NODE*> nnodes{&nnodes_pool};
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                    getNeighborNodes(*nc, nnodes);
                }

                // compute gradient with nodes from all common cells
                if ( rp_method < 2 ) {
                    // find cells common to edge
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                        getNeighborNodes(*nc, nnodes);
                    }
                    // compute gradient with nodes from all common cells
                    g = grad3d->compute(curr_pt, nodes[nodeNo].getTT(threadNo), nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(1);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[nodeNo]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }

                nodeNoPrev = nodeNo;

                // find cell for which gradient intersect opposing face
                bool foundIntersection = false;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                    std::array<T2,3> nb;
                    size_t n=0;
                    for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                        if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                            nb[n++] = *nn;
                        }
                    }
                    std::sort(nb.begin(), nb.end());

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle( nodeNo, g, nb[0], nb[1], nb[2], pt_i);
                    if ( !foundIntersection ) {
                        continue;
                    }

                    //  check if cell is (one of) TxCell(s)
                    for (size_t nt=0; nt<Tx.size(); ++nt) {
                        if ( *nc == txCell[nt] ) {
                            std::array<T2,4> itmp = getPrimary(*nc);

                            if ( processVel )
                                s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                            nodes[itmp[0]],
                                                                            nodes[itmp[1]],
                                                                            nodes[itmp[2]],
                                                                            nodes[itmp[3]]);
                            else
                                s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                         nodes[itmp[0]],
                                                                         nodes[itmp[1]],
                                                                         nodes[itmp[2]],
                                                                         nodes[itmp[3]]);

                            tt += t0[nt] + 0.5*(s1 + s2) * curr_pt.getDistance( Tx[nt] );

                            prev_pt = curr_pt;
                            curr_pt = Tx[nt];
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            ds = curr_pt.getDistance( prev_pt );
                            std::set<T2> allNodes;
                            allNodes.insert(itmp[0]);
                            allNodes.insert(itmp[1]);
                            allNodes.insert(itmp[2]);
                            allNodes.insert(itmp[3]);

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }

                            reachedTx = true;
                            break;
                        }
                    }
                    if ( reachedTx ) {
                        break;
                    }

                    prev_pt = curr_pt;
                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, nb, onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    // compute terms of matrix M
                    mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                    ds = curr_pt.getDistance( prev_pt );

                    tt += 0.5*(s1 + s2) * ds;
                    s1 = s2;

                    std::set<T2> allNodes;
                    allNodes.insert(nodeNoPrev);
                    if ( onNode ) {
                        allNodes.insert(nodeNo);
                    } else if ( onEdge ) {
                        allNodes.insert( edgeNodes[0] );
                        allNodes.insert( edgeNodes[1] );
                    } else { // onFace
                        allNodes.insert( faceNodes[0] );
                        allNodes.insert( faceNodes[1] );
                        allNodes.insert( faceNodes[2] );
                    }
                    if ( processVel ) {
                        T1 s_sq = computeSlowness(mid_pt, true);
                        s_sq *= s_sq;  // slowness squared
                        update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                    } else {
                        update_m_data(m_data, m, allNodes, mid_pt, ds);
                    }

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell1(faceNodes, nodeNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        m_data.resize(0);
                        reachedTx = true;
                    }
                    break;
                }

                if ( foundIntersection == false ) {
                    nodeNoPrev = nodeNo;
                    // projet gradient on face and find intersection on next edge
                    sxyz<T1> pt_i;
                    foundIntersection = projectOnFace(curr_pt, nodeNo, g, edgeNodes, pt_i, threadNo);

                    if ( foundIntersection == false ) {
                        std::cout << "\n\nWarning: finding raypath (onNode) failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        m_data.resize(0);
                        reachedTx = true;
                        break;
                    }
                    // find which cell we are in
                    for ( size_t nc=0; nc<tetrahedra.size(); ++nc ) {
                        std::array<T2,4> itet = {tetrahedra[nc].i[0],
                            tetrahedra[nc].i[1],
                            tetrahedra[nc].i[2],
                            tetrahedra[nc].i[3]};
                        // because we are at the surface, there is only one cell with nodeNo & edgeNodes
                        if (std::find(itet.begin(), itet.end(), nodeNo) != itet.end() &&
                            std::find(itet.begin(), itet.end(), edgeNodes[0]) != itet.end() &&
                            std::find(itet.begin(), itet.end(), edgeNodes[1]) != itet.end()) {
                            cellNo = static_cast<T2>(nc);
                            break;
                        }
                    }

                    //  check if cell is (one of) TxCell(s)
                    for (size_t nt=0; nt<Tx.size(); ++nt) {
                        if ( cellNo == txCell[nt] ) {
                            std::array<T2,4> itmp = getPrimary(cellNo);

                            if ( processVel )
                                s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                            nodes[itmp[0]],
                                                                            nodes[itmp[1]],
                                                                            nodes[itmp[2]],
                                                                            nodes[itmp[3]]);
                            else
                                s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                         nodes[itmp[0]],
                                                                         nodes[itmp[1]],
                                                                         nodes[itmp[2]],
                                                                         nodes[itmp[3]]);

                            mid_pt = static_cast<T1>(0.5)*(curr_pt + Tx[nt]);
                            ds = curr_pt.getDistance( Tx[nt] );

                            tt += t0[nt] + 0.5*(s1 + s2) * ds;

                            std::set<T2> allNodes;
                            allNodes.insert( itmp[0] );
                            allNodes.insert( itmp[1] );
                            allNodes.insert( itmp[2] );
                            allNodes.insert( itmp[3] );
                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }

                            reachedTx = true;
                            break;
                        }
                    }
                    if ( reachedTx ) {
                        break;
                    }

                    // we might be on one of the nodes
                    if ( pt_i.getDistance(nodes[edgeNodes[0]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[0];
                        pt_i = nodes[edgeNodes[0]];
                        onEdge = false;
                    } else if ( pt_i.getDistance(nodes[edgeNodes[1]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[1];
                        pt_i = nodes[edgeNodes[1]];
                        onEdge = false;
                    } else {
                        onEdge = true;
                        onNode = false;
                    }

                    prev_pt = curr_pt;
                    curr_pt = pt_i;

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                    ds = curr_pt.getDistance( prev_pt );
                    tt += 0.5*(s1 + s2) * ds;
                    s1 = s2;

                    std::set<T2> allNodes;
                    allNodes.insert( nodeNoPrev );
                    allNodes.insert( edgeNodes[0] );
                    allNodes.insert( edgeNodes[1] );
                    if ( processVel ) {
                        T1 s_sq = computeSlowness(mid_pt, true);
                        s_sq *= s_sq;  // slowness squared
                        update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                    } else {
                        update_m_data(m_data, m, allNodes, mid_pt, ds);
                    }
                }

            } else if ( onEdge ) {

                // find cells common to edge
                std::vector<T2> cells;
                std::byte nnodes_buf[8192];
                std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                std::pmr::set<NODE*> nnodes{&nnodes_pool};
                for ( auto nc0=nodes[edgeNodes[0]].getOwners().begin(); nc0!=nodes[edgeNodes[0]].getOwners().end(); ++nc0 ) {
                    if ( std::find(nodes[edgeNodes[1]].getOwners().begin(), nodes[edgeNodes[1]].getOwners().end(), *nc0)!=nodes[edgeNodes[1]].getOwners().end() ) {
                        cells.push_back( *nc0 );
                        getNeighborNodes(*nc0, nnodes);
                    }
                }
                if ( rp_method < 2 ) {
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    for (size_t n=0; n<cells.size(); ++n ) {
                        getNeighborNodes(cells[n], nnodes);
                    }
                    T1 d01 = nodes[edgeNodes[0]].getDistance(nodes[edgeNodes[1]]);
                    T1 w0 = curr_pt.getDistance(nodes[edgeNodes[1]]) / d01;
                    T1 w1 = curr_pt.getDistance(nodes[edgeNodes[0]]) / d01;
                    T1 curr_t = nodes[edgeNodes[0]].getTT(threadNo)*w0 + nodes[edgeNodes[1]].getTT(threadNo)*w1;
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(2);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[edgeNodes[0]]);
                    ref_pt[1] = &(nodes[edgeNodes[1]]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, edgeNodes, Tx, txCell);

                edgeNodesPrev = edgeNodes;

                bool foundIntersection = false;
                for (size_t n=0; n<cells.size(); ++n ) {

                    T2 testCellNo = cells[n];

                    // there are 2 faces that might be intersected
                    std::array<T2,2> edgeNodes2;
                    size_t n2=0;
                    for ( auto nn=this->neighbors[testCellNo].begin(); nn!= this->neighbors[testCellNo].end(); ++nn ) {
                        if ( *nn!=edgeNodes[0] && *nn!=edgeNodes[1] && nodes[*nn].isPrimary() ) {
                            edgeNodes2[n2++] = *nn;
                        }
                    }

                    sxyz<T1> pt_i;
                    T2 itmpNode;
                    foundIntersection = intersectVecTriangle(curr_pt, g,
                                                             edgeNodes[0],
                                                             edgeNodes2[0],
                                                             edgeNodes2[1], pt_i);
                    itmpNode = edgeNodes[0];
                    if ( !foundIntersection ) {
                        foundIntersection = intersectVecTriangle(curr_pt, g,
                                                                 edgeNodes[1],
                                                                 edgeNodes2[0],
                                                                 edgeNodes2[1], pt_i);
                        itmpNode = edgeNodes[1];
                    }
                    if ( !foundIntersection ) {
                        continue;
                    }

                    cellNo = testCellNo;
                    prev_pt = curr_pt;
                    curr_pt = pt_i;
                    bool break_flag = check_pt_location(curr_pt, this->neighbors[cellNo],
                                                        {itmpNode, edgeNodes2[0], edgeNodes2[1]},
                                                        onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    // compute terms of matrix M
                    mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                    ds = curr_pt.getDistance( prev_pt );
                    tt += 0.5*(s1 + s2) * ds;
                    s1 = s2;

                    std::set<T2> allNodes;
                    allNodes.insert( edgeNodesPrev[0] );
                    allNodes.insert( edgeNodesPrev[1] );
                    if ( onNode ) {
                        allNodes.insert(nodeNo);
                    } else if ( onEdge ) {
                        allNodes.insert( edgeNodes[0] );
                        allNodes.insert( edgeNodes[1] );
                    } else { // onFace
                        allNodes.insert( faceNodes[0] );
                        allNodes.insert( faceNodes[1] );
                        allNodes.insert( faceNodes[2] );
                    }
                    if ( processVel ) {
                        T1 s_sq = computeSlowness(mid_pt, true);
                        s_sq *= s_sq;  // slowness squared
                        update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                    } else {
                        update_m_data(m_data, m, allNodes, mid_pt, ds);
                    }

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        m_data.resize(0);
                        reachedTx = true;
                    }
                    break;
                }
                if ( foundIntersection == false ) {

                    edgeNodesPrev = edgeNodes;
                    sxyz<T1> pt_i;
                    g = projectOnFace(curr_pt, g, edgeNodes, cells, pt_i, threadNo);
                    if ( g.x==0.0 && g.y==0.0 && g.z==0.0 ) {
                        foundIntersection = false;
                    } else {
                        foundIntersection = true;
                    }

                    if ( foundIntersection == false || curr_pt == pt_i ) {
                        std::cout << "\n\nWarning: finding raypath (onEdge) failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        m_data.resize(0);
                        reachedTx = true;
                        break;
                    }

                    // we might be on one of the nodes
                    if ( pt_i.getDistance(nodes[edgeNodes[0]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[0];
                        pt_i = nodes[edgeNodes[0]];
                        onEdge = false;
                    } else if ( pt_i.getDistance(nodes[edgeNodes[1]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[1];
                        pt_i = nodes[edgeNodes[1]];
                        onEdge = false;
                    }

                    prev_pt = curr_pt;
                    curr_pt = pt_i;

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                    ds = curr_pt.getDistance( prev_pt );
                    tt += 0.5*(s1 + s2) * ds;
                    s1 = s2;

                    std::set<T2> allNodes;
                    allNodes.insert( edgeNodesPrev[0] );
                    allNodes.insert( edgeNodesPrev[1] );
                    allNodes.insert( edgeNodes[0] );
                    allNodes.insert( edgeNodes[1] );

                    if ( processVel ) {
                        T1 s_sq = computeSlowness(mid_pt, true);
                        s_sq *= s_sq;  // slowness squared
                        update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                    } else {
                        update_m_data(m_data, m, allNodes, mid_pt, ds);
                    }
                }

            } else if ( onFace ) { // on Face

                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( rp_method < 2 ) {
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    getNeighborNodes(cellNo, nnodes);
                    T1 curr_t;
                    if ( atRx ) {
                        curr_t = Interpolator<T1>::trilinearTime(curr_pt,
                                                                 nodes[itmp[0]],
                                                                 nodes[itmp[1]],
                                                                 nodes[itmp[2]],
                                                                 nodes[itmp[3]],
                                                                 threadNo);
                        atRx = false;
                    } else {
                        curr_t = Interpolator<T1>::bilinearTime(curr_pt,
                                                                nodes[faceNodes[0]],
                                                                nodes[faceNodes[1]],
                                                                nodes[faceNodes[2]],
                                                                threadNo);
                    }
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(3);
                    if ( atRx ) {
                        ref_pt[0] = &(nodes[itmp[0]]);
                        ref_pt[1] = &(nodes[itmp[1]]);
                        ref_pt[2] = &(nodes[itmp[2]]);
                        ref_pt.push_back( &(nodes[itmp[3]]) );
                        atRx = false;
                    } else {
                        ref_pt[0] = &(nodes[faceNodes[0]]);
                        ref_pt[1] = &(nodes[faceNodes[1]]);
                        ref_pt[2] = &(nodes[faceNodes[2]]);
                    }
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }
                // there are 3 faces that might be intersected

                faceNodesPrev = faceNodes;

                bool foundIntersection = false;
                for ( size_t n=0; n<4; ++n ) {
                    if ( ind[n] == faceNodes ) continue;

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);

                    if ( !foundIntersection )
                        continue;

                    prev_pt = curr_pt;
                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    // compute terms of matrix M
                    mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                    ds = curr_pt.getDistance( prev_pt );

                    tt += 0.5*(s1 + s2) * ds;
                    s1 = s2;

                    std::set<T2> allNodes;
                    allNodes.insert( faceNodesPrev[0] );
                    allNodes.insert( faceNodesPrev[1] );
                    allNodes.insert( faceNodesPrev[2] );
                    if ( onNode ) {
                        allNodes.insert(nodeNo);
                    } else if ( onEdge ) {
                        allNodes.insert( edgeNodes[0] );
                        allNodes.insert( edgeNodes[1] );
                    } else { // onFace
                        allNodes.insert( faceNodes[0] );
                        allNodes.insert( faceNodes[1] );
                        allNodes.insert( faceNodes[2] );
                    }
                    if ( processVel ) {
                        T1 s_sq = computeSlowness(mid_pt, true);
                        s_sq *= s_sq;  // slowness squared
                        update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                    } else {
                        update_m_data(m_data, m, allNodes, mid_pt, ds);
                    }

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        m_data.resize(0);
                        reachedTx = true;
                    }
                    break;
                }

                if ( foundIntersection == false ) {

                    // we must be on an face with gradient pointing slightly outward tetrahedron
                    // return in other cell but keep gradient
                    cellNo = findAdjacentCell2(faceNodes, cellNo);

                    std::array<T2,4> itmp = getPrimary(cellNo);
                    ind[0] = { { itmp[0], itmp[1], itmp[2] } };
                    ind[1] = { { itmp[0], itmp[1], itmp[3] } };
                    ind[2] = { { itmp[0], itmp[2], itmp[3] } };
                    ind[3] = { { itmp[1], itmp[2], itmp[3] } };

                    for ( size_t n=0; n<4; ++n ) {
                        std::sort( ind[n].begin(), ind[n].end() );
                    }

                    for ( size_t n=0; n<4; ++n ) {
                        if ( ind[n] == faceNodes ) continue;

                        sxyz<T1> pt_i;
                        foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                                 ind[n][1], ind[n][2],
                                                                 pt_i);

                        if ( !foundIntersection ) {
                            continue;
                        }
                        prev_pt = curr_pt;
                        curr_pt = pt_i;

                        bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                            nodeNo, onEdge, edgeNodes,
                                                            onFace, faceNodes);

                        s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                             faceNodes);

                        mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                        ds = curr_pt.getDistance( prev_pt );

                        tt += 0.5*(s1 + s2) * ds;
                        s1 = s2;

                        std::set<T2> allNodes;
                        allNodes.insert( faceNodesPrev[0] );
                        allNodes.insert( faceNodesPrev[1] );
                        allNodes.insert( faceNodesPrev[2] );
                        if ( onNode ) {
                            allNodes.insert(nodeNo);
                        } else if ( onEdge ) {
                            allNodes.insert( edgeNodes[0] );
                            allNodes.insert( edgeNodes[1] );
                        } else { // onFace
                            allNodes.insert( faceNodes[0] );
                            allNodes.insert( faceNodes[1] );
                            allNodes.insert( faceNodes[2] );
                        }
                        if ( processVel ) {
                            T1 s_sq = computeSlowness(mid_pt, true);
                            s_sq *= s_sq;  // slowness squared
                            update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                        } else {
                            update_m_data(m_data, m, allNodes, mid_pt, ds);
                        }

                        if ( break_flag ) break;

                        // find next cell
                        cellNo = findAdjacentCell2(faceNodes, cellNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                            << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                            tt = 0.0;
                            m_data.resize(0);
                            reachedTx = true;
                        }
                        break;
                    }
                }
                if ( foundIntersection == false ) {
                    // we must be going outside the mesh
                    // hack: project gradient on face and continue
                    faceNodesPrev = faceNodes;

                    g = projectOnFace(g, faceNodes);

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecEdge(curr_pt, g, faceNodes, pt_i, edgeNodes);

                    if ( foundIntersection == false ) {
                        std::cout << "\n\nWarning: finding raypath (onFace) failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        m_data.resize(0);
                        reachedTx = true;
                        break;
                    }

                    // we might be on one of the nodes
                    if ( pt_i.getDistance(nodes[edgeNodes[0]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[0];
                        pt_i = nodes[edgeNodes[0]];
                        onEdge = false;
                        onFace = false;
                    } else if ( pt_i.getDistance(nodes[edgeNodes[1]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[1];
                        pt_i = nodes[edgeNodes[1]];
                        onEdge = false;
                        onFace = false;
                    } else {
                        onEdge = true;
                        onFace = false;
                    }

                    prev_pt = curr_pt;
                    curr_pt = pt_i;

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                    ds = curr_pt.getDistance( prev_pt );

                    tt += 0.5*(s1 + s2) * ds;
                    s1 = s2;

                    std::set<T2> allNodes;

                    allNodes.insert( faceNodesPrev[0] );
                    allNodes.insert( faceNodesPrev[1] );
                    allNodes.insert( faceNodesPrev[2] );
                    allNodes.insert( edgeNodes[0] );
                    allNodes.insert( edgeNodes[1] );

                    if ( processVel ) {
                        T1 s_sq = computeSlowness(mid_pt, true);
                        s_sq *= s_sq;  // slowness squared
                        update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                    } else {
                        update_m_data(m_data, m, allNodes, mid_pt, ds);
                    }

                }
            } else { // at Rx, somewhere in a tetrahedron

                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( rp_method < 2 ) {
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    getNeighborNodes(cellNo, nnodes);
                    T1 curr_t = Interpolator<T1>::trilinearTime(curr_pt,
                                                                nodes[itmp[0]],
                                                                nodes[itmp[1]],
                                                                nodes[itmp[2]],
                                                                nodes[itmp[3]],
                                                                threadNo);
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(4);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[itmp[0]]);
                    ref_pt[1] = &(nodes[itmp[1]]);
                    ref_pt[2] = &(nodes[itmp[2]]);
                    ref_pt[3] = &(nodes[itmp[3]]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }
                // there are 4 faces that might be intersected

                bool foundIntersection = false;
                for ( size_t n=0; n<4; ++n ) {

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);

                    if ( !foundIntersection )
                        continue;

                    prev_pt = curr_pt;
                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                    ds = curr_pt.getDistance( prev_pt );

                    tt += 0.5*(s1 + s2) * ds;
                    s1 = s2;

                    std::set<T2> allNodes;
                    allNodes.insert( tetrahedra[cellNo].i[0] );
                    allNodes.insert( tetrahedra[cellNo].i[1] );
                    allNodes.insert( tetrahedra[cellNo].i[2] );
                    allNodes.insert( tetrahedra[cellNo].i[3] );

                    if ( processVel ) {
                        T1 s_sq = computeSlowness(mid_pt, true);
                        s_sq *= s_sq;  // slowness squared
                        update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                    } else {
                        update_m_data(m_data, m, allNodes, mid_pt, ds);
                    }

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        m_data.resize(0);
                        reachedTx = true;
                    }
                    break;
                }
                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath within cell failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    tt = 0.0;
                    m_data.resize(0);
                    reachedTx = true;
                }
            }

            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        if ( curr_pt.getDistance( Tx[nt] ) < minDist ) {
                            tt += t0[nt];
                            reachedTx = true;
                            break;
                        }
                    } else if ( txOnEdge[nt] ) {
                        if ( curr_pt.getDistance(nodes[txEdges[nt][0]]) < minDist ||
                            curr_pt.getDistance(nodes[txEdges[nt][1]]) < minDist ) {

                            if ( processVel )
                                s2 = Interpolator<T1>::linearVel(Tx[nt],
                                                                 nodes[txEdges[nt][0]],
                                                                 nodes[txEdges[nt][1]]);
                            else
                                s2 = Interpolator<T1>::linear(Tx[nt],
                                                              nodes[txEdges[nt][0]],
                                                              nodes[txEdges[nt][1]]);

                            reachedTx = true;

                            prev_pt = curr_pt;
                            curr_pt = Tx[nt];
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            ds = curr_pt.getDistance( prev_pt );
                            tt += t0[nt] + 0.5*(s1 + s2) * ds;

                            std::set<T2> allNodes;
                            allNodes.insert(txEdges[nt][0]);
                            allNodes.insert(txEdges[nt][1]);
                            allNodes.insert(nodeNo);

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }

                            break;
                        }
                    } else if ( txOnFace[nt] ) {
                        if ( curr_pt.getDistance(nodes[txFaces[nt][0]]) < minDist ||
                            curr_pt.getDistance(nodes[txFaces[nt][1]]) < minDist ||
                            curr_pt.getDistance(nodes[txFaces[nt][2]]) < minDist ) {

                            if ( processVel )
                                s2 = Interpolator<T1>::bilinearTriangleVel(Tx[nt],
                                                                           nodes[txFaces[nt][0]],
                                                                           nodes[txFaces[nt][1]],
                                                                           nodes[txFaces[nt][2]]);
                            else {
                                s2 = Interpolator<T1>::bilinearTriangle(Tx[nt],
                                                                        nodes[txFaces[nt][0]],
                                                                        nodes[txFaces[nt][1]],
                                                                        nodes[txFaces[nt][2]]);
                            }

                            reachedTx = true;

                            prev_pt = curr_pt;
                            curr_pt = Tx[nt];
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            ds = curr_pt.getDistance( prev_pt );
                            tt += t0[nt] + 0.5*(s1 + s2) * ds;

                            std::set<T2> allNodes;
                            allNodes.insert(txFaces[nt][0]);
                            allNodes.insert(txFaces[nt][1]);
                            allNodes.insert(txFaces[nt][2]);
                            allNodes.insert(nodeNo);

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }

                            break;
                        }
                    }
                }
            } if ( onEdge ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        if ( txNode[nt] == edgeNodes[0] || txNode[nt] == edgeNodes[1] ) {
                            s2 = nodes[txNode[nt]].getNodeSlowness();
                            reachedTx = true;

                            prev_pt = curr_pt;
                            curr_pt = Tx[nt];
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            ds = curr_pt.getDistance( prev_pt );
                            tt += t0[nt] + 0.5*(s1 + s2) * ds;

                            std::set<T2> allNodes;
                            allNodes.insert(edgeNodes[0]);
                            allNodes.insert(edgeNodes[1]);

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }

                            break;
                        }
                    } else {
                        std::array<T2,4> itmp = getPrimary(txCell[nt]);
                        // find adjacent cells
                        const T2 ind[6][2] = {
                            {itmp[0], itmp[1]},
                            {itmp[0], itmp[2]},
                            {itmp[0], itmp[3]},
                            {itmp[1], itmp[2]},
                            {itmp[1], itmp[3]},
                            {itmp[2], itmp[3]} };
                        for ( size_t ne=0; ne<6; ++ne ) {
                            if ( (ind[ne][0] == edgeNodes[0] && ind[ne][1] == edgeNodes[1]) ||
                                (ind[ne][0] == edgeNodes[1] && ind[ne][1] == edgeNodes[0]) ) {
                                if ( processVel )
                                    s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                nodes[itmp[0]],
                                                                                nodes[itmp[1]],
                                                                                nodes[itmp[2]],
                                                                                nodes[itmp[3]]);
                                else
                                    s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                             nodes[itmp[0]],
                                                                             nodes[itmp[1]],
                                                                             nodes[itmp[2]],
                                                                             nodes[itmp[3]]);
                                reachedTx = true;

                                prev_pt = curr_pt;
                                curr_pt = Tx[nt];
                                mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                                ds = curr_pt.getDistance( prev_pt );
                                tt += t0[nt] + 0.5*(s1 + s2) * ds;

                                std::set<T2> allNodes;
                                allNodes.insert(itmp[0]);
                                allNodes.insert(itmp[1]);
                                allNodes.insert(itmp[2]);
                                allNodes.insert(itmp[3]);

                                if ( processVel ) {
                                    T1 s_sq = computeSlowness(mid_pt, true);
                                    s_sq *= s_sq;  // slowness squared
                                    update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                                } else {
                                    update_m_data(m_data, m, allNodes, mid_pt, ds);
                                }

                                break;
                            }
                        }
                    }
                }
            } else {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        for ( auto nc=nodes[txNode[nt]].getOwners().begin(); nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                            if ( cellNo == *nc ) {
                                std::array<T2,4> itmp = getPrimary(cellNo);

                                if ( processVel )
                                    s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                nodes[itmp[0]],
                                                                                nodes[itmp[1]],
                                                                                nodes[itmp[2]],
                                                                                nodes[itmp[3]]);
                                else
                                    s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                             nodes[itmp[0]],
                                                                             nodes[itmp[1]],
                                                                             nodes[itmp[2]],
                                                                             nodes[itmp[3]]);
                                reachedTx = true;

                                prev_pt = curr_pt;
                                curr_pt = Tx[nt];
                                mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                                ds = curr_pt.getDistance( prev_pt );
                                tt += t0[nt] + 0.5*(s1 + s2) * ds;

                                std::set<T2> allNodes;
                                allNodes.insert(itmp[0]);
                                allNodes.insert(itmp[1]);
                                allNodes.insert(itmp[2]);
                                allNodes.insert(itmp[3]);

                                if ( processVel ) {
                                    T1 s_sq = computeSlowness(mid_pt, true);
                                    s_sq *= s_sq;  // slowness squared
                                    update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                                } else {
                                    update_m_data(m_data, m, allNodes, mid_pt, ds);
                                }

                                break;
                            }
                        }
                    } else {
                        if ( cellNo == txCell[nt] ) {
                            std::array<T2,4> itmp = getPrimary(cellNo);
                            if ( processVel )
                                s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                            nodes[itmp[0]],
                                                                            nodes[itmp[1]],
                                                                            nodes[itmp[2]],
                                                                            nodes[itmp[3]]);
                            else
                                s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                         nodes[itmp[0]],
                                                                         nodes[itmp[1]],
                                                                         nodes[itmp[2]],
                                                                         nodes[itmp[3]]);
                            reachedTx = true;

                            prev_pt = curr_pt;
                            curr_pt = Tx[nt];
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            ds = curr_pt.getDistance( prev_pt );
                            tt += t0[nt] + 0.5*(s1 + s2) * ds;

                            std::set<T2> allNodes;
                            allNodes.insert(itmp[0]);
                            allNodes.insert(itmp[1]);
                            allNodes.insert(itmp[2]);
                            allNodes.insert(itmp[3]);

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }

                        } else {
                            for ( size_t nn=0; nn<txNeighborCells[nt].size(); ++nn ) {
                                if ( cellNo == txNeighborCells[nt][nn] ) {
                                    std::array<T2,4> itmp = getPrimary(txCell[nt]);
                                    std::array<T2,3> ind[4] = {
                                        { { itmp[0], itmp[1], itmp[2] } },
                                        { { itmp[0], itmp[1], itmp[3] } },
                                        { { itmp[0], itmp[2], itmp[3] } },
                                        { { itmp[1], itmp[2], itmp[3] } }
                                    };
                                    bool found = false;
                                    for ( size_t n=0; n<4; ++n ) {
                                        std::sort( ind[n].begin(), ind[n].end() );
                                        if ( faceNodes == ind[n] ) {
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found ) {
                                        itmp = getPrimary(cellNo);
                                        if ( processVel )
                                            s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                        nodes[itmp[0]],
                                                                                        nodes[itmp[1]],
                                                                                        nodes[itmp[2]],
                                                                                        nodes[itmp[3]]);
                                        else
                                            s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                                     nodes[itmp[0]],
                                                                                     nodes[itmp[1]],
                                                                                     nodes[itmp[2]],
                                                                                     nodes[itmp[3]]);

                                        reachedTx = true;

                                        prev_pt = curr_pt;
                                        curr_pt = Tx[nt];
                                        mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                                        ds = curr_pt.getDistance( prev_pt );
                                        tt += t0[nt] + 0.5*(s1 + s2) * ds;

                                        std::set<T2> allNodes;
                                        allNodes.insert(itmp[0]);
                                        allNodes.insert(itmp[1]);
                                        allNodes.insert(itmp[2]);
                                        allNodes.insert(itmp[3]);

                                        if ( processVel ) {
                                            T1 s_sq = computeSlowness(mid_pt, true);
                                            s_sq *= s_sq;  // slowness squared
                                            update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                                        } else {
                                            update_m_data(m_data, m, allNodes, mid_pt, ds);
                                        }

                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if ( reachedTx ) break;
                }
            }
            itertions ++;
            if (itertions > max_itertions){
                std::cout << "\n\nWarning: raypath failed to converge for Rx : "
                << Rx.x << ' ' << Rx.y << ' ' << Rx.z << " (infinite loop risk)"<<std::endl;
                tt = 0.0;
                m_data.resize(0);
                break; // break while loop
            }
        }

        delete grad3d;
    }


    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                          const sxyz<T1> &Rx,
                                          std::vector<sxyz<T1>>& r_data,
                                          std::vector<sijv<T1>>& m_data,
                                          const size_t RxNo,
                                          const size_t threadNo) const {
        // Important note: m_data holds terms of traveltime derivative w/r to slowness (not velocity)
        T1 minDist = small;
        std::vector<sxyz<T1>> r_tmp;
        r_tmp.emplace_back( Rx );

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                return;
            }
        }

        const txInfo_t& txi = getTxInfo(Tx, threadNo);
        const std::vector<bool>& txOnNode = txi.txOnNode;
        const std::vector<bool>& txOnEdge = txi.txOnEdge;
        const std::vector<bool>& txOnFace = txi.txOnFace;
        const std::vector<T2>& txNode = txi.txNode;
        const std::vector<T2>& txCell = txi.txCell;
        const std::vector<std::array<T2,2>>& txEdges = txi.txEdges;
        const std::vector<std::array<T2,3>>& txFaces = txi.txFaces;
        const std::vector<std::vector<T2>>& txNeighborCells = txi.txNeighborCells;

        T2 cellNo, nodeNo, nodeNoPrev;
        sxyz<T1> curr_pt( Rx ), mid_pt, prev_pt( Rx );
        sijv<T1> m;
        m.i = RxNo;

        bool onNode = false;
        bool onEdge = false;
        bool onFace = false;
        std::array<T2,2> edgeNodes, edgeNodesPrev;
        std::array<T2,3> faceNodes, faceNodesPrev;
        Grad3D<T1,NODE>* grad3d = nullptr;
        if ( rp_method == 0 ) {
            grad3d = new Grad3D_ls_fo<T1,NODE>();
        } else if ( rp_method == 1 ) {
            grad3d = new Grad3D_ls_so<T1,NODE>();
        } else if ( rp_method == 2 ) {
            grad3d = new Grad3D_ab<T1,NODE>();
        }
        bool reachedTx = false;

        for ( T2 nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == curr_pt ) {
                nodeNo = nn;
                onNode = true;
                break;
            }
        }
        if ( !onNode ) {
            cellNo = getCellNo( curr_pt );

            std::array<T2,4> itmp = getPrimary(cellNo);
            T2 ind[6][2] = {
                {itmp[0], itmp[1]},
                {itmp[0], itmp[2]},
                {itmp[0], itmp[3]},
                {itmp[1], itmp[2]},
                {itmp[1], itmp[3]},
                {itmp[2], itmp[3]} };

            for ( size_t n=0; n<6; ++n ) {
                if ( areCollinear(curr_pt, nodes[ind[n][0]], nodes[ind[n][1]]) ) {
                    onEdge = true;
                    edgeNodes[0] = ind[n][0];
                    edgeNodes[1] = ind[n][1];
                    break;
                }
            }
            if ( !onEdge ) {
                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }

                for ( size_t n=0; n<4; ++n ) {
                    if ( areCoplanar(curr_pt, nodes[ind[n][0]], nodes[ind[n][1]], nodes[ind[n][2]]) ) {
                        onFace = true;
                        //  faceNodes shoud not be assigned, face was not intersected
                        break;
                    }
                }
            }
        }

        for ( auto nt=0; nt<txCell.size(); ++nt ) {
            if ( cellNo == txCell[nt] ) {
                r_tmp.push_back( Tx[nt] );
                reachedTx = true;
                break;
            }
        }

        T1 ds;
        sxyz<T1> g;
        while ( reachedTx == false ) {

            if ( onNode ) {

                // find cells common to edge
                std::byte nnodes_buf[8192];
                std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                std::pmr::set<NODE*> nnodes{&nnodes_pool};
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                    getNeighborNodes(*nc, nnodes);
                }

                // compute gradient with nodes from all common cells
                if ( rp_method < 2 ) {
                    // find cells common to edge
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                        getNeighborNodes(*nc, nnodes);
                    }
                    // compute gradient with nodes from all common cells
                    g = grad3d->compute(curr_pt, nodes[nodeNo].getTT(threadNo), nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(1);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[nodeNo]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }

                nodeNoPrev = nodeNo;

                // find cell for which gradient intersect opposing face
                bool foundIntersection = false;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                    std::array<T2,3> nb;
                    size_t n=0;
                    for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                        if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                            nb[n++] = *nn;
                        }
                    }
                    std::sort(nb.begin(), nb.end());

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle( nodeNo, g, nb[0], nb[1], nb[2], pt_i);
                    if ( !foundIntersection ) {
                        continue;
                    }

                    //  check if cell is (one of) TxCell(s)
                    for (size_t nt=0; nt<Tx.size(); ++nt) {
                        if ( *nc == txCell[nt] ) {
                            std::array<T2,4> itmp = getPrimary(*nc);

                            r_tmp.push_back( Tx[nt] );

                            prev_pt = curr_pt;
                            curr_pt = Tx[nt];
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            ds = curr_pt.getDistance( prev_pt );
                            std::set<T2> allNodes;
                            allNodes.insert(itmp[0]);
                            allNodes.insert(itmp[1]);
                            allNodes.insert(itmp[2]);
                            allNodes.insert(itmp[3]);

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }

                            reachedTx = true;
                            break;
                        }
                    }
                    if ( reachedTx ) {
                        break;
                    }


                    prev_pt = curr_pt;
                    curr_pt = pt_i;
                    r_tmp.push_back( curr_pt );

                    if ( r_tmp.size() > 1 ) {
                        // compute terms of matrix M
                        mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                        ds = curr_pt.getDistance( prev_pt );
                    }

                    bool break_flag=false;
                    for ( n=0; n<3; ++n ) {
                        if ( nodes[ nb[n] ].getDistance( curr_pt ) < small ) {

                            nodeNo = nb[n];
                            onNode = true;
                            onEdge = false;
                            onFace = false;

                            if ( r_tmp.size() > 1) {
                                std::set<T2> allNodes;
                                allNodes.insert(nodeNoPrev);
                                allNodes.insert(nodeNo);

                                if ( processVel ) {
                                    T1 s_sq = computeSlowness(mid_pt, true);
                                    s_sq *= s_sq;  // slowness squared
                                    update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                                } else {
                                    update_m_data(m_data, m, allNodes, mid_pt, ds);
                                }
                            }

                            break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;

                    for ( size_t n1=0; n1<3; ++n1 ) {
                        size_t n2 = (n1+1)%3;
                        if ( areCollinear(curr_pt, nodes[nb[n1]], nodes[nb[n2]]) ) {
                            edgeNodesPrev = edgeNodes;

                            edgeNodes[0] = nb[n1];
                            edgeNodes[1] = nb[n2];
                            onNode = false;
                            onEdge = true;
                            onFace = false;

                            if ( r_tmp.size() > 1) {
                                std::set<T2> allNodes;
                                allNodes.insert(nodeNoPrev);
                                allNodes.insert( edgeNodes[0] );
                                allNodes.insert( edgeNodes[1] );

                                if ( processVel ) {
                                    T1 s_sq = computeSlowness(mid_pt, true);
                                    s_sq *= s_sq;  // slowness squared
                                    update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                                } else {
                                    update_m_data(m_data, m, allNodes, mid_pt, ds);
                                }
                            }

                            break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;

                    onNode = false;
                    onEdge = false;
                    onFace = true;

                    faceNodesPrev = faceNodes;
                    faceNodes = nb;

                    if ( r_tmp.size() > 1) {
                        std::set<T2> allNodes;
                        allNodes.insert(nodeNoPrev);
                        allNodes.insert( faceNodes[0] );
                        allNodes.insert( faceNodes[1] );
                        allNodes.insert( faceNodes[2] );

                        if ( processVel ) {
                            T1 s_sq = computeSlowness(mid_pt, true);
                            s_sq *= s_sq;  // slowness squared
                            update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                        } else {
                            update_m_data(m_data, m, allNodes, mid_pt, ds);
                        }
                    }

                    // find next cell
                    cellNo = findAdjacentCell1(faceNodes, nodeNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        m_data.resize(0);
                        reachedTx = true;
                    }
                    break;
                }

                if ( foundIntersection == false ) {
                    nodeNoPrev = nodeNo;
                    // projet gradient on face and find intersection on next edge
                    sxyz<T1> pt_i;
                    foundIntersection = projectOnFace(curr_pt, nodeNo, g, edgeNodes, pt_i, threadNo);

                    if ( foundIntersection == false ) {
                        std::cout << "\n\nWarning: finding raypath (onNode) failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        m_data.resize(0);
                        reachedTx = true;
                        break;
                    }
                    // find which cell we are in
                    for ( size_t nc=0; nc<tetrahedra.size(); ++nc ) {
                        std::array<T2,4> itet = {tetrahedra[nc].i[0],
                            tetrahedra[nc].i[1],
                            tetrahedra[nc].i[2],
                            tetrahedra[nc].i[3]};
                        // because we are at the surface, there is only one cell with nodeNo & edgeNodes
                        if (std::find(itet.begin(), itet.end(), nodeNo) != itet.end() &&
                            std::find(itet.begin(), itet.end(), edgeNodes[0]) != itet.end() &&
                            std::find(itet.begin(), itet.end(), edgeNodes[1]) != itet.end()) {
                            cellNo = static_cast<T2>(nc);
                            break;
                        }
                    }

                    //  check if cell is (one of) TxCell(s)
                    for (size_t nt=0; nt<Tx.size(); ++nt) {
                        if ( cellNo == txCell[nt] ) {

                            r_tmp.push_back( Tx[nt] );

                            mid_pt = static_cast<T1>(0.5)*(curr_pt + Tx[nt]);
                            ds = curr_pt.getDistance( Tx[nt] );

                            std::set<T2> allNodes;
                            allNodes.insert( tetrahedra[cellNo].i[0] );
                            allNodes.insert( tetrahedra[cellNo].i[1] );
                            allNodes.insert( tetrahedra[cellNo].i[2] );
                            allNodes.insert( tetrahedra[cellNo].i[3] );

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }

                            reachedTx = true;
                            break;
                        }
                    }
                    if ( reachedTx ) {
                        break;
                    }

                    // we might be on one of the nodes
                    if ( pt_i.getDistance(nodes[edgeNodes[0]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[0];
                        pt_i = nodes[edgeNodes[0]];
                        onEdge = false;
                    } else if ( pt_i.getDistance(nodes[edgeNodes[1]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[1];
                        pt_i = nodes[edgeNodes[1]];
                        onEdge = false;
                    } else {
                        onEdge = true;
                        onNode = false;
                    }

                    prev_pt = curr_pt;
                    curr_pt = pt_i;

                    r_tmp.push_back( curr_pt );

                    mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                    ds = curr_pt.getDistance( prev_pt );
                    std::set<T2> allNodes;
                    allNodes.insert( nodeNoPrev );
                    allNodes.insert( edgeNodes[0] );
                    allNodes.insert( edgeNodes[1] );

                    if ( processVel ) {
                        T1 s_sq = computeSlowness(mid_pt, true);
                        s_sq *= s_sq;  // slowness squared
                        update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                    } else {
                        update_m_data(m_data, m, allNodes, mid_pt, ds);
                    }
                }

            } else if ( onEdge ) {

                // find cells common to edge
                std::vector<T2> cells;
                std::byte nnodes_buf[8192];
                std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                std::pmr::set<NODE*> nnodes{&nnodes_pool};
                for ( auto nc0=nodes[edgeNodes[0]].getOwners().begin(); nc0!=nodes[edgeNodes[0]].getOwners().end(); ++nc0 ) {
                    if ( std::find(nodes[edgeNodes[1]].getOwners().begin(), nodes[edgeNodes[1]].getOwners().end(), *nc0)!=nodes[edgeNodes[1]].getOwners().end() ) {
                        cells.push_back( *nc0 );
                        getNeighborNodes(*nc0, nnodes);
                    }
                }
                if ( rp_method < 2 ) {
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    for (size_t n=0; n<cells.size(); ++n ) {
                        getNeighborNodes(cells[n], nnodes);
                    }
                    T1 d01 = nodes[edgeNodes[0]].getDistance(nodes[edgeNodes[1]]);
                    T1 w0 = curr_pt.getDistance(nodes[edgeNodes[1]]) / d01;
                    T1 w1 = curr_pt.getDistance(nodes[edgeNodes[0]]) / d01;
                    T1 curr_t = nodes[edgeNodes[0]].getTT(threadNo)*w0 + nodes[edgeNodes[1]].getTT(threadNo)*w1;
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(2);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[edgeNodes[0]]);
                    ref_pt[1] = &(nodes[edgeNodes[1]]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, edgeNodes, Tx, txCell);

                edgeNodesPrev = edgeNodes;

                bool foundIntersection = false;
                for (size_t n=0; n<cells.size(); ++n ) {

                    T2 testCellNo = cells[n];

                    // there are 2 faces that might be intersected
                    std::array<T2,2> edgeNodes2;
                    size_t n2=0;
                    for ( auto nn=this->neighbors[testCellNo].begin(); nn!= this->neighbors[testCellNo].end(); ++nn ) {
                        if ( *nn!=edgeNodes[0] && *nn!=edgeNodes[1] && nodes[*nn].isPrimary() ) {
                            edgeNodes2[n2++] = *nn;
                        }
                    }

                    sxyz<T1> pt_i;
                    T2 itmpNode;
                    foundIntersection = intersectVecTriangle(curr_pt, g,
                                                             edgeNodes[0],
                                                             edgeNodes2[0],
                                                             edgeNodes2[1], pt_i);
                    itmpNode = edgeNodes[0];
                    if ( !foundIntersection ) {
                        foundIntersection = intersectVecTriangle(curr_pt, g,
                                                                 edgeNodes[1],
                                                                 edgeNodes2[0],
                                                                 edgeNodes2[1], pt_i);
                        itmpNode = edgeNodes[1];
                    }
                    if ( !foundIntersection ) {
                        continue;
                    }

                    cellNo = testCellNo;
                    prev_pt = curr_pt;
                    curr_pt = pt_i;
                    r_tmp.push_back( curr_pt );

                    if (r_tmp.size() > 1 ) {
                        // compute terms of matrix M
                        mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                        ds = curr_pt.getDistance( prev_pt );
                    }

                    bool break_flag = false;
                    for ( size_t n2=0; n2<4; ++n2 ) {
                        if ( nodes[ this->neighbors[cellNo][n2] ].getDistance( curr_pt ) < small ) {

                            nodeNo = this->neighbors[cellNo][n2];
                            onNode = true;
                            onEdge = false;
                            onFace = false;

                            if ( r_tmp.size() > 1) {
                                std::set<T2> allNodes;
                                allNodes.insert(nodeNo);
                                allNodes.insert( edgeNodesPrev[0] );
                                allNodes.insert( edgeNodesPrev[1] );

                                if ( processVel ) {
                                    T1 s_sq = computeSlowness(mid_pt, true);
                                    s_sq *= s_sq;  // slowness squared
                                    update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                                } else {
                                    update_m_data(m_data, m, allNodes, mid_pt, ds);
                                }
                            }

                            break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;

                    if ( areCollinear(curr_pt, nodes[itmpNode], nodes[edgeNodes2[0]]) ) {

                        edgeNodes[0] = itmpNode;
                        edgeNodes[1] = edgeNodes2[0];
                        onNode = false;
                        onEdge = true;
                        onFace = false;

                        if ( r_tmp.size() > 1) {
                            std::set<T2> allNodes;
                            allNodes.insert( edgeNodesPrev[0] );
                            allNodes.insert( edgeNodesPrev[1] );
                            allNodes.insert( edgeNodes[0] );
                            allNodes.insert( edgeNodes[1] );

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }
                        }

                        break_flag = true;
                        break;
                    } else if ( areCollinear(curr_pt, nodes[itmpNode], nodes[edgeNodes2[1]]) ) {

                        edgeNodes[0] = itmpNode;
                        edgeNodes[1] = edgeNodes2[1];
                        onNode = false;
                        onEdge = true;
                        onFace = false;

                        if ( r_tmp.size() > 1) {
                            std::set<T2> allNodes;
                            allNodes.insert( edgeNodesPrev[0] );
                            allNodes.insert( edgeNodesPrev[1] );
                            allNodes.insert( edgeNodes[0] );
                            allNodes.insert( edgeNodes[1] );

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }
                        }

                        break_flag = true;
                        break;
                    } else if ( areCollinear(curr_pt, nodes[edgeNodes2[0]], nodes[edgeNodes2[1]]) ) {

                        edgeNodes[0] = edgeNodes2[0];
                        edgeNodes[1] = edgeNodes2[1];
                        onNode = false;
                        onEdge = true;
                        onFace = false;

                        if ( r_tmp.size() > 1) {
                            std::set<T2> allNodes;
                            allNodes.insert( edgeNodesPrev[0] );
                            allNodes.insert( edgeNodesPrev[1] );
                            allNodes.insert( edgeNodes[0] );
                            allNodes.insert( edgeNodes[1] );

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }
                        }

                        break_flag = true;
                        break;
                    }
                    if ( break_flag ) break;

                    onNode = false;
                    onEdge = false;
                    onFace = true;

                    faceNodesPrev = faceNodes;
                    faceNodes[0] = itmpNode;
                    faceNodes[1] = edgeNodes2[0];
                    faceNodes[2] = edgeNodes2[1];
                    std::sort(faceNodes.begin(), faceNodes.end());

                    if ( r_tmp.size() > 1) {
                        std::set<T2> allNodes;
                        allNodes.insert( edgeNodesPrev[0] );
                        allNodes.insert( edgeNodesPrev[1] );
                        allNodes.insert( faceNodes[0] );
                        allNodes.insert( faceNodes[1] );
                        allNodes.insert( faceNodes[2] );

                        if ( processVel ) {
                            T1 s_sq = computeSlowness(mid_pt, true);
                            s_sq *= s_sq;  // slowness squared
                            update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                        } else {
                            update_m_data(m_data, m, allNodes, mid_pt, ds);
                        }
                    }

                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        m_data.resize(0);
                        reachedTx = true;
                    }
                    break;
                }
                if ( foundIntersection == false ) {

                    edgeNodesPrev = edgeNodes;
                    sxyz<T1> pt_i;
                    g = projectOnFace(curr_pt, g, edgeNodes, cells, pt_i, threadNo);
                    if ( g.x==0.0 && g.y==0.0 && g.z==0.0 ) {
                        foundIntersection = false;
                    } else {
                        foundIntersection = true;
                    }

                    if ( foundIntersection == false || curr_pt == pt_i ) {
                        std::cout << "\n\nWarning: finding raypath (onEdge) failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        m_data.resize(0);
                        reachedTx = true;
                        break;
                    }

                    // we might be on one of the nodes
                    if ( pt_i.getDistance(nodes[edgeNodes[0]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[0];
                        pt_i = nodes[edgeNodes[0]];
                        onEdge = false;
                    } else if ( pt_i.getDistance(nodes[edgeNodes[1]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[1];
                        pt_i = nodes[edgeNodes[1]];
                        onEdge = false;
                    }

                    prev_pt = curr_pt;
                    curr_pt = pt_i;

                    r_tmp.push_back( curr_pt );

                    mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                    ds = curr_pt.getDistance( prev_pt );
                    std::set<T2> allNodes;

                    allNodes.insert( edgeNodesPrev[0] );
                    allNodes.insert( edgeNodesPrev[1] );
                    allNodes.insert( edgeNodes[0] );
                    allNodes.insert( edgeNodes[1] );

                    if ( processVel ) {
                        T1 s_sq = computeSlowness(mid_pt, true);
                        s_sq *= s_sq;  // slowness squared
                        update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                    } else {
                        update_m_data(m_data, m, allNodes, mid_pt, ds);
                    }
                }

            } else if ( onFace ) { // on Face

                std::byte nnodes_buf[8192];
                std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                std::pmr::set<NODE*> nnodes{&nnodes_pool};
                getNeighborNodes(cellNo, nnodes);

                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( rp_method < 2 ) {
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    getNeighborNodes(cellNo, nnodes);
                    T1 curr_t;
                    if ( r_tmp.size() <= 1 ) {
                        curr_t = Interpolator<T1>::trilinearTime(curr_pt,
                                                                 nodes[itmp[0]],
                                                                 nodes[itmp[1]],
                                                                 nodes[itmp[2]],
                                                                 nodes[itmp[3]],
                                                                 threadNo);
                    } else {
                        curr_t = Interpolator<T1>::bilinearTime(curr_pt,
                                                                nodes[faceNodes[0]],
                                                                nodes[faceNodes[1]],
                                                                nodes[faceNodes[2]],
                                                                threadNo);
                    }
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(3);
                    if ( r_tmp.size() <= 1 ) {
                        ref_pt[0] = &(nodes[itmp[0]]);
                        ref_pt[1] = &(nodes[itmp[1]]);
                        ref_pt[2] = &(nodes[itmp[2]]);
                        ref_pt.push_back( &(nodes[itmp[3]]) );
                    } else {
                        ref_pt[0] = &(nodes[faceNodes[0]]);
                        ref_pt[1] = &(nodes[faceNodes[1]]);
                        ref_pt[2] = &(nodes[faceNodes[2]]);
                    }
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }
                // there are 3 faces that might be intersected

                faceNodesPrev = faceNodes;

                bool foundIntersection = false;
                for ( size_t n=0; n<4; ++n ) {
                    if ( ind[n] == faceNodes ) continue;

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);

                    if ( !foundIntersection )
                        continue;

                    prev_pt = curr_pt;
                    curr_pt = pt_i;
                    r_tmp.push_back( curr_pt );

                    if (r_tmp.size() > 1 ) {
                        // compute terms of matrix M
                        mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                        ds = curr_pt.getDistance( prev_pt );
                    }

                    bool break_flag = false;
                    for ( size_t n2=0; n2<3; ++n2 ) {
                        if ( nodes[ ind[n][n2] ].getDistance( curr_pt ) < small ) {
                            nodeNoPrev = nodeNo;

                            nodeNo = ind[n][n2];
                            onNode = true;
                            onEdge = false;
                            onFace = false;

                            if ( r_tmp.size() > 1) {
                                std::set<T2> allNodes;
                                allNodes.insert(nodeNo);
                                allNodes.insert( faceNodesPrev[0] );
                                allNodes.insert( faceNodesPrev[1] );
                                allNodes.insert( faceNodesPrev[2] );

                                if ( processVel ) {
                                    T1 s_sq = computeSlowness(mid_pt, true);
                                    s_sq *= s_sq;  // slowness squared
                                    update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                                } else {
                                    update_m_data(m_data, m, allNodes, mid_pt, ds);
                                }
                            }

                            break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;

                    for ( size_t n1=0; n1<3; ++n1 ) {
                        size_t n2 = (n1+1)%3;
                        if ( areCollinear(curr_pt, nodes[ind[n][n1]], nodes[ind[n][n2]]) ) {
                            edgeNodesPrev = edgeNodes;

                            edgeNodes[0] = ind[n][n1];
                            edgeNodes[1] = ind[n][n2];
                            onNode = false;
                            onEdge = true;
                            onFace = false;

                            if ( r_tmp.size() > 1) {
                                std::set<T2> allNodes;
                                allNodes.insert( edgeNodes[0] );
                                allNodes.insert( edgeNodes[1] );
                                allNodes.insert( faceNodesPrev[0] );
                                allNodes.insert( faceNodesPrev[1] );
                                allNodes.insert( faceNodesPrev[2] );

                                if ( processVel ) {
                                    T1 s_sq = computeSlowness(mid_pt, true);
                                    s_sq *= s_sq;  // slowness squared
                                    update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                                } else {
                                    update_m_data(m_data, m, allNodes, mid_pt, ds);
                                }
                            }

                            break_flag = true;
                            break;
                        }
                    }
                    if ( break_flag ) break;

                    onNode = false;
                    onEdge = false;
                    onFace = true;

                    faceNodes = ind[n];

                    if ( r_tmp.size() > 1) {
                        std::set<T2> allNodes;
                        allNodes.insert( faceNodesPrev[0] );
                        allNodes.insert( faceNodesPrev[1] );
                        allNodes.insert( faceNodesPrev[2] );
                        allNodes.insert( faceNodes[0] );
                        allNodes.insert( faceNodes[1] );
                        allNodes.insert( faceNodes[2] );

                        if ( processVel ) {
                            T1 s_sq = computeSlowness(mid_pt, true);
                            s_sq *= s_sq;  // slowness squared
                            update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                        } else {
                            update_m_data(m_data, m, allNodes, mid_pt, ds);
                        }
                    }

                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        m_data.resize(0);
                        reachedTx = true;
                    }
                    break;
                }

                if ( foundIntersection == false ) {

                    // we must be on an face with gradient pointing slightly outward tetrahedron
                    // return in other cell but keep gradient
                    cellNo = findAdjacentCell2(faceNodes, cellNo);

                    std::array<T2,4> itmp = getPrimary(cellNo);
                    ind[0] = { { itmp[0], itmp[1], itmp[2] } };
                    ind[1] = { { itmp[0], itmp[1], itmp[3] } };
                    ind[2] = { { itmp[0], itmp[2], itmp[3] } };
                    ind[3] = { { itmp[1], itmp[2], itmp[3] } };

                    for ( size_t n=0; n<4; ++n ) {
                        std::sort( ind[n].begin(), ind[n].end() );
                    }

                    for ( size_t n=0; n<4; ++n ) {
                        if ( ind[n] == faceNodes ) continue;

                        sxyz<T1> pt_i;
                        foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                                 ind[n][1], ind[n][2],
                                                                 pt_i);

                        if ( !foundIntersection ) {
                            continue;
                        }
                        prev_pt = curr_pt;
                        curr_pt = pt_i;
                        r_tmp.push_back( curr_pt );

                        if (r_tmp.size() > 1 ) {
                            // compute terms of matrix M
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            ds = curr_pt.getDistance( prev_pt );
                        }

                        bool break_flag = false;
                        for ( size_t n2=0; n2<3; ++n2 ) {
                            if ( nodes[ ind[n][n2] ].getDistance( curr_pt ) < small ) {
                                nodeNoPrev = nodeNo;

                                nodeNo = ind[n][n2];
                                onNode = true;
                                onEdge = false;
                                onFace = false;

                                if ( r_tmp.size() > 1) {
                                    std::set<T2> allNodes;
                                    allNodes.insert(nodeNo);
                                    allNodes.insert( faceNodesPrev[0] );
                                    allNodes.insert( faceNodesPrev[1] );
                                    allNodes.insert( faceNodesPrev[2] );

                                    if ( processVel ) {
                                        T1 s_sq = computeSlowness(mid_pt, true);
                                        s_sq *= s_sq;  // slowness squared
                                        update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                                    } else {
                                        update_m_data(m_data, m, allNodes, mid_pt, ds);
                                    }
                                }

                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;

                        for ( size_t n1=0; n1<3; ++n1 ) {
                            size_t n2 = (n1+1)%3;
                            if ( areCollinear(curr_pt, nodes[ind[n][n1]], nodes[ind[n][n2]]) ) {
                                edgeNodesPrev = edgeNodes;

                                edgeNodes[0] = ind[n][n1];
                                edgeNodes[1] = ind[n][n2];
                                onNode = false;
                                onEdge = true;
                                onFace = false;

                                if ( r_tmp.size() > 1) {
                                    std::set<T2> allNodes;
                                    allNodes.insert( edgeNodes[0] );
                                    allNodes.insert( edgeNodes[1] );
                                    allNodes.insert( faceNodesPrev[0] );
                                    allNodes.insert( faceNodesPrev[1] );
                                    allNodes.insert( faceNodesPrev[2] );

                                    if ( processVel ) {
                                        T1 s_sq = computeSlowness(mid_pt, true);
                                        s_sq *= s_sq;  // slowness squared
                                        update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                                    } else {
                                        update_m_data(m_data, m, allNodes, mid_pt, ds);
                                    }
                                }

                                break_flag = true;
                                break;
                            }
                        }
                        if ( break_flag ) break;

                        onNode = false;
                        onEdge = false;
                        onFace = true;

                        faceNodes = ind[n];

                        if ( r_tmp.size() > 1) {
                            std::set<T2> allNodes;
                            allNodes.insert( faceNodesPrev[0] );
                            allNodes.insert( faceNodesPrev[1] );
                            allNodes.insert( faceNodesPrev[2] );
                            allNodes.insert( faceNodes[0] );
                            allNodes.insert( faceNodes[1] );
                            allNodes.insert( faceNodes[2] );

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }
                        }

                        // find next cell
                        cellNo = findAdjacentCell2(faceNodes, cellNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                            << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                            r_tmp.resize(1);
                            r_tmp[0] = Rx;
                            m_data.resize(0);
                            reachedTx = true;
                        }
                        break;
                    }
                }
                if ( foundIntersection == false ) {
                    faceNodesPrev = faceNodes;

                    g = projectOnFace(g, faceNodes);

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecEdge(curr_pt, g, faceNodes, pt_i, edgeNodes);

                    if ( foundIntersection == false ) {
                        std::cout << "\n\nWarning: finding raypath (onFace) failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        m_data.resize(0);
                        reachedTx = true;
                        break;
                    }

                    // we might be on one of the nodes
                    if ( pt_i.getDistance(nodes[edgeNodes[0]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[0];
                        pt_i = nodes[edgeNodes[0]];
                        onEdge = false;
                        onFace = false;
                    } else if ( pt_i.getDistance(nodes[edgeNodes[1]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[1];
                        pt_i = nodes[edgeNodes[1]];
                        onEdge = false;
                        onFace = false;
                    } else {
                        onEdge = true;
                        onFace = false;
                    }

                    prev_pt = curr_pt;
                    curr_pt = pt_i;

                    r_tmp.push_back( curr_pt );

                    mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                    ds = curr_pt.getDistance( prev_pt );
                    std::set<T2> allNodes;

                    allNodes.insert( faceNodesPrev[0] );
                    allNodes.insert( faceNodesPrev[1] );
                    allNodes.insert( faceNodesPrev[2] );
                    allNodes.insert( edgeNodes[0] );
                    allNodes.insert( edgeNodes[1] );

                    if ( processVel ) {
                        T1 s_sq = computeSlowness(mid_pt, true);
                        s_sq *= s_sq;  // slowness squared
                        update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                    } else {
                        update_m_data(m_data, m, allNodes, mid_pt, ds);
                    }
                }
            } else { // at Rx, somewhere in a tetrahedron

                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( rp_method < 2 ) {
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    getNeighborNodes(cellNo, nnodes);
                    T1 curr_t = Interpolator<T1>::trilinearTime(curr_pt,
                                                                nodes[itmp[0]],
                                                                nodes[itmp[1]],
                                                                nodes[itmp[2]],
                                                                nodes[itmp[3]],
                                                                threadNo);
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(4);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[itmp[0]]);
                    ref_pt[1] = &(nodes[itmp[1]]);
                    ref_pt[2] = &(nodes[itmp[2]]);
                    ref_pt[3] = &(nodes[itmp[3]]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }
                // there are 4 faces that might be intersected

                bool foundIntersection = false;
                for ( size_t n=0; n<4; ++n ) {

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);

                    if ( !foundIntersection )
                        continue;

                    prev_pt = curr_pt;
                    curr_pt = pt_i;
                    r_tmp.push_back( curr_pt );
                    mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                    ds = curr_pt.getDistance( prev_pt );

                    std::set<T2> allNodes;
                    allNodes.insert( tetrahedra[cellNo].i[0] );
                    allNodes.insert( tetrahedra[cellNo].i[1] );
                    allNodes.insert( tetrahedra[cellNo].i[2] );
                    allNodes.insert( tetrahedra[cellNo].i[3] );

                    if ( processVel ) {
                        T1 s_sq = computeSlowness(mid_pt, true);
                        s_sq *= s_sq;  // slowness squared
                        update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                    } else {
                        update_m_data(m_data, m, allNodes, mid_pt, ds);
                    }

                    bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        m_data.resize(0);
                        reachedTx = true;
                    }
                    break;
                }
                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath within cell failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_data.resize(1);
                    r_data[0] = Rx;
                    m_data.resize(0);
                    reachedTx = true;
                }
            }

            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        if ( curr_pt.getDistance( Tx[nt] ) < minDist ) {
                            reachedTx = true;
                            break;
                        }
                    } else if ( txOnEdge[nt] ) {
                        if ( curr_pt.getDistance(nodes[txEdges[nt][0]]) < minDist ||
                            curr_pt.getDistance(nodes[txEdges[nt][1]]) < minDist ) {

                            r_tmp.push_back(Tx[nt]);
                            reachedTx = true;

                            prev_pt = curr_pt;
                            curr_pt = Tx[nt];
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            ds = curr_pt.getDistance( prev_pt );
                            std::set<T2> allNodes;
                            allNodes.insert(txEdges[nt][0]);
                            allNodes.insert(txEdges[nt][1]);
                            allNodes.insert(nodeNo);

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }

                            break;
                        }
                    } else if ( txOnFace[nt] ) {
                        if ( curr_pt.getDistance(nodes[txFaces[nt][0]]) < minDist ||
                            curr_pt.getDistance(nodes[txFaces[nt][1]]) < minDist ||
                            curr_pt.getDistance(nodes[txFaces[nt][2]]) < minDist ) {

                            r_tmp.push_back(Tx[nt]);
                            reachedTx = true;

                            prev_pt = curr_pt;
                            curr_pt = Tx[nt];
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            ds = curr_pt.getDistance( prev_pt );
                            std::set<T2> allNodes;
                            allNodes.insert(txFaces[nt][0]);
                            allNodes.insert(txFaces[nt][1]);
                            allNodes.insert(txFaces[nt][2]);
                            allNodes.insert(nodeNo);

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }

                            break;
                        }
                    }
                }
            } if ( onEdge ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        if ( txNode[nt] == edgeNodes[0] || txNode[nt] == edgeNodes[1] ) {
                            r_tmp.push_back( Tx[nt] );
                            reachedTx = true;

                            prev_pt = curr_pt;
                            curr_pt = Tx[nt];
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            ds = curr_pt.getDistance( prev_pt );
                            std::set<T2> allNodes;
                            allNodes.insert(edgeNodes[0]);
                            allNodes.insert(edgeNodes[1]);

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }

                            break;
                        }
                    } else {
                        std::array<T2,4> itmp = getPrimary(txCell[nt]);
                        // find adjacent cells
                        const T2 ind[6][2] = {
                            {itmp[0], itmp[1]},
                            {itmp[0], itmp[2]},
                            {itmp[0], itmp[3]},
                            {itmp[1], itmp[2]},
                            {itmp[1], itmp[3]},
                            {itmp[2], itmp[3]} };
                        for ( size_t ne=0; ne<6; ++ne ) {
                            if ( (ind[ne][0] == edgeNodes[0] && ind[ne][1] == edgeNodes[1]) ||
                                (ind[ne][0] == edgeNodes[1] && ind[ne][1] == edgeNodes[0]) ) {
                                r_tmp.push_back( Tx[nt] );
                                reachedTx = true;

                                prev_pt = curr_pt;
                                curr_pt = Tx[nt];
                                mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                                ds = curr_pt.getDistance( prev_pt );
                                std::set<T2> allNodes;
                                allNodes.insert(itmp[0]);
                                allNodes.insert(itmp[1]);
                                allNodes.insert(itmp[2]);
                                allNodes.insert(itmp[3]);

                                if ( processVel ) {
                                    T1 s_sq = computeSlowness(mid_pt, true);
                                    s_sq *= s_sq;  // slowness squared
                                    update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                                } else {
                                    update_m_data(m_data, m, allNodes, mid_pt, ds);
                                }

                                break;
                            }
                        }
                    }
                }
            } else {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        for ( auto nc=nodes[txNode[nt]].getOwners().begin(); nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                            if ( cellNo == *nc ) {
                                r_tmp.push_back( Tx[nt] );
                                reachedTx = true;

                                std::array<T2,4> itmp = getPrimary(cellNo);
                                prev_pt = curr_pt;
                                curr_pt = Tx[nt];
                                mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                                ds = curr_pt.getDistance( prev_pt );
                                std::set<T2> allNodes;
                                allNodes.insert(itmp[0]);
                                allNodes.insert(itmp[1]);
                                allNodes.insert(itmp[2]);
                                allNodes.insert(itmp[3]);

                                if ( processVel ) {
                                    T1 s_sq = computeSlowness(mid_pt, true);
                                    s_sq *= s_sq;  // slowness squared
                                    update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                                } else {
                                    update_m_data(m_data, m, allNodes, mid_pt, ds);
                                }

                                break;
                            }
                        }
                    } else {
                        if ( cellNo == txCell[nt] ) {
                            r_tmp.push_back( Tx[nt] );

                            std::array<T2,4> itmp = getPrimary(cellNo);
                            prev_pt = curr_pt;
                            curr_pt = Tx[nt];
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            ds = curr_pt.getDistance( prev_pt );
                            std::set<T2> allNodes;
                            allNodes.insert(itmp[0]);
                            allNodes.insert(itmp[1]);
                            allNodes.insert(itmp[2]);
                            allNodes.insert(itmp[3]);

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }

                            reachedTx = true;
                        } else {
                            for ( size_t nn=0; nn<txNeighborCells[nt].size(); ++nn ) {
                                if ( cellNo == txNeighborCells[nt][nn] ) {
                                    std::array<T2,4> itmp = getPrimary(txCell[nt]);
                                    std::array<T2,3> ind[4] = {
                                        { { itmp[0], itmp[1], itmp[2] } },
                                        { { itmp[0], itmp[1], itmp[3] } },
                                        { { itmp[0], itmp[2], itmp[3] } },
                                        { { itmp[1], itmp[2], itmp[3] } }
                                    };
                                    bool found = false;
                                    for ( size_t n=0; n<4; ++n ) {
                                        std::sort( ind[n].begin(), ind[n].end() );
                                        if ( faceNodes == ind[n] ) {
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found ) {
                                        r_tmp.push_back( Tx[nt] );
                                        reachedTx = true;

                                        prev_pt = curr_pt;
                                        curr_pt = Tx[nt];
                                        mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                                        ds = curr_pt.getDistance( prev_pt );
                                        std::set<T2> allNodes;
                                        allNodes.insert(itmp[0]);
                                        allNodes.insert(itmp[1]);
                                        allNodes.insert(itmp[2]);
                                        allNodes.insert(itmp[3]);

                                        if ( processVel ) {
                                            T1 s_sq = computeSlowness(mid_pt, true);
                                            s_sq *= s_sq;  // slowness squared
                                            update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                                        } else {
                                            update_m_data(m_data, m, allNodes, mid_pt, ds);
                                        }

                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if ( reachedTx ) break;
                }
            }
        }
        // for inversion, the order should be from Tx to Rx, so we reorder...
        size_t npts = r_tmp.size();
        r_data.resize( npts );
        for ( size_t nn=0; nn<npts; ++nn ) {
            r_data[nn] = r_tmp[ npts-1-nn ];
        }

        delete grad3d;
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::getRaypath(const std::vector<sxyz<T1>>& Tx,
                                          const std::vector<T1>& t0,
                                          const sxyz<T1> &Rx,
                                          std::vector<sxyz<T1>>& r_data,
                                          std::vector<sijv<T1>>& m_data,
                                          const size_t RxNo,
                                          T1 &tt,
                                          const size_t threadNo) const {
        // Important note: m_data holds terms of traveltime derivative w/r to slowness (not velocity)
        T1 minDist = small;
        std::vector<sxyz<T1>> r_tmp;
        r_tmp.emplace_back( Rx );
        tt = 0.0;
        T1 s1=0.0, s2=0.0;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        const txInfo_t& txi = getTxInfo(Tx, threadNo);
        const std::vector<bool>& txOnNode = txi.txOnNode;
        const std::vector<bool>& txOnEdge = txi.txOnEdge;
        const std::vector<bool>& txOnFace = txi.txOnFace;
        const std::vector<T2>& txNode = txi.txNode;
        const std::vector<T2>& txCell = txi.txCell;
        const std::vector<std::array<T2,2>>& txEdges = txi.txEdges;
        const std::vector<std::array<T2,3>>& txFaces = txi.txFaces;
        const std::vector<std::vector<T2>>& txNeighborCells = txi.txNeighborCells;

        T2 cellNo=0, nodeNo=0, nodeNoPrev=0;
        sxyz<T1> curr_pt( Rx ), mid_pt, prev_pt( Rx );
        sijv<T1> m;
        m.i = RxNo;

        bool onNode = false;
        bool onEdge = false;
        bool onFace = false;
        std::array<T2,2> edgeNodes{ {0, 0} }, edgeNodesPrev;
        std::array<T2,3> faceNodes{ {0, 0, 0} }, faceNodesPrev;
        Grad3D<T1,NODE>* grad3d = nullptr;
        if ( rp_method == 0 ) {
            grad3d = new Grad3D_ls_fo<T1,NODE>();
        } else if ( rp_method == 1 ) {
            grad3d = new Grad3D_ls_so<T1,NODE>();
        } else if ( rp_method == 2 ) {
            grad3d = new Grad3D_ab<T1,NODE>();
        }
        bool reachedTx = false;

        for ( T2 nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == curr_pt ) {
                nodeNo = nn;
                onNode = true;
                s1 = nodes[nodeNo].getNodeSlowness();
                break;
            }
        }
        if ( !onNode ) {
            cellNo = getCellNo( curr_pt );

            std::array<T2,4> itmp = getPrimary(cellNo);
            T2 ind[6][2] = {
                {itmp[0], itmp[1]},
                {itmp[0], itmp[2]},
                {itmp[0], itmp[3]},
                {itmp[1], itmp[2]},
                {itmp[1], itmp[3]},
                {itmp[2], itmp[3]} };

            for ( size_t n=0; n<6; ++n ) {
                if ( areCollinear(curr_pt, nodes[ind[n][0]], nodes[ind[n][1]]) ) {
                    onEdge = true;
                    edgeNodes[0] = ind[n][0];
                    edgeNodes[1] = ind[n][1];

                    if ( processVel )
                        s1 = Interpolator<T1>::linearVel(r_tmp.back(),
                                                         nodes[edgeNodes[0]],
                                                         nodes[edgeNodes[1]]);
                    else
                        s1 = Interpolator<T1>::linear(r_tmp.back(),
                                                      nodes[edgeNodes[0]],
                                                      nodes[edgeNodes[1]]);

                    break;
                }
            }
            if ( !onEdge ) {
                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }

                for ( size_t n=0; n<4; ++n ) {
                    if ( areCoplanar(curr_pt, nodes[ind[n][0]], nodes[ind[n][1]], nodes[ind[n][2]]) ) {
                        onFace = true;
                        faceNodes[0] = ind[n][0];
                        faceNodes[1] = ind[n][1];
                        faceNodes[2] = ind[n][2];

                        if ( processVel )
                            s1 = Interpolator<T1>::bilinearTriangleVel(curr_pt,
                                                                       nodes[faceNodes[0]],
                                                                       nodes[faceNodes[1]],
                                                                       nodes[faceNodes[2]]);
                        else
                            s1 = Interpolator<T1>::bilinearTriangle(curr_pt,
                                                                    nodes[faceNodes[0]],
                                                                    nodes[faceNodes[1]],
                                                                    nodes[faceNodes[2]]);
                        //  faceNodes shoud not be assigned, face was not intersected
                        faceNodes = {0, 0, 0};
                        break;
                    }
                }
                if ( !onFace ) {
                    if ( processVel )
                        s1 = Interpolator<T1>::trilinearTriangleVel(curr_pt,
                                                                    nodes[itmp[0]],
                                                                    nodes[itmp[1]],
                                                                    nodes[itmp[2]],
                                                                    nodes[itmp[3]]);
                    else
                        s1 = Interpolator<T1>::trilinearTriangle(curr_pt,
                                                                 nodes[itmp[0]],
                                                                 nodes[itmp[1]],
                                                                 nodes[itmp[2]],
                                                                 nodes[itmp[3]]);
                }
            }
        }

        for ( size_t nt=0; nt<txCell.size(); ++nt ) {
            if ( cellNo == txCell[nt] ) {
                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( processVel )
                    s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                nodes[itmp[0]],
                                                                nodes[itmp[1]],
                                                                nodes[itmp[2]],
                                                                nodes[itmp[3]]);
                else
                    s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                             nodes[itmp[0]],
                                                             nodes[itmp[1]],
                                                             nodes[itmp[2]],
                                                             nodes[itmp[3]]);

                mid_pt = static_cast<T1>(0.5)*(curr_pt + Tx[nt]);
                T1 ds = curr_pt.getDistance( Tx[nt] );

                std::set<T2> allNodes;
                allNodes.insert( itmp[0] );
                allNodes.insert( itmp[1] );
                allNodes.insert( itmp[2] );
                allNodes.insert( itmp[3] );

                if ( processVel ) {
                    T1 s_sq = computeSlowness(mid_pt, true);
                    s_sq *= s_sq;  // slowness squared
                    update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                } else {
                    update_m_data(m_data, m, allNodes, mid_pt, ds);
                }

                tt += t0[nt] + 0.5*(s1 + s2) * ds;
                r_tmp.push_back( Tx[nt] );
                reachedTx = true;
                break;
            }
        }

        T1 ds;
        sxyz<T1> g;
        size_t max_itertions = tetrahedra.size();
        size_t itertions = 0;
        while ( reachedTx == false ) {

            if ( onNode ) {

                // find cells common to edge
                std::byte nnodes_buf[8192];
                std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                std::pmr::set<NODE*> nnodes{&nnodes_pool};
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                    getNeighborNodes(*nc, nnodes);
                }

                // compute gradient with nodes from all common cells
                if ( rp_method < 2 ) {
                    // find cells common to edge
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                        getNeighborNodes(*nc, nnodes);
                    }
                    // compute gradient with nodes from all common cells
                    g = grad3d->compute(curr_pt, nodes[nodeNo].getTT(threadNo), nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(1);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[nodeNo]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }

                nodeNoPrev = nodeNo;

                // find cell for which gradient intersect opposing face
                bool foundIntersection = false;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {

                    std::array<T2,3> nb;
                    size_t n=0;
                    for (auto nn=this->neighbors[*nc].begin(); nn!=this->neighbors[*nc].end(); ++nn ) {
                        if ( *nn != nodeNo && nodes[*nn].isPrimary() ) {
                            nb[n++] = *nn;
                        }
                    }
                    std::sort(nb.begin(), nb.end());

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle( nodeNo, g, nb[0], nb[1], nb[2], pt_i);
                    if ( !foundIntersection ) {
                        continue;
                    }

                    //  check if cell is (one of) TxCell(s)
                    for (size_t nt=0; nt<Tx.size(); ++nt) {
                        if ( *nc == txCell[nt] ) {
                            std::array<T2,4> itmp = getPrimary(*nc);

                            if ( processVel )
                                s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                            nodes[itmp[0]],
                                                                            nodes[itmp[1]],
                                                                            nodes[itmp[2]],
                                                                            nodes[itmp[3]]);
                            else
                                s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                         nodes[itmp[0]],
                                                                         nodes[itmp[1]],
                                                                         nodes[itmp[2]],
                                                                         nodes[itmp[3]]);

                            tt += t0[nt] + 0.5*(s1 + s2) * r_tmp.back().getDistance( Tx[nt] );
                            r_tmp.push_back( Tx[nt] );

                            prev_pt = curr_pt;
                            curr_pt = Tx[nt];
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            ds = curr_pt.getDistance( prev_pt );
                            std::set<T2> allNodes;
                            allNodes.insert(itmp[0]);
                            allNodes.insert(itmp[1]);
                            allNodes.insert(itmp[2]);
                            allNodes.insert(itmp[3]);

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }

                            reachedTx = true;
                            break;
                        }
                    }
                    if ( reachedTx ) {
                        break;
                    }

                    prev_pt = curr_pt;
                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, nb, onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    r_tmp.push_back( curr_pt );

                    // compute terms of matrix M
                    mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                    ds = curr_pt.getDistance( prev_pt );

                    tt += 0.5*(s1 + s2) * ds;
                    s1 = s2;

                    std::set<T2> allNodes;
                    allNodes.insert(nodeNoPrev);
                    if ( onNode ) {
                        allNodes.insert(nodeNo);
                    } else if ( onEdge ) {
                        allNodes.insert( edgeNodes[0] );
                        allNodes.insert( edgeNodes[1] );
                    } else { // onFace
                        allNodes.insert( faceNodes[0] );
                        allNodes.insert( faceNodes[1] );
                        allNodes.insert( faceNodes[2] );
                    }

                    if ( processVel ) {
                        T1 s_sq = computeSlowness(mid_pt, true);
                        s_sq *= s_sq;  // slowness squared
                        update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                    } else {
                        update_m_data(m_data, m, allNodes, mid_pt, ds);
                    }

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell1(faceNodes, nodeNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        m_data.resize(0);
                        reachedTx = true;
                    }
                    break;
                }

                if ( foundIntersection == false ) {
                    nodeNoPrev = nodeNo;
                    // projet gradient on face and find intersection on next edge
                    sxyz<T1> pt_i;
                    foundIntersection = projectOnFace(curr_pt, nodeNo, g, edgeNodes, pt_i, threadNo);

                    if ( foundIntersection == false ) {
                        std::cout << "\n\nWarning: finding raypath (onNode) failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        m_data.resize(0);
                        reachedTx = true;
                        break;
                    }
                    // find which cell we are in
                    for ( size_t nc=0; nc<tetrahedra.size(); ++nc ) {
                        std::array<T2,4> itet = {tetrahedra[nc].i[0],
                            tetrahedra[nc].i[1],
                            tetrahedra[nc].i[2],
                            tetrahedra[nc].i[3]};
                        // because we are at the surface, there is only one cell with nodeNo & edgeNodes
                        if (std::find(itet.begin(), itet.end(), nodeNo) != itet.end() &&
                            std::find(itet.begin(), itet.end(), edgeNodes[0]) != itet.end() &&
                            std::find(itet.begin(), itet.end(), edgeNodes[1]) != itet.end()) {
                            cellNo = static_cast<T2>(nc);
                            break;
                        }
                    }

                    //  check if cell is (one of) TxCell(s)
                    for (size_t nt=0; nt<Tx.size(); ++nt) {
                        if ( cellNo == txCell[nt] ) {
                            std::array<T2,4> itmp = getPrimary(cellNo);

                            if ( processVel )
                                s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                            nodes[itmp[0]],
                                                                            nodes[itmp[1]],
                                                                            nodes[itmp[2]],
                                                                            nodes[itmp[3]]);
                            else
                                s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                         nodes[itmp[0]],
                                                                         nodes[itmp[1]],
                                                                         nodes[itmp[2]],
                                                                         nodes[itmp[3]]);

                            r_tmp.push_back( Tx[nt] );

                            mid_pt = static_cast<T1>(0.5)*(curr_pt + Tx[nt]);
                            ds = curr_pt.getDistance( Tx[nt] );

                            tt += t0[nt] + 0.5*(s1 + s2) * curr_pt.getDistance( Tx[nt] );

                            std::set<T2> allNodes;
                            allNodes.insert( itmp[0] );
                            allNodes.insert( itmp[1] );
                            allNodes.insert( itmp[2] );
                            allNodes.insert( itmp[3] );

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }

                            reachedTx = true;
                            break;
                        }
                    }
                    if ( reachedTx ) {
                        break;
                    }

                    // we might be on one of the nodes
                    if ( pt_i.getDistance(nodes[edgeNodes[0]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[0];
                        pt_i = nodes[edgeNodes[0]];
                        onEdge = false;
                    } else if ( pt_i.getDistance(nodes[edgeNodes[1]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[1];
                        pt_i = nodes[edgeNodes[1]];
                        onEdge = false;
                    } else {
                        onEdge = true;
                        onNode = false;
                    }

                    prev_pt = curr_pt;
                    curr_pt = pt_i;

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    r_tmp.push_back( curr_pt );

                    mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                    ds = curr_pt.getDistance( prev_pt );

                    tt += 0.5*(s1 + s2) * ds;
                    s1 = s2;

                    std::set<T2> allNodes;
                    allNodes.insert( nodeNoPrev );
                    allNodes.insert( edgeNodes[0] );
                    allNodes.insert( edgeNodes[1] );

                    if ( processVel ) {
                        T1 s_sq = computeSlowness(mid_pt, true);
                        s_sq *= s_sq;  // slowness squared
                        update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                    } else {
                        update_m_data(m_data, m, allNodes, mid_pt, ds);
                    }
                }

            } else if ( onEdge ) {

                // find cells common to edge
                std::vector<T2> cells;
                std::byte nnodes_buf[8192];
                std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                std::pmr::set<NODE*> nnodes{&nnodes_pool};
                for ( auto nc0=nodes[edgeNodes[0]].getOwners().begin(); nc0!=nodes[edgeNodes[0]].getOwners().end(); ++nc0 ) {
                    if ( std::find(nodes[edgeNodes[1]].getOwners().begin(), nodes[edgeNodes[1]].getOwners().end(), *nc0)!=nodes[edgeNodes[1]].getOwners().end() ) {
                        cells.push_back( *nc0 );
                        getNeighborNodes(*nc0, nnodes);
                    }
                }
                if ( rp_method < 2 ) {
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    for (size_t n=0; n<cells.size(); ++n ) {
                        getNeighborNodes(cells[n], nnodes);
                    }
                    T1 d01 = nodes[edgeNodes[0]].getDistance(nodes[edgeNodes[1]]);
                    T1 w0 = curr_pt.getDistance(nodes[edgeNodes[1]]) / d01;
                    T1 w1 = curr_pt.getDistance(nodes[edgeNodes[0]]) / d01;
                    T1 curr_t = nodes[edgeNodes[0]].getTT(threadNo)*w0 + nodes[edgeNodes[1]].getTT(threadNo)*w1;
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(2);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[edgeNodes[0]]);
                    ref_pt[1] = &(nodes[edgeNodes[1]]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, edgeNodes, Tx, txCell);

                edgeNodesPrev = edgeNodes;

                bool foundIntersection = false;
                for (size_t n=0; n<cells.size(); ++n ) {

                    T2 testCellNo = cells[n];

                    // there are 2 faces that might be intersected
                    std::array<T2,2> edgeNodes2;
                    size_t n2=0;
                    for ( auto nn=this->neighbors[testCellNo].begin(); nn!= this->neighbors[testCellNo].end(); ++nn ) {
                        if ( *nn!=edgeNodes[0] && *nn!=edgeNodes[1] && nodes[*nn].isPrimary() ) {
                            edgeNodes2[n2++] = *nn;
                        }
                    }

                    sxyz<T1> pt_i;
                    T2 itmpNode;
                    foundIntersection = intersectVecTriangle(curr_pt, g,
                                                             edgeNodes[0],
                                                             edgeNodes2[0],
                                                             edgeNodes2[1], pt_i);
                    itmpNode = edgeNodes[0];
                    if ( !foundIntersection ) {
                        foundIntersection = intersectVecTriangle(curr_pt, g,
                                                                 edgeNodes[1],
                                                                 edgeNodes2[0],
                                                                 edgeNodes2[1], pt_i);
                        itmpNode = edgeNodes[1];
                    }
                    if ( !foundIntersection ) {
                        continue;
                    }

                    cellNo = testCellNo;
                    prev_pt = curr_pt;
                    curr_pt = pt_i;
                    bool break_flag = check_pt_location(curr_pt, this->neighbors[cellNo],
                                                        {itmpNode, edgeNodes2[0], edgeNodes2[1]},
                                                        onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    r_tmp.push_back( curr_pt );

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    // compute terms of matrix M
                    mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                    ds = curr_pt.getDistance( prev_pt );
                    tt += 0.5*(s1 + s2) * ds;
                    s1 = s2;

                    std::set<T2> allNodes;
                    allNodes.insert( edgeNodesPrev[0] );
                    allNodes.insert( edgeNodesPrev[1] );
                    if ( onNode ) {
                        allNodes.insert(nodeNo);
                    } else if ( onEdge ) {
                        allNodes.insert( edgeNodes[0] );
                        allNodes.insert( edgeNodes[1] );
                    } else { // onFace
                        allNodes.insert( faceNodes[0] );
                        allNodes.insert( faceNodes[1] );
                        allNodes.insert( faceNodes[2] );
                    }

                    if ( processVel ) {
                        T1 s_sq = computeSlowness(mid_pt, true);
                        s_sq *= s_sq;  // slowness squared
                        update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                    } else {
                        update_m_data(m_data, m, allNodes, mid_pt, ds);
                    }

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        m_data.resize(0);
                        reachedTx = true;
                    }
                    break;
                }
                if ( foundIntersection == false ) {

                    edgeNodesPrev = edgeNodes;
                    sxyz<T1> pt_i;
                    g = projectOnFace(curr_pt, g, edgeNodes, cells, pt_i, threadNo);
                    if ( g.x==0.0 && g.y==0.0 && g.z==0.0 ) {
                        foundIntersection = false;
                    } else {
                        foundIntersection = true;
                    }

                    if ( foundIntersection == false || curr_pt == pt_i ) {
                        std::cout << "\n\nWarning: finding raypath (onEdge) failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        m_data.resize(0);
                        reachedTx = true;
                        break;
                    }

                    // we might be on one of the nodes
                    if ( pt_i.getDistance(nodes[edgeNodes[0]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[0];
                        pt_i = nodes[edgeNodes[0]];
                        onEdge = false;
                    } else if ( pt_i.getDistance(nodes[edgeNodes[1]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[1];
                        pt_i = nodes[edgeNodes[1]];
                        onEdge = false;
                    }

                    prev_pt = curr_pt;
                    curr_pt = pt_i;

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    tt += 0.5*(s1 + s2) * r_tmp.back().getDistance( curr_pt );
                    s1 = s2;
                    r_tmp.push_back( curr_pt );

                    mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                    ds = curr_pt.getDistance( prev_pt );
                    std::set<T2> allNodes;

                    allNodes.insert( edgeNodesPrev[0] );
                    allNodes.insert( edgeNodesPrev[1] );
                    allNodes.insert( edgeNodes[0] );
                    allNodes.insert( edgeNodes[1] );

                    if ( processVel ) {
                        T1 s_sq = computeSlowness(mid_pt, true);
                        s_sq *= s_sq;  // slowness squared
                        update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                    } else {
                        update_m_data(m_data, m, allNodes, mid_pt, ds);
                    }
                }

            } else if ( onFace ) { // on Face

                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( rp_method < 2 ) {
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    getNeighborNodes(cellNo, nnodes);
                    T1 curr_t;
                    if ( r_tmp.size() <= 1 ) {
                        curr_t = Interpolator<T1>::trilinearTime(curr_pt,
                                                                 nodes[itmp[0]],
                                                                 nodes[itmp[1]],
                                                                 nodes[itmp[2]],
                                                                 nodes[itmp[3]],
                                                                 threadNo);
                    } else {
                        curr_t = Interpolator<T1>::bilinearTime(curr_pt,
                                                                nodes[faceNodes[0]],
                                                                nodes[faceNodes[1]],
                                                                nodes[faceNodes[2]],
                                                                threadNo);
                    }
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(3);
                    if ( r_tmp.size() <= 1 ) {
                        ref_pt[0] = &(nodes[itmp[0]]);
                        ref_pt[1] = &(nodes[itmp[1]]);
                        ref_pt[2] = &(nodes[itmp[2]]);
                        ref_pt.push_back( &(nodes[itmp[3]]) );
                    } else {
                        ref_pt[0] = &(nodes[faceNodes[0]]);
                        ref_pt[1] = &(nodes[faceNodes[1]]);
                        ref_pt[2] = &(nodes[faceNodes[2]]);
                    }
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }
                // there are 3 faces that might be intersected

                faceNodesPrev = faceNodes;

                bool foundIntersection = false;
                for ( size_t n=0; n<4; ++n ) {
                    if ( ind[n] == faceNodes ) continue;

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);

                    if ( !foundIntersection )
                        continue;

                    prev_pt = curr_pt;
                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    r_tmp.push_back( curr_pt );

                    // compute terms of matrix M
                    mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                    ds = curr_pt.getDistance( prev_pt );

                    tt += 0.5*(s1 + s2) * ds;
                    s1 = s2;

                    std::set<T2> allNodes;
                    allNodes.insert( faceNodesPrev[0] );
                    allNodes.insert( faceNodesPrev[1] );
                    allNodes.insert( faceNodesPrev[2] );
                    if ( onNode ) {
                        allNodes.insert(nodeNo);
                    } else if ( onEdge ) {
                        allNodes.insert( edgeNodes[0] );
                        allNodes.insert( edgeNodes[1] );
                    } else { // onFace
                        allNodes.insert( faceNodes[0] );
                        allNodes.insert( faceNodes[1] );
                        allNodes.insert( faceNodes[2] );
                    }

                    if ( processVel ) {
                        T1 s_sq = computeSlowness(mid_pt, true);
                        s_sq *= s_sq;  // slowness squared
                        update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                    } else {
                        update_m_data(m_data, m, allNodes, mid_pt, ds);
                    }

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        m_data.resize(0);
                        reachedTx = true;
                    }
                    break;
                }

                if ( foundIntersection == false ) {

                    // we must be on an face with gradient pointing slightly outward tetrahedron
                    // return in other cell but keep gradient
                    cellNo = findAdjacentCell2(faceNodes, cellNo);

                    std::array<T2,4> itmp = getPrimary(cellNo);
                    ind[0] = { { itmp[0], itmp[1], itmp[2] } };
                    ind[1] = { { itmp[0], itmp[1], itmp[3] } };
                    ind[2] = { { itmp[0], itmp[2], itmp[3] } };
                    ind[3] = { { itmp[1], itmp[2], itmp[3] } };

                    for ( size_t n=0; n<4; ++n ) {
                        std::sort( ind[n].begin(), ind[n].end() );
                    }

                    for ( size_t n=0; n<4; ++n ) {
                        if ( ind[n] == faceNodes ) continue;

                        sxyz<T1> pt_i;
                        foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                                 ind[n][1], ind[n][2],
                                                                 pt_i);

                        if ( !foundIntersection ) {
                            continue;
                        }
                        prev_pt = curr_pt;
                        curr_pt = pt_i;

                        bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                            nodeNo, onEdge, edgeNodes,
                                                            onFace, faceNodes);

                        s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                             faceNodes);

                        r_tmp.push_back( curr_pt );

                        mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                        ds = curr_pt.getDistance( prev_pt );

                        tt += 0.5*(s1 + s2) * ds;
                        s1 = s2;

                        std::set<T2> allNodes;
                        allNodes.insert( faceNodesPrev[0] );
                        allNodes.insert( faceNodesPrev[1] );
                        allNodes.insert( faceNodesPrev[2] );
                        if ( onNode ) {
                            allNodes.insert(nodeNo);
                        } else if ( onEdge ) {
                            allNodes.insert( edgeNodes[0] );
                            allNodes.insert( edgeNodes[1] );
                        } else { // onFace
                            allNodes.insert( faceNodes[0] );
                            allNodes.insert( faceNodes[1] );
                            allNodes.insert( faceNodes[2] );
                        }

                        if ( processVel ) {
                            T1 s_sq = computeSlowness(mid_pt, true);
                            s_sq *= s_sq;  // slowness squared
                            update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                        } else {
                            update_m_data(m_data, m, allNodes, mid_pt, ds);
                        }

                        if ( break_flag ) break;

                        // find next cell
                        cellNo = findAdjacentCell2(faceNodes, cellNo);
                        if ( cellNo == std::numeric_limits<T2>::max() ) {
                            std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                            << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                            tt = 0.0;
                            r_tmp.resize(1);
                            r_tmp[0] = Rx;
                            m_data.resize(0);
                            reachedTx = true;
                        }
                        break;
                    }
                }
                if ( foundIntersection == false ) {
                    // we must be going outside the mesh
                    // hack: project gradient on face and continue
                    faceNodesPrev = faceNodes;

                    g = projectOnFace(g, faceNodes);

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecEdge(curr_pt, g, faceNodes, pt_i, edgeNodes);

                    if ( foundIntersection == false ) {
                        std::cout << "\n\nWarning: finding raypath (onFace) failed to converge for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        m_data.resize(0);
                        reachedTx = true;
                        break;
                    }

                    // we might be on one of the nodes
                    if ( pt_i.getDistance(nodes[edgeNodes[0]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[0];
                        pt_i = nodes[edgeNodes[0]];
                        onEdge = false;
                        onFace = false;
                    } else if ( pt_i.getDistance(nodes[edgeNodes[1]]) < min_dist ) {
                        onNode = true;
                        nodeNo = edgeNodes[1];
                        pt_i = nodes[edgeNodes[1]];
                        onEdge = false;
                        onFace = false;
                    } else {
                        onEdge = true;
                        onFace = false;
                    }

                    prev_pt = curr_pt;
                    curr_pt = pt_i;

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    r_tmp.push_back( curr_pt );

                    mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                    ds = curr_pt.getDistance( prev_pt );

                    tt += 0.5*(s1 + s2) * ds;
                    s1 = s2;

                    std::set<T2> allNodes;

                    allNodes.insert( faceNodesPrev[0] );
                    allNodes.insert( faceNodesPrev[1] );
                    allNodes.insert( faceNodesPrev[2] );
                    allNodes.insert( edgeNodes[0] );
                    allNodes.insert( edgeNodes[1] );

                    if ( processVel ) {
                        T1 s_sq = computeSlowness(mid_pt, true);
                        s_sq *= s_sq;  // slowness squared
                        update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                    } else {
                        update_m_data(m_data, m, allNodes, mid_pt, ds);
                    }

                }
            } else { // at Rx, somewhere in a tetrahedron

                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( rp_method < 2 ) {
                    std::byte nnodes_buf[8192];
                    std::pmr::monotonic_buffer_resource nnodes_pool{nnodes_buf, sizeof(nnodes_buf)};
                    std::pmr::set<NODE*> nnodes{&nnodes_pool};
                    getNeighborNodes(cellNo, nnodes);
                    T1 curr_t = Interpolator<T1>::trilinearTime(curr_pt,
                                                                nodes[itmp[0]],
                                                                nodes[itmp[1]],
                                                                nodes[itmp[2]],
                                                                nodes[itmp[3]],
                                                                threadNo);
                    g = grad3d->compute(curr_pt, curr_t, nnodes, threadNo);
                } else {
                    std::vector<NODE*> ref_pt(4);
                    std::vector<std::vector<std::array<NODE*,3>>> opp_pts;
                    ref_pt[0] = &(nodes[itmp[0]]);
                    ref_pt[1] = &(nodes[itmp[1]]);
                    ref_pt[2] = &(nodes[itmp[2]]);
                    ref_pt[3] = &(nodes[itmp[3]]);
                    getNeighborNodesAB(ref_pt, opp_pts);
                    g = dynamic_cast<Grad3D_ab<T1,NODE>*>(grad3d)->compute(curr_pt, ref_pt, opp_pts, threadNo);
                }
                checkCloseToTx(curr_pt, g, cellNo, Tx, txCell);

                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };
                for ( size_t n=0; n<4; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }
                // there are 4 faces that might be intersected

                bool foundIntersection = false;
                for ( size_t n=0; n<4; ++n ) {

                    sxyz<T1> pt_i;
                    foundIntersection = intersectVecTriangle(curr_pt, g, ind[n][0],
                                                             ind[n][1], ind[n][2],
                                                             pt_i);

                    if ( !foundIntersection )
                        continue;

                    prev_pt = curr_pt;
                    curr_pt = pt_i;

                    bool break_flag = check_pt_location(curr_pt, ind[n], onNode,
                                                        nodeNo, onEdge, edgeNodes,
                                                        onFace, faceNodes);

                    s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                         faceNodes);

                    r_tmp.push_back( curr_pt );

                    mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                    ds = curr_pt.getDistance( prev_pt );

                    tt += 0.5*(s1 + s2) * ds;
                    s1 = s2;

                    std::set<T2> allNodes;
                    allNodes.insert( tetrahedra[cellNo].i[0] );
                    allNodes.insert( tetrahedra[cellNo].i[1] );
                    allNodes.insert( tetrahedra[cellNo].i[2] );
                    allNodes.insert( tetrahedra[cellNo].i[3] );

                    if ( processVel ) {
                        T1 s_sq = computeSlowness(mid_pt, true);
                        s_sq *= s_sq;  // slowness squared
                        update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                    } else {
                        update_m_data(m_data, m, allNodes, mid_pt, ds);
                    }

                    if ( break_flag ) break;

                    // find next cell
                    cellNo = findAdjacentCell2(faceNodes, cellNo);
                    if ( cellNo == std::numeric_limits<T2>::max() ) {
                        std::cout << "\n\nWarning: finding raypath failed to converge (cell not found) for Rx "
                        << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                        tt = 0.0;
                        r_tmp.resize(1);
                        r_tmp[0] = Rx;
                        m_data.resize(0);
                        reachedTx = true;
                    }
                    break;
                }
                if ( foundIntersection == false ) {
                    std::cout << "\n\nWarning: finding raypath within cell failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    tt = 0.0;
                    r_data.resize(1);
                    r_data[0] = Rx;
                    m_data.resize(0);
                    reachedTx = true;
                }
            }

            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        if ( curr_pt.getDistance( Tx[nt] ) < minDist ) {
                            tt += t0[nt];
                            reachedTx = true;
                            break;
                        }
                    } else if ( txOnEdge[nt] ) {
                        if ( curr_pt.getDistance(nodes[txEdges[nt][0]]) < minDist ||
                            curr_pt.getDistance(nodes[txEdges[nt][1]]) < minDist ) {

                            if ( processVel )
                                s2 = Interpolator<T1>::linearVel(Tx[nt],
                                                                 nodes[txEdges[nt][0]],
                                                                 nodes[txEdges[nt][1]]);
                            else
                                s2 = Interpolator<T1>::linear(Tx[nt],
                                                              nodes[txEdges[nt][0]],
                                                              nodes[txEdges[nt][1]]);

                            r_tmp.push_back(Tx[nt]);
                            reachedTx = true;

                            prev_pt = curr_pt;
                            curr_pt = Tx[nt];
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            ds = curr_pt.getDistance( prev_pt );
                            tt += t0[nt] + 0.5*(s1 + s2) * ds;

                            std::set<T2> allNodes;
                            allNodes.insert(txEdges[nt][0]);
                            allNodes.insert(txEdges[nt][1]);
                            allNodes.insert(nodeNo);

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }

                            break;
                        }
                    } else if ( txOnFace[nt] ) {
                        if ( curr_pt.getDistance(nodes[txFaces[nt][0]]) < minDist ||
                            curr_pt.getDistance(nodes[txFaces[nt][1]]) < minDist ||
                            curr_pt.getDistance(nodes[txFaces[nt][2]]) < minDist ) {

                            if ( processVel )
                                s2 = Interpolator<T1>::bilinearTriangleVel(Tx[nt],
                                                                           nodes[txFaces[nt][0]],
                                                                           nodes[txFaces[nt][1]],
                                                                           nodes[txFaces[nt][2]]);
                            else
                                s2 = Interpolator<T1>::bilinearTriangle(Tx[nt],
                                                                        nodes[txFaces[nt][0]],
                                                                        nodes[txFaces[nt][1]],
                                                                        nodes[txFaces[nt][2]]);

                            r_tmp.push_back(Tx[nt]);
                            reachedTx = true;

                            prev_pt = curr_pt;
                            curr_pt = Tx[nt];
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            ds = curr_pt.getDistance( prev_pt );
                            tt += t0[nt] + 0.5*(s1 + s2) * ds;

                            std::set<T2> allNodes;
                            allNodes.insert(txFaces[nt][0]);
                            allNodes.insert(txFaces[nt][1]);
                            allNodes.insert(txFaces[nt][2]);
                            allNodes.insert(nodeNo);

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }

                            break;
                        }
                    }
                }
            } if ( onEdge ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        if ( txNode[nt] == edgeNodes[0] || txNode[nt] == edgeNodes[1] ) {
                            s2 = nodes[txNode[nt]].getNodeSlowness();
                            r_tmp.push_back( Tx[nt] );
                            reachedTx = true;

                            prev_pt = curr_pt;
                            curr_pt = Tx[nt];
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            ds = curr_pt.getDistance( prev_pt );
                            tt += t0[nt] + 0.5*(s1 + s2) * ds;

                            std::set<T2> allNodes;
                            allNodes.insert(edgeNodes[0]);
                            allNodes.insert(edgeNodes[1]);

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }

                            break;
                        }
                    } else {
                        std::array<T2,4> itmp = getPrimary(txCell[nt]);
                        // find adjacent cells
                        const T2 ind[6][2] = {
                            {itmp[0], itmp[1]},
                            {itmp[0], itmp[2]},
                            {itmp[0], itmp[3]},
                            {itmp[1], itmp[2]},
                            {itmp[1], itmp[3]},
                            {itmp[2], itmp[3]} };
                        for ( size_t ne=0; ne<6; ++ne ) {
                            if ( (ind[ne][0] == edgeNodes[0] && ind[ne][1] == edgeNodes[1]) ||
                                (ind[ne][0] == edgeNodes[1] && ind[ne][1] == edgeNodes[0]) ) {
                                if ( processVel )
                                    s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                nodes[itmp[0]],
                                                                                nodes[itmp[1]],
                                                                                nodes[itmp[2]],
                                                                                nodes[itmp[3]]);
                                else
                                    s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                             nodes[itmp[0]],
                                                                             nodes[itmp[1]],
                                                                             nodes[itmp[2]],
                                                                             nodes[itmp[3]]);
                                r_tmp.push_back( Tx[nt] );
                                reachedTx = true;

                                prev_pt = curr_pt;
                                curr_pt = Tx[nt];
                                mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                                ds = curr_pt.getDistance( prev_pt );
                                tt += t0[nt] + 0.5*(s1 + s2) * ds;

                                std::set<T2> allNodes;
                                allNodes.insert(itmp[0]);
                                allNodes.insert(itmp[1]);
                                allNodes.insert(itmp[2]);
                                allNodes.insert(itmp[3]);

                                if ( processVel ) {
                                    T1 s_sq = computeSlowness(mid_pt, true);
                                    s_sq *= s_sq;  // slowness squared
                                    update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                                } else {
                                    update_m_data(m_data, m, allNodes, mid_pt, ds);
                                }

                                break;
                            }
                        }
                    }
                }
            } else {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        for ( auto nc=nodes[txNode[nt]].getOwners().begin(); nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                            if ( cellNo == *nc ) {
                                std::array<T2,4> itmp = getPrimary(cellNo);

                                if ( processVel )
                                    s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                nodes[itmp[0]],
                                                                                nodes[itmp[1]],
                                                                                nodes[itmp[2]],
                                                                                nodes[itmp[3]]);
                                else
                                    s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                             nodes[itmp[0]],
                                                                             nodes[itmp[1]],
                                                                             nodes[itmp[2]],
                                                                             nodes[itmp[3]]);
                                r_tmp.push_back( Tx[nt] );
                                reachedTx = true;

                                prev_pt = curr_pt;
                                curr_pt = Tx[nt];
                                mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                                ds = curr_pt.getDistance( prev_pt );
                                tt += t0[nt] + 0.5*(s1 + s2) * ds;

                                std::set<T2> allNodes;
                                allNodes.insert(itmp[0]);
                                allNodes.insert(itmp[1]);
                                allNodes.insert(itmp[2]);
                                allNodes.insert(itmp[3]);

                                if ( processVel ) {
                                    T1 s_sq = computeSlowness(mid_pt, true);
                                    s_sq *= s_sq;  // slowness squared
                                    update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                                } else {
                                    update_m_data(m_data, m, allNodes, mid_pt, ds);
                                }

                                break;
                            }
                        }
                    } else {
                        if ( cellNo == txCell[nt] ) {
                            std::array<T2,4> itmp = getPrimary(cellNo);
                            if ( processVel )
                                s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                            nodes[itmp[0]],
                                                                            nodes[itmp[1]],
                                                                            nodes[itmp[2]],
                                                                            nodes[itmp[3]]);
                            else
                                s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                         nodes[itmp[0]],
                                                                         nodes[itmp[1]],
                                                                         nodes[itmp[2]],
                                                                         nodes[itmp[3]]);
                            r_tmp.push_back( Tx[nt] );
                            reachedTx = true;

                            prev_pt = curr_pt;
                            curr_pt = Tx[nt];
                            mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                            ds = curr_pt.getDistance( prev_pt );
                            tt += t0[nt] + 0.5*(s1 + s2) * ds;

                            std::set<T2> allNodes;
                            allNodes.insert(itmp[0]);
                            allNodes.insert(itmp[1]);
                            allNodes.insert(itmp[2]);
                            allNodes.insert(itmp[3]);

                            if ( processVel ) {
                                T1 s_sq = computeSlowness(mid_pt, true);
                                s_sq *= s_sq;  // slowness squared
                                update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                            } else {
                                update_m_data(m_data, m, allNodes, mid_pt, ds);
                            }

                        } else {
                            for ( size_t nn=0; nn<txNeighborCells[nt].size(); ++nn ) {
                                if ( cellNo == txNeighborCells[nt][nn] ) {
                                    std::array<T2,4> itmp = getPrimary(txCell[nt]);
                                    std::array<T2,3> ind[4] = {
                                        { { itmp[0], itmp[1], itmp[2] } },
                                        { { itmp[0], itmp[1], itmp[3] } },
                                        { { itmp[0], itmp[2], itmp[3] } },
                                        { { itmp[1], itmp[2], itmp[3] } }
                                    };
                                    bool found = false;
                                    for ( size_t n=0; n<4; ++n ) {
                                        std::sort( ind[n].begin(), ind[n].end() );
                                        if ( faceNodes == ind[n] ) {
                                            found = true;
                                            break;
                                        }
                                    }
                                    if ( found ) {
                                        itmp = getPrimary(cellNo);
                                        if ( processVel )
                                            s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                        nodes[itmp[0]],
                                                                                        nodes[itmp[1]],
                                                                                        nodes[itmp[2]],
                                                                                        nodes[itmp[3]]);
                                        else
                                            s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                                     nodes[itmp[0]],
                                                                                     nodes[itmp[1]],
                                                                                     nodes[itmp[2]],
                                                                                     nodes[itmp[3]]);
                                        r_tmp.push_back( Tx[nt] );
                                        reachedTx = true;

                                        prev_pt = curr_pt;
                                        curr_pt = Tx[nt];
                                        mid_pt = static_cast<T1>(0.5)*(curr_pt + prev_pt);
                                        ds = curr_pt.getDistance( prev_pt );
                                        tt += t0[nt] + 0.5*(s1 + s2) * ds;

                                        std::set<T2> allNodes;
                                        allNodes.insert(itmp[0]);
                                        allNodes.insert(itmp[1]);
                                        allNodes.insert(itmp[2]);
                                        allNodes.insert(itmp[3]);

                                        if ( processVel ) {
                                            T1 s_sq = computeSlowness(mid_pt, true);
                                            s_sq *= s_sq;  // slowness squared
                                            update_m_data(m_data, m, allNodes, mid_pt, ds, s_sq);
                                        } else {
                                            update_m_data(m_data, m, allNodes, mid_pt, ds);
                                        }

                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if ( reachedTx ) break;
                }
            }
            itertions ++;
            if (itertions > max_itertions){
                std::cout << "\n\nWarning: raypath failed to converge for Rx : "
                << Rx.x << ' ' << Rx.y << ' ' << Rx.z <<" (infinite loop risk)"<<std::endl;
                tt = 0.0;
                r_tmp.resize(1);
                r_tmp[0] = Rx;
                m_data.resize(0);
                break; // break while loop
            }
        }
        // for inversion, the order should be from Tx to Rx, so we reorder...
        size_t npts = r_tmp.size();
        r_data.resize( npts );
        for ( size_t nn=0; nn<npts; ++nn ) {
            r_data[nn] = r_tmp[ npts-1-nn ];
        }

        delete grad3d;
    }


    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dun<T1,T2,NODE>::getTraveltime_blti(const std::vector<sxyz<T1>>& Tx,
                                                const std::vector<T1>& t0,
                                                const sxyz<T1> &Rx,
                                                const size_t threadNo) const {
        T1 minDist = min_dist;
        T1 tt = 0.0;
        T1 tt_s1, tt_s2;


        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                return t0[ns];
            }
        }

        const txInfo_t& txi = getTxInfoBlti(Tx, threadNo);
        const std::vector<bool>& txOnNode = txi.txOnNode;
        const std::vector<T2>& txNode = txi.txNode;
        const std::vector<T2>& txCell = txi.txCell;
        const std::vector<std::vector<T2>>& txNeighborCells = txi.txNeighborCells;

        T2 cellNo, nodeNo;
        sxyz<T1> curr_pt( Rx ), prev_pt( Rx );
        bool onNode = false;
        bool onEdge = false;
        bool onFace = false;
        std::array<T2,2> edgeNodes;
        std::array<T2,3> faceNodes={{0,0,0}};
        std::array<T2,3> faceNodesPrev={{0,0,0}};
        bool reachedTx = false;
        T1 Xmax=-std::numeric_limits<T1>::max();
        T1 Xmin=std::numeric_limits<T1>::max();
        T1 Ymax=-std::numeric_limits<T1>::max();
        T1 Ymin=std::numeric_limits<T1>::max();
        T1 Zmax=-std::numeric_limits<T1>::max();
        T1 Zmin=std::numeric_limits<T1>::max();
        for ( T2 nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == curr_pt ) {
                nodeNo = nn;
                onNode = true;
                tt_s1 = nodes[nodeNo].getNodeSlowness();
            }
            Xmax=nodes[nn].getX()>Xmax?nodes[nn].getX():Xmax;
            Ymax=nodes[nn].getY()>Ymax?nodes[nn].getY():Ymax;
            Zmax=nodes[nn].getZ()>Zmax?nodes[nn].getZ():Zmax;

            Xmin=nodes[nn].getX()<Xmin?nodes[nn].getX():Xmin;
            Ymin=nodes[nn].getY()<Ymin?nodes[nn].getY():Ymin;
            Zmin=nodes[nn].getZ()<Zmin?nodes[nn].getZ():Zmin;
        }
        if ( !onNode ) {
            cellNo = getCellNo( curr_pt );

            std::array<T2,4> itmp = getPrimary(cellNo);
            T2 ind[6][2] = {
                {itmp[0], itmp[1]},
                {itmp[0], itmp[2]},
                {itmp[0], itmp[3]},
                {itmp[1], itmp[2]},
                {itmp[1], itmp[3]},
                {itmp[2], itmp[3]} };

            for ( size_t n=0; n<6; ++n ) {
                if ( areCollinear(curr_pt, nodes[ind[n][0]], nodes[ind[n][1]]) ) {
                    onEdge = true;
                    edgeNodes[0] = ind[n][0];
                    edgeNodes[1] = ind[n][1];

                    if ( processVel )
                        tt_s1 = Interpolator<T1>::linearVel(curr_pt,
                                                            nodes[edgeNodes[0]],
                                                            nodes[edgeNodes[1]]);
                    else
                        tt_s1 = Interpolator<T1>::linear(curr_pt,
                                                         nodes[edgeNodes[0]],
                                                         nodes[edgeNodes[1]]);

                    break;
                }
            }
            if ( !onEdge ) {
                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };

                for ( size_t n=0; n<4; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }

                for ( size_t n=0; n<4; ++n ) {
                    if ( areCoplanar(curr_pt, nodes[ind[n][0]], nodes[ind[n][1]], nodes[ind[n][2]]) ) {
                        onFace = true;
                        faceNodes[0] = ind[n][0];
                        faceNodes[1] = ind[n][1];
                        faceNodes[2] = ind[n][2];

                        if ( processVel )
                            tt_s1 = Interpolator<T1>::bilinearTriangleVel(curr_pt,
                                                                          nodes[faceNodes[0]],
                                                                          nodes[faceNodes[1]],
                                                                          nodes[faceNodes[2]]);
                        else
                            tt_s1 = Interpolator<T1>::bilinearTriangle(curr_pt,
                                                                       nodes[faceNodes[0]],
                                                                       nodes[faceNodes[1]],
                                                                       nodes[faceNodes[2]]);
                        break;
                    }
                }
            }
        }

        if (!onEdge && ! onNode && ! onFace){
            std::array<T2,4> itmp = getPrimary(cellNo);
            if ( processVel )
                tt_s1 = Interpolator<T1>::trilinearTriangleVel(curr_pt,
                                                               nodes[itmp[0]],
                                                               nodes[itmp[1]],
                                                               nodes[itmp[2]],
                                                               nodes[itmp[3]]);
            else
                tt_s1 = Interpolator<T1>::trilinearTriangle(curr_pt,
                                                            nodes[itmp[0]],
                                                            nodes[itmp[1]],
                                                            nodes[itmp[2]],
                                                            nodes[itmp[3]]);
        }
        for(auto nt=0;nt<txCell.size();++nt){
            if (getCellNo( Rx )==txCell[nt]){
                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( processVel )
                    tt_s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                   nodes[itmp[0]],
                                                                   nodes[itmp[1]],
                                                                   nodes[itmp[2]],
                                                                   nodes[itmp[3]]);
                else
                    tt_s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                nodes[itmp[0]],
                                                                nodes[itmp[1]],
                                                                nodes[itmp[2]],
                                                                nodes[itmp[3]]);

                tt += t0[nt] + 0.5*(tt_s1 + tt_s2) * curr_pt.getDistance( Tx[nt] );
                reachedTx=true;
                break;
            }
        }

        T1 time=std::numeric_limits<T1>::max();
        bool inlimitD=false;
        while ( reachedTx == false ) {
            sxyz<T1> NodeSource;
            bool NearSource=false;
            for(size_t nt=0;nt<Tx.size();++nt){
                for(auto n=this->neighbors[txCell[nt]].begin();n!=this->neighbors[txCell[nt]].begin()+4;++n){
                    for(auto nc=nodes[*n].getOwners().begin();nc!=nodes[*n].getOwners().end();++nc){
                        if(*nc==cellNo){
                            NearSource=true;
                            NodeSource=Tx[nt];
                            break;
                        }
                    }
                    if (NearSource)
                        break;
                }
                if (curr_pt.getDistance(Tx[nt])<=source_radius){
                    NearSource=true;
                    NodeSource=Tx[nt];
                    break;
                }

            }

            if ( onNode ) {

                T1 t_i=std::numeric_limits<T1>::max();
                T1 Slow;
                sxyz<T1> pt_i;
                T2 newNode;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                    T2 cellNoi = *nc;
                    std::array<T2, 3>ind;
                    for (T2 i=0;i<4;++i){
                        if (this->neighbors[cellNoi][i]==nodeNo){
                            ind[0]=this->neighbors[cellNoi][(i+1)%4];
                            ind[1]=this->neighbors[cellNoi][(i+2)%4];
                            ind[2]=this->neighbors[cellNoi][(i+3)%4];
                            break;
                        }
                    }
                    std::sort(ind.begin(), ind.end());
                    if (NearSource){
                        bool flag=false;
                        std::array<T1,3> Barycenter;
                        if (blti_solver_around_source(NodeSource, curr_pt, ind, Barycenter)==true){
                            pt_i.x=Barycenter[0]*nodes[ind[0]].getX()+Barycenter[1]*nodes[ind[1]].getX()+Barycenter[2]*nodes[ind[2]].getX();
                            pt_i.y=Barycenter[0]*nodes[ind[0]].getY()+Barycenter[1]*nodes[ind[1]].getY()+Barycenter[2]*nodes[ind[2]].getY();
                            pt_i.z=Barycenter[0]*nodes[ind[0]].getZ()+Barycenter[1]*nodes[ind[1]].getZ()+Barycenter[2]*nodes[ind[2]].getZ();
                            if(NodeSource.getDistance(pt_i)>NodeSource.getDistance(curr_pt))
                                continue;
                            for(T2 ni=0;ni<3;++ni){
                                if(std::abs(1.0-Barycenter[ni])<minDist*minDist){
                                    onNode = true;
                                    onEdge = false;
                                    onFace = false;
                                    newNode=ind[ni];
                                    cellNo=cellNoi;
                                    flag=true;
                                    break;
                                }
                            }
                            if (flag)
                                break;
                            for(T2 ni=0;ni<3;++ni){
                                if(std::abs(Barycenter[ni])<minDist*minDist){
                                    onNode = false;
                                    onEdge = true;
                                    onFace = false;
                                    edgeNodes[0] = ind[(ni+1)%3];
                                    edgeNodes[1] = ind[(ni+2)%3];
                                    cellNo=cellNoi;
                                    flag=true;
                                    break;
                                }
                            }
                            if (flag)
                                break;

                            onNode = false;
                            onEdge = false;
                            onFace = true;
                            faceNodesPrev=faceNodes;
                            faceNodes=ind;
                            cellNo=cellNoi;
                            break;
                        }
                        continue;
                    }
                    if ( ind == faceNodesPrev) continue;
                    T1 s=0.25*(nodes[this->neighbors[cellNoi][0]].getNodeSlowness()+
                               nodes[this->neighbors[cellNoi][1]].getNodeSlowness()+
                               nodes[this->neighbors[cellNoi][2]].getNodeSlowness()+
                               nodes[this->neighbors[cellNoi][3]].getNodeSlowness());

                    sxyz<T1> pt;
                    if(blti_raytrace(curr_pt, ind, pt, threadNo, s)){
                        T1 t=Interpolator<T1>::bilinearTime(pt, nodes[ind[0]], nodes[ind[1]], nodes[ind[2]], threadNo);
                        if (t>=time)
                            continue;
                        t+=curr_pt.getDistance(pt)*s;
                        if (t<t_i){
                            t_i=t;
                            Slow=s;
                            pt_i=pt;;
                            onNode = false;
                            onEdge = false;
                            onFace = true;
                            faceNodesPrev=faceNodes;
                            faceNodes = ind;
                            cellNo=cellNoi;
                        }
                    }
                    /////// minimum at each edges
                    std::array<T2,2> indEdges[3] = {{{ ind[0],ind[1]}},{{ind[0],ind[2]}},{{ ind[1],ind[2]}}};
                    for(T2 n=0;n<3;++n){

                        sxyz<T1> pt;
                        if(blti2D_raytrace(curr_pt, indEdges[n][0], indEdges[n][1], pt, threadNo, s) ==false) continue;
                        T1 dist0=pt.getDistance(nodes[indEdges[n][0]]);
                        T1 dist1=pt.getDistance(nodes[indEdges[n][1]]);
                        T1 t= (dist1*nodes[indEdges[n][0]].getTT(threadNo)+dist0*nodes[indEdges[n][1]].getTT(threadNo))/(dist0+dist1);
                        if (t>=time)
                            continue;
                        t+=curr_pt.getDistance(pt)*s;
                        if (t<t_i ){
                            t_i=t;
                            Slow=s;
                            pt_i=pt;
                            edgeNodes[0] = indEdges[n][0];
                            edgeNodes[1] = indEdges[n][1];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            cellNo=cellNoi;
                        }

                    }
                    ////// nodes
                    for (T2 i=0;i<4;++i){
                        if (this->neighbors[cellNoi][i]==nodeNo)
                            continue;
                        T1 t=nodes[this->neighbors[cellNoi][i]].getTT(threadNo);
                        if (t>=time)
                            continue;
                        t+=s*nodes[this->neighbors[cellNoi][i]].getDistance(curr_pt);
                        if (t<t_i){
                            t_i=t;
                            Slow=s;
                            pt_i=sxyz<T1>(nodes[this->neighbors[cellNoi][i]]);
                            onNode = true;
                            onEdge = false;
                            onFace = false;
                            newNode=this->neighbors[cellNoi][i];
                            cellNo=cellNoi;
                        }
                    }

                }

                prev_pt = curr_pt;
                curr_pt = pt_i;
                time=t_i-curr_pt.getDistance(prev_pt)*Slow;
                if (onNode){
                    nodeNo=newNode;
                }
                tt_s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                        faceNodes);

                tt += 0.5*(tt_s1 + tt_s2) * prev_pt.getDistance( curr_pt );
                tt_s1 = tt_s2;
                if (!NearSource && t_i==std::numeric_limits<T1>::max()){
                    std::cout << "\n\nWarning: finding raypath on node failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    reachedTx = true;

                }

            } else if ( onEdge ) {

                // find cells common to edge
                std::vector<T2> cells;
                T1 Slow;
                for ( auto nc0=nodes[edgeNodes[0]].getOwners().begin(); nc0!=nodes[edgeNodes[0]].getOwners().end(); ++nc0 ) {
                    if ( std::find(nodes[edgeNodes[1]].getOwners().begin(), nodes[edgeNodes[1]].getOwners().end(), *nc0)!=nodes[edgeNodes[1]].getOwners().end() ) {
                        cells.push_back( *nc0 );
                    }
                }
                std::array<T2,2> edgeNodestmp;
                T1 t_i=std::numeric_limits<T1>::max();
                sxyz<T1> pt_i;
                for (T2 nc=0; nc<cells.size(); ++nc ) {
                    T2 cellNoi = cells[nc];
                    // there are 2 faces that might be intersected
                    std::array<T2,2> edgeNodes2;
                    size_t n2=0;
                    for ( auto nn=this->neighbors[cellNoi].begin(); nn!= this->neighbors[cellNoi].begin()+4; ++nn ) {
                        if ( *nn!=edgeNodes[0] && *nn!=edgeNodes[1] ) {
                            edgeNodes2[n2++] = *nn;
                        }
                    }
                    std::array<T2,3> ind[2] = {
                        { { edgeNodes[0], edgeNodes2[0], edgeNodes2[1]}},
                        { { edgeNodes[1], edgeNodes2[0], edgeNodes2[1]}}};
                    std::sort( ind[0].begin(), ind[0].end());
                    std::sort( ind[1].begin(), ind[1].end());
                    std::array<T2,2> indEdges[5] = {
                        {{ edgeNodes[0],edgeNodes2[0]}},
                        {{ edgeNodes[0],edgeNodes2[1]}},
                        {{ edgeNodes[1],edgeNodes2[0]}},
                        {{ edgeNodes[1],edgeNodes2[1]}},
                        {{ edgeNodes2[0],edgeNodes2[1]}}};
                    T1 s=0.25*(nodes[edgeNodes[0]].getNodeSlowness()+
                               nodes[edgeNodes[1]].getNodeSlowness()+
                               nodes[edgeNodes2[0]].getNodeSlowness()+
                               nodes[edgeNodes2[1]].getNodeSlowness());
                    bool flag=false;
                    for ( size_t n=0; n<2; ++n ) {
                        if (NearSource){
                            std::array<T1,3> Barycenter;
                            if (blti_solver_around_source(NodeSource, curr_pt, ind[n], Barycenter)==true){
                                pt_i.x=Barycenter[0]*nodes[ind[n][0]].getX()+Barycenter[1]*nodes[ind[n][1]].getX()+Barycenter[2]*nodes[ind[n][2]].getX();
                                pt_i.y=Barycenter[0]*nodes[ind[n][0]].getY()+Barycenter[1]*nodes[ind[n][1]].getY()+Barycenter[2]*nodes[ind[n][2]].getY();
                                pt_i.z=Barycenter[0]*nodes[ind[n][0]].getZ()+Barycenter[1]*nodes[ind[n][1]].getZ()+Barycenter[2]*nodes[ind[n][2]].getZ();
                                if(NodeSource.getDistance(pt_i)>NodeSource.getDistance(curr_pt))
                                    continue;
                                for(T2 ni=0;ni<3;++ni){
                                    if(std::abs(1.0-Barycenter[ni])<minDist*minDist){
                                        onNode = true;
                                        onEdge = false;
                                        onFace = false;
                                        nodeNo=ind[n][ni];
                                        cellNo=cellNoi;
                                        flag=true;
                                        break;
                                    }
                                }
                                if (flag)
                                    break;
                                for(T2 ni=0;ni<3;++ni){
                                    if(std::abs(Barycenter[ni])<minDist*minDist){
                                        onNode = false;
                                        onEdge = true;
                                        onFace = false;
                                        edgeNodes[0] = ind[n][(ni+1)%3];
                                        edgeNodes[1] = ind[n][(ni+2)%3];
                                        cellNo=cellNoi;
                                        flag=true;
                                        break;
                                    }
                                }
                                if (flag)
                                    break;
                                onNode = false;
                                onEdge = false;
                                onFace = true;
                                faceNodesPrev=faceNodes;
                                faceNodes=ind[n];
                                cellNo=cellNoi;
                                flag=true;
                                break;
                            }
                            continue;
                        }
                        if ( ind[n] == faceNodesPrev) continue;
                        sxyz<T1> pt;
                        if(blti_raytrace(curr_pt, ind[n], pt, threadNo, s)==false) continue;
                        T1 t=Interpolator<T1>::bilinearTime(pt, nodes[ind[n][0]], nodes[ind[n][1]], nodes[ind[n][2]], threadNo);
                        if (t>=time)
                            continue;
                        t+=curr_pt.getDistance(pt)*s;

                        if (t<t_i){
                            t_i=t;
                            Slow=s;
                            pt_i=pt;
                            onNode = false;
                            onEdge = false;
                            onFace = true;
                            faceNodesPrev=faceNodes;
                            faceNodes = ind[n];
                            cellNo=cellNoi;
                        }
                    }
                    if (flag)
                        break;
                    if (NearSource)
                        continue;
                    for(T2 n=0;n<5;++n){
                        sxyz<T1> pt;
                        if(blti2D_raytrace(curr_pt, indEdges[n][0], indEdges[n][1], pt, threadNo, s) ==false) continue;
                        T1 dist0=pt.getDistance(nodes[indEdges[n][0]]);
                        T1 dist1=pt.getDistance(nodes[indEdges[n][1]]);
                        T1 t= (dist1*nodes[indEdges[n][0]].getTT(threadNo)+dist0*nodes[indEdges[n][1]].getTT(threadNo))/(dist0+dist1);
                        if (t>=time)
                            continue;
                        t+=curr_pt.getDistance(pt)*s;
                        if (t<t_i ){
                            t_i=t;
                            Slow=s;
                            pt_i=pt;
                            edgeNodestmp[0] = indEdges[n][0];
                            edgeNodestmp[1] = indEdges[n][1];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            cellNo=cellNoi;

                        }

                    }
                    for (T2 i=0;i<4;++i){
                        T1 t=nodes[this->neighbors[cellNoi][i]].getTT(threadNo);
                        if (t>=time)
                            continue;
                        t+=s*nodes[this->neighbors[cellNoi][i]].getDistance(curr_pt);
                        if (t<t_i){
                            t_i=t;
                            Slow=s;
                            pt_i=sxyz<T1>(nodes[this->neighbors[cellNoi][i]]);
                            onNode = true;
                            onEdge = false;
                            onFace = false;
                            nodeNo=this->neighbors[cellNoi][i];
                            cellNo=cellNoi;
                        }
                    }

                    // find next cell
                }
                if (!NearSource &&  t_i==std::numeric_limits<T1>::max()){
                    std::cout << "\n\nWarning: finding raypath on edge failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    reachedTx = true;

                }
                if (onEdge){
                    edgeNodes[0] = edgeNodestmp[0];
                    edgeNodes[1] = edgeNodestmp[1];
                }
                prev_pt = curr_pt;
                curr_pt = pt_i;
                tt_s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                        faceNodes);

                tt += 0.5*(tt_s1 + tt_s2) * prev_pt.getDistance( curr_pt );
                tt_s1 = tt_s2;
                time=t_i-curr_pt.getDistance(prev_pt)*Slow;
            } else{ // on Face


                T2 cellNo1;
                if (!onFace || inlimitD)
                    cellNo1=cellNo;
                else
                    cellNo1=findAdjacentCell2(faceNodes, cellNo,curr_pt);

                sxyz<T1> pt_i;
                T1 Slow;
                std::array<T2,3> ind[8] = {
                    { { this->neighbors[cellNo][0], this->neighbors[cellNo][1], this->neighbors[cellNo][2] } },
                    { { this->neighbors[cellNo][0], this->neighbors[cellNo][1], this->neighbors[cellNo][3] } },
                    { { this->neighbors[cellNo][0], this->neighbors[cellNo][2], this->neighbors[cellNo][3] } },
                    { { this->neighbors[cellNo][1], this->neighbors[cellNo][2], this->neighbors[cellNo][3] } },
                    { { this->neighbors[cellNo1][0], this->neighbors[cellNo1][1], this->neighbors[cellNo1][2] } },
                    { { this->neighbors[cellNo1][0], this->neighbors[cellNo1][1], this->neighbors[cellNo1][3] } },
                    { { this->neighbors[cellNo1][0], this->neighbors[cellNo1][2], this->neighbors[cellNo1][3] } },
                    { { this->neighbors[cellNo1][1], this->neighbors[cellNo1][2], this->neighbors[cellNo1][3] } }
                };
                for ( size_t n=0; n<8; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }
                // there are 3 faces that might be intersected

                T1 t_i=std::numeric_limits<T1>::max();
                T2 face;
                T1 s1=0.25*(nodes[this->neighbors[cellNo][0]].getNodeSlowness()+
                            nodes[this->neighbors[cellNo][1]].getNodeSlowness()+
                            nodes[this->neighbors[cellNo][2]].getNodeSlowness()+
                            nodes[this->neighbors[cellNo][3]].getNodeSlowness());
                T1 s2=0.25*(nodes[this->neighbors[cellNo1][0]].getNodeSlowness()+
                            nodes[this->neighbors[cellNo1][1]].getNodeSlowness()+
                            nodes[this->neighbors[cellNo1][2]].getNodeSlowness()+
                            nodes[this->neighbors[cellNo1][3]].getNodeSlowness());

                for ( size_t n=0; n<8; ++n ) {
                    if ( ind[n] == faceNodes ||   ind[n] == faceNodesPrev ) continue;
                    if (NearSource){
                        bool flag=false;
                        std::array<T1,3> Barycenter;
                        if (blti_solver_around_source(NodeSource, curr_pt, ind[n], Barycenter)==true){
                            pt_i.x=Barycenter[0]*nodes[ind[n][0]].getX()+Barycenter[1]*nodes[ind[n][1]].getX()+Barycenter[2]*nodes[ind[n][2]].getX();
                            pt_i.y=Barycenter[0]*nodes[ind[n][0]].getY()+Barycenter[1]*nodes[ind[n][1]].getY()+Barycenter[2]*nodes[ind[n][2]].getY();
                            pt_i.z=Barycenter[0]*nodes[ind[n][0]].getZ()+Barycenter[1]*nodes[ind[n][1]].getZ()+Barycenter[2]*nodes[ind[n][2]].getZ();

                            if(NodeSource.getDistance(pt_i)>NodeSource.getDistance(curr_pt))
                                continue;
                            for(T2 ni=0;ni<3;++ni){
                                if(std::abs(1.0-Barycenter[ni])<minDist*minDist){
                                    onNode = true;
                                    onEdge = false;
                                    onFace = false;
                                    nodeNo=ind[n][ni];
                                    faceNodesPrev=faceNodes;
                                    flag=true;
                                    break;
                                }
                            }
                            if (flag)
                                break;
                            for(T2 ni=0;ni<3;++ni){
                                if(std::abs(Barycenter[ni])<minDist*minDist){
                                    onNode = false;
                                    onEdge = true;
                                    onFace = false;
                                    edgeNodes[0] = ind[n][(ni+1)%3];
                                    edgeNodes[1] = ind[n][(ni+2)%3];
                                    faceNodesPrev=faceNodes;
                                    break;
                                }
                            }
                            if (flag)
                                break;
                            onNode = false;
                            onEdge = false;
                            onFace = true;
                            face=static_cast<T1>(n);
                            break;
                        }
                        continue;
                    }
                    T1 s;
                    sxyz<T1> pt;
                    s=(n<4)?s1:s2;

                    if(blti_raytrace(curr_pt, ind[n], pt,threadNo, s)==false) continue;
                    T1 t=Interpolator<T1>::bilinearTime(pt, nodes[ind[n][0]], nodes[ind[n][1]], nodes[ind[n][2]], threadNo);
                    if (t>=time)
                        continue;
                    t+=curr_pt.getDistance(pt)*s;
                    if (t<t_i){
                        t_i=t;
                        Slow=s;
                        pt_i=pt;
                        onNode = false;
                        onEdge = false;
                        onFace = true;
                        face=static_cast<T1>(n);
                    }
                }
                if (! NearSource){
                    std::array<T2,2> indEdges[9] = {
                        {{ this->neighbors[cellNo][0],this->neighbors[cellNo][1]}},
                        {{ this->neighbors[cellNo][0],this->neighbors[cellNo][2]}},
                        {{ this->neighbors[cellNo][0],this->neighbors[cellNo][3]}},
                        {{ this->neighbors[cellNo][1],this->neighbors[cellNo][2]}},
                        {{ this->neighbors[cellNo][1],this->neighbors[cellNo][3]}},
                        {{ this->neighbors[cellNo][2],this->neighbors[cellNo][3]}}
                    };
                    T2 forthNode=0;
                    for (T2 i=0; i<4;++i){
                        if (this->neighbors[cellNo1][i]!=this->neighbors[cellNo][0] &&
                            this->neighbors[cellNo1][i]!=this->neighbors[cellNo][1]&&
                            this->neighbors[cellNo1][i]!=this->neighbors[cellNo][2]&&
                            this->neighbors[cellNo1][i]!=this->neighbors[cellNo][3]){
                            forthNode=i;
                            break;
                        }
                    }
                    indEdges[6]={this->neighbors[cellNo1][forthNode],this->neighbors[cellNo1][(forthNode+1)%4]};
                    indEdges[7]={this->neighbors[cellNo1][forthNode],this->neighbors[cellNo1][(forthNode+2)%4]};
                    indEdges[8]={this->neighbors[cellNo1][forthNode],this->neighbors[cellNo1][(forthNode+3)%4]};
                    for(T2 n=0;n<9;++n){
                        T1 s;
                        sxyz<T1> pt;
                        s=(n<6)?s1:s2;
                        if(blti2D_raytrace(curr_pt, indEdges[n][0], indEdges[n][1], pt, threadNo, s) ==false) continue;
                        T1 dist0=pt.getDistance(nodes[indEdges[n][0]]);
                        T1 dist1=pt.getDistance(nodes[indEdges[n][1]]);
                        T1 t= (dist1*nodes[indEdges[n][0]].getTT(threadNo)+dist0*nodes[indEdges[n][1]].getTT(threadNo))/(dist0+dist1);
                        if (t>=time)
                            continue;
                        t+=curr_pt.getDistance(pt)*s;
                        if (t<t_i ){
                            t_i=t;
                            Slow=s;
                            pt_i=pt;
                            edgeNodes[0] = indEdges[n][0];
                            edgeNodes[1] = indEdges[n][1];
                            faceNodesPrev=faceNodes;
                            onNode = false;
                            onEdge = true;
                            onFace = false;;

                        }

                    }
                    for (T2 n=0;n<4;++n){
                        T1 t=nodes[this->neighbors[cellNo][n]].getTT(threadNo);
                        if(t<time){
                            sxyz<T1> pt;
                            t+=s1*nodes[this->neighbors[cellNo][n]].getDistance(curr_pt);
                            if(t<t_i){
                                t_i=t;
                                Slow=s1;
                                pt_i=sxyz<T1>(nodes[this->neighbors[cellNo][n]]);
                                onNode = true;
                                onEdge = false;
                                onFace = false;
                                faceNodesPrev=faceNodes;
                                nodeNo=this->neighbors[cellNo][n];
                            }
                        }
                    }
                    if(nodes[this->neighbors[cellNo1][forthNode]].getTT(threadNo)<time){
                        sxyz<T1> pt;
                        if(nodes[this->neighbors[cellNo1][forthNode]].getTT(threadNo)+s2*nodes[this->neighbors[cellNo1][forthNode]].getDistance(curr_pt)<t_i){
                            t_i=nodes[this->neighbors[cellNo1][forthNode]].getTT(threadNo)+s2*nodes[this->neighbors[cellNo1][forthNode]].getDistance(curr_pt);
                            Slow=s2;
                            pt_i=sxyz<T1>(nodes[this->neighbors[cellNo1][forthNode]]);
                            onNode = true;
                            onEdge = false;
                            onFace = false;
                            faceNodesPrev=faceNodes;
                            nodeNo=this->neighbors[cellNo1][forthNode];
                        }
                    }
                }
                if (t_i==std::numeric_limits<T1>::max() && ! NearSource){
                    std::cout << "\n\nWarning: finding raypath on face failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    reachedTx = true;
                    continue;
                }

                prev_pt = curr_pt;
                curr_pt = pt_i;
                inlimitD=false;
                if (std::abs(curr_pt.x-Xmax)<=minDist*minDist ||std::abs(curr_pt.x-Xmin)<=minDist*minDist )
                    inlimitD=true;
                if (std::abs(curr_pt.y-Ymax)<=minDist*minDist ||std::abs(curr_pt.y-Ymin)<=minDist*minDist )
                    inlimitD=true;
                if (std::abs(curr_pt.z-Zmax)<=minDist*minDist ||std::abs(curr_pt.z-Zmin)<=minDist*minDist )
                    inlimitD=true;

                if (onFace){
                    faceNodesPrev=faceNodes;
                    faceNodes = ind[face];
                    if (!inlimitD){
                        if (face<4)
                            cellNo=findAdjacentCell2(faceNodes, cellNo,curr_pt);

                        else
                            cellNo=findAdjacentCell2(faceNodes, cellNo1,curr_pt);
                    }
                }

                tt_s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                        faceNodes);
                tt += 0.5*(tt_s1 + tt_s2) * prev_pt.getDistance( curr_pt );
                tt_s1 = tt_s2;
                time=t_i-curr_pt.getDistance(prev_pt)*Slow;

            }
            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( curr_pt.getDistance( Tx[nt] ) < minDist*minDist ) {
                        std::array<T2,4> itmp = getPrimary(txCell[nt]);
                        if ( processVel )

                            tt_s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                           nodes[itmp[0]],
                                                                           nodes[itmp[1]],
                                                                           nodes[itmp[2]],
                                                                           nodes[itmp[3]]);
                        else
                            tt_s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                        nodes[itmp[0]],
                                                                        nodes[itmp[1]],
                                                                        nodes[itmp[2]],
                                                                        nodes[itmp[3]]);

                        tt += t0[nt] + 0.5*(tt_s1 + tt_s2) * curr_pt.getDistance( Tx[nt] );
                        reachedTx = true;
                        break;
                    }
                    for(size_t n=0;n<4;++n){
                        sxyz<T1> NearTx={nodes[this->neighbors[txCell[nt]][n]]};
                        if (curr_pt.getDistance(NearTx)< minDist*minDist){
                            std::array<T2,4> itmp = getPrimary(txCell[nt]);
                            if ( processVel )

                                tt_s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                               nodes[itmp[0]],
                                                                               nodes[itmp[1]],
                                                                               nodes[itmp[2]],
                                                                               nodes[itmp[3]]);
                            else
                                tt_s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                            nodes[itmp[0]],
                                                                            nodes[itmp[1]],
                                                                            nodes[itmp[2]],
                                                                            nodes[itmp[3]]);

                            tt += t0[nt] + 0.5*(tt_s1 + tt_s2) * curr_pt.getDistance( Tx[nt] );
                            reachedTx = true;
                            break;
                        }
                    }

                }
            } else {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        for ( auto nc=nodes[txNode[nt]].getOwners().begin(); nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                            if ( cellNo == *nc ) {
                                std::array<T2,4> itmp = getPrimary(cellNo);

                                if ( processVel )
                                    tt_s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                   nodes[itmp[0]],
                                                                                   nodes[itmp[1]],
                                                                                   nodes[itmp[2]],
                                                                                   nodes[itmp[3]]);
                                else
                                    tt_s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                                nodes[itmp[0]],
                                                                                nodes[itmp[1]],
                                                                                nodes[itmp[2]],
                                                                                nodes[itmp[3]]);

                                tt += t0[nt] + 0.5*(tt_s1 + tt_s2) * curr_pt.getDistance( Tx[nt] );
                                reachedTx = true;
                                break;
                            }
                        }
                    } else {
                        if ( cellNo == txCell[nt] ) {
                            std::array<T2,4> itmp = getPrimary(cellNo);
                            if ( processVel )
                                tt_s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                               nodes[itmp[0]],
                                                                               nodes[itmp[1]],
                                                                               nodes[itmp[2]],
                                                                               nodes[itmp[3]]);
                            else
                                tt_s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                            nodes[itmp[0]],
                                                                            nodes[itmp[1]],
                                                                            nodes[itmp[2]],
                                                                            nodes[itmp[3]]);

                            tt += t0[nt] + 0.5*(tt_s1 + tt_s2) * curr_pt.getDistance( Tx[nt] );
                            reachedTx = true;
                        } else {
                            for ( size_t nn=0; nn<txNeighborCells[nt].size(); ++nn ) {
                                if ( cellNo == txNeighborCells[nt][nn] ) {
                                    std::array<T2,4> itmp = getPrimary(txCell[nt]);
                                    if ( processVel )
                                        tt_s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                       nodes[itmp[0]],
                                                                                       nodes[itmp[1]],
                                                                                       nodes[itmp[2]],
                                                                                       nodes[itmp[3]]);
                                    else
                                        tt_s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                                    nodes[itmp[0]],
                                                                                    nodes[itmp[1]],
                                                                                    nodes[itmp[2]],
                                                                                    nodes[itmp[3]]);

                                    tt += t0[nt] + 0.5*(tt_s1 + tt_s2) * curr_pt.getDistance( Tx[nt] );
                                    reachedTx = true;
                                    break;
                                }
                            }
                        }
                    }
                    if ( reachedTx ) break;
                }
            }
        }
        return tt;
    }


    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::getRaypath_blti(const std::vector<sxyz<T1>>& Tx,
                                               const std::vector<T1>& t0,
                                               const sxyz<T1> &Rx,
                                               std::vector<sxyz<T1>> &r_data,
                                               T1& tt,
                                               const size_t threadNo) const {

        T1 minDist = min_dist;
        std::vector<sxyz<T1>> r_tmp;
        r_tmp.emplace_back( Rx );
        tt = 0.0;
        T1 tt_s1, tt_s2;

        for ( size_t ns=0; ns<Tx.size(); ++ns ) {
            if ( Rx == Tx[ns] ) {
                tt = t0[ns];
                return;
            }
        }

        const txInfo_t& txi = getTxInfoBlti(Tx, threadNo);
        const std::vector<bool>& txOnNode = txi.txOnNode;
        const std::vector<T2>& txNode = txi.txNode;
        const std::vector<T2>& txCell = txi.txCell;
        const std::vector<std::vector<T2>>& txNeighborCells = txi.txNeighborCells;

        T2 cellNo, nodeNo;
        sxyz<T1> curr_pt( Rx ), prev_pt( Rx );
        bool onNode = false;
        bool onEdge = false;
        bool onFace = false;
        std::array<T2,2> edgeNodes;
        std::array<T2,3> faceNodes={{0,0,0}};
        std::array<T2,3> faceNodesPrev={{0,0,0}};
        bool reachedTx = false;
        T1 Xmax=-std::numeric_limits<T1>::max();
        T1 Xmin=std::numeric_limits<T1>::max();
        T1 Ymax=-std::numeric_limits<T1>::max();
        T1 Ymin=std::numeric_limits<T1>::max();
        T1 Zmax=-std::numeric_limits<T1>::max();
        T1 Zmin=std::numeric_limits<T1>::max();
        for ( T2 nn=0; nn<nodes.size(); ++nn ) {
            if ( nodes[nn] == curr_pt ) {
                nodeNo = nn;
                onNode = true;
                tt_s1 = nodes[nodeNo].getNodeSlowness();
            }
            Xmax=nodes[nn].getX()>Xmax?nodes[nn].getX():Xmax;
            Ymax=nodes[nn].getY()>Ymax?nodes[nn].getY():Ymax;
            Zmax=nodes[nn].getZ()>Zmax?nodes[nn].getZ():Zmax;

            Xmin=nodes[nn].getX()<Xmin?nodes[nn].getX():Xmin;
            Ymin=nodes[nn].getY()<Ymin?nodes[nn].getY():Ymin;
            Zmin=nodes[nn].getZ()<Zmin?nodes[nn].getZ():Zmin;
        }
        if ( !onNode ) {
            cellNo = getCellNo( curr_pt );

            std::array<T2,4> itmp = getPrimary(cellNo);
            T2 ind[6][2] = {
                {itmp[0], itmp[1]},
                {itmp[0], itmp[2]},
                {itmp[0], itmp[3]},
                {itmp[1], itmp[2]},
                {itmp[1], itmp[3]},
                {itmp[2], itmp[3]} };

            for ( size_t n=0; n<6; ++n ) {
                if ( areCollinear(curr_pt, nodes[ind[n][0]], nodes[ind[n][1]]) ) {
                    onEdge = true;
                    edgeNodes[0] = ind[n][0];
                    edgeNodes[1] = ind[n][1];

                    if ( processVel )
                        tt_s1 = Interpolator<T1>::linearVel(r_tmp.back(),
                                                            nodes[edgeNodes[0]],
                                                            nodes[edgeNodes[1]]);
                    else
                        tt_s1 = Interpolator<T1>::linear(r_tmp.back(),
                                                         nodes[edgeNodes[0]],
                                                         nodes[edgeNodes[1]]);

                    break;
                }
            }
            if ( !onEdge ) {
                std::array<T2,3> ind[4] = {
                    { { itmp[0], itmp[1], itmp[2] } },
                    { { itmp[0], itmp[1], itmp[3] } },
                    { { itmp[0], itmp[2], itmp[3] } },
                    { { itmp[1], itmp[2], itmp[3] } }
                };

                for ( size_t n=0; n<4; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }

                for ( size_t n=0; n<4; ++n ) {
                    if ( areCoplanar(curr_pt, nodes[ind[n][0]], nodes[ind[n][1]], nodes[ind[n][2]]) ) {
                        onFace = true;
                        faceNodes[0] = ind[n][0];
                        faceNodes[1] = ind[n][1];
                        faceNodes[2] = ind[n][2];

                        if ( processVel )
                            tt_s1 = Interpolator<T1>::bilinearTriangleVel(curr_pt,
                                                                          nodes[faceNodes[0]],
                                                                          nodes[faceNodes[1]],
                                                                          nodes[faceNodes[2]]);
                        else
                            tt_s1 = Interpolator<T1>::bilinearTriangle(curr_pt,
                                                                       nodes[faceNodes[0]],
                                                                       nodes[faceNodes[1]],
                                                                       nodes[faceNodes[2]]);
                        break;
                    }
                }
            }
        }

        if (!onEdge && ! onNode && ! onFace){
            std::array<T2,4> itmp = getPrimary(cellNo);
            if ( processVel )
                tt_s1 = Interpolator<T1>::trilinearTriangleVel(curr_pt,
                                                               nodes[itmp[0]],
                                                               nodes[itmp[1]],
                                                               nodes[itmp[2]],
                                                               nodes[itmp[3]]);
            else
                tt_s1 = Interpolator<T1>::trilinearTriangle(curr_pt,
                                                            nodes[itmp[0]],
                                                            nodes[itmp[1]],
                                                            nodes[itmp[2]],
                                                            nodes[itmp[3]]);
        }
        for(auto nt=0;nt<txCell.size();++nt){
            if (getCellNo( Rx )==txCell[nt]){
                std::array<T2,4> itmp = getPrimary(cellNo);
                if ( processVel )
                    tt_s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                   nodes[itmp[0]],
                                                                   nodes[itmp[1]],
                                                                   nodes[itmp[2]],
                                                                   nodes[itmp[3]]);
                else
                    tt_s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                nodes[itmp[0]],
                                                                nodes[itmp[1]],
                                                                nodes[itmp[2]],
                                                                nodes[itmp[3]]);

                tt += t0[nt] + 0.5*(tt_s1 + tt_s2) * r_tmp.back().getDistance( Tx[nt] );
                r_tmp.push_back(Tx[nt]);
                reachedTx=true;
                break;
            }
        }

        T1 time=std::numeric_limits<T1>::max();
        bool inlimitD=false;
        while ( reachedTx == false ) {
            sxyz<T1> NodeSource;
            bool NearSource=false;
            for(size_t nt=0;nt<Tx.size();++nt){
                for(auto n=this->neighbors[txCell[nt]].begin();n!=this->neighbors[txCell[nt]].begin()+4;++n){
                    for(auto nc=nodes[*n].getOwners().begin();nc!=nodes[*n].getOwners().end();++nc){
                        if(*nc==cellNo){
                            NearSource=true;
                            NodeSource=Tx[nt];
                            break;
                        }
                    }
                    if (NearSource)
                        break;
                }
                if (curr_pt.getDistance(Tx[nt])<50.0*1.e-3){
                    NearSource=true;
                    NodeSource=Tx[nt];
                    break;
                }

            }

            if ( onNode ) {

                T1 t_i=std::numeric_limits<T1>::max();
                T1 Slow;
                sxyz<T1> pt_i;
                T2 newNode;
                for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                    T2 cellNoi = *nc;
                    std::array<T2, 3>ind;
                    for (T2 i=0;i<4;++i){
                        if (this->neighbors[cellNoi][i]==nodeNo){
                            ind[0]=this->neighbors[cellNoi][(i+1)%4];
                            ind[1]=this->neighbors[cellNoi][(i+2)%4];
                            ind[2]=this->neighbors[cellNoi][(i+3)%4];
                            break;
                        }
                    }
                    std::sort(ind.begin(), ind.end());
                    if (NearSource){
                        bool flag=false;
                        std::array<T1,3> Barycenter;
                        if (blti_solver_around_source(NodeSource, curr_pt, ind, Barycenter)==true){
                            pt_i.x=Barycenter[0]*nodes[ind[0]].getX()+Barycenter[1]*nodes[ind[1]].getX()+Barycenter[2]*nodes[ind[2]].getX();
                            pt_i.y=Barycenter[0]*nodes[ind[0]].getY()+Barycenter[1]*nodes[ind[1]].getY()+Barycenter[2]*nodes[ind[2]].getY();
                            pt_i.z=Barycenter[0]*nodes[ind[0]].getZ()+Barycenter[1]*nodes[ind[1]].getZ()+Barycenter[2]*nodes[ind[2]].getZ();
                            if(NodeSource.getDistance(pt_i)>NodeSource.getDistance(curr_pt))
                                continue;
                            for(T2 ni=0;ni<3;++ni){
                                if(std::abs(1.0-Barycenter[ni])<minDist*minDist){
                                    onNode = true;
                                    onEdge = false;
                                    onFace = false;
                                    newNode=ind[ni];
                                    cellNo=cellNoi;
                                    flag=true;
                                    break;
                                }
                            }
                            if (flag)
                                break;
                            for(T2 ni=0;ni<3;++ni){
                                if(std::abs(Barycenter[ni])<minDist*minDist){
                                    onNode = false;
                                    onEdge = true;
                                    onFace = false;
                                    edgeNodes[0] = ind[(ni+1)%3];
                                    edgeNodes[1] = ind[(ni+2)%3];
                                    cellNo=cellNoi;
                                    flag=true;
                                    break;
                                }
                            }
                            if (flag)
                                break;

                            onNode = false;
                            onEdge = false;
                            onFace = true;
                            faceNodesPrev=faceNodes;
                            faceNodes=ind;
                            cellNo=cellNoi;
                            break;
                        }
                        continue;
                    }
                    if ( ind == faceNodesPrev) continue;
                    T1 s=0.25*(nodes[this->neighbors[cellNoi][0]].getNodeSlowness()+
                               nodes[this->neighbors[cellNoi][1]].getNodeSlowness()+
                               nodes[this->neighbors[cellNoi][2]].getNodeSlowness()+
                               nodes[this->neighbors[cellNoi][3]].getNodeSlowness());

                    sxyz<T1> pt;
                    if(blti_raytrace(curr_pt, ind, pt, threadNo, s)){
                        T1 t=Interpolator<T1>::bilinearTime(pt, nodes[ind[0]], nodes[ind[1]], nodes[ind[2]], threadNo);
                        if (t>=time)
                            continue;
                        t+=curr_pt.getDistance(pt)*s;
                        if (t<t_i){
                            t_i=t;
                            Slow=s;
                            pt_i=pt;;
                            onNode = false;
                            onEdge = false;
                            onFace = true;
                            faceNodesPrev=faceNodes;
                            faceNodes = ind;
                            cellNo=cellNoi;
                        }
                    }
                    /////// minimum at each edges
                    std::array<T2,2> indEdges[3] = {{{ ind[0],ind[1]}},{{ind[0],ind[2]}},{{ ind[1],ind[2]}}};
                    for(T2 n=0;n<3;++n){

                        sxyz<T1> pt;
                        if(blti2D_raytrace(curr_pt, indEdges[n][0], indEdges[n][1], pt, threadNo, s) ==false) continue;
                        T1 dist0=pt.getDistance(nodes[indEdges[n][0]]);
                        T1 dist1=pt.getDistance(nodes[indEdges[n][1]]);
                        T1 t= (dist1*nodes[indEdges[n][0]].getTT(threadNo)+dist0*nodes[indEdges[n][1]].getTT(threadNo))/(dist0+dist1);
                        if (t>time)
                            continue;
                        t+=curr_pt.getDistance(pt)*s;
                        if (t<t_i ){
                            t_i=t;
                            Slow=s;
                            pt_i=pt;
                            edgeNodes[0] = indEdges[n][0];
                            edgeNodes[1] = indEdges[n][1];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            cellNo=cellNoi;
                        }

                    }
                    ////// nodes
                    for (T2 i=0;i<4;++i){
                        if (this->neighbors[cellNoi][i]==nodeNo)
                            continue;
                        T1 t=nodes[this->neighbors[cellNoi][i]].getTT(threadNo);
                        if (t>=time)
                            continue;
                        t+=s*nodes[this->neighbors[cellNoi][i]].getDistance(curr_pt);
                        if (t<t_i){
                            t_i=t;
                            Slow=s;
                            pt_i=sxyz<T1>(nodes[this->neighbors[cellNoi][i]]);
                            onNode = true;
                            onEdge = false;
                            onFace = false;
                            newNode=this->neighbors[cellNoi][i];
                            cellNo=cellNoi;
                        }
                    }

                }

                prev_pt = curr_pt;
                curr_pt = pt_i;
                time=t_i-curr_pt.getDistance(prev_pt)*Slow;
                if (onNode){
                    nodeNo=newNode;
                }
                tt_s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                        faceNodes);
                tt += 0.5*(tt_s1 + tt_s2) * r_tmp.back().getDistance( curr_pt );
                tt_s1 = tt_s2;
                r_tmp.push_back( curr_pt );
                if (!NearSource && t_i==std::numeric_limits<T1>::max()){
                    std::cout << "\n\nWarning: finding raypath on edge failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    tt = 0.0;
                    reachedTx = true;
                }

            } else if ( onEdge ) {

                // find cells common to edge
                std::vector<T2> cells;
                T1 Slow;
                for ( auto nc0=nodes[edgeNodes[0]].getOwners().begin(); nc0!=nodes[edgeNodes[0]].getOwners().end(); ++nc0 ) {
                    if ( std::find(nodes[edgeNodes[1]].getOwners().begin(), nodes[edgeNodes[1]].getOwners().end(), *nc0)!=nodes[edgeNodes[1]].getOwners().end() ) {
                        cells.push_back( *nc0 );
                    }
                }
                std::array<T2,2> edgeNodestmp;
                T1 t_i=std::numeric_limits<T1>::max();
                sxyz<T1> pt_i;
                for (T2 nc=0; nc<cells.size(); ++nc ) {
                    T2 cellNoi = cells[nc];
                    // there are 2 faces that might be intersected
                    std::array<T2,2> edgeNodes2;
                    size_t n2=0;
                    for ( auto nn=this->neighbors[cellNoi].begin(); nn!= this->neighbors[cellNoi].begin()+4; ++nn ) {
                        if ( *nn!=edgeNodes[0] && *nn!=edgeNodes[1] ) {
                            edgeNodes2[n2++] = *nn;
                        }
                    }
                    std::array<T2,3> ind[2] = {
                        { { edgeNodes[0], edgeNodes2[0], edgeNodes2[1]}},
                        { { edgeNodes[1], edgeNodes2[0], edgeNodes2[1]}}};
                    std::sort( ind[0].begin(), ind[0].end());
                    std::sort( ind[1].begin(), ind[1].end());
                    std::array<T2,2> indEdges[5] = {
                        {{ edgeNodes[0],edgeNodes2[0]}},
                        {{ edgeNodes[0],edgeNodes2[1]}},
                        {{ edgeNodes[1],edgeNodes2[0]}},
                        {{ edgeNodes[1],edgeNodes2[1]}},
                        {{ edgeNodes2[0],edgeNodes2[1]}}};
                    T1 s=0.25*(nodes[edgeNodes[0]].getNodeSlowness()+
                               nodes[edgeNodes[1]].getNodeSlowness()+
                               nodes[edgeNodes2[0]].getNodeSlowness()+
                               nodes[edgeNodes2[1]].getNodeSlowness());
                    bool flag=false;
                    for ( size_t n=0; n<2; ++n ) {
                        if (NearSource){
                            std::array<T1,3> Barycenter;
                            if (blti_solver_around_source(NodeSource, curr_pt, ind[n], Barycenter)==true){
                                pt_i.x=Barycenter[0]*nodes[ind[n][0]].getX()+Barycenter[1]*nodes[ind[n][1]].getX()+Barycenter[2]*nodes[ind[n][2]].getX();
                                pt_i.y=Barycenter[0]*nodes[ind[n][0]].getY()+Barycenter[1]*nodes[ind[n][1]].getY()+Barycenter[2]*nodes[ind[n][2]].getY();
                                pt_i.z=Barycenter[0]*nodes[ind[n][0]].getZ()+Barycenter[1]*nodes[ind[n][1]].getZ()+Barycenter[2]*nodes[ind[n][2]].getZ();
                                if(NodeSource.getDistance(pt_i)>NodeSource.getDistance(curr_pt))
                                    continue;
                                for(T2 ni=0;ni<3;++ni){
                                    if(std::abs(1.0-Barycenter[ni])<minDist*minDist){
                                        onNode = true;
                                        onEdge = false;
                                        onFace = false;
                                        nodeNo=ind[n][ni];
                                        cellNo=cellNoi;
                                        flag=true;
                                        break;
                                    }
                                }
                                if (flag)
                                    break;
                                for(T2 ni=0;ni<3;++ni){
                                    if(std::abs(Barycenter[ni])<minDist*minDist){
                                        onNode = false;
                                        onEdge = true;
                                        onFace = false;
                                        edgeNodes[0] = ind[n][(ni+1)%3];
                                        edgeNodes[1] = ind[n][(ni+2)%3];
                                        cellNo=cellNoi;
                                        flag=true;
                                        break;
                                    }
                                }
                                if (flag)
                                    break;
                                onNode = false;
                                onEdge = false;
                                onFace = true;
                                faceNodesPrev=faceNodes;
                                faceNodes=ind[n];
                                cellNo=cellNoi;
                                flag=true;
                                break;
                            }
                            continue;
                        }
                        if ( ind[n] == faceNodesPrev) continue;
                        sxyz<T1> pt;
                        if(blti_raytrace(curr_pt, ind[n], pt, threadNo, s)==false) continue;
                        T1 t=Interpolator<T1>::bilinearTime(pt, nodes[ind[n][0]], nodes[ind[n][1]], nodes[ind[n][2]], threadNo);
                        if (t>time)
                            continue;
                        t+=curr_pt.getDistance(pt)*s;

                        if (t<t_i){
                            t_i=t;
                            Slow=s;
                            pt_i=pt;
                            onNode = false;
                            onEdge = false;
                            onFace = true;
                            faceNodesPrev=faceNodes;
                            faceNodes = ind[n];
                            cellNo=cellNoi;
                        }
                    }
                    if (flag)
                        break;
                    if (NearSource)
                        continue;
                    for(T2 n=0;n<5;++n){
                        sxyz<T1> pt;
                        if(blti2D_raytrace(curr_pt, indEdges[n][0], indEdges[n][1], pt, threadNo, s) ==false) continue;
                        T1 dist0=pt.getDistance(nodes[indEdges[n][0]]);
                        T1 dist1=pt.getDistance(nodes[indEdges[n][1]]);
                        T1 t= (dist1*nodes[indEdges[n][0]].getTT(threadNo)+dist0*nodes[indEdges[n][1]].getTT(threadNo))/(dist0+dist1);
                        if (t>time)
                            continue;
                        t+=curr_pt.getDistance(pt)*s;
                        if (t<t_i ){
                            t_i=t;
                            Slow=s;
                            pt_i=pt;
                            edgeNodestmp[0] = indEdges[n][0];
                            edgeNodestmp[1] = indEdges[n][1];
                            onNode = false;
                            onEdge = true;
                            onFace = false;
                            cellNo=cellNoi;

                        }

                    }
                    for (T2 i=0;i<4;++i){
                        T1 t=nodes[this->neighbors[cellNoi][i]].getTT(threadNo);
                        if (t>time)
                            continue;
                        t+=s*nodes[this->neighbors[cellNoi][i]].getDistance(curr_pt);
                        if (t<t_i){
                            t_i=t;
                            Slow=s;
                            pt_i=sxyz<T1>(nodes[this->neighbors[cellNoi][i]]);
                            onNode = true;
                            onEdge = false;
                            onFace = false;
                            nodeNo=this->neighbors[cellNoi][i];
                            cellNo=cellNoi;
                        }
                    }

                    // find next cell
                }
                if (!NearSource && t_i==std::numeric_limits<T1>::max()) {
                    std::cout << "\n\nWarning: finding raypath on edge failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    tt = 0.0;
                    reachedTx = true;
                }
                if (onEdge){
                    edgeNodes[0] = edgeNodestmp[0];
                    edgeNodes[1] = edgeNodestmp[1];
                }
                prev_pt = curr_pt;
                curr_pt = pt_i;
                tt_s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                        faceNodes);
                tt += 0.5*(tt_s1 + tt_s2) * r_tmp.back().getDistance( curr_pt );
                tt_s1 = tt_s2;
                r_tmp.push_back( curr_pt );
                time=t_i-curr_pt.getDistance(prev_pt)*Slow;
            } else{ // on Face


                T2 cellNo1;
                if (!onFace || inlimitD)
                    cellNo1=cellNo;
                else
                    cellNo1=findAdjacentCell2(faceNodes, cellNo,curr_pt);

                sxyz<T1> pt_i;
                T1 Slow;
                std::array<T2,3> ind[8] = {
                    { { this->neighbors[cellNo][0], this->neighbors[cellNo][1], this->neighbors[cellNo][2] } },
                    { { this->neighbors[cellNo][0], this->neighbors[cellNo][1], this->neighbors[cellNo][3] } },
                    { { this->neighbors[cellNo][0], this->neighbors[cellNo][2], this->neighbors[cellNo][3] } },
                    { { this->neighbors[cellNo][1], this->neighbors[cellNo][2], this->neighbors[cellNo][3] } },
                    {{ this->neighbors[cellNo1][0], this->neighbors[cellNo1][1], this->neighbors[cellNo1][2] }},
                    { { this->neighbors[cellNo1][0], this->neighbors[cellNo1][1], this->neighbors[cellNo1][3]}},
                    { { this->neighbors[cellNo1][0], this->neighbors[cellNo1][2], this->neighbors[cellNo1][3] } },
                    { { this->neighbors[cellNo1][1], this->neighbors[cellNo1][2], this->neighbors[cellNo1][3] } }
                };
                for ( size_t n=0; n<8; ++n ) {
                    std::sort( ind[n].begin(), ind[n].end() );
                }
                // there are 3 faces that might be intersected

                T1 t_i=std::numeric_limits<T1>::max();
                T2 face;
                T1 s1=0.25*(nodes[this->neighbors[cellNo][0]].getNodeSlowness()+
                            nodes[this->neighbors[cellNo][1]].getNodeSlowness()+
                            nodes[this->neighbors[cellNo][2]].getNodeSlowness()+
                            nodes[this->neighbors[cellNo][3]].getNodeSlowness());
                T1 s2=0.25*(nodes[this->neighbors[cellNo1][0]].getNodeSlowness()+
                            nodes[this->neighbors[cellNo1][1]].getNodeSlowness()+
                            nodes[this->neighbors[cellNo1][2]].getNodeSlowness()+
                            nodes[this->neighbors[cellNo1][3]].getNodeSlowness());

                for ( size_t n=0; n<8; ++n ) {
                    if ( ind[n] == faceNodes) continue;
                    if (NearSource){
                        bool flag=false;
                        std::array<T1,3> Barycenter;
                        if (blti_solver_around_source(NodeSource, curr_pt, ind[n], Barycenter)==true){
                            pt_i.x=Barycenter[0]*nodes[ind[n][0]].getX()+Barycenter[1]*nodes[ind[n][1]].getX()+Barycenter[2]*nodes[ind[n][2]].getX();
                            pt_i.y=Barycenter[0]*nodes[ind[n][0]].getY()+Barycenter[1]*nodes[ind[n][1]].getY()+Barycenter[2]*nodes[ind[n][2]].getY();
                            pt_i.z=Barycenter[0]*nodes[ind[n][0]].getZ()+Barycenter[1]*nodes[ind[n][1]].getZ()+Barycenter[2]*nodes[ind[n][2]].getZ();

                            if(NodeSource.getDistance(pt_i)>NodeSource.getDistance(curr_pt))
                                continue;
                            for(T2 ni=0;ni<3;++ni){
                                if(std::abs(1.0-Barycenter[ni])<minDist*minDist){
                                    onNode = true;
                                    onEdge = false;
                                    onFace = false;
                                    faceNodesPrev=faceNodes;
                                    nodeNo=ind[n][ni];
                                    flag=true;
                                    break;
                                }
                            }
                            if (flag)
                                break;
                            for(T2 ni=0;ni<3;++ni){
                                if(std::abs(Barycenter[ni])<minDist*minDist){
                                    onNode = false;
                                    onEdge = true;
                                    onFace = false;
                                    edgeNodes[0] = ind[n][(ni+1)%3];
                                    edgeNodes[1] = ind[n][(ni+2)%3];
                                    faceNodesPrev=faceNodes;
                                    break;
                                }
                            }
                            if (flag)
                                break;
                            onNode = false;
                            onEdge = false;
                            onFace = true;
                            face=static_cast<T1>(n);
                            break;
                        }
                        continue;
                    }
                    T1 s;
                    sxyz<T1> pt;
                    s=(n<4)?s1:s2;

                    if(blti_raytrace(curr_pt, ind[n], pt, threadNo, s)==false) continue;
                    T1 t=Interpolator<T1>::bilinearTime(pt, nodes[ind[n][0]], nodes[ind[n][1]], nodes[ind[n][2]], threadNo);
                    if (t>time)
                        continue;
                    t+=curr_pt.getDistance(pt)*s;
                    if (t<t_i){
                        t_i=t;
                        Slow=s;
                        pt_i=pt;
                        onNode = false;
                        onEdge = false;
                        onFace = true;
                        face=static_cast<T1>(n);
                    }
                }
                if (!NearSource){
                    std::array<T2,2> indEdges[9] = {
                        {{ this->neighbors[cellNo][0],this->neighbors[cellNo][1]}},
                        {{ this->neighbors[cellNo][0],this->neighbors[cellNo][2]}},
                        {{ this->neighbors[cellNo][0],this->neighbors[cellNo][3]}},
                        {{ this->neighbors[cellNo][1],this->neighbors[cellNo][2]}},
                        {{ this->neighbors[cellNo][1],this->neighbors[cellNo][3]}},
                        {{ this->neighbors[cellNo][2],this->neighbors[cellNo][3]}}
                    };
                    T2 forthNode=0;
                    for (T2 i=0; i<4;++i){
                        if (this->neighbors[cellNo1][i]!=this->neighbors[cellNo][0] &&
                            this->neighbors[cellNo1][i]!=this->neighbors[cellNo][1]&&
                            this->neighbors[cellNo1][i]!=this->neighbors[cellNo][2]&&
                            this->neighbors[cellNo1][i]!=this->neighbors[cellNo][3]){
                            forthNode=i;
                            break;
                        }
                    }
                    indEdges[6]={this->neighbors[cellNo1][forthNode],this->neighbors[cellNo1][(forthNode+1)%4]};
                    indEdges[7]={this->neighbors[cellNo1][forthNode],this->neighbors[cellNo1][(forthNode+2)%4]};
                    indEdges[8]={this->neighbors[cellNo1][forthNode],this->neighbors[cellNo1][(forthNode+3)%4]};
                    for(T2 n=0;n<9;++n){
                        T1 s;
                        sxyz<T1> pt;
                        s=(n<6)?s1:s2;
                        if(blti2D_raytrace(curr_pt, indEdges[n][0], indEdges[n][1], pt, threadNo, s) ==false) continue;
                        T1 dist0=pt.getDistance(nodes[indEdges[n][0]]);
                        T1 dist1=pt.getDistance(nodes[indEdges[n][1]]);
                        T1 t= (dist1*nodes[indEdges[n][0]].getTT(threadNo)+dist0*nodes[indEdges[n][1]].getTT(threadNo))/(dist0+dist1);
                        if (t>time)
                            continue;
                        t+=curr_pt.getDistance(pt)*s;
                        if (t<t_i ){
                            t_i=t;
                            Slow=s;
                            pt_i=pt;
                            edgeNodes[0] = indEdges[n][0];
                            edgeNodes[1] = indEdges[n][1];
                            faceNodesPrev=faceNodes;
                            onNode = false;
                            onEdge = true;
                            onFace = false;

                        }

                    }
                    for (T2 n=0;n<4;++n){
                        T1 t=nodes[this->neighbors[cellNo][n]].getTT(threadNo);
                        if(t<time){
                            sxyz<T1> pt;
                            t+=s1*nodes[this->neighbors[cellNo][n]].getDistance(curr_pt);
                            if(t<t_i){
                                t_i=t;
                                Slow=s1;
                                pt_i=sxyz<T1>(nodes[this->neighbors[cellNo][n]]);
                                onNode = true;
                                onEdge = false;
                                onFace = false;
                                faceNodesPrev=faceNodes;
                                nodeNo=this->neighbors[cellNo][n];
                            }
                        }
                    }
                    if(nodes[this->neighbors[cellNo1][forthNode]].getTT(threadNo)<time){
                        sxyz<T1> pt;
                        if(nodes[this->neighbors[cellNo1][forthNode]].getTT(threadNo)+s2*nodes[this->neighbors[cellNo1][forthNode]].getDistance(curr_pt)<t_i){
                            t_i=nodes[this->neighbors[cellNo1][forthNode]].getTT(threadNo)+s2*nodes[this->neighbors[cellNo1][forthNode]].getDistance(curr_pt);
                            Slow=s2;
                            pt_i=sxyz<T1>(nodes[this->neighbors[cellNo1][forthNode]]);
                            onNode = true;
                            onEdge = false;
                            onFace = false;
                            faceNodesPrev=faceNodes;
                            nodeNo=this->neighbors[cellNo1][forthNode];
                        }
                    }
                }
                if (t_i==std::numeric_limits<T1>::max()){
                    std::cout << "\n\nWarning: finding raypath on face failed to converge for Rx "
                    << Rx.x << ' ' << Rx.y << ' ' << Rx.z << std::endl;
                    r_tmp.resize(1);
                    r_tmp[0] = Rx;
                    tt = 0.0;
                    reachedTx = true;
                    continue;
                }

                prev_pt = curr_pt;
                curr_pt = pt_i;
                inlimitD=false;
                if (std::abs(curr_pt.x-Xmax)<=minDist*minDist ||std::abs(curr_pt.x-Xmin)<=minDist*minDist )
                    inlimitD=true;
                if (std::abs(curr_pt.y-Ymax)<=minDist*minDist ||std::abs(curr_pt.y-Ymin)<=minDist*minDist )
                    inlimitD=true;
                if (std::abs(curr_pt.z-Zmax)<=minDist*minDist ||std::abs(curr_pt.z-Zmin)<=minDist*minDist )
                    inlimitD=true;

                if (onFace){
                    faceNodesPrev=faceNodes;
                    faceNodes = ind[face];
                    if (!inlimitD){
                        if (face<4)
                            cellNo=findAdjacentCell2(faceNodes, cellNo,curr_pt);
                        else
                            cellNo=findAdjacentCell2(faceNodes, cellNo1,curr_pt);
                    }
                }

                tt_s2 = computeSlowness(curr_pt, onNode, nodeNo, onEdge, edgeNodes,
                                        faceNodes);
                tt += 0.5*(tt_s1 + tt_s2) * r_tmp.back().getDistance( curr_pt );
                tt_s1 = tt_s2;
                r_tmp.push_back( curr_pt );
                time=t_i-curr_pt.getDistance(prev_pt)*Slow;
            }
            if ( onNode ) {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( curr_pt.getDistance( Tx[nt] ) < minDist*minDist ) {
                        tt += t0[nt];
                        reachedTx = true;
                        break;
                    }
                    for(size_t n=0;n<4;++n){
                        sxyz<T1> NearTx={nodes[this->neighbors[txCell[nt]][n]]};
                        if (curr_pt.getDistance(NearTx)< minDist*minDist){
                            tt += t0[nt];
                            reachedTx = true;
                            r_tmp.push_back( Tx[nt] );
                            break;
                        }
                    }

                }
            } else {
                for ( size_t nt=0; nt<Tx.size(); ++nt ) {
                    if ( txOnNode[nt] ) {
                        for ( auto nc=nodes[txNode[nt]].getOwners().begin(); nc!=nodes[txNode[nt]].getOwners().end(); ++nc ) {
                            if ( cellNo == *nc ) {
                                std::array<T2,4> itmp = getPrimary(cellNo);

                                if ( processVel )
                                    tt_s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                   nodes[itmp[0]],
                                                                                   nodes[itmp[1]],
                                                                                   nodes[itmp[2]],
                                                                                   nodes[itmp[3]]);
                                else
                                    tt_s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                                nodes[itmp[0]],
                                                                                nodes[itmp[1]],
                                                                                nodes[itmp[2]],
                                                                                nodes[itmp[3]]);

                                tt += t0[nt] + 0.5*(tt_s1 + tt_s2) * r_tmp.back().getDistance( Tx[nt] );
                                r_tmp.push_back( Tx[nt] );
                                reachedTx = true;
                                break;
                            }
                        }
                    } else {
                        if ( cellNo == txCell[nt] ) {
                            std::array<T2,4> itmp = getPrimary(cellNo);
                            if ( processVel )
                                tt_s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                               nodes[itmp[0]],
                                                                               nodes[itmp[1]],
                                                                               nodes[itmp[2]],
                                                                               nodes[itmp[3]]);
                            else
                                tt_s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                            nodes[itmp[0]],
                                                                            nodes[itmp[1]],
                                                                            nodes[itmp[2]],
                                                                            nodes[itmp[3]]);

                            tt += t0[nt] + 0.5*(tt_s1 + tt_s2) * r_tmp.back().getDistance( Tx[nt] );
                            r_tmp.push_back( Tx[nt] );
                            reachedTx = true;
                        } else {
                            for ( size_t nn=0; nn<txNeighborCells[nt].size(); ++nn ) {
                                if ( cellNo == txNeighborCells[nt][nn] ) {
                                    std::array<T2,4> itmp = getPrimary(cellNo);
                                    if ( processVel )
                                        tt_s2 = Interpolator<T1>::trilinearTriangleVel(Tx[nt],
                                                                                       nodes[itmp[0]],
                                                                                       nodes[itmp[1]],
                                                                                       nodes[itmp[2]],
                                                                                       nodes[itmp[3]]);
                                    else
                                        tt_s2 = Interpolator<T1>::trilinearTriangle(Tx[nt],
                                                                                    nodes[itmp[0]],
                                                                                    nodes[itmp[1]],
                                                                                    nodes[itmp[2]],
                                                                                    nodes[itmp[3]]);

                                    tt += t0[nt] + 0.5*(tt_s1 + tt_s2) * r_tmp.back().getDistance( Tx[nt] );
                                    r_tmp.push_back( Tx[nt] );
                                    reachedTx = true;
                                    break;
                                }
                            }
                        }
                    }
                    if ( reachedTx ) break;
                }
            }
        }
        // for inversion, the order should be from Tx to Rx, so we reorder...
        size_t npts = r_tmp.size();
        r_data.resize( npts );
        for ( size_t nn=0; nn<npts; ++nn ) {
            r_data[nn] = r_tmp[ npts-1-nn ];
        }
    }

    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::blti_raytrace(const sxyz<T1>& curr_pt,
                                             const std::array<T2,3>&face,
                                             sxyz<T1>& next_pt,
                                             const size_t threadNo,
                                             const T1& s)const{

        NODE *vertexA=&(nodes[face[0]]);
        NODE *vertexB=&(nodes[face[1]]);
        NODE *vertexC=&(nodes[face[2]]);

        if ( vertexA->getTT(threadNo) > vertexB->getTT(threadNo) )
            std::swap(vertexA, vertexB);
        if ( vertexA->getTT(threadNo) > vertexC->getTT(threadNo) )
            std::swap(vertexA, vertexC);

        T1 u = vertexB->getTT(threadNo) - vertexA->getTT(threadNo);
        T1 v = vertexC->getTT(threadNo) - vertexA->getTT(threadNo);
        sxyz<T1> v_b = { vertexC->getX() - vertexA->getX(),
            vertexC->getY() - vertexA->getY(),
            vertexC->getZ() - vertexA->getZ() };
        sxyz<T1> v_c = { vertexB->getX() - vertexA->getX(),
            vertexB->getY() - vertexA->getY(),
            vertexB->getZ() - vertexA->getZ() };

        sxyz<T1> v_n = cross(v_b, v_c);

        T1 b = norm( v_b );
        T1 c = norm( v_c );
        T1 d2 = dot(v_b, v_c);


        T1 phi=norm(v_n);
        T1 d_tmp = -vertexA->getX()*v_n.x - vertexA->getY()*v_n.y - vertexA->getZ()*v_n.z;

        T1 k = -(d_tmp + v_n.x*curr_pt.x + v_n.y*curr_pt.y+ v_n.z*curr_pt.z)/norm2(v_n);

        sxyz<T1> pt;
        pt.x = curr_pt.x+ k*v_n.x;
        pt.y = curr_pt.y + k*v_n.y;
        pt.z = curr_pt.z + k*v_n.z;

        T1 rho0 = curr_pt.getDistance( pt );

        // project point on AB
        sxyz<T1> v_pt = {pt.x-vertexA->getX(), pt.y-vertexA->getY(), pt.z-vertexA->getZ()};
        //// decomposition of Ap
        sxz<T1> AtA_Vect1={b*b,d2};
        sxz<T1> AtA_Vect2={d2,c*c};
        sxz<T1> Atb={dot(v_b,v_pt),dot(v_c,v_pt)};
        T1 DeT=det(AtA_Vect1,AtA_Vect2);
        T1 xi0=det(AtA_Vect1,Atb)/DeT;
        T1 zeta0=det(Atb,AtA_Vect2)/DeT;

        T1 beta = u*b*b - v*d2;
        T1 gamma = v*c*c - u*d2;
        T1 w_tilde2 = (s*s*phi*phi-u*u*b*b-v*v*c*c+2.0*u*v*d2 );
        if (w_tilde2>0.0){
            T1 xi_tilde = -std::abs(beta)*rho0/(phi*sqrt(w_tilde2));
            T1 zeta_tilde = -std::abs(gamma)*rho0/(phi*sqrt(w_tilde2));
            T1 xi = xi_tilde + xi0;
            T1 zeta = zeta_tilde + zeta0;
            if ( 0.<xi && xi<1. && 0.<zeta && zeta<1. && 0.<(xi+zeta) && (xi+zeta)<1. ){
                next_pt.x=xi*vertexB->getX()+zeta*vertexC->getX()+(1-xi-zeta)*vertexA->getX();
                next_pt.y=xi*vertexB->getY()+zeta*vertexC->getY()+(1-xi-zeta)*vertexA->getY();
                next_pt.z=xi*vertexB->getZ()+zeta*vertexC->getZ()+(1-xi-zeta)*vertexA->getZ();
                return true;
            }else
                return false;

        }else
            return false;
    }

    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::blti2D_raytrace(const sxyz<T1>& curr_pt,
                                               const T2& node1,
                                               const T2& node2,
                                               sxyz<T1>& next_pt,
                                               const size_t threadNo,
                                               const T1& s) const{

        NODE *vertexA=&(nodes[node1]);
        NODE *vertexB=&(nodes[node2]);
        if ( vertexA->getTT(threadNo) > vertexB->getTT(threadNo) ) {
            std::swap(vertexA, vertexB);
        }
        T1 u = vertexB->getTT(threadNo) - vertexA->getTT(threadNo);

        sxyz<T1> v_b = { curr_pt.x- vertexA->getX(),
            curr_pt.y- vertexA->getY(),
            curr_pt.z- vertexA->getZ() };
        sxyz<T1> v_c = { vertexB->getX() - vertexA->getX(),
            vertexB->getY() - vertexA->getY(),
            vertexB->getZ() - vertexA->getZ() };

        T1 c = norm( v_c );

        T1 w2 = (s*s*c*c - u*u);
        if (w2<0.0)
            return false;
        T1 w = sqrt( w2 );
        T1 k = dot(v_b,v_c)/dot(v_c,v_c);
        sxyz<T1> pt;
        pt.x = vertexA->getX() + k*v_c.x;
        pt.y = vertexA->getY() + k*v_c.y;
        pt.z = vertexA->getZ() + k*v_c.z;

        T1 rho0 = curr_pt.getDistance( pt );

        T1 xi0 = k;
        T1 xi = xi0 -u*rho0/(w*c);

        if ( 0.0<xi && xi<1.0 ) {
            next_pt.x=(1-xi)*vertexA->getX()+xi*vertexB->getX();
            next_pt.y=(1-xi)*vertexA->getY()+xi*vertexB->getY();
            next_pt.z=(1-xi)*vertexA->getZ()+xi*vertexB->getZ();
            return true;
        }
        return false;
    }

    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::blti_solver_around_source(const sxyz<T1>& source,
                                                         const sxyz<T1>& curr_pt,
                                                         const std::array<T2,3>& face,
                                                         std::array<T1,3>& barycenters) const {
        T2 node1=face[0];
        T2 node2=face[1];
        T2 node3=face[2];
        sxyz<T1> PQ=source-curr_pt;
        sxyz<T1>PA={nodes[node1].getX()-curr_pt.x,nodes[node1].getY()-curr_pt.y,nodes[node1].getZ()-curr_pt.z};
        sxyz<T1>PB={nodes[node2].getX()-curr_pt.x,nodes[node2].getY()-curr_pt.y,nodes[node2].getZ()-curr_pt.z};
        sxyz<T1>PC={nodes[node3].getX()-curr_pt.x,nodes[node3].getY()-curr_pt.y,nodes[node3].getZ()-curr_pt.z};
        sxyz<T1> m=cross(PQ,PC);
        T1 u=dot(PB,m);
        T1 v=-dot(PA,m);
        T1 relativeZero=min_dist;
        if(signum(u)!=signum(v) && fabs(v)>relativeZero && fabs(u)>relativeZero)
            return false;
        T1 w=tripleScalar(PQ, PB, PA);
        if(signum(u)!=signum(w) && abs(w)>relativeZero && abs(u)>relativeZero)
            return false;
        T1 denom=1.0/(u+v+w);
        barycenters[0]=denom*u;
        barycenters[1]=denom*v;
        barycenters[2]=denom*w;
        return true;
    }


    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::check_pt_location(sxyz<T1> &curr_pt,
                                                 const std::array<T2,3> &ind,
                                                 bool &onNode,
                                                 T2 &nodeNo,
                                                 bool &onEdge,
                                                 std::array<T2,2> &edgeNodes,
                                                 bool &onFace,
                                                 std::array<T2,3> &faceNodes) const {
        // check if point is on vertex, edge or face

        bool break_flag = false;

        for ( size_t n=0; n<3; ++n ) {
            if ( nodes[ ind[n] ].getDistance( curr_pt ) < min_dist) {
                curr_pt = nodes[ ind[n] ];
                nodeNo = ind[n];
                onNode = true;
                onEdge = false;
                onFace = false;
                break_flag = true;
                return break_flag;
            }
        }

        for ( size_t n1=0; n1<3; ++n1 ) {
            size_t n2 = (n1+1)%3;
            if ( areCollinear(curr_pt, nodes[ind[n1]], nodes[ind[n2]]) ) {
                edgeNodes[0] = ind[n1];
                edgeNodes[1] = ind[n2];
                onNode = false;
                onEdge = true;
                onFace = false;
                break_flag = true;
                return break_flag;
            }
        }

        onNode = false;
        onEdge = false;
        onFace = true;

        faceNodes = ind;

        return break_flag;
    }


    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::check_pt_location(sxyz<T1> &curr_pt,
                                                 const std::vector<T2> &ind1,
                                                 const std::array<T2,3> &ind2,
                                                 bool &onNode,
                                                 T2 &nodeNo,
                                                 bool &onEdge,
                                                 std::array<T2,2> &edgeNodes,
                                                 bool &onFace,
                                                 std::array<T2,3> &faceNodes) const {
        // chech if point is on vertex, edge or face

        bool break_flag = false;

        for ( size_t n=0; n<4; ++n ) {
            if ( nodes[ ind1[n] ].getDistance( curr_pt ) < min_dist ) {
                curr_pt = nodes[ ind1[n] ];
                nodeNo = ind1[n];
                onNode = true;
                onEdge = false;
                onFace = false;
                break_flag = true;
                return break_flag;
            }
        }

        for ( size_t n1=0; n1<3; ++n1 ) {
            size_t n2 = (n1+1)%3;
            if ( areCollinear(curr_pt, nodes[ind2[n1]], nodes[ind2[n2]]) ) {
                edgeNodes[0] = ind2[n1];
                edgeNodes[1] = ind2[n2];
                onNode = false;
                onEdge = true;
                onFace = false;
                break_flag = true;
                return break_flag;
            }
        }

        onNode = false;
        onEdge = false;
        onFace = true;

        faceNodes = ind2;
        std::sort(faceNodes.begin(), faceNodes.end());

        return break_flag;
    }


    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::update_m_data(std::vector<sijv<T1>>& m_data,
                                             sijv<T1>& m,
                                             const std::set<T2>& allNodes,
                                             const sxyz<T1>& mid_pt,
                                             const T1 ds) const {
        // valid for slowness
        std::vector<T1> w;
        T1 sum_w = 0.0;
        for ( auto it=allNodes.begin(); it!=allNodes.end(); ++it ) {
            w.push_back( 1./nodes[*it].getDistance( mid_pt ) );
            sum_w += w.back();
        }
        size_t nn=0;
        for ( auto it=allNodes.begin(); it!=allNodes.end(); ++it ) {
            m.j = *it;
            m.v = ds * w[nn++]/sum_w;
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


    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::update_m_data(std::vector<sijv<T1>>& m_data,
                                             sijv<T1>& m,
                                             const std::set<T2>& allNodes,
                                             const sxyz<T1>& mid_pt,
                                             const T1 ds,
                                             const T1 s_sq) const {
        // valid for velocity
        std::vector<T1> w;
        T1 sum_w = 0.0;
        for ( auto it=allNodes.begin(); it!=allNodes.end(); ++it ) {
            w.push_back( 1./nodes[*it].getDistance( mid_pt ) );
            sum_w += w.back();
        }
        size_t nn=0;
        for ( auto it=allNodes.begin(); it!=allNodes.end(); ++it ) {
            m.j = *it;
            m.v = -s_sq * ds * w[nn++]/sum_w;
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


    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::intersectVecTriangle(const T2 iO, const sxyz<T1> &vec,
                                                    const T2 iA, T2 iB, T2 iC,
                                                    sxyz<T1> &pt_i) const {

        sxyz<T1> OA = {nodes[iA].getX()-nodes[iO].getX(), nodes[iA].getY()-nodes[iO].getY(), nodes[iA].getZ()-nodes[iO].getZ()};
        // check if counterclockwise
        sxyz<T1> AB = {nodes[iB].getX()-nodes[iA].getX(),
            nodes[iB].getY()-nodes[iA].getY(),
            nodes[iB].getZ()-nodes[iA].getZ()};
        sxyz<T1> AC = {nodes[iC].getX()-nodes[iA].getX(),
            nodes[iC].getY()-nodes[iA].getY(),
            nodes[iC].getZ()-nodes[iA].getZ()};
        sxyz<T1> n = cross(AB, AC);
        if ( dot(OA, n) > 0. ) std::swap(iB, iC);

        sxyz<T1> OB = {nodes[iB].getX()-nodes[iO].getX(), nodes[iB].getY()-nodes[iO].getY(), nodes[iB].getZ()-nodes[iO].getZ()};
        sxyz<T1> OC = {nodes[iC].getX()-nodes[iO].getX(), nodes[iC].getY()-nodes[iO].getY(), nodes[iC].getZ()-nodes[iO].getZ()};

        T1 u, v, w;
        u = tripleScalar(vec, OC, OB);
        if ( u<0.0 ) return false;
        v = tripleScalar(vec, OA, OC);
        if ( v<0.0 ) return false;
        w = tripleScalar(vec, OB, OA);
        if ( w<0.0 ) return false;

        T1 den = 1./(u+v+w);
        u *= den;
        v *= den;
        w *= den;

        pt_i.x = u*nodes[iA].getX() + v*nodes[iB].getX() + w*nodes[iC].getX();
        pt_i.y = u*nodes[iA].getY() + v*nodes[iB].getY() + w*nodes[iC].getY();
        pt_i.z = u*nodes[iA].getZ() + v*nodes[iB].getZ() + w*nodes[iC].getZ();

        return true;
    }

    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::intersectVecTriangle(const sxyz<T1> &O, const sxyz<T1> &vec,
                                                    const T2 iA, T2 iB, T2 iC,
                                                    sxyz<T1> &pt_i) const {

        sxyz<T1> OA = {nodes[iA].getX()-O.x, nodes[iA].getY()-O.y, nodes[iA].getZ()-O.z};
        // check if counterclockwise
        sxyz<T1> AB = {nodes[iB].getX()-nodes[iA].getX(),
            nodes[iB].getY()-nodes[iA].getY(),
            nodes[iB].getZ()-nodes[iA].getZ()};
        sxyz<T1> AC = {nodes[iC].getX()-nodes[iA].getX(),
            nodes[iC].getY()-nodes[iA].getY(),
            nodes[iC].getZ()-nodes[iA].getZ()};
        sxyz<T1> n = cross(AB, AC);
        if ( dot(OA, n) > 0. ) std::swap(iB, iC);

        sxyz<T1> OB = {nodes[iB].getX()-O.x, nodes[iB].getY()-O.y, nodes[iB].getZ()-O.z};
        sxyz<T1> OC = {nodes[iC].getX()-O.x, nodes[iC].getY()-O.y, nodes[iC].getZ()-O.z};

        T1 u, v, w;
        u = tripleScalar(vec, OC, OB);
        if ( u<0.0 ) return false;
        v = tripleScalar(vec, OA, OC);
        if ( v<0.0 ) return false;
        w = tripleScalar(vec, OB, OA);
        if ( w<0.0 ) return false;

        T1 den = 1./(u+v+w);
        u *= den;
        v *= den;
        w *= den;

        pt_i.x = u*nodes[iA].getX() + v*nodes[iB].getX() + w*nodes[iC].getX();
        pt_i.y = u*nodes[iA].getY() + v*nodes[iB].getY() + w*nodes[iC].getY();
        pt_i.z = u*nodes[iA].getZ() + v*nodes[iB].getZ() + w*nodes[iC].getZ();

        return true;
    }

    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::intersectVecEdge(const sxyz<T1>& curr_pt,
                                                const sxyz<T1>& g,
                                                std::array<T2,3>& faceNodes,
                                                sxyz<T1>&  pt_i,
                                                std::array<T2,2>& edgeNodes) const {

        // from http://mathworld.wolfram.com/Line-LineIntersection.html

#ifdef DEBUG_RP
        if ( verbose > 1 ) {
            std::cout << "\n\n\n\n\nIn intersectVecEdge (face)\n\n";

            std::cout << "fig = plt.figure()\nax = fig.add_subplot(111, projection='3d')\n";
            std::cout << "cpt = np.array([" << curr_pt.x << ", " << curr_pt.y << ", " << curr_pt.z << "])\n";
            std::cout << "g = np.array([" << g.x << ", " << g.y << ", " << g.z << "])\n";
            std::cout << "pt2 = cpt + 0.05*g\n";
            std::cout << "c1 = np.array([" << nodes[faceNodes[0]].getX() << ", " << nodes[faceNodes[0]].getY() << ", " << nodes[faceNodes[0]].getZ() << "])\n";
            std::cout << "c2 = np.array([" << nodes[faceNodes[1]].getX() << ", " << nodes[faceNodes[1]].getY() << ", " << nodes[faceNodes[1]].getZ() << "])\n";
            std::cout << "c3 = np.array([" << nodes[faceNodes[2]].getX() << ", " << nodes[faceNodes[2]].getY() << ", " << nodes[faceNodes[2]].getZ() << "])\n";

            std::cout << "ax.plot([c1[0], c2[0], c3[0], c1[0]], [c1[1], c2[1], c3[1], c1[1]], [c1[2], c2[2], c3[2], c1[2]])\n";
            std::cout << "ax.plot([cpt[0], pt2[0]], [cpt[1], pt2[1]], [cpt[2], pt2[2]], c='r')\n";
            std::cout << "ax.scatter(cpt[0], cpt[1], cpt[2], c='r')\n";
        }
#endif

        sxyz<T1> x1 = {nodes[faceNodes[0]].getX(),
            nodes[faceNodes[0]].getY(),
            nodes[faceNodes[0]].getZ()};
        sxyz<T1> x2 = {nodes[faceNodes[1]].getX(),
            nodes[faceNodes[1]].getY(),
            nodes[faceNodes[1]].getZ()};
        sxyz<T1> x4 = curr_pt + static_cast<T1>(10.0)*x1.getDistance(x2) * g;

        sxyz<T1> a = x2 - x1;
        sxyz<T1> b = x4 - curr_pt;
        sxyz<T1> c = curr_pt - x1;

        sxyz<T1> ab = cross(a, b);

        pt_i = x1 + (dot(cross(c, b), ab) / norm2(ab)) * a;

#ifdef DEBUG_RP
        if ( verbose > 1 ) {
            std::cout << "pti = np.array([" << pt_i.x << ", " << pt_i.y << ", " << pt_i.z << "])\n";
            std::cout << "ax.scatter(pti[0], pti[1], pti[2], c='g')\n";
        }
#endif

        sxyz<T1> mid_pt = static_cast<T1>(0.5) * ( curr_pt + pt_i );

        bool test1 = testInTriangle(&(nodes[ faceNodes[0] ]),
                                    &(nodes[ faceNodes[1] ]),
                                    &(nodes[ faceNodes[2] ]), mid_pt);

        // check if we are between x1 & x2

        b = pt_i - x1;
        T1 dab = dot(a, b);
        bool test2 = dab > 0.0 && dab <= norm2(a);

        // check if going in the same direction as g

        b = pt_i - curr_pt;
        bool test3 = dot(b, g) > 0.0;

        if ( test1 && test2 && test3) {
            edgeNodes[0] = faceNodes[0];
            edgeNodes[1] = faceNodes[1];
            return true;
        }

        x1 = {nodes[faceNodes[0]].getX(),
            nodes[faceNodes[0]].getY(),
            nodes[faceNodes[0]].getZ()};
        x2 = {nodes[faceNodes[2]].getX(),
            nodes[faceNodes[2]].getY(),
            nodes[faceNodes[2]].getZ()};

        a = x2 - x1;
        b = x4 - curr_pt;
        c = curr_pt - x1;

        ab = cross(a, b);

        pt_i = x1 + (dot(cross(c, b), ab) / norm2(ab)) * a;

#ifdef DEBUG_RP
        if ( verbose > 1 ) {
            std::cout << "pti = np.array([" << pt_i.x << ", " << pt_i.y << ", " << pt_i.z << "])\n";
            std::cout << "ax.scatter(pti[0], pti[1], pti[2], c='k')\n";
        }
#endif

        mid_pt = static_cast<T1>(0.5) * ( curr_pt + pt_i );

        test1 = testInTriangle(&(nodes[ faceNodes[0] ]),
                               &(nodes[ faceNodes[1] ]),
                               &(nodes[ faceNodes[2] ]), mid_pt);
        // check if we are between x1 & x2

        b = pt_i - x1;
        dab = dot(a, b);
        test2 = dab > 0.0 && dab <= norm2(a);

        b = pt_i - curr_pt;
        test3 = dot(b, g) > 0.0;

        if ( test1 && test2 && test3 ) {
            edgeNodes[0] = faceNodes[0];
            edgeNodes[1] = faceNodes[2];
            return true;
        }

        x1 = {nodes[faceNodes[1]].getX(),
            nodes[faceNodes[1]].getY(),
            nodes[faceNodes[1]].getZ()};
        x2 = {nodes[faceNodes[2]].getX(),
            nodes[faceNodes[2]].getY(),
            nodes[faceNodes[2]].getZ()};

        a = x2 - x1;
        b = x4 - curr_pt;
        c = curr_pt - x1;

        ab = cross(a, b);

        pt_i = x1 + (dot(cross(c, b), ab) / norm2(ab)) * a;

#ifdef DEBUG_RP
        if ( verbose > 1 ) {
            std::cout << "pti = np.array([" << pt_i.x << ", " << pt_i.y << ", " << pt_i.z << "])\n";
            std::cout << "ax.scatter(pti[0], pti[1], pti[2], c='b')\n";
            std::cout << "plt.show()\n\n";
        }
#endif

        mid_pt = static_cast<T1>(0.5) * ( curr_pt + pt_i );

        test1 = testInTriangle(&(nodes[ faceNodes[0] ]),
                               &(nodes[ faceNodes[1] ]),
                               &(nodes[ faceNodes[2] ]), mid_pt);
        b = pt_i - x1;
        dab = dot(a, b);
        test2 = dab > 0.0 && dab <= norm2(a);

        b = pt_i - curr_pt;
        test3 = dot(b, g) > 0.0;

        if ( test1 && test2 && test3 ) {
            edgeNodes[0] = faceNodes[1];
            edgeNodes[1] = faceNodes[2];
            return true;
        }

        return false;
    }

    template<typename T1, typename T2, typename NODE>
    T2 Grid3Dun<T1,T2,NODE>::findAdjacentCell1(const std::array<T2,3> &faceNodes,
                                               const T2 nodeNo) const {

        std::vector<T2> cells;
        for ( auto nc0=nodes[faceNodes[0]].getOwners().begin(); nc0!=nodes[faceNodes[0]].getOwners().end(); ++nc0 ) {
            if ( std::find(nodes[faceNodes[1]].getOwners().begin(), nodes[faceNodes[1]].getOwners().end(), *nc0)!=nodes[faceNodes[1]].getOwners().end() &&
                std::find(nodes[faceNodes[2]].getOwners().begin(), nodes[faceNodes[2]].getOwners().end(), *nc0)!=nodes[faceNodes[2]].getOwners().end() )
                cells.push_back( *nc0 );
        }
        if ( cells.size() == 1 ) {
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

    template<typename T1, typename T2, typename NODE>
    T2 Grid3Dun<T1,T2,NODE>::findAdjacentCell2(const std::array<T2,3> &faceNodes,
                                               const T2 cellNo) const {

        std::vector<T2> cells;
        for ( auto nc0=nodes[faceNodes[0]].getOwners().begin(); nc0!=nodes[faceNodes[0]].getOwners().end(); ++nc0 ) {
            if ( std::find(nodes[faceNodes[1]].getOwners().begin(), nodes[faceNodes[1]].getOwners().end(), *nc0)!=nodes[faceNodes[1]].getOwners().end() &&
                std::find(nodes[faceNodes[2]].getOwners().begin(), nodes[faceNodes[2]].getOwners().end(), *nc0)!=nodes[faceNodes[2]].getOwners().end() )
                cells.push_back( *nc0 );
        }
        if ( cells.size() == 1 ) {
            return cells[0];
        }
        if ( cellNo == cells[0] ) {
            return cells[1];
        } else if ( cellNo == cells[1] ) {
            return cells[0];
        }
        return std::numeric_limits<T2>::max();
    }

    template<typename T1, typename T2, typename NODE>
    T2 Grid3Dun<T1,T2,NODE>::findAdjacentCell2(const std::array<T2,3> &faceNodes,
                                               const T2 & cellNo,
                                               const sxyz<T1>& curr_pt ) const {

        T2 ac = findAdjacentCell2(faceNodes, cellNo);
        if ( ac!=cellNo && ac != std::numeric_limits<T2>::max() )
            return ac;
        std::set<T2> AdjacentCells;
        for ( size_t n1=0; n1<3; ++n1 ) {
            size_t n2((n1+1)%3), n3((n1+2)%3);
            for ( auto nc0=nodes[faceNodes[n1]].getOwners().begin(); nc0!=nodes[faceNodes[n1]].getOwners().end(); ++nc0 ) {
                for (size_t N=0; N<4; ++N ) {
                    for ( auto nc=nodes[this->neighbors[(*nc0)][N]].getOwners().begin(); nc!=nodes[this->neighbors[(*nc0)][N]].getOwners().end(); ++nc ) {
                        for (size_t iD=0; iD<4; ++iD ) {
                            size_t iA((iD+1)%4), iB((iD+2)%4), iC((iD+3)%4);
                            bool coarsetofine=(testInTriangle(&nodes[faceNodes[n1]],
                                                              &nodes[faceNodes[n2]],
                                                              &nodes[faceNodes[n3]],
                                                              sxyz<T1>(nodes[this->neighbors[(*nc)][iA]]))) &&
                            (testInTriangle(&nodes[faceNodes[n1]],
                                            &nodes[faceNodes[n2]],
                                            &nodes[faceNodes[n3]],
                                            sxyz<T1>(nodes[this->neighbors[(*nc)][iB]]))) &&
                            (testInTriangle(&nodes[faceNodes[n1]],
                                            &nodes[faceNodes[n2]],
                                            &nodes[faceNodes[n3]],
                                            sxyz<T1>(nodes[this->neighbors[(*nc)][iC]]))) && (*nc!=cellNo);

                            bool fintocoarse=(testInTriangle(&nodes[this->neighbors[(*nc)][iA]],
                                                             &nodes[this->neighbors[(*nc)][iB]],
                                                             &nodes[this->neighbors[(*nc)][iC]],
                                                             sxyz<T1>(nodes[faceNodes[n1]]))) &&
                            (testInTriangle(&nodes[this->neighbors[(*nc)][iA]],
                                            &nodes[this->neighbors[(*nc)][iB]],
                                            &nodes[this->neighbors[(*nc)][iC]],
                                            sxyz<T1>(nodes[faceNodes[n2]]))) &&
                            (testInTriangle(&nodes[this->neighbors[(*nc)][iA]],
                                            &nodes[this->neighbors[(*nc)][iB]],
                                            &nodes[this->neighbors[(*nc)][iC]],
                                            sxyz<T1>(nodes[faceNodes[n3]]))) && (*nc!=cellNo);

                            if(fintocoarse || coarsetofine){
                                AdjacentCells.insert(*(nc));
                                break;
                            }
                        }
                    }
                }
            }
        }



        for ( auto nc=AdjacentCells.begin();nc!=AdjacentCells.end();++nc ) {
            for ( T2 iD=0; iD<4; ++iD ) {
                T2 iA((iD+1)%4), iB((iD+2)%4), iC((iD+3)%4);
                if (testInTriangle(&nodes[this->neighbors[*nc][iA]],
                                   &nodes[this->neighbors[*nc][iB]],
                                   &nodes[this->neighbors[*nc][iC]],
                                   curr_pt)){
                    return (*nc);
                }
            }
        }
        return std::numeric_limits<T2>::max();
    }


    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::plotCell(const T2 cellNo, const sxyz<T1> &pt, const sxyz<T1> &g) const {


        if ( true ) {
            std::array<T2,4> itmp = getPrimary(cellNo);
            T2 i0 = itmp[0];
            T2 i1 = itmp[1];
            T2 i2 = itmp[2];
            T2 i3 = itmp[3];

            std::cout << "\nplot3(["<<nodes[ i0 ].getX()<<' ' << nodes[ i1 ].getX() <<"],["
            <<nodes[ i0 ].getY()<<' ' << nodes[ i1 ].getY() <<"],["
            <<nodes[ i0 ].getZ()<<' ' << nodes[ i1 ].getZ() <<"]); hold on;\n";
            std::cout << "plot3(["<<nodes[ i0 ].getX()<<' ' << nodes[ i2 ].getX() <<"],["
            <<nodes[ i0 ].getY()<<' ' << nodes[ i2 ].getY() <<"],["
            <<nodes[ i0 ].getZ()<<' ' << nodes[ i2 ].getZ() <<"])\n";
            std::cout << "plot3(["<<nodes[ i0 ].getX()<<' ' << nodes[ i3 ].getX() <<"],["
            <<nodes[ i0 ].getY()<<' ' << nodes[ i3 ].getY() <<"],["
            <<nodes[ i0 ].getZ()<<' ' << nodes[ i3 ].getZ() <<"])\n";
            std::cout << "plot3(["<<nodes[ i1 ].getX()<<' '<<nodes[ i2 ].getX()<<' '<<nodes[ i3 ].getX()<<' '<<nodes[ i1 ].getX()<<"],["
            <<nodes[ i1 ].getY()<<' '<<nodes[ i2 ].getY()<<' '<<nodes[ i3 ].getY()<<' '<<nodes[ i1 ].getY()<<"],["
            <<nodes[ i1 ].getZ()<<' '<<nodes[ i2 ].getZ()<<' '<<nodes[ i3 ].getZ()<<' '<<nodes[ i1 ].getZ()<<"])\n";
            std::cout << "plot3(["<<pt.x<< ' ' << pt.x+g.x<<"],["<<pt.y<< ' ' << pt.y+g.y<<"],["<<pt.z<< ' ' << pt.z+g.z<<"],'r')\naxis equal\n\n";
        }
    }

    template<typename T1, typename T2, typename NODE>
    template<typename SetT>
    void Grid3Dun<T1,T2,NODE>::getNeighborNodes(const T2 cellNo,
                                                SetT &nnodes) const {

        for ( size_t n=0; n<this->neighbors[cellNo].size(); ++n ) {
            if ( nodes[this->neighbors[cellNo][n]].isPrimary() ) {
                T2 nodeNo = this->neighbors[cellNo][n];
                nnodes.insert( &(nodes[nodeNo]) );
                if ( rp_method == 1 ) {  // second-order
                    for ( auto nc=nodes[nodeNo].getOwners().cbegin(); nc!=nodes[nodeNo].getOwners().cend(); ++nc ) {
                        for ( size_t nn=0; nn<this->neighbors[*nc].size(); ++nn ) {
                            if ( nodes[ this->neighbors[*nc][nn] ].isPrimary() ) {
                                nnodes.insert( &(nodes[ this->neighbors[*nc][nn] ]) );
                            }
                        }
                    }
                }
            }
        }
        assert(nnodes.size()>0);
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::getNeighborNodesAB(const std::vector<NODE*>& ref_pt,
                                                  std::vector<std::vector<std::array<NODE*,3>>>& opp_pts) const {
        opp_pts.resize( ref_pt.size() );
        for ( size_t nr=0; nr<ref_pt.size(); ++nr ) {
            opp_pts[nr].resize( ref_pt[nr]->getOwners().size() );
            for ( size_t nc=0; nc<ref_pt[nr]->getOwners().size(); ++nc ) {
                T2 cellNo = ref_pt[nr]->getOwners()[nc];
                size_t ind=0;
                for ( size_t nn=0; nn<this->neighbors[cellNo].size(); ++nn ) {
                    if ( &(nodes[ this->neighbors[cellNo][nn] ]) != ref_pt[nr] && nodes[ this->neighbors[cellNo][nn] ].isPrimary() ) {
                        opp_pts[nr][nc][ind++] = &(nodes[ this->neighbors[cellNo][nn] ]);
                    }
                }
            }
        }
    }


    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dun<T1,T2,NODE>::computeSlowness(sxyz<T1> pt, const bool isTranslated) const {

        // Calculate the slowness of any point that is not on a node

        if (this->translateOrigin == true && isTranslated == false) {
            pt -= this->origin;
            checkPts(std::vector<sxyz<T1>> {pt}, true);
        }
        T2 cellNo = this->getCellNo( pt );

        std::vector<NODE*> interpNodes;

        for (size_t n=0; n < this->neighbors[ cellNo ].size(); n++){
            if ( nodes[this->neighbors[ cellNo ][n] ].isPrimary() ){
                interpNodes.push_back( &(nodes[this->neighbors[ cellNo ][n] ]) );
            }
        }
        if ( processVel )
            return Interpolator<T1>::trilinearTriangleVel( pt, interpNodes );
        else
            return Interpolator<T1>::trilinearTriangle( pt, interpNodes );
    }


    template<typename T1, typename T2, typename NODE>
    T1 Grid3Dun<T1,T2,NODE>::computeSlowness(const sxyz<T1>& curr_pt,
                                             const bool onNode,
                                             const T2 nodeNo,
                                             const bool onEdge,
                                             const std::array<T2,2>& edgeNodes,
                                             const std::array<T2,3>& faceNodes) const {
        if ( onNode ) {
            return nodes[nodeNo].getNodeSlowness();
        } else {
            if ( processVel ) {
                if ( onEdge ) {
                    return Interpolator<T1>::linearVel(curr_pt,
                                                       nodes[edgeNodes[0]],
                                                       nodes[edgeNodes[1]]);
                } else {
                    return Interpolator<T1>::bilinearTriangleVel(curr_pt,
                                                                 nodes[faceNodes[0]],
                                                                 nodes[faceNodes[1]],
                                                                 nodes[faceNodes[2]]);
                }
            } else {
                if ( onEdge ) {
                    return Interpolator<T1>::linear(curr_pt,
                                                    nodes[edgeNodes[0]],
                                                    nodes[edgeNodes[1]]);
                } else {
                    return Interpolator<T1>::bilinearTriangle(curr_pt,
                                                              nodes[faceNodes[0]],
                                                              nodes[faceNodes[1]],
                                                              nodes[faceNodes[2]]);
                }
            }
        }
    }

    template<typename T1, typename T2, typename NODE>
    sxyz<T1> Grid3Dun<T1,T2,NODE>::projectOnFace(const sxyz<T1>& curr_pt,
                                             const sxyz<T1>& g,
                                             std::array<T2,2>& edgeNodes,
                                             const std::vector<T2> cells,
                                             sxyz<T1>&  pt_i,
                                             const size_t & threadNo) const {

        // find all primary nodes of common cells
        std::vector<std::array<T2,4>> primary(cells.size());
        for (size_t nc=0; nc<cells.size(); ++nc) {
            size_t np = 0;
            for ( auto nn=this->neighbors[cells[nc]].begin(); nn!= this->neighbors[cells[nc]].end(); ++nn ) {
                if (nodes[*nn].isPrimary()) {
                    primary[nc][np++] = *nn;
                }
            }
        }
        // at the surface, the edge is common to two triangles
        // we must have only two other nodes at the surface, which are not common to the cells
        std::array<T2,2> edgeNodes2;
        size_t ne = 0;
        for (size_t nc=0; nc<cells.size(); ++nc) {
            for (size_t np=0; np<4; ++np) {
                if (primary[nc][np] == edgeNodes[0] || primary[nc][np] == edgeNodes[1]) {
                    continue;
                }
                bool found=false;
                for (size_t nc2=0; nc2<cells.size(); ++nc2) {
                    if ( nc == nc2 ) {
                        continue;
                    }
                    if ( std::find(primary[nc2].begin(), primary[nc2].end(), primary[nc][np])!=primary[nc2].end() ) {
                        found = true;  // primary node is common to two cells
                        break;
                    }
                }
                if ( found == false ) {
                    edgeNodes2[ne++] = primary[nc][np];
                }
            }
        }
        sxyz<T1> g_proj;
        T2 thirdnode;//the third node of the face on which g is to be projected

        if (nodes[edgeNodes2[0]].getTT(threadNo) < nodes[edgeNodes2[1]].getTT(threadNo)) {
            thirdnode = edgeNodes2[0];
        } else {
            thirdnode = edgeNodes2[1];
        }
        std::array<T2,2> oldEdgeNodes = edgeNodes;
        sxyz<T1> tmp = projectOnFace(g, {edgeNodes[0], edgeNodes[1], thirdnode});
        // find the intersection between tmp and the face edges
        for (size_t i = 0; i < 2; i++) {
            sxyz<T1> x1 = {nodes[edgeNodes[i]].getX(),
                nodes[edgeNodes[i]].getY(),
                nodes[edgeNodes[i]].getZ()};
            sxyz<T1> x2 = {nodes[thirdnode].getX(),
                nodes[thirdnode].getY(),
                nodes[thirdnode].getZ()};
            sxyz<T1> x4 = curr_pt + static_cast<T1>(10.0)*x1.getDistance(x2) * tmp;

            sxyz<T1> a = x2 - x1;
            sxyz<T1> b = x4 - curr_pt;
            sxyz<T1> c = curr_pt - x1;

            sxyz<T1> ab = cross(a, b);
            T1 p = dot(cross(c, b), ab) / norm2(ab);
            if(p >= 0. && p <= 1.){// 0 <= p <= 1 is a sufficient condition to check if pt_i0 is between x1 and x2 and that the midpoint of the line segment [curr_p, pt_i0] is inside the face
                sxyz<T1> pt_i0 = p * x2 + static_cast<T1>(1.0 -p) * x1;
                pt_i = pt_i0;
                edgeNodes[0] = edgeNodes[i];
                edgeNodes[1] = thirdnode;
                g_proj = pt_i0 - curr_pt;
                break;

            }
        }
        // if tmp points outside the tested edges or g_proj and g have opposite direction, we use the head wave traveltimes to find g_proj
        if ((g_proj.x == 0. && g_proj.y == 0. && g_proj.z == 0.) ||
            (dot(g,g_proj) < 0 )){
            T2 minttNode;
            minttNode = nodes[oldEdgeNodes[0]].getTT(threadNo) < nodes[oldEdgeNodes[1]].getTT(threadNo)? oldEdgeNodes[0]:oldEdgeNodes[1];
            minttNode = nodes[minttNode].getTT(threadNo) < nodes[thirdnode].getTT(threadNo)? minttNode:thirdnode;
            pt_i = nodes[minttNode];
            edgeNodes[0] = minttNode;
            g_proj = pt_i - curr_pt;

        }
        return g_proj;
    }
    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::projectOnFace(const sxyz<T1>& curr_pt,
                                             const T2 nodeNo,
                                             sxyz<T1>& g,
                                             std::array<T2,2>& edgeNodes,
                                             sxyz<T1>& pt_i,
                                             const size_t & threadNo) const {
        std::set<std::array<T2,3>> faces;  // for some reason, unordered_set does not compile

        // loop over cells to find faces connected to current point
        for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
            // find nodes other that nodeNo
            std::array<T2,3> tmpnodes;
            size_t n = 0;
            for ( size_t nn=0; nn<4; ++nn ) {
                if ( nodeNo != tetrahedra[*nc].i[nn] ) {
                    tmpnodes[n++] = tetrahedra[*nc].i[nn];
                }
            }

            std::sort(tmpnodes.begin(), tmpnodes.end());
            faces.insert({nodeNo, tmpnodes[0], tmpnodes[1]});
            faces.insert({nodeNo, tmpnodes[0], tmpnodes[2]});
            faces.insert({nodeNo, tmpnodes[1], tmpnodes[2]});
        }

        bool found = false;
        sxyz<T1> g_old = g;
        T1 minAngle = std::numeric_limits<T1>::max();
        T1 ng = norm(g);
        // find projection that is closest to current gradient
        for ( auto fn=faces.begin(); fn!=faces.end(); ++fn ) {
            sxyz<T1> gtmp = projectOnFace(g_old, *fn);

            // find pt of intersection with opposing edge

            sxyz<T1> x1 = {nodes[(*fn)[1]].getX(),
                nodes[(*fn)[1]].getY(),
                nodes[(*fn)[1]].getZ()};
            sxyz<T1> x2 = {nodes[(*fn)[2]].getX(),
                nodes[(*fn)[2]].getY(),
                nodes[(*fn)[2]].getZ()};
            sxyz<T1> x4 = curr_pt + static_cast<T1>(100.0)*x1.getDistance(x2) * gtmp;

            sxyz<T1> a = x2 - x1;
            sxyz<T1> b = x4 - curr_pt;
            sxyz<T1> c = curr_pt - x1;

            sxyz<T1> ab = cross(a, b);
            T1 p = dot(cross(c, b), ab) / norm2(ab);

            if(p >= 0. && p <= 1.){// 0 <= p <= 1 is a sufficient condition to check if pt_i0 is between x1 and x2 and that the midpoint of the line segment [curr_p, pt_i0] is inside the face

                sxyz<T1> pt_i2 = p * x2 + static_cast<T1>(1.0 -p) * x1;

                c = pt_i2 - curr_pt;
                T1 dcg = dot(c, g_old);
                T1 angle = std::acos(dcg/(norm(c)*ng));

                if (  angle < minAngle ) {
                    g = pt_i2 - curr_pt;
                    pt_i = pt_i2;
                    edgeNodes[0] = (*fn)[1];
                    edgeNodes[1] = (*fn)[2];
                    minAngle = angle;
                    found = true;
                }

            }

        }
        if (!found || dot(g,g_old) < 0 ){
            T1 tmin = std::numeric_limits<T1>::max();
            T2 nodeTmin=0;
            for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
                for ( size_t nn=0; nn<4; ++nn ) {
                    T2 node_i = tetrahedra[*nc].i[nn];
                    if ( nodeNo != node_i && nodes[node_i].getTT(threadNo) < tmin ) {
                        nodeTmin = node_i;
                        tmin = nodes[node_i].getTT(threadNo);
                    }
                }
            }
            pt_i = nodes[nodeTmin];
            g = pt_i - curr_pt;
            edgeNodes[0] = nodeTmin;
            found = true;

        }
        return found;
    }
    template<typename T1, typename T2, typename NODE>
    sxyz<T1> Grid3Dun<T1,T2,NODE>::projectOnFace(const sxyz<T1>& curr_pt,
                                                 const sxyz<T1>& g,
                                                 std::array<T2,2>& edgeNodes,
                                                 const std::vector<T2> cells,
                                                 sxyz<T1>&  pt_i) const {

        // find all primary nodes of common cells
        std::vector<std::array<T2,4>> primary(cells.size());
        for (size_t nc=0; nc<cells.size(); ++nc) {
            size_t np = 0;
            for ( auto nn=this->neighbors[cells[nc]].begin(); nn!= this->neighbors[cells[nc]].end(); ++nn ) {
                if (nodes[*nn].isPrimary()) {
                    primary[nc][np++] = *nn;
                }
            }
        }
        // at the surface, the edge is common to two triangles
        // we must have only two other nodes at the surface, which are not common to the cells
        std::array<T2,2> edgeNodes2;
        size_t ne = 0;
        for (size_t nc=0; nc<cells.size(); ++nc) {
            for (size_t np=0; np<4; ++np) {
                if (primary[nc][np] == edgeNodes[0] || primary[nc][np] == edgeNodes[1]) {
                    continue;
                }
                bool found=false;
                for (size_t nc2=0; nc2<cells.size(); ++nc2) {
                    if ( nc == nc2 ) {
                        continue;
                    }
                    if ( std::find(primary[nc2].begin(), primary[nc2].end(), primary[nc][np])!=primary[nc2].end() ) {
                        found = true;  // primary node is common to two cells
                        break;
                    }
                }
                if ( found == false ) {
                    edgeNodes2[ne++] = primary[nc][np];
                }
            }
        }

        sxyz<T1> g_proj;
        T1 minAngle = std::numeric_limits<T1>::max();
        std::array<T2,2> oldEdgeNodes = edgeNodes;
        // find projection that is closest to current gradient
        for (auto en=edgeNodes2.begin(); en!=edgeNodes2.end(); ++en) {

            std::array<T2,2> tmpEdgeNodes = oldEdgeNodes;

            sxyz<T1> tmp = projectOnFace(g, {tmpEdgeNodes[0], tmpEdgeNodes[1], *en});
            T1 a2 = acos(dot(tmp, g)/(norm(tmp)*norm(g)));

#ifdef DEBUG_RP
            if ( verbose > 1 ) {
                std::cout << "\n\n\n\n\n# In projectOnFace\n\n";
                std::cout << "fig = plt.figure()\nax = fig.add_subplot(111, projection='3d')\n";
                std::cout << "cpt = np.array([" << curr_pt.x << ", " << curr_pt.y << ", " << curr_pt.z << "])\n";
                std::cout << "g = np.array([" << g.x << ", " << g.y << ", " << g.z << "])\n";
                std::cout << "pt2 = cpt + 0.05*g\n";
                std::cout << "c1 = np.array([" << nodes[tmpEdgeNodes[0]].getX() << ", " << nodes[tmpEdgeNodes[0]].getY() << ", " << nodes[tmpEdgeNodes[0]].getZ() << "])\n";
                std::cout << "c2 = np.array([" << nodes[tmpEdgeNodes[1]].getX() << ", " << nodes[tmpEdgeNodes[1]].getY() << ", " << nodes[tmpEdgeNodes[1]].getZ() << "])\n";
                std::cout << "g2 = np.array([" << tmp.x << ", " << tmp.y << ", " << tmp.z << "])\n";
                std::cout << "c3 = np.array([" << nodes[*en].getX() << ", " << nodes[*en].getY() << ", " << nodes[*en].getZ() << "])\n";
                std::cout << "a2 = " << a2 << "\n";
                std::cout << "pt3 = cpt + 0.05*g2\n";
                std::cout << "ax.plot([c1[0], c2[0], c3[0], c1[0]], [c1[1], c2[1], c3[1], c1[1]], [c1[2], c2[2], c3[2], c1[2]])\n";
                std::cout << "ax.plot([cpt[0], pt2[0]], [cpt[1], pt2[1]], [cpt[2], pt2[2]], c='g')\n";
                std::cout << "ax.plot([cpt[0], pt3[0]], [cpt[1], pt3[1]], [cpt[2], pt3[2]], c='k')\n";
                std::cout << "ax.scatter(cpt[0], cpt[1], cpt[2], c='r')\n";
                std::cout << "ax.scatter(c1[0], c1[1], c1[2], c='k')\n";
                std::cout << "ax.scatter(c2[0], c2[1], c2[2], c='k')\n";
                std::cout << "ax.scatter(pt2[0], pt2[1], pt2[2], c='g')\n";
            }
#endif
            if ( a2 < minAngle ) {
                // compute intersection point and check if within edge nodes

                sxyz<T1> x1 = {nodes[tmpEdgeNodes[0]].getX(),
                    nodes[tmpEdgeNodes[0]].getY(),
                    nodes[tmpEdgeNodes[0]].getZ()};
                sxyz<T1> x2 = {nodes[*en].getX(),
                    nodes[*en].getY(),
                    nodes[*en].getZ()};
                sxyz<T1> x4 = curr_pt + static_cast<T1>(10.0)*x1.getDistance(x2) * tmp;

                sxyz<T1> a = x2 - x1;
                sxyz<T1> b = x4 - curr_pt;
                sxyz<T1> c = curr_pt - x1;

                sxyz<T1> ab = cross(a, b);

                sxyz<T1> pt_i0 = x1 + (dot(cross(c, b), ab) / norm2(ab)) * a;

                sxyz<T1> mid_pt = static_cast<T1>(0.5) * ( curr_pt + pt_i0 );

                bool test1 = testInTriangle(&(nodes[ tmpEdgeNodes[0] ]),
                                            &(nodes[ tmpEdgeNodes[1] ]),
                                            &(nodes[ *en ]), mid_pt);

                // check if we are between x1 & x2

                b = pt_i0 - x1;
                T1 dab = dot(a, b);
                bool test2 = dab > 0.0 && norm2(b) <= norm2(a);

                // check if going in the same direction as g

                b = pt_i0 - curr_pt;
                T1 theta = acos(dot(b, g) / (norm(b)*norm(g)));
#ifdef DEBUG_RP
                if ( verbose > 1 ) {
                    std::cout << "pti = np.array([" << pt_i0.x << ", " << pt_i0.y << ", " << pt_i0.z << "])\n";
                    std::cout << "ax.scatter(pti[0], pti[1], pti[2], c='g')\n";
                    std::cout << "th = " << theta*180.0/3.14159 << '\n';
                }
#endif
                bool test3 = theta < theta_cut;  // defined in ttcr_t.h

                if ( test1 && test2 && test3 ) {
                    pt_i = pt_i0;
                    edgeNodes[0] = tmpEdgeNodes[0];
                    edgeNodes[1] = *en;
                    g_proj = tmp;
                    minAngle = a2;
#ifdef DEBUG_RP
                    if ( verbose > 1 ) {
                        std::cout << "\n\n\n\n\n# found\n\nfig = plt.figure()\nax = fig.add_subplot(111, projection='3d')\n";
                        std::cout << "pti = np.array([" << pt_i0.x << ", " << pt_i0.y << ", " << pt_i0.z << "])\n";
                        std::cout << "c1 = np.array([" << nodes[edgeNodes[0]].getX() << ", " << nodes[edgeNodes[0]].getY() << ", " << nodes[edgeNodes[0]].getZ() << "])\n";
                        std::cout << "c2 = np.array([" << nodes[edgeNodes[1]].getX() << ", " << nodes[edgeNodes[1]].getY() << ", " << nodes[edgeNodes[1]].getZ() << "])\n";
                        std::cout << "ax.plot([c1[0], c2[0]], [c1[1], c2[1]], [c1[2], c2[2]], 'k')\n";

                        std::cout << "ax.scatter(pti[0], pti[1], pti[2], c='r')\n";
                        std::cout << "ax.scatter(c1[0], c1[1], c1[2], c='k')\n";
                        std::cout << "ax.scatter(c2[0], c2[1], c2[2], c='k')\n";
                        std::cout << "plt.show()\n\n\n\n\n\n";
                    }
#endif
                    continue;
                }

                x1 = {nodes[tmpEdgeNodes[1]].getX(),
                    nodes[tmpEdgeNodes[1]].getY(),
                    nodes[tmpEdgeNodes[1]].getZ()};
                x2 = {nodes[*en].getX(),
                    nodes[*en].getY(),
                    nodes[*en].getZ()};

                a = x2 - x1;
                b = x4 - curr_pt;
                c = curr_pt - x1;

                ab = cross(a, b);

                pt_i0 = x1 + (dot(cross(c, b), ab) / norm2(ab)) * a;

                mid_pt = static_cast<T1>(0.5) * ( curr_pt + pt_i0 );

                test1 = testInTriangle(&(nodes[ tmpEdgeNodes[0] ]),
                                       &(nodes[ tmpEdgeNodes[1] ]),
                                       &(nodes[ *en ]), mid_pt);

                b = pt_i0 - x1;
                dab = dot(a, b);
                test2 = dab > 0.0 && dab <= norm2(a);

                b = pt_i0 - curr_pt;
                theta = acos(dot(b, g) / (norm(b)*norm(g)));
#ifdef DEBUG_RP
                if ( verbose > 1 ) {
                    std::cout << "pti = np.array([" << pt_i0.x << ", " << pt_i0.y << ", " << pt_i0.z << "])\n";
                    std::cout << "ax.scatter(pti[0], pti[1], pti[2], c='k')\n";
                    std::cout << "plt.show()\n\n";
                    std::cout << "th = " << theta*180.0/3.14159 << '\n';
                }
#endif
                test3 = theta < theta_cut;  // defined in ttcr_t.h

                if ( test1 && test2 && test3 ) {
                    pt_i = pt_i0;
                    edgeNodes[0] = *en;
                    edgeNodes[1] = tmpEdgeNodes[1];
                    g_proj = tmp;
                    minAngle = a2;
#ifdef DEBUG_RP
                    if ( verbose > 1 ) {
                        std::cout << "\n\n\n\n\n# found\n\nfig = plt.figure()\nax = fig.add_subplot(111, projection='3d')\n";
                        std::cout << "pti = np.array([" << pt_i0.x << ", " << pt_i0.y << ", " << pt_i0.z << "])\n";
                        std::cout << "c1 = np.array([" << nodes[edgeNodes[0]].getX() << ", " << nodes[edgeNodes[0]].getY() << ", " << nodes[edgeNodes[0]].getZ() << "])\n";
                        std::cout << "c2 = np.array([" << nodes[edgeNodes[1]].getX() << ", " << nodes[edgeNodes[1]].getY() << ", " << nodes[edgeNodes[1]].getZ() << "])\n";
                        std::cout << "ax.plot([c1[0], c2[0]], [c1[1], c2[1]], [c1[2], c2[2]], 'k')\n";
                        std::cout << "ax.scatter(pti[0], pti[1], pti[2], c='r')\n";
                        std::cout << "ax.scatter(c1[0], c1[1], c1[2], c='k')\n";
                        std::cout << "ax.scatter(c2[0], c2[1], c2[2], c='k')\n";
                        std::cout << "plt.show()\n\n\n\n\n\n";
                    }
#endif
                    continue;
                }
            }
        }
        return g_proj;
    }

    template<typename T1, typename T2, typename NODE>
    bool Grid3Dun<T1,T2,NODE>::projectOnFace(const sxyz<T1>& curr_pt,
                                             const T2 nodeNo,
                                             sxyz<T1>& g,
                                             std::array<T2,2>& edgeNodes,
                                             sxyz<T1>& pt_i) const {

        std::set<std::array<T2,3>> faces;  // for some reason, unordered_set does not compile

        // loop over cells to find faces connected to current point
        for ( auto nc=nodes[nodeNo].getOwners().begin(); nc!=nodes[nodeNo].getOwners().end(); ++nc ) {
            // find nodes other that nodeNo
            std::array<T2,3> tmpnodes;
            size_t n = 0;
            for ( size_t nn=0; nn<4; ++nn ) {
                if ( nodeNo != tetrahedra[*nc].i[nn] ) {
                    tmpnodes[n++] = tetrahedra[*nc].i[nn];
                }
            }

            std::sort(tmpnodes.begin(), tmpnodes.end());
            faces.insert({nodeNo, tmpnodes[0], tmpnodes[1]});
            faces.insert({nodeNo, tmpnodes[0], tmpnodes[2]});
            faces.insert({nodeNo, tmpnodes[1], tmpnodes[2]});
        }

        bool found = false;
        sxyz<T1> g_old = g;
        T1 minAngle = 9999.9;
        T1 ng = norm(g);
        // find projection that is closest to current gradient
        for ( auto fn=faces.begin(); fn!=faces.end(); ++fn ) {
            sxyz<T1> gtmp = projectOnFace(g_old, *fn);

            // find pt of intersection with opposing edge

            sxyz<T1> x1 = {nodes[(*fn)[1]].getX(),
                nodes[(*fn)[1]].getY(),
                nodes[(*fn)[1]].getZ()};
            sxyz<T1> x2 = {nodes[(*fn)[2]].getX(),
                nodes[(*fn)[2]].getY(),
                nodes[(*fn)[2]].getZ()};
            sxyz<T1> x4 = curr_pt + static_cast<T1>(100.0)*x1.getDistance(x2) * gtmp;

            sxyz<T1> a = x2 - x1;
            sxyz<T1> b = x4 - curr_pt;
            sxyz<T1> c = curr_pt - x1;

            sxyz<T1> ab = cross(a, b);

            sxyz<T1> pt_i2 = x1 + (dot(cross(c, b), ab) / norm2(ab)) * a;

            // check if pt_i is between x1 and x2
            b = pt_i2 - x1;
            T1 dab = dot(a, b);

            // check if going in the same direction as g
            c = pt_i2 - curr_pt;

            T1 dcg = dot(c, g_old);
            T1 angle = std::acos(dcg/(norm(c)*ng));

            if ( dab > 0.0 && norm2(b) <= norm2(a) && dcg > 0.0 && angle < minAngle ) {
                g = gtmp;
                pt_i = pt_i2;
                edgeNodes[0] = (*fn)[1];
                edgeNodes[1] = (*fn)[2];

                minAngle = angle;
                found = true;
            }
        }
        return found;
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::computeD(std::vector<sxyz<T1>> pts,
                                        std::vector<std::vector<sijv<T1>>> &d_data,
                                        const bool translated) const {

        if (this->translateOrigin == true && translated == false) {
            for ( size_t n=0; n<pts.size(); ++n ) {
                pts[n] -= this->origin;
            }
        }

        if ( d_data.size() != pts.size() ) {
            d_data.resize(pts.size());
        }
        for ( size_t i=0; i<pts.size(); ++i ) {
            d_data[i].resize(0);
        }

        for ( size_t np=0; np<pts.size(); ++np ) {
            bool found = false;
            for ( size_t nn=0; nn<nodes.size(); ++nn ) {
                if ( nodes[nn].getDistance(pts[np])<small ) {
                    found = true;
                    d_data[np].push_back( {np, nn, 1.0} );
                    break;
                }
            }
            if ( !found ) {
                T2 cellNO = this->getCellNo(pts[np]);
                if ( cellNO==std::numeric_limits<T2>::max() ) {
                    std::ostringstream msg;
                    msg << "Error: Point (" << pts[np] << ") outside grid.";
                    throw std::runtime_error(msg.str());
                }
                std::array<T1,4> weights;
                T1 sum (0.0);
                for ( size_t n=0; n<4; ++n ) {
                    weights[n] = 1.0/nodes[this->neighbors[cellNO][n]].getDistance(pts[np]);
                    sum += weights[n];
                }
                for ( size_t n=0; n<4; ++n ) {
                    weights[n] /= sum;
                    d_data[np].push_back( {np, this->neighbors[cellNO][n], weights[n]});
                }
            }
        }
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::getSurroundingNodes(const T2 nodeNumber,
                                                   const T2 minNbrPoints,
                                                   std::set<T2>& surroundingNodes) const {
        std::set<T2> layer;
        layer.insert(nodeNumber);
        T1 dx, dy, dz;
        int nzx, nzy, nzz;
        nzx = nzy = nzz = 0;

        while ( (surroundingNodes.size()+layer.size()-1) < minNbrPoints ) {
            std::copy(layer.begin(), layer.end(), std::inserter(surroundingNodes, surroundingNodes.end()));
            std::vector<T2> nextlayer;
            for ( auto nn=layer.begin(); nn!=layer.end(); ++nn ) {
                for ( auto cel=nodes[*nn].getOwners().begin(); cel!=nodes[*nn].getOwners().end(); cel++ ) {
                    for ( size_t i=0; i<4; ++i ) {
                        // first 4 nodes are the primary nodes
                        if ( surroundingNodes.find( this->neighbors[*cel][i] ) != surroundingNodes.end() )
                            continue;
                        dx = (nodes[nodeNumber].getX() - nodes[this->neighbors[*cel][i]].getX());
                        dy = (nodes[nodeNumber].getY() - nodes[this->neighbors[*cel][i]].getY());
                        dz = (nodes[nodeNumber].getZ() - nodes[this->neighbors[*cel][i]].getZ());
                        if ( dx == 0.0 )
                            nzx++;
                        if ( dy == 0.0 )
                            nzy++;
                        if ( dz == 0.0 )
                            nzz++;
                        if ( ( dx == 0.0 && nzx > 2 ) ||
                            ( dy == 0.0 && nzy > 2 ) ||
                            ( dz == 0.0 && nzz > 2 ) ) {
                            // allow only 2 nodes on X, Y or Z planes (typically external faces of domain)
                            //   as more lead to poorly conditionned system
                            continue;
                        }
                        nextlayer.push_back(this->neighbors[*cel][i]);

                    }
                }
            }
            if ( nextlayer.size() == 0 ) {
                throw std::runtime_error("Problem finding surrounding nodes");
            }
            layer.clear();
            std::copy( nextlayer.begin(), nextlayer.end(), std::inserter(layer,layer.end()));
        }
        std::copy(layer.begin(), layer.end(), std::inserter(surroundingNodes,surroundingNodes.end()));
        surroundingNodes.erase(nodeNumber);
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::buildA(const T2 nodeNumber,
                                      const std::set<T2>& surroundingNodes,
                                      const bool weighting,
                                      const int order,
                                      Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>& A,
                                      Eigen::Matrix<T1,Eigen::Dynamic, Eigen::Dynamic>& W) const {

        size_t npt = surroundingNodes.size();
        if ( order == 2 ) {
            A.resize(npt, 9);
        } else {
            A.resize(npt, 3);
        }

        if ( weighting ) {
            W.setZero(npt, npt);
            size_t i = 0;
            for ( auto nd=surroundingNodes.begin();nd!=surroundingNodes.end(); ++nd ) {
                A(i,0) = nodes[*nd].getX() - nodes[nodeNumber].getX();
                A(i,1) = nodes[*nd].getY() - nodes[nodeNumber].getY();
                A(i,2) = nodes[*nd].getZ() - nodes[nodeNumber].getZ();

                if ( order == 2 ) {
                    A(i,3) = 0.5*A(i,0)*A(i,0);
                    A(i,4) = 0.5*A(i,1)*A(i,1);
                    A(i,5) = 0.5*A(i,2)*A(i,2);

                    A(i,6) = A(i,0)*A(i,1);
                    A(i,7) = A(i,0)*A(i,2);
                    A(i,8) = A(i,1)*A(i,2);
                }

                W(i,i) = sqrt(1./(A(i,0)*A(i,0) + A(i,1)*A(i,1) + A(i,2)*A(i,2)));
                i++;
            }
            A = W*A;
        } else {
            size_t i = 0;
            for ( auto nd=surroundingNodes.begin(); nd!=surroundingNodes.end(); ++nd ) {
                A(i,0) = nodes[*nd].getX() - nodes[nodeNumber].getX();
                A(i,1) = nodes[*nd].getY() - nodes[nodeNumber].getY();
                A(i,2) = nodes[*nd].getZ() - nodes[nodeNumber].getZ();

                if ( order == 2 ) {
                    A(i,3) = 0.5*A(i,0)*A(i,0);
                    A(i,4) = 0.5*A(i,1)*A(i,1);
                    A(i,5) = 0.5*A(i,2)*A(i,2);

                    A(i,6) = A(i,0)*A(i,1);
                    A(i,7) = A(i,0)*A(i,2);
                    A(i,8) = A(i,1)*A(i,2);
                }

                i++;
            }
        }
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::buildA2(const T2 nodeNumber,
                                       const std::set<T2>& surroundingNodes,
                                       const bool weighting,
                                       const int order,
                                       Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>& A,
                                       Eigen::Matrix<T1,Eigen::Dynamic, Eigen::Dynamic>& W) const {

        size_t npt = surroundingNodes.size();
        if ( order == 2 ) {
            A.resize(npt, 10);
        } else {
            A.resize(npt, 4);
        }

        if ( weighting ) {
            W.setZero(npt, npt);
            size_t i = 0;
            for ( auto nd=surroundingNodes.begin();nd!=surroundingNodes.end(); ++nd ) {
                A(i,0) = nodes[*nd].getX() - nodes[nodeNumber].getX();
                A(i,1) = nodes[*nd].getY() - nodes[nodeNumber].getY();
                A(i,2) = nodes[*nd].getZ() - nodes[nodeNumber].getZ();

                if ( order == 1 ) {
                    A(i,3) = 1.0;
                } else {  // order == 2
                    A(i,3) = 0.5*A(i,0)*A(i,0);
                    A(i,4) = 0.5*A(i,1)*A(i,1);
                    A(i,5) = 0.5*A(i,2)*A(i,2);

                    A(i,6) = A(i,0)*A(i,1);
                    A(i,7) = A(i,0)*A(i,2);
                    A(i,8) = A(i,1)*A(i,2);

                    A(i,9) = 1.0;
                }

                W(i,i) = sqrt(1./(A(i,0)*A(i,0) + A(i,1)*A(i,1) + A(i,2)*A(i,2)));
                i++;
            }
            A = W*A;
        } else {
            size_t i = 0;
            for ( auto nd=surroundingNodes.begin(); nd!=surroundingNodes.end(); ++nd ) {
                A(i,0) = nodes[*nd].getX() - nodes[nodeNumber].getX();
                A(i,1) = nodes[*nd].getY() - nodes[nodeNumber].getY();
                A(i,2) = nodes[*nd].getZ() - nodes[nodeNumber].getZ();

                if ( order == 1 ) {
                    A(i,3) = 1.0;
                } else {  // order == 2
                    A(i,3) = 0.5*A(i,0)*A(i,0);
                    A(i,4) = 0.5*A(i,1)*A(i,1);
                    A(i,5) = 0.5*A(i,2)*A(i,2);

                    A(i,6) = A(i,0)*A(i,1);
                    A(i,7) = A(i,0)*A(i,2);
                    A(i,8) = A(i,1)*A(i,2);

                    A(i,9) = 1.0;
                }

                i++;
            }
        }
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::fill_k_data(const T2 nodeNo,
                                           const std::set<T2>& surroundingNodes,
                                           const int i,
                                           const int j,
                                           const int k,
                                           const Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>& Acoefs,
                                           std::vector<std::vector<std::vector<siv<T1>>>>& k_data) const {
        k_data[0][nodeNo].resize(0);
        k_data[1][nodeNo].resize(0);
        k_data[2][nodeNo].resize(0);
        std::array<T1,3> sum{0., 0., 0.};
        siv<T1> coef;

        size_t c = 0;
        for ( auto nd=surroundingNodes.begin(); nd!=surroundingNodes.end(); ++nd ) {
            coef.i = *nd;
            coef.v = Acoefs(i, c);
            sum[0] += Acoefs(i, c);
            k_data[0][nodeNo].push_back(coef);
            coef.v = Acoefs(j, c);
            sum[1] += Acoefs(j, c);
            k_data[1][nodeNo].push_back(coef);
            coef.v = Acoefs(k, c);
            sum[2] += Acoefs(k, c);
            k_data[2][nodeNo].push_back(coef);
            c++;
        }
        coef.i = nodeNo;
        coef.v = -sum[0];
        k_data[0][nodeNo].push_back(coef);
        coef.v = -sum[1];
        k_data[1][nodeNo].push_back(coef);
        coef.v = -sum[2];
        k_data[2][nodeNo].push_back(coef);
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::fill_k_data2(const T2 nodeNo,
                                            const std::set<T2>& surroundingNodes,
                                            const int i,
                                            const int j,
                                            const int k,
                                            const Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>& Acoefs,
                                            std::vector<std::vector<std::vector<siv<T1>>>>& k_data) const {
        k_data[0][nodeNo].resize(0);
        k_data[1][nodeNo].resize(0);
        k_data[2][nodeNo].resize(0);
        siv<T1> coef;

        size_t c = 0;
        for ( auto nd=surroundingNodes.begin(); nd!=surroundingNodes.end(); ++nd ) {
            coef.i = *nd;
            coef.v = Acoefs(i, c);
            k_data[0][nodeNo].push_back(coef);
            coef.v = Acoefs(j, c);
            k_data[1][nodeNo].push_back(coef);
            coef.v = Acoefs(k, c);
            k_data[2][nodeNo].push_back(coef);
            c++;
        }
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::computeK(std::vector<std::vector<std::vector<siv<T1>>>>& k_data,
                                        const int order,
                                        const int taylorSeriesOrder,
                                        const bool weighting,
                                        const bool s0inside,
                                        const int additionnalPoints) const {
        if ( order!=1 && order!=2 ) {
            throw std::runtime_error("order in computeK should be 1 or 2");
        }
        if ( taylorSeriesOrder!=1 && taylorSeriesOrder!=2 ) {
            throw std::runtime_error("taylorSeriesOrder in computeK should be 1 or 2");
        }
        if ( order == 2 && taylorSeriesOrder == 1 ) {
            throw std::runtime_error("2nd order derivative operator requires 2nd order Taylor series expansion");
        }

        if ( k_data.size() != 3 ) {
            k_data.resize(3);
        }

        k_data[0].resize(nPrimary);
        k_data[1].resize(nPrimary);
        k_data[2].resize(nPrimary);

        T2 minNbrPoints = 4;
        if ( taylorSeriesOrder == 2 ) {
            minNbrPoints = 10 + additionnalPoints;
        }

        int neededRank = 9;
        if ( taylorSeriesOrder == 1 ) {
            neededRank = 3;
        }
        if ( s0inside ) {
            neededRank += 1;
        }

        for ( T2 n=0; n<nPrimary; ++n ) {

            // find surrounding nodes
            std::set<T2> surroundingNodes;
            getSurroundingNodes(n, minNbrPoints, surroundingNodes);

            Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> Acoefs;
            Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> A;
            Eigen::Matrix<T1,Eigen::Dynamic, Eigen::Dynamic> W;

            if ( s0inside ) {
                buildA2(n, surroundingNodes, weighting, taylorSeriesOrder, A, W);
            } else {
                buildA(n, surroundingNodes, weighting, taylorSeriesOrder, A, W);
            }

            Eigen::Index rank = pseudoInverse(A, Acoefs);

            if ( rank < neededRank ) {
                surroundingNodes.clear();
                getSurroundingNodes(n, 2*minNbrPoints, surroundingNodes);
                if ( s0inside ) {
                    buildA2(n, surroundingNodes, weighting, taylorSeriesOrder, A, W);
                } else {
                    buildA(n, surroundingNodes, weighting, taylorSeriesOrder, A, W);
                }
                rank = pseudoInverse(A, Acoefs);
                if ( rank < neededRank ) {
                    throw std::runtime_error("Mesh appears poorly conditionned, unable to compute matrix K");
                }
            }

            if ( weighting ) {
                Acoefs *= W;
            }

            if ( order == 1 ) {
                if ( s0inside ) {
                    fill_k_data2(n, surroundingNodes, 0, 1, 2, Acoefs, k_data);
                } else {
                    fill_k_data(n, surroundingNodes, 0, 1, 2, Acoefs, k_data);
                }
            } else {  // order == 2
                if ( s0inside ) {
                    fill_k_data2(n, surroundingNodes, 3, 4, 5, Acoefs, k_data);
                } else {
                    fill_k_data(n, surroundingNodes, 3, 4, 5, Acoefs, k_data);
                }
            }

        }  // end of for nnodes
    }

    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::interpVelocitySecondary(const T2 nSecondary) {

        T2 nNodes = nPrimary;

        std::map<std::array<T2,2>,std::vector<T2>> lineMap;
        std::array<T2,2> lineKey;
        typename std::map<std::array<T2,2>,std::vector<T2>>::iterator lineIt;

        size_t nFaceNodes = 0;
        for ( size_t n=1; n<=(nSecondary-1); ++n ) {
            nFaceNodes += n;
        }

        for ( T2 ntet=0; ntet<tetrahedra.size(); ++ntet ) {

            // for each triangle
            for ( T2 ntri=0; ntri<4; ++ntri ) {

                // start from ntri to avoid redundancy
                for ( size_t nl=ntri; nl<3; ++nl ) {

                    lineKey = {tetrahedra[ntet].i[ iNodes[ntri][nl] ],
                        tetrahedra[ntet].i[ iNodes[ntri][(nl+1)%3] ]};
                    std::sort(lineKey.begin(), lineKey.end());

                    lineIt = lineMap.find( lineKey );
                    if ( lineIt == lineMap.end() ) {
                        // not found, insert new pair
                        lineMap[ lineKey ] = std::vector<T2>(nSecondary);
                    } else {
                        continue;
                    }

                    T1 slope = (1.0/nodes[lineKey[1]].getNodeSlowness() - 1.0/nodes[lineKey[0]].getNodeSlowness())/
                    nodes[lineKey[1]].getDistance(nodes[lineKey[0]]);

                    for ( size_t n2=0; n2<nSecondary; ++n2 ) {
                        T1 s = 1.0/(1.0/nodes[lineKey[0]].getNodeSlowness() + slope * nodes[nNodes].getDistance(nodes[lineKey[0]]));
                        nodes[nNodes].setNodeSlowness( s );
                        lineMap[lineKey][n2] = nNodes++;
                    }
                }
            }
        }


        if ( nSecondary > 1 ) {

            std::map<std::array<T2,3>,std::vector<T2>> faceMap;
            std::array<T2,3> faceKey;
            typename std::map<std::array<T2,3>,std::vector<T2>>::iterator faceIt;

            ptrdiff_t ncut = nSecondary - 1;

            for ( T2 ntet=0; ntet<tetrahedra.size(); ++ntet ) {

                // for each triangle
                for ( T2 ntri=0; ntri<4; ++ntri ) {

                    faceKey = {tetrahedra[ntet].i[ iNodes[ntri][0] ],
                        tetrahedra[ntet].i[ iNodes[ntri][1] ],
                        tetrahedra[ntet].i[ iNodes[ntri][2] ]};
                    std::sort(faceKey.begin(), faceKey.end());


                    faceIt = faceMap.find( faceKey );
                    if ( faceIt == faceMap.end() ) {
                        // not found, insert new pair
                        faceMap[ faceKey ] = std::vector<T2>(nFaceNodes);
                    } else {
                        continue;
                    }

                    std::vector<NODE*> inodes;
                    inodes.push_back( &(nodes[faceKey[0]]) );
                    inodes.push_back( &(nodes[faceKey[1]]) );
                    inodes.push_back( &(nodes[faceKey[2]]) );

                    size_t ifn = 0;
                    for ( ptrdiff_t n=0; n<ncut; ++n ) {
                        size_t nseg = ncut+1-n;
                        for ( size_t n2=0; n2<nseg-1; ++n2 ) {

                            T1 s = Interpolator<T1>::bilinearTriangleVel(nodes[nNodes], inodes);
                            nodes[nNodes].setNodeSlowness( s );

                            faceMap[faceKey][ifn++] = nNodes++;

                        }
                    }
                }
            }
        }
    }


    template<typename T1, typename T2, typename NODE>
    void Grid3Dun<T1,T2,NODE>::interpSlownessSecondary(const T2 nSecondary) {

        T2 nNodes = this->nPrimary;

        std::map<std::array<T2,2>,std::vector<T2>> lineMap;
        std::array<T2,2> lineKey;
        typename std::map<std::array<T2,2>,std::vector<T2>>::iterator lineIt;

        size_t nFaceNodes = 0;
        for ( size_t n=1; n<=(nSecondary-1); ++n ) {
            nFaceNodes += n;
        }

        for ( T2 ntet=0; ntet<tetrahedra.size(); ++ntet ) {

            // for each triangle
            for ( T2 ntri=0; ntri<4; ++ntri ) {

                // start from ntri to avoid redundancy
                for ( size_t nl=ntri; nl<3; ++nl ) {

                    lineKey = {tetrahedra[ntet].i[ iNodes[ntri][nl] ],
                        tetrahedra[ntet].i[ iNodes[ntri][(nl+1)%3] ]};
                    std::sort(lineKey.begin(), lineKey.end());

                    lineIt = lineMap.find( lineKey );
                    if ( lineIt == lineMap.end() ) {
                        // not found, insert new pair
                        lineMap[ lineKey ] = std::vector<T2>(nSecondary);
                    } else {
                        continue;
                    }

                    T1 slope = (nodes[lineKey[1]].getNodeSlowness() - nodes[lineKey[0]].getNodeSlowness())/
                    nodes[lineKey[1]].getDistance(nodes[lineKey[0]]);

                    for ( size_t n2=0; n2<nSecondary; ++n2 ) {
                        T1 s = nodes[lineKey[0]].getNodeSlowness() + slope * nodes[nNodes].getDistance(nodes[lineKey[0]]);
                        nodes[nNodes].setNodeSlowness( s );
                        lineMap[lineKey][n2] = nNodes++;
                    }
                }
            }
        }


        if ( nSecondary > 1 ) {

            std::map<std::array<T2,3>,std::vector<T2>> faceMap;
            std::array<T2,3> faceKey;
            typename std::map<std::array<T2,3>,std::vector<T2>>::iterator faceIt;

            ptrdiff_t ncut = nSecondary - 1;

            for ( T2 ntet=0; ntet<tetrahedra.size(); ++ntet ) {

                // for each triangle
                for ( T2 ntri=0; ntri<4; ++ntri ) {

                    faceKey = {tetrahedra[ntet].i[ iNodes[ntri][0] ],
                        tetrahedra[ntet].i[ iNodes[ntri][1] ],
                        tetrahedra[ntet].i[ iNodes[ntri][2] ]};
                    std::sort(faceKey.begin(), faceKey.end());


                    faceIt = faceMap.find( faceKey );
                    if ( faceIt == faceMap.end() ) {
                        // not found, insert new pair
                        faceMap[ faceKey ] = std::vector<T2>(nFaceNodes);
                    } else {
                        continue;
                    }

                    std::vector<NODE*> inodes;
                    inodes.push_back( &(nodes[faceKey[0]]) );
                    inodes.push_back( &(nodes[faceKey[1]]) );
                    inodes.push_back( &(nodes[faceKey[2]]) );

                    size_t ifn = 0;
                    for ( ptrdiff_t n=0; n<ncut; ++n ) {
                        size_t nseg = ncut+1-n;
                        for ( size_t n2=0; n2<nseg-1; ++n2 ) {

                            T1 s = Interpolator<T1>::bilinearTriangle(nodes[nNodes], inodes);
                            nodes[nNodes].setNodeSlowness( s );

                            faceMap[faceKey][ifn++] = nNodes++;

                        }
                    }
                }
            }
        }
    }

    template<typename T1, typename T2, typename NODE>
    const T1 Grid3Dun<T1,T2,NODE>::getAverageEdgeLength() const {
        std::set<std::array<T2,2>> edges;
        typename std::set<std::array<T2,2>>::iterator edgIt;
        T2 iNodes[6][2] = {
            {0,1},
            {0,2},
            {0,3},
            {1,2},
            {1,3},
            {2,3}
        };
        T1 sum = 0.0;
        for (size_t ntet=0; ntet<tetrahedra.size(); ++ntet) {
            for (size_t n=0; n<6; ++n) {
                std::array<T2, 2> edgei = {tetrahedra[ntet].i[iNodes[n][0]],
                    tetrahedra[ntet].i[iNodes[n][1]]};
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

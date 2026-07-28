//
//  NodeKDTree2D.h
//  ttcr
//

/**
 * @file NodeKDTree2D.h
 * @brief Static 2-D (x-z plane) nearest-neighbour index over the nodes of an
 *        unstructured or rectilinear grid.
 *
 * The two-dimensional counterpart of NodeKDTree.h: a nanoflann k-d tree used to
 * accelerate point location and on-node lookups. Only the x and z coordinates
 * are indexed, matching the x-z convention of the 2-D grids.
 *
 * The node coordinates are **copied** at construction, so the index is
 * self-contained and independent of any later growth of the node vector — such
 * as secondary nodes being appended after the grid is built. The flip side is
 * that it is a snapshot: if node coordinates change, the index must be rebuilt.
 *
 * Build it once; queries are read-only and hence safe to call concurrently from
 * several threads.
 *
 * @sa NodeKDTree.h, Grid2Dun.h, Grid2Duc.h
 */

#ifndef ttcr_NodeKDTree2D_h
#define ttcr_NodeKDTree2D_h

#include <array>
#include <cstddef>
#include <vector>

#include "nanoflann.hpp"

namespace ttcr {

    /**
     * @brief k-d tree over the nodes of a 2-D grid, in the x-z plane.
     *
     * @tparam T1 floating-point type of the node coordinates.
     * @tparam T2 integer type used for node indices, normally @c uint32_t.
     *
     * Non-copyable: the nanoflann index holds a reference to the point cloud
     * member, so a copy would leave the new index pointing at the original's
     * cloud. Hold it by pointer or construct it in place.
     */
    template<typename T1, typename T2>
    class NodeKDTree2D {
    public:
        /**
         * @brief Build the index from the first @p nNodes entries of @p nodes.
         * @tparam NODE node type exposing @c getX() and @c getZ().
         * @param[in] nodes   node vector; only the leading @p nNodes are indexed.
         * @param[in] nNodes  number of nodes to index. Unlike the 3-D version this
         *                    is not restricted to primary nodes — pass whichever
         *                    leading range the caller wants searchable.
         * @pre @p nodes holds at least @p nNodes entries.
         */
        template<typename NODE>
        NodeKDTree2D(const std::vector<NODE>& nodes, const size_t nNodes)
        : cloud(makeCloud(nodes, nNodes)),
          index(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10)) {}

        NodeKDTree2D(const NodeKDTree2D&) = delete;             ///< Non-copyable; see the class note.
        NodeKDTree2D& operator=(const NodeKDTree2D&) = delete;  ///< Non-copy-assignable; see the class note.

        /**
         * @brief Find the node closest to a point in the x-z plane.
         * @param[in] x query x coordinate.
         * @param[in] z query z coordinate.
         * @return Grid index of the nearest node.
         * @note Const and free of shared mutable state, hence callable
         *       concurrently from several threads.
         */
        T2 findNearest(const T1 x, const T1 z) const {
            const T1 query[2] = { x, z };
            T2 idx = 0;
            T1 d2 = 0;
            index.knnSearch(query, 1, &idx, &d2);
            return idx;
        }

    private:
        /// Node coordinates in the layout nanoflann expects, satisfying its
        /// dataset adaptor requirements.
        struct PointCloud {
            std::vector<std::array<T1,2>> pts;  ///< One (x,z) pair per indexed node.
            /// @return Number of indexed points (nanoflann adaptor hook).
            inline size_t kdtree_get_point_count() const { return pts.size(); }
            /// @param idx point number. @param dim 0 for x, 1 for z. @return That coordinate (nanoflann adaptor hook).
            inline T1 kdtree_get_pt(const size_t idx, const size_t dim) const {
                return pts[idx][dim];
            }
            /// @return False, declining to supply a precomputed bounding box (nanoflann adaptor hook).
            template<class BBOX>
            bool kdtree_get_bbox(BBOX&) const { return false; }
        };

        /// Concrete nanoflann index: single-index adaptor, squared-L2 metric, 2 dimensions.
        using index_t = nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L2_Simple_Adaptor<T1, PointCloud>, PointCloud, 2, T2>;

        /**
         * @brief Copy the x and z node coordinates into a @ref PointCloud.
         * @tparam NODE node type exposing @c getX() and @c getZ().
         * @param[in] nodes   source node vector.
         * @param[in] nNodes  number of leading nodes to copy.
         * @return Populated point cloud, moved into the @ref cloud member.
         */
        template<typename NODE>
        static PointCloud makeCloud(const std::vector<NODE>& nodes,
                                    const size_t nNodes) {
            PointCloud c;
            c.pts.resize(nNodes);
            for ( size_t n=0; n<nNodes; ++n ) {
                c.pts[n][0] = nodes[n].getX();
                c.pts[n][1] = nodes[n].getZ();
            }
            return c;
        }

        PointCloud cloud;  ///< Owned copy of the coordinates; declared before @ref index, which references it.
        index_t index;     ///< nanoflann index built over @ref cloud.
    };
}

#endif

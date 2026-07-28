//
//  NodeKDTree.h
//  ttcr
//

/**
 * @file NodeKDTree.h
 * @brief Static 3-D nearest-neighbour index over the primary nodes of an
 *        unstructured mesh.
 *
 * Wraps a nanoflann k-d tree so that point location (@c getCellNo and friends)
 * can start from the closest primary node instead of scanning every cell,
 * turning a linear search into a logarithmic one.
 *
 * The node coordinates are **copied** at construction, so the index is
 * self-contained and independent of any later growth of the node vector — such
 * as secondary nodes being appended after the mesh is built. The flip side is
 * that it is a snapshot: if node coordinates change, the index must be rebuilt.
 *
 * Build it once; queries are read-only and hence safe to call concurrently from
 * several threads.
 *
 * @sa NodeKDTree2D.h, Grid3Dun.h, Grid3Duc.h
 */

#ifndef ttcr_NodeKDTree_h
#define ttcr_NodeKDTree_h

#include <array>
#include <cstddef>
#include <vector>

#include "nanoflann.hpp"

namespace ttcr {

    /**
     * @brief k-d tree over the primary nodes of a 3-D mesh.
     *
     * @tparam T1 floating-point type of the node coordinates.
     * @tparam T2 integer type used for node indices, normally @c uint32_t.
     *
     * Non-copyable: the nanoflann index holds a reference to the point cloud
     * member, so a copy would leave the new index pointing at the original's
     * cloud. Hold it by pointer or construct it in place.
     */
    template<typename T1, typename T2>
    class NodeKDTree {
    public:
        /**
         * @brief Build the index from the first @p nPrimary entries of @p nodes.
         * @tparam NODE node type exposing @c getX(), @c getY() and @c getZ().
         * @param[in] nodes     node vector; only the leading primary nodes are indexed.
         * @param[in] nPrimary  number of primary nodes, i.e. the mesh vertices,
         *                      excluding any secondary nodes stored after them.
         * @pre @p nodes holds at least @p nPrimary entries.
         */
        template<typename NODE>
        NodeKDTree(const std::vector<NODE>& nodes, const size_t nPrimary)
        : cloud(makeCloud(nodes, nPrimary)),
          index(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10)) {}

        NodeKDTree(const NodeKDTree&) = delete;             ///< Non-copyable; see the class note.
        NodeKDTree& operator=(const NodeKDTree&) = delete;  ///< Non-copy-assignable; see the class note.

        /**
         * @brief Find the primary node closest to a point.
         * @param[in] x query x coordinate.
         * @param[in] y query y coordinate.
         * @param[in] z query z coordinate.
         * @return Grid index of the nearest primary node.
         * @note Const and free of shared mutable state, hence callable
         *       concurrently from several threads.
         */
        T2 findNearest(const T1 x, const T1 y, const T1 z) const {
            const T1 query[3] = { x, y, z };
            T2 idx = 0;
            T1 d2 = 0;
            index.knnSearch(query, 1, &idx, &d2);
            return idx;
        }

    private:
        /// Node coordinates in the layout nanoflann expects, satisfying its
        /// dataset adaptor requirements.
        struct PointCloud {
            std::vector<std::array<T1,3>> pts;  ///< One (x,y,z) triple per indexed node.
            /// @return Number of indexed points (nanoflann adaptor hook).
            inline size_t kdtree_get_point_count() const { return pts.size(); }
            /// @param idx point number. @param dim 0, 1 or 2. @return That coordinate (nanoflann adaptor hook).
            inline T1 kdtree_get_pt(const size_t idx, const size_t dim) const {
                return pts[idx][dim];
            }
            /// @return False, declining to supply a precomputed bounding box (nanoflann adaptor hook).
            template<class BBOX>
            bool kdtree_get_bbox(BBOX&) const { return false; }
        };

        /// Concrete nanoflann index: single-index adaptor, squared-L2 metric, 3 dimensions.
        using index_t = nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L2_Simple_Adaptor<T1, PointCloud>, PointCloud, 3, T2>;

        /**
         * @brief Copy the primary node coordinates into a @ref PointCloud.
         * @tparam NODE node type exposing @c getX(), @c getY() and @c getZ().
         * @param[in] nodes     source node vector.
         * @param[in] nPrimary  number of leading nodes to copy.
         * @return Populated point cloud, moved into the @ref cloud member.
         */
        template<typename NODE>
        static PointCloud makeCloud(const std::vector<NODE>& nodes,
                                    const size_t nPrimary) {
            PointCloud c;
            c.pts.resize(nPrimary);
            for ( size_t n=0; n<nPrimary; ++n ) {
                c.pts[n][0] = nodes[n].getX();
                c.pts[n][1] = nodes[n].getY();
                c.pts[n][2] = nodes[n].getZ();
            }
            return c;
        }

        PointCloud cloud;  ///< Owned copy of the coordinates; declared before @ref index, which references it.
        index_t index;     ///< nanoflann index built over @ref cloud.
    };
}

#endif

//
//  NodeKDTree.h
//  ttcr
//
//  Static 3-D nearest-neighbour index over the primary nodes of an
//  unstructured mesh, used to accelerate point location (getCellNo).
//
//  The node coordinates are copied at construction, so the index is
//  self-contained and independent of any later growth of the node vector
//  (e.g. appended secondary nodes).  Build it once; queries are read-only and
//  hence safe to call concurrently from several threads.
//

#ifndef ttcr_NodeKDTree_h
#define ttcr_NodeKDTree_h

#include <array>
#include <cstddef>
#include <vector>

#include "nanoflann.hpp"

namespace ttcr {

    template<typename T1, typename T2>
    class NodeKDTree {
    public:
        template<typename NODE>
        NodeKDTree(const std::vector<NODE>& nodes, const size_t nPrimary)
        : cloud(makeCloud(nodes, nPrimary)),
          index(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10)) {}

        NodeKDTree(const NodeKDTree&) = delete;
        NodeKDTree& operator=(const NodeKDTree&) = delete;

        // Grid index of the primary node closest to (x, y, z).
        T2 findNearest(const T1 x, const T1 y, const T1 z) const {
            const T1 query[3] = { x, y, z };
            T2 idx = 0;
            T1 d2 = 0;
            index.knnSearch(query, 1, &idx, &d2);
            return idx;
        }

    private:
        struct PointCloud {
            std::vector<std::array<T1,3>> pts;
            inline size_t kdtree_get_point_count() const { return pts.size(); }
            inline T1 kdtree_get_pt(const size_t idx, const size_t dim) const {
                return pts[idx][dim];
            }
            template<class BBOX>
            bool kdtree_get_bbox(BBOX&) const { return false; }
        };

        using index_t = nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L2_Simple_Adaptor<T1, PointCloud>, PointCloud, 3, T2>;

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

        PointCloud cloud;
        index_t index;
    };
}

#endif

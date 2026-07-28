//
//  Node2Dnsp.h
//  ttcr
//
//  Created by Bernard Giroux on 2012-09-25.
//  Copyright (c) 2012 Bernard Giroux. All rights reserved.
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
 * @file Node2Dnsp.h
 * @brief Node of a 2-D shortest-path grid with node-based slowness.
 *
 * Declares ttcr::Node2Dnsp, the node type of Grid2Drnsp. It is to
 * ttcr::Node2Dcsp what ttcr::Node2Dn is to ttcr::Node2Dc: the same
 * shortest-path bookkeeping, plus a slowness value carried on the node.
 *
 * See @ref n2dcsp_parents for what the parent arrays are for and
 * @ref n2dc_threads for the per-thread convention. Compared with
 * ttcr::Node2Dcsp this class adds @ref ttcr::Node2Dnsp::getNodeSlowness, drops
 * @c getSize(), and offers fewer constructors.
 *
 * @note Uses @ref ttcr::sxz and @ref ttcr::small but includes only Node.h,
 *       relying on includers to have pulled in ttcr_t.h first — as Node2Dn.h
 *       does, and unlike Node2Dcsp.h which includes it explicitly.
 *
 * @sa Node2Dn.h, Node2Dcsp.h, Grid2Drnsp.h
 */

#ifndef ttcr_Node2Dnsp_h
#define ttcr_Node2Dnsp_h

#include <cmath>
#include <limits>

#include "Node.h"

namespace ttcr {

    /**
     * @brief Node of a 2-D shortest-path grid, slowness held at the node.
     *
     * Lives in the x-z plane; @ref getY always returns 0 so the class still
     * satisfies the 3-D ttcr::Node interface.
     *
     * @tparam T1 floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2 integer type of grid, cell and parent indices, normally
     *            @c uint32_t.
     *
     * @warning Owns three raw @c new[] arrays and declares a copy constructor
     *          and a destructor but **no copy-assignment operator** — see the
     *          same warning on ttcr::Node2Dcsp. Copy-construct, never assign.
     */
    template<typename T1, typename T2>
    class Node2Dnsp : public Node<T1> {
    public:
        /**
         * @brief Construct an unpositioned node with everything unset.
         * @param nt number of threads, i.e. the size of each of the three arrays.
         * @post Coordinates and slowness are 0; traveltimes and both parent
         *       indices hold their type's maximum, the "not yet reached" marker.
         */
        Node2Dnsp(const size_t nt=1) :
        nThreads(nt),
        tt(0),
        x(0.0), z(0.0), slowness(0.0),
        gridIndex(std::numeric_limits<T2>::max()),
        nodeParent(0),
        cellParent(0),
        owners(0),
        primary(false)
        {
            tt = new T1[nt];
            nodeParent = new T2[nt];
            cellParent = new T2[nt];

            for ( size_t n=0; n<nt; ++n ) {
                tt[n] = std::numeric_limits<T1>::max();
                nodeParent[n] = std::numeric_limits<T2>::max();
                cellParent[n] = std::numeric_limits<T2>::max();
            }
        }


        /**
         * @brief Construct from a point, with a known traveltime on one thread.
         *
         * Used to seed a source position; the parent indices stay unset, which
         * is what marks the node as a raypath origin.
         *
         * @param t  traveltime to store.
         * @param s  point supplying the coordinates.
         * @param nt number of threads.
         * @param i  thread number the traveltime belongs to.
         * @pre @p i is less than @p nt; the write is unchecked.
         * @note Slowness is left at 0 and must be set separately with
         *       @ref setNodeSlowness.
         */
        Node2Dnsp(const T1 t, const sxz<T1>& s, const size_t nt, const size_t i) :
        nThreads(nt),
        tt(0),
        x(s.x), z(s.z), slowness(0.0),
        gridIndex(std::numeric_limits<T2>::max()),
        nodeParent(0),
        cellParent(0),
        owners(std::vector<T2>(0)),
        primary(false)
        {
            tt = new T1[nt];
            nodeParent = new T2[nt];
            cellParent = new T2[nt];

            for ( size_t n=0; n<nt; ++n ) {
                tt[n] = std::numeric_limits<T1>::max();
                nodeParent[n] = std::numeric_limits<T2>::max();
                cellParent[n] = std::numeric_limits<T2>::max();
            }
            tt[i] = t;
        }

        /**
         * @brief Copy constructor; deep-copies all three per-thread arrays.
         * @param node node to copy.
         */
        Node2Dnsp(const Node2Dnsp<T1,T2>& node) :
        nThreads(node.nThreads),
        tt(0),
        x(node.x), z(node.z), slowness(node.slowness),
        gridIndex(node.gridIndex),
        nodeParent(0),
        cellParent(0),
        owners(node.owners),
        primary(node.primary)
        {
            tt = new T1[nThreads];
            nodeParent = new T2[nThreads];
            cellParent = new T2[nThreads];

            for ( size_t n=0; n<nThreads; ++n ) {
                tt[n] = node.tt[n];
                nodeParent[n] = node.nodeParent[n];
                cellParent[n] = node.cellParent[n];
            }
        }


        /// Destructor; frees all three per-thread arrays.
        virtual ~Node2Dnsp() {
            delete [] tt;
            delete [] nodeParent;
            delete [] cellParent;
        }

        /**
         * @brief Reset one thread's traveltime and both parent indices.
         * @param thread_no thread whose slots are reset.
         * @note Resets all three arrays together, so a stale raypath cannot
         *       survive into the next source's solve.
         */
        void reinit(const size_t thread_no) { //=0) {
            tt[thread_no] = std::numeric_limits<T1>::max();
            nodeParent[thread_no] = std::numeric_limits<T2>::max();
            cellParent[thread_no] = std::numeric_limits<T2>::max();
        }

        /// @param i thread number. @return Traveltime stored for that thread.
        T1 getTT(const size_t i) const { return tt[i]; }
        /// @param t traveltime. @param i thread number. Stores @p t for that thread.
        void setTT(const T1 t, const size_t i) { tt[i] = t; }

        /**
         * @brief Set position and grid index in one call.
         * @param xx    x coordinate.
         * @param zz    z coordinate.
         * @param index index of this node in the grid's node list.
         */
        void setXZindex(const T1 xx, const T1 zz, const T2 index) {
            x=xx; z=zz; gridIndex = index;  }

        /**
         * @brief Set position from a point object, plus the grid index.
         * @tparam SXZ any type with @c x and @c z members.
         * @param s     point supplying the coordinates.
         * @param index index of this node in the grid's node list.
         */
        template<typename SXZ>
        void setXYZindex(const SXZ& s, const T2 index) {
            x=s.x; z=s.z; gridIndex = index;  }

        T1 getX() const {              ///< @return x coordinate.
            return x;
        }
        void setX(const T1 xx) { x = xx; }  ///< @param xx new x coordinate.

        /// @return Always 0 — the node is planar; present to satisfy ttcr::Node.
        T1 getY() const { return 0.0; }

        T1 getZ() const { return z; }       ///< @return z coordinate.
        void setZ(const T1 zz) { z = zz; }  ///< @param zz new z coordinate.

        /**
         * @brief Slowness at this node.
         * @return Slowness, i.e. reciprocal velocity.
         * @note This accessor is what distinguishes ttcr::Node2Dnsp from
         *       ttcr::Node2Dcsp.
         */
        T1 getNodeSlowness() const { return slowness; }
        /**
         * @brief Set the slowness at this node.
         * @param s slowness value.
         * @warning Unvalidated: a zero or negative slowness is stored as given.
         */
        void setNodeSlowness(const T1 s) { slowness = s; }

        /// @return Index of this node in the grid's node list.
        T2 getGridIndex() const { return gridIndex; }
        /// @param index new index of this node in the grid's node list.
        void setGridIndex(const T2 index) { gridIndex = index; }

        /**
         * @brief Node the ray arrived from, on a given thread.
         * @param i thread number.
         * @return Index of the parent node, or @c T2's maximum if this node was
         *         not reached — which is also how a source node reads.
         * @sa @ref n2dcsp_parents
         */
        T2 getNodeParent(const size_t i) const { return nodeParent[i]; }
        /**
         * @brief Record the node the ray arrived from.
         * @param index parent node index.
         * @param i     thread number.
         * @note Spelled with a lower-case @c n, unlike @ref setCellParent.
         */
        void setnodeParent(const T2 index, const size_t i) { nodeParent[i] = index; }

        /**
         * @brief Cell the ray crossed to reach this node, on a given thread.
         * @param i thread number.
         * @return Index of the traversed cell, or @c T2's maximum if unset.
         * @sa @ref n2dcsp_parents
         */
        T2 getCellParent(const size_t i) const { return cellParent[i]; }
        /**
         * @brief Record the cell the ray crossed to reach this node.
         * @param index cell index.
         * @param i     thread number.
         */
        void setCellParent(const T2 index, const size_t i) { cellParent[i] = index; }

        /**
         * @brief Mark this node as primary or secondary.
         * @param o true for a primary (corner) node, false for a secondary one.
         */
        void setPrimary( const bool o ) { primary = o; }

        /**
         * @brief Record that a cell touches this node.
         * @param o index of the owning cell.
         * @note Appends without checking for duplicates.
         */
        void pushOwner(const T2 o) { owners.push_back(o); }
        /// @return Indices of the cells touching this node.
        const std::vector<T2>& getOwners() const { return owners; }

        /// @param node other node. @return Euclidean distance to it.
        T1 getDistance( const Node2Dnsp<T1,T2>& node ) const {
            return sqrt( (x-node.x)*(x-node.x) + (z-node.z)*(z-node.z) );
        }

        /// @param node point. @return Euclidean distance to it.
        T1 getDistance( const sxz<T1>& node ) const {
            return sqrt( (x-node.x)*(x-node.x) + (z-node.z)*(z-node.z) );
        }

        /// @param node point. @return Absolute difference in x only.
        T1 getDistanceX( const sxz<T1>& node ) const {
            return std::abs( x-node.x );
        }

        /// @param node point. @return Absolute difference in z only.
        T1 getDistanceZ( const sxz<T1>& node ) const {
            return std::abs( z-node.z );
        }

        /**
         * @brief Test whether the node sits at a given location.
         * @param node point to compare against.
         * @return True if both coordinates agree to within @ref small.
         * @note Compares positions only; the tolerance is absolute.
         */
        // operator to test if same location
        bool operator==( const sxz<T1>& node ) const {
            return std::abs(x-node.x)<small && std::abs(z-node.z)<small;
        }

        /// @return Always 2 — the spatial dimension of this node.
        int getDimension() const { return 2; }

        /// @return True for a primary (corner) node, false for a secondary one.
        const bool isPrimary() const { return primary; }

    private:
        size_t nThreads;               ///< Number of slots in each per-thread array.
        T1 *tt;                        ///< travel time, one slot per thread; owned
        T1 x;                          ///< x coordinate
        T1 z;                          ///< z coordinate
        T1 slowness;                   ///< Slowness at this node; the member ttcr::Node2Dcsp lacks.
        T2 gridIndex;                  ///< index of this node in the list of the grid
        T2 *nodeParent;                ///< index of parent node of the ray, per thread; owned
        T2 *cellParent;                ///< index of cell traversed by the ray, per thread; owned
        std::vector<T2> owners;        ///< indices of cells touching the node
        /// True for a primary (corner) node, false for a secondary one.
        /// @note The trailing source comment "5= primary" is stale — this is a
        ///       @c bool, not an order code.
        bool primary;				   // indicate the order of the node: 5= primary,

    };

    /**
     * @brief Add two nodes' positions.
     * @return The sum, as a plain point — arithmetic on nodes yields
     *         @ref sxz, since a sum has no traveltime or grid identity.
     */
    template<typename T1, typename T2>
    sxz<T1> operator+(const Node2Dnsp<T1,T2>& lhs, const Node2Dnsp<T1,T2>& rhs) {
        return sxz<T1>( lhs.getX()+rhs.getX(), lhs.getZ()+rhs.getZ() );
    }

    /// @brief Vector from a node to a point. @return The difference, as a point.
    template <typename T1, typename T2>
    sxz<T1> operator-(const sxz<T1>& lhs, const Node2Dnsp<T1,T2>& rhs) {
        return sxz<T1>(lhs.x-rhs.getX(), lhs.z-rhs.getZ());
    }

    /**
     * @brief Write a node's coordinates to a stream as "x z".
     * @param os stream to write to.
     * @param s  node to write.
     * @return @p os, for chaining.
     * @note The @c struct elaborated-type-specifier is vestigial — this is a class.
     */
    template <typename T1, typename T2>
    std::ostream& operator<< (std::ostream& os, const struct Node2Dnsp<T1,T2>& s) {
        os << s.getX() << ' ' << s.getZ();
        return os;
    }

}

#endif

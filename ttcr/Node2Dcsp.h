//
//  Node2Dcsp.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-03-03.
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
 * @file Node2Dcsp.h
 * @brief Node of a 2-D shortest-path grid with cell-based slowness.
 *
 * Declares ttcr::Node2Dcsp, the node type of Grid2Drcsp. It is
 * ttcr::Node2Dc plus the bookkeeping the shortest-path method needs.
 *
 * @section n2dcsp_parents Raypath bookkeeping
 * The shortest-path method treats the grid as a graph and relaxes edges until
 * every node holds the least traveltime from the source. Recovering the ray
 * afterwards requires knowing, for each node, which neighbour it was reached
 * *from* — so alongside the traveltime this class stores two parallel arrays:
 *
 * - @ref ttcr::Node2Dcsp::nodeParent — index of the node the ray arrived from,
 * - @ref ttcr::Node2Dcsp::cellParent — index of the cell it crossed to get here.
 *
 * Walking @c nodeParent from a receiver back to the source reconstructs the
 * raypath; @c cellParent identifies which cell each segment traversed, which is
 * what the sensitivity-matrix code accumulates into.
 *
 * All three arrays are per-thread, for the reason given in @ref n2dc_threads,
 * and @ref ttcr::Node2Dcsp::reinit resets all three together.
 *
 * @sa Node2Dc.h, Node2Dnsp.h, Grid2Drcsp.h
 */

#ifndef ttcr_Node2Dcsp_h
#define ttcr_Node2Dcsp_h

#include <cmath>
#include <limits>

#include "ttcr_t.h"
#include "Node.h"

namespace ttcr {

    /**
     * @brief Node of a 2-D shortest-path grid, slowness held on the cells.
     *
     * Lives in the x-z plane; @ref getY always returns 0 so the class still
     * satisfies the 3-D ttcr::Node interface.
     *
     * @tparam T1 floating-point type of coordinates and traveltimes.
     * @tparam T2 integer type of grid, cell and parent indices, normally
     *            @c uint32_t.
     *
     * @warning Owns three raw @c new[] arrays and declares a copy constructor
     *          and a destructor but **no copy-assignment operator**, so the
     *          implicit one would shallow-copy all three pointers and lead to a
     *          triple double-@c delete[]. Nothing currently assigns nodes;
     *          treat the type as copy-constructible only.
     * @note Members are @c private here, unlike ttcr::Node2Dc's @c protected —
     *       this class is a leaf and is not derived from.
     */
    template<typename T1, typename T2>
    class Node2Dcsp : public Node<T1> {
    public:
        /**
         * @brief Construct an unpositioned node with everything unset.
         * @param nt number of threads, i.e. the size of each of the three arrays.
         * @post Traveltimes are @c T1's maximum and both parent indices are
         *       @c T2's maximum — the "not yet reached" markers.
         */
        Node2Dcsp(const size_t nt=1) :
        nThreads(nt),
        x(0.0), z(0.0),
        gridIndex(std::numeric_limits<T2>::max()),
        tt(nullptr),
        nodeParent(nullptr),
        cellParent(nullptr),
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
         * @brief Construct a positioned node with everything else unset.
         * @param xx    x coordinate.
         * @param zz    z coordinate.
         * @param index index of this node in the grid's node list.
         * @param nt    number of threads.
         */
        Node2Dcsp(const T1 xx, const T1 zz, const T2 index, const size_t nt=1) :
        nThreads(nt),
        x(xx), z(zz),
        gridIndex(index),
        tt(nullptr),
        nodeParent(nullptr),
        cellParent(nullptr),
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
         * @brief Construct a positioned node with a known traveltime on one thread.
         *
         * Used to seed a source position; the parent indices stay unset, which
         * is what marks the node as a raypath origin.
         *
         * @param t  traveltime to store.
         * @param xx x coordinate.
         * @param zz z coordinate.
         * @param nt number of threads.
         * @param i  thread number the traveltime belongs to.
         * @pre @p i is less than @p nt; the write is unchecked.
         */
        Node2Dcsp(const T1 t, const T1 xx, const T1 zz, const size_t nt, const size_t i) :
        nThreads(nt),
        x(xx), z(zz),
        gridIndex(std::numeric_limits<T2>::max()),
        tt(nullptr),
        nodeParent(nullptr),
        cellParent(nullptr),
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
         * @brief Construct from a point, with a known traveltime on one thread.
         * @tparam SXZ any type with @c x and @c z members.
         * @param t  traveltime to store.
         * @param s  point supplying the coordinates.
         * @param nt number of threads.
         * @param i  thread number the traveltime belongs to.
         * @pre @p i is less than @p nt.
         */
        template<typename SXZ>
        Node2Dcsp(const T1 t, const SXZ &s, const size_t nt, const size_t i) :
        nThreads(nt),
        x(s.x), z(s.z),
        gridIndex(std::numeric_limits<T2>::max()),
        tt(nullptr),
        nodeParent(nullptr),
        cellParent(nullptr),
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
        Node2Dcsp(const Node2Dcsp<T1,T2>& node) :
        nThreads(node.nThreads),
        x(node.x), z(node.z),
        gridIndex(node.gridIndex),
        tt(nullptr),
        nodeParent(nullptr),
        cellParent(nullptr),
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
        virtual ~Node2Dcsp() {
            delete [] tt;
            delete [] nodeParent;
            delete [] cellParent;
        }

        /**
         * @brief Reset one thread's traveltime and both parent indices.
         * @param thread_no thread whose slots are reset.
         * @note Resets all three arrays together, so a stale raypath cannot
         *       survive into the next source's solve. Other threads are
         *       untouched.
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
         * @brief Record that a cell touches this node.
         * @param o index of the owning cell.
         * @note Appends without checking for duplicates.
         */
        void pushOwner(const T2 o) { owners.push_back(o); }
        /// @return Indices of the cells touching this node.
        const std::vector<T2>& getOwners() const { return owners; }

        /// @param node other node. @return Euclidean distance to it.
        T1 getDistance( const Node2Dcsp<T1,T2>& node ) const {
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

        /**
         * @brief Approximate memory footprint of this node, in bytes.
         * @return Estimated size, including the three heap arrays and the
         *         owners vector.
         * @note The @c 2*nThreads*sizeof(T2) term accounts for @ref nodeParent
         *       and @ref cellParent, so the formula fits this class — it is
         *       ttcr::Node2Dc that inherited it without having those arrays.
         *       @ref primary is still unaccounted for.
         */
        size_t getSize() const {
            return sizeof(size_t) + nThreads*sizeof(T1) + 2*sizeof(T1) +
            (1+2*nThreads)*sizeof(T2) + owners.size() * sizeof(T2);
        }

        /// @return Always 2 — the spatial dimension of this node.
        int getDimension() const { return 2; }

        /**
         * @brief Mark this node as primary or secondary.
         * @param p true for a primary (corner) node, false for a secondary one.
         */
        void setPrimary(const bool p) { primary = p; }
        /// @return True for a primary (corner) node, false for a secondary one.
        const bool isPrimary() const { return primary; }

    private:
        size_t nThreads;               ///< Number of slots in each per-thread array.
        T1 x;                          ///< x coordinate
        T1 z;                          ///< z coordinate
        T2 gridIndex;                  ///< index of this node in the list of the grid
        T1 *tt;                        ///< travel time, one slot per thread; owned
        T2 *nodeParent;                ///< index of parent node of the ray, per thread; owned
        T2 *cellParent;                ///< index of cell traversed by the ray, per thread; owned
        std::vector<T2> owners;        ///< indices of cells touching the node
        bool primary;                  ///< True for a primary (corner) node, false for a secondary one.

    };

    /**
     * @brief Write a node's coordinates to a stream as "x z".
     * @param os stream to write to.
     * @param n  node to write.
     * @return @p os, for chaining.
     */
    template<typename T1, typename T2>
    std::ostream& operator<< (std::ostream& os, const Node2Dcsp<T1, T2> &n) {
        os << n.getX() << ' ' << n.getZ();
        return os;
    }

    /**
     * @brief Vector from a node to a point.
     * @return The difference, as a point.
     * @note Only this one subtraction overload is provided, unlike
     *       ttcr::Node2Dc which also offers node-node and node-point forms.
     */
    template <typename T1, typename T2>
    sxz<T1> operator-(const sxz<T1>& lhs, const Node2Dcsp<T1,T2>& rhs) {
        return sxz<T1>(lhs.x-rhs.getX(), lhs.z-rhs.getZ());
    }
}

#endif

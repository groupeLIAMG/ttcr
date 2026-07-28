//
//  Node3Dnsp.h
//  ttcr
//
//  Created by Bernard Giroux on 12-08-14.
//  Copyright (c) 2012 Bernard Giroux. All rights reserved.
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

/**
 * @file Node3Dnsp.h
 * @brief Node of a 3-D shortest-path grid with node-based slowness.
 *
 * Declares ttcr::Node3Dnsp, the node type of Grid3Drnsp and Grid3Dunsp. It is
 * to ttcr::Node3Dcsp what ttcr::Node3Dn is to ttcr::Node3Dc: the same
 * shortest-path raypath bookkeeping (see @ref n2dcsp_parents), plus a slowness
 * value carried on the node.
 *
 * @sa Node3Dn.h, Node3Dcsp.h, Node2Dnsp.h, Grid3Drnsp.h, Grid3Dunsp.h
 */

#ifndef ttcr_Node3Dnsp_h
#define ttcr_Node3Dnsp_h

#include <cmath>
#include <limits>

#include "Node.h"

namespace ttcr {

    template<typename T1, typename T2>
    /**
     * @brief Node of a 3-D shortest-path grid, slowness held at the node.
     *
     * @tparam T1 floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2 integer type of grid, cell and parent indices, normally
     *            @c uint32_t.
     *
     * @warning Owns three raw @c new[] arrays and declares a copy constructor
     *          and a destructor but **no copy-assignment operator** — see the
     *          same warning on ttcr::Node3Dcsp. Copy-construct, never assign.
     * @note Unlike ttcr::Node3Dcsp, the first constructor takes no default
     *       argument: a thread count must always be supplied.
     */
    class Node3Dnsp : public Node<T1> {
    public:
        /**
         * @brief Construct an unpositioned node with everything unset.
         * @param nt number of threads, i.e. the size of each of the three arrays.
         * @post Coordinates and slowness are 0; traveltimes and both parent
         *       indices hold their type's maximum, the "not yet reached" marker.
         */
        Node3Dnsp(const size_t nt) :
        nThreads(nt),
        tt(new T1[nt]),
        x(0.0f), y(0.0f), z(0.0f),
        gridIndex(std::numeric_limits<T2>::max()),
        nodeParent(new T2[nt]),
        cellParent(new T2[nt]),
        owners(std::vector<T2>(0)),
        slowness(0),
        primary(0)
        {
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
         * @param yy y coordinate.
         * @param zz z coordinate.
         * @param nt number of threads.
         * @param i  thread number the traveltime belongs to.
         * @pre @p i is less than @p nt; the write is unchecked.
         * @note Slowness is left at 0 and must be set separately with
         *       @ref setNodeSlowness.
         */
        Node3Dnsp(const T1 t, const T1 xx, const T1 yy, const T1 zz, const size_t nt,
                  const size_t i) :
        nThreads(nt),
        tt(new T1[nt]),
        x(xx), y(yy), z(zz),
        gridIndex(std::numeric_limits<T2>::max()),
        nodeParent(new T2[nt]),
        cellParent(new T2[nt]),
        owners(std::vector<T2>(0)),
        slowness(0),
        primary(0)
        {
            for ( size_t n=0; n<nt; ++n ) {
                tt[n] = std::numeric_limits<T1>::max();
                nodeParent[n] = std::numeric_limits<T2>::max();
                cellParent[n] = std::numeric_limits<T2>::max();
            }
            tt[i]=t;
        }

        /**
         * @brief Construct from a point, with a known traveltime on one thread.
         * @param t  traveltime to store.
         * @param s  point supplying the coordinates.
         * @param nt number of threads.
         * @param i  thread number the traveltime belongs to.
         * @pre @p i is less than @p nt.
         */
        Node3Dnsp(const T1 t, const sxyz<T1>& s, const size_t nt,
                  const size_t i) :
        nThreads(nt),
        tt(new T1[nt]),
        x(s.x), y(s.y), z(s.z),
        gridIndex(std::numeric_limits<T2>::max()),
        nodeParent(new T2[nt]),
        cellParent(new T2[nt]),
        owners(std::vector<T2>(0)),
        slowness(0),
        primary(0)
        {
            for ( size_t n=0; n<nt; ++n ) {
                tt[n] = std::numeric_limits<T1>::max();
                nodeParent[n] = std::numeric_limits<T2>::max();
                cellParent[n] = std::numeric_limits<T2>::max();
            }
            tt[i]=t;
        }

        /**
         * @brief Copy constructor; deep-copies all three per-thread arrays.
         * @param node node to copy.
         */
        Node3Dnsp(const Node3Dnsp<T1,T2>& node) :
        nThreads(node.nThreads),
        tt(0),
        x(node.x), y(node.y), z(node.z),
        gridIndex(node.gridIndex),
        nodeParent(0),
        cellParent(0),
        owners(node.owners),
        slowness(node.slowness),
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
        virtual ~Node3Dnsp() {
            delete [] tt;
            delete [] nodeParent;
            delete [] cellParent;
        }

        /**
         * @brief Reset one thread's traveltime and both parent indices.
         * @param n thread whose slots are reset.
         * @note Resets all three arrays together, so a stale raypath cannot
         *       survive into the next source's solve.
         */
        // Sets the vectors to the right size of threads and initialize it
        void reinit(const size_t n) {
            tt[n] = std::numeric_limits<T1>::max();
            nodeParent[n] = std::numeric_limits<T2>::max();
            cellParent[n] = std::numeric_limits<T2>::max();
        }

        /// @param n thread number. @return Traveltime stored for that thread.
        T1 getTT(const size_t n) const { return tt[n]; }
        /// @param t traveltime. @param n thread number. Stores @p t for that thread.
        void setTT(const T1 t, const size_t n ) { tt[n] = t; }

        /**
         * @brief Set position and grid index in one call.
         * @param xx    x coordinate.
         * @param yy    y coordinate.
         * @param zz    z coordinate.
         * @param index index of this node in the grid's node list.
         */
        void setXYZindex(const T1 xx, const T1 yy, const T1 zz, const T2 index) {
            x=xx; y=yy; z=zz; gridIndex = index;  }
        /**
         * @brief Set position from a point, plus the grid index.
         * @param s     point supplying the coordinates.
         * @param index index of this node in the grid's node list.
         */
        void setXYZindex(const sxyz<T1>& s, const T2 index) {
            x=s.x; y=s.y; z=s.z; gridIndex = index;  }

        T1 getX() const { return x; }       ///< @return x coordinate.
        void setX(const T1 xx) { x = xx; }  ///< @param xx new x coordinate.

        T1 getY() const { return y; }       ///< @return y coordinate.
        void setY(const T1 yy) { y = yy; }  ///< @param yy new y coordinate.

        T1 getZ() const { return z; }       ///< @return z coordinate.
        void setZ(const T1 zz) { z = zz; }  ///< @param zz new z coordinate.

        /// @return Index of this node in the grid's node list.
        T2 getGridIndex() const { return gridIndex; }
        /// @param index new index of this node in the grid's node list.
        void setGridIndex(const T2 index) { gridIndex = index; }

        /**
         * @brief Node the ray arrived from, on a given thread.
         * @param n thread number.
         * @return Index of the parent node, or @c T2's maximum if this node was
         *         not reached — which is also how a source node reads.
         * @sa @ref n2dcsp_parents
         */
        T2 getNodeParent(const size_t n) const { return nodeParent[n]; }
        /**
         * @brief Record the node the ray arrived from.
         * @param index parent node index.
         * @param n     thread number.
         * @note Spelled with a lower-case @c n, unlike @ref setCellParent.
         */
        void setnodeParent(const T2 index, const size_t n) { nodeParent[n] = index; }

        /**
         * @brief Cell the ray crossed to reach this node, on a given thread.
         * @param n thread number.
         * @return Index of the traversed cell, or @c T2's maximum if unset.
         * @sa @ref n2dcsp_parents
         */
        T2 getCellParent(const size_t n) const { return cellParent[n]; }
        /**
         * @brief Record the cell the ray crossed to reach this node.
         * @param index cell index.
         * @param n     thread number.
         */
        void setCellParent(const T2 index, const size_t n) { cellParent[n] = index; }

        /**
         * @brief Mark this node as primary or secondary.
         * @param o nonzero for a primary (corner) node, 0 for a secondary one.
         */
        void setPrimary( const bool o ) { primary = o; }

        /**
         * @brief Slowness at this node.
         * @return Slowness, i.e. reciprocal velocity.
         * @note This accessor is what distinguishes ttcr::Node3Dnsp from
         *       ttcr::Node3Dcsp.
         */
        T1 getNodeSlowness() const { return slowness; }
        /**
         * @brief Set the slowness at this node.
         * @param s slowness value.
         * @warning Unvalidated: a zero or negative slowness is stored as given.
         */
        void setNodeSlowness(const T1 s) { slowness = s; }

        /**
         * @brief Record that a cell touches this node.
         * @param o index of the owning cell.
         * @note Appends without checking for duplicates.
         */
        void pushOwner(const T2 o) { owners.push_back(o); }
        /// @return Indices of the cells touching this node.
        const std::vector<T2>& getOwners() const { return owners; }

        /// @param node other node. @return Euclidean distance to it.
        T1 getDistance( const Node3Dnsp<T1,T2>& node ) const {
            return sqrt( (x-node.x)*(x-node.x) + (y-node.y)*(y-node.y) + (z-node.z)*(z-node.z) );
        }

        /// @param node point. @return Euclidean distance to it.
        T1 getDistance( const sxyz<T1>& node ) const {
            return sqrt( (x-node.x)*(x-node.x) + (y-node.y)*(y-node.y) + (z-node.z)*(z-node.z) );
        }

        /**
         * @brief Test whether the node sits at a given location.
         * @param node point to compare against.
         * @return True if all three coordinates agree to within @ref small.
         * @note Compares positions only; the tolerance is absolute.
         */
        // operator to test if same location
        bool operator==( const sxyz<T1>& node ) const {
            return std::abs(x-node.x)<small && std::abs(y-node.y)<small && std::abs(z-node.z)<small;
        }

        /**
         * @brief Translate the node by subtracting a vector, in place.
         * @param node offset to subtract.
         * @return Reference to this node, for chaining.
         * @note Used to shift every node when
         *       ttcr::input_parameters::translateOrigin is set.
         */
        Node3Dnsp<T1, T2>& operator-=(const sxyz<T1>& node) {
            this->x -= node.x;
            this->y -= node.y;
            this->z -= node.z;
            return *this;
        }

        /**
         * @brief Approximate memory footprint of this node, in bytes.
         * @return Estimated size, including the three heap arrays and the
         *         owners vector.
         * @note The @c 4*sizeof(T1) term covers x, y, z and slowness, and the
         *       @c 2*nThreads*sizeof(T2) term the two parent arrays, so the
         *       formula fits this class. @ref primary is unaccounted for.
         */
        size_t getSize() const {
            return 2*sizeof(size_t) + nThreads*sizeof(T1) + 4*sizeof(T1) +
            (1+2*nThreads)*sizeof(T2) + owners.size() * sizeof(T2);
        }

        /// @return Always 3 — the spatial dimension of this node.
        int getDimension() const { return 3; }

        /// @return True for a primary (corner) node, false for a secondary one.
        bool const isPrimary() const { return primary; }

    private:
        size_t nThreads;                ///< Number of slots in each per-thread array.
        T1 *tt;                         ///< travel time for the multiple source points; owned
        T1 x;                           ///< x coordinate [km]
        T1 y;							///< y coordinate [km]
        T1 z;                           ///< z coordinate [km]
        T2 gridIndex;                   ///< index of this node in the list of the grid
        T2 *nodeParent;                 ///< index of parent node of the ray for each thread; owned
        T2 *cellParent;                 ///< index of cell traversed by the ray for each thread; owned
        std::vector<T2> owners;         ///< indices of cells touching the node
        /// Slowness at the node [s/km].
        T1 slowness;					// slowness at the node [s/km]
        /// True for a primary (corner) node, false for a secondary one.
        /// @note The trailing comment "indicate the order of the node" predates
        ///       this being a @c bool.
        bool primary;

    };

    /**
     * @brief Add two nodes' positions.
     * @return The sum, as a plain point — arithmetic on nodes yields
     *         @ref sxyz, since a sum has no traveltime or grid identity.
     */
    template<typename T1, typename T2>
    sxyz<T1> operator+(const Node3Dnsp<T1,T2>& lhs, const Node3Dnsp<T1,T2>& rhs) {
        return sxyz<T1>( lhs.getX()+rhs.getX(), lhs.getY()+rhs.getY(), lhs.getZ()+rhs.getZ() );
    }
    /**
     * @brief Write a node's coordinates to a stream as "x y z".
     * @param os stream to write to.
     * @param n  node to write.
     * @return @p os, for chaining.
     */
    template<typename T1, typename T2>
    std::ostream& operator<< (std::ostream& os, const Node3Dnsp<T1, T2> &n) {
        os << n.getX() << ' ' << n.getY() << ' ' << n.getZ();
        return os;
    }

    /// @brief Vector from @p rhs to @p lhs. @return The difference, as a point.
    template<typename T1, typename T2>
    sxyz<T1> operator-(const Node3Dnsp<T1,T2>& lhs, const Node3Dnsp<T1,T2>& rhs) {
        return sxyz<T1>( lhs.getX()-rhs.getX(), lhs.getY()-rhs.getY(), lhs.getZ()-rhs.getZ() );
    }

    /// @brief Vector from a node to a point. @return The difference, as a point.
    template<typename T1, typename T2>
    sxyz<T1> operator-(const sxyz<T1>& lhs, const Node3Dnsp<T1,T2>& rhs) {
        return sxyz<T1>( lhs.x-rhs.getX(), lhs.y-rhs.getY(), lhs.z-rhs.getZ() );
    }

    /// @brief Vector from a point to a node. @return The difference, as a point.
    template<typename T1, typename T2>
    sxyz<T1> operator-(const Node3Dnsp<T1,T2>& lhs, const sxyz<T1>& rhs) {
        return sxyz<T1>( lhs.getX()-rhs.x, lhs.getY()-rhs.y, lhs.getZ()-rhs.z );
    }

}

#endif

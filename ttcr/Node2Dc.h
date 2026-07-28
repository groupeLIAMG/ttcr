//
// Node2Dc.h
// ttcr
//
// Created by Bernard Giroux on 08-04-24.
// Copyright 2008 Bernard Giroux.
//
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
 * @file Node2Dc.h
 * @brief Node of a 2-D grid whose slowness is constant within each cell.
 *
 * Declares ttcr::Node2Dc, the node type used by the 2-D grid classes that store
 * slowness **per cell** rather than per node — the @c uc family
 * (Grid2Ducfs, Grid2Ducfm, Grid2Ducdsp) and the cell-based rectilinear
 * Grid2Drcdsp. The node itself therefore carries no slowness value: it holds a
 * position, a traveltime, and the indices of the cells it belongs to.
 *
 * @section n2dc_family The Node2D family
 * Six closely-related classes, differing along two axes:
 *
 * | | cell slowness | node slowness |
 * |---|---|---|
 * | plain | ttcr::Node2Dc | ttcr::Node2Dn |
 * | + shortest-path bookkeeping | ttcr::Node2Dcsp | ttcr::Node2Dnsp |
 * | + dynamic (temporary) | ttcr::Node2Dcd | ttcr::Node2Dnd |
 *
 * The @c n variants add a @c slowness member; the @c sp variants add the parent
 * node and cell arrays the shortest-path method needs to walk a raypath back to
 * its source; the @c d variants are lightweight temporaries used by the dynamic
 * shortest-path solvers.
 *
 * @section n2dc_threads Per-thread traveltimes
 * The traveltime is not a scalar but an array of @c nThreads entries. Each
 * worker thread propagates a different source through the *same* node objects
 * and writes only its own slot, so the grid is shared without locking. Every
 * accessor that touches a traveltime therefore takes a thread number, and
 * passing the wrong one silently reads another source's solution.
 *
 * @sa Node.h, Node2Dn.h, Node2Dcsp.h, Node2Dcd.h, Grid2Duc.h
 */

#ifndef ttcr_NodeDc_h
#define ttcr_NodeDc_h

#include <cmath>
#include <limits>

#include "ttcr_t.h"
#include "Node.h"

namespace ttcr {

    /**
     * @brief Node of a 2-D grid with cell-based slowness.
     *
     * Lives in the x-z plane; @ref getY always returns 0 so the class still
     * satisfies the 3-D ttcr::Node interface.
     *
     * @tparam T1 floating-point type of coordinates and traveltimes.
     * @tparam T2 integer type of grid and cell indices, normally @c uint32_t.
     *
     * @warning Owns a raw @c new[] array and declares a copy constructor and a
     *          destructor but **no copy-assignment operator**, so the implicitly
     *          generated one copies the pointer rather than the array —
     *          assigning one node to another would give two objects owning the
     *          same buffer and a double @c delete[]. Nothing in the project
     *          currently assigns nodes, so this is latent, but treat the type as
     *          copy-constructible only. Note that because the copy constructor
     *          is user-declared no move operations are generated either, so
     *          @c std::vector reallocation goes through the (correct, deep)
     *          copy constructor.
     */
    template<typename T1, typename T2>
    class Node2Dc : public Node<T1> {
    public:
        /**
         * @brief Construct an unpositioned node with unset traveltimes.
         * @param nt number of threads, i.e. the size of the traveltime array.
         * @post Coordinates are 0, the grid index is @c T2's maximum, and every
         *       traveltime is @c T1's maximum — the "not yet reached" marker the
         *       solvers test against.
         */
        Node2Dc(const size_t nt=1) :
        nThreads(nt),
        x(0.0), z(0.0),
        gridIndex(std::numeric_limits<T2>::max()),
        tt(nullptr),
        owners(0),
        primary(false)
        {
            tt = new T1[nt];

            for ( size_t n=0; n<nt; ++n ) {
                tt[n] = std::numeric_limits<T1>::max();
            }
        }

        /**
         * @brief Construct a positioned node with unset traveltimes.
         * @param xx    x coordinate.
         * @param zz    z coordinate.
         * @param index index of this node in the grid's node list.
         * @param nt    number of threads.
         */
        Node2Dc(const T1 xx, const T1 zz, const T2 index, const size_t nt=1) :
        nThreads(nt),
        x(xx), z(zz),
        gridIndex(index),
        tt(nullptr),
        owners(0),
        primary(false)
        {
            tt = new T1[nt];

            for ( size_t n=0; n<nt; ++n ) {
                tt[n] = std::numeric_limits<T1>::max();
            }
        }


        /**
         * @brief Construct a positioned node with a known traveltime on one thread.
         *
         * Used to seed a source position: every slot is left at "not reached"
         * except slot @p i, which receives @p t.
         *
         * @param t  traveltime to store.
         * @param xx x coordinate.
         * @param zz z coordinate.
         * @param nt number of threads.
         * @param i  thread number the traveltime belongs to.
         * @pre @p i is less than @p nt; the write is unchecked.
         * @note Leaves the grid index unset (@c T2's maximum) — unlike the
         *       positioned constructor above, this one is for nodes that are not
         *       part of the grid's node list.
         */
        Node2Dc(const T1 t, const T1 xx, const T1 zz, const size_t nt, const size_t i) :
        nThreads(nt),
        x(xx), z(zz),
        gridIndex(std::numeric_limits<T2>::max()),
        tt(nullptr),
        owners(std::vector<T2>(0)),
        primary(false)
        {
            tt = new T1[nt];

            for ( size_t n=0; n<nt; ++n ) {
                tt[n] = std::numeric_limits<T1>::max();
            }
            tt[i] = t;
        }

        /**
         * @brief Copy constructor; deep-copies the traveltime array.
         * @param node node to copy.
         */
        Node2Dc(const Node2Dc<T1,T2>& node) :
        nThreads(node.nThreads),
        x(node.x), z(node.z),
        gridIndex(node.gridIndex),
        tt(nullptr),
        owners(node.owners),
        primary(node.primary)
        {
            tt = new T1[nThreads];

            for ( size_t n=0; n<nThreads; ++n ) {
                tt[n] = node.tt[n];
            }
        }


        /// Destructor; frees the traveltime array.
        virtual ~Node2Dc() {
            delete [] tt;
        }

        /**
         * @brief Reset one thread's traveltime to "not reached".
         * @param thread_no thread whose slot is reset.
         * @note Called before each new source is propagated. Only the given
         *       slot is touched, so other threads keep their solutions.
         */
        void reinit(const size_t thread_no) { //=0) {
            tt[thread_no] = std::numeric_limits<T1>::max();
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
            x=xx; z=zz; gridIndex = index; }

        /**
         * @brief Set position from a point object, plus the grid index.
         * @tparam SXZ any type with @c x and @c z members.
         * @param s     point supplying the coordinates.
         * @param index index of this node in the grid's node list.
         * @note Named @c setXYZindex for symmetry with the 3-D nodes; any @c y
         *       the argument carries is ignored.
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
         * @brief Record that a cell touches this node.
         * @param o index of the owning cell.
         * @note Appends without checking for duplicates; the grid builder is
         *       responsible for calling this once per adjacent cell.
         */
        void pushOwner(const T2 o) { owners.push_back(o); }
        /// @return Indices of the cells touching this node.
        const std::vector<T2>& getOwners() const { return owners; }

        /// @param node other node. @return Euclidean distance to it.
        T1 getDistance( const Node2Dc<T1,T2>& node ) const {
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
         * @note Compares positions only — traveltimes, index and owners are
         *       ignored. The tolerance is absolute (@f$10^{-4}@f$ model units),
         *       not relative.
         */
        // operator to test if same location
        bool operator==( const sxz<T1>& node ) const {
            return std::abs(x-node.x)<small && std::abs(z-node.z)<small;
        }

        /**
         * @brief Approximate memory footprint of this node, in bytes.
         * @return Estimated size, including the heap-allocated traveltime array
         *         and the owners vector.
         * @warning The formula is inherited from ttcr::Node2Dcsp and
         *          **overestimates** here: its @c 2*nThreads*sizeof(T2) term
         *          accounts for that class's @c nodeParent and @c cellParent
         *          arrays, which this class does not have. It also omits
         *          @ref primary. Treat the value as indicative only.
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
        void setPrimary(const bool p=true) {
            primary = p;
        }
        /// @return True for a primary (corner) node, false for a secondary one.
        const bool isPrimary() const { return primary; }

    protected:
        size_t nThreads;               ///< Number of traveltime slots, i.e. concurrent sources.
        T1 x;                          ///< x coordinate
        T1 z;                          ///< z coordinate
        T2 gridIndex;                  ///< index of this node in the list of the grid
        T1 *tt;                        ///< travel time, one slot per thread; owned, @c new[] allocated
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
    std::ostream& operator<< (std::ostream& os, const Node2Dc<T1, T2> &n) {
        os << n.getX() << ' ' << n.getZ();
        return os;
    }

    /**
     * @brief Add two nodes' positions.
     * @return The sum, as a plain point — arithmetic on nodes yields
     *         @ref sxz, since a sum has no traveltime or grid identity.
     */
    template<typename T1, typename T2>
    sxz<T1> operator+(const Node2Dc<T1,T2>& lhs, const Node2Dc<T1,T2>& rhs) {
        return sxz<T1>( lhs.getX()+rhs.getX(), lhs.getZ()+rhs.getZ() );
    }

    /// @brief Vector from @p rhs to @p lhs. @return The difference, as a point.
    template<typename T1, typename T2>
    sxz<T1> operator-(const Node2Dc<T1,T2>& lhs, const Node2Dc<T1,T2>& rhs) {
        return sxz<T1>( lhs.getX()-rhs.getX(), lhs.getZ()-rhs.getZ() );
    }

    /// @brief Vector from a node to a point. @return The difference, as a point.
    template<typename T1, typename T2>
    sxz<T1> operator-(const sxz<T1>& lhs, const Node2Dc<T1,T2>& rhs) {
        return sxz<T1>( lhs.x-rhs.getX(), lhs.z-rhs.getZ() );
    }

    /// @brief Vector from a point to a node. @return The difference, as a point.
    template<typename T1, typename T2>
    sxz<T1> operator-(const Node2Dc<T1,T2>& lhs, const sxz<T1>& rhs) {
        return sxz<T1>( lhs.getX()-rhs.x, lhs.getZ()-rhs.z);
    }

}

#endif

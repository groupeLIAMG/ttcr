//
//  Node3Dc.h
//  ttcr
//
//  Created by Bernard Giroux on 08-04-24.
//
//  Modified by Benoit Larouche on 12-07-12
//  	: now support parallel raytracing from many source points
//  	  on the same 3D grid simultaneously, using OpenMP.
//  	  Velocity model is sampled at the primary nodes
//  	  Secondary nodes velocity is linearly interpolated
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
 * @file Node3Dc.h
 * @brief Node of a 3-D grid whose slowness is constant within each cell.
 *
 * Declares ttcr::Node3Dc, the node type used by the 3-D grid classes that store
 * slowness **per cell** — Grid3Ducfs, Grid3Ducfm, Grid3Ducdsp and the
 * cell-based rectilinear Grid3Drcdsp.
 *
 * @section n3dc_family The Node3D family
 * The 3-D counterpart of the Node2D family, with the same 2x3 structure:
 *
 * | | cell slowness | node slowness |
 * |---|---|---|
 * | plain | ttcr::Node3Dc | ttcr::Node3Dn |
 * | + shortest-path bookkeeping | ttcr::Node3Dcsp | ttcr::Node3Dnsp |
 * | + dynamic (temporary) | ttcr::Node3Dcd | ttcr::Node3Dnd |
 *
 * Differences from the 2-D classes, beyond carrying a @c y coordinate:
 * - ttcr::Node3Dc::getDimension returns 3 and ttcr::Node3Dc::getY returns a
 *   real coordinate rather than 0;
 * - a mutating @ref ttcr::Node3Dc::operator-=() is provided, used to shift every
 *   node when the grid origin is translated;
 * - there are no @c getDistanceX / @c getDistanceZ helpers.
 *
 * Per-thread traveltimes work exactly as in the 2-D family: @c tt is an array
 * with one slot per worker thread, so several sources can be propagated through
 * the same shared nodes without locking. See @ref n2dc_threads.
 *
 * @sa Node.h, Node2Dc.h, Node3Dn.h, Node3Dcsp.h, Node3Dcd.h, Grid3Duc.h
 */

#ifndef ttcr_Node3Dc_h
#define ttcr_Node3Dc_h

#include <cmath>
#include <limits>

#include "Node.h"

namespace ttcr {

    /**
     * @brief Node of a 3-D grid with cell-based slowness.
     *
     * @tparam T1 floating-point type of coordinates and traveltimes.
     * @tparam T2 integer type of grid and cell indices, normally @c uint32_t.
     *
     * @warning Owns a raw @c new[] array and declares a copy constructor and a
     *          destructor but **no copy-assignment operator**, so the implicitly
     *          generated one copies the pointer rather than the array — two
     *          objects would own one buffer and double-@c delete[] it. Nothing
     *          in the project assigns nodes, so this is latent; treat the type
     *          as copy-constructible only.
     */
    template<typename T1, typename T2>
    class Node3Dc : public Node<T1> {
    public:
        /**
         * @brief Construct an unpositioned node with unset traveltimes.
         * @param nt number of threads, i.e. the size of the traveltime array.
         * @post Coordinates are 0, the grid index is @c T2's maximum, every
         *       traveltime is @c T1's maximum (the "not yet reached" marker),
         *       and @ref primary is false.
         */
        Node3Dc(const size_t nt=1) :
        nThreads(nt),
        x(0.0f), y(0.0f), z(0.0f),
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
         * @param yy    y coordinate.
         * @param zz    z coordinate.
         * @param index index of this node in the grid's node list.
         * @param nt    number of threads.
         */
        Node3Dc(const T1 xx, const T1 yy, const T1 zz, const T2 index,
                const size_t nt) :
        nThreads(nt),
        x(xx), y(yy), z(zz),
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
         * @param yy y coordinate.
         * @param zz z coordinate.
         * @param nt number of threads.
         * @param i  thread number the traveltime belongs to.
         * @pre @p i is less than @p nt; the write is unchecked.
         * @note Leaves the grid index unset — this form is for nodes that are
         *       not part of the grid's node list.
         */
        Node3Dc(const T1 t, const T1 xx, const T1 yy, const T1 zz, const size_t nt,
                const size_t i) :
        nThreads(nt),
        x(xx), y(yy), z(zz),
        gridIndex(std::numeric_limits<T2>::max()),
        tt(nullptr),
        owners(0),
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
        Node3Dc(const Node3Dc<T1,T2>& node) :
        nThreads(node.nThreads),
        x(node.x), y(node.y), z(node.z),
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
        virtual ~Node3Dc() {
            delete [] tt;
        }

        /**
         * @brief Reset one thread's traveltime to "not reached".
         * @param n thread whose slot is reset.
         * @note Called before each new source is propagated; other threads keep
         *       their solutions.
         */
        // Sets the vectors to the right size of threads and initialize it
        void reinit(const size_t n) {
            tt[n] = std::numeric_limits<T1>::max();
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
         * @brief Record that a cell touches this node.
         * @param o index of the owning cell.
         * @note Appends without checking for duplicates.
         */
        void pushOwner(const T2 o) { owners.push_back(o); }
        /// @return Indices of the cells touching this node.
        const std::vector<T2>& getOwners() const { return owners; }

        /// @param node other node. @return Euclidean distance to it.
        T1 getDistance( const Node3Dc<T1,T2>& node ) const {
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
         * @note Compares positions only — traveltimes, index and owners are
         *       ignored. The tolerance is absolute, in model units.
         */
        // operator to test if same location
        bool operator==( const sxyz<T1>& node ) const {
            return std::abs(x-node.x)<small && std::abs(y-node.y)<small && std::abs(z-node.z)<small;
        }

        /**
         * @brief Translate the node by subtracting a vector, in place.
         * @param node offset to subtract.
         * @return Reference to this node, for chaining.
         * @note This is how a grid shifts its nodes when
         *       ttcr::input_parameters::translateOrigin is set — the grid
         *       classes apply @c *node -= origin to every node. It mutates the
         *       node, unlike the free @c operator- overloads, which return a
         *       plain @ref sxyz.
         */
        Node3Dc<T1, T2>& operator-=(const sxyz<T1>& node) {
            this->x -= node.x;
            this->y -= node.y;
            this->z -= node.z;
            return *this;
        }

        /**
         * @brief Approximate memory footprint of this node, in bytes.
         * @return Estimated size, including the heap-allocated traveltime array
         *         and the owners vector.
         * @warning **Overestimates**: the @c 2*nThreads*sizeof(T2) term accounts
         *          for the @c nodeParent and @c cellParent arrays of
         *          ttcr::Node3Dcsp, which this class does not have. It also
         *          omits @ref primary. Indicative only.
         */
        size_t getSize() const {
            return sizeof(size_t) + nThreads*sizeof(T1) + 3*sizeof(T1) +
            (1+2*nThreads)*sizeof(T2) + owners.size() * sizeof(T2);
        }

        /// @return Always 3 — the spatial dimension of this node.
        int getDimension() const { return 3; }

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
        size_t nThreads;            ///< Number of traveltime slots, i.e. concurrent sources.
        T1 x;                       ///< x coordinate [km]
        T1 y;						///< y coordinate [km]
        T1 z;                       ///< z coordinate [km]
        T2 gridIndex;               ///< index of this node in the list of the grid
        T1 *tt;                     ///< travel time for the multiple source points; owned, @c new[] allocated
        std::vector<T2> owners;     ///< indices of cells touching the node
        bool primary;               ///< True for a primary (corner) node; defaults to true in this class.
    };

    /**
     * @brief Add two nodes' positions.
     * @return The sum, as a plain point — arithmetic on nodes yields
     *         @ref sxyz, since a sum has no traveltime or grid identity.
     */
    template<typename T1, typename T2>
    sxyz<T1> operator+(const Node3Dc<T1,T2>& lhs, const Node3Dc<T1,T2>& rhs) {
        return sxyz<T1>( lhs.getX()+rhs.getX(), lhs.getY()+rhs.getY(), lhs.getZ()+rhs.getZ() );
    }

    /// @brief Vector from @p rhs to @p lhs. @return The difference, as a point.
    template<typename T1, typename T2>
    sxyz<T1> operator-(const Node3Dc<T1,T2>& lhs, const Node3Dc<T1,T2>& rhs) {
        return sxyz<T1>( lhs.getX()-rhs.getX(), lhs.getY()-rhs.getY(), lhs.getZ()-rhs.getZ() );
    }

    /// @brief Vector from a node to a point. @return The difference, as a point.
    template<typename T1, typename T2>
    sxyz<T1> operator-(const sxyz<T1>& lhs, const Node3Dc<T1,T2>& rhs) {
        return sxyz<T1>( lhs.x-rhs.getX(), lhs.y-rhs.getY(), lhs.z-rhs.getZ() );
    }

    /// @brief Vector from a point to a node. @return The difference, as a point.
    template<typename T1, typename T2>
    sxyz<T1> operator-(const Node3Dc<T1,T2>& lhs, const sxyz<T1>& rhs) {
        return sxyz<T1>( lhs.getX()-rhs.x, lhs.getY()-rhs.y, lhs.getZ()-rhs.z);
    }

    /**
     * @brief Write a node's coordinates to a stream as "x y z".
     * @param os stream to write to.
     * @param n  node to write.
     * @return @p os, for chaining.
     */
    template<typename T1, typename T2>
    std::ostream& operator<< (std::ostream& os, const Node3Dc<T1, T2> &n) {
        os << n.getX() << ' ' << n.getY() << ' ' << n.getZ();
        return os;
    }

}

#endif

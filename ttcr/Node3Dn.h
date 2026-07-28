//
//  Node3Dn.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-04-21.
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
 * @file Node3Dn.h
 * @brief Node of a 3-D grid whose slowness is defined at the nodes.
 *
 * Declares ttcr::Node3Dn, the counterpart of ttcr::Node3Dc for grids that store
 * slowness **per node** — the @c rn and @c un families (Grid3Drnfs,
 * Grid3Drndsp, Grid3Dunfs, Grid3Dunfm, Grid3Dundsp and the OpenCL variants),
 * plus the node-based rectilinear Grid3Drcfs.
 *
 * It differs from ttcr::Node3Dc by carrying a @ref ttcr::Node3Dn::slowness
 * member, reachable through @ref ttcr::Node3Dn::getNodeSlowness, and by
 * defaulting @c primary to false rather than true. See @ref n3dc_family for how
 * the six Node3D classes relate and @ref n2dc_threads for the per-thread
 * traveltime convention.
 *
 * @sa Node.h, Node3Dc.h, Node3Dnsp.h, Node3Dnd.h, Node2Dn.h, Grid3Dun.h
 */

#ifndef ttcr_Node3Dn_h
#define ttcr_Node3Dn_h

#include <cmath>
#include <limits>
#include <vector>

#include "Node.h"

namespace ttcr {

    /**
     * @brief Node of a 3-D grid with node-based slowness.
     *
     * @tparam T1 floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2 integer type of grid and cell indices, normally @c uint32_t.
     *
     * @warning Owns a raw @c new[] array and declares a copy constructor and a
     *          destructor but **no copy-assignment operator** — see the same
     *          warning on ttcr::Node3Dc. Copy-construct, never assign.
     * @note Unlike ttcr::Node3Dc, the default constructor here takes no default
     *       argument: a thread count must always be supplied.
     */
    template<typename T1, typename T2>
    class Node3Dn : public Node<T1> {
    public:
        /**
         * @brief Construct an unpositioned node with unset traveltimes.
         * @param nt number of threads, i.e. the size of the traveltime array.
         * @post Coordinates and slowness are 0, the grid index is @c T2's
         *       maximum, every traveltime is @c T1's maximum, and @ref primary
         *       is false.
         */
        Node3Dn(const size_t nt) :
        nThreads(nt),
        tt(new T1[nt]),
        x(0.0f), y(0.0f), z(0.0f),
        gridIndex(std::numeric_limits<T2>::max()),
        owners(std::vector<T2>(0)),
        slowness(0), primary(0)
        {
            for ( size_t n=0; n<nt; ++n ) {
                tt[n] = std::numeric_limits<T1>::max();
            }
        }

        /**
         * @brief Construct a positioned node with a known traveltime on one thread.
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
        Node3Dn(const T1 t, const T1 xx, const T1 yy, const T1 zz, const size_t nt,
                const size_t i) :
        nThreads(nt),
        tt(new T1[nt]),
        x(xx), y(yy), z(zz),
        gridIndex(std::numeric_limits<T2>::max()),
        owners(std::vector<T2>(0)),
        slowness(0), primary(0)
        {
            for ( size_t n=0; n<nt; ++n ) {
                tt[n] = std::numeric_limits<T1>::max();
            }
            tt[i]=t;
        }

        /**
         * @brief Construct from a point, with a known traveltime on one thread.
         * @param t traveltime to store.
         * @param s point supplying the coordinates.
         * @param nt number of threads.
         * @param i thread number the traveltime belongs to.
         * @pre @p i is less than @p nt.
         */
        Node3Dn(const T1 t, const sxyz<T1>& s, const size_t nt,
                const size_t i) :
        nThreads(nt),
        tt(new T1[nt]),
        x(s.x), y(s.y), z(s.z),
        gridIndex(std::numeric_limits<T2>::max()),
        owners(std::vector<T2>(0)),
        slowness(0), primary(0)
        {
            for ( size_t n=0; n<nt; ++n ) {
                tt[n] = std::numeric_limits<T1>::max();
            }
            tt[i]=t;
        }

        /**
         * @brief Copy constructor; deep-copies the traveltime array.
         * @param node node to copy.
         */
        Node3Dn(const Node3Dn<T1,T2>& node) :
        nThreads(node.nThreads),
        tt(0),
        x(node.x), y(node.y), z(node.z),
        gridIndex(node.gridIndex),
        owners(node.owners),
        slowness(node.slowness),
        primary(node.primary)
        {
            tt = new T1[nThreads];

            for ( size_t n=0; n<nThreads; ++n ) {
                tt[n] = node.tt[n];
            }
        }

        /// Destructor; frees the traveltime array.
        virtual ~Node3Dn() {
            delete [] tt;
        }

        /**
         * @brief Reset one thread's traveltime to "not reached".
         * @param n thread whose slot is reset.
         */
        // Sets the vectors to the right size of threads and initialize it
        void reinit(const size_t n) {
            tt[n] = std::numeric_limits<T1>::max();
        }

        /// @return Number of traveltime slots this node was built with.
        const size_t getNThreads() const { return nThreads; }

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

        /**
         * @brief Mark this node as primary or secondary.
         * @param o true for a primary (corner) node, false for a secondary one.
         */
        void setPrimary(const bool o) { primary = o; }

        /// @return Index of this node in the grid's node list.
        T2 getGridIndex() const { return gridIndex; }
        /// @param index new index of this node in the grid's node list.
        void setGridIndex(const T2 index) { gridIndex = index; }

        /**
         * @brief Slowness at this node.
         * @return Slowness, i.e. reciprocal velocity.
         * @note This accessor is what distinguishes ttcr::Node3Dn from
         *       ttcr::Node3Dc, whose slowness lives on the cells instead.
         */
        T1 getNodeSlowness() const { return slowness; }
        /**
         * @brief Set the slowness at this node.
         * @param s slowness value.
         * @warning Unvalidated: a zero or negative slowness is stored as given
         *          and will produce infinite or nonsensical traveltimes.
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
        T1 getDistance( const Node3Dn<T1,T2>& node ) const {
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
         *       ttcr::input_parameters::translateOrigin is set. Mutating,
         *       unlike the free @c operator- overloads.
         */
        Node3Dn<T1, T2>& operator-=(const sxyz<T1>& node) {
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
         *          for the parent arrays of ttcr::Node3Dnsp, which this class
         *          does not have. The @c 4*sizeof(T1) term does correctly cover
         *          x, y, z and slowness. @ref primary is unaccounted for.
         */
        size_t getSize() const {
            return 2*sizeof(size_t) + nThreads*sizeof(T1) + 4*sizeof(T1) +
            (1+2*nThreads)*sizeof(T2) + owners.size() * sizeof(T2);
        }

        /// @return Always 3 — the spatial dimension of this node.
        int getDimension() const { return 3; }

        /// @return True for a primary (corner) node, false for a secondary one.
        const bool isPrimary() const { return primary; }

    protected:
        size_t nThreads;                ///< Number of traveltime slots, i.e. concurrent sources.
        T1 *tt;                         ///< travel time for the multiple source points; owned
        T1 x;                           ///< x coordinate [km]
        T1 y;							///< y coordinate [km]
        T1 z;                           ///< z coordinate [km]
        T2 gridIndex;                   ///< index of this node in the list of the grid
        std::vector<T2> owners;         ///< indices of cells touching the node
        /// Slowness at the node [s/km].
        /// @note The slowness is read by every @c rn / @c un grid class.
        T1 slowness;					// slowness at the node [s/km]
        bool primary;                   ///< True for a primary (corner) node, false for a secondary one.
    };

    /**
     * @brief Add two nodes' positions.
     * @return The sum, as a plain point — arithmetic on nodes yields
     *         @ref sxyz, since a sum has no traveltime or grid identity.
     */
    template<typename T1, typename T2>
    sxyz<T1> operator+(const Node3Dn<T1,T2>& lhs, const Node3Dn<T1,T2>& rhs) {
        return sxyz<T1>( lhs.getX()+rhs.getX(), lhs.getY()+rhs.getY(), lhs.getZ()+rhs.getZ() );
    }

    /// @brief Vector from @p rhs to @p lhs. @return The difference, as a point.
    template<typename T1, typename T2>
    sxyz<T1> operator-(const Node3Dn<T1,T2>& lhs, const Node3Dn<T1,T2>& rhs) {
        return sxyz<T1>( lhs.getX()-rhs.getX(), lhs.getY()-rhs.getY(), lhs.getZ()-rhs.getZ() );
    }

    /// @brief Vector from a node to a point. @return The difference, as a point.
    template<typename T1, typename T2>
    sxyz<T1> operator-(const sxyz<T1>& lhs, const Node3Dn<T1,T2>& rhs) {
        return sxyz<T1>( lhs.x-rhs.getX(), lhs.y-rhs.getY(), lhs.z-rhs.getZ() );
    }

    /// @brief Vector from a point to a node. @return The difference, as a point.
    template<typename T1, typename T2>
    sxyz<T1> operator-(const Node3Dn<T1,T2>& lhs, const sxyz<T1>& rhs) {
        return sxyz<T1>( lhs.getX()-rhs.x, lhs.getY()-rhs.y, lhs.getZ()-rhs.z );
    }

    /**
     * @brief Write a node's coordinates to a stream as "x y z".
     * @param os stream to write to.
     * @param n  node to write.
     * @return @p os, for chaining.
     */
    template<typename T1, typename T2>
    std::ostream& operator<< (std::ostream& os, const Node3Dn<T1, T2> &n) {
        os << n.getX() << ' ' << n.getY() << ' ' << n.getZ();
        return os;
    }
}

#endif

//
//  Node2Dn.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-04-15.
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
 * @file Node2Dn.h
 * @brief Node of a 2-D grid whose slowness is defined at the nodes.
 *
 * Declares ttcr::Node2Dn, the counterpart of ttcr::Node2Dc for grids that store
 * slowness **per node** rather than per cell — the @c rn and @c un families
 * (Grid2Drnfs, Grid2Drndsp, Grid2Dunfs, Grid2Dunfm, Grid2Dundsp and the OpenCL
 * variants), plus the node-based rectilinear Grid2Drcfs.
 *
 * It differs from ttcr::Node2Dc in exactly one respect: it carries a
 * @ref ttcr::Node2Dn::slowness member, reachable through
 * @ref ttcr::Node2Dn::getNodeSlowness. Everything else — the per-thread
 * traveltime array, the owners list, the primary flag — is the same. See
 * @ref n2dc_family for how the six Node2D classes relate, and
 * @ref n2dc_threads for the per-thread traveltime convention.
 *
 * @note Uses @ref ttcr::sxz and @ref ttcr::small but includes only Node.h; it
 *       compiles because every includer pulls in ttcr_t.h first. Unlike
 *       Node2Dc.h, which includes ttcr_t.h explicitly.
 *
 * @sa Node.h, Node2Dc.h, Node2Dnsp.h, Node2Dnd.h, Grid2Dun.h
 */

#ifndef ttcr_Node2Dn_h
#define ttcr_Node2Dn_h

#include <cmath>
#include <limits>

#include "Node.h"

namespace ttcr {

    /**
     * @brief Node of a 2-D grid with node-based slowness.
     *
     * Lives in the x-z plane; @ref getY always returns 0 so the class still
     * satisfies the 3-D ttcr::Node interface.
     *
     * @tparam T1 floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2 integer type of grid and cell indices, normally @c uint32_t.
     *
     * @warning Owns a raw @c new[] array and declares a copy constructor and a
     *          destructor but **no copy-assignment operator** — see the same
     *          warning on ttcr::Node2Dc. Copy-construct, never assign.
     */
    template<typename T1, typename T2>
    class Node2Dn : public Node<T1> {
    public:
        /**
         * @brief Construct an unpositioned node with unset traveltimes.
         * @param nt number of threads, i.e. the size of the traveltime array.
         * @post Coordinates and slowness are 0, the grid index is @c T2's
         *       maximum, and every traveltime is @c T1's maximum — the "not yet
         *       reached" marker.
         */
        Node2Dn(const size_t nt=1) :
        nThreads(nt), tt(0),
        x(0.0), z(0.0), slowness(0.0),
        gridIndex(std::numeric_limits<T2>::max()),
        owners(0), primary(false)
        {
            tt = new T1[nt];

            for ( size_t n=0; n<nt; ++n ) {
                tt[n] = std::numeric_limits<T1>::max();
            }
        }

        /**
         * @brief Construct a positioned node with a known traveltime on one thread.
         * @param t  traveltime to store.
         * @param _x x coordinate.
         * @param _z z coordinate.
         * @param nt number of threads.
         * @param i  thread number the traveltime belongs to.
         * @pre @p i is less than @p nt; the write is unchecked.
         * @note Slowness is left at 0 and must be set separately with
         *       @ref setNodeSlowness.
         */
        Node2Dn(const T1 t, const T1 _x, const T1 _z, const size_t nt, const size_t i) :
        nThreads(nt), tt(0),
        x(_x), z(_z), slowness(0.0),
        gridIndex(std::numeric_limits<T2>::max()),
        owners(std::vector<T2>(0)), primary(false)
        {
            tt = new T1[nt];

            for ( size_t n=0; n<nt; ++n ) {
                tt[n] = std::numeric_limits<T1>::max();
            }
            tt[i] = t;
        }

        /**
         * @brief Construct from a point, with a known traveltime on one thread.
         * @param t traveltime to store.
         * @param s point supplying the coordinates.
         * @param nt number of threads.
         * @param i thread number the traveltime belongs to.
         * @pre @p i is less than @p nt.
         */
        Node2Dn(const T1 t, const sxz<T1>& s, const size_t nt, const size_t i) :
        nThreads(nt), tt(0),
        x(s.x), z(s.z), slowness(0.0),
        gridIndex(std::numeric_limits<T2>::max()),
        owners(std::vector<T2>(0)), primary(false)
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
        Node2Dn(const Node2Dn<T1,T2>& node) :
        nThreads(node.nThreads), tt(0),
        x(node.x), z(node.z), slowness(node.slowness),
        gridIndex(node.gridIndex),
        owners(node.owners), primary(node.primary)
        {
            tt = new T1[nThreads];

            for ( size_t n=0; n<nThreads; ++n ) {
                tt[n] = node.tt[n];
            }
        }


        /// Destructor; frees the traveltime array.
        virtual ~Node2Dn() {
            delete [] tt;
        }

        /**
         * @brief Reset one thread's traveltime to "not reached".
         * @param thread_no thread whose slot is reset.
         * @note Called before each new source is propagated; other threads keep
         *       their solutions.
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
            x=xx; z=zz; gridIndex = index;  }

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

        /**
         * @brief Slowness at this node.
         * @return Slowness, i.e. reciprocal velocity.
         * @note This accessor is what distinguishes ttcr::Node2Dn from
         *       ttcr::Node2Dc, whose slowness lives on the cells instead.
         */
        T1 getNodeSlowness() const { return slowness; }
        /**
         * @brief Set the slowness at this node.
         * @param s slowness value.
         * @warning Unvalidated: a zero or negative slowness is stored as given
         *          and will produce infinite or nonsensical traveltimes.
         */
        void setNodeSlowness(const T1 s) { slowness = s; }

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
        T1 getDistance( const Node2Dn<T1,T2>& node ) const {
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

        /**
         * @brief Mark this node as primary or secondary.
         * @param p true for a primary (corner) node, false for a secondary one.
         * @note Unlike ttcr::Node2Dc::setPrimary this has no default argument.
         */
        void setPrimary(const bool p) { primary = p; }
        /// @return True for a primary (corner) node, false for a secondary one.
        const bool isPrimary() const { return primary; }

    protected:
        size_t nThreads;               ///< Number of traveltime slots, i.e. concurrent sources.
        T1 *tt;                        ///< travel time, one slot per thread; owned, @c new[] allocated
        T1 x;                          ///< x coordinate
        T1 z;                          ///< z coordinate
        T1 slowness;                   ///< Slowness at this node; the member ttcr::Node2Dc lacks.
        T2 gridIndex;                  ///< index of this node in the list of the grid
        std::vector<T2> owners;        ///< indices of cells touching the node
        bool primary;                  ///< True for a primary (corner) node, false for a secondary one.
    };

    /**
     * @brief Add two nodes' positions.
     * @return The sum, as a plain point — arithmetic on nodes yields
     *         @ref sxz, since a sum has no traveltime or grid identity.
     */
    template <typename T1, typename T2>
    sxz<T1> operator+(const Node2Dn<T1,T2>& lhs, const Node2Dn<T1,T2>& rhs) {
        return sxz<T1>(lhs.getX()+rhs.getX(), lhs.getZ()+rhs.getZ());
    }
    /// @brief Vector from @p rhs to @p lhs. @return The difference, as a point.
    template <typename T1, typename T2>
    sxz<T1> operator-(const Node2Dn<T1,T2>& lhs, const Node2Dn<T1,T2>& rhs) {
        return sxz<T1>(lhs.getX()-rhs.getX(), lhs.getZ()-rhs.getZ());
    }
    /// @brief Vector from a node to a point. @return The difference, as a point.
    template <typename T1, typename T2>
    sxz<T1> operator-(const sxz<T1>& lhs, const Node2Dn<T1,T2>& rhs) {
        return sxz<T1>(lhs.x-rhs.getX(), lhs.z-rhs.getZ());
    }

    /**
     * @brief Write a node's coordinates to a stream as "x z".
     * @param os stream to write to.
     * @param s  node to write.
     * @return @p os, for chaining.
     * @note The @c struct elaborated-type-specifier in the signature is
     *       vestigial — ttcr::Node2Dn is a class.
     */
    template <typename T1, typename T2>
    std::ostream& operator<< (std::ostream& os, const struct Node2Dn<T1,T2>& s) {
        os << s.getX() << ' ' << s.getZ();
        return os;
    }
}
#endif

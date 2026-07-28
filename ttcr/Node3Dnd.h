//
//  Node3Dnd.h.h
//  ttcr
//
//  Created by Bernard Giroux on 18-10-20.
//  Copyright © 2018 Bernard Giroux. All rights reserved.
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
 * @file Node3Dnd.h
 * @brief Temporary "dynamic" node for the node-based 3-D shortest-path solvers.
 *
 * Declares ttcr::Node3Dnd, the node-slowness counterpart of ttcr::Node3Dcd,
 * used by the dynamic shortest-path grids Grid3Drndsp and Grid3Dundsp. Dynamic
 * nodes are created in separate threads, so the base class is constructed with
 * 1 as the input argument for @c nThreads.
 *
 * @sa Node3Dn.h, Node3Dcd.h, Node2Dnd.h, Grid3Drndsp.h, Grid3Dundsp.h
 */

#ifndef ttcr_Node3Dnd_h
#define ttcr_Node3Dnd_h

#include "Node3Dn.h"

// for temporary dynamic nodes
// dynamic nodes are created in separate threads, so base class is created with
// 1 as input argument for nThreads

namespace ttcr {

    /**
     * @brief Single-threaded temporary node with node-based slowness.
     *
     * Always constructs its base with a thread count of 1 and overrides the
     * traveltime accessors to ignore the thread number.
     *
     * @tparam T1 floating-point type of coordinates, slowness and traveltimes.
     * @tparam T2 integer type of grid and cell indices.
     *
     * @warning @ref getTT and @ref setTT **discard** their thread argument and
     *          always use slot 0 — see ttcr::Node3Dcd for the rationale and the
     *          consequences.
     */
    template<typename T1, typename T2>
    class Node3Dnd : public Node3Dn<T1,T2> {
    public:
        /// Construct an unpositioned dynamic node with a single traveltime slot.
        Node3Dnd() : Node3Dn<T1,T2>(1) {}

        /// Destructor.
        virtual ~Node3Dnd() {}

        /**
         * @brief Construct a positioned dynamic node with a known traveltime.
         * @param t  traveltime to store in the single slot.
         * @param xx x coordinate.
         * @param yy y coordinate.
         * @param zz z coordinate.
         * @param nt ignored — the base is always built with 1 slot.
         * @param i  ignored — the traveltime always goes to slot 0.
         * @note @p nt and @p i exist only to match ttcr::Node3Dn's signature.
         */
        Node3Dnd(const T1 t, const T1 xx, const T1 yy, const T1 zz, const size_t nt,
                 const size_t i) : Node3Dn<T1,T2>(t, xx, yy, zz, 1, 0) {}

        /**
         * @brief Construct from a point, with a known traveltime.
         * @param t traveltime to store in the single slot.
         * @param s point supplying the coordinates.
         * @param nt ignored — the base is always built with 1 slot.
         * @param i  ignored — the traveltime always goes to slot 0.
         */
        Node3Dnd(const T1 t, const sxyz<T1>& s, const size_t nt,
                 const size_t i) : Node3Dn<T1,T2>(t, s, 1, 0) {}

        /**
         * @brief Copy constructor.
         * @param node node to copy.
         * @note Carries over the grid index, owners, slowness and primary flag
         *       on top of the base construction. ttcr::Node3Dcd's copy
         *       constructor omits the primary flag; this one does not.
         */
        Node3Dnd(const Node3Dnd<T1,T2>& node) :
        Node3Dn<T1,T2>(node.tt[0], node.x, node.y, node.z, 1, 0)
        {
            this->gridIndex = node.gridIndex;
            this->owners = node.owners;
            this->slowness = node.slowness;
            this->primary = node.primary;
        }

        /// @param n ignored. @return The single stored traveltime.
        T1 getTT(const size_t n) const { return this->tt[0]; }
        /// @param t traveltime to store. @param n ignored.
        void setTT(const T1 t, const size_t n ) { this->tt[0] = t; }

    };

    /**
     * @brief Write a node's coordinates to a stream as "x y z".
     * @param os stream to write to.
     * @param n  node to write.
     * @return @p os, for chaining.
     */
    template<typename T1, typename T2>
    std::ostream& operator<< (std::ostream& os, const Node3Dnd<T1, T2> &n) {
        os << n.getX() << ' ' << n.getY() << ' ' << n.getZ();
        return os;
    }

    /// @brief Vector from @p rhs to @p lhs. @return The difference, as a point.
    template<typename T1, typename T2>
    sxyz<T1> operator-(const Node3Dnd<T1,T2>& lhs, const Node3Dnd<T1,T2>& rhs) {
        return sxyz<T1>( lhs.getX()-rhs.getX(), lhs.getY()-rhs.getY(), lhs.getZ()-rhs.getZ() );
    }

    /// @brief Vector from a node to a point. @return The difference, as a point.
    template<typename T1, typename T2>
    sxyz<T1> operator-(const sxyz<T1>& lhs, const Node3Dnd<T1,T2>& rhs) {
        return sxyz<T1>( lhs.x-rhs.getX(), lhs.y-rhs.getY(), lhs.z-rhs.getZ() );
    }

    /// @brief Vector from a point to a node. @return The difference, as a point.
    template<typename T1, typename T2>
    sxyz<T1> operator-(const Node3Dnd<T1,T2>& lhs, const sxyz<T1>& rhs) {
        return sxyz<T1>( lhs.getX()-rhs.x, lhs.getY()-rhs.y, lhs.getZ()-rhs.z );
    }

}

#endif /* Node3Dnd_h_h */

//
//  Node2Dcd.h
//  ttcr
//
//  Created by Bernard Giroux on 2021-02-23.
//  Copyright © 2021 Bernard Giroux. All rights reserved.
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
 * @file Node2Dcd.h
 * @brief Temporary "dynamic" node for the cell-based 2-D shortest-path solvers.
 *
 * Declares ttcr::Node2Dcd, a lightweight ttcr::Node2Dc used by the dynamic
 * shortest-path grids (Grid2Drcdsp, Grid2Ducdsp). Those solvers insert extra
 * nodes near the source at run time to refine the traveltime there; such nodes
 * are created, used and discarded within a single worker thread, so unlike the
 * permanent grid nodes they need only **one** traveltime slot rather than one
 * per thread.
 *
 * @sa Node2Dc.h, Node2Dnd.h, Grid2Drcdsp.h, Grid2Ducdsp.h
 */

#ifndef ttcr_Node2Dcd_h
#define ttcr_Node2Dcd_h

#include "Node2Dc.h"

// for temporary dynamic nodes

namespace ttcr {

    /**
     * @brief Single-threaded temporary node with cell-based slowness.
     *
     * Always constructs its base with a thread count of 1 and overrides the
     * traveltime accessors to ignore the thread number, so a dynamic node can be
     * passed to code written against the multi-threaded node interface while
     * carrying only one traveltime.
     *
     * @tparam T1 floating-point type of coordinates and traveltimes.
     * @tparam T2 integer type of grid and cell indices.
     *
     * @warning @ref getTT and @ref setTT **discard** their thread argument and
     *          always use slot 0. That is deliberate — these nodes are
     *          thread-local — but it means an instance must never be shared
     *          between threads, and calling with a nonzero thread number
     *          silently reads or writes slot 0 instead of failing.
     */
    template<typename T1, typename T2>
    class Node2Dcd : public Node2Dc<T1,T2> {
    public:
        /// Construct an unpositioned dynamic node with a single traveltime slot.
        Node2Dcd() : Node2Dc<T1,T2>(1) {}

        /**
         * @brief Construct a positioned dynamic node with a known traveltime.
         * @param t  traveltime to store in the single slot.
         * @param xx x coordinate.
         * @param zz z coordinate.
         * @param nt ignored — the base is always built with 1 slot.
         * @param i  ignored — the traveltime always goes to slot 0.
         * @note @p nt and @p i exist only so the signature matches
         *       ttcr::Node2Dc's, letting the two be used interchangeably in
         *       templates. They have no effect.
         */
        Node2Dcd(const T1 t, const T1 xx, const T1 zz, const size_t nt,
                 const size_t i) : Node2Dc<T1,T2>(t, xx, zz, 1, 0) {}

        /**
         * @brief Copy constructor.
         * @param node node to copy.
         * @note Copies the grid index and owners explicitly on top of the base
         *       construction. It does **not** copy the @c primary flag, so a
         *       copy is always non-primary regardless of the original — unlike
         *       ttcr::Node2Dnd, which does copy it.
         */
        Node2Dcd(const Node2Dcd<T1,T2>& node) :
        Node2Dc<T1,T2>(node.tt[0], node.x, node.z, 1, 0)
        {
            this->gridIndex = node.gridIndex;
            this->owners = node.owners;
        }

        /// @param n ignored. @return The single stored traveltime.
        T1 getTT(const size_t n) const { return this->tt[0]; }
        /// @param t traveltime to store. @param n ignored.
        void setTT(const T1 t, const size_t n ) { this->tt[0] = t; }
    };

    /**
     * @brief Write a node's coordinates to a stream as "x z".
     * @param os stream to write to.
     * @param n  node to write.
     * @return @p os, for chaining.
     */
    template<typename T1, typename T2>
    std::ostream& operator<< (std::ostream& os, const Node2Dcd<T1, T2> &n) {
        os << n.getX() << ' ' << n.getZ();
        return os;
    }

    /// @brief Vector from @p rhs to @p lhs. @return The difference, as a point.
    template<typename T1, typename T2>
    sxz<T1> operator-(const Node2Dcd<T1,T2>& lhs, const Node2Dcd<T1,T2>& rhs) {
        return sxz<T1>( lhs.getX()-rhs.getX(), lhs.getZ()-rhs.getZ() );
    }

    /// @brief Vector from a node to a point. @return The difference, as a point.
    template<typename T1, typename T2>
    sxz<T1> operator-(const sxz<T1>& lhs, const Node2Dcd<T1,T2>& rhs) {
        return sxz<T1>( lhs.x-rhs.getX(), lhs.z-rhs.getZ() );
    }

    /// @brief Vector from a point to a node. @return The difference, as a point.
    template<typename T1, typename T2>
    sxz<T1> operator-(const Node2Dcd<T1,T2>& lhs, const sxz<T1>& rhs) {
        return sxz<T1>( lhs.getX()-rhs.x, lhs.getZ()-rhs.z );
    }
}
#endif /* Node2Dcd_h */

//
//  Node.h
//  ttcr
//
//  Created by Bernard Giroux on 12-08-14.
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
 * @file Node.h
 * @brief Minimal abstract node interface and a traveltime-based comparator.
 *
 * Declares @ref ttcr::Node, the small polymorphic base exposing the coordinate
 * and per-thread traveltime accessors that grid/mesh algorithms rely on, and
 * @ref ttcr::CompareNodePtr, a comparator used to order node pointers by
 * traveltime in a priority queue during the shortest-path solve.
 */

#ifndef ttcr_Node_h
#define ttcr_Node_h

#include <cstddef>

namespace ttcr {

    /**
     * @brief Abstract interface for a grid/mesh node.
     *
     * Defines the read-only accessors the traveltime solvers use, without
     * committing to a storage layout; concrete node classes supply the
     * implementation.
     *
     * @tparam T floating-point type for coordinates and traveltimes.
     */
    template<typename T>
    class Node {
    public:
        /// @param n thread number. @return Traveltime stored for thread @p n.
        virtual T getTT(const size_t n) const = 0;
        /// @return Spatial dimension of the node (2 or 3).
        virtual int getDimension() const = 0;
        virtual T getX() const = 0;  ///< @return x coordinate of the node.
        virtual T getY() const = 0;  ///< @return y coordinate of the node.
        virtual T getZ() const = 0;  ///< @return z coordinate of the node.
        /// @return True if this is a primary (corner) node, false for a secondary node.
        virtual const bool isPrimary() const = 0;
    };


    /**
     * @brief Comparator ordering node pointers by traveltime for a min-priority queue.
     *
     * Compares the traveltimes of two nodes for a chosen thread. The comparison
     * is reversed (`>`) so that a `std::priority_queue` — which pops its largest
     * element — yields the node with the *smallest* traveltime first, as the
     * shortest-path propagation requires.
     *
     * @tparam T floating-point type for traveltimes.
     */
    template<typename T>
    class CompareNodePtr {
        // Overloaded operator for the priority queue, compare the "n"th traveltimes of two nodes.
    private:
        size_t n;  ///< Thread number whose traveltimes are compared.
    public:
        /// @param nn thread number to compare traveltimes for.
        CompareNodePtr(const size_t nn) : n(nn) {}
        /**
         * @param n1 first node.
         * @param n2 second node.
         * @return True if @p n1 has a larger traveltime than @p n2 (min-heap ordering).
         */
        bool operator()(const Node<T>* n1, const Node<T>* n2) const {
            //  The priority_queue must return the minimum time!!!
            return n1->getTT(n) > n2->getTT(n);
        }
    };

}

#endif

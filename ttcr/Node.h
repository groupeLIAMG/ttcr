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

#ifndef ttcr_Node_h
#define ttcr_Node_h

namespace ttcr {

    template<typename T>
    class Node {
    public:
        virtual T getTT(const size_t n) const = 0;
        virtual int getDimension() const = 0;
        virtual T getX() const = 0;
        virtual T getY() const = 0;
        virtual T getZ() const = 0;
        virtual const bool isPrimary() const = 0;
    };


    template<typename T>
    class CompareNodePtr {
        // Overloaded operator for the priority queue, compare the "n"th traveltimes of two nodes.
    private:
        size_t n;
    public:
        CompareNodePtr(const size_t nn) : n(nn) {}
        bool operator()(const Node<T>* n1, const Node<T>* n2) const {
            //  The priority_queue must return the minimum time!!!
            return n1->getTT(n) > n2->getTT(n);
        }
    };

}

#endif

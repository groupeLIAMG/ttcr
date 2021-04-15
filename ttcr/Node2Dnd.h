//
//  Node2Dnd.h
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

#ifndef ttcr_Node2Dnd_h
#define ttcr_Node2Dnd_h

#include "Node2Dn.h"

// for temporary dynamic nodes
// dynamic nodes are created in separate threads, so base class is created with
// 1 as input argument for nThreads

namespace ttcr {

    template<typename T1, typename T2>
    class Node2Dnd : public Node2Dn<T1,T2> {
    public:
        Node2Dnd() : Node2Dn<T1,T2>(1) {}

        virtual ~Node2Dnd() {}

        Node2Dnd(const T1 t, const T1 xx, const T1 zz, const size_t nt,
                 const size_t i) : Node2Dn<T1,T2>(t, xx, zz, 1, 0) {}

        Node2Dnd(const T1 t, const sxz<T1>& s, const size_t nt,
                 const size_t i) : Node2Dn<T1,T2>(t, s, 1, 0) {}

        Node2Dnd(const Node2Dnd<T1,T2>& node) :
        Node2Dn<T1,T2>(node.tt[0], node.x, node.z, 1, 0)
        {
            this->gridIndex = node.gridIndex;
            this->owners = node.owners;
            this->slowness = node.slowness;
            this->primary = node.primary;
        }

        T1 getTT(const size_t n) const { return this->tt[0]; }
        void setTT(const T1 t, const size_t n ) { this->tt[0] = t; }

    };

    template<typename T1, typename T2>
    std::ostream& operator<< (std::ostream& os, const Node2Dnd<T1, T2> &n) {
        os << n.getX() << ' ' << n.getZ();
        return os;
    }

    template<typename T1, typename T2>
    sxz<T1> operator-(const Node2Dnd<T1,T2>& lhs, const Node2Dnd<T1,T2>& rhs) {
        return sxz<T1>( lhs.getX()-rhs.getX(), lhs.getZ()-rhs.getZ() );
    }

    template<typename T1, typename T2>
    sxz<T1> operator-(const sxz<T1>& lhs, const Node2Dnd<T1,T2>& rhs) {
        return sxz<T1>( lhs.x-rhs.getX(), lhs.z-rhs.getZ() );
    }

    template<typename T1, typename T2>
    sxz<T1> operator-(const Node2Dnd<T1,T2>& lhs, const sxz<T1>& rhs) {
        return sxz<T1>( lhs.getX()-rhs.x, lhs.getZ()-rhs.z );
    }

}

#endif /* Node2Dnd_h */

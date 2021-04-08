//
//  Node2Dc.h
// ttRaisCourbes
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

#ifndef ttcr_NodeDc_h
#define ttcr_NodeDc_h

#include <cmath>
#include <limits>

#include "ttcr_t.h"
#include "Node.h"

namespace ttcr {

    template<typename T1, typename T2>
    class Node2Dc : public Node<T1> {
    public:
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


        virtual ~Node2Dc() {
            delete [] tt;
        }

        void reinit(const size_t thread_no) { //=0) {
            tt[thread_no] = std::numeric_limits<T1>::max();
        }

        T1 getTT(const size_t i) const { return tt[i]; }
        void setTT(const T1 t, const size_t i) { tt[i] = t; }

        void setXZindex(const T1 xx, const T1 zz, const T2 index) {
            x=xx; z=zz; gridIndex = index; }

        template<typename SXZ>
        void setXYZindex(const SXZ& s, const T2 index) {
            x=s.x; z=s.z; gridIndex = index;  }

        T1 getX() const {
            return x;
        }
        void setX(const T1 xx) { x = xx; }

        T1 getY() const { return 0.0; }

        T1 getZ() const { return z; }
        void setZ(const T1 zz) { z = zz; }

        T2 getGridIndex() const { return gridIndex; }
        void setGridIndex(const T2 index) { gridIndex = index; }

        void pushOwner(const T2 o) { owners.push_back(o); }
        const std::vector<T2>& getOwners() const { return owners; }

        T1 getDistance( const Node2Dc<T1,T2>& node ) const {
            return sqrt( (x-node.x)*(x-node.x) + (z-node.z)*(z-node.z) );
        }

        T1 getDistance( const sxz<T1>& node ) const {
            return sqrt( (x-node.x)*(x-node.x) + (z-node.z)*(z-node.z) );
        }

        T1 getDistanceX( const sxz<T1>& node ) const {
            return std::abs( x-node.x );
        }

        T1 getDistanceZ( const sxz<T1>& node ) const {
            return std::abs( z-node.z );
        }

        // operator to test if same location
        bool operator==( const sxz<T1>& node ) const {
            return std::abs(x-node.x)<small && std::abs(z-node.z)<small;
        }

        size_t getSize() const {
            return sizeof(size_t) + nThreads*sizeof(T1) + 2*sizeof(T1) +
            (1+2*nThreads)*sizeof(T2) + owners.size() * sizeof(T2);
        }

        int getDimension() const { return 2; }

        void setPrimary(const bool p=true) {
            primary = p;
        }
        const bool isPrimary() const { return primary; }

    protected:
        size_t nThreads;
        T1 x;                          // x coordinate
        T1 z;                          // z coordinate
        T2 gridIndex;                  // index of this node in the list of the grid
        T1 *tt;                        // travel time
        std::vector<T2> owners;        // indices of cells touching the node
        bool primary;
    };

    template<typename T1, typename T2>
    std::ostream& operator<< (std::ostream& os, const Node2Dc<T1, T2> &n) {
        os << n.getX() << ' ' << n.getZ();
        return os;
    }

    template<typename T1, typename T2>
    sxz<T1> operator+(const Node2Dc<T1,T2>& lhs, const Node2Dc<T1,T2>& rhs) {
        return sxz<T1>( lhs.getX()+rhs.getX(), lhs.getZ()+rhs.getZ() );
    }

    template<typename T1, typename T2>
    sxz<T1> operator-(const Node2Dc<T1,T2>& lhs, const Node2Dc<T1,T2>& rhs) {
        return sxz<T1>( lhs.getX()-rhs.getX(), lhs.getZ()-rhs.getZ() );
    }

    template<typename T1, typename T2>
    sxz<T1> operator-(const sxz<T1>& lhs, const Node2Dc<T1,T2>& rhs) {
        return sxz<T1>( lhs.x-rhs.getX(), lhs.z-rhs.getZ() );
    }

    template<typename T1, typename T2>
    sxz<T1> operator-(const Node2Dc<T1,T2>& lhs, const sxz<T1>& rhs) {
        return sxz<T1>( lhs.getX()-rhs.x, lhs.getZ()-rhs.z);
    }

}

#endif

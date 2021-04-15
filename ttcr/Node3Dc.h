//
//  Node3Dc.h
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

#ifndef ttcr_Node3Dc_h
#define ttcr_Node3Dc_h

#include <cmath>
#include <limits>

#include "Node.h"

namespace ttcr {

    template<typename T1, typename T2>
    class Node3Dc : public Node<T1> {
    public:
        Node3Dc(const size_t nt=1) :
        nThreads(nt),
        x(0.0f), y(0.0f), z(0.0f),
        gridIndex(std::numeric_limits<T2>::max()),
        tt(nullptr),
        owners(0),
        primary(true)
        {
            tt = new T1[nt];

            for ( size_t n=0; n<nt; ++n ) {
                tt[n] = std::numeric_limits<T1>::max();
            }
        }

        Node3Dc(const T1 xx, const T1 yy, const T1 zz, const T2 index,
                const size_t nt) :
        nThreads(nt),
        x(xx), y(yy), z(zz),
        gridIndex(index),
        tt(nullptr),
        owners(0),
        primary(true)
        {
            tt = new T1[nt];

            for ( size_t n=0; n<nt; ++n ) {
                tt[n] = std::numeric_limits<T1>::max();
            }
        }

        Node3Dc(const T1 t, const T1 xx, const T1 yy, const T1 zz, const size_t nt,
                const size_t i) :
        nThreads(nt),
        x(xx), y(yy), z(zz),
        gridIndex(std::numeric_limits<T2>::max()),
        tt(nullptr),
        owners(0),
        primary(true)
        {
            tt = new T1[nt];

            for ( size_t n=0; n<nt; ++n ) {
                tt[n] = std::numeric_limits<T1>::max();
            }
            tt[i] = t;
        }

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

        virtual ~Node3Dc() {
            delete [] tt;
        }

        // Sets the vectors to the right size of threads and initialize it
        void reinit(const size_t n) {
            tt[n] = std::numeric_limits<T1>::max();
        }

        T1 getTT(const size_t n) const { return tt[n]; }
        void setTT(const T1 t, const size_t n ) { tt[n] = t; }

        void setXYZindex(const T1 xx, const T1 yy, const T1 zz, const T2 index) {
            x=xx; y=yy; z=zz; gridIndex = index;  }

        T1 getX() const { return x; }
        void setX(const T1 xx) { x = xx; }

        T1 getY() const { return y; }
        void setY(const T1 yy) { y = yy; }

        T1 getZ() const { return z; }
        void setZ(const T1 zz) { z = zz; }

        T2 getGridIndex() const { return gridIndex; }
        void setGridIndex(const T2 index) { gridIndex = index; }

        void pushOwner(const T2 o) { owners.push_back(o); }
        const std::vector<T2>& getOwners() const { return owners; }

        T1 getDistance( const Node3Dc<T1,T2>& node ) const {
            return sqrt( (x-node.x)*(x-node.x) + (y-node.y)*(y-node.y) + (z-node.z)*(z-node.z) );
        }

        T1 getDistance( const sxyz<T1>& node ) const {
            return sqrt( (x-node.x)*(x-node.x) + (y-node.y)*(y-node.y) + (z-node.z)*(z-node.z) );
        }

        // operator to test if same location
        bool operator==( const sxyz<T1>& node ) const {
            return std::abs(x-node.x)<small && std::abs(y-node.y)<small && std::abs(z-node.z)<small;
        }

        Node3Dc<T1, T2>& operator-=(const sxyz<T1>& node) {
            this->x -= node.x;
            this->y -= node.y;
            this->z -= node.z;
            return *this;
        }

        size_t getSize() const {
            return sizeof(size_t) + nThreads*sizeof(T1) + 3*sizeof(T1) +
            (1+2*nThreads)*sizeof(T2) + owners.size() * sizeof(T2);
        }

        int getDimension() const { return 3; }

        void setPrimary(const bool p=true) {
            primary = p;
        }
        const bool isPrimary() const { return primary; }

    protected:
        size_t nThreads;
        T1 x;                       // x coordinate [km]
        T1 y;						// y coordinate [km]
        T1 z;                       // z coordinate [km]
        T2 gridIndex;               // index of this node in the list of the grid
        T1 *tt;                     // travel time for the multiple source points
        std::vector<T2> owners;     // indices of cells touching the node
        bool primary;
    };

    template<typename T1, typename T2>
    sxyz<T1> operator+(const Node3Dc<T1,T2>& lhs, const Node3Dc<T1,T2>& rhs) {
        return sxyz<T1>( lhs.getX()+rhs.getX(), lhs.getY()+rhs.getY(), lhs.getZ()+rhs.getZ() );
    }

    template<typename T1, typename T2>
    sxyz<T1> operator-(const Node3Dc<T1,T2>& lhs, const Node3Dc<T1,T2>& rhs) {
        return sxyz<T1>( lhs.getX()-rhs.getX(), lhs.getY()-rhs.getY(), lhs.getZ()-rhs.getZ() );
    }

    template<typename T1, typename T2>
    sxyz<T1> operator-(const sxyz<T1>& lhs, const Node3Dc<T1,T2>& rhs) {
        return sxyz<T1>( lhs.x-rhs.getX(), lhs.y-rhs.getY(), lhs.z-rhs.getZ() );
    }

    template<typename T1, typename T2>
    sxyz<T1> operator-(const Node3Dc<T1,T2>& lhs, const sxyz<T1>& rhs) {
        return sxyz<T1>( lhs.getX()-rhs.x, lhs.getY()-rhs.y, lhs.getZ()-rhs.z);
    }

    template<typename T1, typename T2>
    std::ostream& operator<< (std::ostream& os, const Node3Dc<T1, T2> &n) {
        os << n.getX() << ' ' << n.getY() << ' ' << n.getZ();
        return os;
    }

}

#endif

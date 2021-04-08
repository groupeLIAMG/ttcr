//
//  Node2Dcsp.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-03-03.
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

#ifndef ttcr_Node2Dcsp_h
#define ttcr_Node2Dcsp_h

#include <cmath>
#include <limits>

#include "ttcr_t.h"
#include "Node.h"

namespace ttcr {

    template<typename T1, typename T2>
    class Node2Dcsp : public Node<T1> {
    public:
        Node2Dcsp(const size_t nt=1) :
        nThreads(nt),
        x(0.0), z(0.0),
        gridIndex(std::numeric_limits<T2>::max()),
        tt(nullptr),
        nodeParent(nullptr),
        cellParent(nullptr),
        owners(0),
        primary(false)
        {
            tt = new T1[nt];
            nodeParent = new T2[nt];
            cellParent = new T2[nt];

            for ( size_t n=0; n<nt; ++n ) {
                tt[n] = std::numeric_limits<T1>::max();
                nodeParent[n] = std::numeric_limits<T2>::max();
                cellParent[n] = std::numeric_limits<T2>::max();
            }
        }

        Node2Dcsp(const T1 xx, const T1 zz, const T2 index, const size_t nt=1) :
        nThreads(nt),
        x(xx), z(zz),
        gridIndex(index),
        tt(nullptr),
        nodeParent(nullptr),
        cellParent(nullptr),
        owners(0),
        primary(false)
        {
            tt = new T1[nt];
            nodeParent = new T2[nt];
            cellParent = new T2[nt];

            for ( size_t n=0; n<nt; ++n ) {
                tt[n] = std::numeric_limits<T1>::max();
                nodeParent[n] = std::numeric_limits<T2>::max();
                cellParent[n] = std::numeric_limits<T2>::max();
            }
        }


        Node2Dcsp(const T1 t, const T1 xx, const T1 zz, const size_t nt, const size_t i) :
        nThreads(nt),
        x(xx), z(zz),
        gridIndex(std::numeric_limits<T2>::max()),
        tt(nullptr),
        nodeParent(nullptr),
        cellParent(nullptr),
        owners(std::vector<T2>(0)),
        primary(false)
        {
            tt = new T1[nt];
            nodeParent = new T2[nt];
            cellParent = new T2[nt];

            for ( size_t n=0; n<nt; ++n ) {
                tt[n] = std::numeric_limits<T1>::max();
                nodeParent[n] = std::numeric_limits<T2>::max();
                cellParent[n] = std::numeric_limits<T2>::max();
            }
            tt[i] = t;
        }

        template<typename SXZ>
        Node2Dcsp(const T1 t, const SXZ &s, const size_t nt, const size_t i) :
        nThreads(nt),
        x(s.x), z(s.z),
        gridIndex(std::numeric_limits<T2>::max()),
        tt(nullptr),
        nodeParent(nullptr),
        cellParent(nullptr),
        owners(std::vector<T2>(0)),
        primary(false)
        {
            tt = new T1[nt];
            nodeParent = new T2[nt];
            cellParent = new T2[nt];

            for ( size_t n=0; n<nt; ++n ) {
                tt[n] = std::numeric_limits<T1>::max();
                nodeParent[n] = std::numeric_limits<T2>::max();
                cellParent[n] = std::numeric_limits<T2>::max();
            }
            tt[i] = t;
        }

        Node2Dcsp(const Node2Dcsp<T1,T2>& node) :
        nThreads(node.nThreads),
        x(node.x), z(node.z),
        gridIndex(node.gridIndex),
        tt(nullptr),
        nodeParent(nullptr),
        cellParent(nullptr),
        owners(node.owners),
        primary(node.primary)
        {
            tt = new T1[nThreads];
            nodeParent = new T2[nThreads];
            cellParent = new T2[nThreads];

            for ( size_t n=0; n<nThreads; ++n ) {
                tt[n] = node.tt[n];
                nodeParent[n] = node.nodeParent[n];
                cellParent[n] = node.cellParent[n];
            }
        }

        virtual ~Node2Dcsp() {
            delete [] tt;
            delete [] nodeParent;
            delete [] cellParent;
        }

        void reinit(const size_t thread_no) { //=0) {
            tt[thread_no] = std::numeric_limits<T1>::max();
            nodeParent[thread_no] = std::numeric_limits<T2>::max();
            cellParent[thread_no] = std::numeric_limits<T2>::max();
        }

        T1 getTT(const size_t i) const { return tt[i]; }
        void setTT(const T1 t, const size_t i) { tt[i] = t; }

        void setXZindex(const T1 xx, const T1 zz, const T2 index) {
            x=xx; z=zz; gridIndex = index;  }

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

        T2 getNodeParent(const size_t i) const { return nodeParent[i]; }
        void setnodeParent(const T2 index, const size_t i) { nodeParent[i] = index; }

        T2 getCellParent(const size_t i) const { return cellParent[i]; }
        void setCellParent(const T2 index, const size_t i) { cellParent[i] = index; }

        void pushOwner(const T2 o) { owners.push_back(o); }
        const std::vector<T2>& getOwners() const { return owners; }

        T1 getDistance( const Node2Dcsp<T1,T2>& node ) const {
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

        void setPrimary(const bool p) { primary = p; }
        const bool isPrimary() const { return primary; }

    private:
        size_t nThreads;
        T1 x;                          // x coordinate
        T1 z;                          // z coordinate
        T2 gridIndex;                  // index of this node in the list of the grid
        T1 *tt;                        // travel time
        T2 *nodeParent;                // index of parent node of the ray
        T2 *cellParent;                // index of cell traversed by the ray
        std::vector<T2> owners;        // indices of cells touching the node
        bool primary;

    };

    template<typename T1, typename T2>
    std::ostream& operator<< (std::ostream& os, const Node2Dcsp<T1, T2> &n) {
        os << n.getX() << ' ' << n.getZ();
        return os;
    }

    template <typename T1, typename T2>
    sxz<T1> operator-(const sxz<T1>& lhs, const Node2Dcsp<T1,T2>& rhs) {
        return sxz<T1>(lhs.x-rhs.getX(), lhs.z-rhs.getZ());
    }
}

#endif

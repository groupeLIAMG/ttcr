//
//  Node2Dnsp.h
//  ttcr
//
//  Created by Bernard Giroux on 2012-09-25.
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

#ifndef ttcr_Node2Dnsp_h
#define ttcr_Node2Dnsp_h

#include <cmath>
#include <limits>

#include "Node.h"

namespace ttcr {

    template<typename T1, typename T2>
    class Node2Dnsp : public Node<T1> {
    public:
        Node2Dnsp(const size_t nt=1) :
        nThreads(nt),
        tt(0),
        x(0.0), z(0.0), slowness(0.0),
        gridIndex(std::numeric_limits<T2>::max()),
        nodeParent(0),
        cellParent(0),
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


        Node2Dnsp(const T1 t, const sxz<T1>& s, const size_t nt, const size_t i) :
        nThreads(nt),
        tt(0),
        x(s.x), z(s.z), slowness(0.0),
        gridIndex(std::numeric_limits<T2>::max()),
        nodeParent(0),
        cellParent(0),
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

        Node2Dnsp(const Node2Dnsp<T1,T2>& node) :
        nThreads(node.nThreads),
        tt(0),
        x(node.x), z(node.z), slowness(node.slowness),
        gridIndex(node.gridIndex),
        nodeParent(0),
        cellParent(0),
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


        virtual ~Node2Dnsp() {
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

        T1 getNodeSlowness() const { return slowness; }
        void setNodeSlowness(const T1 s) { slowness = s; }

        T2 getGridIndex() const { return gridIndex; }
        void setGridIndex(const T2 index) { gridIndex = index; }

        T2 getNodeParent(const size_t i) const { return nodeParent[i]; }
        void setnodeParent(const T2 index, const size_t i) { nodeParent[i] = index; }

        T2 getCellParent(const size_t i) const { return cellParent[i]; }
        void setCellParent(const T2 index, const size_t i) { cellParent[i] = index; }

        void setPrimary( const bool o ) { primary = o; }

        void pushOwner(const T2 o) { owners.push_back(o); }
        const std::vector<T2>& getOwners() const { return owners; }

        T1 getDistance( const Node2Dnsp<T1,T2>& node ) const {
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

        int getDimension() const { return 2; }

        const bool isPrimary() const { return primary; }

    private:
        size_t nThreads;
        T1 *tt;                        // travel time
        T1 x;                          // x coordinate
        T1 z;                          // z coordinate
        T1 slowness;
        T2 gridIndex;                  // index of this node in the list of the grid
        T2 *nodeParent;                // index of parent node of the ray
        T2 *cellParent;                // index of cell traversed by the ray
        std::vector<T2> owners;        // indices of cells touching the node
        bool primary;				   // indicate the order of the node: 5= primary,

    };

    template<typename T1, typename T2>
    sxz<T1> operator+(const Node2Dnsp<T1,T2>& lhs, const Node2Dnsp<T1,T2>& rhs) {
        return sxz<T1>( lhs.getX()+rhs.getX(), lhs.getZ()+rhs.getZ() );
    }

    template <typename T1, typename T2>
    sxz<T1> operator-(const sxz<T1>& lhs, const Node2Dnsp<T1,T2>& rhs) {
        return sxz<T1>(lhs.x-rhs.getX(), lhs.z-rhs.getZ());
    }

    template <typename T1, typename T2>
    std::ostream& operator<< (std::ostream& os, const struct Node2Dnsp<T1,T2>& s) {
        os << s.getX() << ' ' << s.getZ();
        return os;
    }

}

#endif

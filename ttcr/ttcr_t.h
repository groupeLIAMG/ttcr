//
//  ttcr_t.h
//  ttcr_t
//
//  Created by Bernard Giroux on 08-07-01.
//
//

//
// Copyright (C) 2008 Bernard Giroux.
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

#ifndef __TTCR_T_H__
#define __TTCR_T_H__

#include <cmath>
#include <type_traits>
#include <vector>

namespace ttcr {
    
    const double small = 1.e-5;
    const double pi = 4.0*atan(1.0);
    
    const size_t iLength[4][3]={{0,1,2},{1,3,4},{2,3,5},{0,4,5}};
    const size_t iNodes[4][3] = {
        {0,1,2},  // (relative) indices of nodes of 1st triangle
        {1,2,3},  // (relative) indices of nodes of 2nd triangle
        {0,2,3},  // (relative) indices of nodes of 3rd triangle
        {0,1,3}   // (relative) indices of nodes of 4th triangle
    };
    
    //
    //
    //              1
    //            ,/|`\
    //          ,/  |  `\
    //        ,0    '.   `4
    //      ,/       1     `\
    //    ,/         |       `\
    //   0-----5-----'.--------3
    //    `\.         |      ,/
    //       `\.      |     3
    //          `2.   '. ,/
    //             `\. |/
    //                `2
    //
    //
    //  triangle 0:  0-1  1-2  2-0     (first occurence of edge underlined)
    //               ---  ---  ---
    //  triangle 1:  1-2  2-3  3-1
    //                    ---  ---
    //  triangle 2:  0-2  2-3  3-0
    //                         ---
    //  triangle 3:  0-1  1-3  3-0
    
    // for triangle "itri", indices of length for edge "iedge" are
    // iLength[itri][iedge] = {{0,1,2},{1,3,4},{2,3,5},{0,4,5}}
    
    
    template<typename T>
    struct sxz {
        T x;
        T z;
        
        sxz() : x(0), z(0) {}
        sxz(const T x_, const T z_) : x(x_), z(z_) {}
        
        bool operator==(const sxz<T>& s) const {
            //        return (fabs(x-a.x)<small && fabs(z-a.z)<small);
            return x==s.x && z==s.z;
        }
        
        sxz<T>& operator+=(const sxz<T>& s)
        {
            x += s.x;
            z += s.z;
            return *this;
        }
        
        sxz<T>& operator-=(const sxz<T>& s) {
            x -= s.x;
            z -= s.z;
            return *this;
        }
        
        sxz<T>& operator/=(const T s) {
            x /= s;
            z /= s;
            return *this;
        }
        
        sxz<T>& operator=(const T s) {
            x = s;
            z = s;
            return *this;
        }
        
        sxz<T> & operator=(const sxz<T>& s) {
            if (this != &s) { // protect against invalid self-assignment
                x = s.x;
                z = s.z;
            }
            return *this;
        }
        
        template<typename NODE>
        sxz<T>& operator=(const NODE& n) {
            x = n.getX();
            z = n.getZ();
            return *this;
        }
        
        sxz<T>& operator *=(const T value) {
            x *= value;
            z *= value;
            return *this;
        }
        
        bool operator<(const sxz<T>& s) const {
            return x<s.x && z<s.z;
        }
        
        T getDistance(const sxz<T>& s) const {
            return sqrt( (x-s.x)*(x-s.x) + (z-s.z)*(z-s.z) );
        }
        
        void normalize() {
            T n = sqrt( x*x + z*z );
            x /= n;
            z /= n;
        }
    };
    
    template<typename T>
    std::ostream& operator<< (std::ostream& os, const struct sxz<T> &s) {
        os << s.x << ' ' << s.z;
        return os;
    }
    
    template<typename T>
    std::istream& operator>> (std::istream& is, struct sxz<T> &s) {
        is >> s.x >> s.z;
        return is;
    }
    
    template <typename T>
    sxz<T> operator+(const sxz<T>& lhs, const sxz<T>& rhs)
    {
        return sxz<T>(lhs.x+rhs.x, lhs.z+rhs.z);
    }
    
    template <typename T>
    sxz<T> operator-(const sxz<T>& lhs, const sxz<T>& rhs)
    {
        return sxz<T>(lhs.x-rhs.x, lhs.z-rhs.z);
    }
    
    template <typename T>
    sxz<T> operator*(const T& scalar, const sxz<T>& v)
    {
        return sxz<T>(scalar*v.x, scalar*v.z);
    }
    
    template <typename T>
    sxz<T> operator/(const sxz<T>& v, const T& scalar)
    {
        return sxz<T>(v.x/scalar, v.z/scalar);
    }
    
    template<typename T>
    sxz<T> normalize(const sxz<T>& s) {
        sxz<T> s2 = s;
        s2.normalize();
        return s2;
    }
    
    
    
    template<typename T>
    struct sxyz {
        T x;
        T y;
        T z;
        
        sxyz() : x(0), y(0), z(0) {}
        sxyz(const T x_, const T y_, const T z_) : x(x_), y(y_), z(z_) {}
        
        bool operator==(const sxyz<T>& s) const {
            //        return (fabs(x-a.x)<small && fabs(z-a.z)<small);
            return x==s.x && y==s.y && z==s.z;
        }
        
        sxyz<T>& operator=(const sxyz<T>& s) {
            if (this != &s) { // protect against invalid self-assignment
                x = s.x;
                y = s.y;
                z = s.z;
            }
            return *this;
        }
        
        sxyz<T>& operator+=(const sxyz<T>& s)
        {
            x += s.x;
            y += s.y;
            z += s.z;
            return *this;
        }
        
        sxyz<T>& operator-=(const sxyz<T>& s) {
            x -= s.x;
            y -= s.y;
            z -= s.z;
            return *this;
        }
        
        sxyz<T>& operator=(const T s) {
            x = s;
            y = s;
            z = s;
            return *this;
        }
        
        sxyz<T>& operator/=(const T s) {
            x /= s;
            y /= s;
            z /= s;
            return *this;
        }
        
        template<typename NODE>
        sxyz<T>& operator=(const NODE& n) {
            x = n.getX();
            y = n.getY();
            z = n.getZ();
            return *this;
        }
        
        sxyz<T>& operator *=(const T value) {
            x *= value;
            y *= value;
            z *= value;
            return *this;
        }
        
        bool operator<(const sxyz<T>& s) const {
            return x<s.x && y<s.y && z<s.z;
        }
        
        T getDistance(const sxyz<T>& s) const {
            return sqrt( (x-s.x)*(x-s.x) + (y-s.y)*(y-s.y) + (z-s.z)*(z-s.z) );
        }
        
        void normalize() {
            T n = sqrt( x*x + y*y + z*z );
            x /= n;
            y /= n;
            z /= n;
        }
        
    };
    
    template<typename T>
    std::ostream& operator<< (std::ostream& os, const struct sxyz<T> &s) {
        os << s.x << ' ' << s.y << ' ' << s.z;
        return os;
    }
    
    template<typename T>
    std::istream& operator>> (std::istream& is, struct sxyz<T> &s) {
        is >> s.x >> s.y >> s.z;
        return is;
    }
    
    template <typename T>
    sxyz<T> operator+(const sxyz<T>& lhs, const sxyz<T>& rhs)
    {
        return sxyz<T>(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z);
    }
    
    template <typename T>
    sxyz<T> operator-(const sxyz<T>& lhs, const sxyz<T>& rhs)
    {
        return sxyz<T>(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z);
    }
    
    template <typename T>
    sxyz<T> operator*(const T& scalar, const sxyz<T>& v)
    {
        return sxyz<T>(scalar*v.x, scalar*v.y, scalar*v.z);
    }
    
    template <typename T>
    sxyz<T> operator/(const sxyz<T>& v, const T& scalar)
    {
        return sxyz<T>(v.x/scalar, v.y/scalar, v.z/scalar);
    }
    
    template<typename T>
    sxyz<T> normalize(const sxyz<T>& s) {
        sxyz<T> s2 = s;
        s2.normalize();
        return s2;
    }
    
    
    
    
    
    template<typename T>
    struct sij {
        T i;
        T j;
        
        sij() : i(0), j(0) {}
        sij(const T i_, const T j_) : i(i_), j(j_) {}
    };
    
    template<typename T>
    struct sijk {
        T i;
        T j;
        T k;
        
        sijk() : i(0), j(0), k(0) {}
        sijk(const T i_, const T j_, const T k_) : i(i_), j(j_), k(k_) {}
    };
    
    template<typename T>
    struct siv {
        size_t i;   // index
        T v;        // value
        
        siv<T>& operator+=(const siv<T>& s) {
            v += s.v;
            return *this;
        }
    };
    
    template<typename T>
    struct siv2 {
        size_t i;    // index
        T v;         // first value
        T v2;        // second value
        
        siv2<T>& operator+=(const siv2<T>& s) {
            v += s.v;
            v2 += s.v2;
            return *this;
        }
    };
    
    template<typename T>
    class CompareSiv_i {
    public:
        bool operator()(const siv<T> n1, const siv<T> n2) const {
            return n1.i < n2.i;
        }
    };
    
    template<typename T>
    class CompareSiv_v {
    public:
        bool operator()(const siv<T> n1, const siv<T> n2) const {
            return n1.v < n2.v;
        }
    };
    
    template<typename T>
    class CompareSiv_vr {
    public:
        bool operator()(const siv<T> n1, const siv<T> n2) const {
            return n1.v > n2.v;
        }
    };
    
    template<typename T>
    class CompareSiv2_i {
    public:
        bool operator()(const siv2<T> n1, const siv2<T> n2) const {
            return n1.i < n2.i;
        }
    };
    
    
    template<typename T>
    struct txPar {
        sxz<T> pt;
        T t0;
        T theta;
        T diam;
        bool inWater;
    };
    
    template<typename T>
    struct rxPar {
        std::vector<sxz<T>> pts;
        std::vector<T> theta;
        std::vector<T> diam;
        std::vector<bool> inWater;
    };
    
    template<typename T>
    struct lineElem {
        T i[2];
        T physical_entity;
    };
    
    template<typename T>
    struct triangleElem {
        T i[3];
        T physical_entity;
    };
    
    template<typename T>
    struct tetrahedronElem {
        T i[4];
        T physical_entity;
    };
    
    template<typename T1, typename T2>
    struct triangleElemAngle : triangleElem<T2>{
        //	T2 i[3];    // indices of nodes
        T1 a[3];    // angles at nodes
        T1 l[3];    // length of opposite edges
        
        triangleElemAngle(const triangleElem<T2>& t) {
            this->i[0] = t.i[0];
            this->i[1] = t.i[1];
            this->i[2] = t.i[2];
            this->physical_entity = t.physical_entity;
        }
    };
    
    template<typename T1, typename T2>
    struct tetrahedronElemAngle {
        T2 i[4];
        T1 a[4][3];    // angles at nodes
        T1 l[6];       // length of edges
        
        tetrahedronElemAngle(const tetrahedronElem<T2>& t) {
            i[0] = t.i[0];
            i[1] = t.i[1];
            i[2] = t.i[2];
            i[3] = t.i[3];
        }
    };
    
    template<typename T, typename N>
    struct virtualNode {
        N *node1;
        N *node2;
        T a[3];
        T e[3];
    };
    
    
    template<typename T>
    T det(const sxz<T>& u, const sxz<T>& v) {
        return u.x*v.z - u.z*v.x;
    };
    
    template<typename T>
    T det(const sxyz<T>& v1, const sxyz<T>& v2, const sxyz<T>& v3) {
        return -v1.z*v2.y*v3.x + v1.y*v2.z*v3.x + v1.z*v2.x*v3.y -
        v1.x*v2.z*v3.y - v1.y*v2.x*v3.z + v1.x*v2.y*v3.z;
    };
    
    template<typename T>
    T tripleScalar(const sxyz<T>& v1, const sxyz<T>& v2, const sxyz<T>& v3) {
        return v1.x*(v2.y*v3.z-v3.y*v2.z) - v1.y*(v2.x*v3.z-v3.x*v2.z) + v1.z*(v2.x*v3.y-v3.x*v2.y);
    };
    
    template<typename T>
    T det4(const sxyz<T>& v1, const sxyz<T>& v2,
           const sxyz<T>& v3, const sxyz<T>& v4) {
        return -v1.z*v2.y*v3.x + v1.y*v2.z*v3.x + v1.z*v2.x*v3.y - v1.x*v2.z*v3.y -
        v1.y*v2.x*v3.z + v1.x*v2.y*v3.z + v1.z*v2.y*v4.x - v1.y*v2.z*v4.x -
        v1.z*v3.y*v4.x + v2.z*v3.y*v4.x + v1.y*v3.z*v4.x - v2.y*v3.z*v4.x -
        v1.z*v2.x*v4.y + v1.x*v2.z*v4.y + v1.z*v3.x*v4.y - v2.z*v3.x*v4.y -
        v1.x*v3.z*v4.y + v2.x*v3.z*v4.y + v1.y*v2.x*v4.z - v1.x*v2.y*v4.z -
        v1.y*v3.x*v4.z + v2.y*v3.x*v4.z + v1.x*v3.y*v4.z - v2.x*v3.y*v4.z;
    };
    
    template<typename T>
    T dot(const sxz<T>& v1, const sxz<T>& v2) {
        return v1.x*v2.x + v1.z*v2.z;
    }
    
    template<typename T>
    T cross(const sxz<T>& v1, const sxz<T>& v2) {
        return v1.z*v2.x - v1.x*v2.z;
    }
    
    template<typename T>
    T dot(const sxyz<T>& v1, const sxyz<T>& v2) {
        return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
    }
    
    template<typename T>
    sxyz<T> cross(const sxyz<T>& v1, const sxyz<T>& v2) {
        sxyz<T> v3;
        v3.x = v1.y*v2.z - v1.z*v2.y;
        v3.y = v1.z*v2.x - v1.x*v2.z;
        v3.z = v1.x*v2.y - v1.y*v2.x;
        return v3;
    }
    
    
    template<typename T>
    T norm(const sxz<T>& v) {
        return sqrt( v.x*v.x + v.z*v.z );
    }
    
    template<typename T>
    T norm(const sxyz<T>& v) {
        return sqrt( v.x*v.x + v.y*v.y + v.z*v.z );
    }
    
    template<typename T>
    T norm2(const sxyz<T>& v) {
        return v.x*v.x + v.y*v.y + v.z*v.z;
    }
    
    
    // following 3 fct from
    // http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
    template <typename T> inline constexpr
    int signum(T x, std::false_type is_signed) {
        return T(0) < x;
    }
    
    template <typename T> inline constexpr
    int signum(T x, std::true_type is_signed) {
        return (T(0) < x) - (x < T(0));
    }
    
    template <typename T> inline constexpr
    int signum(T x) {
        return signum(x, std::is_signed<T>());
    }
    
}

#endif

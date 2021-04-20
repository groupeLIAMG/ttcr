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

#ifndef ttcr_ttrc_t_h
#define ttcr_ttrc_t_h

#include <cmath>
#include <type_traits>
#include <vector>

namespace ttcr {

    extern int verbose;

    const double small = 1.e-4;
    const double small2 = small*small;
    const double small3 = small*small*small;
    const double pi = 4.0*atan(1.0);
    const double theta_cut = 65. * pi / 180.;  // for raytracing -> unstructured meshes

    const size_t iLength[4][3]={{0,1,2},{1,3,4},{2,3,5},{0,4,5}};
    const size_t iNodes[4][3] = {
        {0,1,2},  // (relative) indices of nodes of 1st triangle
        {1,2,3},  // (relative) indices of nodes of 2nd triangle
        {0,2,3},  // (relative) indices of nodes of 3rd triangle
        {0,1,3}   // (relative) indices of nodes of 4th triangle
    };

    /*

     1
     ,/|`\
     ,/  |  `\
     ,0    '.   `4
     ,/       1     `\
     ,/         |       `\
     0-----5-----'.--------3
     `\.         |      ,/
     `\.      |     3
     `2.   '. ,/
     `\. |/
     `2


     triangle 0:  0-1  1-2  2-0     (first occurence of edge underlined)
     ---  ---  ---
     triangle 1:  1-2  2-3  3-1
     ---  ---
     triangle 2:  0-2  2-3  3-0
     ---
     triangle 3:  0-1  1-3  3-0

     for triangle "itri", indices of length for edge "iedge" are
     iLength[itri][iedge] = {{0,1,2},{1,3,4},{2,3,5},{0,4,5}}

     */

    template<typename T>
    struct sxz {
        T x;
        T z;

        sxz() : x(0), z(0) {}
        sxz(const T x_, const T z_) : x(x_), z(z_) {}
        sxz(const T v) : x(v), z(v) {}
        template<typename NODE>
        sxz(const NODE& n) : x(n.getX()), z(n.getZ()) {}

        bool operator==(const sxz<T>& s) const {
            //        return (std::abs(x-a.x)<small && std::abs(z-a.z)<small);
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

    template <typename T, typename NODE>
    sxz<T> operator+(const NODE& lhs, const sxz<T>& rhs)
    {
        return sxz<T>(lhs.getX()+rhs.x, lhs.getZ()+rhs.z);
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
        sxyz(const T v) : x(v), y(v), z(v) {}
        template<typename NODE>
        sxyz(const NODE& n) : x(n.getX()), y(n.getY()), z(n.getZ()) {}

        bool operator==(const sxyz<T>& s) const {
            //        return (std::abs(x-a.x)<small && std::abs(z-a.z)<small);
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

    template <typename T, typename NODE>
    sxyz<T> operator+(const NODE& lhs, const sxyz<T>& rhs)
    {
        return sxyz<T>(lhs.getX()+rhs.x, lhs.getY()+rhs.y, lhs.getZ()+rhs.z);
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

        siv() : i(0), v(0) {}
        siv(const size_t i_, const T v_) : i(i_), v(v_) {}

        siv<T>& operator+=(const siv<T>& s) {
            v += s.v;
            return *this;
        }
    };

    template<typename T>
    struct sijv {
        size_t i;  // index in 1st dim
        size_t j;  // index in 2nd dim
        T v;       // value

        sijv() : i(0), j(0), v(0) {}
        sijv(const size_t i_, const size_t j_, const T v_) : i(i_), j(j_), v(v_) {}

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

        triangleElem() : i{ 0,0,0 }, physical_entity(0) {}
        triangleElem(const T i0, const T i1, const T i2, const T pe=0) : i{ i0,i1,i2 }, physical_entity(pe) {}
    };

    template<typename T>
    struct tetrahedronElem {
        T i[4];
        T physical_entity;

        tetrahedronElem() : i{ 0,0,0,0 }, physical_entity(0) {}
        tetrahedronElem(const T i0, const T i1, const T i2, const T i3, const T pe=0) : i{ i0,i1,i2,i3 }, physical_entity(pe) {}

        tetrahedronElem(const tetrahedronElem& tet) : i{ tet.i[0],tet.i[1],tet.i[2],tet.i[3] }, physical_entity(tet.physical_entity) {}
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
    }

    template<typename T>
    T det(const sxyz<T>& v1, const sxyz<T>& v2, const sxyz<T>& v3) {
        return -v1.z*v2.y*v3.x + v1.y*v2.z*v3.x + v1.z*v2.x*v3.y -
        v1.x*v2.z*v3.y - v1.y*v2.x*v3.z + v1.x*v2.y*v3.z;
    }

    template<typename T>
    T tripleScalar(const sxyz<T>& v1, const sxyz<T>& v2, const sxyz<T>& v3) {
        return v1.x*(v2.y*v3.z-v3.y*v2.z) - v1.y*(v2.x*v3.z-v3.x*v2.z) + v1.z*(v2.x*v3.y-v3.x*v2.y);
    }

    template<typename T>
    T det4(const sxyz<T>& v1, const sxyz<T>& v2,
           const sxyz<T>& v3, const sxyz<T>& v4) {
        return -v1.z*v2.y*v3.x + v1.y*v2.z*v3.x + v1.z*v2.x*v3.y - v1.x*v2.z*v3.y -
        v1.y*v2.x*v3.z + v1.x*v2.y*v3.z + v1.z*v2.y*v4.x - v1.y*v2.z*v4.x -
        v1.z*v3.y*v4.x + v2.z*v3.y*v4.x + v1.y*v3.z*v4.x - v2.y*v3.z*v4.x -
        v1.z*v2.x*v4.y + v1.x*v2.z*v4.y + v1.z*v3.x*v4.y - v2.z*v3.x*v4.y -
        v1.x*v3.z*v4.y + v2.x*v3.z*v4.y + v1.y*v2.x*v4.z - v1.x*v2.y*v4.z -
        v1.y*v3.x*v4.z + v2.y*v3.x*v4.z + v1.x*v3.y*v4.z - v2.x*v3.y*v4.z;
    }

    template<typename T>
    T dot(const sxz<T>& v1, const sxz<T>& v2) {
        return v1.x*v2.x + v1.z*v2.z;
    }

    //    template<typename T>
    //    T dot(const T v1, const sxz<T>& v2) {
    //        // v1 considered to be a vector along y
    //        return 0.0;
    //    }

    template<typename T>
    T dot(const T v1, const T v2) {
        // both v1 & v2 considered to be a vector along y
        return v1*v2;
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
    T norm2(const sxz<T>& v) {
        return v.x*v.x + v.z*v.z;
    }

    template<typename T>
    T norm2(const sxyz<T>& v) {
        return v.x*v.x + v.y*v.y + v.z*v.z;
    }

    template<typename T>
    void projNorm(const sxyz<T>& b, const sxyz<T>& c, const sxyz<T>& p, T& xi0, T& zeta0) {
        /*
         B
         o
         /   \
         /          \
         /       o P       \
         A o----------------------o  C

         b is _unit_ vector AC
         c is _unit_ vector AB
         p is vector AP
         xi0 & zeta0 are normalized coordinates of P in plane ABC (Lelièvre et al, 2011, GJI 184, 885-896)

         solved using xi*c + zeta*b = p

         | c.x b.x |                        | p.x |
         A =  | c.y b.y |   x = |  xi0  |    b = | p.y |
         | c.z b.z |       | zeta0 |        | p.z |

         solve AT A x = AT b
         */

        T ata11 = norm2(c);
        T ata12 = b.x*c.x + b.y*c.y + b.z*c.z;
        T ata21 = ata12;
        T ata22 = norm2(b);
        T atb1 = c.x*p.x + c.y*p.y + c.z*p.z;
        T atb2 = b.x*p.x + b.y*p.y + b.z*p.z;

        T det = ata11*ata22 - ata12*ata21;
        if ( det == 0.0 ) {
            xi0 = -1.0;
            zeta0 = -1.0;
        } else {
            xi0   = (atb1*ata22 - ata12*atb2) / det;
            zeta0 = (ata11*atb2 - atb1*ata21) / det;
        }
    }

#ifndef _MSC_VER
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
#else
    template <typename T>
    int signum(T x, std::false_type is_signed) {
        return T(0) < x;
    }

    template <typename T>
    int signum(T x, std::true_type is_signed) {
        return (T(0) < x) - (x < T(0));
    }

    template <typename T>
    int signum(T x) {
        return signum(x, std::is_signed<T>());
    }
#endif

    template<typename T>
    bool sameSide(const sxyz<T>& v1, const sxyz<T>& v2, const sxyz<T>& v3,
                  const sxyz<T>& v4,  const sxyz<T>& p) {
        sxyz<T> normal = cross(v2 - v1, v3 - v1);
        T dotV4 = dot(normal, v4 - v1);
        T dotP = dot(normal, p - v1);
        return signum(dotV4) == signum(dotP);
    }

}

#endif

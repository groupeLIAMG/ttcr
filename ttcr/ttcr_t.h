//
//  ttcr_t.h
//  ttcr
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

/**
 * @file ttcr_t.h
 * @brief Core types and geometric helpers for the ttcr traveltime library.
 *
 * This header defines the small, dependency-free building blocks shared across
 * the ttcr grids and meshes: 2D/3D coordinate structs (@ref ttcr::sxz,
 * @ref ttcr::sxyz), index/value pairs used to assemble sparse sensitivity
 * matrices (@ref ttcr::siv, @ref ttcr::sijv, ...), mesh element descriptors
 * (@ref ttcr::triangleElem, @ref ttcr::tetrahedronElem, ...), and a collection
 * of free functions implementing common vector algebra (dot, cross,
 * determinants, norms, projections).
 *
 * Most types are templated on the numeric type @c T (typically @c double or
 * @c float); index types are templated separately where the two differ.
 */

#ifndef ttcr_ttrc_t_h
#define ttcr_ttrc_t_h

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

/// Namespace grouping all ttcr traveltime-computation types and helpers.
namespace ttcr {

    /// Global verbosity level; higher values emit more diagnostic output.
    extern int verbose;

    /**
     * @brief Non-zero enables the OpenCL GPU profiling breakdown.
     *
     * Set from the "profile" parameter-file keyword and read directly by the
     * OpenCL grids, mirroring how @ref verbose is consumed.
     */
    extern int gpu_profile;

    /**
     * @brief Accumulate the sweep-to-sweep change and measure the traveltime scale.
     *
     * Walks the nodes once, summing @f$|\Delta t|@f$ against the previous sweep
     * and refreshing @p times, and at the same time measures the **range** of
     * the traveltime field.  That range is what @ref fsmTolerance turns the
     * dimensionless tolerance into an absolute one, so the stop test means the
     * same thing whatever units the model is expressed in.
     *
     * The range is used rather than the largest traveltime because the field
     * carries the source origin time: @c t0 may be an absolute timestamp, and
     * @f$t_{max}@f$ alone would then be dominated by it and the tolerance
     * loosened by orders of magnitude.  Nodes the sweeps have not reached yet
     * still hold @c std::numeric_limits<T1>::max() and are left out of both
     * ends; before the first sweep every node is unreached, so the range comes
     * back as zero and the tolerance with it, which keeps the loop running.
     *
     * @param[in]     nodes    grid nodes, queried through @c getTT.
     * @param[in,out] times    previous sweep's traveltimes; overwritten with
     *                         the current ones.
     * @param[in]     threadNo thread whose traveltimes to read.
     * @param[out]    tref     traveltime range of the reached nodes.
     * @return The summed absolute change over all nodes (an L1 sum, so it
     *         grows with the node count; @ref fsmTolerance scales to match).
     */
    template<typename T1, typename NODES>
    T1 fsmChange(const NODES& nodes, std::vector<T1>& times,
                 const size_t threadNo, T1& tref) {
        const T1 unreached = std::numeric_limits<T1>::max();
        T1 change = 0.0;
        T1 tmin = unreached;
        T1 tmax = 0.0;
        for ( size_t n=0; n<nodes.size(); ++n ) {
            const T1 t = nodes[n].getTT(threadNo);
            change += std::abs( times[n] - t );
            times[n] = t;
            if ( t < unreached ) {
                if ( t > tmax ) tmax = t;
                if ( t < tmin ) tmin = t;
            }
        }
        tref = tmax > tmin ? tmax - tmin : static_cast<T1>(0);
        return change;
    }

    /**
     * @brief @ref fsmChange over a plain traveltime vector.
     *
     * The OpenCL grids hold the field in a buffer downloaded from the device
     * rather than in the nodes, so they need the same reduction over a
     * @c std::vector.
     *
     * @param[in]     tt    current traveltimes.
     * @param[in,out] times previous sweep's traveltimes; overwritten with @p tt.
     * @param[out]    tref  traveltime range of the reached entries.
     * @return The summed absolute change.
     */
    template<typename T1>
    T1 fsmChange(const std::vector<T1>& tt, std::vector<T1>& times, T1& tref) {
        const T1 unreached = std::numeric_limits<T1>::max();
        T1 change = 0.0;
        T1 tmin = unreached;
        T1 tmax = 0.0;
        for ( size_t n=0; n<tt.size(); ++n ) {
            change += std::abs( times[n] - tt[n] );
            times[n] = tt[n];
            if ( tt[n] < unreached ) {
                if ( tt[n] > tmax ) tmax = tt[n];
                if ( tt[n] < tmin ) tmin = tt[n];
            }
        }
        tref = tmax > tmin ? tmax - tmin : static_cast<T1>(0);
        return change;
    }

    /**
     * @brief Absolute stop threshold for a sweep loop.
     *
     * @c epsilon (ttcr::input_parameters::epsilon) is a **dimensionless**
     * tolerance: the sweeps stop once the mean change per node has fallen below
     * that fraction of the traveltime range.  The loops compare against an L1
     * sum, so the node count is folded in here rather than into @c epsilon at
     * construction.
     *
     * @param epsilon relative tolerance, as supplied by the caller.
     * @param tref    traveltime range, from @ref fsmChange.
     * @param nnodes  number of grid nodes.
     * @return The value the summed change must fall below.
     */
    template<typename T1>
    inline T1 fsmTolerance(const T1 epsilon, const T1 tref, const size_t nnodes) {
        return epsilon * tref * static_cast<T1>(nnodes);
    }

    /**
     * @brief Report a fast-sweeping pass that stopped on its iteration cap.
     *
     * A sweep loop ends either because the change in nodal traveltime dropped
     * below the tolerance, or because it ran out of iterations.  The second
     * case is not a converged solve, and it used to be silent: the traveltimes
     * were returned as if they were final.  Callers of the fast-sweeping grids
     * must call this when the loop exits with @c niter equal to @c nitermax.
     *
     * @param pass    name of the pass, e.g. @c "WENO3", used in the message.
     * @param niter   number of iterations run, i.e. the cap.
     * @param change  summed absolute change over all nodes in the last sweep.
     * @param epsilon relative tolerance, as supplied by the caller.
     * @param tref    traveltime range, from @ref fsmChange.
     * @param nnodes  number of grid nodes, used to turn the L1 sum back into
     *                the per-node figures the caller reasons about.
     *
     * @note Written to @c std::cerr regardless of @ref verbose: a capped solve
     *       is a wrong answer, not a diagnostic detail.
     * @sa ttcr::input_parameters::epsilon, ttcr::input_parameters::nitermax
     */
    inline void warnFSMnotConverged(const std::string& pass, const int niter,
                                    const double change, const double epsilon,
                                    const double tref, const size_t nnodes) {
        const double scale = nnodes > 0 ? static_cast<double>(nnodes) : 1.0;
        const double mean = change/scale;
        std::cerr << "Warning: " << pass << " fast sweeping stopped after "
                  << niter << " iterations (nitermax) without converging: "
                  << "mean change per node is " << mean;
        if ( tref > 0.0 ) {
            std::cerr << ", i.e. " << mean/tref << " of the traveltime range ("
                      << tref << "), above epsilon = " << epsilon;
        } else {
            std::cerr << ", and the traveltime range is not yet established";
        }
        std::cerr << ".  The traveltimes are not converged; raise nitermax "
                  << "and/or lower epsilon.\n";
    }

    const double small = 1.e-4;              ///< Small tolerance used for geometric comparisons.
    const double small2 = small*small;       ///< @ref small squared.
    const double small3 = small*small*small; ///< @ref small cubed.
    const double pi = 4.0*atan(1.0);         ///< The mathematical constant &pi;.
    const double theta_cut = 65. * pi / 180.;  ///< Angle threshold (radians) for raytracing on unstructured meshes.

    /// Per-triangle edge-length indices; see the diagram below.
    const size_t iLength[4][3]={{0,1,2},{1,3,4},{2,3,5},{0,4,5}};
    /// Relative node indices of the four triangular faces of a tetrahedron.
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

    /**
     * @brief A point or vector in the 2D (x, z) plane.
     * @tparam T Numeric type of the coordinates (e.g. @c double or @c float).
     */
    template<typename T>
    struct sxz {
        T x;  ///< Horizontal coordinate.
        T z;  ///< Vertical (depth) coordinate.

        /// Constructs the origin (0, 0).
        sxz() : x(0), z(0) {}
        /// Constructs from explicit @p x_ and @p z_ coordinates.
        sxz(const T x_, const T z_) : x(x_), z(z_) {}
        /// Constructs with both coordinates set to @p v.
        sxz(const T v) : x(v), z(v) {}
        /// Constructs from any node type exposing @c getX() and @c getZ().
        template<typename NODE>
        sxz(const NODE& n) : x(n.getX()), z(n.getZ()) {}

        /// Equality: true when both coordinates match exactly.
        bool operator==(const sxz<T>& s) const {
            return x==s.x && z==s.z;
        }

        /// Component-wise addition.
        sxz<T>& operator+=(const sxz<T>& s)
        {
            x += s.x;
            z += s.z;
            return *this;
        }

        /// Component-wise subtraction.
        sxz<T>& operator-=(const sxz<T>& s) {
            x -= s.x;
            z -= s.z;
            return *this;
        }

        /// Divides both coordinates by the scalar @p s.
        sxz<T>& operator/=(const T s) {
            x /= s;
            z /= s;
            return *this;
        }

        /// Assigns the scalar @p s to both coordinates.
        sxz<T>& operator=(const T s) {
            x = s;
            z = s;
            return *this;
        }

        /// Copy assignment (self-assignment safe).
        sxz<T> & operator=(const sxz<T>& s) {
            if (this != &s) { // protect against invalid self-assignment
                x = s.x;
                z = s.z;
            }
            return *this;
        }

        /// Assigns from any node type exposing @c getX() and @c getZ().
        template<typename NODE>
        sxz<T>& operator=(const NODE& n) {
            x = n.getX();
            z = n.getZ();
            return *this;
        }

        /// Multiplies both coordinates by the scalar @p value.
        sxz<T>& operator *=(const T value) {
            x *= value;
            z *= value;
            return *this;
        }

        /// Strict ordering: true only when @e both coordinates are less than @p s's.
        bool operator<(const sxz<T>& s) const {
            return x<s.x && z<s.z;
        }

        /// Returns the Euclidean distance to point @p s.
        T getDistance(const sxz<T>& s) const {
            return sqrt( (x-s.x)*(x-s.x) + (z-s.z)*(z-s.z) );
        }

        /// Scales this vector in place to unit length.
        void normalize() {
            T n = sqrt( x*x + z*z );
            x /= n;
            z /= n;
        }
    };

    /// Stream-inserts an @ref sxz as space-separated "x z".
    template<typename T>
    std::ostream& operator<< (std::ostream& os, const struct sxz<T> &s) {
        os << s.x << ' ' << s.z;
        return os;
    }

    /// Stream-extracts an @ref sxz from two whitespace-separated values.
    template<typename T>
    std::istream& operator>> (std::istream& is, struct sxz<T> &s) {
        is >> s.x >> s.z;
        return is;
    }

    /// Adds a node (via @c getX()/@c getZ()) and an @ref sxz.
    template <typename T, typename NODE>
    sxz<T> operator+(const NODE& lhs, const sxz<T>& rhs)
    {
        return sxz<T>(lhs.getX()+rhs.x, lhs.getZ()+rhs.z);
    }

    /// Component-wise sum of two @ref sxz vectors.
    template <typename T>
    sxz<T> operator+(const sxz<T>& lhs, const sxz<T>& rhs)
    {
        return sxz<T>(lhs.x+rhs.x, lhs.z+rhs.z);
    }

    /// Component-wise difference of two @ref sxz vectors.
    template <typename T>
    sxz<T> operator-(const sxz<T>& lhs, const sxz<T>& rhs)
    {
        return sxz<T>(lhs.x-rhs.x, lhs.z-rhs.z);
    }

    /// Scalar-vector product.
    template <typename T>
    sxz<T> operator*(const T& scalar, const sxz<T>& v)
    {
        return sxz<T>(scalar*v.x, scalar*v.z);
    }

    /// Divides an @ref sxz vector by a scalar.
    template <typename T>
    sxz<T> operator/(const sxz<T>& v, const T& scalar)
    {
        return sxz<T>(v.x/scalar, v.z/scalar);
    }

    /// Returns a unit-length copy of @p s.
    template<typename T>
    sxz<T> normalize(const sxz<T>& s) {
        sxz<T> s2 = s;
        s2.normalize();
        return s2;
    }



    /**
     * @brief A point or vector in 3D (x, y, z) space.
     * @tparam T Numeric type of the coordinates (e.g. @c double or @c float).
     */
    template<typename T>
    struct sxyz {
        T x;  ///< x coordinate.
        T y;  ///< y coordinate.
        T z;  ///< z (depth) coordinate.

        /// Constructs the origin (0, 0, 0).
        sxyz() : x(0), y(0), z(0) {}
        /// Constructs from explicit @p x_, @p y_ and @p z_ coordinates.
        sxyz(const T x_, const T y_, const T z_) : x(x_), y(y_), z(z_) {}
        /// Constructs with all three coordinates set to @p v.
        sxyz(const T v) : x(v), y(v), z(v) {}
        /// Constructs from any node type exposing @c getX(), @c getY() and @c getZ().
        template<typename NODE>
        sxyz(const NODE& n) : x(n.getX()), y(n.getY()), z(n.getZ()) {}

        T getX() const { return x; }  ///< Returns the x coordinate.
        T getY() const { return y; }  ///< Returns the y coordinate.
        T getZ() const { return z; }  ///< Returns the z coordinate.

        /// Equality: true when all three coordinates match exactly.
        bool operator==(const sxyz<T>& s) const {
            return x==s.x && y==s.y && z==s.z;
        }

        /// Copy assignment (self-assignment safe).
        sxyz<T>& operator=(const sxyz<T>& s) {
            if (this != &s) { // protect against invalid self-assignment
                x = s.x;
                y = s.y;
                z = s.z;
            }
            return *this;
        }

        /// Component-wise addition.
        sxyz<T>& operator+=(const sxyz<T>& s)
        {
            x += s.x;
            y += s.y;
            z += s.z;
            return *this;
        }

        /// Component-wise subtraction.
        sxyz<T>& operator-=(const sxyz<T>& s) {
            x -= s.x;
            y -= s.y;
            z -= s.z;
            return *this;
        }

        /// Assigns the scalar @p s to all three coordinates.
        sxyz<T>& operator=(const T s) {
            x = s;
            y = s;
            z = s;
            return *this;
        }

        /// Divides all three coordinates by the scalar @p s.
        sxyz<T>& operator/=(const T s) {
            x /= s;
            y /= s;
            z /= s;
            return *this;
        }

        /// Assigns from any node type exposing @c getX(), @c getY() and @c getZ().
        template<typename NODE>
        sxyz<T>& operator=(const NODE& n) {
            x = n.getX();
            y = n.getY();
            z = n.getZ();
            return *this;
        }

        /// Multiplies all three coordinates by the scalar @p value.
        sxyz<T>& operator *=(const T value) {
            x *= value;
            y *= value;
            z *= value;
            return *this;
        }

        /// Strict ordering: true only when @e all coordinates are less than @p s's.
        bool operator<(const sxyz<T>& s) const {
            return x<s.x && y<s.y && z<s.z;
        }

        /// Returns the Euclidean distance to point @p s.
        T getDistance(const sxyz<T>& s) const {
            return sqrt( (x-s.x)*(x-s.x) + (y-s.y)*(y-s.y) + (z-s.z)*(z-s.z) );
        }

        /// Scales this vector in place to unit length.
        void normalize() {
            T n = sqrt( x*x + y*y + z*z );
            x /= n;
            y /= n;
            z /= n;
        }

    };

    /// Stream-inserts an @ref sxyz as space-separated "x y z".
    template<typename T>
    std::ostream& operator<< (std::ostream& os, const struct sxyz<T> &s) {
        os << s.x << ' ' << s.y << ' ' << s.z;
        return os;
    }

    /// Stream-extracts an @ref sxyz from three whitespace-separated values.
    template<typename T>
    std::istream& operator>> (std::istream& is, struct sxyz<T> &s) {
        is >> s.x >> s.y >> s.z;
        return is;
    }

    /// Adds a node (via @c getX()/@c getY()/@c getZ()) and an @ref sxyz.
    template <typename T, typename NODE>
    sxyz<T> operator+(const NODE& lhs, const sxyz<T>& rhs)
    {
        return sxyz<T>(lhs.getX()+rhs.x, lhs.getY()+rhs.y, lhs.getZ()+rhs.z);
    }

    /// Component-wise sum of two @ref sxyz vectors.
    template <typename T>
    sxyz<T> operator+(const sxyz<T>& lhs, const sxyz<T>& rhs)
    {
        return sxyz<T>(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z);
    }

    /// Component-wise difference of two @ref sxyz vectors.
    template <typename T>
    sxyz<T> operator-(const sxyz<T>& lhs, const sxyz<T>& rhs)
    {
        return sxyz<T>(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z);
    }

    /// Scalar-vector product.
    template <typename T>
    sxyz<T> operator*(const T& scalar, const sxyz<T>& v)
    {
        return sxyz<T>(scalar*v.x, scalar*v.y, scalar*v.z);
    }

    /// Divides an @ref sxyz vector by a scalar.
    template <typename T>
    sxyz<T> operator/(const sxyz<T>& v, const T& scalar)
    {
        return sxyz<T>(v.x/scalar, v.y/scalar, v.z/scalar);
    }

    /// Returns a unit-length copy of @p s.
    template<typename T>
    sxyz<T> normalize(const sxyz<T>& s) {
        sxyz<T> s2 = s;
        s2.normalize();
        return s2;
    }



    /**
     * @brief A pair of integer indices (i, j).
     * @tparam T Integer index type.
     */
    template<typename T>
    struct sij {
        T i;  ///< First index.
        T j;  ///< Second index.

        /// Constructs (0, 0).
        sij() : i(0), j(0) {}
        /// Constructs from explicit indices @p i_ and @p j_.
        sij(const T i_, const T j_) : i(i_), j(j_) {}
    };

    /**
     * @brief A triple of integer indices (i, j, k).
     * @tparam T Integer index type.
     */
    template<typename T>
    struct sijk {
        T i;  ///< First index.
        T j;  ///< Second index.
        T k;  ///< Third index.

        /// Constructs (0, 0, 0).
        sijk() : i(0), j(0), k(0) {}
        /// Constructs from explicit indices @p i_, @p j_ and @p k_.
        sijk(const T i_, const T j_, const T k_) : i(i_), j(j_), k(k_) {}
    };

    /**
     * @brief An (index, value) pair, e.g. one nonzero entry of a sparse vector.
     * @tparam T Numeric type of the value.
     */
    template<typename T>
    struct siv {
        size_t i;   ///< Index.
        T v;        ///< Value.

        /// Constructs (index 0, value 0).
        siv() : i(0), v(0) {}
        /// Constructs from an explicit index @p i_ and value @p v_.
        siv(const size_t i_, const T v_) : i(i_), v(v_) {}

        /// Accumulates the value of @p s (indices are left unchanged).
        siv<T>& operator+=(const siv<T>& s) {
            v += s.v;
            return *this;
        }
    };

    /// Stream-inserts an @ref siv as space-separated "i v".
    template<typename T>
    std::ostream& operator<< (std::ostream& os, const struct siv<T> &s) {
        os << s.i << ' ' << s.v;
        return os;
    }

    /**
     * @brief An (i, j, value) triple, e.g. one nonzero entry of a sparse matrix.
     * @tparam T Numeric type of the value.
     */
    template<typename T>
    struct sijv {
        size_t i;  ///< Index in the 1st dimension (row).
        size_t j;  ///< Index in the 2nd dimension (column).
        T v;       ///< Value.

        /// Constructs (0, 0, value 0).
        sijv() : i(0), j(0), v(0) {}
        /// Constructs from explicit indices @p i_, @p j_ and value @p v_.
        sijv(const size_t i_, const size_t j_, const T v_) : i(i_), j(j_), v(v_) {}

    };

    /**
     * @brief An index paired with two values.
     * @tparam T Numeric type of the values.
     */
    template<typename T>
    struct siv2 {
        size_t i;    ///< Index.
        T v;         ///< First value.
        T v2;        ///< Second value.

        /// Accumulates both values of @p s (index is left unchanged).
        siv2<T>& operator+=(const siv2<T>& s) {
            v += s.v;
            v2 += s.v2;
            return *this;
        }
    };

    /**
     * @brief An index paired with four values.
     *
     * Used to hold the sensitivity of the traveltime to the medium parameters
     * of one cell, one value per parameter.  The meaning and the order of the
     * values are those documented by the cell class filling the struct; a class
     * describing fewer than four parameters leaves the trailing values at zero,
     * the number of meaningful values being given by its @c nParams member.
     *
     * @tparam T Numeric type of the values.
     */
    template<typename T>
    struct siv4 {
        size_t i;    ///< Index.
        T v;         ///< First value.
        T v2;        ///< Second value.
        T v3;        ///< Third value.
        T v4;        ///< Fourth value.

        /// Constructs (index 0, values 0).
        siv4() : i(0), v(0), v2(0), v3(0), v4(0) {}

        /// Accumulates the four values of @p s (index is left unchanged).
        siv4<T>& operator+=(const siv4<T>& s) {
            v += s.v;
            v2 += s.v2;
            v3 += s.v3;
            v4 += s.v4;
            return *this;
        }
    };

    /**
     * @brief An index paired with five values.
     *
     * Five-parameter counterpart of @ref siv4, needed by the tilted
     * transversely isotropic cells, whose medium parameters are @f$V_{P0}@f$,
     * @f$V_{S0}@f$, @f$\epsilon@f$, @f$\delta@f$ and the tilt angle.
     *
     * @tparam T Numeric type of the values.
     */
    template<typename T>
    struct siv5 {
        size_t i;    ///< Index.
        T v;         ///< First value.
        T v2;        ///< Second value.
        T v3;        ///< Third value.
        T v4;        ///< Fourth value.
        T v5;        ///< Fifth value.

        /// Constructs (index 0, values 0).
        siv5() : i(0), v(0), v2(0), v3(0), v4(0), v5(0) {}

        /// Accumulates the five values of @p s (index is left unchanged).
        siv5<T>& operator+=(const siv5<T>& s) {
            v += s.v;
            v2 += s.v2;
            v3 += s.v3;
            v4 += s.v4;
            v5 += s.v5;
            return *this;
        }
    };

    /// Comparator ordering @ref siv by ascending index.
    template<typename T>
    class CompareSiv_i {
    public:
        bool operator()(const siv<T> n1, const siv<T> n2) const {
            return n1.i < n2.i;
        }
    };

    /// Comparator ordering @ref siv by ascending value.
    template<typename T>
    class CompareSiv_v {
    public:
        bool operator()(const siv<T> n1, const siv<T> n2) const {
            return n1.v < n2.v;
        }
    };

    /// Comparator ordering @ref siv by descending value (reverse).
    template<typename T>
    class CompareSiv_vr {
    public:
        bool operator()(const siv<T> n1, const siv<T> n2) const {
            return n1.v > n2.v;
        }
    };

    /// Comparator ordering @ref siv2 by ascending index.
    template<typename T>
    class CompareSiv2_i {
    public:
        bool operator()(const siv2<T> n1, const siv2<T> n2) const {
            return n1.i < n2.i;
        }
    };

    /// Comparator ordering @ref siv4 by ascending index.
    template<typename T>
    class CompareSiv4_i {
    public:
        bool operator()(const siv4<T> n1, const siv4<T> n2) const {
            return n1.i < n2.i;
        }
    };

    /// Comparator ordering @ref siv5 by ascending index.
    template<typename T>
    class CompareSiv5_i {
    public:
        bool operator()(const siv5<T> n1, const siv5<T> n2) const {
            return n1.i < n2.i;
        }
    };

    /**
     * @brief Whether a cell class reports its sensitivity into a given container
     *
     * The number of values a cell reports is the number of medium parameters it
     * describes, so each cell class provides computeDistance() for one container
     * only.  A grid, however, has to override every raytrace() of Grid2D,
     * whichever cells it holds, and the overrides are emitted with the vtable
     * whether or not they are ever called.  This tells them apart, so that the
     * ones the cells cannot serve fail when called rather than when compiled.
     */
    template<typename CELL, typename SIV, typename NODE, typename S,
             typename = void>
    struct cell_reports_into : std::false_type {};

    template<typename CELL, typename SIV, typename NODE, typename S>
    struct cell_reports_into<CELL, SIV, NODE, S,
        std::void_t<decltype(std::declval<const CELL&>().computeDistance(
            std::declval<const NODE&>(), std::declval<const S&>(),
            std::declval<SIV&>()))>> : std::true_type {};



    /**
     * @brief Parameters of a single transmitter (source).
     * @tparam T Numeric type.
     */
    template<typename T>
    struct txPar {
        sxz<T> pt;      ///< Source location.
        T t0;           ///< Source excitation time.
        T theta;        ///< Emission half-angle (aperture).
        T diam;         ///< Source diameter.
        bool inWater;   ///< True if the source lies in water.
    };

    /**
     * @brief Parameters of a set of receivers.
     * @tparam T Numeric type.
     */
    template<typename T>
    struct rxPar {
        std::vector<sxz<T>> pts;      ///< Receiver locations.
        std::vector<T> theta;         ///< Per-receiver acceptance half-angle.
        std::vector<T> diam;          ///< Per-receiver diameter.
        std::vector<bool> inWater;    ///< Per-receiver in-water flag.
    };

    /**
     * @brief A line (2-node) mesh element.
     * @tparam T Integer node-index type.
     */
    template<typename T>
    struct lineElem {
        T i[2];             ///< Indices of the two end nodes.
        T physical_entity;  ///< Physical-entity/material tag.
    };

    /**
     * @brief A triangular (3-node) mesh element.
     * @tparam T Integer node-index type.
     */
    template<typename T>
    struct triangleElem {
        T i[3];             ///< Indices of the three corner nodes.
        T physical_entity;  ///< Physical-entity/material tag.

        /// Constructs a zero-initialized element.
        triangleElem() : i{ 0,0,0 }, physical_entity(0) {}
        /// Constructs from three node indices and an optional physical-entity tag @p pe.
        triangleElem(const T i0, const T i1, const T i2, const T pe=0) : i{ i0,i1,i2 }, physical_entity(pe) {}
    };

    /**
     * @brief A tetrahedral (4-node) mesh element.
     * @tparam T Integer node-index type.
     */
    template<typename T>
    struct tetrahedronElem {
        T i[4];             ///< Indices of the four corner nodes.
        T physical_entity;  ///< Physical-entity/material tag.

        /// Constructs a zero-initialized element.
        tetrahedronElem() : i{ 0,0,0,0 }, physical_entity(0) {}
        /// Constructs from four node indices and an optional physical-entity tag @p pe.
        tetrahedronElem(const T i0, const T i1, const T i2, const T i3, const T pe=0) : i{ i0,i1,i2,i3 }, physical_entity(pe) {}

        /// Copy constructor.
        tetrahedronElem(const tetrahedronElem& tet) : i{ tet.i[0],tet.i[1],tet.i[2],tet.i[3] }, physical_entity(tet.physical_entity) {}
    };

    /**
     * @brief A triangular element augmented with cached geometry.
     * @tparam T1 Numeric type for angles and lengths.
     * @tparam T2 Integer node-index type (inherited @ref triangleElem parameter).
     */
    template<typename T1, typename T2>
    struct triangleElemAngle : triangleElem<T2>{
        //	T2 i[3];    // indices of nodes
        T1 a[3];    ///< Interior angles at the three nodes.
        T1 l[3];    ///< Lengths of the edges opposite each node.

        /// Constructs from a plain @ref triangleElem, copying node indices and tag.
        triangleElemAngle(const triangleElem<T2>& t) {
            this->i[0] = t.i[0];
            this->i[1] = t.i[1];
            this->i[2] = t.i[2];
            this->physical_entity = t.physical_entity;
        }
    };

    /**
     * @brief A tetrahedral element augmented with cached geometry.
     * @tparam T1 Numeric type for angles and lengths.
     * @tparam T2 Integer node-index type.
     */
    template<typename T1, typename T2>
    struct tetrahedronElemAngle {
        T2 i[4];       ///< Indices of the four corner nodes.
        T1 a[4][3];    ///< Angles at each node (3 per face-node).
        T1 l[6];       ///< Lengths of the six edges.

        /// Constructs from a plain @ref tetrahedronElem, copying node indices.
        tetrahedronElemAngle(const tetrahedronElem<T2>& t) {
            i[0] = t.i[0];
            i[1] = t.i[1];
            i[2] = t.i[2];
            i[3] = t.i[3];
        }
    };

    /**
     * @brief A virtual node interpolated between two real nodes.
     * @tparam T Numeric type.
     * @tparam N Node type.
     */
    template<typename T, typename N>
    struct virtualNode {
        N *node1;   ///< First anchoring node.
        N *node2;   ///< Second anchoring node.
        T a[3];     ///< Interpolation coefficients/angles.
        T e[3];     ///< Edge-related terms.
    };


    /// 2x2 determinant of two @ref sxz vectors (the 2D cross product).
    template<typename T>
    T det(const sxz<T>& u, const sxz<T>& v) {
        return u.x*v.z - u.z*v.x;
    }

    /// 3x3 determinant of three @ref sxyz vectors taken as columns/rows.
    template<typename T>
    T det(const sxyz<T>& v1, const sxyz<T>& v2, const sxyz<T>& v3) {
        return -v1.z*v2.y*v3.x + v1.y*v2.z*v3.x + v1.z*v2.x*v3.y -
        v1.x*v2.z*v3.y - v1.y*v2.x*v3.z + v1.x*v2.y*v3.z;
    }

    /// Scalar triple product v1 &middot; (v2 &times; v3).
    template<typename T>
    T tripleScalar(const sxyz<T>& v1, const sxyz<T>& v2, const sxyz<T>& v3) {
        return v1.x*(v2.y*v3.z-v3.y*v2.z) - v1.y*(v2.x*v3.z-v3.x*v2.z) + v1.z*(v2.x*v3.y-v3.x*v2.y);
    }

    /// Determinant used for tetrahedron orientation/volume from four vertices.
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

    /// Dot product of two 2D @ref sxz vectors.
    template<typename T>
    T dot(const sxz<T>& v1, const sxz<T>& v2) {
        return v1.x*v2.x + v1.z*v2.z;
    }


    /// Dot product of two scalars, each treated as a vector along y.
    template<typename T>
    T dot(const T v1, const T v2) {
        // both v1 & v2 considered to be a vector along y
        return v1*v2;
    }

    /// 2D cross product (z-component) of two @ref sxz vectors.
    template<typename T>
    T cross(const sxz<T>& v1, const sxz<T>& v2) {
        return v1.z*v2.x - v1.x*v2.z;
    }

    /// Dot product of two 3D @ref sxyz vectors.
    template<typename T>
    T dot(const sxyz<T>& v1, const sxyz<T>& v2) {
        return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
    }

    /// Cross product v1 &times; v2 of two @ref sxyz vectors.
    template<typename T>
    sxyz<T> cross(const sxyz<T>& v1, const sxyz<T>& v2) {
        sxyz<T> v3;
        v3.x = v1.y*v2.z - v1.z*v2.y;
        v3.y = v1.z*v2.x - v1.x*v2.z;
        v3.z = v1.x*v2.y - v1.y*v2.x;
        return v3;
    }

    /// Euclidean (L2) norm of a 2D @ref sxz vector.
    template<typename T>
    T norm(const sxz<T>& v) {
        return sqrt( v.x*v.x + v.z*v.z );
    }

    /// Euclidean (L2) norm of a 3D @ref sxyz vector.
    template<typename T>
    T norm(const sxyz<T>& v) {
        return sqrt( v.x*v.x + v.y*v.y + v.z*v.z );
    }

    /// Squared Euclidean norm of a 2D @ref sxz vector (avoids the sqrt).
    template<typename T>
    T norm2(const sxz<T>& v) {
        return v.x*v.x + v.z*v.z;
    }

    /// Squared Euclidean norm of a 3D @ref sxyz vector (avoids the sqrt).
    template<typename T>
    T norm2(const sxyz<T>& v) {
        return v.x*v.x + v.y*v.y + v.z*v.z;
    }

    /**
     * @brief Projects a point onto a triangle plane in normalized coordinates.
     *
     * Computes the normalized coordinates (@p xi0, @p zeta0) of point P within
     * the plane of triangle ABC, following Lelièvre et al. (2011), GJI 184,
     * 885-896. Solves @f$ \xi\,c + \zeta\,b = p @f$ in a least-squares sense
     * via the normal equations @f$ A^T A x = A^T b @f$. If the system is
     * singular, both outputs are set to -1.
     *
     * @param b Unit vector along edge AC.
     * @param c Unit vector along edge AB.
     * @param p Vector AP from vertex A to the point P.
     * @param[out] xi0   Normalized coordinate along @p c.
     * @param[out] zeta0 Normalized coordinate along @p b.
     */
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
    /// Sign implementation for unsigned types: returns 1 if @p x > 0, else 0.
    template <typename T> inline constexpr
    int signum(T x, std::false_type is_signed) {
        return T(0) < x;
    }

    /// Sign implementation for signed types: returns -1, 0 or 1.
    template <typename T> inline constexpr
    int signum(T x, std::true_type is_signed) {
        return (T(0) < x) - (x < T(0));
    }

    /// Returns the sign of @p x (-1, 0 or 1), dispatching on signedness of @c T.
    template <typename T> inline constexpr
    int signum(T x) {
        return signum(x, std::is_signed<T>());
    }
#else
    /// Sign implementation for unsigned types (MSVC): returns 1 if @p x > 0, else 0.
    template <typename T>
    int signum(T x, std::false_type is_signed) {
        return T(0) < x;
    }

    /// Sign implementation for signed types (MSVC): returns -1, 0 or 1.
    template <typename T>
    int signum(T x, std::true_type is_signed) {
        return (T(0) < x) - (x < T(0));
    }

    /// Returns the sign of @p x (-1, 0 or 1), dispatching on signedness of @c T.
    template <typename T>
    int signum(T x) {
        return signum(x, std::is_signed<T>());
    }
#endif

    /**
     * @brief Tests whether points @p v4 and @p p lie on the same side of a plane.
     *
     * The plane is defined by the triangle (@p v1, @p v2, @p v3). Returns true
     * when @p v4 and @p p fall on the same side (or on the plane), determined
     * by comparing the sign of their projections onto the triangle normal.
     *
     * @param v1 First vertex defining the plane.
     * @param v2 Second vertex defining the plane.
     * @param v3 Third vertex defining the plane.
     * @param v4 Reference point.
     * @param p  Query point.
     * @return @c true if @p p is on the same side as @p v4.
     */
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

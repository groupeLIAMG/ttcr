//
//  Cell.h
//  ttcr
//
//  Created by Bernard Giroux on 16-02-27.
//  Copyright © 2016 Bernard Giroux. All rights reserved.
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
 * @file Cell.h
 * @brief Classes holding the physical properties of the cells of a grid
 *
 * The classes defined in this file are policy classes, passed as the `CELL`
 * template argument of the grid classes that assume a constant velocity per
 * cell (Grid2Drc*, Grid3Drc*, Grid2Duc*, Grid3Duc*).  Each class stores the
 * medium parameters of every cell of a grid, and provides the operations
 * needed by the raytracing algorithms to convert a ray segment into a
 * traveltime increment.
 *
 * The classes share an implicit (duck-typed) interface made of:
 * - a constructor taking the number of cells of the grid;
 * - a set of `setXxx()` methods used to assign the medium parameters.  The
 *   setters that are meaningless for a given class throw a std::logic_error,
 *   so that a mismatch between the requested physics and the parameters
 *   supplied by the caller is reported at run time;
 * - `computeDt()` overloads returning the traveltime needed to travel between
 *   two points located in a given cell;
 * - `computeDistance()` overloads storing, in a siv or siv2 struct, the length
 *   of a ray segment (used when building raypaths and sensitivity matrices).
 *
 * @note The interface is not implemented uniformly.  A setter that is simply
 * absent from a class yields a compilation error rather than the run-time
 * std::logic_error described above; `setChi()` and `setPsi()`, for instance,
 * exist only in Cell, CellWeaklyAnelliptical and CellElliptical3D.  The 3D
 * classes are the least complete: none of them provides `computeDistance()`,
 * and CellVTI_PSV3D, CellVTI_SH3D and CellWeaklyAnelliptical3D provide neither
 * `setSlowness()` nor `getSlowness()`.
 *
 * Anisotropy is described either with anisotropy ratios (elliptical classes)
 * or with Thomsen's (1986) parameters (VTI classes).  Classes whose name does
 * not end with `3D` are meant for 2D grids lying in the x-z plane; the y
 * dimension is ignored.
 *
 * @sa Grid2Drc, Grid3Drc, Grid2Duc, Grid3Duc
 */

#ifndef ttcr_Cell_h
#define ttcr_Cell_h

#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

#include "ttcr_t.h"

namespace ttcr {

    /**
     * @brief Isotropic cells, defined by their slowness
     *
     * Simplest of the cell classes: the traveltime across a cell is the
     * straight-ray distance times the slowness of the cell.
     *
     * @tparam T    type of the medium parameters and traveltimes (float or double)
     * @tparam NODE type of the grid nodes (Node2Dc, Node2Dcsp, Node3Dc, ...)
     * @tparam S    type of the coordinate struct (sxz<T> in 2D, sxyz<T> in 3D)
     */
    template <typename T, typename NODE, typename S>
    class Cell {
    public:
        /// Number of medium parameters of the cells: the slowness
        static constexpr size_t nParams = 1;
        /**
         * @brief Constructor
         * @param n number of cells of the grid
         */
        Cell(const size_t n) : slowness(std::vector<T>(n)) { }

        /**
         * @brief Assign the slowness of the cells
         * @param s slowness values, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setSlowness(const std::vector<T>& s) {
            if ( slowness.size() != s.size() ) {
                throw std::length_error("Error: slowness vectors of incompatible size.");
            }
            slowness = s;
        }

        /**
         * @brief Get the slowness of a cell
         * @param i index of the cell
         * @return slowness of cell @p i
         * @throws std::out_of_range if @p i is not a valid cell index
         */
        const T getSlowness(const size_t i) const {
            return slowness.at(i);
        }

        /// @brief Not applicable: always throws std::logic_error
        void setXi(const std::vector<T>& s) {
            throw std::logic_error("Error: xi not defined for Cell.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setChi(const std::vector<T>& s) {
            throw std::logic_error("Error: chi not defined for Cell.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setPsi(const std::vector<T>& s) {
            throw std::logic_error("Error: psi not defined for Cell.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setTiltAngle(const std::vector<T>& s) {
            throw std::logic_error("Error: TiltAngle not defined for Cell.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setVp0(const std::vector<T>& s) {
            throw std::logic_error("Error: Vp0 not defined for Cell.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setVs0(const std::vector<T>& s) {
            throw std::logic_error("Error: Vs0 not defined for Cell.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setDelta(const std::vector<T>& s) {
            throw std::logic_error("Error: delta not defined for Cell.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setEpsilon(const std::vector<T>& s) {
            throw std::logic_error("Error: epsilon not defined for Cell.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setGamma(const std::vector<T>& s) {
            throw std::logic_error("Error: gamma not defined for Cell.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setS2(const std::vector<T>& r) {
            throw std::logic_error("Error: s2 not defined for Cell.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setS4(const std::vector<T>& r) {
            throw std::logic_error("Error: s4 not defined for Cell.");
        }

        /**
         * @brief Traveltime between two points of a cell
         * @param source first point
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime, i.e. the distance times the slowness of the cell
         */
        T computeDt(const S& source, const S& node,
                    const size_t cellNo) const {
            return slowness[cellNo] * source.getDistance( node );
        }

        /**
         * @brief Traveltime between a node and a point of a cell
         * @param source node from which the ray originates
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime, i.e. the distance times the slowness of the cell
         */
        T computeDt(const NODE& source, const S& node,
                    const size_t cellNo) const {
            return slowness[cellNo] * source.getDistance( node );
        }

        /**
         * @brief Traveltime between two nodes of a cell
         * @param source node from which the ray originates
         * @param node   second node
         * @param cellNo index of the cell holding the two nodes
         * @return the traveltime, i.e. the distance times the slowness of the cell
         */
        T computeDt(const NODE& source, const NODE& node,
                    const size_t cellNo) const {
            return slowness[cellNo] * source.getDistance( node );
        }

        /**
         * @brief Length of a ray segment, stored in a siv struct
         * @param[in]  source node from which the ray segment originates
         * @param[in]  node   end point of the ray segment
         * @param[out] cell   struct in which the length is stored (member `v`)
         */
        void computeDistance(const NODE& source, const S& node,
                             siv<T>& cell) const {
            cell.v = source.getDistance( node );
        }
        /**
         * @brief Length of a ray segment, stored in a siv2 struct
         *
         * The medium being isotropic, the whole length is stored in member `v`
         * and member `v2` is set to zero.  Zeroing it is required: siv2 has no
         * default constructor, so callers that declare it without an
         * initializer (as the raytracing routines do) would otherwise
         * accumulate an indeterminate value through siv2::operator+=.
         *
         * @param[in]  source node from which the ray segment originates
         * @param[in]  node   end point of the ray segment
         * @param[out] cell   struct in which the length is stored (member `v`),
         *                    member `v2` being set to zero
         */
        void computeDistance(const NODE& source, const S& node,
                             siv2<T>& cell) const {
            cell.v = source.getDistance( node );
            cell.v2 = 0;
        }

    private:
        std::vector<T> slowness;  ///< slowness of the cells
    };


    /**
     * @brief Cells with elliptical anisotropy, in 2D (y dimension ignored)
     *
     * The medium is described by the horizontal slowness @f$s_x@f$ of the cell
     * and by the anisotropy ratio @f$\xi = s_z / s_x@f$, so that the traveltime
     * along a segment of components @f$(l_x, l_z)@f$ is
     * @f[ dt = s_x \sqrt{l_x^2 + \xi^2 l_z^2} @f]
     * The axes of the ellipse are aligned with the axes of the grid; see
     * CellTiltedElliptical for the tilted case.
     *
     * @tparam T    type of the medium parameters and traveltimes (float or double)
     * @tparam NODE type of the grid nodes (Node2Dc, Node2Dcsp, ...)
     * @tparam S    type of the coordinate struct (sxz<T>)
     */
    template <typename T, typename NODE, typename S>
    class CellElliptical {
    public:
        /// Number of medium parameters of the cells: the slowness and the anisotropy ratio
        static constexpr size_t nParams = 2;
        /**
         * @brief Constructor
         * @param n number of cells of the grid
         */
        CellElliptical(const size_t n) :
        slowness(std::vector<T>(n)),
        xi(std::vector<T>(n,0.0)) {
        }

        /**
         * @brief Assign the horizontal slowness of the cells
         * @param s slowness values, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setSlowness(const std::vector<T>& s) {
            if ( slowness.size() != s.size() ) {
                throw std::length_error("Error: slowness vectors of incompatible size.");
            }
            slowness = s;
        }

        /**
         * @brief Get the horizontal slowness of a cell
         * @param i index of the cell
         * @return slowness of cell @p i
         * @throws std::out_of_range if @p i is not a valid cell index
         */
        const T getSlowness(const size_t i) const {
            return slowness.at(i);
        }

        /**
         * @brief Assign the anisotropy ratio of the cells
         * @param s values of @f$\xi = s_z / s_x@f$, one per cell; they are
         *          squared and stored as such
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setXi(const std::vector<T>& s) {
            if ( xi.size() != s.size() ) {
                throw std::length_error("Error: xi vectors of incompatible size.");
            }
            for ( size_t n=0; n<xi.size(); ++n ) {
                xi[n] = s[n]*s[n];
            }
        }

        /// @brief Not applicable: always throws std::logic_error
        void setTiltAngle(const std::vector<T>& s) {
            throw std::logic_error("Error: TiltAngle not defined for CellElliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setVp0(const std::vector<T>& s) {
            throw std::logic_error("Error: Vp0 not defined for CellElliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setVs0(const std::vector<T>& s) {
            throw std::logic_error("Error: Vs0 not defined for CellElliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setDelta(const std::vector<T>& s) {
            throw std::logic_error("Error: delta not defined for CellElliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setEpsilon(const std::vector<T>& s) {
            throw std::logic_error("Error: epsilon not defined for CellElliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setGamma(const std::vector<T>& s) {
            throw std::logic_error("Error: gamma not defined for CellElliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setS2(const std::vector<T>& r) {
            throw std::logic_error("Error: s2 not defined for CellElliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setS4(const std::vector<T>& r) {
            throw std::logic_error("Error: s4 not defined for CellElliptical.");
        }

        /**
         * @brief Traveltime between two points of a cell
         * @param source first point
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime along the segment joining the two points
         */
        T computeDt(const S& source, const S& node,
                    const size_t cellNo) const {
            T lx = node.x - source.x;
            T lz = node.z - source.z;
            return slowness[cellNo] * std::sqrt( lx*lx + xi[cellNo]*lz*lz );
        }

        /**
         * @brief Traveltime between a node and a point of a cell
         * @param source node from which the ray originates
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime along the segment joining the two points
         */
        T computeDt(const NODE& source, const S& node,
                    const size_t cellNo) const {
            T lx = node.x - source.getX();
            T lz = node.z - source.getZ();
            return slowness[cellNo] * std::sqrt( lx*lx + xi[cellNo]*lz*lz );
        }

        /**
         * @brief Traveltime between two nodes of a cell
         * @param source node from which the ray originates
         * @param node   second node
         * @param cellNo index of the cell holding the two nodes
         * @return the traveltime along the segment joining the two nodes
         */
        T computeDt(const NODE& source, const NODE& node,
                    const size_t cellNo) const {
            T lx = node.getX() - source.getX();
            T lz = node.getZ() - source.getZ();
            return slowness[cellNo] * std::sqrt( lx*lx + xi[cellNo]*lz*lz );
        }

        /**
         * @brief Not applicable: always throws std::logic_error
         *
         * A single length cannot describe the sensitivity of an anisotropic
         * cell; the siv2 overload must be used instead.  This condition was
         * formerly reported with assert(), which is compiled out whenever
         * NDEBUG is defined.
         *
         * @param[in]  source unused
         * @param[in]  node   unused
         * @param[out] cell   unused
         * @throws std::logic_error always
         */
        void computeDistance(const NODE& source, const S& node,
                             siv<T>& cell) const {
            throw std::logic_error("Error: computeDistance with siv not defined for CellElliptical.");
        }
        /**
         * @brief Sensitivity of a segment traveltime to the medium parameters
         *
         * With @f$ dt = s\sqrt{l_x^2 + \xi^2 l_z^2} @f$,
         * @f[ \frac{\partial dt}{\partial s} = \frac{dt}{s}, \qquad
         *     \frac{\partial dt}{\partial \xi} = \frac{s^2 \xi l_z^2}{dt} @f]
         *
         * @param[in]  source node from which the ray segment originates
         * @param[in]  node   end point of the ray segment
         * @param[out] cell   derivatives w/r to the slowness along the fast
         *                    axis in member `v`, and to the anisotropy ratio
         *                    in member `v2`
         */
        void computeDistance(const NODE& source, const S& node,
                             siv2<T>& cell) const {
            const T lx = node.x - source.getX();
            const T lz = node.z - source.getZ();
            const T q  = std::sqrt( lx*lx + xi[cell.i]*lz*lz );
            const T s  = slowness[cell.i];
            const T dt = s*q;
            cell.v  = q;
            cell.v2 = ( dt > 0. ) ?
                      s*s*std::sqrt(xi[cell.i])*lz*lz/dt : T(0);
        }

    private:
        std::vector<T> slowness;  ///< horizontal slowness of the cells
        std::vector<T> xi;        ///< anisotropy ratio, xi = sz / sx, *** squared ***
    };


    /**
     * @brief Cells with tilted elliptical anisotropy, in 2D (y dimension ignored)
     *
     * Same parametrization as CellElliptical, with the axes of the ellipse
     * rotated by the tilt angle @f$\theta_t@f$ of the cell.  The components
     * @f$(l_x, l_z)@f$ of a ray segment are first rotated in the frame of the
     * ellipse,
     * @f[ t_1 = l_x \cos\theta_t + l_z \sin\theta_t, \quad
     *     t_2 = l_z \cos\theta_t - l_x \sin\theta_t @f]
     * and the traveltime is then @f$ dt = s \sqrt{t_1^2 + \xi^2 t_2^2} @f$.
     *
     * @tparam T    type of the medium parameters and traveltimes (float or double)
     * @tparam NODE type of the grid nodes (Node2Dc, Node2Dcsp, ...)
     * @tparam S    type of the coordinate struct (sxz<T>)
     */
    template <typename T, typename NODE, typename S>
    class CellTiltedElliptical {
    public:
        /// Number of medium parameters of the cells: the slowness, the anisotropy ratio and the tilt angle
        static constexpr size_t nParams = 3;
        /**
         * @brief Constructor
         * @param n number of cells of the grid
         */
        CellTiltedElliptical(const size_t n) :
        slowness(std::vector<T>(n)),
        xi(std::vector<T>(n,0.0)),
        tAngle(std::vector<T>(n,0.0)),
        ca(std::vector<T>(n,1.0)),
        sa(std::vector<T>(n,0.0)) {
        }

        /**
         * @brief Assign the slowness of the cells, along the fast axis of the ellipse
         * @param s slowness values, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setSlowness(const std::vector<T>& s) {
            if ( slowness.size() != s.size() ) {
                throw std::length_error("Error: slowness vectors of incompatible size.");
            }
            slowness = s;
        }

        /**
         * @brief Get the slowness of a cell
         * @param i index of the cell
         * @return slowness of cell @p i
         * @throws std::out_of_range if @p i is not a valid cell index
         */
        const T getSlowness(const size_t i) const {
            return slowness.at(i);
        }

        /**
         * @brief Assign the anisotropy ratio of the cells
         * @param s values of @f$\xi = s_z / s_x@f$, one per cell; they are
         *          squared and stored as such
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setXi(const std::vector<T>& s) {
            if ( xi.size() != s.size() ) {
                throw std::length_error("Error: xi vectors of incompatible size.");
            }
            for ( size_t n=0; n<xi.size(); ++n ) {
                xi[n] = s[n]*s[n];
            }
        }

        /**
         * @brief Assign the tilt angle of the cells
         *
         * The sine and cosine of the angles are precomputed and cached.
         *
         * @param s tilt angles in radians, measured from the z axis, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setTiltAngle(const std::vector<T>& s) {
            if ( tAngle.size() != s.size() ) {
                throw std::length_error("Error: angle vectors of incompatible size.");
            }
            for ( size_t n=0; n<tAngle.size(); ++n ) {
                tAngle[n] = s[n];
                ca[n] = std::cos(s[n]);
                sa[n] = std::sin(s[n]);
            }
        }

        /// @brief Not applicable: always throws std::logic_error
        void setVp0(const std::vector<T>& s) {
            throw std::logic_error("Error: Vp0 not defined for CellTiltedElliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setVs0(const std::vector<T>& s) {
            throw std::logic_error("Error: Vs0 not defined for CellTiltedElliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setDelta(const std::vector<T>& s) {
            throw std::logic_error("Error: delta not defined for CellTiltedElliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setEpsilon(const std::vector<T>& s) {
            throw std::logic_error("Error: epsilon not defined for CellTiltedElliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setGamma(const std::vector<T>& s) {
            throw std::logic_error("Error: gamma not defined for CellTiltedElliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setS2(const std::vector<T>& r) {
            throw std::logic_error("Error: s2 not defined for CellTiltedElliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setS4(const std::vector<T>& r) {
            throw std::logic_error("Error: s4 not defined for CellTiltedElliptical.");
        }

        /**
         * @brief Traveltime between two points of a cell
         * @param source first point
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime along the segment joining the two points
         */
        T computeDt(const S& source, const S& node,
                    const size_t cellNo) const {
            T lx = node.x - source.x;
            T lz = node.z - source.z;
            T t1 = lx * ca[cellNo] + lz * sa[cellNo];
            T t2 = lz * ca[cellNo] - lx * sa[cellNo];

            return slowness[cellNo] * std::sqrt( t1*t1 + xi[cellNo]*t2*t2 );
        }

        /**
         * @brief Traveltime between a node and a point of a cell
         * @param source node from which the ray originates
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime along the segment joining the two points
         */
        T computeDt(const NODE& source, const S& node,
                    const size_t cellNo) const {
            T lx = node.x - source.getX();
            T lz = node.z - source.getZ();
            T t1 = lx * ca[cellNo] + lz * sa[cellNo];
            T t2 = lz * ca[cellNo] - lx * sa[cellNo];

            return slowness[cellNo] * std::sqrt( t1*t1 + xi[cellNo]*t2*t2 );
        }

        /**
         * @brief Traveltime between two nodes of a cell
         * @param source node from which the ray originates
         * @param node   second node
         * @param cellNo index of the cell holding the two nodes
         * @return the traveltime along the segment joining the two nodes
         */
        T computeDt(const NODE& source, const NODE& node,
                    const size_t cellNo) const {
            T lx = node.getX() - source.getX();
            T lz = node.getZ() - source.getZ();
            T t1 = lx * ca[cellNo] + lz * sa[cellNo];
            T t2 = lz * ca[cellNo] - lx * sa[cellNo];

            return slowness[cellNo] * std::sqrt( t1*t1 + xi[cellNo]*t2*t2 );
        }

        /**
         * @brief Not applicable: always throws std::logic_error
         *
         * A single length cannot describe the sensitivity of an anisotropic
         * cell; the siv2 overload must be used instead.  This condition was
         * formerly reported with assert(), which is compiled out whenever
         * NDEBUG is defined.
         *
         * @param[in]  source unused
         * @param[in]  node   unused
         * @param[out] cell   unused
         * @throws std::logic_error always
         */
        void computeDistance(const NODE& source, const S& node,
                             siv<T>& cell) const {
            throw std::logic_error("Error: computeDistance with siv not defined for CellTiltedElliptical.");
        }
        /**
         * @brief Components of a ray segment, stored in a siv2 struct
         * @param[in]  source node from which the ray segment originates
         * @param[in]  node   end point of the ray segment
         * @param[out] cell   struct holding the horizontal component of the
         *                    segment in member `v`, and its vertical component
         *                    in member `v2`
         */
        void computeDistance(const NODE& source, const S& node,
                             siv2<T>& cell) const {
            cell.v  = std::abs(node.x - source.getX());
            cell.v2 = std::abs(node.z - source.getZ());
        }

        /**
         * @brief Sensitivity of a segment traveltime to the medium parameters
         *
         * With @f$ t_1 = l_x\cos\theta + l_z\sin\theta @f$,
         * @f$ t_2 = l_z\cos\theta - l_x\sin\theta @f$ and
         * @f$ dt = s\sqrt{t_1^2 + \xi^2 t_2^2} @f$,
         * @f[ \frac{\partial dt}{\partial s} = \frac{dt}{s}, \quad
         *     \frac{\partial dt}{\partial \xi} = \frac{s^2 \xi t_2^2}{dt},
         *     \quad
         *     \frac{\partial dt}{\partial \theta}
         *       = \frac{s^2 t_1 t_2 (1-\xi^2)}{dt} @f]
         *
         * @param[in]  source node from which the ray segment originates
         * @param[in]  node   end point of the ray segment
         * @param[out] cell   derivatives w/r to the slowness, the anisotropy
         *                    ratio and the tilt angle, in members `v`, `v2`
         *                    and `v3`; member `v4` is left at zero
         */
        void computeDistance(const NODE& source, const S& node,
                             siv4<T>& cell) const {
            const T lx = node.x - source.getX();
            const T lz = node.z - source.getZ();
            const T t1 = lx*ca[cell.i] + lz*sa[cell.i];
            const T t2 = lz*ca[cell.i] - lx*sa[cell.i];
            const T x2 = xi[cell.i];              // xi is stored squared
            const T q  = std::sqrt( t1*t1 + x2*t2*t2 );
            const T s  = slowness[cell.i];
            const T dt = s*q;
            cell.v  = q;
            cell.v2 = ( dt > 0. ) ? s*s*std::sqrt(x2)*t2*t2/dt : T(0);
            cell.v3 = ( dt > 0. ) ? s*s*t1*t2*(1.-x2)/dt : T(0);
            cell.v4 = T(0);
        }

    private:
        std::vector<T> slowness;  ///< slowness of the cells along the fast axis
        std::vector<T> xi;        ///< anisotropy ratio, xi = sz / sx, *** squared ***
        std::vector<T> tAngle;    ///< column-wise (z axis) anisotropy angle of the cells, in radians
        std::vector<T> ca;        ///< cosine of tAngle
        std::vector<T> sa;        ///< sine of tAngle
    };


    /**
     * @brief Group (energy) velocity tables for VTI qP/qSV media
     *
     * In anisotropic media the traveltime along a ray segment is its length
     * divided by the **group** (energy) velocity taken in the direction of the
     * segment, not by the phase velocity: the two coincide only along the
     * symmetry directions.  Unlike the elliptical and weakly anelliptical
     * cases, the group velocity of the coupled qP/qSV system has no closed
     * form, so it is tabulated.
     *
     * For a phase angle @f$\theta@f$ measured from the vertical symmetry axis,
     * the group velocity magnitude and the ray (group) angle are
     * @f[ V_g = \sqrt{v^2 + \left(\frac{dv}{d\theta}\right)^2},
     *     \qquad
     *     \psi = \theta + \arctan\left(\frac{1}{v}\frac{dv}{d\theta}\right) @f]
     * The parametric curve @f$(\psi, V_g)@f$ is sampled densely in
     * @f$\theta@f$ and scattered into a table indexed by @f$\psi@f$, keeping
     * the *largest* @f$V_g@f$ of every bin.  Where the qSV group-velocity
     * surface is multivalued (triplications) this automatically selects the
     * first arrival, without any explicit branch handling.
     *
     * Tables are shared by all the cells having the same
     * @f$(V_{P0}, V_{S0}, \epsilon, \delta)@f$, so layered and blocky models
     * hold only a handful of them.  A model in which every cell has distinct
     * parameters holds one table per cell, costing 901 values per cell.
     *
     * @tparam T type of the medium parameters and traveltimes (float or double)
     */
    template <typename T>
    class VTI_PSV_GroupVel {
    public:
        /**
         * @brief Build one table per distinct set of medium parameters
         * @param Vp0     vertical P-wave velocity of the cells
         * @param Vs0     vertical S-wave velocity of the cells
         * @param epsilon Thomsen's parameter epsilon of the cells
         * @param delta   Thomsen's parameter delta of the cells
         * @param sign    +1 for the qP phase, -1 for the qSV phase
         */
        void build(const std::vector<T>& Vp0, const std::vector<T>& Vs0,
                   const std::vector<T>& epsilon, const std::vector<T>& delta,
                   const T sign) {
            const size_t n = Vp0.size();
            sgn = sign;
            tables.clear();
            index.assign(n, 0);
            std::map<std::array<T,4>, size_t> seen;
            for ( size_t i=0; i<n; ++i ) {
                const std::array<T,4> key = {{ Vp0[i], Vs0[i], epsilon[i], delta[i] }};
                const typename std::map<std::array<T,4>, size_t>::const_iterator it
                    = seen.find(key);
                if ( it == seen.end() ) {
                    index[i] = tables.size();
                    seen[key] = index[i];
                    tables.push_back( tabulate(Vp0[i], Vs0[i], epsilon[i],
                                               delta[i], sign) );
                } else {
                    index[i] = it->second;
                }
            }
        }

        /// @brief True if no table has been built yet
        bool empty() const { return tables.empty(); }

        /// @brief Number of distinct tables held
        size_t size() const { return tables.size(); }

        /**
         * @brief Group velocity in the direction of a ray segment
         * @param lh     horizontal component of the segment
         * @param lz     vertical component of the segment
         * @param cellNo index of the cell holding the segment
         * @return the group (energy) velocity of the first arrival
         * @throws std::logic_error if the medium parameters were not all set
         */
        T velocity(const T lh, const T lz, const size_t cellNo) const {
            const Table& t = table(cellNo);
            size_t i;
            T w;
            locate(t, lh, lz, i, w);
            return (1.-w)*t.vg[i] + w*t.vg[i+1];
        }

        /**
         * @brief Sensitivity of a segment traveltime to the medium parameters
         *
         * The traveltime across a segment of length @f$\ell@f$ making an angle
         * @f$\psi@f$ with the symmetry axis is @f$\ell / V_g(\psi)@f$.  Writing
         * the group velocity through the phase velocity of the branch carrying
         * the arrival, @f$V_g = v / \cos(\psi - \theta_s)@f$, and noting that
         * @f$\theta_s@f$ is stationary for @f$\cos(\psi-\theta)/v(\theta)@f$,
         * the envelope theorem removes the derivative of @f$\theta_s@f$ with
         * respect to @f$p@f$ and leaves
         * @f[ \frac{\partial}{\partial p}\frac{\ell}{V_g}
         *     = -\frac{\ell}{V_g\,v}\,\frac{\partial v}{\partial p}
         *     \qquad \text{evaluated at } \theta_s @f]
         * which holds on whichever branch is selected and is therefore valid
         * through the triplications of the qSV wave.
         *
         * @param[in]  lh     horizontal component of the segment
         * @param[in]  lz     vertical component of the segment
         * @param[in]  cellNo index of the cell holding the segment
         * @param[out] out    the four derivatives, w/r to @f$V_{P0}@f$,
         *                    @f$V_{S0}@f$, @f$\epsilon@f$ and @f$\delta@f$
         * @throws std::logic_error if the medium parameters were not all set
         */
        void sensitivity(const T lh, const T lz, const size_t cellNo,
                         T* out) const {
            const Table& t = table(cellNo);
            size_t i;
            T w;
            locate(t, lh, lz, i, w);
            const T vg = (1.-w)*t.vg[i] + w*t.vg[i+1];
            // the branch may change from one sample to the next at a
            // triplication cusp; interpolating across such a step would
            // describe neither branch, so the nearer sample is taken instead
            const T th = ( std::abs(t.th[i+1]-t.th[i]) < branchStep() ) ?
                         (1.-w)*t.th[i] + w*t.th[i+1] :
                         ( w < 0.5 ? t.th[i] : t.th[i+1] );
            T v, dvdVp0, dvdVs0, dvdEps, dvdDlt;
            phaseVelocityDerivatives(th, t.Vp0, t.Vs0, t.eps, t.dlt, sgn,
                                     v, dvdVp0, dvdVs0, dvdEps, dvdDlt);
            const T ell = std::sqrt(lh*lh + lz*lz);
            const T c = ( vg > 0. && v > 0. ) ? -ell/(vg*v) : T(0);
            out[0] = c*dvdVp0;
            out[1] = c*dvdVs0;
            out[2] = c*dvdEps;
            out[3] = c*dvdDlt;
        }

        /**
         * @brief Sensitivity of a segment traveltime to the ray angle
         *
         * @f$\partial(\ell/V_g)/\partial\psi = -\ell V_g^{-2}\,dV_g/d\psi@f$.
         * Tilting the medium by @f$\theta_t@f$ shifts the angle the segment
         * makes with the symmetry axis by the same amount, so this is also the
         * derivative with respect to the tilt angle of a tilted cell.
         *
         * @param lh     horizontal component of the segment
         * @param lz     vertical component of the segment
         * @param cellNo index of the cell holding the segment
         * @return the derivative of the segment traveltime w/r to the angle
         *         @f$\psi@f$ the segment makes with the symmetry axis, signed
         *         as @f$\mathrm{atan2}(l_h, l_z)@f$
         * @throws std::logic_error if the medium parameters were not all set
         *
         * @note @f$dV_g/d\psi@f$ is obtained by differencing the tabulated
         * group velocity, so its accuracy degrades at the cusps bounding a qSV
         * triplication, where @f$V_g@f$ has a corner.
         */
        T dt_dpsi(const T lh, const T lz, const size_t cellNo) const {
            const Table& t = table(cellNo);
            size_t i;
            T w;
            locate(t, lh, lz, i, w);
            const T vg  = (1.-w)*t.vg[i]  + w*t.vg[i+1];
            const T dvg = (1.-w)*t.dvg[i] + w*t.dvg[i+1];
            if ( vg <= 0. ) return T(0);
            const T ell = std::sqrt(lh*lh + lz*lz);
            // the tables span [0, pi/2]; folding the angle of the segment into
            // that interval reverses the sense in which it varies whenever the
            // two components have opposite signs
            const T fold = ( (lh < 0.) == (lz < 0.) ) ? T(1) : T(-1);
            return -fold*ell*dvg/(vg*vg);
        }

        /**
         * @brief Phase velocity and its derivative w/r to the phase angle
         * @param[in]  theta phase angle measured from the vertical axis
         * @param[in]  Vp0   vertical P-wave velocity
         * @param[in]  Vs0   vertical S-wave velocity
         * @param[in]  eps   Thomsen's parameter epsilon
         * @param[in]  dlt   Thomsen's parameter delta
         * @param[in]  sign  +1 for the qP phase, -1 for the qSV phase
         * @param[out] v     phase velocity
         * @param[out] dv    derivative of the phase velocity w/r to theta
         */
        static void phaseVelocity(const T theta, const T Vp0, const T Vs0,
                                  const T eps, const T dlt, const T sign,
                                  T& v, T& dv) {
            const T f  = 1. - (Vs0*Vs0)/(Vp0*Vp0);
            const T st = std::sin(theta);
            const T s  = st*st;
            const T A  = 1. + 2.*eps*s/f;
            const T R  = A*A - 8.*(eps-dlt)*s*(1.-s)/f;
            const T sqR = std::sqrt( R > 0. ? R : 0. );
            const T G  = 1. + eps*s - f/2. + sign*(f/2.)*sqR;
            const T dRds = 2.*A*(2.*eps/f) - 8.*(eps-dlt)*(1.-2.*s)/f;
            const T dGds = eps + ( sqR > 0. ? sign*(f/4.)*dRds/sqR : 0. );
            const T sqG = std::sqrt( G > 0. ? G : 0. );
            v  = Vp0*sqG;
            dv = ( sqG > 0. ? Vp0*dGds*std::sin(2.*theta)/(2.*sqG) : 0. );
        }

        /**
         * @brief Phase velocity and its derivatives w/r to the medium parameters
         *
         * Obtained by differentiating the same expression as phaseVelocity(),
         * @f$v = V_{P0}\sqrt{G}@f$, the dependence on the vertical velocities
         * running through @f$f = 1 - V_{S0}^2/V_{P0}^2@f$.
         *
         * @param[in]  theta  phase angle measured from the vertical axis
         * @param[in]  Vp0    vertical P-wave velocity
         * @param[in]  Vs0    vertical S-wave velocity
         * @param[in]  eps    Thomsen's parameter epsilon
         * @param[in]  dlt    Thomsen's parameter delta
         * @param[in]  sign   +1 for the qP phase, -1 for the qSV phase
         * @param[out] v      phase velocity
         * @param[out] dvdVp0 derivative of the phase velocity w/r to Vp0
         * @param[out] dvdVs0 derivative of the phase velocity w/r to Vs0
         * @param[out] dvdEps derivative of the phase velocity w/r to epsilon
         * @param[out] dvdDlt derivative of the phase velocity w/r to delta
         */
        static void phaseVelocityDerivatives(const T theta, const T Vp0,
                                             const T Vs0, const T eps,
                                             const T dlt, const T sign,
                                             T& v, T& dvdVp0, T& dvdVs0,
                                             T& dvdEps, T& dvdDlt) {
            const T f  = 1. - (Vs0*Vs0)/(Vp0*Vp0);
            const T st = std::sin(theta);
            const T s  = st*st;
            const T A  = 1. + 2.*eps*s/f;
            const T R  = A*A - 8.*(eps-dlt)*s*(1.-s)/f;
            const T sqR = std::sqrt( R > 0. ? R : 0. );
            const T G  = 1. + eps*s - f/2. + sign*(f/2.)*sqR;
            const T sqG = std::sqrt( G > 0. ? G : 0. );
            v = Vp0*sqG;
            if ( sqG <= 0. || sqR <= 0. ) {
                dvdVp0 = dvdVs0 = dvdEps = dvdDlt = T(0);
                return;
            }
            const T dRdEps = 2.*A*(2.*s/f) - 8.*s*(1.-s)/f;
            const T dRdDlt = 8.*s*(1.-s)/f;
            const T dRdf   = 2.*A*(-2.*eps*s/(f*f))
                           + 8.*(eps-dlt)*s*(1.-s)/(f*f);
            const T dGdEps = s + sign*(f/4.)*dRdEps/sqR;
            const T dGdDlt =     sign*(f/4.)*dRdDlt/sqR;
            const T dGdf   = -0.5 + sign*0.5*sqR + sign*(f/4.)*dRdf/sqR;
            const T dfdVp0 =  2.*Vs0*Vs0/(Vp0*Vp0*Vp0);
            const T dfdVs0 = -2.*Vs0/(Vp0*Vp0);
            dvdVp0 = sqG + Vp0*dGdf*dfdVp0/(2.*sqG);
            dvdVs0 =       Vp0*dGdf*dfdVs0/(2.*sqG);
            dvdEps =       Vp0*dGdEps/(2.*sqG);
            dvdDlt =       Vp0*dGdDlt/(2.*sqG);
        }

    private:
        /// Number of samples of a table, spanning [0, pi/2]
        static constexpr size_t nSamples = 901;
        /// Oversampling of the phase angle when building a table
        static constexpr size_t oversampling = 16;

        /**
         * @brief Group velocity of one medium, tabulated against the ray angle
         *
         * The three vectors span @f$[0, \pi/2]@f$, the VTI symmetries giving
         * the remaining quadrants.  The medium parameters are kept because the
         * sensitivities are formed from the phase velocity of the branch
         * carrying the arrival.
         */
        struct Table {
            std::vector<T> ps;   ///< ray angle of the sample kept in each bin
            std::vector<T> vg;   ///< group velocity of the first arrival
            std::vector<T> th;   ///< phase angle of the branch, in [0, pi/2]
            std::vector<T> dvg;  ///< derivative of vg w/r to the ray angle
            T Vp0;               ///< vertical P-wave velocity
            T Vs0;               ///< vertical S-wave velocity
            T eps;               ///< Thomsen's parameter epsilon
            T dlt;               ///< Thomsen's parameter delta
        };

        std::vector<Table> tables;  ///< one table per distinct medium
        std::vector<size_t> index;  ///< table used by each cell
        T sgn;                      ///< +1 for the qP phase, -1 for the qSV phase

        static T halfPi() { return static_cast<T>(1.57079632679489661923); }
        static T pi()     { return static_cast<T>(3.14159265358979323846); }

        /// @brief Largest step of the phase angle treated as a single branch
        static T branchStep() { return static_cast<T>(0.05); }

        /// @brief Table of a cell
        /// @throws std::logic_error if the medium parameters were not all set
        const Table& table(const size_t cellNo) const {
            if ( tables.empty() ) {
                throw std::logic_error("Error: medium parameters of VTI_PSV cells not set.");
            }
            return tables[ index[cellNo] ];
        }

        /**
         * @brief Sample bracketing the direction of a segment, and its weight
         *
         * Interpolation uses the ray angles actually sampled rather than the
         * centres of the bins holding them: the sample kept in a bin is the one
         * of largest group velocity, which lies towards the edge of the bin the
         * group velocity increases to, and assuming it sits at the centre
         * distorts the interpolated slope.
         */
        static void locate(const Table& t, const T lh, const T lz,
                           size_t& i, T& w) {
            // VTI symmetry: the tables span [0, pi/2] only
            const T psi = std::atan2( std::abs(lh), std::abs(lz) );
            const T x = psi / halfPi() * (nSamples-1);
            if ( !(x > 0.) ) {
                i = 0;
            } else if ( x >= nSamples-1 ) {
                i = nSamples-2;
            } else {
                i = static_cast<size_t>(x);
            }
            // the sampled angles are ordered but need not bracket psi from the
            // bin index alone
            while ( i+1 < nSamples-1 && t.ps[i+1] < psi ) ++i;
            while ( i > 0 && t.ps[i] > psi ) --i;
            const T d = t.ps[i+1] - t.ps[i];
            w = ( d > 0. ) ? (psi - t.ps[i])/d : T(0);
            if ( w < 0. ) w = T(0);
            if ( w > 1. ) w = T(1);
        }

        /// @brief Tabulate the group velocity over [0, pi/2] for one medium
        static Table tabulate(const T Vp0, const T Vs0, const T eps,
                              const T dlt, const T sign) {
            Table t;
            t.Vp0 = Vp0;
            t.Vs0 = Vs0;
            t.eps = eps;
            t.dlt = dlt;
            t.vg.assign(nSamples, T(0));
            t.th.assign(nSamples, T(0));
            t.ps.assign(nSamples, T(0));
            const size_t m = (nSamples-1)*oversampling + 1;
            for ( size_t k=0; k<m; ++k ) {
                const T theta = pi() * static_cast<T>(k) / static_cast<T>(m-1);
                T v, dv;
                phaseVelocity(theta, Vp0, Vs0, eps, dlt, sign, v, dv);
                if ( v <= 0. ) continue;
                const T vg = std::sqrt(v*v + dv*dv);
                // fold the ray angle into [0, pi/2] using the VTI symmetries
                T psi = std::abs( theta + std::atan2(dv, v) );
                if ( psi > halfPi() ) psi = pi() - psi;
                psi = std::abs(psi);
                size_t j = static_cast<size_t>( psi/halfPi()*(nSamples-1) + 0.5 );
                if ( j >= nSamples ) j = nSamples-1;
                if ( vg > t.vg[j] ) {          // keep the first arrival
                    t.vg[j] = vg;
                    t.ps[j] = psi;
                    // the phase velocity and its derivatives depend on the
                    // phase angle only through sin^2, so folding it likewise
                    // loses nothing and keeps neighbouring samples comparable
                    t.th[j] = ( theta > halfPi() ) ? pi() - theta : theta;
                }
            }
            fillGaps(t.vg, t.th, t.ps);
            // derivative of the group velocity w/r to the ray angle; Vg is
            // symmetric about both ends of the interval, so it vanishes there
            t.dvg.assign(nSamples, T(0));
            for ( size_t i=1; i+1<nSamples; ++i ) {
                const T d = t.ps[i+1] - t.ps[i-1];
                if ( d > 0. ) t.dvg[i] = (t.vg[i+1] - t.vg[i-1])/d;
            }
            return t;
        }

        /**
         * @brief Fill the bins that no sample reached
         *
         * The group velocity is interpolated linearly; the phase angle is
         * copied from the nearer filled bin, a gap being as likely as not to
         * straddle a change of branch.
         */
        static void fillGaps(std::vector<T>& vg, std::vector<T>& th,
                             std::vector<T>& ps) {
            const T dpsi = halfPi()/static_cast<T>(nSamples-1);
            size_t lo = 0;
            while ( lo < vg.size() && vg[lo] == T(0) ) ++lo;
            if ( lo == vg.size() ) {
                throw std::runtime_error("Error: could not tabulate the group "
                                         "velocity of a VTI_PSV cell.");
            }
            for ( size_t i=0; i<lo; ++i ) {
                vg[i] = vg[lo];
                th[i] = th[lo];
                ps[i] = static_cast<T>(i)*dpsi;
            }
            size_t hi = vg.size()-1;
            while ( vg[hi] == T(0) ) --hi;
            for ( size_t i=hi+1; i<vg.size(); ++i ) {
                vg[i] = vg[hi];
                th[i] = th[hi];
                ps[i] = static_cast<T>(i)*dpsi;
            }
            size_t i = lo;
            while ( i < hi ) {
                if ( vg[i+1] != T(0) ) { ++i; continue; }
                size_t j = i+1;
                while ( vg[j] == T(0) ) ++j;
                for ( size_t k=i+1; k<j; ++k ) {
                    const T w = static_cast<T>(k-i)/static_cast<T>(j-i);
                    vg[k] = (1.-w)*vg[i] + w*vg[j];
                    th[k] = ( w < 0.5 ) ? th[i] : th[j];
                    ps[k] = static_cast<T>(k)*dpsi;
                }
                i = j;
            }
        }
    };


    /**
     * @brief Cells with VTI anisotropy, P or SV phase, in 2D (y dimension ignored)
     *
     * The medium is described by the vertical velocities @f$V_{P0}@f$ and
     * @f$V_{S0}@f$ and by Thomsen's (1986) parameters @f$\epsilon@f$ and
     * @f$\delta@f$.  The phase velocity of the qP (@f$\pm = +@f$) and qSV
     * (@f$\pm = -@f$) waves, for a phase angle @f$\theta@f$ measured from the
     * vertical (symmetry) axis, is
     * @f[ v(\theta) = V_{P0} \left[ 1 + \epsilon \sin^2\theta - \frac{f}{2}
     *     \pm \frac{f}{2} \sqrt{ \left( 1 + \frac{2 \epsilon \sin^2\theta}{f}
     *     \right)^2 - \frac{2 (\epsilon - \delta) \sin^2 2\theta}{f} }
     *     \right]^{1/2} @f]
     * with @f$ f = 1 - V_{S0}^2 / V_{P0}^2 @f$.
     *
     * @pre Every cell satisfies @f$V_{P0} > 0@f$ and @f$V_{S0} \neq V_{P0}@f$.  The
     * values are not validated: @f$V_{P0} = 0@f$ or @f$V_{S0} = V_{P0}@f$
     * (which makes @f$f = 0@f$) produces a division by zero.
     *
     * The traveltime along a ray segment is its length divided by the
     * **group** (energy) velocity taken in the direction of the segment.  That
     * velocity has no closed form here and is tabulated by VTI_PSV_GroupVel,
     * which also selects the first arrival where the qSV group-velocity surface
     * is multivalued.  The phase velocity must not be used in its place: the
     * two coincide only along the symmetry directions, and using the phase
     * velocity underestimates the traveltime by up to about 9% for qSV in
     * strongly anisotropic media.
     *
     * The tables are (re)built as soon as the four medium parameters have been
     * assigned, and again whenever one of them or the phase is changed.
     *
     * The phase to model is selected with setPhase().
     *
     * @tparam T    type of the medium parameters and traveltimes (float or double)
     * @tparam NODE type of the grid nodes (Node2Dc, Node2Dcsp, ...)
     * @tparam S    type of the coordinate struct (sxz<T>)
     */
    template <typename T, typename NODE, typename S>
    class CellVTI_PSV {
    public:
        /// Number of medium parameters of the cells: the two vertical velocities and Thomsen's epsilon and delta
        static constexpr size_t nParams = 4;
        /**
         * @brief Constructor
         *
         * The qP phase is selected by default.
         *
         * @param n number of cells of the grid
         */
        CellVTI_PSV(const size_t n) :
        sign(1.0),
        Vp0(std::vector<T>(n)),
        Vs0(std::vector<T>(n)),
        epsilon(std::vector<T>(n)),
        delta(std::vector<T>(n)),
        paramsSet(0) {
        }

        /**
         * @brief Assign the vertical P-wave velocity of the cells
         * @param s values of @f$V_{P0}@f$, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setVp0(const std::vector<T>& s) {
            if ( Vp0.size() != s.size() ) {
                throw std::length_error("Error: Vp0 vectors of incompatible size.");
            }
            Vp0 = s;
            paramsSet |= 1;
            buildTables();
        }

        /**
         * @brief Assign the vertical S-wave velocity of the cells
         * @param s values of @f$V_{S0}@f$, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setVs0(const std::vector<T>& s) {
            if ( Vs0.size() != s.size() ) {
                throw std::length_error("Error: Vs0 vectors of incompatible size.");
            }
            Vs0 = s;
            paramsSet |= 2;
            buildTables();
        }

        /**
         * @brief Assign Thomsen's parameter @f$\epsilon@f$ of the cells
         * @param s values of @f$\epsilon@f$, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setEpsilon(const std::vector<T>& s) {
            if ( epsilon.size() != s.size() ) {
                throw std::length_error("Error: epsilon vectors of incompatible size.");
            }
            epsilon = s;
            paramsSet |= 4;
            buildTables();
        }

        /**
         * @brief Assign Thomsen's parameter @f$\delta@f$ of the cells
         * @param s values of @f$\delta@f$, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setDelta(const std::vector<T>& s) {
            if ( delta.size() != s.size() ) {
                throw std::length_error("Error: delta vectors of incompatible size.");
            }
            delta = s;
            paramsSet |= 8;
            buildTables();
        }

        /**
         * @brief Select the phase to model
         * @param p 1 for the qP wave, any other value for the qSV wave
         */
        void setPhase(const int p) {
            if ( p==1 ) sign = 1.;  // P wave
            else sign = -1.;        // SV wave
            buildTables();
        }

        /// @brief Not applicable: always throws std::logic_error
        void setSlowness(const std::vector<T>& s) {
            throw std::logic_error("Error: slowness not defined for CellVTI_PSV.");
        }

        /// @brief Not applicable: always throws std::logic_error
        const T getSlowness(const size_t i) const {
            throw std::logic_error("Error: slowness not defined for CellVTI_PSV.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setXi(const std::vector<T>& s) {
            throw std::logic_error("Error: xi not defined for CellVTI_PSV.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setTiltAngle(const std::vector<T>& s) {
            throw std::logic_error("Error: TiltAngle not defined for CellVTI_PSV.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setGamma(const std::vector<T>& s) {
            throw std::logic_error("Error: gamma not defined for CellVTI_PSV.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setS2(const std::vector<T>& r) {
            throw std::logic_error("Error: s2 not defined for CellVTI_PSV.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setS4(const std::vector<T>& r) {
            throw std::logic_error("Error: s4 not defined for CellVTI_PSV.");
        }

        /**
         * @brief Traveltime between two points of a cell
         * @param source first point
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime along the segment joining the two points
         */
        T computeDt(const S& source, const S& node,
                    const size_t cellNo) const {
            return source.getDistance( node ) /
                   gv.velocity(node.x - source.x, node.z - source.z, cellNo);
        }

        /**
         * @brief Traveltime between a node and a point of a cell
         * @param source node from which the ray originates
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime along the segment joining the two points
         */
        T computeDt(const NODE& source, const S& node,
                    const size_t cellNo) const {
            return source.getDistance( node ) /
                   gv.velocity(node.x - source.getX(), node.z - source.getZ(), cellNo);
        }

        /**
         * @brief Traveltime between two nodes of a cell
         * @param source node from which the ray originates
         * @param node   second node
         * @param cellNo index of the cell holding the two nodes
         * @return the traveltime along the segment joining the two nodes
         */
        T computeDt(const NODE& source, const NODE& node,
                    const size_t cellNo) const {
            return source.getDistance( node ) /
                   gv.velocity(node.getX() - source.getX(), node.getZ() - source.getZ(), cellNo);
        }

        /**
         * @brief Not applicable: always throws std::logic_error
         *
         * A single length cannot describe the sensitivity of an anisotropic
         * cell; the siv2 overload must be used instead.  This condition was
         * formerly reported with assert(), which is compiled out whenever
         * NDEBUG is defined.
         *
         * @param[in]  source unused
         * @param[in]  node   unused
         * @param[out] cell   unused
         * @throws std::logic_error always
         */
        void computeDistance(const NODE& source, const S& node,
                             siv<T>& cell) const {
            throw std::logic_error("Error: computeDistance with siv not defined for CellVTI_PSV.");
        }
        /**
         * @brief Components of a ray segment, stored in a siv2 struct
         * @param[in]  source node from which the ray segment originates
         * @param[in]  node   end point of the ray segment
         * @param[out] cell   struct holding the horizontal component of the
         *                    segment in member `v`, and its vertical component
         *                    in member `v2`
         */
        void computeDistance(const NODE& source, const S& node,
                             siv2<T>& cell) const {
            cell.v  = std::abs(node.x - source.getX());
            cell.v2 = std::abs(node.z - source.getZ());
        }

        /**
         * @brief Sensitivity of a segment traveltime to the medium parameters
         *
         * Formed by VTI_PSV_GroupVel::sensitivity(), which applies the envelope
         * theorem on the branch of the slowness surface carrying the arrival.
         *
         * @param[in]  source node from which the ray segment originates
         * @param[in]  node   end point of the ray segment
         * @param[out] cell   derivatives w/r to Vp0, Vs0, epsilon and delta,
         *                    in members `v`, `v2`, `v3` and `v4`
         */
        void computeDistance(const NODE& source, const S& node,
                             siv4<T>& cell) const {
            T s[4];
            gv.sensitivity(node.x - source.getX(), node.z - source.getZ(),
                           cell.i, s);
            cell.v  = s[0];
            cell.v2 = s[1];
            cell.v3 = s[2];
            cell.v4 = s[3];
        }

    private:
        /**
         * @brief (Re)build the group-velocity tables
         *
         * Does nothing until the four medium parameters have all been assigned,
         * so that the setters may be called in any order.
         */
        void buildTables() {
            if ( paramsSet == 15 ) {
                gv.build(Vp0, Vs0, epsilon, delta, sign);
            }
        }

        T sign;                  ///< +1 for the qP phase, -1 for the qSV phase
        std::vector<T> Vp0;      ///< vertical P-wave velocity of the cells
        std::vector<T> Vs0;      ///< vertical S-wave velocity of the cells
        std::vector<T> epsilon;  ///< Thomsen's parameter epsilon of the cells
        std::vector<T> delta;    ///< Thomsen's parameter delta of the cells
        unsigned char paramsSet; ///< bitmask of the medium parameters assigned
        VTI_PSV_GroupVel<T> gv;  ///< group-velocity tables
    };


    /**
     * @brief Cells with TTI anisotropy, qP or qSV phase, in 2D (y dimension ignored)
     *
     * Tilted counterpart of CellVTI_PSV: the same transversely isotropic medium,
     * with its symmetry axis rotated by the tilt angle @f$\theta_t@f$ of the
     * cell, assigned with setTiltAngle().
     *
     * Tilting the medium rotates its slowness and group-velocity surfaces
     * rigidly, so the group velocity of a segment is that of the underlying VTI
     * medium taken in the frame of the symmetry axis.  The components
     * @f$(l_x, l_z)@f$ of the segment are rotated with the same convention as
     * CellTiltedElliptical,
     * @f[ t_1 = l_x \cos\theta_t + l_z \sin\theta_t, \quad
     *     t_2 = l_z \cos\theta_t - l_x \sin\theta_t @f]
     * @f$t_2@f$ being the component along the symmetry axis.  A segment making
     * an angle @f$\psi@f$ with the vertical therefore makes an angle
     * @f$\psi + \theta_t@f$ with the symmetry axis, which lies at
     * @f$-\theta_t@f$ from the vertical.
     *
     * The tables of VTI_PSV_GroupVel are reused unchanged, and changing the
     * tilt angle does not rebuild them: they are expressed in the frame of the
     * symmetry axis, which the tilt does not alter.
     *
     * With @f$\theta_t = 0@f$ the class is equivalent to CellVTI_PSV.
     *
     * @pre Every cell satisfies @f$V_{P0} > 0@f$ and @f$V_{S0} \neq V_{P0}@f$.
     *
     * The phase to model is selected with setPhase().
     *
     * @tparam T    type of the medium parameters and traveltimes (float or double)
     * @tparam NODE type of the grid nodes (Node2Dc, Node2Dcsp, ...)
     * @tparam S    type of the coordinate struct (sxz<T>)
     *
     * @sa CellVTI_PSV, VTI_PSV_GroupVel, CellTTI_SH
     */
    template <typename T, typename NODE, typename S>
    class CellTTI_PSV {
    public:
        /// Number of medium parameters of the cells: the two vertical velocities, Thomsen's epsilon and delta, and the tilt angle
        static constexpr size_t nParams = 5;
        /**
         * @brief Constructor
         *
         * The qP phase and a vertical symmetry axis are selected by default.
         *
         * @param n number of cells of the grid
         */
        CellTTI_PSV(const size_t n) :
        sign(1.0),
        Vp0(std::vector<T>(n)),
        Vs0(std::vector<T>(n)),
        epsilon(std::vector<T>(n)),
        delta(std::vector<T>(n)),
        tAngle(std::vector<T>(n, 0.0)),
        ca(std::vector<T>(n, 1.0)),
        sa(std::vector<T>(n, 0.0)),
        paramsSet(0) {
        }

        /**
         * @brief Assign the vertical P-wave velocity of the cells
         * @param s values of @f$V_{P0}@f$, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setVp0(const std::vector<T>& s) {
            if ( Vp0.size() != s.size() ) {
                throw std::length_error("Error: Vp0 vectors of incompatible size.");
            }
            Vp0 = s;
            paramsSet |= 1;
            buildTables();
        }

        /**
         * @brief Assign the vertical S-wave velocity of the cells
         * @param s values of @f$V_{S0}@f$, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setVs0(const std::vector<T>& s) {
            if ( Vs0.size() != s.size() ) {
                throw std::length_error("Error: Vs0 vectors of incompatible size.");
            }
            Vs0 = s;
            paramsSet |= 2;
            buildTables();
        }

        /**
         * @brief Assign Thomsen's parameter @f$\epsilon@f$ of the cells
         * @param s values of @f$\epsilon@f$, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setEpsilon(const std::vector<T>& s) {
            if ( epsilon.size() != s.size() ) {
                throw std::length_error("Error: epsilon vectors of incompatible size.");
            }
            epsilon = s;
            paramsSet |= 4;
            buildTables();
        }

        /**
         * @brief Assign Thomsen's parameter @f$\delta@f$ of the cells
         * @param s values of @f$\delta@f$, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setDelta(const std::vector<T>& s) {
            if ( delta.size() != s.size() ) {
                throw std::length_error("Error: delta vectors of incompatible size.");
            }
            delta = s;
            paramsSet |= 8;
            buildTables();
        }

        /**
         * @brief Assign the angle between the symmetry axis and the vertical
         * @param s values of @f$\theta_t@f$, in radians, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setTiltAngle(const std::vector<T>& s) {
            if ( tAngle.size() != s.size() ) {
                throw std::length_error("Error: angle vectors of incompatible size.");
            }
            for ( size_t n=0; n<tAngle.size(); ++n ) {
                tAngle[n] = s[n];
                ca[n] = std::cos(s[n]);
                sa[n] = std::sin(s[n]);
            }
            // the group-velocity tables live in the frame of the symmetry axis
            // and are unaffected by the tilt: no rebuild needed
        }

        /**
         * @brief Select the phase to model
         * @param p 1 for the qP wave, any other value for the qSV wave
         */
        void setPhase(const int p) {
            if ( p==1 ) sign = 1.;  // P wave
            else sign = -1.;        // SV wave
            buildTables();
        }

        /// @brief Not applicable: always throws std::logic_error
        void setSlowness(const std::vector<T>& s) {
            throw std::logic_error("Error: slowness not defined for CellTTI_PSV.");
        }

        /// @brief Not applicable: always throws std::logic_error
        const T getSlowness(const size_t i) const {
            throw std::logic_error("Error: slowness not defined for CellTTI_PSV.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setXi(const std::vector<T>& s) {
            throw std::logic_error("Error: xi not defined for CellTTI_PSV.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setGamma(const std::vector<T>& s) {
            throw std::logic_error("Error: gamma not defined for CellTTI_PSV.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setS2(const std::vector<T>& r) {
            throw std::logic_error("Error: s2 not defined for CellTTI_PSV.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setS4(const std::vector<T>& r) {
            throw std::logic_error("Error: s4 not defined for CellTTI_PSV.");
        }

        /**
         * @brief Traveltime between two points of a cell
         * @param source first point
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime along the segment joining the two points
         */
        T computeDt(const S& source, const S& node,
                    const size_t cellNo) const {
            const T lx = node.x - source.x;
            const T lz = node.z - source.z;
            return source.getDistance( node ) /
                   gv.velocity(lx*ca[cellNo] + lz*sa[cellNo],
                               lz*ca[cellNo] - lx*sa[cellNo], cellNo);
        }

        /**
         * @brief Traveltime between a node and a point of a cell
         * @param source node from which the ray originates
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime along the segment joining the two points
         */
        T computeDt(const NODE& source, const S& node,
                    const size_t cellNo) const {
            const T lx = node.x - source.getX();
            const T lz = node.z - source.getZ();
            return source.getDistance( node ) /
                   gv.velocity(lx*ca[cellNo] + lz*sa[cellNo],
                               lz*ca[cellNo] - lx*sa[cellNo], cellNo);
        }

        /**
         * @brief Traveltime between two nodes of a cell
         * @param source node from which the ray originates
         * @param node   second node
         * @param cellNo index of the cell holding the two nodes
         * @return the traveltime along the segment joining the two nodes
         */
        T computeDt(const NODE& source, const NODE& node,
                    const size_t cellNo) const {
            const T lx = node.getX() - source.getX();
            const T lz = node.getZ() - source.getZ();
            return source.getDistance( node ) /
                   gv.velocity(lx*ca[cellNo] + lz*sa[cellNo],
                               lz*ca[cellNo] - lx*sa[cellNo], cellNo);
        }

        /**
         * @brief Not applicable: always throws std::logic_error
         * @param[in]  source unused
         * @param[in]  node   unused
         * @param[out] cell   unused
         * @throws std::logic_error always
         */
        void computeDistance(const NODE& source, const S& node,
                             siv<T>& cell) const {
            throw std::logic_error("Error: computeDistance with siv not defined for CellTTI_PSV.");
        }

        /**
         * @brief Components of a ray segment, stored in a siv2 struct
         * @param[in]  source node from which the ray segment originates
         * @param[in]  node   end point of the ray segment
         * @param[out] cell   struct holding the horizontal component of the
         *                    segment in member `v`, and its vertical component
         *                    in member `v2`
         *
         * @note As for the other anisotropic classes, these two components do
         * not suffice to form the sensitivity of the traveltime to the medium
         * parameters of a tilted transversely isotropic cell.
         */
        void computeDistance(const NODE& source, const S& node,
                             siv2<T>& cell) const {
            cell.v  = std::abs(node.x - source.getX());
            cell.v2 = std::abs(node.z - source.getZ());
        }

        /**
         * @brief Sensitivity of a segment traveltime to the medium parameters
         *
         * The segment is expressed in the frame of the symmetry axis and passed
         * to VTI_PSV_GroupVel.  Tilting the medium shifts the angle the segment
         * makes with that axis by the tilt angle, so the derivative with respect
         * to the tilt is the derivative with respect to the ray angle.
         *
         * @param[in]  source node from which the ray segment originates
         * @param[in]  node   end point of the ray segment
         * @param[out] cell   derivatives w/r to Vp0, Vs0, epsilon, delta and
         *                    the tilt angle, in members `v` to `v5`
         */
        void computeDistance(const NODE& source, const S& node,
                             siv5<T>& cell) const {
            const T lx = node.x - source.getX();
            const T lz = node.z - source.getZ();
            const T t1 = lx*ca[cell.i] + lz*sa[cell.i];
            const T t2 = lz*ca[cell.i] - lx*sa[cell.i];
            T s[4];
            gv.sensitivity(t1, t2, cell.i, s);
            cell.v  = s[0];
            cell.v2 = s[1];
            cell.v3 = s[2];
            cell.v4 = s[3];
            cell.v5 = gv.dt_dpsi(t1, t2, cell.i);
        }

    private:
        /**
         * @brief (Re)build the group-velocity tables
         *
         * Does nothing until the four medium parameters have all been assigned,
         * so that the setters may be called in any order.
         */
        void buildTables() {
            if ( paramsSet == 15 ) {
                gv.build(Vp0, Vs0, epsilon, delta, sign);
            }
        }

        T sign;                  ///< +1 for the qP phase, -1 for the qSV phase
        std::vector<T> Vp0;      ///< vertical P-wave velocity of the cells
        std::vector<T> Vs0;      ///< vertical S-wave velocity of the cells
        std::vector<T> epsilon;  ///< Thomsen's parameter epsilon of the cells
        std::vector<T> delta;    ///< Thomsen's parameter delta of the cells
        std::vector<T> tAngle;   ///< angle of the symmetry axis, in radians
        std::vector<T> ca;       ///< cosine of tAngle
        std::vector<T> sa;       ///< sine of tAngle
        unsigned char paramsSet; ///< bitmask of the medium parameters assigned
        VTI_PSV_GroupVel<T> gv;  ///< group-velocity tables
    };





    /**
     * @brief Cells with VTI anisotropy, SH phase, in 2D (y dimension ignored)
     *
     * The medium is described by the vertical S-wave velocity @f$V_{S0}@f$ and
     * by Thomsen's (1986) parameter @f$\gamma@f$.  The phase velocity of the SH
     * wave, for a phase angle @f$\theta@f$ measured from the vertical
     * (symmetry) axis, is
     * @f[ v(\theta) = V_{S0} \sqrt{1 + 2 \gamma \sin^2\theta} @f]
     * This phase velocity is elliptical, with semi-axis @f$V_{S0}@f$ along the
     * vertical (symmetry) axis and @f$V_{S0}\sqrt{1 + 2\gamma}@f$ along the
     * horizontal one.  The traveltime along a ray segment of components
     * @f$(l_x, l_z)@f$ therefore has the closed form
     * @f[ dt = \frac{1}{V_{S0}}
     *         \sqrt{\frac{l_x^2}{1 + 2\gamma} + l_z^2} @f]
     * which is the group (energy) velocity result.  Dividing the segment length
     * by the phase velocity taken in the direction of the segment would
     * underestimate the traveltime: phase and group velocities coincide only
     * along the symmetry directions.
     *
     * @tparam T    type of the medium parameters and traveltimes (float or double)
     * @tparam NODE type of the grid nodes (Node2Dc, Node2Dcsp, ...)
     * @tparam S    type of the coordinate struct (sxz<T>)
     */
    template <typename T, typename NODE, typename S>
    class CellVTI_SH {
    public:
        /// Number of medium parameters of the cells: the vertical S-wave velocity and Thomsen's gamma
        static constexpr size_t nParams = 2;
        /**
         * @brief Constructor
         * @param n number of cells of the grid
         */
        CellVTI_SH(const size_t n) :
        Vs0(std::vector<T>(n)),
        gamma(std::vector<T>(n)) {
        }

        /**
         * @brief Assign the vertical S-wave velocity of the cells
         * @param s values of @f$V_{S0}@f$, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setVs0(const std::vector<T>& s) {
            if ( Vs0.size() != s.size() ) {
                throw std::length_error("Error: Vs0 vectors of incompatible size.");
            }
            Vs0 = s;
        }

        /**
         * @brief Assign Thomsen's parameter @f$\gamma@f$ of the cells
         * @param s values of @f$\gamma@f$, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setGamma(const std::vector<T>& s) {
            if ( gamma.size() != s.size() ) {
                throw std::length_error("Error: gamma vectors of incompatible size.");
            }
            gamma = s;
        }

        /// @brief Not applicable: always throws std::logic_error
        void setSlowness(const std::vector<T>& s) {
            throw std::logic_error("Error: slowness not defined for CellVTI_SH.");
        }

        /// @brief Not applicable: always throws std::logic_error
        const T getSlowness(const size_t i) const {
            throw std::logic_error("Error: slowness not defined for CellVTI_SH.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setXi(const std::vector<T>& s) {
            throw std::logic_error("Error: xi not defined for CellVTI_SH.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setTiltAngle(const std::vector<T>& s) {
            throw std::logic_error("Error: TiltAngle not defined for CellVTI_SH.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setVp0(const std::vector<T>& s) {
            throw std::logic_error("Error: Vp0 not defined for CellVTI_SH.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setDelta(const std::vector<T>& s) {
            throw std::logic_error("Error: delta not defined for CellVTI_SH.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setEpsilon(const std::vector<T>& s) {
            throw std::logic_error("Error: epsilon not defined for CellVTI_SH.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setS2(const std::vector<T>& r) {
            throw std::logic_error("Error: s2 not defined for CellVTI_SH.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setS4(const std::vector<T>& r) {
            throw std::logic_error("Error: s4 not defined for CellVTI_SH.");
        }

        /**
         * @brief Traveltime between two points of a cell
         * @param source first point
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime along the segment joining the two points
         */
        T computeDt(const S& source, const S& node,
                    const size_t cellNo) const {
            T lx = node.x - source.x;
            T lz = node.z - source.z;
            return std::sqrt( lx*lx / (1. + 2.*gamma[cellNo]) + lz*lz ) / Vs0[cellNo];
        }

        /**
         * @brief Traveltime between a node and a point of a cell
         * @param source node from which the ray originates
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime along the segment joining the two points
         */
        T computeDt(const NODE& source, const S& node,
                    const size_t cellNo) const {
            T lx = node.x - source.getX();
            T lz = node.z - source.getZ();
            return std::sqrt( lx*lx / (1. + 2.*gamma[cellNo]) + lz*lz ) / Vs0[cellNo];
        }

        /**
         * @brief Traveltime between two nodes of a cell
         * @param source node from which the ray originates
         * @param node   second node
         * @param cellNo index of the cell holding the two nodes
         * @return the traveltime along the segment joining the two nodes
         */
        T computeDt(const NODE& source, const NODE& node,
                    const size_t cellNo) const {
            T lx = node.getX() - source.getX();
            T lz = node.getZ() - source.getZ();
            return std::sqrt( lx*lx / (1. + 2.*gamma[cellNo]) + lz*lz ) / Vs0[cellNo];
        }

        /**
         * @brief Not applicable: always throws std::logic_error
         *
         * A single length cannot describe the sensitivity of an anisotropic
         * cell; the siv2 overload must be used instead.  This condition was
         * formerly reported with assert(), which is compiled out whenever
         * NDEBUG is defined.
         *
         * @param[in]  source unused
         * @param[in]  node   unused
         * @param[out] cell   unused
         * @throws std::logic_error always
         */
        void computeDistance(const NODE& source, const S& node,
                             siv<T>& cell) const {
            throw std::logic_error("Error: computeDistance with siv not defined for CellVTI_SH.");
        }
        /**
         * @brief Components of a ray segment, stored in a siv2 struct
         * @param[in]  source node from which the ray segment originates
         * @param[in]  node   end point of the ray segment
         * @param[out] cell   struct holding the horizontal component of the
         *                    segment in member `v`, and its vertical component
         *                    in member `v2`
         */
        void computeDistance(const NODE& source, const S& node,
                             siv2<T>& cell) const {
            const T lx = node.x - source.getX();
            const T lz = node.z - source.getZ();
            const T g  = 1. + 2.*gamma[cell.i];
            const T dt = std::sqrt( lx*lx/g + lz*lz )/Vs0[cell.i];
            cell.v  = -dt/Vs0[cell.i];
            cell.v2 = ( dt > 0. ) ?
                      -lx*lx/(Vs0[cell.i]*Vs0[cell.i]*dt*g*g) : T(0);
        }

    private:
        std::vector<T> Vs0;    ///< vertical S-wave velocity of the cells
        std::vector<T> gamma;  ///< Thomsen's parameter gamma of the cells
    };


    /**
     * @brief Cells with TTI anisotropy, SH phase, in 2D (y dimension ignored)
     *
     * Tilted counterpart of CellVTI_SH: the symmetry axis is rotated by the tilt
     * angle @f$\theta_t@f$ of the cell, assigned with setTiltAngle() and using
     * the same convention as CellTiltedElliptical.
     *
     * The SH phase velocity being elliptical, the traveltime keeps a closed
     * form once the segment is expressed in the frame of the symmetry axis,
     * @f[ t_1 = l_x \cos\theta_t + l_z \sin\theta_t, \quad
     *     t_2 = l_z \cos\theta_t - l_x \sin\theta_t @f]
     * @f[ dt = \frac{1}{V_{S0}}
     *         \sqrt{\frac{t_1^2}{1 + 2\gamma} + t_2^2} @f]
     * @f$t_2@f$ being the component along the symmetry axis, which lies at
     * @f$-\theta_t@f$ from the vertical.
     *
     * With @f$\theta_t = 0@f$ the class is equivalent to CellVTI_SH.
     *
     * @tparam T    type of the medium parameters and traveltimes (float or double)
     * @tparam NODE type of the grid nodes (Node2Dc, Node2Dcsp, ...)
     * @tparam S    type of the coordinate struct (sxz<T>)
     *
     * @sa CellVTI_SH, CellTTI_PSV
     */
    template <typename T, typename NODE, typename S>
    class CellTTI_SH {
    public:
        /// Number of medium parameters of the cells: the vertical S-wave velocity, Thomsen's gamma and the tilt angle
        static constexpr size_t nParams = 3;
        /**
         * @brief Constructor
         *
         * A vertical symmetry axis is selected by default.
         *
         * @param n number of cells of the grid
         */
        CellTTI_SH(const size_t n) :
        Vs0(std::vector<T>(n)),
        gamma(std::vector<T>(n)),
        tAngle(std::vector<T>(n, 0.0)),
        ca(std::vector<T>(n, 1.0)),
        sa(std::vector<T>(n, 0.0)) {
        }

        /**
         * @brief Assign the vertical S-wave velocity of the cells
         * @param s values of @f$V_{S0}@f$, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setVs0(const std::vector<T>& s) {
            if ( Vs0.size() != s.size() ) {
                throw std::length_error("Error: Vs0 vectors of incompatible size.");
            }
            Vs0 = s;
        }

        /**
         * @brief Assign Thomsen's parameter @f$\gamma@f$ of the cells
         * @param s values of @f$\gamma@f$, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setGamma(const std::vector<T>& s) {
            if ( gamma.size() != s.size() ) {
                throw std::length_error("Error: gamma vectors of incompatible size.");
            }
            gamma = s;
        }

        /**
         * @brief Assign the angle between the symmetry axis and the vertical
         * @param s values of @f$\theta_t@f$, in radians, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setTiltAngle(const std::vector<T>& s) {
            if ( tAngle.size() != s.size() ) {
                throw std::length_error("Error: angle vectors of incompatible size.");
            }
            for ( size_t n=0; n<tAngle.size(); ++n ) {
                tAngle[n] = s[n];
                ca[n] = std::cos(s[n]);
                sa[n] = std::sin(s[n]);
            }
        }

        /// @brief Not applicable: always throws std::logic_error
        void setSlowness(const std::vector<T>& s) {
            throw std::logic_error("Error: slowness not defined for CellTTI_SH.");
        }

        /// @brief Not applicable: always throws std::logic_error
        const T getSlowness(const size_t i) const {
            throw std::logic_error("Error: slowness not defined for CellTTI_SH.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setXi(const std::vector<T>& s) {
            throw std::logic_error("Error: xi not defined for CellTTI_SH.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setVp0(const std::vector<T>& s) {
            throw std::logic_error("Error: Vp0 not defined for CellTTI_SH.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setDelta(const std::vector<T>& s) {
            throw std::logic_error("Error: delta not defined for CellTTI_SH.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setEpsilon(const std::vector<T>& s) {
            throw std::logic_error("Error: epsilon not defined for CellTTI_SH.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setS2(const std::vector<T>& r) {
            throw std::logic_error("Error: s2 not defined for CellTTI_SH.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setS4(const std::vector<T>& r) {
            throw std::logic_error("Error: s4 not defined for CellTTI_SH.");
        }

        /**
         * @brief Traveltime between two points of a cell
         * @param source first point
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime along the segment joining the two points
         */
        T computeDt(const S& source, const S& node,
                    const size_t cellNo) const {
            return dt(node.x - source.x, node.z - source.z, cellNo);
        }

        /**
         * @brief Traveltime between a node and a point of a cell
         * @param source node from which the ray originates
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime along the segment joining the two points
         */
        T computeDt(const NODE& source, const S& node,
                    const size_t cellNo) const {
            return dt(node.x - source.getX(), node.z - source.getZ(), cellNo);
        }

        /**
         * @brief Traveltime between two nodes of a cell
         * @param source node from which the ray originates
         * @param node   second node
         * @param cellNo index of the cell holding the two nodes
         * @return the traveltime along the segment joining the two nodes
         */
        T computeDt(const NODE& source, const NODE& node,
                    const size_t cellNo) const {
            return dt(node.getX() - source.getX(),
                      node.getZ() - source.getZ(), cellNo);
        }

        /**
         * @brief Not applicable: always throws std::logic_error
         * @param[in]  source unused
         * @param[in]  node   unused
         * @param[out] cell   unused
         * @throws std::logic_error always
         */
        void computeDistance(const NODE& source, const S& node,
                             siv<T>& cell) const {
            throw std::logic_error("Error: computeDistance with siv not defined for CellTTI_SH.");
        }

        /**
         * @brief Components of a ray segment, stored in a siv2 struct
         * @param[in]  source node from which the ray segment originates
         * @param[in]  node   end point of the ray segment
         * @param[out] cell   struct holding the horizontal component of the
         *                    segment in member `v`, and its vertical component
         *                    in member `v2`
         *
         * @note As for the other anisotropic classes, these two components do
         * not suffice to form the sensitivity of the traveltime to the medium
         * parameters of a tilted transversely isotropic cell.
         */
        void computeDistance(const NODE& source, const S& node,
                             siv2<T>& cell) const {
            cell.v  = std::abs(node.x - source.getX());
            cell.v2 = std::abs(node.z - source.getZ());
        }

        /**
         * @brief Sensitivity of a segment traveltime to the medium parameters
         *
         * With @f$ g = 1+2\gamma @f$ and
         * @f$ dt = V_{S0}^{-1}\sqrt{t_1^2/g + t_2^2} @f$,
         * @f[ \frac{\partial dt}{\partial V_{S0}} = -\frac{dt}{V_{S0}}, \quad
         *     \frac{\partial dt}{\partial \gamma}
         *       = -\frac{t_1^2}{V_{S0}^2\,dt\,g^2}, \quad
         *     \frac{\partial dt}{\partial \theta}
         *       = \frac{t_1 t_2 (g^{-1}-1)}{V_{S0}^2\,dt} @f]
         *
         * @param[in]  source node from which the ray segment originates
         * @param[in]  node   end point of the ray segment
         * @param[out] cell   derivatives w/r to the vertical S-wave velocity,
         *                    Thomsen's gamma and the tilt angle, in members
         *                    `v`, `v2` and `v3`; member `v4` is left at zero
         */
        void computeDistance(const NODE& source, const S& node,
                             siv4<T>& cell) const {
            const T lx = node.x - source.getX();
            const T lz = node.z - source.getZ();
            const T t1 = lx*ca[cell.i] + lz*sa[cell.i];
            const T t2 = lz*ca[cell.i] - lx*sa[cell.i];
            const T g  = 1. + 2.*gamma[cell.i];
            const T vs = Vs0[cell.i];
            const T dt = std::sqrt( t1*t1/g + t2*t2 )/vs;
            cell.v  = -dt/vs;
            cell.v2 = ( dt > 0. ) ? -t1*t1/(vs*vs*dt*g*g) : T(0);
            cell.v3 = ( dt > 0. ) ? t1*t2*(1./g - 1.)/(vs*vs*dt) : T(0);
            cell.v4 = T(0);
        }

    private:
        /// @brief Traveltime along a segment of components (lx, lz)
        T dt(const T lx, const T lz, const size_t cellNo) const {
            const T t1 = lx*ca[cellNo] + lz*sa[cellNo];
            const T t2 = lz*ca[cellNo] - lx*sa[cellNo];
            return std::sqrt( t1*t1 / (1. + 2.*gamma[cellNo]) + t2*t2 ) / Vs0[cellNo];
        }

        std::vector<T> Vs0;    ///< vertical S-wave velocity of the cells
        std::vector<T> gamma;  ///< Thomsen's parameter gamma of the cells
        std::vector<T> tAngle; ///< angle of the symmetry axis, in radians
        std::vector<T> ca;     ///< cosine of tAngle
        std::vector<T> sa;     ///< sine of tAngle
    };




    /**
     * @brief Cells with weakly anelliptical anisotropy, in 2D (y dimension ignored)
     *
     * The energy (group) velocity is expanded in powers of @f$\sin^2\theta@f$,
     * @f$\theta@f$ being the angle measured from the vertical axis,
     * @f[ v(\theta) = v_0 \left[ 1 + \left( s_2 + s_4 \sin^2\theta \right)
     *     \sin^2\theta \right] @f]
     * where @f$v_0@f$ is the vertical velocity of the cell and @f$s_2@f$ and
     * @f$s_4@f$ are the second- and fourth-order anisotropy coefficients.  As
     * the expression yields the energy velocity directly, the traveltime is
     * simply the length of the ray segment divided by @f$v(\theta)@f$.
     *
     * @tparam T    type of the medium parameters and traveltimes (float or double)
     * @tparam NODE type of the grid nodes (Node2Dc, Node2Dcsp, ...)
     * @tparam S    type of the coordinate struct (sxz<T>)
     */
    template <typename T, typename NODE, typename S>
    class CellWeaklyAnelliptical {
    public:
        /// Number of medium parameters of the cells: the vertical slowness and the two anisotropy coefficients
        static constexpr size_t nParams = 3;
        /**
         * @brief Constructor
         * @param n number of cells of the grid
         */
        CellWeaklyAnelliptical(const size_t n) :
        v0(std::vector<T>(n)),
        s2(std::vector<T>(n)),
        s4(std::vector<T>(n)) {
        }

        /**
         * @brief Get the vertical slowness of a cell
         *
         * @pre the vertical velocity of cell @p i is non-zero; it is not
         * validated, and a zero velocity produces a division by zero.
         *
         * @param i index of the cell
         * @return the reciprocal of the vertical velocity of cell @p i
         * @throws std::out_of_range if @p i is not a valid cell index
         */
        const T getSlowness(const size_t i) const {
            return 1. / v0.at(i);
        }

        /**
         * @brief Assign the vertical slowness of the cells
         *
         * The values are inverted and stored as vertical velocities.
         *
         * @pre every value of @p s is non-zero; this is not validated, and a
         * zero slowness produces a division by zero.
         *
         * @param s slowness values, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setSlowness(const std::vector<T>& s) {
            if ( v0.size() != s.size() ) {
                throw std::length_error("Error: slowness vectors of incompatible size.");
            }
            for ( size_t n=0; n<v0.size(); ++n ) {
                v0[n] = 1. / s[n];
            }
        }

        /**
         * @brief Assign the second-order anisotropy coefficient of the cells
         * @param r values of @f$s_2@f$, one per cell
         * @throws std::length_error if @p r does not hold one value per cell
         */
        void setS2(const std::vector<T>& r) {
            if ( s2.size() != r.size() ) {
                throw std::length_error("Error: s2 vectors of incompatible size.");
            }
            s2 = r;
        }

        /**
         * @brief Assign the fourth-order anisotropy coefficient of the cells
         * @param r values of @f$s_4@f$, one per cell
         * @throws std::length_error if @p r does not hold one value per cell
         */
        void setS4(const std::vector<T>& r) {
            if ( s4.size() != r.size() ) {
                throw std::length_error("Error: s4 vectors of incompatible size.");
            }
            s4 = r;
        }

        /// @brief Not applicable: always throws std::logic_error
        void setXi(const std::vector<T>& s) {
            throw std::logic_error("Error: xi not defined for CellWeaklyAnelliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setChi(const std::vector<T>& s) {
            throw std::logic_error("Error: chi not defined for CellWeaklyAnelliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setPsi(const std::vector<T>& s) {
            throw std::logic_error("Error: psi not defined for CellWeaklyAnelliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setTiltAngle(const std::vector<T>& s) {
            throw std::logic_error("Error: TiltAngle not defined for CellWeaklyAnelliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setVp0(const std::vector<T>& s) {
            throw std::logic_error("Error: Vp0 not defined for CellWeaklyAnelliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setVs0(const std::vector<T>& s) {
            throw std::logic_error("Error: Vs0 not defined for CellWeaklyAnelliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setDelta(const std::vector<T>& s) {
            throw std::logic_error("Error: delta not defined for CellWeaklyAnelliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setEpsilon(const std::vector<T>& s) {
            throw std::logic_error("Error: epsilon not defined for CellWeaklyAnelliptical.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setGamma(const std::vector<T>& s) {
            throw std::logic_error("Error: gamma not defined for CellWeaklyAnelliptical.");
        }

        /**
         * @brief Traveltime between two points of a cell
         * @param source first point
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime along the segment joining the two points
         */
        T computeDt(const S& source, const S& node,
                    const size_t cellNo) const {
            T v = get_energy_vel(node.x - source.x, node.z - source.z, cellNo);
            return source.getDistance( node ) / v;
        }

        /**
         * @brief Traveltime between a node and a point of a cell
         * @param source node from which the ray originates
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime along the segment joining the two points
         */
        T computeDt(const NODE& source, const S& node,
                    const size_t cellNo) const {
            T v = get_energy_vel(node.x - source.getX(), node.z - source.getZ(), cellNo);
            return source.getDistance( node ) / v;
        }

        /**
         * @brief Traveltime between two nodes of a cell
         * @param source node from which the ray originates
         * @param node   second node
         * @param cellNo index of the cell holding the two nodes
         * @return the traveltime along the segment joining the two nodes
         */
        T computeDt(const NODE& source, const NODE& node,
                    const size_t cellNo) const {
            T v = get_energy_vel(node.getX() - source.getX(), node.getZ() - source.getZ(), cellNo);
            return source.getDistance( node ) / v;
        }

        /**
         * @brief Not applicable: always throws std::logic_error
         *
         * A single length cannot describe the sensitivity of an anisotropic
         * cell; the siv2 overload must be used instead.  This condition was
         * formerly reported with assert(), which is compiled out whenever
         * NDEBUG is defined.
         *
         * @param[in]  source unused
         * @param[in]  node   unused
         * @param[out] cell   unused
         * @throws std::logic_error always
         */
        void computeDistance(const NODE& source, const S& node,
                             siv<T>& cell) const {
            throw std::logic_error("Error: computeDistance with siv not defined for CellWeaklyAnelliptical.");
        }
        /**
         * @brief Components of a ray segment, stored in a siv2 struct
         * @param[in]  source node from which the ray segment originates
         * @param[in]  node   end point of the ray segment
         * @param[out] cell   struct holding the horizontal component of the
         *                    segment in member `v`, and its vertical component
         *                    in member `v2`
         */
        void computeDistance(const NODE& source, const S& node,
                             siv2<T>& cell) const {
            cell.v  = std::abs(node.x - source.getX());
            cell.v2 = std::abs(node.z - source.getZ());
        }

        /**
         * @brief Sensitivity of a segment traveltime to the medium parameters
         *
         * With @f$ u = \sin^2\theta @f$, @f$ P = 1 + (s_2 + s_4 u)u @f$ and
         * @f$ dt = \ell s / P @f$, @f$s@f$ being the vertical slowness,
         * @f[ \frac{\partial dt}{\partial s} = \frac{dt}{s}, \quad
         *     \frac{\partial dt}{\partial s_2} = -\frac{dt\,u}{P}, \quad
         *     \frac{\partial dt}{\partial s_4} = -\frac{dt\,u^2}{P} @f]
         *
         * @param[in]  source node from which the ray segment originates
         * @param[in]  node   end point of the ray segment
         * @param[out] cell   derivatives w/r to the vertical slowness and to
         *                    the two anisotropy coefficients, in members `v`,
         *                    `v2` and `v3`; member `v4` is left at zero
         */
        void computeDistance(const NODE& source, const S& node,
                             siv4<T>& cell) const {
            const T lx = node.x - source.getX();
            const T lz = node.z - source.getZ();
            const T st = std::sin( std::atan2(lx, lz) );
            const T u  = st*st;
            const T P  = 1. + (s2[cell.i] + s4[cell.i]*u)*u;
            const T dt = std::sqrt(lx*lx + lz*lz)/(v0[cell.i]*P);
            cell.v  = dt*v0[cell.i];          // d(dt)/d(slowness), s = 1/v0
            cell.v2 = -dt*u/P;
            cell.v3 = -dt*u*u/P;
            cell.v4 = T(0);
        }

    private:
        std::vector<T> v0;  ///< vertical velocity of the cells
        std::vector<T> s2;  ///< second-order anisotropy coefficient of the cells
        std::vector<T> s4;  ///< fourth-order anisotropy coefficient of the cells

        /**
         * @brief Energy velocity in a given direction
         * @param dx     horizontal component of the propagation direction
         * @param dz     vertical component of the propagation direction
         * @param cellNo index of the cell
         * @return the energy velocity of cell @p cellNo in direction (@p dx, @p dz)
         */
        T get_energy_vel(const T dx, const T dz, const size_t cellNo) const {
            // theta: angle w/r to vertical axis
            T sin_theta = sin( atan2(dx, dz) );
            sin_theta *= sin_theta;  // square
            return v0[cellNo] * (1. + (s2[cellNo] + s4[cellNo] * sin_theta) * sin_theta);
        }
    };


    /**
     * @brief Cells with elliptical anisotropy, in 3D
     *
     * The medium is described by the vertical slowness @f$s_z@f$ of the cell
     * and by the two anisotropy ratios @f$\chi = s_x / s_z@f$ and
     * @f$\psi = s_y / s_z@f$, so that the traveltime along a segment of
     * components @f$(l_x, l_y, l_z)@f$ is
     * @f[ dt = s_z \sqrt{\chi^2 l_x^2 + \psi^2 l_y^2 + l_z^2} @f]
     * The axes of the ellipsoid are aligned with the axes of the grid.
     *
     * @tparam T    type of the medium parameters and traveltimes (float or double)
     * @tparam NODE type of the grid nodes (Node3Dc, Node3Dcsp, ...)
     * @tparam S    type of the coordinate struct (sxyz<T>)
     */
    template <typename T, typename NODE, typename S>
    class CellElliptical3D {
    public:
        /**
         * @brief Constructor
         * @param n number of cells of the grid
         */
        CellElliptical3D(const size_t n) :
        slowness(std::vector<T>(n)),
        chi(std::vector<T>(n)),
        psi(std::vector<T>(n)) {
        }

        /**
         * @brief Assign the vertical slowness of the cells
         * @param s values of @f$s_z@f$, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setSlowness(const std::vector<T>& s) {
            if ( slowness.size() != s.size() ) {
                throw std::length_error("Error: slowness vectors of incompatible size.");
            }
            slowness = s;
        }

        /**
         * @brief Get the vertical slowness of a cell
         * @param i index of the cell
         * @return slowness of cell @p i
         * @throws std::out_of_range if @p i is not a valid cell index
         */
        const T getSlowness(const size_t i) const {
            return slowness.at(i);
        }

        /**
         * @brief Assign the anisotropy ratio along x of the cells
         * @param s values of @f$\chi = s_x / s_z@f$, one per cell; they are
         *          squared and stored as such
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setChi(const std::vector<T>& s) {
            if ( chi.size() != s.size() ) {
                throw std::length_error("Error: chi vectors of incompatible size.");
            }
            for ( size_t n=0; n<chi.size(); ++n ) {
                chi[n] = s[n]*s[n];
            }
        }

        /**
         * @brief Assign the anisotropy ratio along y of the cells
         * @param s values of @f$\psi = s_y / s_z@f$, one per cell; they are
         *          squared and stored as such
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setPsi(const std::vector<T>& s) {
            if ( psi.size() != s.size() ) {
                throw std::length_error("Error: psi vectors of incompatible size.");
            }
            for ( size_t n=0; n<psi.size(); ++n ) {
                psi[n] = s[n]*s[n];
            }
        }

        /// @brief Not applicable: always throws std::logic_error
        void setTiltAngle(const std::vector<T>& s) {
            throw std::logic_error("Error: TiltAngle not defined for CellElliptical3D.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setVp0(const std::vector<T>& s) {
            throw std::logic_error("Error: Vp0 not defined for CellElliptical3D.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setVs0(const std::vector<T>& s) {
            throw std::logic_error("Error: Vs0 not defined for CellElliptical3D.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setDelta(const std::vector<T>& s) {
            throw std::logic_error("Error: delta not defined for CellElliptical3D.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setEpsilon(const std::vector<T>& s) {
            throw std::logic_error("Error: epsilon not defined for CellElliptical3D.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setGamma(const std::vector<T>& s) {
            throw std::logic_error("Error: gamma not defined for CellElliptical3D.");
        }

        /**
         * @brief Traveltime between two points of a cell
         * @param source first point
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime along the segment joining the two points
         */
        T computeDt(const S& source, const S& node,
                    const size_t cellNo) const {
            T lx = node.x - source.x;
            T ly = node.y - source.y;
            T lz = node.z - source.z;
            return slowness[cellNo] * std::sqrt( chi[cellNo]*lx*lx + psi[cellNo]*ly*ly + lz*lz );
        }

        /**
         * @brief Traveltime between a node and a point of a cell
         * @param source node from which the ray originates
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime along the segment joining the two points
         */
        T computeDt(const NODE& source, const S& node,
                    const size_t cellNo) const {
            T lx = node.x - source.getX();
            T ly = node.y - source.getY();
            T lz = node.z - source.getZ();
            return slowness[cellNo] * std::sqrt( chi[cellNo]*lx*lx + psi[cellNo]*ly*ly + lz*lz );
        }

        /**
         * @brief Traveltime between two nodes of a cell
         * @param source node from which the ray originates
         * @param node   second node
         * @param cellNo index of the cell holding the two nodes
         * @return the traveltime along the segment joining the two nodes
         */
        T computeDt(const NODE& source, const NODE& node,
                    const size_t cellNo) const {
            T lx = node.getX() - source.getX();
            T ly = node.getY() - source.getY();
            T lz = node.getZ() - source.getZ();
            return slowness[cellNo] * std::sqrt( chi[cellNo]*lx*lx + psi[cellNo]*ly*ly + lz*lz );
        }

    private:
        std::vector<T> slowness;  ///< this vector contains sz
        std::vector<T> chi;       ///< anisotropy ratio, chi = sx / sz, *** squared ***
        std::vector<T> psi;       ///< anisotropy ratio, psi = sy / sz, *** squared ***
    };



    /**
     * @brief Cells with VTI anisotropy, P or SV phase, in 3D
     *
     * Three-dimensional counterpart of CellVTI_PSV: the medium being
     * transversely isotropic about the vertical axis, the phase velocity
     * depends only on the angle @f$\theta@f$ between the ray segment and the
     * vertical axis, which is here computed from the horizontal offset
     * @f$\sqrt{l_x^2 + l_y^2}@f$ and the vertical offset @f$l_z@f$.  See
     * CellVTI_PSV for the expression of @f$v(\theta)@f$.
     *
     * @pre Every cell satisfies @f$V_{P0} > 0@f$ and @f$V_{S0} \neq V_{P0}@f$.  The
     * values are not validated: @f$V_{P0} = 0@f$ or @f$V_{S0} = V_{P0}@f$
     * (which makes @f$f = 0@f$) produces a division by zero.
     *
     * The phase to model is selected with setPhase().
     *
     * @tparam T    type of the medium parameters and traveltimes (float or double)
     * @tparam NODE type of the grid nodes (Node3Dc, Node3Dcsp, ...)
     * @tparam S    type of the coordinate struct (sxyz<T>)
     */
    template <typename T, typename NODE, typename S>
    class CellVTI_PSV3D {
    public:
        /**
         * @brief Constructor
         *
         * The qP phase is selected by default.
         *
         * @param n number of cells of the grid
         */
        CellVTI_PSV3D(const size_t n) :
        sign(1.0),
        Vp0(std::vector<T>(n)),
        Vs0(std::vector<T>(n)),
        epsilon(std::vector<T>(n)),
        delta(std::vector<T>(n)),
        paramsSet(0) {
        }

        /**
         * @brief Assign the vertical P-wave velocity of the cells
         * @param s values of @f$V_{P0}@f$, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setVp0(const std::vector<T>& s) {
            if ( Vp0.size() != s.size() ) {
                throw std::length_error("Error: Vp0 vectors of incompatible size.");
            }
            Vp0 = s;
            paramsSet |= 1;
            buildTables();
        }

        /**
         * @brief Assign the vertical S-wave velocity of the cells
         * @param s values of @f$V_{S0}@f$, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setVs0(const std::vector<T>& s) {
            if ( Vs0.size() != s.size() ) {
                throw std::length_error("Error: Vs0 vectors of incompatible size.");
            }
            Vs0 = s;
            paramsSet |= 2;
            buildTables();
        }

        /**
         * @brief Assign Thomsen's parameter @f$\epsilon@f$ of the cells
         * @param s values of @f$\epsilon@f$, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setEpsilon(const std::vector<T>& s) {
            if ( epsilon.size() != s.size() ) {
                throw std::length_error("Error: epsilon vectors of incompatible size.");
            }
            epsilon = s;
            paramsSet |= 4;
            buildTables();
        }

        /**
         * @brief Assign Thomsen's parameter @f$\delta@f$ of the cells
         * @param s values of @f$\delta@f$, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setDelta(const std::vector<T>& s) {
            if ( delta.size() != s.size() ) {
                throw std::length_error("Error: delta vectors of incompatible size.");
            }
            delta = s;
            paramsSet |= 8;
            buildTables();
        }

        /**
         * @brief Select the phase to model
         * @param p 1 for the qP wave, any other value for the qSV wave
         */
        void setPhase(const int p) {
            if ( p==1 ) sign = 1.;  // P wave
            else sign = -1.;        // SV wave
            buildTables();
        }

        /// @brief Not applicable: always throws std::logic_error
        void setXi(const std::vector<T>& s) {
            throw std::logic_error("Error: xi not defined for CellVTI_PSV3D.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setTiltAngle(const std::vector<T>& s) {
            throw std::logic_error("Error: TiltAngle not defined for CellVTI_PSV3D.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setGamma(const std::vector<T>& s) {
            throw std::logic_error("Error: gamma not defined for CellVTI_PSV3D.");
        }

        /**
         * @brief Traveltime between a node and a point of a cell
         * @param source node from which the ray originates
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime along the segment joining the two points
         */
        T computeDt(const NODE& source, const S& node,
                    const size_t cellNo) const {
            T lx = node.x - source.getX();
            T ly = node.y - source.getY();
            return source.getDistance( node ) /
                   gv.velocity(std::sqrt(lx*lx + ly*ly), node.z - source.getZ(), cellNo);
        }

        /**
         * @brief Traveltime between two nodes of a cell
         * @param source node from which the ray originates
         * @param node   second node
         * @param cellNo index of the cell holding the two nodes
         * @return the traveltime along the segment joining the two nodes
         */
        T computeDt(const NODE& source, const NODE& node,
                    const size_t cellNo) const {
            T lx = node.getX() - source.getX();
            T ly = node.getY() - source.getY();
            return source.getDistance( node ) /
                   gv.velocity(std::sqrt(lx*lx + ly*ly), node.getZ() - source.getZ(), cellNo);
        }

    private:
        /**
         * @brief (Re)build the group-velocity tables
         *
         * Does nothing until the four medium parameters have all been assigned,
         * so that the setters may be called in any order.
         */
        void buildTables() {
            if ( paramsSet == 15 ) {
                gv.build(Vp0, Vs0, epsilon, delta, sign);
            }
        }

        T sign;                  ///< +1 for the qP phase, -1 for the qSV phase
        std::vector<T> Vp0;      ///< vertical P-wave velocity of the cells
        std::vector<T> Vs0;      ///< vertical S-wave velocity of the cells
        std::vector<T> epsilon;  ///< Thomsen's parameter epsilon of the cells
        std::vector<T> delta;    ///< Thomsen's parameter delta of the cells
        unsigned char paramsSet; ///< bitmask of the medium parameters assigned
        VTI_PSV_GroupVel<T> gv;  ///< group-velocity tables
    };



    /**
     * @brief Cells with VTI anisotropy, SH phase, in 3D
     *
     * Three-dimensional counterpart of CellVTI_SH.  The phase velocity for a
     * phase angle @f$\theta@f$ measured from the vertical (symmetry) axis is
     * @f[ v(\theta) = V_{S0} \sqrt{1 + 2 \gamma \sin^2\theta} @f]
     * and, being elliptical, yields the closed-form traveltime along a ray
     * segment of components @f$(l_x, l_y, l_z)@f$
     * @f[ dt = \frac{1}{V_{S0}}
     *         \sqrt{\frac{l_x^2 + l_y^2}{1 + 2\gamma} + l_z^2} @f]
     * which is the group (energy) velocity result.  See CellVTI_SH for why the
     * phase velocity must not be used directly.
     *
     * @tparam T    type of the medium parameters and traveltimes (float or double)
     * @tparam NODE type of the grid nodes (Node3Dc, Node3Dcsp, ...)
     * @tparam S    type of the coordinate struct (sxyz<T>)
     */
    template <typename T, typename NODE, typename S>
    class CellVTI_SH3D {
    public:
        /**
         * @brief Constructor
         * @param n number of cells of the grid
         */
        CellVTI_SH3D(const size_t n) :
        Vs0(std::vector<T>(n)),
        gamma(std::vector<T>(n)) {
        }

        /**
         * @brief Assign the vertical S-wave velocity of the cells
         * @param s values of @f$V_{S0}@f$, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setVs0(const std::vector<T>& s) {
            if ( Vs0.size() != s.size() ) {
                throw std::length_error("Error: Vs0 vectors of incompatible size.");
            }
            Vs0 = s;
        }

        /**
         * @brief Assign Thomsen's parameter @f$\gamma@f$ of the cells
         * @param s values of @f$\gamma@f$, one per cell
         * @throws std::length_error if @p s does not hold one value per cell
         */
        void setGamma(const std::vector<T>& s) {
            if ( gamma.size() != s.size() ) {
                throw std::length_error("Error: gamma vectors of incompatible size.");
            }
            gamma = s;
        }

        /// @brief Not applicable: always throws std::logic_error
        void setXi(const std::vector<T>& s) {
            throw std::logic_error("Error: xi not defined for CellVTI_SH3D.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setTiltAngle(const std::vector<T>& s) {
            throw std::logic_error("Error: TiltAngle not defined for CellVTI_SH3D.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setVp0(const std::vector<T>& s) {
            throw std::logic_error("Error: Vp0 not defined for CellVTI_SH3D.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setDelta(const std::vector<T>& s) {
            throw std::logic_error("Error: delta not defined for CellVTI_SH3D.");
        }

        /// @brief Not applicable: always throws std::logic_error
        void setEpsilon(const std::vector<T>& s) {
            throw std::logic_error("Error: epsilon not defined for CellVTI_SH3D.");
        }

        /**
         * @brief Traveltime between a node and a point of a cell
         * @param source node from which the ray originates
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime along the segment joining the two points
         */
        T computeDt(const NODE& source, const S& node,
                    const size_t cellNo) const {
            T lx = node.x - source.getX();
            T ly = node.y - source.getY();
            T lz = node.z - source.getZ();
            return std::sqrt( (lx*lx + ly*ly) / (1. + 2.*gamma[cellNo]) + lz*lz ) / Vs0[cellNo];
        }

        /**
         * @brief Traveltime between two nodes of a cell
         * @param source node from which the ray originates
         * @param node   second node
         * @param cellNo index of the cell holding the two nodes
         * @return the traveltime along the segment joining the two nodes
         */
        T computeDt(const NODE& source, const NODE& node,
                    const size_t cellNo) const {
            T lx = node.getX() - source.getX();
            T ly = node.getY() - source.getY();
            T lz = node.getZ() - source.getZ();
            return std::sqrt( (lx*lx + ly*ly) / (1. + 2.*gamma[cellNo]) + lz*lz ) / Vs0[cellNo];
        }

    private:
        std::vector<T> Vs0;    ///< vertical S-wave velocity of the cells
        std::vector<T> gamma;  ///< Thomsen's parameter gamma of the cells
    };

    /**
     * @brief Cells with weakly anelliptical anisotropy, in 3D
     *
     * Three-dimensional counterpart of CellWeaklyAnelliptical: the angle
     * @f$\theta@f$ from the vertical axis is computed from the horizontal
     * offset @f$\sqrt{l_x^2 + l_y^2}@f$ and the vertical offset @f$l_z@f$, and
     * the energy velocity is
     * @f[ v(\theta) = v_0 \left[ 1 + \left( s_2 + s_4 \sin^2\theta \right)
     *     \sin^2\theta \right] @f]
     *
     * @warning Unlike its 2D counterpart, this class takes vertical velocities
     * rather than slownesses (setV0()), and provides neither getSlowness() nor
     * the setters that throw for the parameters it does not support.
     *
     * @tparam T    type of the medium parameters and traveltimes (float or double)
     * @tparam NODE type of the grid nodes (Node3Dc, Node3Dcsp, ...)
     * @tparam S    type of the coordinate struct (sxyz<T>)
     */
    template <typename T, typename NODE, typename S>
    class CellWeaklyAnelliptical3D {
    public:
        /**
         * @brief Constructor
         * @param n number of cells of the grid
         */
        CellWeaklyAnelliptical3D(const size_t n) :
        v0(std::vector<T>(n)),
        s2(std::vector<T>(n)),
        s4(std::vector<T>(n)) {
        }

        /**
         * @brief Assign the vertical velocity of the cells
         * @param v values of @f$v_0@f$, one per cell
         * @throws std::length_error if @p v does not hold one value per cell
         */
        void setV0(const std::vector<T>& v) {
            if ( v0.size() != v.size() ) {
                throw std::length_error("Error: velocity vectors of incompatible size.");
            }
            v0 = v;
        }

        /**
         * @brief Assign the second-order anisotropy coefficient of the cells
         * @param r values of @f$s_2@f$, one per cell
         * @throws std::length_error if @p r does not hold one value per cell
         */
        void setS2(const std::vector<T>& r) {
            if ( s2.size() != r.size() ) {
                throw std::length_error("Error: s2 vectors of incompatible size.");
            }
            s2 = r;
        }

        /**
         * @brief Assign the fourth-order anisotropy coefficient of the cells
         * @param r values of @f$s_4@f$, one per cell
         * @throws std::length_error if @p r does not hold one value per cell
         */
        void setS4(const std::vector<T>& r) {
            if ( s4.size() != r.size() ) {
                throw std::length_error("Error: s4 vectors of incompatible size.");
            }
            s4 = r;
        }

        /**
         * @brief Traveltime between a node and a point of a cell
         * @param source node from which the ray originates
         * @param node   second point
         * @param cellNo index of the cell holding the two points
         * @return the traveltime along the segment joining the two points
         */
        T computeDt(const NODE& source, const S& node,
                    const size_t cellNo) const {
            // theta: angle w/r to vertical axis
            T lx = node.x - source.getX();
            T ly = node.y - source.getY();
            lx = std::sqrt( lx*lx + ly*ly ); // horizontal distance
            T v = get_energy_vel(lx, node.z - source.getZ(), cellNo);
            return source.getDistance( node ) / v;
        }

        /**
         * @brief Traveltime between two nodes of a cell
         * @param source node from which the ray originates
         * @param node   second node
         * @param cellNo index of the cell holding the two nodes
         * @return the traveltime along the segment joining the two nodes
         */
        T computeDt(const NODE& source, const NODE& node,
                    const size_t cellNo) const {
            // theta: angle w/r to vertical axis
            T lx = node.getX() - source.getX();
            T ly = node.getY() - source.getY();
            lx = std::sqrt( lx*lx + ly*ly ); // horizontal distance
            T v = get_energy_vel(lx, node.getZ() - source.getZ(), cellNo);
            return source.getDistance( node ) / v;
        }

    private:
        std::vector<T> v0;  ///< vertical velocity of the cells
        std::vector<T> s2;  ///< second-order anisotropy coefficient of the cells
        std::vector<T> s4;  ///< fourth-order anisotropy coefficient of the cells

        /**
         * @brief Energy velocity in a given direction
         * @param dx     horizontal component of the propagation direction
         * @param dz     vertical component of the propagation direction
         * @param cellNo index of the cell
         * @return the energy velocity of cell @p cellNo in direction (@p dx, @p dz)
         */
        T get_energy_vel(const T dx, const T dz, const size_t cellNo) const {
            // theta: angle w/r to vertical axis
            T sin_theta = sin( atan2(dx, dz) );
            sin_theta *= sin_theta;  // square
            return v0[cellNo] * (1. + (s2[cellNo] + s4[cellNo] * sin_theta) * sin_theta);
        }

    };

}

#endif /* Cell_h */

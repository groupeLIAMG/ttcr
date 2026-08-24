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
            if ( tables.empty() ) {
                throw std::logic_error("Error: medium parameters of VTI_PSV cells not set.");
            }
            const std::vector<T>& tab = tables[ index[cellNo] ];
            // VTI symmetry: the table spans [0, pi/2] only
            const T psi = std::atan2( std::abs(lh), std::abs(lz) );
            const T x = psi / halfPi() * (nSamples-1);
            size_t i = static_cast<size_t>(x);
            if ( i >= nSamples-1 ) {
                return tab[nSamples-1];
            }
            const T w = x - static_cast<T>(i);
            return (1.-w)*tab[i] + w*tab[i+1];
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

    private:
        /// Number of samples of a table, spanning [0, pi/2]
        static const size_t nSamples = 901;
        /// Oversampling of the phase angle when building a table
        static const size_t oversampling = 16;

        std::vector<std::vector<T>> tables;  ///< one table per distinct medium
        std::vector<size_t> index;           ///< table used by each cell

        static T halfPi() { return static_cast<T>(1.57079632679489661923); }
        static T pi()     { return static_cast<T>(3.14159265358979323846); }

        /// @brief Tabulate the group velocity over [0, pi/2] for one medium
        static std::vector<T> tabulate(const T Vp0, const T Vs0, const T eps,
                                       const T dlt, const T sign) {
            std::vector<T> tab(nSamples, T(0));
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
                if ( vg > tab[j] ) tab[j] = vg;   // keep the first arrival
            }
            fillGaps(tab);
            return tab;
        }

        /// @brief Linearly interpolate the bins that no sample reached
        static void fillGaps(std::vector<T>& tab) {
            size_t lo = 0;
            while ( lo < tab.size() && tab[lo] == T(0) ) ++lo;
            if ( lo == tab.size() ) {
                throw std::runtime_error("Error: could not tabulate the group "
                                         "velocity of a VTI_PSV cell.");
            }
            for ( size_t i=0; i<lo; ++i ) tab[i] = tab[lo];
            size_t hi = tab.size()-1;
            while ( tab[hi] == T(0) ) --hi;
            for ( size_t i=hi+1; i<tab.size(); ++i ) tab[i] = tab[hi];
            size_t i = lo;
            while ( i < hi ) {
                if ( tab[i+1] != T(0) ) { ++i; continue; }
                size_t j = i+1;
                while ( tab[j] == T(0) ) ++j;
                for ( size_t k=i+1; k<j; ++k ) {
                    const T w = static_cast<T>(k-i)/static_cast<T>(j-i);
                    tab[k] = (1.-w)*tab[i] + w*tab[j];
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
            cell.v  = std::abs(node.x - source.getX());
            cell.v2 = std::abs(node.z - source.getZ());
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

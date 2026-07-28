//
//  utils.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-02-15.
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

/**
 * @file utils.h
 * @brief Free-function toolbox: computational geometry, reflector construction
 *        and VTK output helpers.
 *
 * A grab-bag of small templates shared by the grid and mesh classes. Three
 * loose groups:
 *
 * - **Geometric predicates and measures** — @ref ttcr::areCoplanar,
 *   @ref ttcr::areCollinear, @ref ttcr::testInTriangle,
 *   @ref ttcr::testInTriangleBoundingBox, @ref ttcr::barycentric,
 *   @ref ttcr::triangleArea2D, @ref ttcr::distSqPointToSegment,
 *   @ref ttcr::distPointToPlane, @ref ttcr::distPointToLine. These underpin
 *   point location on unstructured meshes and the raypath/interface
 *   intersection tests.
 * - **Model and result plumbing** — @ref ttcr::buildReflectors turns the named
 *   surfaces of a Gmsh mesh into receiver arrays; @ref ttcr::saveRayPaths
 *   writes traced rays to VTK.
 * - **Odds and ends** — @ref ttcr::factorial, @ref ttcr::pseudoInverse,
 *   @ref ttcr::to_string.
 *
 * @section utils_tol Tolerances
 * The predicates compare against @ref ttcr::small2 (that is, @ref ttcr::small
 * squared, @f$10^{-8}@f$) rather than testing for exact equality, so they are
 * tolerant of round-off but **not scale-invariant**: the tolerance is absolute,
 * in model units. On a model whose coordinates are in metres this is a tenth of
 * a millimetre; on one in kilometres it is a tenth of a micron.
 *
 * @section utils_overloads Overload pattern
 * Most predicates come in two or four flavours, differing only in how the
 * points arrive: as `NODE*` (reading coordinates through `getX()`/`getY()`/
 * `getZ()`) or as plain @ref ttcr::sxyz / @ref ttcr::sxz values, and in 3-D or
 * 2-D (x-z). The bodies are otherwise identical.
 *
 * @note Several functions here are compiled to no-ops unless @c VTK is defined;
 *       each says so.
 *
 * @sa ttcr_t.h, Interpolator.h, MSHReader.h, Rcv.h
 */

#ifndef ttcr_utils_h
#define ttcr_utils_h

#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>

#ifdef VTK
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyLine.h"
#include "vtkSmartPointer.h"
#include "vtkXMLPolyDataWriter.h"
#endif

#include <Eigen/Dense>

#include "MSHReader.h"
#include "Rcv.h"

namespace ttcr {

    /**
     * @brief Compute the factorial of a number.
     *
     * Recursive and @c constexpr, so a call with a literal argument is folded
     * at compile time.
     *
     * @tparam T integer type of the argument and result.
     * @param n number to compute
     * @returns value of factorial
     * @warning No overflow check. @f$13!@f$ already exceeds 32 bits and
     *          @f$21!@f$ exceeds 64, and the result silently wraps.
     */
    template<typename T>
    constexpr T factorial(T n)
    {
        return n <= 1 ? 1 : (n * factorial(n-1));
    }

    /**
     * @brief Test if four points are coplanar.
     *
     * Points are sxyz objects. Coplanarity is judged by the scalar triple
     * product, @f$|(x_3-x_1)\cdot[(x_2-x_1)\times(x_4-x_3)]|@f$, which is
     * (twice) the volume of the tetrahedron they span: zero exactly when the
     * four points share a plane.
     *
     * @tparam T underlying type of sxyz objects
     * @param x1 first point
     * @param x2 second point
     * @param x3 third point
     * @param x4 fourth point
     * @returns result of test
     * @note The triple product carries the cube of the coordinate scale while
     *       @ref small2 is a fixed absolute tolerance, so this test grows
     *       stricter as the model shrinks and looser as it grows. See
     *       @ref utils_tol.
     */
    template<typename T>
    bool areCoplanar(const sxyz<T> &x1, const sxyz<T> &x2,
                     const sxyz<T> &x3, const sxyz<T> &x4) {
        return (std::abs( dot( x3-x1, cross(x2-x1, x4-x3) ) )<small2);
    }

    /**
     * @brief Test if four points are coplanar.
     *
     * Overload taking the last three points as Node objects; otherwise
     * identical to @ref areCoplanar(const sxyz<T>&, const sxyz<T>&, const sxyz<T>&, const sxyz<T>&).
     *
     * @tparam T underlying type of sxyz object
     * @tparam NODE type of Node objects
     * @param x1 first point
     * @param x2 second point
     * @param x3 third point
     * @param x4 fourth point
     * @returns result of test
     */
    template<typename T, typename NODE>
    bool areCoplanar(const sxyz<T> &x1, const NODE &x2,
                     const NODE &x3, const NODE &x4) {
        return (std::abs( dot( x3-x1, cross(x2-x1, x4-x3) ) )<small2);
    }

    /**
     * @brief Test if three points are collinear.
     *
     * Points are sxyz and Node objects. Collinearity is judged from the cross
     * product @f$(pt-n_0)\times(pt-n_1)@f$, whose magnitude is twice the area
     * of the triangle they span and so vanishes exactly when they share a line.
     *
     * @tparam T underlying type of sxyz object
     * @tparam NODE type of Node objects
     * @param pt first point
     * @param n0 second point
     * @param n1 third point
     * @returns result of test
     * @note Compares an unsquared norm against @ref small2, which is @ref small
     *       *squared* — so the effective tolerance here is @f$10^{-8}@f$, not
     *       @f$10^{-4}@f$, making this markedly stricter than the name suggests.
     * @sa http://mathworld.wolfram.com/Collinear.html
     */
    template<typename T, typename NODE>
    bool areCollinear(const sxyz<T> &pt, const NODE &n0, const NODE &n1) {

        // http://mathworld.wolfram.com/Collinear.html
        //
        sxyz<T> v = cross(pt-n0, pt-n1);
        return norm(v)<small2;
    }

    /**
     * @brief Test if three points are collinear, in the x-z plane.
     *
     * The 2-D counterpart of
     * @ref areCollinear(const sxyz<T>&, const NODE&, const NODE&). In 2-D the
     * cross product is the scalar @f$(pt-n_0)\times(pt-n_1)@f$, so its absolute
     * value is compared directly.
     *
     * @tparam T underlying type of sxz object
     * @tparam NODE type of Node objects
     * @param pt first point
     * @param n0 second point
     * @param n1 third point
     * @returns result of test
     * @note As in the 3-D overload the comparison is against @ref small2, i.e.
     *       @f$10^{-8}@f$.
     * @warning Calls @c abs unqualified. With @c \<cstdlib\> in scope ahead of
     *          the floating-point overloads this can select the integer
     *          @c abs(int) and truncate the argument to zero, making every
     *          triple look collinear.
     */
    template<typename T, typename NODE>
    bool areCollinear(const sxz<T> &pt, const NODE &n0, const NODE &n1) {

        // http://mathworld.wolfram.com/Collinear.html
        //
        T v = cross(pt-n0, pt-n1);
        return abs(v)<small2;
    }

    /**
     * @brief Compute twice the signed area of a 2-D triangle.
     *
     * Evaluates
     * @f[ 2A = x_1(y_2-y_3) + x_2(y_3-y_1) + x_3(y_1-y_2) @f]
     *
     * @tparam T type of coordinate values
     * @param x1 x coordinate of first point
     * @param y1 y coordinate of first point
     * @param x2 x coordinate of second point
     * @param y2 y coordinate of second point
     * @param x3 x coordinate of third point
     * @param y3 y coordinate of third point
     * @returns value of area
     *
     * @warning Despite the name this is **twice** the area, and it is
     *          **signed** — negative for a clockwise triple. Neither is an
     *          oversight: @ref barycentric divides these values by an
     *          equally-unnormalised denominator, so the factor of two cancels,
     *          and the sign is what lets it detect points outside the triangle.
     *          Callers wanting a true area must halve it and take the absolute
     *          value.
     */
    template<typename T>
    T triangleArea2D(T x1, T y1, T x2, T y2, T x3, T y3) {
        return (x1-x2)*(y2-y3) - (x2-x3)*(y1-y2);
    }

    /**
     * @brief Compute the barycentric coordinates of a point given a triangle.
     *
     * The triangle lives in 3-D, so the areas are computed after projecting
     * onto whichever coordinate plane the triangle is *least* edge-on to —
     * chosen from the largest component of the unnormalised normal. Projecting
     * onto a plane the triangle is nearly perpendicular to would collapse it
     * and destroy precision; picking the largest component maximises the
     * projected area and keeps the ratios well conditioned.
     *
     * The returned weights satisfy @f$u+v+w=1@f$ and reproduce the point as
     * @f$p = ua + vb + wc@f$. All three are non-negative exactly when @p p lies
     * inside the triangle, which is how @ref testInTriangle uses them.
     *
     * @tparam T underlying type of sxyz object
     * @tparam NODE type of Node objects
     * @param[in] a 1st node of triangle
     * @param[in] b 2nd node of triangle
     * @param[in] c 3rd node of triangle
     * @param[in] p point for which barycentric coordinate are computed
     * @param[out] u 1st barycentric coordinate
     * @param[out] v 2nd barycentric coordinate
     * @param[out] w 3rd barycentric coordinate
     *
     * @pre The three nodes are not collinear. A degenerate triangle has a zero
     *      normal, so the reciprocal is not finite and the weights come back as
     *      infinities or NaNs — the division is not guarded.
     * @note @p p is not required to lie in the plane of the triangle; it is
     *       effectively projected onto it, and no residual is reported.
     */
    template<typename T, typename NODE>
    void barycentric(const NODE *a,
                     const NODE *b,
                     const NODE *c,
                     const sxyz<T> &p,
                     T &u, T &v, T &w) {

        sxyz<T> ab = {b->getX()-a->getX(), b->getY()-a->getY(), b->getZ()-a->getZ()};
        sxyz<T> ac = {c->getX()-a->getX(), c->getY()-a->getY(), c->getZ()-a->getZ()};

        // Unnormalized triangle normal
        sxyz<T> m = cross(ab, ac);

        // Nominators and one-over-denominator for u and v ratios
        T nu, nv, ood;

        // Absolute components for determining projection plane
        T x = std::abs(m.x), y = std::abs(m.y), z = std::abs(m.z);

        // Compute areas in plane of largest projection
        if (x >= y && x >= z) {
            // x is largest, project to the yz plane
            nu = triangleArea2D(p.y, p.z, b->getY(), b->getZ(), c->getY(), c->getZ()); // Area of PBC in yz plane
            nv = triangleArea2D(p.y, p.z, c->getY(), c->getZ(), a->getY(), a->getZ()); // Area of PCA in yz plane
            ood = 1.0 / m.x; // 1/(2*area of ABC in yz plane)
        } else if (y >= x && y >= z) {
            // y is largest, project to the xz plane
            nu = triangleArea2D(p.x, p.z, b->getX(), b->getZ(), c->getX(), c->getZ());
            nv = triangleArea2D(p.x, p.z, c->getX(), c->getZ(), a->getX(), a->getZ());
            ood = 1.0 / -m.y;
        } else {
            // z is largest, project to the xy plane
            nu = triangleArea2D(p.x, p.y, b->getX(), b->getY(), c->getX(), c->getY());
            nv = triangleArea2D(p.x, p.y, c->getX(), c->getY(), a->getX(), a->getY());
            ood = 1.0 / m.z;
        }
        u = nu * ood;
        v = nv * ood;
        w = 1.0 - u - v;
    }

    /**
     * @brief Compute the barycentric coordinates of a point given a triangle.
     *
     * Overload taking the triangle as plain @ref sxyz values rather than Node
     * pointers; the algorithm is identical to
     * @ref barycentric(const NODE*, const NODE*, const NODE*, const sxyz<T>&, T&, T&, T&),
     * whose documentation describes the projection scheme and the degenerate
     * case.
     *
     * @tparam T underlying type of sxyz object
     * @param[in] a 1st node of triangle
     * @param[in] b 2nd node of triangle
     * @param[in] c 3rd node of triangle
     * @param[in] p point for which barycentric coordinate are computed
     * @param[out] u 1st barycentric coordinate
     * @param[out] v 2nd barycentric coordinate
     * @param[out] w 3rd barycentric coordinate
     *
     * @pre The three vertices are not collinear.
     */
    template<typename T>
    void barycentric(const sxyz<T> &a,
                     const sxyz<T> &b,
                     const sxyz<T> &c,
                     const sxyz<T> &p,
                     T &u, T &v, T &w) {

        sxyz<T> ab = b - a;
        sxyz<T> ac = c - a;

        // Unnormalized triangle normal
        sxyz<T> m = cross(ab, ac);

        // Nominators and one-over-denominator for u and v ratios
        T nu, nv, ood;

        // Absolute components for determining projection plane
        T x = std::abs(m.x), y = std::abs(m.y), z = std::abs(m.z);

        // Compute areas in plane of largest projection
        if (x >= y && x >= z) {
            // x is largest, project to the yz plane
            nu = triangleArea2D(p.y, p.z, b.y, b.z, c.y, c.z); // Area of PBC in yz plane
            nv = triangleArea2D(p.y, p.z, c.y, c.z, a.y, a.z); // Area of PCA in yz plane
            ood = 1.0 / m.x; // 1/(2*area of ABC in yz plane)
        } else if (y >= x && y >= z) {
            // y is largest, project to the xz plane
            nu = triangleArea2D(p.x, p.z, b.x, b.z, c.x, c.z);
            nv = triangleArea2D(p.x, p.z, c.x, c.z, a.x, a.z);
            ood = 1.0 / -m.y;
        } else {
            // z is largest, project to the xy plane
            nu = triangleArea2D(p.x, p.y, b.x, b.y, c.x, c.y);
            nv = triangleArea2D(p.x, p.y, c.x, c.y, a.x, a.y);
            ood = 1.0 / m.z;
        }
        u = nu * ood;
        v = nv * ood;
        w = 1.0 - u - v;
    }

    /**
     * @brief Test if a point is within the axis-aligned box bounding a triangle.
     *
     * A cheap conservative reject, used to skip the much costlier
     * @ref barycentric call in @ref testInTriangle. The box is grown by
     * @ref small2 on every side so points marginally outside still pass.
     *
     * @tparam T underlying type of sxyz object
     * @tparam NODE type of Node objects
     * @param vertexA 1st point defining triangle
     * @param vertexB 2nd point defining triangle
     * @param vertexC 3rd point defining triangle
     * @param E point to test
     * @returns test value
     *
     * @note One-sided: false is conclusive (the point cannot be in the
     *       triangle), true is not (the box is strictly larger than the
     *       triangle). Only ever use it as a pre-filter.
     */
    template<typename T, typename NODE>
    bool testInTriangleBoundingBox(const NODE *vertexA,
                                   const NODE *vertexB,
                                   const NODE *vertexC,
                                   const sxyz<T> &E) {
        T xMin = vertexA->getX() < vertexB->getX() ? vertexA->getX() : vertexB->getX();
        xMin = xMin < vertexC->getX() ? xMin : vertexC->getX();
        T xMax = vertexA->getX() > vertexB->getX() ? vertexA->getX() : vertexB->getX();
        xMax = xMax > vertexC->getX() ? xMax : vertexC->getX();

        T yMin = vertexA->getY() < vertexB->getY() ? vertexA->getY() : vertexB->getY();
        yMin = yMin < vertexC->getY() ? yMin : vertexC->getY();
        T yMax = vertexA->getY() > vertexB->getY() ? vertexA->getY() : vertexB->getY();
        yMax = yMax > vertexC->getY() ? yMax : vertexC->getY();

        T zMin = vertexA->getZ() < vertexB->getZ() ? vertexA->getZ() : vertexB->getZ();
        zMin = zMin < vertexC->getZ() ? zMin : vertexC->getZ();
        T zMax = vertexA->getZ() > vertexB->getZ() ? vertexA->getZ() : vertexB->getZ();
        zMax = zMax > vertexC->getZ() ? zMax : vertexC->getZ();

        if (E.x < xMin-small2 || xMax+small2 < E.x ||
            E.y < yMin-small2 || yMax+small2 < E.y ||
            E.z < zMin-small2 || zMax+small2 < E.z)
            return false;
        else
            return true;
    }

    /**
     * @brief Test if a point is within the box bounding a triangle, in the x-z plane.
     *
     * As
     * @ref testInTriangleBoundingBox(const NODE*, const NODE*, const NODE*, const sxyz<T>&),
     * but only x and z are bounded — the triangle's y extent is ignored.
     *
     * @tparam T underlying type of sxyz object
     * @tparam NODE type of Node objects
     * @param vertexA 1st point defining triangle
     * @param vertexB 2nd point defining triangle
     * @param vertexC 3rd point defining triangle
     * @param E point to test
     * @returns test value
     *
     * @note Conservative in the same one-sided sense as the 3-D overload.
     */
    template<typename T, typename NODE>
    bool testInTriangleBoundingBox(const NODE *vertexA,
                                   const NODE *vertexB,
                                   const NODE *vertexC,
                                   const sxz<T> &E) {
        T xMin = vertexA->getX() < vertexB->getX() ? vertexA->getX() : vertexB->getX();
        xMin = xMin < vertexC->getX() ? xMin : vertexC->getX();
        T xMax = vertexA->getX() > vertexB->getX() ? vertexA->getX() : vertexB->getX();
        xMax = xMax > vertexC->getX() ? xMax : vertexC->getX();

        T zMin = vertexA->getZ() < vertexB->getZ() ? vertexA->getZ() : vertexB->getZ();
        zMin = zMin < vertexC->getZ() ? zMin : vertexC->getZ();
        T zMax = vertexA->getZ() > vertexB->getZ() ? vertexA->getZ() : vertexB->getZ();
        zMax = zMax > vertexC->getZ() ? zMax : vertexC->getZ();

        if (E.x < xMin-small2 || xMax+small2 < E.x ||
            E.z < zMin-small2 || zMax+small2 < E.z)
            return false;
        else
            return true;
    }

    /**
     * @brief Test if a point is within the axis-aligned box bounding a triangle.
     *
     * Overload taking plain @ref sxyz values instead of Node pointers.
     *
     * @tparam T underlying type of sxyz objects
     * @param vertexA 1st point defining triangle
     * @param vertexB 2nd point defining triangle
     * @param vertexC 3rd point defining triangle
     * @param E point to test
     * @returns test value
     *
     * @note Conservative pre-filter: false is conclusive, true is not.
     */
    template<typename T>
    bool testInTriangleBoundingBox(const sxyz<T> &vertexA,
                                   const sxyz<T> &vertexB,
                                   const sxyz<T> &vertexC,
                                   const sxyz<T> &E) {
        T xMin = vertexA.x < vertexB.x ? vertexA.x : vertexB.x;
        xMin = xMin < vertexC.x ? xMin : vertexC.x;
        T xMax = vertexA.x > vertexB.x ? vertexA.x : vertexB.x;
        xMax = xMax > vertexC.x ? xMax : vertexC.x;

        T yMin = vertexA.y < vertexB.y ? vertexA.y : vertexB.y;
        yMin = yMin < vertexC.y ? yMin : vertexC.y;
        T yMax = vertexA.y > vertexB.y ? vertexA.y : vertexB.y;
        yMax = yMax > vertexC.y ? yMax : vertexC.y;

        T zMin = vertexA.z < vertexB.z ? vertexA.z : vertexB.z;
        zMin = zMin < vertexC.z ? zMin : vertexC.z;
        T zMax = vertexA.z > vertexB.z ? vertexA.z : vertexB.z;
        zMax = zMax > vertexC.z ? zMax : vertexC.z;

        if (E.x < xMin-small2 || xMax+small2 < E.x ||
            E.y < yMin-small2 || yMax+small2 < E.y ||
            E.z < zMin-small2 || zMax+small2 < E.z)
            return false;
        else
            return true;
    }

    /**
     * @brief Test if a point is within the box bounding a triangle, in the x-z plane.
     *
     * Overload taking plain @ref sxz values; bounds x and z only.
     *
     * @tparam T underlying type of sxyz objects
     * @param vertexA 1st point defining triangle
     * @param vertexB 2nd point defining triangle
     * @param vertexC 3rd point defining triangle
     * @param E point to test
     * @returns test value
     *
     * @note Conservative pre-filter: false is conclusive, true is not.
     */
    template<typename T>
    bool testInTriangleBoundingBox(const sxz<T> &vertexA,
                                   const sxz<T> &vertexB,
                                   const sxz<T> &vertexC,
                                   const sxz<T> &E) {
        T xMin = vertexA.x < vertexB.x ? vertexA.x : vertexB.x;
        xMin = xMin < vertexC.x ? xMin : vertexC.x;
        T xMax = vertexA.x > vertexB.x ? vertexA.x : vertexB.x;
        xMax = xMax > vertexC.x ? xMax : vertexC.x;
        T zMin = vertexA.z < vertexB.z ? vertexA.z : vertexB.z;
        zMin = zMin < vertexC.z ? zMin : vertexC.z;
        T zMax = vertexA.z > vertexB.z ? vertexA.z : vertexB.z;
        zMax = zMax > vertexC.z ? zMax : vertexC.z;

        if (E.x < xMin-small2 || xMax+small2 < E.x ||
            E.z < zMin-small2 || zMax+small2 < E.z)
            return false;
        else
            return true;
    }

    /**
     * @brief Compute the closest **squared** distance between a point and a line segment.
     *
     * Projects @p E onto the infinite line through the segment and clamps the
     * projection parameter to @f$[0,1]@f$, so a point beyond either end
     * measures to that endpoint rather than to the line.
     *
     * @tparam T underlying type of sxyz object
     * @tparam NODE type of Node objects
     * @param vertexA 1st point defining segment
     * @param vertexB 2nd point defining segment
     * @param E point to consider
     * @returns distance
     *
     * @warning The result is the **squared** distance — as the name says, but
     *          not as the @c \@returns line above suggests. Compare it against a
     *          squared tolerance such as @ref small2, which is what
     *          @ref testInTriangle does. No square root is taken.
     * @pre The segment has non-zero length; a degenerate segment divides by
     *      zero. Not guarded.
     */
    template<typename T, typename NODE>
    T distSqPointToSegment(const NODE *vertexA,
                           const NODE *vertexB,
                           const sxyz<T> &E) {

        T p1_p2_squareLength = (vertexB->getX() - vertexA->getX())*(vertexB->getX() - vertexA->getX()) +
        (vertexB->getY() - vertexA->getY())*(vertexB->getY() - vertexA->getY()) +
        (vertexB->getZ() - vertexA->getZ())*(vertexB->getZ() - vertexA->getZ());

        T dotProd = ((E.x - vertexA->getX())*(vertexB->getX() - vertexA->getX()) +
                     (E.y - vertexA->getY())*(vertexB->getY() - vertexA->getY()) +
                     (E.z - vertexA->getZ())*(vertexB->getZ() - vertexA->getZ())) / p1_p2_squareLength;

        if ( dotProd < 0.0 ) {
            return (E.x-vertexA->getX())*(E.x-vertexA->getX()) +
            (E.y-vertexA->getY())*(E.y-vertexA->getY()) +
            (E.z-vertexA->getZ())*(E.z-vertexA->getZ());
        } else if ( dotProd > 1.0 ) {
            return (E.x-vertexB->getX())*(E.x-vertexB->getX()) +
            (E.y-vertexB->getY())*(E.y-vertexB->getY()) +
            (E.z-vertexB->getZ())*(E.z-vertexB->getZ());
        } else {
            T p_p1_squareLength = (vertexA->getX() - E.x)*(vertexA->getX() - E.x) +
            (vertexA->getY() - E.y)*(vertexA->getY() - E.y) +
            (vertexA->getZ() - E.z)*(vertexA->getZ() - E.z);
            return p_p1_squareLength - dotProd * dotProd * p1_p2_squareLength;
        }
    }

    /**
     * @brief Compute the closest **squared** distance between a point and a line segment.
     *
     * Overload taking plain @ref sxyz values instead of Node pointers; same
     * clamped projection as
     * @ref distSqPointToSegment(const NODE*, const NODE*, const sxyz<T>&).
     *
     * @tparam T underlying type of sxyz object
     * @param vertexA 1st point defining segment
     * @param vertexB 2nd point defining segment
     * @param E point to consider
     * @returns distance
     *
     * @warning Returns the **squared** distance; compare against a squared
     *          tolerance.
     * @pre The segment has non-zero length.
     */
    template<typename T>
    T distSqPointToSegment(const sxyz<T> &vertexA,
                           const sxyz<T> &vertexB,
                           const sxyz<T> &E) {

        T p1_p2_squareLength = (vertexB.x - vertexA.x)*(vertexB.x - vertexA.x) +
        (vertexB.y - vertexA.y)*(vertexB.y - vertexA.y) +
        (vertexB.z - vertexA.z)*(vertexB.z - vertexA.z);

        T dotProd = ((E.x - vertexA.x)*(vertexB.x - vertexA.x) +
                     (E.y - vertexA.y)*(vertexB.y - vertexA.y) +
                     (E.z - vertexA.z)*(vertexB.z - vertexA.z)) / p1_p2_squareLength;

        if ( dotProd < 0.0 ) {
            return (E.x-vertexA.x)*(E.x-vertexA.x) +
            (E.y-vertexA.y)*(E.y-vertexA.y) +
            (E.z-vertexA.z)*(E.z-vertexA.z);
        } else if ( dotProd > 1.0 ) {
            return (E.x-vertexB.x)*(E.x-vertexB.x) +
            (E.y-vertexB.y)*(E.y-vertexB.y) +
            (E.z-vertexB.z)*(E.z-vertexB.z);
        } else {
            T p_p1_squareLength = (vertexA.x - E.x)*(vertexA.x - E.x) +
            (vertexA.y - E.y)*(vertexA.y - E.y) +
            (vertexA.z - E.z)*(vertexA.z - E.z);
            return p_p1_squareLength - dotProd * dotProd * p1_p2_squareLength;
        }
    }

    /**
     * @brief Compute the closest **signed** distance between a point and a plane.
     *
     * The plane is given by three points; its unit normal comes from
     * @f$(b-a)\times(c-a)@f$ and the distance is @f$n\cdot(a-pt)@f$.
     *
     * @tparam T underlying type of sxyz object
     * @tparam NODE type of Node objects
     * @param a 1st point defining plane
     * @param b 2nd point defining plane
     * @param c 3rd point defining plane
     * @param pt point to consider
     * @returns distance
     *
     * @note The result is **signed** — its sign tells you which side of the
     *       plane @p pt is on, which is what makes it usable for
     *       above/below tests. Take the absolute value for a magnitude.
     *       The sign follows the winding of @p a, @p b, @p c.
     * @pre The three points are not collinear; otherwise the normal is zero and
     *      normalising it is undefined.
     */
    template<typename T, typename NODE>
    T distPointToPlane(const NODE *a,
                       const NODE *b,
                       const NODE *c,
                       const sxyz<T> &pt) {
        sxyz<T> ab = {b->getX()-a->getX(), b->getY()-a->getY(), b->getZ()-a->getZ()};
        sxyz<T> ac = {c->getX()-a->getX(), c->getY()-a->getY(), c->getZ()-a->getZ()};

        // triangle normal
        sxyz<T> n = cross(ab, ac);
        n.normalize();
        return dot(n, {a->getX() - pt.x, a->getY() - pt.y, a->getZ() - pt.z});
    }

    /**
     * @brief Compute the closest distance between a point and a line, in the x-z plane.
     *
     * Forms the implicit line equation @f$ax+bz+c=0@f$ through @p pt1 and
     * @p pt2 and evaluates @f$|ax+bz+c|/\sqrt{a^2+b^2}@f$.
     *
     * @tparam T underlying type of sxyz object
     * @param pt1 1st point defining line
     * @param pt2 2nd point defining line
     * @param pt point to consider
     * @returns distance
     *
     * @note Measures to the **infinite line**, not to the segment — unlike
     *       @ref distSqPointToSegment, which clamps to the endpoints. Returns
     *       a true distance, not a squared one.
     * @pre @p pt1 and @p pt2 are distinct; coincident points give
     *      @f$a=b=0@f$ and a division by zero.
     * @warning Calls @c abs and @c sqrt unqualified; see the note on
     *          @ref areCollinear(const sxz<T>&, const NODE&, const NODE&).
     */
    template<typename T>
    T distPointToLine(const sxz<T>& pt1, const sxz<T>& pt2, const sxz<T>& pt) {
        // get eq of line
        T a = pt2.z - pt1.z;
        T b = pt1.x - pt2.x;
        T c = pt2.x*pt1.z - pt1.x*pt2.z;
        // compute distance
        return abs(a*pt.x + b*pt.z + c) / sqrt(a*a + b*b);
    }

    /**
     * @brief Test if a point lies within a triangle.
     *
     * Three stages, cheapest first: reject on the bounding box
     * (@ref testInTriangleBoundingBox), then accept if the barycentric weights
     * are all non-negative, and finally accept if the point is within
     * @ref small2 of any of the three edges. That last stage is what makes a
     * point sitting exactly on an edge or vertex — the common case when a
     * raypath crosses from one mesh cell into the next — count as inside.
     *
     * @tparam T underlying type of sxyz object
     * @tparam NODE type of Node objects
     * @param vertexA 1st point defining triangle
     * @param vertexB 2nd point defining triangle
     * @param vertexC 3rd point defining triangle
     * @param E point to consider
     * @returns test results
     *
     * @note Because of the edge tolerance, a point on a shared edge tests
     *       inside **both** adjoining triangles. Callers scanning for the
     *       containing cell should take the first hit rather than assume it is
     *       unique.
     * @note @p E is not required to be coplanar with the triangle; it is
     *       effectively projected, so a point well off the plane can test
     *       inside. Pair with @ref distPointToPlane where that matters.
     */
    template<typename T, typename NODE>
    bool testInTriangle(const NODE *vertexA,
                        const NODE *vertexB,
                        const NODE *vertexC,
                        const sxyz<T> &E) {

        if ( ! testInTriangleBoundingBox(vertexA, vertexB, vertexC, E) )
            return false;

        T u, v, w;
        barycentric(vertexA, vertexB, vertexC, E, u, v, w);
        if ( v >= 0.0 && w >= 0.0 && (v + w) <= 1.0 )
            return true;

        if (distSqPointToSegment(vertexA, vertexB, E) <= small2)
            return true;
        if (distSqPointToSegment(vertexA, vertexC, E) <= small2)
            return true;
        if (distSqPointToSegment(vertexB, vertexC, E) <= small2)
            return true;

        return false;
    }

    /**
     * @brief Test if a point lies within a triangle.
     *
     * Overload taking plain @ref sxyz values instead of Node pointers; the
     * three-stage test is identical to
     * @ref testInTriangle(const NODE*, const NODE*, const NODE*, const sxyz<T>&),
     * whose documentation describes the edge tolerance and its consequences.
     *
     * @tparam T underlying type of sxyz object
     * @param vertexA 1st point defining triangle
     * @param vertexB 2nd point defining triangle
     * @param vertexC 3rd point defining triangle
     * @param E point to consider
     * @returns test results
     */
    template<typename T>
    bool testInTriangle(const sxyz<T> &vertexA,
                        const sxyz<T> &vertexB,
                        const sxyz<T> &vertexC,
                        const sxyz<T> &E) {

        if ( ! testInTriangleBoundingBox(vertexA, vertexB, vertexC, E) )
            return false;

        T u, v, w;
        barycentric(vertexA, vertexB, vertexC, E, u, v, w);
        if ( v >= 0.0 && w >= 0.0 && (v + w) <= 1.0 )
            return true;

        if (distSqPointToSegment(vertexA, vertexB, E) <= small2)
            return true;
        if (distSqPointToSegment(vertexA, vertexC, E) <= small2)
            return true;
        if (distSqPointToSegment(vertexB, vertexC, E) <= small2)
            return true;

        return false;
    }

    /**
     * @brief Compute the Moore-Penrose pseudo-inverse of a matrix by SVD.
     *
     * Decomposes @f$A = U\Sigma V^{T}@f$ and forms
     * @f$A^{+} = V\Sigma^{+}U^{T}@f$, where @f$\Sigma^{+}@f$ inverts the
     * non-zero singular values. Used to solve the over- or under-determined
     * least-squares systems that arise when fitting a traveltime gradient.
     *
     * @tparam T scalar type of the matrices.
     * @param[in] A matrix to invert; any shape, square or not.
     * @param[out] pi the pseudo-inverse, resized as needed.
     * @returns Rank of @p A, as reported by the decomposition.
     *
     * @warning The @f$\Sigma^{+}@f$ diagonal is built from
     *          @c nonzeroSingularValues() with **no relative threshold**, so a
     *          singular value that is merely tiny rather than exactly zero is
     *          still inverted, and its reciprocal dominates the result. For a
     *          near-rank-deficient @p A this is numerically unstable; check the
     *          returned rank.
     * @note Uses @c JacobiSVD, which is accurate but @f$O(n^3)@f$ with a large
     *       constant — fine for the small systems here, poor for big ones.
     * @sa Grad.h
     */
    template<typename T>
    Eigen::Index pseudoInverse(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A,
                               Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& pi) {

        Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> S;
        S.resize(svd.nonzeroSingularValues(), svd.nonzeroSingularValues());
        S.setZero();

        for ( ptrdiff_t n=0; n<svd.nonzeroSingularValues(); ++n ) {
            S(n, n) = 1.0 / svd.singularValues()(n);
        }
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> U = svd.matrixU().block(0, 0, svd.matrixU().rows(), svd.nonzeroSingularValues());
        pi = svd.matrixV() * S * U.transpose();

        return svd.rank();
    }

    /**
     * @brief Build reflectors from interfaces between two lithologies.
     *
     * Each named 2-D physical entity of the Gmsh mesh becomes one reflector.
     * Its triangles are discretised into a point set — the three vertices, then
     * @p nsecondary evenly spaced points along each edge, then a triangular
     * lattice of interior face points — and that set is handed to a
     * ttcr::Rcv, so a reflector is modelled as a dense receiver array and
     * reflected arrivals fall out of the ordinary traveltime computation.
     *
     * Points are accumulated in a `std::set`, so vertices and edge points
     * shared between adjacent triangles are stored once rather than duplicated.
     *
     * @tparam T underlying type of sxyz & Rcv objects
     * @param[in] reader reader used to extract indicides of nodes makign the reflectors
     * @param[in] nodes grid nodes
     * @param[in] nsrc number of sources to model
     * @param[in] nsecondary number of secondary nodes
     * @param[out] reflectors vector of Rcv objects making the reflectors
     *
     * @pre @p nsecondary is at least 1.
     * @throws std::invalid_argument if @p nsecondary is 0. The face-point loop
     *         is bounded by `ncut = nsecondary - 1` on an unsigned type, so 0
     *         would wrap it to `SIZE_MAX`, divide by a zero segment count and
     *         try to insert an unbounded number of points. This is reachable in
     *         normal use — ttcr::input_parameters::nn defaults to zero and the
     *         caller in grids.h passes `par.nn[0]` whenever
     *         ttcr::input_parameters::processReflectors is set — so it is
     *         rejected up front rather than left to hang.
     * @warning Each ttcr::Rcv is constructed from the reflector's **name**,
     *          which lands in its filename field. ttcr::Rcv::save_rcvfile would
     *          therefore write to a file named after the reflector, in the
     *          working directory.
     * @warning Reports a mismatch between reflector names and indices on
     *          @c std::cerr and calls @c exit(1) rather than throwing.
     * @note Reflectors are appended to @p reflectors, which is not cleared
     *       first.
     * @sa MSHReader::getPhysicalNames, Rcv
     */
    template<typename T>
    void buildReflectors(const MSHReader &reader,
                         const std::vector<sxyz<T>> &nodes,
                         const size_t nsrc,
                         const size_t nsecondary,
                         std::vector<Rcv<T>> &reflectors) {

        if ( nsecondary == 0 ) {
            // ncut below would wrap to SIZE_MAX on this unsigned type, and the
            // face-point loop would then divide by a zero segment count and try
            // to insert an unbounded number of points.
            throw std::invalid_argument("Error: building reflectors requires at "
                                        "least one secondary node (set \"secondary "
                                        "nodes\" in the parameter file).");
        }

        std::vector<std::string> reflector_names = reader.getPhysicalNames(2);
        std::vector<size_t> indices = reader.getPhysicalIndices(2);

        if ( reflector_names.size() != indices.size() ) {
            std::cerr << "Error - definition of reflectors\n";
            exit(1);
        }

        std::vector<triangleElem<uint32_t>> triangles;
        reader.readTriangleElements(triangles);

        sxyz<T> pt1, pt2, pt3, pt4, pt5, pt6, d1, d2;
        size_t ncut = nsecondary - 1;

        for ( size_t ni=0; ni<indices.size(); ++ni ) {

            reflectors.push_back( reflector_names[ni] );

            std::set<sxyz<T>> refl_pts;  // use set to avoid duplicate points
            typename std::set<sxyz<T>>::iterator it;

            for ( size_t nt=0; nt<triangles.size(); ++nt ) {
                if ( indices[ni] == triangles[nt].physical_entity ) {
                    pt1 = nodes[ triangles[nt].i[0] ];
                    pt2 = nodes[ triangles[nt].i[1] ];
                    pt3 = nodes[ triangles[nt].i[2] ];

                    // edge nodes
                    // 1st edge
                    d1.x = (pt2.x-pt1.x)/(nsecondary+1);
                    d1.y = (pt2.y-pt1.y)/(nsecondary+1);
                    d1.z = (pt2.z-pt1.z)/(nsecondary+1);

                    refl_pts.insert( pt1 );
                    for ( size_t n2=0; n2<nsecondary; ++n2 ) {
                        pt4.x = pt1.x+(1+n2)*d1.x;
                        pt4.y = pt1.y+(1+n2)*d1.y;
                        pt4.z = pt1.z+(1+n2)*d1.z;
                        refl_pts.insert( pt4 );
                    }
                    refl_pts.insert( pt2 );

                    // 2nd edge
                    d1.x = (pt3.x-pt2.x)/(nsecondary+1);
                    d1.y = (pt3.y-pt2.y)/(nsecondary+1);
                    d1.z = (pt3.z-pt2.z)/(nsecondary+1);

                    for ( size_t n2=0; n2<nsecondary; ++n2 ) {
                        pt4.x = pt2.x+(1+n2)*d1.x;
                        pt4.y = pt2.y+(1+n2)*d1.y;
                        pt4.z = pt2.z+(1+n2)*d1.z;
                        refl_pts.insert( pt4 );
                    }
                    refl_pts.insert( pt3 );

                    // 3rd edge
                    d1.x = (pt1.x-pt3.x)/(nsecondary+1);
                    d1.y = (pt1.y-pt3.y)/(nsecondary+1);
                    d1.z = (pt1.z-pt3.z)/(nsecondary+1);

                    for ( size_t n2=0; n2<nsecondary; ++n2 ) {
                        pt4.x = pt3.x+(1+n2)*d1.x;
                        pt4.y = pt3.y+(1+n2)*d1.y;
                        pt4.z = pt3.z+(1+n2)*d1.z;
                        refl_pts.insert( pt4 );
                    }

                    // face nodes
                    d2.x = (pt1.x-pt2.x)/(nsecondary+1);
                    d2.y = (pt1.y-pt2.y)/(nsecondary+1);
                    d2.z = (pt1.z-pt2.z)/(nsecondary+1);

                    for ( size_t n=0; n<ncut; ++n ) {

                        pt4.x = pt3.x+(1+n)*d1.x;
                        pt4.y = pt3.y+(1+n)*d1.y;
                        pt4.z = pt3.z+(1+n)*d1.z;

                        pt5.x = pt2.x+(1+n)*d2.x;
                        pt5.y = pt2.y+(1+n)*d2.y;
                        pt5.z = pt2.z+(1+n)*d2.z;

                        size_t nseg = ncut+1-n;

                        sxyz<T> d = { (pt5.x-pt4.x)/nseg,
                            (pt5.y-pt4.y)/nseg,
                            (pt5.z-pt4.z)/nseg };

                        for ( size_t n2=0; n2<nseg-1; ++n2 ) {
                            pt6.x = pt1.x+(1+n2)*d.x;
                            pt6.y = pt1.y+(1+n2)*d.y;
                            pt6.z = pt1.z+(1+n2)*d.z;
                            refl_pts.insert( pt6 );
                        }
                    }
                }
            }
            for (it=refl_pts.begin(); it!=refl_pts.end(); ++it) {
                reflectors.back().add_coord( *it );
            }
            reflectors.back().init_tt( nsrc );
        }
    }

    /**
     * @brief Save 3D raypaths in a vtk file.
     *
     * Writes the rays as a VTK PolyData of polylines, one line per raypath,
     * in binary XML (@c .vtp).
     *
     * @tparam T underlying type of sxyz objects
     * @param fname name of file for saving paths
     * @param r_data raypath coordinates
     *
     * @note Compiled to an empty function unless @c VTK is defined — without
     *       VTK support the call silently writes nothing.
     */
    template<typename T>
    void saveRayPaths(const std::string &fname,
                      const std::vector<std::vector<sxyz<T>>> &r_data) {

#ifdef VTK
        vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();
        vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();

        size_t npts=0;
        for ( size_t n=0; n<r_data.size(); ++n ) {
            npts += r_data[n].size();
        }
        pts->SetNumberOfPoints(npts);
        for ( size_t n=0, npts=0; n<r_data.size(); ++n ) {
            for ( size_t np=0; np<r_data[n].size(); ++np, ++npts ) {
                pts->InsertPoint(npts, r_data[n][np].x, r_data[n][np].y, r_data[n][np].z);
            }
        }
        polydata->SetPoints(pts);

        for ( size_t n=0, npts=0; n<r_data.size(); ++n ) {
            vtkSmartPointer<vtkPolyLine> line = vtkSmartPointer<vtkPolyLine>::New();
            line->GetPointIds()->SetNumberOfIds( r_data[n].size() );
            for ( size_t np=0; np<r_data[n].size(); ++np, ++npts ) {
                line->GetPointIds()->SetId(np, npts);
            }
            cellarray->InsertNextCell(line);
        }
        polydata->SetLines(cellarray);

        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer->SetFileName( fname.c_str() );
        writer->SetInputData( polydata );
        writer->SetDataModeToBinary();
        writer->Update();
#endif

    }

    /**
     * @brief Save 2D raypaths in a vtk file.
     *
     * As the 3-D overload, but points are emitted as (x, 0, z) — the missing y
     * is written as zero — and degenerate raypaths of fewer than two points are
     * skipped, since they cannot form a drawable polyline.
     *
     * @tparam T underlying type of sxz objects
     * @param fname name of file for saving paths
     * @param r_data raypath coordinates
     *
     * @note Compiled to an empty function unless @c VTK is defined.
     * @note A skipped path still contributes its points to the point array —
     *       they are written, just not joined into a line — so the file may
     *       hold points that belong to no polyline.
     */
    template<typename T>
    void saveRayPaths(const std::string &fname,
                      const std::vector<std::vector<sxz<T>>> &r_data) {

#ifdef VTK
        vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();
        vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();

        size_t npts=0;
        for ( size_t n=0; n<r_data.size(); ++n ) {
            npts += r_data[n].size();
        }
        pts->SetNumberOfPoints(npts);
        for ( size_t n=0, npts=0; n<r_data.size(); ++n ) {
            for ( size_t np=0; np<r_data[n].size(); ++np, ++npts ) {
                pts->InsertPoint(npts, r_data[n][np].x, 0.0, r_data[n][np].z);
            }
        }
        polydata->SetPoints(pts);

        for ( size_t n=0, npts=0; n<r_data.size(); ++n ) {
            if ( r_data[n].size() < 2 ) {
                // A path of fewer than 2 points cannot form a polyline, but its
                // points were still written to the point array above: step over
                // them so the following lines keep referring to the right ones.
                npts += r_data[n].size();
                continue;
            }
            vtkSmartPointer<vtkPolyLine> line = vtkSmartPointer<vtkPolyLine>::New();
            line->GetPointIds()->SetNumberOfIds( r_data[n].size() );
            for ( size_t np=0; np<r_data[n].size(); ++np, ++npts ) {
                line->GetPointIds()->SetId(np, npts);
            }
            cellarray->InsertNextCell(line);
        }
        polydata->SetLines(cellarray);

        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer->SetFileName( fname.c_str() );
        writer->SetInputData( polydata );
        writer->SetDataModeToBinary();
        writer->Update();
#endif

    }

    /**
     * @brief Create a string from any streamable value.
     *
     * Formats @p value through a @c std::ostringstream.
     *
     * @tparam T type of object for which string is created
     * @param value value to convert to string
     * @returns string representation of value
     *
     * @note Predates the availability of @c std::to_string here and, unlike it,
     *       works for any type with an @c operator<<. Being in namespace
     *       @c ttcr it does not conflict, but an unqualified @c to_string call
     *       inside the namespace will find this one by ordinary lookup.
     * @note Uses the stream's default formatting — six significant digits for
     *       floating-point values, which loses precision. Use an explicit
     *       stream with @c setprecision where the full value matters.
     */
    template<typename T>
    std::string to_string( const T & value )
    {
        // utiliser un flux de sortie pour créer la chaîne
        std::ostringstream oss;
        // écrire la valeur dans le flux
        oss << value;
        // renvoyer une string
        return oss.str();
    }

}

#endif

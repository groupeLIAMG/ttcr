//
//  Interpolator.h
//  ttcr
//
//  Created by Bernard Giroux on 2012-11-23.
//  Copyright (c) 2012 Bernard Giroux. All rights reserved.
//  Copyright (c) 2018 Bernard Giroux, Maher Nasr. All rights reserved.
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
 * \file Interpolator.h
 * \brief Static interpolation routines used to evaluate fields inside grid cells
 *
 * ttcr::Interpolator is a stateless collection of interpolation formulas.  The
 * grid classes call them to evaluate slowness, velocity or traveltime at a
 * point that does not coincide with a node, typically while propagating a
 * raypath through a cell.
 *
 * \sa Grid2Drn, Grid3Drn, Grid2Dun, Grid3Dun, Grad
 */

#ifndef ttcr_Interpolator_h
#define ttcr_Interpolator_h

#include <array>
#include <cassert>
#include <cmath>
#include <vector>

#include "ttcr_t.h"

namespace ttcr {

    /**
     * \brief Collection of static interpolation routines
     *
     * The class is never instantiated; it only groups static member functions
     * under a common numeric type.  Three families of routines are provided.
     *
     * **Interpolation on a regular grid** — linear(), bilinear() and
     * trilinear() take plain C arrays: the coordinate arrays hold the
     * evaluation point at index 0 and the two bracketing grid coordinates at
     * indices 1 and 2, while the value array holds the field at the corners of
     * the cell.
     *
     * **Interpolation from a set of nodes** — the remaining routines take node
     * objects, and weight them according to the position of the evaluation
     * point.  The naming reflects the geometry rather than the mathematics:
     * despite their names, `barycentricTriangle*` (triangle in the x-z plane),
     * `bilinearTriangle*` (triangle in 3D) and `trilinearTriangle*`
     * (tetrahedron) all compute barycentric weights, i.e. ratios of lengths,
     * areas and volumes respectively.  inverseDistance() instead weights each
     * node by the reciprocal of its distance to the evaluation point.
     *
     * **Naming suffixes** — for most geometries three variants exist:
     * - the plain name interpolates *slowness* linearly, returning
     *   \f$\sum_i w_i s_i\f$;
     * - the `...Vel` variant interpolates *velocity* linearly and returns the
     *   corresponding slowness, \f$1 / \sum_i w_i v_i\f$, computed as
     *   \f$1 / \sum_i w_i / s_i\f$.  The two give different results, and the
     *   grid classes choose between them at run time (their `processVel` flag);
     * - the `...Time` variant applies the same weights to the traveltimes
     *   stored in the nodes, and the `...Weight` variant merely returns the
     *   weights.
     *
     * \tparam T numeric type of the coordinates and of the interpolated field
     *           (float or double)
     *
     * \note NODE template parameters are duck-typed: depending on the routine,
     * a node must provide `getX()`, `getY()`, `getZ()`, `getNodeSlowness()`,
     * `getTT()` and/or `getDistance()`.
     */
    template<class T> class Interpolator
    {
    public:
        /**
         * \brief Linear interpolation between two samples
         *
         * Evaluates the field at `x[0]`, knowing `s[0]` at `x[1]` and `s[1]`
         * at `x[2]`:
         * \f[ s = \frac{s_0 (x_2 - x_0) + s_1 (x_0 - x_1)}{x_2 - x_1} \f]
         *
         * \pre `x[1] != x[2]`; this is not checked.
         *
         * \param x array of 3 abscissas: evaluation point, then the two
         *          bracketing sample abscissas
         * \param s array of 2 values, at `x[1]` and `x[2]`
         * \return the interpolated value at `x[0]`
         */
        inline static T linear(const T x[], const T s[]) {

            // evaluate s @ x[0]
            // with
            // s[0] @ x[1]
            // s[1] @ x[2]

            return (s[0]*(x[2]-x[0]) + s[1]*(x[0]-x[1]))/(x[2]-x[1]);
        }

        /**
         * \brief Bilinear interpolation over a rectangular cell
         *
         * Evaluates the field at `(x[0], y[0])` from its values at the four
         * corners of the rectangle \f$[x_1, x_2] \times [y_1, y_2]\f$, given
         * in the order (x1,y1), (x1,y2), (x2,y1), (x2,y2).
         *
         * \pre `x[1] != x[2]` and `y[1] != y[2]`; this is not checked.
         *
         * \param x array of 3 abscissas: evaluation point, then the two cell edges
         * \param y array of 3 ordinates: evaluation point, then the two cell edges
         * \param s array of 4 corner values, ordered as described above
         * \return the interpolated value at `(x[0], y[0])`
         */
        inline static T bilinear(const T x[], const T y[], const T s[]) {

            // evaluate s @ (x[0], y[0])
            // with
            // s[0] @ (x[1], y[1])
            // s[1] @ (x[1], y[2])
            // s[2] @ (x[2], y[1])
            // s[3] @ (x[2], y[2])

            return (s[0]*(x[2]-x[0])*(y[2]-y[0]) +
                    s[1]*(x[2]-x[0])*(y[0]-y[1]) +
                    s[2]*(x[0]-x[1])*(y[2]-y[0]) +
                    s[3]*(x[0]-x[1])*(y[0]-y[1]))
            /((x[2]-x[1])*(y[2]-y[1]));
        }

        /**
         * \brief Trilinear interpolation over a rectangular cuboid cell
         *
         * Evaluates the field at `(x[0], y[0], z[0])` from its values at the
         * eight corners of the cuboid \f$[x_1, x_2] \times [y_1, y_2] \times
         * [z_1, z_2]\f$.  The corner values are ordered with z varying
         * fastest, then y, then x: (x1,y1,z1), (x1,y1,z2), (x1,y2,z1),
         * (x1,y2,z2), (x2,y1,z1), (x2,y1,z2), (x2,y2,z1), (x2,y2,z2).
         *
         * \pre `x[1] != x[2]`, `y[1] != y[2]` and `z[1] != z[2]`; this is not
         * checked.
         *
         * \param x array of 3 abscissas: evaluation point, then the two cell faces
         * \param y array of 3 ordinates: evaluation point, then the two cell faces
         * \param z array of 3 elevations: evaluation point, then the two cell faces
         * \param s array of 8 corner values, ordered as described above
         * \return the interpolated value at `(x[0], y[0], z[0])`
         */
        inline static T trilinear(const T x[], const T y[], const T z[], const T s[]) {

            // evaluate s @ (x[0], y[0], z[0])
            // with
            // s[0] @ (x[1], y[1], z[1])
            // s[1] @ (x[1], y[1], z[2])
            // s[2] @ (x[1], y[2], z[1])
            // s[3] @ (x[1], y[2], z[2])
            // s[4] @ (x[2], y[1], z[1])
            // s[5] @ (x[2], y[1], z[2])
            // s[6] @ (x[2], y[2], z[1])
            // s[7] @ (x[2], y[2], z[2])

            return (s[0]*(x[2]-x[0])*(y[2]-y[0])*(z[2]-z[0]) +
                    s[1]*(x[2]-x[0])*(y[2]-y[0])*(z[0]-z[1]) +
                    s[2]*(x[2]-x[0])*(y[0]-y[1])*(z[2]-z[0]) +
                    s[3]*(x[2]-x[0])*(y[0]-y[1])*(z[0]-z[1]) +
                    s[4]*(x[0]-x[1])*(y[2]-y[0])*(z[2]-z[0]) +
                    s[5]*(x[0]-x[1])*(y[2]-y[0])*(z[0]-z[1]) +
                    s[6]*(x[0]-x[1])*(y[0]-y[1])*(z[2]-z[0]) +
                    s[7]*(x[0]-x[1])*(y[0]-y[1])*(z[0]-z[1]))
            /((x[2]-x[1])*(y[2]-y[1])*(z[2]-z[1]));
        }


        /**
         * \brief Inverse-distance (Shepard) interpolation of slowness, node form
         *
         * Weights each interpolating node by the reciprocal of its distance to
         * \p node:
         * \f[ s = \frac{\sum_i s_i / d_i}{\sum_i 1 / d_i} \f]
         *
         * \pre \p inodes is not empty and \p node coincides with none of its
         * entries; neither is checked, and a zero distance produces a division
         * by zero.
         *
         * \tparam NODE node type, providing `getDistance()` and `getNodeSlowness()`
         * \param node   point at which the slowness is evaluated
         * \param inodes nodes used for the interpolation
         * \return the interpolated slowness
         */
        template<typename NODE>
        static T inverseDistance(const NODE &node,
                                 const std::vector<NODE*> &inodes) {

            T num=0.;
            T den=0.;
            T w;

            for ( size_t n=0; n<inodes.size(); ++n ) {
                w = 1./inodes[n]->getDistance( node );
                num += w*inodes[n]->getNodeSlowness();
                den += w;
            }

            return num/den;
        }

        /**
         * \brief Inverse-distance (Shepard) interpolation of slowness, coordinate form
         *
         * Identical to the node overload, but the evaluation point is given as
         * a plain coordinate struct rather than as a node.
         *
         * \pre \p inodes is not empty and \p node coincides with none of its
         * entries; neither is checked, and a zero distance produces a division
         * by zero.
         *
         * \tparam NODE node type, providing `getDistance()` and `getNodeSlowness()`
         * \tparam S    coordinate struct type (sxz<T> or sxyz<T>)
         * \param node   point at which the slowness is evaluated
         * \param inodes nodes used for the interpolation
         * \return the interpolated slowness
         */
        template<typename NODE, typename S>
        static T inverseDistance(const S& node,
                                 const std::vector<NODE*> &inodes) {

            T num=0.;
            T den=0.;
            T w;

            for ( size_t n=0; n<inodes.size(); ++n ) {
                w = 1./inodes[n]->getDistance( node );
                num += w*inodes[n]->getNodeSlowness();
                den += w;
            }

            return num/den;
        }

        /**
         * \brief Linear interpolation of slowness along a segment in 3D
         *
         * The weight of each end node is the normalized distance from the
         * evaluation point to the *other* end, so that the result reduces to
         * the slowness of an end node when \p node coincides with it.
         *
         * \pre \p node lies on the segment joining \p node1 and \p node2, and
         * the two nodes are distinct; neither is checked.
         *
         * \tparam NODE node type, providing `getX()`, `getY()`, `getZ()` and
         *              `getNodeSlowness()`
         * \param node  point at which the slowness is evaluated
         * \param node1 first end of the segment
         * \param node2 second end of the segment
         * \return the interpolated slowness
         * \sa linearVel() to interpolate velocity instead
         */
        template<typename NODE>
        static T linear(const sxyz<T>& node,
                        const NODE& node1,
                        const NODE& node2) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> IA = {node1.getX()-node.x, node1.getY()-node.y, node1.getZ()-node.z};
            sxyz<T> IB = {node2.getX()-node.x, node2.getY()-node.y, node2.getZ()-node.z};

            T nAB = norm(AB);
            T w1 = norm(IA)/nAB;
            T w2 = norm(IB)/nAB;

            return w2*node1.getNodeSlowness() + w1*node2.getNodeSlowness();
        }

        /**
         * \brief Linear interpolation of velocity along a segment in 3D
         *
         * Uses the same weights as linear(), but interpolates the velocities
         * of the end nodes and returns the slowness of the result:
         * \f[ s = \left( w_2 / s_1 + w_1 / s_2 \right)^{-1} \f]
         *
         * \pre \p node lies on the segment joining \p node1 and \p node2, the
         * two nodes are distinct, and their slownesses are non-zero; none of
         * this is checked.
         *
         * \tparam NODE node type, providing `getX()`, `getY()`, `getZ()` and
         *              `getNodeSlowness()`
         * \param node  point at which the slowness is evaluated
         * \param node1 first end of the segment
         * \param node2 second end of the segment
         * \return the slowness corresponding to the interpolated velocity
         * \sa linear() to interpolate slowness instead
         */
        template<typename NODE>
        static T linearVel(const sxyz<T>& node,
                           const NODE& node1,
                           const NODE& node2) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> IA = {node1.getX()-node.x, node1.getY()-node.y, node1.getZ()-node.z};
            sxyz<T> IB = {node2.getX()-node.x, node2.getY()-node.y, node2.getZ()-node.z};

            T nAB = norm(AB);
            T w1 = norm(IA)/nAB;
            T w2 = norm(IB)/nAB;

            return 1.0/(w2/node1.getNodeSlowness() + w1/node2.getNodeSlowness());
        }

        /**
         * \brief Barycentric interpolation of slowness in a triangle of the x-z plane
         *
         * The weights are the barycentric coordinates of \p node in the
         * triangle, obtained from the ratios of the sub-triangle areas; the
         * third one is derived as \f$w_3 = 1 - w_1 - w_2\f$.
         *
         * \pre the three nodes are not collinear (the denominator, twice the
         * signed area of the triangle, would vanish) and \p node lies inside
         * the triangle; neither is checked.
         *
         * \tparam NODE node type, providing `getX()`, `getZ()` and `getNodeSlowness()`
         * \param node  point at which the slowness is evaluated
         * \param node1 first vertex of the triangle
         * \param node2 second vertex of the triangle
         * \param node3 third vertex of the triangle
         * \return the interpolated slowness
         * \sa barycentricTriangleVel() to interpolate velocity instead
         */
        template<typename NODE>
        static T barycentricTriangle(const sxz<T>& node,
                                     const NODE& node1,
                                     const NODE& node2,
                                     const NODE& node3) {
            T den = (node2.getZ()-node3.getZ())*(node1.getX()-node3.getX()) +
            (node3.getX()-node2.getX())*(node1.getZ()-node3.getZ());
            T w1 = ((node2.getZ()-node3.getZ())*(node.x-node3.getX()) +
                    (node3.getX()-node2.getX())*(node.z-node3.getZ())) / den;
            T w2 = ((node3.getZ()-node1.getZ())*(node.x-node3.getX()) +
                    (node1.getX()-node3.getX())*(node.z-node3.getZ())) / den;
            T w3 = 1. - w1 - w2;

            return w1*node1.getNodeSlowness() + w2*node2.getNodeSlowness() + w3*node3.getNodeSlowness();
        }

        /**
         * \brief Barycentric interpolation of slowness in a triangle, vector form
         *
         * Identical to the three-node overload, the vertices being taken from
         * the first three entries of \p inodes.
         *
         * \pre \p inodes holds at least 3 nodes, they are not collinear, and
         * \p node lies inside the triangle; none of this is checked.
         *
         * \tparam NODE node type, providing `getX()`, `getZ()` and `getNodeSlowness()`
         * \param node   point at which the slowness is evaluated
         * \param inodes vertices of the triangle
         * \return the interpolated slowness
         */
        template<typename NODE>
        static T barycentricTriangle(const sxz<T>& node,
                                     const std::vector<NODE*> &inodes) {
            T den = (inodes[1]->getZ()-inodes[2]->getZ())*(inodes[0]->getX()-inodes[2]->getX()) +
            (inodes[2]->getX()-inodes[1]->getX())*(inodes[0]->getZ()-inodes[2]->getZ());
            T w1 = ((inodes[1]->getZ()-inodes[2]->getZ())*(node.x-inodes[2]->getX()) +
                    (inodes[2]->getX()-inodes[1]->getX())*(node.z-inodes[2]->getZ())) / den;
            T w2 = ((inodes[2]->getZ()-inodes[0]->getZ())*(node.x-inodes[2]->getX()) +
                    (inodes[0]->getX()-inodes[2]->getX())*(node.z-inodes[2]->getZ())) / den;
            T w3 = 1. - w1 - w2;

            return w1*inodes[0]->getNodeSlowness() + w2*inodes[1]->getNodeSlowness() + w3*inodes[2]->getNodeSlowness();
        }

        /**
         * \brief Barycentric interpolation of velocity in a triangle of the x-z plane
         *
         * Uses the same barycentric weights as barycentricTriangle(), but
         * interpolates the velocities of the vertices and returns the slowness
         * of the result, \f$s = \left( \sum_i w_i / s_i \right)^{-1}\f$.
         *
         * \pre the three nodes are not collinear, \p node lies inside the
         * triangle, and the vertex slownesses are non-zero; none of this is
         * checked.
         *
         * \tparam NODE node type, providing `getX()`, `getZ()` and `getNodeSlowness()`
         * \param node  point at which the slowness is evaluated
         * \param node1 first vertex of the triangle
         * \param node2 second vertex of the triangle
         * \param node3 third vertex of the triangle
         * \return the slowness corresponding to the interpolated velocity
         * \sa barycentricTriangle() to interpolate slowness instead
         */
        template<typename NODE>
        static T barycentricTriangleVel(const sxz<T>& node,
                                        const NODE& node1,
                                        const NODE& node2,
                                        const NODE& node3) {
            T den = (node2.getZ()-node3.getZ())*(node1.getX()-node3.getX()) +
            (node3.getX()-node2.getX())*(node1.getZ()-node3.getZ());
            T w1 = ((node2.getZ()-node3.getZ())*(node.x-node3.getX()) +
                    (node3.getX()-node2.getX())*(node.z-node3.getZ())) / den;
            T w2 = ((node3.getZ()-node1.getZ())*(node.x-node3.getX()) +
                    (node1.getX()-node3.getX())*(node.z-node3.getZ())) / den;
            T w3 = 1. - w1 - w2;

            return (1.0/(w1/node1.getNodeSlowness() + w2/node2.getNodeSlowness() + w3/node3.getNodeSlowness()));
        }

        /**
         * \brief Areal interpolation of velocity in a triangle in 3D
         *
         * The weight of each vertex is the area of the sub-triangle opposite
         * to it, divided by the area of the whole triangle, both obtained from
         * cross products.  The velocities of the vertices are interpolated and
         * the slowness of the result is returned.
         *
         * \pre the three nodes are not collinear, \p node lies in the plane of
         * the triangle and inside it, and the vertex slownesses are non-zero;
         * none of this is checked.
         *
         * \tparam NODE node type, providing `getX()`, `getY()`, `getZ()` and
         *              `getNodeSlowness()`
         * \param node  point at which the slowness is evaluated
         * \param node1 first vertex of the triangle
         * \param node2 second vertex of the triangle
         * \param node3 third vertex of the triangle
         * \return the slowness corresponding to the interpolated velocity
         * \sa bilinearTriangle() to interpolate slowness instead
         */
        template<typename NODE>
        static T bilinearTriangleVel(const sxyz<T>& node,
                                     const NODE& node1,
                                     const NODE& node2,
                                     const NODE& node3) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> AC = {node3.getX()-node1.getX(), node3.getY()-node1.getY(), node3.getZ()-node1.getZ()};

            T S = norm(cross(AB,AC));

            sxyz<T> IB = {node2.getX()-node.x, node2.getY()-node.y, node2.getZ()-node.z};
            sxyz<T> IC = {node3.getX()-node.x, node3.getY()-node.y, node3.getZ()-node.z};

            T w1 = norm(cross(IB,IC))/S;
            T w2 = norm(cross(AC,IC))/S;
            T w3 = norm(cross(IB,AB))/S;

            return (1.0/(w1/node1.getNodeSlowness() + w2/node2.getNodeSlowness() + w3/node3.getNodeSlowness()));
        }

        /**
         * \brief Areal interpolation of slowness in a triangle in 3D
         *
         * Same area-ratio weights as bilinearTriangleVel(), applied directly
         * to the slownesses of the vertices, \f$s = \sum_i w_i s_i\f$.
         *
         * \pre the three nodes are not collinear, and \p node lies in the
         * plane of the triangle and inside it; neither is checked.
         *
         * \tparam NODE node type, providing `getX()`, `getY()`, `getZ()` and
         *              `getNodeSlowness()`
         * \param node  point at which the slowness is evaluated
         * \param node1 first vertex of the triangle
         * \param node2 second vertex of the triangle
         * \param node3 third vertex of the triangle
         * \return the interpolated slowness
         * \sa bilinearTriangleVel() to interpolate velocity instead
         */
        template<typename NODE>
        static T bilinearTriangle(const sxyz<T>& node,
                                  const NODE& node1,
                                  const NODE& node2,
                                  const NODE& node3) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> AC = {node3.getX()-node1.getX(), node3.getY()-node1.getY(), node3.getZ()-node1.getZ()};

            T S = norm(cross(AB,AC));

            sxyz<T> IB = {node2.getX()-node.x, node2.getY()-node.y, node2.getZ()-node.z};
            sxyz<T> IC = {node3.getX()-node.x, node3.getY()-node.y, node3.getZ()-node.z};

            T w1 = norm(cross(IB,IC))/S;
            T w2 = norm(cross(AC,IC))/S;
            T w3 = norm(cross(IB,AB))/S;

            return w1*node1.getNodeSlowness() + w2*node2.getNodeSlowness() + w3*node3.getNodeSlowness();
        }

        /**
         * \brief Areal interpolation of velocity in a triangle in 3D, node form
         *
         * Identical to the sxyz overload, the evaluation point being given as
         * a node rather than as a coordinate struct.
         *
         * \pre the three vertices are not collinear, \p node lies in the plane
         * of the triangle and inside it, and the vertex slownesses are
         * non-zero; none of this is checked.
         *
         * \tparam NODE node type, providing `getX()`, `getY()`, `getZ()` and
         *              `getNodeSlowness()`
         * \param node  node at which the slowness is evaluated
         * \param node1 first vertex of the triangle
         * \param node2 second vertex of the triangle
         * \param node3 third vertex of the triangle
         * \return the slowness corresponding to the interpolated velocity
         */
        template<typename NODE>
        static T bilinearTriangleVel(const NODE& node,
                                     const NODE& node1,
                                     const NODE& node2,
                                     const NODE& node3) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> AC = {node3.getX()-node1.getX(), node3.getY()-node1.getY(), node3.getZ()-node1.getZ()};

            T S = norm(cross(AB,AC));

            sxyz<T> IB = {node2.getX()-node.getX(), node2.getY()-node.getY(), node2.getZ()-node.getZ()};
            sxyz<T> IC = {node3.getX()-node.getX(), node3.getY()-node.getY(), node3.getZ()-node.getZ()};

            T w1 = norm(cross(IB,IC))/S;
            T w2 = norm(cross(AC,IC))/S;
            T w3 = norm(cross(IB,AB))/S;

            return (1.0/(w1/node1.getNodeSlowness() + w2/node2.getNodeSlowness() + w3/node3.getNodeSlowness()));
        }

        /**
         * \brief Areal interpolation of velocity in a triangle in 3D, vector form
         *
         * Identical to the three-vertex overloads, the vertices being taken
         * from the first three entries of \p inodes.  The evaluation point and
         * the vertices may be of different node types.
         *
         * \pre \p inodes holds at least 3 nodes, they are not collinear,
         * \p node lies in the plane of the triangle and inside it, and the
         * vertex slownesses are non-zero; none of this is checked.
         *
         * \tparam NODE1 type of the evaluation node, providing `getX()`,
         *               `getY()` and `getZ()`
         * \tparam NODE2 type of the vertices, additionally providing
         *               `getNodeSlowness()`
         * \param node   node at which the slowness is evaluated
         * \param inodes vertices of the triangle
         * \return the slowness corresponding to the interpolated velocity
         */
        template<typename NODE1, typename NODE2>
        static T bilinearTriangleVel(const NODE1& node, const std::vector<NODE2*> &inodes) {
            sxyz<T> AB = {inodes[1]->getX()-inodes[0]->getX(), inodes[1]->getY()-inodes[0]->getY(), inodes[1]->getZ()-inodes[0]->getZ()};
            sxyz<T> AC = {inodes[2]->getX()-inodes[0]->getX(), inodes[2]->getY()-inodes[0]->getY(), inodes[2]->getZ()-inodes[0]->getZ()};

            T S = norm(cross(AB,AC));

            sxyz<T> IB = {inodes[1]->getX()-node.getX(), inodes[1]->getY()-node.getY(), inodes[1]->getZ()-node.getZ()};
            sxyz<T> IC = {inodes[2]->getX()-node.getX(), inodes[2]->getY()-node.getY(), inodes[2]->getZ()-node.getZ()};

            T w1 = norm(cross(IB,IC))/S;
            T w2 = norm(cross(AC,IC))/S;
            T w3 = norm(cross(IB,AB))/S;

            return (1.0/(w1/inodes[0]->getNodeSlowness() + w2/inodes[1]->getNodeSlowness() + w3/inodes[2]->getNodeSlowness()));
        }

        /**
         * \brief Areal interpolation of slowness in a triangle in 3D, vector form
         *
         * Same area-ratio weights as the bilinearTriangleVel() vector form,
         * applied directly to the slownesses of the vertices.
         *
         * \pre \p inodes holds at least 3 nodes, they are not collinear, and
         * \p node lies in the plane of the triangle and inside it; none of
         * this is checked.
         *
         * \tparam NODE1 type of the evaluation node, providing `getX()`,
         *               `getY()` and `getZ()`
         * \tparam NODE2 type of the vertices, additionally providing
         *               `getNodeSlowness()`
         * \param node   node at which the slowness is evaluated
         * \param inodes vertices of the triangle
         * \return the interpolated slowness
         */
        template<typename NODE1, typename NODE2>
        static T bilinearTriangle(const NODE1& node, const std::vector<NODE2*> &inodes) {
            sxyz<T> AB = {inodes[1]->getX()-inodes[0]->getX(), inodes[1]->getY()-inodes[0]->getY(), inodes[1]->getZ()-inodes[0]->getZ()};
            sxyz<T> AC = {inodes[2]->getX()-inodes[0]->getX(), inodes[2]->getY()-inodes[0]->getY(), inodes[2]->getZ()-inodes[0]->getZ()};

            T S = norm(cross(AB,AC));

            sxyz<T> IB = {inodes[1]->getX()-node.getX(), inodes[1]->getY()-node.getY(), inodes[1]->getZ()-node.getZ()};
            sxyz<T> IC = {inodes[2]->getX()-node.getX(), inodes[2]->getY()-node.getY(), inodes[2]->getZ()-node.getZ()};

            T w1 = norm(cross(IB,IC))/S;
            T w2 = norm(cross(AC,IC))/S;
            T w3 = norm(cross(IB,AB))/S;

            return w1*inodes[0]->getNodeSlowness() + w2*inodes[1]->getNodeSlowness() + w3*inodes[2]->getNodeSlowness();
        }

        /**
         * \brief Areal (barycentric) weights of a point in a triangle in 3D
         *
         * Computes the weights used by bilinearTriangle() and
         * bilinearTriangleVel() without applying them, so that the caller can
         * combine any node quantity.  `w[i]` is the weight of vertex `i+1`.
         *
         * \pre the three vertices are not collinear, and \p node lies in the
         * plane of the triangle and inside it; neither is checked.
         *
         * \tparam NODE node type, providing `getX()`, `getY()` and `getZ()`
         * \param[in]  node  point at which the weights are evaluated
         * \param[in]  node1 first vertex of the triangle
         * \param[in]  node2 second vertex of the triangle
         * \param[in]  node3 third vertex of the triangle
         * \param[out] w     the three weights, summing to one
         */
        template<typename NODE>
        static void bilinearTriangleWeight(const sxyz<T>& node, const NODE* node1,
                                           const NODE* node2, const NODE* node3,
                                           std::array<T,3>& w) {
            sxyz<T> AB = {node2->getX()-node1->getX(), node2->getY()-node1->getY(), node2->getZ()-node1->getZ()};
            sxyz<T> AC = {node3->getX()-node1->getX(), node3->getY()-node1->getY(), node3->getZ()-node1->getZ()};

            T S = norm(cross(AB,AC));

            sxyz<T> IB = {node2->getX()-node.x, node2->getY()-node.y, node2->getZ()-node.z};
            sxyz<T> IC = {node3->getX()-node.x, node3->getY()-node.y, node3->getZ()-node.z};

            w[0] = norm(cross(IB,IC))/S;
            w[1] = norm(cross(AC,IC))/S;
            w[2] = norm(cross(IB,AB))/S;
        }

        /**
         * \brief Areal interpolation of traveltime in a triangle in 3D
         *
         * Same area-ratio weights as bilinearTriangle(), applied to the
         * traveltimes stored in the vertices rather than to their slownesses.
         *
         * \pre the three vertices are not collinear, and \p node lies in the
         * plane of the triangle and inside it; neither is checked.
         *
         * \tparam NODE node type, providing `getX()`, `getY()`, `getZ()` and
         *              `getTT()`
         * \param node  point at which the traveltime is evaluated
         * \param node1 first vertex of the triangle
         * \param node2 second vertex of the triangle
         * \param node3 third vertex of the triangle
         * \param nt    thread number, used to select the traveltime field
         * \return the interpolated traveltime
         */
        template<typename NODE>
        static T bilinearTime(const sxyz<T> &node, const NODE &node1,
                              const NODE &node2, const NODE &node3, const size_t nt) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> AC = {node3.getX()-node1.getX(), node3.getY()-node1.getY(), node3.getZ()-node1.getZ()};

            T S = norm(cross(AB,AC));

            sxyz<T> IB = {node2.getX()-node.x, node2.getY()-node.y, node2.getZ()-node.z};
            sxyz<T> IC = {node3.getX()-node.x, node3.getY()-node.y, node3.getZ()-node.z};

            T w1 = norm(cross(IB,IC))/S;
            T w2 = norm(cross(AC,IC))/S;
            T w3 = norm(cross(IB,AB))/S;

            return w1*node1.getTT(nt) + w2*node2.getTT(nt) + w3*node3.getTT(nt);
        }


        /**
         * \brief Volumetric (barycentric) weights of a point in a tetrahedron
         *
         * Computes the weights used by the trilinearTriangle() family without
         * applying them.  The weight of a vertex is the volume of the
         * sub-tetrahedron opposite to it divided by the volume of the whole
         * tetrahedron, both obtained from determinants.  `w[i]` is the weight
         * of vertex `i+1`.
         *
         * \pre the four vertices are not coplanar (the determinant, six times
         * the volume, would vanish) and \p node lies inside the tetrahedron;
         * neither is checked.
         *
         * \tparam NODE node type, providing `getX()`, `getY()` and `getZ()`
         * \param[in]  node  point at which the weights are evaluated
         * \param[in]  node1 first vertex of the tetrahedron
         * \param[in]  node2 second vertex of the tetrahedron
         * \param[in]  node3 third vertex of the tetrahedron
         * \param[in]  node4 fourth vertex of the tetrahedron
         * \param[out] w     the four weights, summing to one
         */
        template<typename NODE>
        static void trilinearTriangleWeight(const sxyz<T>& node, const NODE* node1,
                                            const NODE* node2, const NODE* node3,
                                            const NODE* node4, std::array<T,4>& w) {
            sxyz<T> AB = {node2->getX()-node1->getX(), node2->getY()-node1->getY(), node2->getZ()-node1->getZ()};
            sxyz<T> AC = {node3->getX()-node1->getX(), node3->getY()-node1->getY(), node3->getZ()-node1->getZ()};
            sxyz<T> AD = {node4->getX()-node1->getX(), node4->getY()-node1->getY(), node4->getZ()-node1->getZ()};

            T V = std::abs(det(AB, AC, AD));

            sxyz<T> IA = {node1->getX()-node.x, node1->getY()-node.y, node1->getZ()-node.z};
            sxyz<T> IB = {node2->getX()-node.x, node2->getY()-node.y, node2->getZ()-node.z};
            sxyz<T> IC = {node3->getX()-node.x, node3->getY()-node.y, node3->getZ()-node.z};
            sxyz<T> ID = {node4->getX()-node.x, node4->getY()-node.y, node4->getZ()-node.z};

            w[0] = std::abs(det(IB, IC, ID))/V;
            w[1] = std::abs(det(IC, IA, ID))/V;
            w[2] = std::abs(det(IB, IA, ID))/V;
            w[3] = std::abs(det(IB, IA, IC))/V;
        }

        /**
         * \brief Volumetric interpolation of slowness in a tetrahedron
         *
         * Applies the volume-ratio weights of trilinearTriangleWeight()
         * directly to the slownesses of the vertices,
         * \f$s = \sum_i w_i s_i\f$.
         *
         * \pre the four vertices are not coplanar and \p node lies inside the
         * tetrahedron; neither is checked.
         *
         * \tparam NODE node type, providing `getX()`, `getY()`, `getZ()` and
         *              `getNodeSlowness()`
         * \param node  point at which the slowness is evaluated
         * \param node1 first vertex of the tetrahedron
         * \param node2 second vertex of the tetrahedron
         * \param node3 third vertex of the tetrahedron
         * \param node4 fourth vertex of the tetrahedron
         * \return the interpolated slowness
         * \sa trilinearTriangleVel() to interpolate velocity instead
         */
        template<typename NODE>
        static T trilinearTriangle(const sxyz<T>& node, const NODE& node1,
                                   const NODE& node2, const NODE& node3,
                                   const NODE& node4) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> AC = {node3.getX()-node1.getX(), node3.getY()-node1.getY(), node3.getZ()-node1.getZ()};
            sxyz<T> AD = {node4.getX()-node1.getX(), node4.getY()-node1.getY(), node4.getZ()-node1.getZ()};

            T V = std::abs(det(AB, AC, AD));

            sxyz<T> IA = {node1.getX()-node.x, node1.getY()-node.y, node1.getZ()-node.z};
            sxyz<T> IB = {node2.getX()-node.x, node2.getY()-node.y, node2.getZ()-node.z};
            sxyz<T> IC = {node3.getX()-node.x, node3.getY()-node.y, node3.getZ()-node.z};
            sxyz<T> ID = {node4.getX()-node.x, node4.getY()-node.y, node4.getZ()-node.z};

            T w1 = std::abs(det(IB, IC, ID))/V;
            T w2 = std::abs(det(IC, IA, ID))/V;
            T w3 = std::abs(det(IB, IA, ID))/V;
            T w4 = std::abs(det(IB, IA, IC))/V;

            return (w1*node1.getNodeSlowness() + w2*node2.getNodeSlowness() + w3*node3.getNodeSlowness() + w4*node4.getNodeSlowness());
        }

        /**
         * \brief Volumetric interpolation of velocity in a tetrahedron
         *
         * Uses the same volume-ratio weights as trilinearTriangle(), but
         * interpolates the velocities of the vertices and returns the slowness
         * of the result, \f$s = \left( \sum_i w_i / s_i \right)^{-1}\f$.
         *
         * \pre the four vertices are not coplanar, \p node lies inside the
         * tetrahedron, and the vertex slownesses are non-zero; none of this is
         * checked.
         *
         * \tparam NODE node type, providing `getX()`, `getY()`, `getZ()` and
         *              `getNodeSlowness()`
         * \param node  point at which the slowness is evaluated
         * \param node1 first vertex of the tetrahedron
         * \param node2 second vertex of the tetrahedron
         * \param node3 third vertex of the tetrahedron
         * \param node4 fourth vertex of the tetrahedron
         * \return the slowness corresponding to the interpolated velocity
         * \sa trilinearTriangle() to interpolate slowness instead
         */
        template<typename NODE>
        static T trilinearTriangleVel(const sxyz<T>& node, const NODE& node1,
                                      const NODE& node2, const NODE& node3,
                                      const NODE& node4) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> AC = {node3.getX()-node1.getX(), node3.getY()-node1.getY(), node3.getZ()-node1.getZ()};
            sxyz<T> AD = {node4.getX()-node1.getX(), node4.getY()-node1.getY(), node4.getZ()-node1.getZ()};

            T V = std::abs(det(AB, AC, AD));

            sxyz<T> IA = {node1.getX()-node.x, node1.getY()-node.y, node1.getZ()-node.z};
            sxyz<T> IB = {node2.getX()-node.x, node2.getY()-node.y, node2.getZ()-node.z};
            sxyz<T> IC = {node3.getX()-node.x, node3.getY()-node.y, node3.getZ()-node.z};
            sxyz<T> ID = {node4.getX()-node.x, node4.getY()-node.y, node4.getZ()-node.z};

            T w1 = std::abs(det(IB, IC, ID))/V;
            T w2 = std::abs(det(IC, IA, ID))/V;
            T w3 = std::abs(det(IB, IA, ID))/V;
            T w4 = std::abs(det(IB, IA, IC))/V;

            return 1.0/(w1/node1.getNodeSlowness() + w2/node2.getNodeSlowness() + w3/node3.getNodeSlowness() + w4/node4.getNodeSlowness());
        }

        /**
         * \brief Volumetric interpolation of slowness in a tetrahedron, node form
         *
         * Identical to the sxyz overload, the evaluation point being given as
         * a node rather than as a coordinate struct.
         *
         * \pre the four vertices are not coplanar and \p node lies inside the
         * tetrahedron; neither is checked.
         *
         * \tparam NODE node type, providing `getX()`, `getY()`, `getZ()` and
         *              `getNodeSlowness()`
         * \param node  node at which the slowness is evaluated
         * \param node1 first vertex of the tetrahedron
         * \param node2 second vertex of the tetrahedron
         * \param node3 third vertex of the tetrahedron
         * \param node4 fourth vertex of the tetrahedron
         * \return the interpolated slowness
         */
        template<typename NODE>
        static T trilinearTriangle(const NODE& node, const NODE& node1,
                                   const NODE& node2, const NODE& node3,
                                   const NODE& node4) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> AC = {node3.getX()-node1.getX(), node3.getY()-node1.getY(), node3.getZ()-node1.getZ()};
            sxyz<T> AD = {node4.getX()-node1.getX(), node4.getY()-node1.getY(), node4.getZ()-node1.getZ()};

            T V = std::abs(det(AB, AC, AD));

            sxyz<T> IA = {node1.getX()-node.getX(), node1.getY()-node.getY(), node1.getZ()-node.getZ()};
            sxyz<T> IB = {node2.getX()-node.getX(), node2.getY()-node.getY(), node2.getZ()-node.getZ()};
            sxyz<T> IC = {node3.getX()-node.getX(), node3.getY()-node.getY(), node3.getZ()-node.getZ()};
            sxyz<T> ID = {node4.getX()-node.getX(), node4.getY()-node.getY(), node4.getZ()-node.getZ()};

            T w1 = std::abs(det(IB, IC, ID))/V;
            T w2 = std::abs(det(IC, IA, ID))/V;
            T w3 = std::abs(det(IB, IA, ID))/V;
            T w4 = std::abs(det(IB, IA, IC))/V;

            return (w1*node1.getNodeSlowness() + w2*node2.getNodeSlowness() + w3*node3.getNodeSlowness() + w4*node4.getNodeSlowness());
        }

        /**
         * \brief Volumetric interpolation of slowness in a tetrahedron, vector form
         *
         * Identical to the four-vertex overloads, the vertices being taken
         * from \p nodes.
         *
         * \pre the four vertices are not coplanar and \p node lies inside the
         * tetrahedron; neither is checked.  The size of \p nodes is checked
         * with assert(), hence only when NDEBUG is not defined.
         *
         * \tparam NODE node type, providing `getX()`, `getY()`, `getZ()` and
         *              `getNodeSlowness()`
         * \param node  point at which the slowness is evaluated
         * \param nodes the 4 vertices of the tetrahedron, passed by value
         * \return the interpolated slowness
         */
        template<typename NODE>
        static T trilinearTriangle(const sxyz<T>& node,
                                   const std::vector<NODE*> nodes) {

            assert(nodes.size()==4);
            sxyz<T> AB = {nodes[1]->getX()-nodes[0]->getX(), nodes[1]->getY()-nodes[0]->getY(), nodes[1]->getZ()-nodes[0]->getZ()};
            sxyz<T> AC = {nodes[2]->getX()-nodes[0]->getX(), nodes[2]->getY()-nodes[0]->getY(), nodes[2]->getZ()-nodes[0]->getZ()};
            sxyz<T> AD = {nodes[3]->getX()-nodes[0]->getX(), nodes[3]->getY()-nodes[0]->getY(), nodes[3]->getZ()-nodes[0]->getZ()};

            T V = std::abs(det(AB, AC, AD));

            sxyz<T> IA = {nodes[0]->getX()-node.x, nodes[0]->getY()-node.y, nodes[0]->getZ()-node.z};
            sxyz<T> IB = {nodes[1]->getX()-node.x, nodes[1]->getY()-node.y, nodes[1]->getZ()-node.z};
            sxyz<T> IC = {nodes[2]->getX()-node.x, nodes[2]->getY()-node.y, nodes[2]->getZ()-node.z};
            sxyz<T> ID = {nodes[3]->getX()-node.x, nodes[3]->getY()-node.y, nodes[3]->getZ()-node.z};

            T w1 = std::abs(det(IB, IC, ID))/V;
            T w2 = std::abs(det(IC, IA, ID))/V;
            T w3 = std::abs(det(IB, IA, ID))/V;
            T w4 = std::abs(det(IB, IA, IC))/V;

            return (w1*nodes[0]->getNodeSlowness() + w2*nodes[1]->getNodeSlowness() + w3*nodes[2]->getNodeSlowness() + w4*nodes[3]->getNodeSlowness());
        }

        /**
         * \brief Volumetric interpolation of velocity in a tetrahedron, vector form
         *
         * Uses the same volume-ratio weights as the trilinearTriangle() vector
         * form, but interpolates the velocities of the vertices and returns
         * the slowness of the result.
         *
         * \pre the four vertices are not coplanar, \p node lies inside the
         * tetrahedron, and the vertex slownesses are non-zero; none of this is
         * checked.  The size of \p nodes is checked with assert(), hence only
         * when NDEBUG is not defined.
         *
         * \tparam NODE node type, providing `getX()`, `getY()`, `getZ()` and
         *              `getNodeSlowness()`
         * \param node  point at which the slowness is evaluated
         * \param nodes the 4 vertices of the tetrahedron, passed by value
         * \return the slowness corresponding to the interpolated velocity
         */
        template<typename NODE>
        static T trilinearTriangleVel(const sxyz<T>& node,
                                      const std::vector<NODE*> nodes) {

            assert(nodes.size()==4);
            sxyz<T> AB = {nodes[1]->getX()-nodes[0]->getX(), nodes[1]->getY()-nodes[0]->getY(), nodes[1]->getZ()-nodes[0]->getZ()};
            sxyz<T> AC = {nodes[2]->getX()-nodes[0]->getX(), nodes[2]->getY()-nodes[0]->getY(), nodes[2]->getZ()-nodes[0]->getZ()};
            sxyz<T> AD = {nodes[3]->getX()-nodes[0]->getX(), nodes[3]->getY()-nodes[0]->getY(), nodes[3]->getZ()-nodes[0]->getZ()};

            T V = std::abs(det(AB, AC, AD));

            sxyz<T> IA = {nodes[0]->getX()-node.x, nodes[0]->getY()-node.y, nodes[0]->getZ()-node.z};
            sxyz<T> IB = {nodes[1]->getX()-node.x, nodes[1]->getY()-node.y, nodes[1]->getZ()-node.z};
            sxyz<T> IC = {nodes[2]->getX()-node.x, nodes[2]->getY()-node.y, nodes[2]->getZ()-node.z};
            sxyz<T> ID = {nodes[3]->getX()-node.x, nodes[3]->getY()-node.y, nodes[3]->getZ()-node.z};

            T w1 = std::abs(det(IB, IC, ID))/V;
            T w2 = std::abs(det(IC, IA, ID))/V;
            T w3 = std::abs(det(IB, IA, ID))/V;
            T w4 = std::abs(det(IB, IA, IC))/V;

            return 1.0/(w1/nodes[0]->getNodeSlowness() + w2/nodes[1]->getNodeSlowness() + w3/nodes[2]->getNodeSlowness() + w4/nodes[3]->getNodeSlowness());
        }

        /**
         * \brief Volumetric interpolation of velocity in a tetrahedron, node form
         *
         * Identical to the sxyz overload, the evaluation point being given as
         * a node rather than as a coordinate struct.
         *
         * \pre the four vertices are not coplanar, \p node lies inside the
         * tetrahedron, and the vertex slownesses are non-zero; none of this is
         * checked.
         *
         * \tparam NODE node type, providing `getX()`, `getY()`, `getZ()` and
         *              `getNodeSlowness()`
         * \param node  node at which the slowness is evaluated
         * \param node1 first vertex of the tetrahedron
         * \param node2 second vertex of the tetrahedron
         * \param node3 third vertex of the tetrahedron
         * \param node4 fourth vertex of the tetrahedron
         * \return the slowness corresponding to the interpolated velocity
         */
        template<typename NODE>
        static T trilinearTriangleVel(const NODE& node, const NODE& node1,
                                      const NODE& node2, const NODE& node3,
                                      const NODE& node4) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> AC = {node3.getX()-node1.getX(), node3.getY()-node1.getY(), node3.getZ()-node1.getZ()};
            sxyz<T> AD = {node4.getX()-node1.getX(), node4.getY()-node1.getY(), node4.getZ()-node1.getZ()};

            T V = std::abs(det(AB, AC, AD));

            sxyz<T> IA = {node1.getX()-node.getX(), node1.getY()-node.getY(), node1.getZ()-node.getZ()};
            sxyz<T> IB = {node2.getX()-node.getX(), node2.getY()-node.getY(), node2.getZ()-node.getZ()};
            sxyz<T> IC = {node3.getX()-node.getX(), node3.getY()-node.getY(), node3.getZ()-node.getZ()};
            sxyz<T> ID = {node4.getX()-node.getX(), node4.getY()-node.getY(), node4.getZ()-node.getZ()};

            T w1 = std::abs(det(IB, IC, ID))/V;
            T w2 = std::abs(det(IC, IA, ID))/V;
            T w3 = std::abs(det(IB, IA, ID))/V;
            T w4 = std::abs(det(IB, IA, IC))/V;

            return 1.0/(w1/node1.getNodeSlowness() + w2/node2.getNodeSlowness() + w3/node3.getNodeSlowness() + w4/node4.getNodeSlowness());
        }

        /**
         * \brief Volumetric interpolation of traveltime in a tetrahedron
         *
         * Same volume-ratio weights as trilinearTriangle(), applied to the
         * traveltimes stored in the vertices rather than to their slownesses.
         *
         * \pre the four vertices are not coplanar and \p node lies inside the
         * tetrahedron; neither is checked.
         *
         * \tparam NODE node type, providing `getX()`, `getY()`, `getZ()` and
         *              `getTT()`
         * \param node  point at which the traveltime is evaluated
         * \param node1 first vertex of the tetrahedron
         * \param node2 second vertex of the tetrahedron
         * \param node3 third vertex of the tetrahedron
         * \param node4 fourth vertex of the tetrahedron
         * \param nt    thread number, used to select the traveltime field
         * \return the interpolated traveltime
         */
        template<typename NODE>
        static T trilinearTime(const sxyz<T> &node, const NODE &node1,
                               const NODE &node2, const NODE &node3,
                               const NODE &node4, const size_t nt) {
            sxyz<T> AB = {node2.getX()-node1.getX(), node2.getY()-node1.getY(), node2.getZ()-node1.getZ()};
            sxyz<T> AC = {node3.getX()-node1.getX(), node3.getY()-node1.getY(), node3.getZ()-node1.getZ()};
            sxyz<T> AD = {node4.getX()-node1.getX(), node4.getY()-node1.getY(), node4.getZ()-node1.getZ()};

            T V = std::abs(det(AB,AC,AD));

            sxyz<T> IA = {node1.getX()-node.x, node1.getY()-node.y, node1.getZ()-node.z};
            sxyz<T> IB = {node2.getX()-node.x, node2.getY()-node.y, node2.getZ()-node.z};
            sxyz<T> IC = {node3.getX()-node.x, node3.getY()-node.y, node3.getZ()-node.z};
            sxyz<T> ID = {node4.getX()-node.x, node4.getY()-node.y, node4.getZ()-node.z};

            T w1 = std::abs(det(IB,IC,ID))/V;
            T w2 = std::abs(det(IC,IA,ID))/V;
            T w3 = std::abs(det(IB,IA,ID))/V;
            T w4 = std::abs(det(IB,IA,IC))/V;

            return w1*node1.getTT(nt) + w2*node2.getTT(nt) + w3*node3.getTT(nt) + w4*node4.getTT(nt);
        }
    };

}

#endif

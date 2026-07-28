//
//  Metric.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-02-01.
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
 * @file Metric.h
 * @brief Distance metrics used to build the node sweeping orders of the
 *        unstructured fast sweeping methods.
 *
 * Declares @ref ttcr::Metric, an abstract functor returning the distance
 * between a node and a point, together with its two concrete implementations
 * @ref ttcr::Metric1 (@f$\ell^1@f$, "taxicab") and @ref ttcr::Metric2
 * (@f$\ell^2@f$, Euclidean).
 *
 * The fast sweeping solvers (`Grid2Dunfs`, `Grid2Ducfs`, `Grid3Dunfs`,
 * `Grid3Ducfs`) cannot sweep along grid axes as the rectilinear variants do,
 * because an unstructured mesh has no natural index ordering.  Their
 * `initOrdering()` instead builds one node ordering per user-supplied reference
 * point, sorting all mesh nodes by their distance to that reference point; the
 * sweeps then run up and down each of those orderings in turn.  The metric
 * chosen here is what defines "distance" for that sort, and is selected at run
 * time by the `order` argument of `initOrdering()` (1 → @ref ttcr::Metric1,
 * anything else → @ref ttcr::Metric2), which in the command-line driver comes
 * from the `basename`/`Params` field `par.order`.
 *
 * @note The metric only fixes the *order* in which nodes are relaxed.  It does
 *       not enter the local traveltime solver, so it affects the convergence
 *       rate of the sweeps, not the fixed point they converge to.
 *
 * @sa Node.h, Grid2Dunfs.h, Grid2Ducfs.h, Grid3Dunfs.h, Grid3Ducfs.h
 */

#ifndef ttcr_Metric_h
#define ttcr_Metric_h

#include <cmath>

#include "Node.h"
#include "ttcr_t.h"

namespace ttcr {

    /**
     * @brief Abstract distance metric between a node and a point.
     *
     * Pure interface with one overload per spatial dimension, so that the same
     * polymorphic pointer serves both the 2D and 3D fast sweeping grids.  A
     * concrete metric is normally created on the heap, used to sort the nodes,
     * and deleted once the orderings are built.
     *
     * @tparam T floating-point type for coordinates and distances.
     */
    template <typename T>
    class Metric {
    public:
        /**
         * @brief Distance between a node and a point in the x-z plane.
         * @param[in] n node whose coordinates are read via Node::getX and
         *              Node::getZ.
         * @param[in] s reference point.
         * @return Distance between @p n and @p s in the implemented metric.
         */
        virtual T l(const Node<T>& n, const sxz<T>& s) const = 0;
        /**
         * @brief Distance between a node and a point in 3D.
         * @param[in] n node whose coordinates are read via Node::getX,
         *              Node::getY and Node::getZ.
         * @param[in] s reference point.
         * @return Distance between @p n and @p s in the implemented metric.
         */
        virtual T l(const Node<T>& n, const sxyz<T>& s) const = 0;
        /// Virtual destructor, so that derived metrics may be deleted through a
        /// @ref Metric pointer.
        virtual ~Metric() {}
    };


    /**
     * @brief Taxicab (Manhattan) @f$\ell^1@f$ metric.
     *
     * Sums the absolute coordinate differences,
     * @f[ l = \sum_i \left| x_i^{node} - x_i^{pt} \right| @f]
     *
     * Cheaper than @ref Metric2 as it needs no square root, and its level sets
     * are diamonds rather than circles, which yields sweeping orders aligned
     * with the coordinate directions.
     *
     * @tparam T floating-point type for coordinates and distances.
     */
    template <typename T>
    class Metric1 : public Metric<T> {
    public:
        /// @copydoc Metric::l(const Node<T>&, const sxz<T>&) const
        T l(const Node<T>& n, const sxz<T>& s) const {
            return std::abs(n.getX()-s.x) + std::abs(n.getZ()-s.z);
        }
        /// @copydoc Metric::l(const Node<T>&, const sxyz<T>&) const
        T l(const Node<T>& n, const sxyz<T>& s) const {
            return std::abs(n.getX()-s.x) + std::abs(n.getY()-s.y) + std::abs(n.getZ()-s.z);
        }
        ~Metric1() {}
    };


    /**
     * @brief Euclidean @f$\ell^2@f$ metric.
     *
     * Straight-line distance,
     * @f[ l = \sqrt{ \sum_i \left( x_i^{node} - x_i^{pt} \right)^2 } @f]
     *
     * Its level sets are circles (2D) or spheres (3D), so nodes are relaxed in
     * order of increasing radial distance from the reference point.  This is
     * the metric used whenever `order != 1`.
     *
     * @tparam T floating-point type for coordinates and distances.
     */
    template <typename T>
    class Metric2 : public Metric<T> {
    public:
        /// @copydoc Metric::l(const Node<T>&, const sxz<T>&) const
        T l(const Node<T>& n, const sxz<T>& s) const {
            return std::sqrt( (n.getX()-s.x)*(n.getX()-s.x) +
                             (n.getZ()-s.z)*(n.getZ()-s.z) );
        }
        /// @copydoc Metric::l(const Node<T>&, const sxyz<T>&) const
        T l(const Node<T>& n, const sxyz<T>& s) const {
            return std::sqrt( (n.getX()-s.x)*(n.getX()-s.x) +
                             (n.getY()-s.y)*(n.getY()-s.y)  +
                             (n.getZ()-s.z)*(n.getZ()-s.z) );
        }
        ~Metric2() {}
    };

}

#endif

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

#ifndef ttcr_Metric_h
#define ttcr_Metric_h

#include <cmath>

#include "Node.h"

namespace ttcr {

    template <typename T>
    class Metric {
    public:
        virtual T l(const Node<T>&, const sxz<T>&) const = 0;
        virtual T l(const Node<T>&, const sxyz<T>&) const = 0;
        virtual ~Metric() {}
    };


    template <typename T>
    class Metric1 : public Metric<T> {
    public:
        T l(const Node<T>& n, const sxz<T>& s) const {
            return std::abs(n.getX()-s.x) + std::abs(n.getZ()-s.z);
        }
        T l(const Node<T>& n, const sxyz<T>& s) const {
            return std::abs(n.getX()-s.x) + std::abs(n.getY()-s.y) + std::abs(n.getZ()-s.z);
        }
        ~Metric1() {}
    };


    template <typename T>
    class Metric2 : public Metric<T> {
    public:
        T l(const Node<T>& n, const sxz<T>& s) const {
            return sqrt( (n.getX()-s.x)*(n.getX()-s.x) +
                        (n.getZ()-s.z)*(n.getZ()-s.z) );
        }
        T l(const Node<T>& n, const sxyz<T>& s) const {
            return sqrt( (n.getX()-s.x)*(n.getX()-s.x) +
                        (n.getY()-s.y)*(n.getY()-s.y)  +
                        (n.getZ()-s.z)*(n.getZ()-s.z) );
        }
        ~Metric2() {}
    };

}

#endif

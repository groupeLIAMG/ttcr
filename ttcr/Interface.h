//
//  Interface.h
//  ttcr
//
//  Created by Bernard Giroux on 2013-01-19.
//  Copyright (c) 2013 Bernard Giroux. All rights reserved.
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
 * @file Interface.h
 * @brief Container for a reflecting interface: its points, their traveltimes
 *        and the raypaths reaching them.
 *
 * Declares ttcr::Interface, which groups the discretised points of a reflector
 * with the traveltime computed at each one and, optionally, the raypath that
 * reached it. It plays the same role for a reflector that Rcv does for a
 * receiver array: a set of positions plus the results attached to them.
 *
 * @warning **Currently unused.** No other file in the project references
 *          ttcr::Interface; reflected arrivals are presently handled inside the
 *          grid classes rather than through this container. It is documented
 *          here as-is, but treat it as unexercised code — nothing compiles it,
 *          so it has never been instantiated (see the note on
 *          ttcr::Interface::addPoint).
 *
 * @sa Rcv.h, Src.h
 */

#ifndef ttcr_Interface_h
#define ttcr_Interface_h

#include <string>
#include <vector>

#include "ttcr_t.h"

namespace ttcr {

    /**
     * @brief Points of a reflecting interface with their traveltimes and raypaths.
     *
     * @tparam P point type used for coordinates, e.g. @ref sxz or @ref sxyz.
     * @tparam T floating-point type of the traveltimes.
     *
     * @warning Unused throughout the project — see the file-level warning.
     */
    template<typename P, class T>
    class Interface {
    public:
        /// Construct an empty interface, to be filled with @ref addPoint.
        Interface() {}
        /**
         * @brief Construct from existing points and traveltimes.
         * @param[in] c interface point coordinates.
         * @param[in] t traveltime at each point; expected to match @p c in size.
         */
        Interface(const std::vector<P> &c, const std::vector<T> &t) :
        coord(c), tt(t) {}

        /// @return Coordinates of the interface points.
        const std::vector<P>& get_coord() const { return coord; }

        /// @return Traveltime at each interface point.
        const std::vector<T>& get_tt() const { return tt; }
        /// @return Writable traveltimes, for a solver to fill in.
        std::vector<T>& get_tt() { return tt; }

        /// @return Raypath reaching each interface point, as a list of points.
        const std::vector<std::vector<P>>& get_r_data() const { return r_data; }
        /// @return Writable raypaths, for a solver to fill in.
        std::vector<std::vector<P>>& get_r_data() { return r_data; }

        /**
         * @brief Append a point to the interface, with traveltime initialised to 0.
         * @param[in] pt coordinates of the new point.
         * @note Grows @ref coord and @ref tt but not @ref r_data, so after using
         *       this the three vectors are no longer parallel. Any consumer must
         *       size @ref r_data itself before writing raypaths.
         */
        void addPoint(const P& pt) {
            coord.push_back( pt );
            tt.push_back( 0 );
        }

    private:
        std::vector<P> coord;               ///< Coordinates of the interface points.
        std::vector<T> tt;                  ///< Traveltime at each point, parallel to @ref coord.
        std::vector<std::vector<P>> r_data; ///< Raypath reaching each point; not resized by @ref addPoint.
    };

}

#endif

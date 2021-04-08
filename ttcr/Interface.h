//
//  SrcRcv.h
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

#ifndef ttcr_Interface_h
#define ttcr_Interface_h

#include <string>
#include <vector>

#include "ttcr_t.h"

namespace ttcr {

    template<typename P, class T>
    class Interface {
    public:
        Interface() {}
        Interface(const std::vector<P> &c, const std::vector<T> &t) :
        coord(c), tt(t) {}

        const std::vector<P>& get_coord() const { return coord; }

        const std::vector<T>& get_tt() const { return tt; }
        std::vector<T>& get_tt() { return tt; }

        const std::vector<std::vector<P>>& get_r_data() const { return r_data; }
        std::vector<std::vector<P>>& get_r_data() { return r_data; }

        void addPoint(const P& pt) {
            coord.push_back( pt );
            tt.push_back( 0 );
        }

    private:
        std::vector<P> coord;
        std::vector<T> tt;
        std::vector<std::vector<P>> r_data;
    };

}

#endif

//
//  Grid3Dttcr.cpp
//  ttcr
//
//  Created by Bernard Giroux on 16-10-11.
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

#include <thread>

#include "Grid3Dttcr.h"

using namespace std;


namespace ttcr {
    
    Grid3Dttcr::Grid3Dttcr(const uint32_t nx, const uint32_t ny, const uint32_t nz, const double ddx,
                           const double minx, const double miny, const double minz,
                           const double eps, const int maxit, const bool w,
                           const size_t nt=1) {
        grid_instance = new grid(ny, ny, nz, ddx, minx, miny, minz, eps, maxit, w, nt);
    }
    
    Grid3Dttcr::setSlowness(const std::vector<double>& slowness) {
        if ( grid_instance->setSlowness(slowness) == 1 ) {
            throw out_of_range("Slowness values must be defined for each grid cell.");
        }
    }
    
    
}

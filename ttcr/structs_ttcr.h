//
//  structs_ttcr.h
//  ttcr_u
//
//  Created by Bernard Giroux on 2012-11-19.
//  Copyright (c) 2012 Bernard Giroux. All rights reserved.
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

#ifndef ttcr_u_structs_ttcr_h
#define ttcr_u_structs_ttcr_h

#include <map>
#include <string>

namespace ttcr {
    
    enum raytracing_method { SHORTEST_PATH, FAST_MARCHING, FAST_SWEEPING };
    
    struct input_parameters {
        uint32_t nn[3];
        int nt;
        int verbose;
        int order;                    // order of l metric
        int nitermax;
        bool inverseDistance;
        bool singlePrecision;
        bool saveRaypaths;
        bool saveModelVTK;
        bool saveGridTT;
        bool time;
        bool processReflectors;
        bool raypath_high_order;
        bool rotated_template;
        bool weno3;
        double epsilon;
        double source_radius;
        raytracing_method method;
        std::string basename;
        std::string modelfile;
        std::string velfile;
        std::string slofile;
        std::string rcvfile;
        std::vector<std::string> srcfiles;
        
        input_parameters() : nn(), nt(0), verbose(0), order(2), nitermax(20),
        inverseDistance(false),	singlePrecision(false), saveRaypaths(false),
        saveModelVTK(false), saveGridTT(false), time(false), processReflectors(false),
        raypath_high_order(false), rotated_template(false), weno3(false),
        epsilon(1.e-15), source_radius(0.0), method(SHORTEST_PATH), basename(),
        modelfile(), velfile(), slofile(), rcvfile(), srcfiles() {}
        
    };
    
}

#endif

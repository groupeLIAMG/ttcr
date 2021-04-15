//
//  structs_ttcr.h
//  ttcr
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

#ifndef ttcr_structs_ttcr_h
#define ttcr_structs_ttcr_h

#include <map>
#include <string>

namespace ttcr {

    enum raytracing_method { SHORTEST_PATH, FAST_MARCHING, FAST_SWEEPING, DYNAMIC_SHORTEST_PATH };
    enum gradient_method : int { LS_FO=0, LS_SO=1, AB=2 };

    struct input_parameters {
        uint32_t nn[3];
        int nt;
        int order;                    // order of l metric
        int nitermax;
        int nTertiary;
        int raypath_method;
        int saveGridTT;
        int min_per_thread;
        bool inverseDistance;
        bool singlePrecision;
        bool saveRaypaths;
        bool saveModelVTK;
        bool saveM;
        bool time;
        bool processReflectors;
        bool projectTxRx;
        bool processVel;
        bool rotated_template;
        bool weno3;
        bool dump_secondary;
        bool tt_from_rp;
        bool useEdgeLength;
        bool translateOrigin;
        double epsilon;
        double source_radius;
        double min_distance_rp;
        double radius_tertiary_nodes;
        raytracing_method method;
        std::string basename;
        std::string modelfile;
        std::string velfile;
        std::string slofile;
        std::string rcvfile;
        std::vector<std::string> srcfiles;

        input_parameters() : nn(), nt(0), order(2), nitermax(20),
        nTertiary(3), raypath_method(LS_SO), saveGridTT(0), min_per_thread(5),
        inverseDistance(false), singlePrecision(false), saveRaypaths(false),
        saveModelVTK(false), saveM(false), time(false), processReflectors(false),
        projectTxRx(false), processVel(false), rotated_template(false),
        weno3(false), dump_secondary(false), tt_from_rp(false),
        useEdgeLength(true), translateOrigin(false),
        epsilon(1.e-15), source_radius(0.0), min_distance_rp(1.e-5),
        radius_tertiary_nodes(0.0), method(SHORTEST_PATH), basename(),
        modelfile(), velfile(), slofile(), rcvfile(), srcfiles() {}

    };

}

#endif

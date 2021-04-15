//
//  structs_msh2vtk.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-10-18.
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

#ifndef ttcr_structs_msh2vtk_h
#define ttcr_structs_msh2vtk_h

#include <string>

struct input_parameters {
    bool rectilinear;
    bool crt;
    bool saveSlowness;
    bool saveReflectors;
    bool saveXYZ;
    int nSecondary;
    double d;
    std::string mshFile;
    std::string vtkFile;
    std::string velFile;
    std::string sloFile;

    input_parameters() : rectilinear(false), crt(false), saveSlowness(false),
    saveReflectors(false), saveXYZ(false), nSecondary(0), d(0.0),
    mshFile(""), vtkFile(""), velFile(""), sloFile("") {}
};

#endif

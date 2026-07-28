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

/**
 * @file structs_msh2vtk.h
 * @brief Run-time options of the @c msh2vtk conversion utility.
 *
 * Holds the single ::input_parameters aggregate that carries the parsed command
 * line of @c msh2vtk from msh2vtk_io.h into the conversion code.
 *
 * @warning This struct is declared in the **global namespace** and shares its
 *          name with ttcr::input_parameters (structs_ttcr.h), which is an
 *          unrelated type with different members.  The two are never included
 *          in the same translation unit; keep it that way, as doing so would
 *          make every unqualified use of @c input_parameters ambiguous.
 *
 * @sa msh2vtk_io.h, structs_ttcr.h
 */

#ifndef ttcr_structs_msh2vtk_h
#define ttcr_structs_msh2vtk_h

#include <string>

/**
 * @brief Options controlling a Gmsh-to-VTK mesh conversion.
 *
 * Default-constructed to a "convert the mesh, write nothing extra" state: no
 * rectilinear resampling, no secondary nodes, no slowness or reflector output.
 */
struct input_parameters {
    bool rectilinear;      ///< Resample the mesh onto a rectilinear grid of spacing @ref d instead of copying it verbatim.
    bool crt;              ///< Read/write point coordinates in CRT format (one record per line, terminated by @c /).
    bool saveSlowness;     ///< Attach the slowness field to the output as a data array.
    bool saveReflectors;   ///< Write the physical entities tagged as reflectors as separate output files.
    bool saveXYZ;          ///< Also dump node coordinates to a plain XYZ text file.
    int nSecondary;        ///< Number of secondary nodes to insert per edge; 0 leaves only the primary mesh vertices.
    double d;              ///< Cell size of the rectilinear grid used when @ref rectilinear is true.
    std::string mshFile;   ///< Input Gmsh @c .msh file.
    std::string vtkFile;   ///< Output VTK file.
    std::string velFile;   ///< Optional velocity file, applied to the mesh before writing.
    std::string sloFile;   ///< Optional slowness file; alternative to @ref velFile.

    /// Construct with conversion disabled everywhere and all filenames empty.
    input_parameters() : rectilinear(false), crt(false), saveSlowness(false),
    saveReflectors(false), saveXYZ(false), nSecondary(0), d(0.0),
    mshFile(""), vtkFile(""), velFile(""), sloFile("") {}
};

#endif

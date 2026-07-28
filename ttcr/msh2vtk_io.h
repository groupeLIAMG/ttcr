//
//  msh2vtk_io.h
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
 * @file msh2vtk_io.h
 * @brief Command-line front end of the @c msh2vtk conversion utility.
 *
 * Declares the two functions that turn an @c argv vector into a populated
 * ::input_parameters record for @c msh2vtk, the standalone tool that converts a
 * Gmsh @c .msh mesh into a VTK unstructured grid.  These are free functions in
 * the global namespace, not in @c ttcr, because @c msh2vtk is built as a
 * separate program rather than as part of the library.
 *
 * @note This header is not self-contained: it names @c std::ostream but pulls in
 *       only structs_msh2vtk.h.  It compiles because every current includer has
 *       already included @c \<iostream\>.
 *
 * @sa structs_msh2vtk.h, MSHReader.h
 */

#ifndef ttcr_msh2vtk_io_h
#define ttcr_msh2vtk_io_h

#include "structs_msh2vtk.h"

/**
 * @brief Print the @c msh2vtk usage message and terminate.
 * @param[in,out] stream    destination stream (@c std::cout for @c -h,
 *                          @c std::cerr for a usage error).
 * @param[in] exit_code     status handed to @c exit(); 0 for a requested help
 *                          listing, nonzero for a bad command line.
 * @warning Does not return — it calls @c exit().
 */
void print_usage (std::ostream& stream, int exit_code);

/**
 * @brief Parse the command line into an ::input_parameters record.
 * @param[in] argc      argument count as received by @c main.
 * @param[in] argv      argument vector as received by @c main.
 * @param[out] params   record filled from the options; fields left unset keep
 *                      the defaults from the ::input_parameters constructor.
 */
void parse_input(int argc, char * argv[], input_parameters &params);

#endif /* defined(__ttcr__msh2vtk_io__) */

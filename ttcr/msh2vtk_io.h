//
//  msh2vtk_io.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-10-18.
//  Copyright (c) 2014 Bernard Giroux. All rights reserved.
//

#ifndef __ttcr__msh2vtk_io__
#define __ttcr__msh2vtk_io__

#include "structs_msh2vtk.h"

void print_usage (std::ostream&, int);
void parse_input(int argc, char * argv[], input_parameters &);

#endif /* defined(__ttcr__msh2vtk_io__) */

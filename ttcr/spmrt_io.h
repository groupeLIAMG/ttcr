//
//  spmrt_io.h
//  ttcr_u
//
//  Created by Bernard Giroux on 2012-11-19.
//  Copyright (c) 2012 Bernard Giroux. All rights reserved.
//

#ifndef __ttcr_u__spmrt_io__
#define __ttcr_u__spmrt_io__

#include <iostream>
#include <string>
#include <vector>

#include "ttcr_t.h"
#include "structs_spmrt.h"

void print_usage (std::ostream&, char *, int);
std::string parse_input(int argc, char * argv[], input_parameters &);
void get_params(const std::string &, input_parameters &);

#endif /* defined(__ttcr_u__spmrt_io__) */

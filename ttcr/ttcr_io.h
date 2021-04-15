//
//  ttcr_io.h
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

#ifndef ttcr_ttcr_io_h
#define ttcr_ttcr_io_h

#include <iostream>
#include <string>
#include <vector>

#include "ttcr_t.h"
#include "structs_ttcr.h"

namespace ttcr {

    void print_usage (std::ostream&, char *, int);
    std::string parse_input(int argc, char * argv[], input_parameters &);
    void get_params(const std::string &, input_parameters &);

    class AtomicWriter {
        std::ostringstream st;
        std::ostream &stream;
    public:
        AtomicWriter(std::ostream &s=std::cout):stream(s) { }
        template <typename T>
        AtomicWriter& operator<<(T const& t) {
            st << t;
            return *this;
        }
        AtomicWriter& operator<<( std::ostream&(*f)(std::ostream&) ) {
            st << f;
            return *this;
        }
        ~AtomicWriter() { stream << st.str(); }
    };

}

#endif /* defined(__ttcr__spmrt_io__) */

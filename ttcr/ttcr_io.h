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

/**
 * @file ttcr_io.h
 * @brief Command-line and parameter-file front end of the @c ttcr programs.
 *
 * Declares the three functions that populate a ttcr::input_parameters record —
 * @ref ttcr::print_usage, @ref ttcr::parse_input for the command-line switches
 * and @ref ttcr::get_params for the parameter file — plus
 * @ref ttcr::AtomicWriter, a small helper that keeps output from interleaving
 * when several solver threads write at once.
 *
 * The two parsers are complementary and run in that order: @c parse_input
 * handles the switches and returns the parameter file name, which is then fed
 * to @c get_params.
 *
 * @sa structs_ttcr.h, grids.h
 */

#ifndef ttcr_ttcr_io_h
#define ttcr_ttcr_io_h

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "ttcr_t.h"
#include "structs_ttcr.h"

namespace ttcr {

    /**
     * @brief Print the usage message and terminate the program.
     * @param[in,out] stream    destination (@c std::cout for @c -h,
     *                          @c std::cerr for a bad option).
     * @param[in] progname      program name, normally @c argv[0].
     * @param[in] exit_code     status handed to @c exit(): 0 for a requested
     *                          help listing, 1 for a command-line error.
     * @warning Does not return — it calls @c exit().
     */
    void print_usage (std::ostream& stream, char *progname, int exit_code);

    /**
     * @brief Parse the command-line switches.
     *
     * Recognises @c -h (usage), @c -k (@ref input_parameters::saveModelVTK),
     * @c -v (increase the global verbosity counter), @c -t
     * (@ref input_parameters::time), @c -s
     * (@ref input_parameters::dump_secondary) and the mandatory @c -p, which
     * names the parameter file.
     *
     * @param[in] argc      argument count as received by @c main.
     * @param[in] argv      argument vector as received by @c main.
     * @param[in,out] params  record whose flag fields are set from the switches.
     * @return Name of the parameter file given with @c -p, to be passed to
     *         @ref get_params.
     * @warning If no option is supplied, or @c -p is missing, this prints the
     *          usage message and exits rather than returning.
     */
    std::string parse_input(int argc, char * argv[], input_parameters &params);

    /**
     * @brief Read the parameter file into @p params.
     *
     * The file is a list of @c "value # keyword, comment" records: the value
     * comes first, then a @c \# separator, then a keyword matched by substring.
     * Unrecognised keywords are skipped silently, and any field the file does
     * not mention keeps its constructor default — so a typo in a keyword
     * degrades to "option ignored" rather than to an error.
     *
     * Each documented field of ttcr::input_parameters names its keyword.
     * @c "srcfile" is the one key that may appear repeatedly, appending to
     * ttcr::input_parameters::srcfiles.
     *
     * @param[in] filename  parameter file, normally the return of
     *                      @ref parse_input.
     * @param[in,out] params  record to fill.
     * @warning Exits the program if @p filename cannot be opened.
     */
    void get_params(const std::string &filename, input_parameters &params);

    /**
     * @brief Buffers a sequence of stream insertions and emits them in one write.
     *
     * Accumulates everything inserted into an internal @c std::ostringstream and
     * flushes it to the target stream in the destructor, as a single @c operator<<
     * call. Because the solvers run several threads that all report progress,
     * writing directly to @c std::cout would interleave characters from different
     * threads mid-line; routing a whole message through a temporary
     * @c AtomicWriter keeps each message contiguous.
     *
     * Typical use is a temporary, so the flush happens at the end of the statement:
     * @code
     * AtomicWriter() << "thread " << n << " done\n";
     * @endcode
     *
     * @note This makes each message contiguous, not the stream thread-safe: the
     *       final write is one insertion, but it is not otherwise synchronised.
     */
    class AtomicWriter {
        std::ostringstream st;   ///< Buffer holding the message until destruction.
        std::ostream &stream;    ///< Destination the buffer is flushed to.
    public:
        /// @param s stream to flush to when this object is destroyed.
        AtomicWriter(std::ostream &s=std::cout):stream(s) { }
        /**
         * @brief Buffer any streamable value.
         * @tparam T type of the value; anything with an @c operator<< overload.
         * @param[in] t value to append to the buffer.
         * @return Reference to this writer, so insertions chain.
         */
        template <typename T>
        AtomicWriter& operator<<(T const& t) {
            st << t;
            return *this;
        }
        /**
         * @brief Buffer a stream manipulator such as @c std::endl.
         * @param[in] f manipulator to apply to the buffer.
         * @return Reference to this writer, so insertions chain.
         */
        AtomicWriter& operator<<( std::ostream&(*f)(std::ostream&) ) {
            st << f;
            return *this;
        }
        /// Flush the accumulated message to the target stream in a single write.
        ~AtomicWriter() { stream << st.str(); }
    };

}

#endif /* defined(__ttcr__spmrt_io__) */

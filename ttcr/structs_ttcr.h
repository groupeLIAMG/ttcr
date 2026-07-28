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

/**
 * @file structs_ttcr.h
 * @brief Run-time configuration of the @c ttcr command-line programs.
 *
 * Defines the solver and gradient-method enumerations together with
 * ttcr::input_parameters, the single record that carries every option from the
 * command line and the parameter file through to grid construction in grids.h.
 *
 * The record is populated by ttcr_io.h (@c parse_input for the switches,
 * @c get_params for the parameter file) and then read by the grid factories,
 * which dispatch on ttcr::input_parameters::method to pick a solver class.
 *
 * @sa ttcr_io.h, grids.h, structs_msh2vtk.h
 */

#ifndef ttcr_structs_ttcr_h
#define ttcr_structs_ttcr_h

#include <cstdint>
#include <map>
#include <string>

namespace ttcr {

    /**
     * @brief Eikonal solver selected for the traveltime computation.
     *
     * Chooses which family of grid classes the factories in grids.h instantiate;
     * the concrete class also depends on the grid type and on whether slowness
     * is stored per cell or per node.
     */
    enum raytracing_method {
        SHORTEST_PATH,          ///< Shortest-path method on a graph of primary and secondary nodes (@c *sp classes).
        FAST_MARCHING,          ///< Fast marching, single-pass narrow band ordered by a priority queue (@c *fm classes).
        FAST_SWEEPING,          ///< Fast sweeping, Gauss-Seidel iterations along alternating orderings (@c *fs classes).
        DYNAMIC_SHORTEST_PATH,  ///< Shortest path with secondary nodes inserted adaptively near the source (@c *dsp classes).
        FAST_SWEEPING_OPENCL    ///< Fast sweeping offloaded to an OpenCL device (@c *fs_OpenCL classes).
    };

    /**
     * @brief Scheme used to estimate the traveltime gradient when tracing rays.
     *
     * The raypath is integrated by stepping down the traveltime gradient, so
     * this selects how that gradient is recovered from the nodal traveltimes.
     * Stored in ttcr::input_parameters::raypath_method.
     *
     * @sa Grad.h
     */
    enum gradient_method : int {
        LS_FO=0,  ///< Least-squares fit, first order.
        LS_SO=1,  ///< Least-squares fit, second order (the default).
        AB=2      ///< Averaged-gradient (Aki-Backus style) estimate.
    };

    /**
     * @brief Every run-time option of a @c ttcr run, in one record.
     *
     * Default-constructed to a usable configuration — second-order least-squares
     * raypaths, the shortest-path solver, double precision, 50 sweep iterations —
     * so a caller only needs to override what it cares about plus the filenames.
     *
     * @note Not all fields apply to all solvers. @ref nitermax, @ref epsilon and
     *       @ref weno3 are read only by the fast sweeping classes;
     *       @ref nTertiary and @ref radius_tertiary_nodes only by the dynamic
     *       shortest-path classes; @ref gpu_max_threads and @ref profile only by
     *       the OpenCL ones.
     */
    struct input_parameters {
        uint32_t nn[3];               ///< Secondary nodes per cell edge, per axis (@c "secondary nodes"). A single value in the file is copied to all three.
        int nt;                       ///< Number of worker threads; 0 lets the program choose (@c "number of threads").
        int order;                    ///< Order of the @f$\ell@f$ metric ordering the fast sweeping sweeps: 1 or 2 (@c "metric order"). @sa Metric.h
        int nitermax;                 ///< Iteration cap for the fast sweeping solvers (@c "max number of iteration").
        int nTertiary;                ///< Tertiary (dynamic) nodes added per edge near the source by the @c dsp solvers (@c "tertiary nodes").
        int raypath_method;           ///< Gradient scheme used to trace rays; values of ::gradient_method (@c "gradient method").
        int saveGridTT;               ///< If > 0, dump the traveltime field over the whole grid to @c \<basename\>_all_tt; the value is passed on to @c saveTT as its format selector.
        int min_per_thread;           ///< Minimum sources per thread before another thread is spawned (@c "min nb Tx per thread").
        int gpu_max_threads;          ///< Cap on concurrent GPU solver streams (OpenCL FSM); 0 = no cap (@c "gpu max threads").
        bool inverseDistance;         ///< Use inverse-distance rather than linear weighting when interpolating onto Tx/Rx positions (@c "inverse distance").
        bool singlePrecision;         ///< Build the grid with @c float instead of @c double (@c "single precision").
        bool saveRaypaths;            ///< Write the traced raypaths alongside the traveltimes (@c "saveRayPaths").
        bool saveModelVTK;            ///< Dump the velocity/slowness model as VTK (command-line @c -k).
        bool saveM;                   ///< Also compute and save the sensitivity (Fréchet derivative) matrix; selects the @c raytrace overload taking @c m_data (@c "save M").
        bool time;                    ///< Report wall-clock time for grid build and raytracing (command-line @c -t).
        bool processReflectors;       ///< Trace reflected arrivals off the interfaces named in the model (@c "process reflectors").
        bool projectTxRx;             ///< Project sources and receivers onto the nearest grid entity instead of rejecting off-grid points (@c "project Tx Rx").
        bool processVel;              ///< Interpolate velocity rather than slowness, returning @f$1/\sum_i w_i/s_i@f$ (@c "interpolate velocity"). @sa Interpolator.h
        bool rotated_template;        ///< Use rotated stencils in the fast sweeping update (@c "rotated template").
        bool weno3;                   ///< Use the 3rd-order WENO stencil in the fast sweeping solvers (@c "fsm high order").
        bool dump_secondary;          ///< Write the generated secondary nodes to an ASCII file (command-line @c -s).
        bool tt_from_rp;              ///< Recompute the receiver traveltime by integrating slowness along the traced raypath, rather than reading the interpolated nodal value (@c "traveltime from raypath").
        bool useEdgeLength;           ///< Interpret @ref radius_tertiary_nodes as a multiple of the mean cell edge length rather than as an absolute distance.
        bool translateOrigin;         ///< Shift the grid so its origin sits at (0,0,0) before solving, improving conditioning for grids far from the origin (@c "translate grid origin").
        bool profile;                 ///< Emit a GPU profiling breakdown (OpenCL) (@c "profile").
        double epsilon;               ///< Fast sweeping convergence tolerance: the **mean** per-node @f$|\Delta t|@f$ between successive sweeps (@c "epsilon").
        double source_radius;         ///< If nonzero, treat the source as a sphere of this radius rather than a point (@c "source radius").
        double min_distance_rp;       ///< Minimum step retained when integrating a raypath; guards against stalled steps (@c "raypath minimum distance").
        double radius_tertiary_nodes; ///< Radius around the source within which tertiary nodes are inserted; scaled by edge length when @ref useEdgeLength is set (@c "radius dynamic nodes").
        raytracing_method method;     ///< Eikonal solver to use. @sa ::raytracing_method
        std::string basename;         ///< Prefix for every output file (@c "basename").
        std::string modelfile;        ///< Grid/mesh file defining the model geometry (@c "modelfile").
        std::string velfile;          ///< Velocity file; mutually exclusive with @ref slofile (@c "velfile").
        std::string slofile;          ///< Slowness file; mutually exclusive with @ref velfile (@c "slofile").
        std::string rcvfile;          ///< Receiver coordinate file (@c "rcvfile"). @sa Rcv, Rcv2D
        std::vector<std::string> srcfiles;  ///< Source files; the @c "srcfile" key may appear repeatedly and each occurrence appends. @sa Src, Src2D

        /// Construct with the documented defaults: shortest path, double
        /// precision, second-order least-squares raypaths, 50 sweep iterations.
        input_parameters() : nn(), nt(0), order(2), nitermax(50),
        nTertiary(3), raypath_method(LS_SO), saveGridTT(0), min_per_thread(5),
        gpu_max_threads(4),
        inverseDistance(false), singlePrecision(false), saveRaypaths(false),
        saveModelVTK(false), saveM(false), time(false), processReflectors(false),
        projectTxRx(false), processVel(false), rotated_template(false),
        weno3(false), dump_secondary(false), tt_from_rp(false),
        useEdgeLength(true), translateOrigin(false), profile(false),
        epsilon(1.e-5), source_radius(0.0), min_distance_rp(1.e-5),
        radius_tertiary_nodes(0.0), method(SHORTEST_PATH), basename(),
        modelfile(), velfile(), slofile(), rcvfile(), srcfiles() {}

    };

}

#endif

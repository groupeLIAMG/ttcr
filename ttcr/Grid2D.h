//
//  Grid2D.h
//  ttcr
//
//  Created by Bernard Giroux on 2014-01-21.
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
 * @file Grid2D.h
 * @brief Abstract base class for 2D traveltime grids and meshes.
 *
 * Defines @ref ttcr::Grid2D, the common interface implemented by every 2D
 * rectilinear-grid and unstructured-mesh raytracing class in ttcr.  It provides
 * the public entry points for single-source and multi-source (threaded)
 * traveltime computation and raypath extraction, dispatching to protected
 * virtual hooks that concrete subclasses override.  The base class itself only
 * implements the buffer sizing / receiver iteration / thread scheduling logic;
 * the actual eikonal solve and raypath tracing live in the derived classes.
 */

#ifndef ttcr_Grid2D_h
#define ttcr_Grid2D_h

#include <algorithm>
#include <array>
#include <fstream>
#include <functional>
#include <future>
#include <thread>
#include <vector>

#include "ctpl_stl.h"

#include "ttcr_t.h"

namespace ttcr {

    /**
     * @brief Abstract base for 2D traveltime grids and raytracing.
     *
     * Concrete subclasses model a 2D velocity/slowness medium (rectilinear grid
     * or triangular mesh) and solve the eikonal equation to obtain first-arrival
     * traveltimes and raypaths.  This base class defines the public API shared by
     * all of them and coordinates receiver iteration and multi-threaded execution
     * over source points, delegating the numerical work to virtual hooks.
     *
     * @tparam T1 Floating-point type for coordinates, slowness and traveltimes
     *            (typically `double` or `float`).
     * @tparam T2 Integer type used for node and cell indices (typically
     *            `uint32_t`).
     * @tparam S  Point type used for source/receiver/raypath coordinates
     *            (e.g. `sxz<double>` for a vertical 2D plane).
     */
    template<typename T1 = double, typename T2 = uint32_t, typename S = sxz<double>>
    class Grid2D {
    public:
        /**
         * @brief Constructor.
         *
         * @param ncells Number of cells in the grid/mesh; sizes the per-cell
         *               neighbour lookup table.
         * @param ttrp   If true, receiver traveltimes are recomputed by
         *               integrating slowness along the extracted raypath rather
         *               than read from the traveltime field.
         * @param nt     Number of worker threads to support.
         * @param up     If true (and `nt > 1`), a persistent thread pool is
         *               created and reused across raytrace calls; otherwise
         *               transient `std::thread`s are spawned per call.
         */
        Grid2D(const size_t ncells, const bool ttrp, const size_t nt=1, const bool up=1) :
        nThreads(nt), tt_from_rp(ttrp), neighbors(std::vector<std::vector<T2>>(ncells)),
        usePool(up)
        {
            if ( nThreads > 1 && usePool ) {
                pool.resize(static_cast<int>(nThreads));
            }
        }

        /// Virtual destructor.
        virtual ~Grid2D() {}

        /// @return Number of worker threads this grid was configured with.
        const size_t getNthreads() const { return nThreads; }

        /**
         * @name Thread-pool call operators
         *
         * Entry points invoked by the thread pool (`ctpl::thread_pool`).  Each
         * simply forwards to the matching single-source @ref raytrace overload,
         * passing the pool-assigned worker index as the thread number.  The set
         * of trailing output containers selects which per-receiver data are
         * produced: `r_data` = raypath points, `l_data` = ray-length
         * sensitivities (`siv`/`siv2`), `m_data` = model-parameter sensitivities
         * (`sijv`).
         * @{
         */
        // operators: for usage with thread pool
        void operator()(int id, const std::vector<S>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<S>& Rx,
                        std::vector<T1>& traveltimes) const {
            this->raytrace(Tx, t0, Rx, traveltimes, id);
        }

        void operator()(int id, const std::vector<S>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<S>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<S>>& r_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, r_data, id);
        }

        void operator()(int id, const std::vector<S>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<S>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<siv2<T1>>>& l_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, l_data, id);
        }

        void operator()(int id, const std::vector<S>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<S>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<siv<T1>>>& l_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, l_data, id);
        }

        void operator()(int id, const std::vector<S>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<S>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<siv4<T1>>>& l_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, l_data, id);
        }

        void operator()(int id, const std::vector<S>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<S>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<siv5<T1>>>& l_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, l_data, id);
        }

        void operator()(int id, const std::vector<S>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<S>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<sijv<T1>>>& m_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, m_data, id);
        }

        void operator()(int id, const std::vector<S>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<S>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<S>>& r_data,
                        std::vector<std::vector<siv2<T1>>>& l_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, r_data, l_data, id);
        }

        void operator()(int id, const std::vector<S>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<S>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<S>>& r_data,
                        std::vector<std::vector<siv<T1>>>& l_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, r_data, l_data, id);
        }

        void operator()(int id, const std::vector<S>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<S>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<S>>& r_data,
                        std::vector<std::vector<siv4<T1>>>& l_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, r_data, l_data, id);
        }

        void operator()(int id, const std::vector<S>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<S>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<S>>& r_data,
                        std::vector<std::vector<siv5<T1>>>& l_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, r_data, l_data, id);
        }

        void operator()(int id, const std::vector<S>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<S>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<S>>& r_data,
                        std::vector<std::vector<sijv<T1>>>& m_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, r_data, m_data, id);
        }
        /// @}

        /**
         * @brief Compute first-arrival traveltimes from a set of sources to a
         *        set of receivers.
         *
         * @param Tx          Source (transmitter) positions.
         * @param t0          Initial time at each source in @p Tx.
         * @param Rx          Receiver positions.
         * @param[out] traveltimes  Computed traveltime at each receiver.
         * @param threadNo    Worker-thread index used for scratch storage.
         */
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              std::vector<T1>& traveltimes,
                              const size_t threadNo=0) const;

        /**
         * @brief Traveltimes to receivers grouped per source-gather.
         *
         * Variant taking receivers as a vector of pointers to per-source
         * receiver lists, filling the corresponding per-source traveltime lists.
         */

        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<const std::vector<S>*>& Rx,
                              std::vector<std::vector<T1>*>& traveltimes,
                              const size_t threadNo=0) const;

        /// @brief As above, also returning the raypath (`r_data`) to each receiver.
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<S>>& r_data,
                              const size_t threadNo=0) const;

        /// @brief Per-source-gather traveltimes and raypaths (pointer-grouped receivers).
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<const std::vector<S>*>& Rx,
                              std::vector<std::vector<T1>*>& traveltimes,
                              std::vector<std::vector<std::vector<S>>*>& r_data,
                              const size_t threadNo=0) const;

        /// @brief Traveltimes, raypaths and per-cell ray-length sensitivities (`siv`).
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<S>>& r_data,
                              std::vector<std::vector<siv<T1>>>& l_data,
                              const size_t threadNo=0) const;

        /// @brief Traveltimes, raypaths and ray-length sensitivities (`siv2`). Subclass must override.
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<S>>& r_data,
                              std::vector<std::vector<siv2<T1>>>& l_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method raytrace should be implemented in subclass");
        }

        /// @brief Traveltimes, raypaths and four-component sensitivities (`siv4`). Subclass must override.
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<S>>& r_data,
                              std::vector<std::vector<siv4<T1>>>& l_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method raytrace should be implemented in subclass");
        }

        /// @brief Traveltimes, raypaths and five-component sensitivities (`siv5`). Subclass must override.
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<S>>& r_data,
                              std::vector<std::vector<siv5<T1>>>& l_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method raytrace should be implemented in subclass");
        }

        /// @brief Traveltimes, raypaths and model-parameter sensitivities (`sijv`). Subclass must override.
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<S>>& r_data,
                              std::vector<std::vector<sijv<T1>>>& m_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method raytrace should be implemented in subclass");
        }

        /// @brief Traveltimes and ray-length sensitivities only (`siv2`, no raypath). Subclass must override.
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<siv2<T1>>>& l_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method raytrace should be implemented in subclass");
        }

        /// @brief Traveltimes and four-component sensitivities only (`siv4`, no raypath). Subclass must override.
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<siv4<T1>>>& l_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method raytrace should be implemented in subclass");
        }

        /// @brief Traveltimes and five-component sensitivities only (`siv5`, no raypath). Subclass must override.
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<siv5<T1>>>& l_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method raytrace should be implemented in subclass");
        }

        /// @brief Traveltimes and model-parameter sensitivities only (`sijv`, no raypath). Subclass must override.
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<sijv<T1>>>& m_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method raytrace should be implemented in subclass");
        }

        /// @brief Traveltimes and per-cell ray-length sensitivities only (`siv`, no raypath).
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<siv<T1>>>& l_data,
                              const size_t threadNo=0) const;

        /// @brief Traveltimes and raypaths, also returning the source velocity @p v0. Subclass must override.
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<S>>& r_data,
                              T1& v0,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method raytrace should be implemented in subclass");
        }

        /// @brief Traveltimes, raypaths, source velocity @p v0 and model sensitivities (`sijv`). Subclass must override.
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<S>>& r_data,
                              T1& v0,
                              std::vector<std::vector<sijv<T1>>>& m_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method raytrace should be implemented in subclass");
        }

        /**
         * @name Multi-source threaded raytracing
         *
         * Overloads taking one source list, initial-time list and receiver list
         * per shot.  Each shot is dispatched to a worker thread (via the thread
         * pool when enabled, otherwise via a block partition over transient
         * threads); a single shot or a single-threaded grid runs inline.  The
         * trailing output containers select which per-receiver data are produced
         * (`r_data`, `l_data`, `m_data`), matching the single-source overloads.
         * @{
         */
        // methods for threaded raytracing
        void raytrace(const std::vector<std::vector<S>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<S>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes) const;

        void raytrace(const std::vector<std::vector<S>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<S>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<S>>>& r_data) const;

        void raytrace(const std::vector<std::vector<S>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<S>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<siv2<T1>>>>& l_data) const;

        void raytrace(const std::vector<std::vector<S>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<S>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<siv<T1>>>>& l_data) const;

        /// @brief Several sources, four-component sensitivities (@ref siv4).
        void raytrace(const std::vector<std::vector<S>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<S>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<siv4<T1>>>>& l_data) const {
            raytraceMulti(Tx, t0, Rx, traveltimes, l_data);
        }

        /// @brief Several sources, five-component sensitivities (@ref siv5).
        void raytrace(const std::vector<std::vector<S>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<S>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<siv5<T1>>>>& l_data) const {
            raytraceMulti(Tx, t0, Rx, traveltimes, l_data);
        }

        void raytrace(const std::vector<std::vector<S>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<S>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<S>>>& r_data,
                      std::vector<std::vector<std::vector<siv2<T1>>>>& l_data) const;

        void raytrace(const std::vector<std::vector<S>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<S>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<S>>>& r_data,
                      std::vector<std::vector<std::vector<sijv<T1>>>>& m_data) const;

        void raytrace(const std::vector<std::vector<S>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<S>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<S>>>& r_data,
                      std::vector<std::vector<std::vector<siv<T1>>>>& l_data) const;

        /// @brief Several sources, raypaths and four-component sensitivities.
        void raytrace(const std::vector<std::vector<S>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<S>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<S>>>& r_data,
                      std::vector<std::vector<std::vector<siv4<T1>>>>& l_data) const {
            raytraceMulti(Tx, t0, Rx, traveltimes, r_data, l_data);
        }

        /// @brief Several sources, raypaths and five-component sensitivities.
        void raytrace(const std::vector<std::vector<S>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<S>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<S>>>& r_data,
                      std::vector<std::vector<std::vector<siv5<T1>>>>& l_data) const {
            raytraceMulti(Tx, t0, Rx, traveltimes, r_data, l_data);
        }

        void raytrace(const std::vector<std::vector<S>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<S>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<sijv<T1>>>>& m_data) const;
        /// @}

        /// @brief Set the per-cell (or per-node) slowness model.
        virtual void setSlowness(const std::vector<T1>& s) {}
        /// @brief Set the slowness model from a raw array of @p ns values. Subclass must override.
        virtual void setSlowness(const T1 *s, const size_t ns) {
            throw std::runtime_error("Method setSlowness should be implemented in subclass");
        }
        /// @brief Set the anisotropy parameter &xi; (weak elliptical anisotropy). Subclass must override.
        virtual void setXi(const std::vector<T1>& x) {
            throw std::runtime_error("Method setXi should be implemented in subclass");
        }
        /// @brief Set the symmetry-axis tilt angle of the anisotropy. Subclass must override.
        virtual void setTiltAngle(const std::vector<T1>& x) {
            throw std::runtime_error("Method setTiltAngle should be implemented in subclass");
        }
        /// @brief Set the P-wave axial velocity Vp0. Subclass must override.
        virtual void setVp0(const std::vector<T1>& s) {
            throw std::runtime_error("Method setVp0 should be implemented in subclass");
        }
        /// @brief Set the S-wave axial velocity Vs0. Subclass must override.
        virtual void setVs0(const std::vector<T1>& s) {
            throw std::runtime_error("Method setVs0 should be implemented in subclass");
        }
        /// @brief Set the Thomsen anisotropy parameter &delta;. Subclass must override.
        virtual void setDelta(const std::vector<T1>& s) {
            throw std::runtime_error("Method setDelta should be implemented in subclass");
        }
        /// @brief Set the Thomsen anisotropy parameter &epsilon;. Subclass must override.
        virtual void setEpsilon(const std::vector<T1>& s) {
            throw std::runtime_error("Method setEpsilon should be implemented in subclass");
        }
        /// @brief Set the Thomsen anisotropy parameter &gamma;. Subclass must override.
        virtual void setGamma(const std::vector<T1>& s) {
            throw std::runtime_error("Method setGamma should be implemented in subclass");
        }
        /// @brief Set the 2nd-order slowness coefficient (`s2`). Subclass must override.
        virtual void setS2(const std::vector<T1>& s) {
            throw std::runtime_error("Method setS2 should be implemented in subclass");
        }
        /// @brief Set the 4th-order slowness coefficient (`s4`). Subclass must override.
        virtual void setS4(const std::vector<T1>& s) {
            throw std::runtime_error("Method setS4 should be implemented in subclass");
        }

        /// @brief Enable/disable recomputing receiver traveltimes from the raypath.
        void setTraveltimeFromRaypath(const bool ttrp) { tt_from_rp = ttrp; }

        /// @param primary If true, count only primary (corner) nodes. @return Number of nodes.
        virtual size_t getNumberOfNodes(const bool primary=false) const { return 1; }
        /// @return Number of cells in the grid/mesh.
        virtual size_t getNumberOfCells() const { return 1; }
        /// @brief Copy the traveltime field for thread @p threadNo into @p tt. Subclass must override.
        virtual void getTT(std::vector<T1>& tt, const size_t threadNo=0) const {
            throw std::runtime_error("Method getTT should be implemented in subclass");
        }
        /// @brief Copy the slowness model into the supplied vector. Subclass must override.
        virtual void getSlowness(std::vector<T1>&) const {
            throw std::runtime_error("Method getSlowness should be implemented in subclass");
        }

        /// @brief Save the traveltime field to a file (format-dependent).
        virtual void saveTT(const std::string &, const int, const size_t nt=0,
                            const int format=1) const {}

        /// @brief Save the traveltime gradient field to a file.
        virtual void saveTTgrad(const std::string &, const size_t nt=0,
                                const bool vtkFormat=0) const {}

        virtual const T1 getXmin() const { return 1; }  ///< @return Minimum x coordinate of the grid.
        virtual const T1 getXmax() const { return 1; }  ///< @return Maximum x coordinate of the grid.
        virtual const T1 getYmin() const { return 1; }  ///< @return Minimum y coordinate of the grid.
        virtual const T1 getYmax() const { return 1; }  ///< @return Maximum y coordinate of the grid.
        virtual const T1 getZmin() const { return 1; }  ///< @return Minimum z coordinate of the grid.
        virtual const T1 getZmax() const { return 1; }  ///< @return Maximum z coordinate of the grid.
        virtual const T1 getDx() const { return 1; }    ///< @return Cell size along x.
        virtual const T1 getDz() const { return 1; }    ///< @return Cell size along z.
        virtual const T2 getNcx() const { return 1; }   ///< @return Number of cells along x.
        virtual const T2 getNcz() const { return 1; }   ///< @return Number of cells along z.
        virtual const T2 getNsnx() const { return 1; }  ///< @return Number of secondary nodes per cell edge along x.
        virtual const T2 getNsnz() const { return 1; }  ///< @return Number of secondary nodes per cell edge along z.

        virtual const int get_niter() const { return 0; }   ///< @return Number of solver iterations performed.
        virtual const int get_niterw() const { return 0; }  ///< @return Number of iterations in the wave-front refinement pass.

        /// @brief Project points onto the grid (e.g. snap to the medium). @return status code.
        virtual int projectPts(std::vector<S>&) const { return 1; }

        /// @brief Interpolate the given nodal field onto secondary nodes, in place.
        virtual void interpolateAtNodes(std::vector<T1> &) const {}
        /// @brief Interpolate the given nodal field, writing the result to the second argument.
        virtual void interpolateAtNodes(const std::vector<T1> &,
                                        std::vector<T1> &) const {}

        /// @brief Dump secondary-node data to a stream (diagnostics).
        virtual void dump_secondary(std::ofstream&) const {};

        /// @brief Retrieve the grid/mesh node coordinates. Subclass must override.
        virtual void getNodes(std::vector<S>& nodes) const {
            throw std::runtime_error("Method getNodes should be implemented in subclass");
        }
        /// @brief Retrieve triangle connectivity as fixed-size index arrays. Subclass must override.
        virtual void getTriangles(std::vector<std::array<T2, 3>>&) const {
            throw std::runtime_error("Method getTriangles should be implemented in subclass");
        }
        /// @brief Retrieve triangle connectivity as index vectors (Cython-friendly). Subclass must override.
        // keep next method until array is supported by cython
        virtual void getTriangles(std::vector<std::vector<T2>>&) const {
            throw std::runtime_error("Method getTriangles should be implemented in subclass");
        }
        
        /**
         * compute slowness at given point
         *
         * @param pt point to consider
         * @returns slowness value
         */
        virtual T1 computeSlowness(const S& pt) const {
            throw std::runtime_error("Method computeSlowness should be implemented in subclass");
        }

        /**
         * @brief Enable or disable use of the persistent thread pool.
         *
         * When enabling with more than one thread, the pool is (re)sized to the
         * configured thread count if needed.
         * @param up New use-pool flag.
         */
        void setUsePool(const bool up) {
            usePool = up;
            if ( nThreads > 1 && usePool && pool.size() != nThreads ) {
                pool.resize(static_cast<int>(nThreads));
            }
        }

        /**
         * @brief Build the interpolation weights mapping model nodes/cells onto
         *        the given points.
         * @param pts        Points to interpolate at.
         * @param[out] d_data Per-point interpolation weights (`siv`).
         * @return status code. Subclass must override.
         */
        virtual int computeD(const std::vector<S>& pts,
                             std::vector<std::vector<siv<T1>>>& d_data) const {
            throw std::runtime_error("Method computeD should be implemented in subclass");
        }

        /// @return Average edge length of the grid/mesh. Subclass must override.
        virtual const T1 getAverageEdgeLength() const {
            throw std::runtime_error("Method getAverageEdgeLength should be implemented in subclass");
        }

#ifdef VTK
        /// @brief Save the model as an unstructured-grid VTK file (`.vtu`).
        virtual void saveModelVTU(const std::string &, const bool saveSlowness=true,
                                  const bool savePhysicalEntity=false) const {}
        /// @brief Save the model as a rectilinear-grid VTK file (`.vtr`).
        virtual void saveModelVTR(const std::string &, const double*,
                                  const bool saveSlowness=true) const {}
#endif

        /**
         * @brief Spread several sources over the threads, whatever the container
         *
         * The work is distributed as in the overloads taking @ref siv or
         * @ref siv2, which predate the wider containers and keep their own
         * copy of it.
         */
        template<typename SIV>
        void raytraceMulti(const std::vector<std::vector<S>>& Tx,
                           const std::vector<std::vector<T1>>& t0,
                           const std::vector<std::vector<S>>& Rx,
                           std::vector<std::vector<T1>>& traveltimes,
                           std::vector<std::vector<std::vector<SIV>>>& l_data) const {
            if ( Tx.size() == 1 ) {
                this->raytrace(Tx[0], t0[0], Rx[0], traveltimes[0], l_data[0], 0);
            } else if ( nThreads == 1 ) {
                for ( size_t n=0; n<Tx.size(); ++n ) {
                    this->raytrace(Tx[n], t0[n], Rx[n], traveltimes[n], l_data[n], 0);
                }
            } else if ( usePool ) {
                std::vector<std::future<void>> results(Tx.size());
                for ( size_t n=0; n<Tx.size(); ++n ) {
                    results[n] = pool.push(std::ref(*this),
                                           std::ref(Tx[n]),
                                           std::ref(t0[n]),
                                           std::ref(Rx[n]),
                                           std::ref(traveltimes[n]),
                                           std::ref(l_data[n]));
                }
                for ( size_t n=0; n<Tx.size(); ++n ) {
                    results[n].get();
                }
            } else {
                std::vector<size_t> blk_size = get_blk_size(Tx.size());
                std::vector<std::thread> threads(blk_size.size());
                size_t blk_start = 0;
                for ( size_t i=0; i<blk_size.size(); ++i ) {
                    size_t blk_end = blk_start + blk_size[i];
                    threads[i]=std::thread( [this,&Tx,&t0,&Rx,&traveltimes,&l_data,blk_start,blk_end,i]{
                        for ( size_t n=blk_start; n<blk_end; ++n ) {
                            this->raytrace(Tx[n], t0[n], Rx[n], traveltimes[n], l_data[n], i);
                        }
                    });
                    blk_start = blk_end;
                }
                std::for_each(threads.begin(),threads.end(), std::mem_fn(&std::thread::join));
            }
        }

        /// @brief The same, also returning the raypaths
        template<typename SIV>
        void raytraceMulti(const std::vector<std::vector<S>>& Tx,
                           const std::vector<std::vector<T1>>& t0,
                           const std::vector<std::vector<S>>& Rx,
                           std::vector<std::vector<T1>>& traveltimes,
                           std::vector<std::vector<std::vector<S>>>& r_data,
                           std::vector<std::vector<std::vector<SIV>>>& l_data) const {
            if ( Tx.size() == 1 ) {
                this->raytrace(Tx[0], t0[0], Rx[0], traveltimes[0], r_data[0], l_data[0], 0);
            } else if ( nThreads == 1 ) {
                for ( size_t n=0; n<Tx.size(); ++n ) {
                    this->raytrace(Tx[n], t0[n], Rx[n], traveltimes[n], r_data[n], l_data[n], 0);
                }
            } else if ( usePool ) {
                std::vector<std::future<void>> results(Tx.size());
                for ( size_t n=0; n<Tx.size(); ++n ) {
                    results[n] = pool.push(std::ref(*this),
                                           std::ref(Tx[n]),
                                           std::ref(t0[n]),
                                           std::ref(Rx[n]),
                                           std::ref(traveltimes[n]),
                                           std::ref(r_data[n]),
                                           std::ref(l_data[n]));
                }
                for ( size_t n=0; n<Tx.size(); ++n ) {
                    results[n].get();
                }
            } else {
                std::vector<size_t> blk_size = get_blk_size(Tx.size());
                std::vector<std::thread> threads(blk_size.size());
                size_t blk_start = 0;
                for ( size_t i=0; i<blk_size.size(); ++i ) {
                    size_t blk_end = blk_start + blk_size[i];
                    threads[i]=std::thread( [this,&Tx,&t0,&Rx,&traveltimes,&r_data,&l_data,blk_start,blk_end,i]{
                        for ( size_t n=blk_start; n<blk_end; ++n ) {
                            this->raytrace(Tx[n], t0[n], Rx[n], traveltimes[n], r_data[n], l_data[n], i);
                        }
                    });
                    blk_start = blk_end;
                }
                std::for_each(threads.begin(),threads.end(), std::mem_fn(&std::thread::join));
            }
        }


    protected:
        size_t nThreads;                         ///< Number of worker threads.
        bool tt_from_rp;                         ///< If true, receiver traveltimes are integrated along the raypath.
        std::vector<std::vector<T2>> neighbors;  ///< For each cell, the indices of nodes belonging to it.

        /**
         * @brief Populate @ref neighbors from a node collection.
         *
         * Inverts each node's owner-cell list so that every cell records the
         * nodes attached to it.
         * @tparam N    Node type exposing `getOwners()`.
         * @param nodes Nodes whose owner information is used.
         */
        template<typename N>
        void buildGridNeighbors(std::vector<N>& nodes) {
            for ( T2 n=0; n<nodes.size(); ++n ) {
                for ( size_t n2=0; n2<nodes[n].getOwners().size(); ++n2) {
                    neighbors[ nodes[n].getOwners()[n2] ].push_back(n);
                }
            }
        }

        /**
         * @brief Core eikonal solve for one source: propagate traveltimes to all
         *        nodes (receiver traveltimes are read out separately).
         *
         * Central virtual hook the public overloads delegate to. Subclass must
         * override.
         * @param Tx       Source positions.
         * @param t0       Initial times at the sources.
         * @param Rx       Receiver positions (used to bound the computation).
         * @param threadNo Worker-thread index.
         */
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method raytrace should be implemented in subclass");
        }

        /// @brief Core eikonal solve, receivers grouped per source-gather. Subclass must override.
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<const std::vector<S>*>& Rx,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method raytrace should be implemented in subclass");
        }

        /// @brief Sample the computed traveltime field at point @p pt. Subclass must override.
        virtual T1 getTraveltime(const S& pt,
                                 const size_t threadNo) const {
            throw std::runtime_error("Method getTraveltime should be implemented in subclass");
        }

        /// @brief Traveltime at @p Rx obtained by integrating slowness along its raypath. Subclass must override.
        virtual T1 getTraveltimeFromRaypath(const std::vector<S>& Tx,
                                            const std::vector<T1>& t0,
                                            const S& Rx,
                                            const size_t threadNo) const {
            throw std::runtime_error("Method getTraveltimeFromRaypath should be implemented in subclass");
        }

        /**
         * @brief Trace the raypath from receiver @p Rx back to the nearest source.
         *
         * @param Tx          Source positions.
         * @param t0          Initial times at the sources.
         * @param Rx          Receiver position.
         * @param[out] r_data Raypath as a sequence of points.
         * @param[out] tt     Traveltime along the extracted raypath.
         * @param threadNo    Worker-thread index. Subclass must override.
         */
        virtual void getRaypath(const std::vector<S>& Tx,
                                const std::vector<T1>& t0,
                                const S& Rx,
                                std::vector<S>& r_data,
                                T1 &tt,
                                const size_t threadNo) const {
            throw std::runtime_error("Method getRaypath should be implemented in subclass");
        }

        /// @brief Trace raypath, also accumulating per-cell ray-length sensitivities (`l_data`). Subclass must override.
        virtual void getRaypath(const std::vector<S>& Tx,
                                const std::vector<T1>& t0,
                                const S& Rx,
                                std::vector<S>& r_data,
                                std::vector<siv<T1>> &l_data,
                                T1 &tt,
                                const size_t threadNo) const {
            throw std::runtime_error("Method getRaypath should be implemented in subclass");
        }

        /// @brief Trace raypath, returning only ray-length sensitivities (`l_data`, no point list). Subclass must override.
        virtual void getRaypath(const std::vector<S>& Tx,
                                const std::vector<T1>& t0,
                                const S& Rx,
                                std::vector<siv<T1>> &l_data,
                                T1 &tt,
                                const size_t threadNo) const {
            throw std::runtime_error("Method getRaypath should be implemented in subclass");
        }

    private:
        bool usePool;                    ///< Whether the persistent thread pool is used for multi-source runs.
        mutable ctpl::thread_pool pool;  ///< Reusable worker-thread pool.

        /**
         * @brief Partition @p nTx sources into per-thread block sizes.
         *
         * Distributes the sources as evenly as possible round-robin across at
         * most @ref nThreads blocks.
         * @param nTx Number of sources to distribute.
         * @return Vector of block sizes summing to @p nTx.
         */
        const std::vector<size_t> get_blk_size(const size_t nTx) const {
            size_t n_blk = nThreads < nTx ? nThreads : nTx;
            std::vector<size_t> blk_size ( n_blk, 0 );
            size_t nj = nTx;
            while ( nj > 0 ) {
                for ( size_t n=0; n<n_blk; ++n ) {
                    blk_size[n] += 1;
                    nj -= 1;
                    if ( nj == 0 ) {
                        break;
                    }
                }
            }
            return blk_size;
        }
    };

    template<typename T1, typename T2, typename S>
    void Grid2D<T1,T2,S>::raytrace(const std::vector<S>& Tx,
                                   const std::vector<T1>& t0,
                                   const std::vector<S>& Rx,
                                   std::vector<T1>& traveltimes,
                                   const size_t threadNo) const {
        this->raytrace(Tx, t0, Rx, threadNo);

        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }

        if ( tt_from_rp ) {
            for (size_t n=0; n<Rx.size(); ++n) {
                traveltimes[n] = this->getTraveltimeFromRaypath(Tx, t0, Rx[n], threadNo);
            }
        } else {
            for (size_t n=0; n<Rx.size(); ++n) {
                traveltimes[n] = this->getTraveltime(Rx[n], threadNo);
            }
        }
    }

    template<typename T1, typename T2, typename S>
    void Grid2D<T1,T2,S>::raytrace(const std::vector<S>& Tx,
                                   const std::vector<T1>& t0,
                                   const std::vector<const std::vector<S>*>& Rx,
                                   std::vector<std::vector<T1>*>& traveltimes,
                                   const size_t threadNo) const {
        this->raytrace(Tx, t0, Rx, threadNo);

        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }

        if ( tt_from_rp ) {
            for (size_t nr=0; nr<Rx.size(); ++nr) {
                traveltimes[nr]->resize( Rx[nr]->size() );
                for (size_t n=0; n<Rx[nr]->size(); ++n) {
                    (*traveltimes[nr])[n] = this->getTraveltimeFromRaypath(Tx, t0, (*Rx[nr])[n], threadNo);
                }
            }
        } else {
            for (size_t nr=0; nr<Rx.size(); ++nr) {
                traveltimes[nr]->resize( Rx[nr]->size() );
                for (size_t n=0; n<Rx[nr]->size(); ++n) {
                    (*traveltimes[nr])[n] = this->getTraveltime((*Rx[nr])[n], threadNo);
                }
            }
        }
    }

    template<typename T1, typename T2, typename S>
    void Grid2D<T1,T2,S>::raytrace(const std::vector<S>& Tx,
                                   const std::vector<T1>& t0,
                                   const std::vector<S>& Rx,
                                   std::vector<T1>& traveltimes,
                                   std::vector<std::vector<S>>& r_data,
                                   const size_t threadNo) const {
        this->raytrace(Tx, t0, Rx, threadNo);

        if ( r_data.size() != Rx.size() ) {
            r_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<r_data.size(); ++ni ) {
            r_data[ni].resize( 0 );
        }
        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        for (size_t n=0; n<Rx.size(); ++n) {
            this->getRaypath(Tx, t0, Rx[n], r_data[n], traveltimes[n], threadNo);
        }
    }

    template<typename T1, typename T2, typename S>
    void Grid2D<T1,T2,S>::raytrace(const std::vector<S>& Tx,
                                   const std::vector<T1>& t0,
                                   const std::vector<const std::vector<S>*>& Rx,
                                   std::vector<std::vector<T1>*>& traveltimes,
                                   std::vector<std::vector<std::vector<S>>*>& r_data,
                                   const size_t threadNo) const {
        this->raytrace(Tx, t0, Rx, threadNo);

        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }

        for (size_t nr=0; nr<Rx.size(); ++nr) {
            r_data[nr]->resize( Rx[nr]->size() );
            for ( size_t ni=0; ni<r_data[nr]->size(); ++ni ) {
                (*r_data[nr])[ni].resize( 0 );
            }
            traveltimes[nr]->resize( Rx[nr]->size() );

            for (size_t n=0; n<Rx[nr]->size(); ++n) {
                this->getRaypath(Tx, t0, (*Rx[nr])[n], (*r_data[nr])[n],
                                 (*traveltimes[nr])[n], threadNo);
            }
        }
    }

    template<typename T1, typename T2, typename S>
    void Grid2D<T1,T2,S>::raytrace(const std::vector<S>& Tx,
                                   const std::vector<T1>& t0,
                                   const std::vector<S>& Rx,
                                   std::vector<T1>& traveltimes,
                                   std::vector<std::vector<S>>& r_data,
                                   std::vector<std::vector<siv<T1>>>& l_data,
                                   const size_t threadNo) const {
        this->raytrace(Tx, t0, Rx, threadNo);

        if ( r_data.size() != Rx.size() ) {
            r_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<r_data.size(); ++ni ) {
            r_data[ni].resize( 0 );
        }
        if ( l_data.size() != Rx.size() ) {
            l_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<l_data.size(); ++ni ) {
            l_data[ni].resize( 0 );
        }
        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        for (size_t n=0; n<Rx.size(); ++n) {
            this->getRaypath(Tx, t0, Rx[n], r_data[n], l_data[n], traveltimes[n], threadNo);
            //  must be sorted to build matrix L
            sort(l_data[n].begin(), l_data[n].end(), CompareSiv_i<T1>());
        }
    }

    template<typename T1, typename T2, typename S>
    void Grid2D<T1,T2,S>::raytrace(const std::vector<S>& Tx,
                                   const std::vector<T1>& t0,
                                   const std::vector<S>& Rx,
                                   std::vector<T1>& traveltimes,
                                   std::vector<std::vector<siv<T1>>>& l_data,
                                   const size_t threadNo) const {
        this->raytrace(Tx, t0, Rx, threadNo);

        if ( l_data.size() != Rx.size() ) {
            l_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<l_data.size(); ++ni ) {
            l_data[ni].resize( 0 );
        }
        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }
        for (size_t n=0; n<Rx.size(); ++n) {
            this->getRaypath(Tx, t0, Rx[n], l_data[n], traveltimes[n], threadNo);
            //  must be sorted to build matrix L
            sort(l_data[n].begin(), l_data[n].end(), CompareSiv_i<T1>());
        }
    }

    template<typename T1, typename T2, typename S>
    void Grid2D<T1,T2,S>::raytrace(const std::vector<std::vector<S>>& Tx,
                                   const std::vector<std::vector<T1>>& t0,
                                   const std::vector<std::vector<S>>& Rx,
                                   std::vector<std::vector<T1>>& traveltimes) const {

        if ( Tx.size() == 1 ) {
            this->raytrace(Tx[0], t0[0], Rx[0], traveltimes[0], 0);
        } else if ( nThreads == 1 ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                this->raytrace(Tx[n], t0[n], Rx[n], traveltimes[n], 0);
            }
        } else if ( usePool ) {
            std::vector<std::future<void>> results(Tx.size());
            for ( size_t n=0; n<Tx.size(); ++n ) {
                results[n] = pool.push(std::ref(*this),
                                       std::ref(Tx[n]),
                                       std::ref(t0[n]),
                                       std::ref(Rx[n]),
                                       std::ref(traveltimes[n]));
            }
            for ( size_t n=0; n<Tx.size(); ++n ) {
                results[n].get();
            }
        } else {
            std::vector<size_t> blk_size = get_blk_size(Tx.size());

            std::vector<std::thread> threads(blk_size.size());
            size_t blk_start = 0;
            for ( size_t i=0; i<blk_size.size(); ++i ) {

                size_t blk_end = blk_start + blk_size[i];
                threads[i]=std::thread( [this,&Tx,&t0,&Rx,&traveltimes,blk_start,blk_end,i]{

                    for ( size_t n=blk_start; n<blk_end; ++n ) {
                        this->raytrace(Tx[n], t0[n], Rx[n], traveltimes[n], i);
                    }
                });

                blk_start = blk_end;
            }
            std::for_each(threads.begin(),threads.end(), std::mem_fn(&std::thread::join));
        }
    }

    template<typename T1, typename T2, typename S>
    void Grid2D<T1,T2,S>::raytrace(const std::vector<std::vector<S>>& Tx,
                                   const std::vector<std::vector<T1>>& t0,
                                   const std::vector<std::vector<S>>& Rx,
                                   std::vector<std::vector<T1>>& traveltimes,
                                   std::vector<std::vector<std::vector<S>>>& r_data) const {

        if ( Tx.size() == 1 ) {
            this->raytrace(Tx[0], t0[0], Rx[0], traveltimes[0], r_data[0], 0);
        } else if ( nThreads == 1 ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                this->raytrace(Tx[n], t0[n], Rx[n], traveltimes[n], r_data[n], 0);
            }
        } else if ( usePool ) {
            std::vector<std::future<void>> results(Tx.size());
            for ( size_t n=0; n<Tx.size(); ++n ) {
                results[n] = pool.push(std::ref(*this),
                                       std::ref(Tx[n]),
                                       std::ref(t0[n]),
                                       std::ref(Rx[n]),
                                       std::ref(traveltimes[n]),
                                       std::ref(r_data[n]));
            }
            for ( size_t n=0; n<Tx.size(); ++n ) {
                results[n].get();
            }
        } else {
            std::vector<size_t> blk_size = get_blk_size(Tx.size());

            std::vector<std::thread> threads(blk_size.size());
            size_t blk_start = 0;
            for ( size_t i=0; i<blk_size.size(); ++i ) {

                size_t blk_end = blk_start + blk_size[i];
                threads[i]=std::thread( [this,&Tx,&t0,&Rx,&traveltimes,&r_data,blk_start,blk_end,i]{

                    for ( size_t n=blk_start; n<blk_end; ++n ) {
                        this->raytrace(Tx[n], t0[n], Rx[n], traveltimes[n], r_data[n], i);
                    }
                });

                blk_start = blk_end;
            }
            std::for_each(threads.begin(),threads.end(), std::mem_fn(&std::thread::join));
        }
    }

    template<typename T1, typename T2, typename S>
    void Grid2D<T1,T2,S>::raytrace(const std::vector<std::vector<S>>& Tx,
                                   const std::vector<std::vector<T1>>& t0,
                                   const std::vector<std::vector<S>>& Rx,
                                   std::vector<std::vector<T1>>& traveltimes,
                                   std::vector<std::vector<std::vector<siv2<T1>>>>& l_data) const {

        if ( Tx.size() == 1 ) {
            this->raytrace(Tx[0], t0[0], Rx[0], traveltimes[0], l_data[0], 0);
        } else if ( nThreads == 1 ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                this->raytrace(Tx[n], t0[n], Rx[n], traveltimes[n], l_data[n], 0);
            }
        } else if ( usePool ) {
            std::vector<std::future<void>> results(Tx.size());
            for ( size_t n=0; n<Tx.size(); ++n ) {
                results[n] = pool.push(std::ref(*this),
                                       std::ref(Tx[n]),
                                       std::ref(t0[n]),
                                       std::ref(Rx[n]),
                                       std::ref(traveltimes[n]),
                                       std::ref(l_data[n]));
            }
            for ( size_t n=0; n<Tx.size(); ++n ) {
                results[n].get();
            }
        } else {
            std::vector<size_t> blk_size = get_blk_size(Tx.size());

            std::vector<std::thread> threads(blk_size.size());
            size_t blk_start = 0;
            for ( size_t i=0; i<blk_size.size(); ++i ) {

                size_t blk_end = blk_start + blk_size[i];
                threads[i]=std::thread( [this,&Tx,&t0,&Rx,&traveltimes,&l_data,blk_start,blk_end,i]{

                    for ( size_t n=blk_start; n<blk_end; ++n ) {
                        this->raytrace(Tx[n], t0[n], Rx[n], traveltimes[n], l_data[n], i);
                    }
                });

                blk_start = blk_end;
            }
            std::for_each(threads.begin(),threads.end(), std::mem_fn(&std::thread::join));
        }
    }

    template<typename T1, typename T2, typename S>
    void Grid2D<T1,T2,S>::raytrace(const std::vector<std::vector<S>>& Tx,
                                   const std::vector<std::vector<T1>>& t0,
                                   const std::vector<std::vector<S>>& Rx,
                                   std::vector<std::vector<T1>>& traveltimes,
                                   std::vector<std::vector<std::vector<siv<T1>>>>& l_data) const {

        if ( Tx.size() == 1 ) {
            this->raytrace(Tx[0], t0[0], Rx[0], traveltimes[0], l_data[0], 0);
        } else if ( nThreads == 1 ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                this->raytrace(Tx[n], t0[n], Rx[n], traveltimes[n], l_data[n], 0);
            }
        } else if ( usePool ) {
            std::vector<std::future<void>> results(Tx.size());
            for ( size_t n=0; n<Tx.size(); ++n ) {
                results[n] = pool.push(std::ref(*this),
                                       std::ref(Tx[n]),
                                       std::ref(t0[n]),
                                       std::ref(Rx[n]),
                                       std::ref(traveltimes[n]),
                                       std::ref(l_data[n]));
            }
            for ( size_t n=0; n<Tx.size(); ++n ) {
                results[n].get();
            }
        } else {
            std::vector<size_t> blk_size = get_blk_size(Tx.size());

            std::vector<std::thread> threads(blk_size.size());
            size_t blk_start = 0;
            for ( size_t i=0; i<blk_size.size(); ++i ) {

                size_t blk_end = blk_start + blk_size[i];
                threads[i]=std::thread( [this,&Tx,&t0,&Rx,&traveltimes,&l_data,blk_start,blk_end,i]{

                    for ( size_t n=blk_start; n<blk_end; ++n ) {
                        this->raytrace(Tx[n], t0[n], Rx[n], traveltimes[n], l_data[n], i);
                    }
                });

                blk_start = blk_end;
            }
            std::for_each(threads.begin(),threads.end(), std::mem_fn(&std::thread::join));
        }
    }

    template<typename T1, typename T2, typename S>
    void Grid2D<T1,T2,S>::raytrace(const std::vector<std::vector<S>>& Tx,
                                   const std::vector<std::vector<T1>>& t0,
                                   const std::vector<std::vector<S>>& Rx,
                                   std::vector<std::vector<T1>>& traveltimes,
                                   std::vector<std::vector<std::vector<sijv<T1>>>>& m_data) const {

        if ( Tx.size() == 1 ) {
            this->raytrace(Tx[0], t0[0], Rx[0], traveltimes[0], m_data[0], 0);
        } else if ( nThreads == 1 ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                this->raytrace(Tx[n], t0[n], Rx[n], traveltimes[n], m_data[n], 0);
            }
        } else if ( usePool ) {
            std::vector<std::future<void>> results(Tx.size());
            for ( size_t n=0; n<Tx.size(); ++n ) {
                results[n] = pool.push(std::ref(*this),
                                       std::ref(Tx[n]),
                                       std::ref(t0[n]),
                                       std::ref(Rx[n]),
                                       std::ref(traveltimes[n]),
                                       std::ref(m_data[n]));
            }
            for ( size_t n=0; n<Tx.size(); ++n ) {
                results[n].get();
            }
        } else {
            std::vector<size_t> blk_size = get_blk_size(Tx.size());

            std::vector<std::thread> threads(blk_size.size());
            size_t blk_start = 0;
            for ( size_t i=0; i<blk_size.size(); ++i ) {

                size_t blk_end = blk_start + blk_size[i];
                threads[i]=std::thread( [this,&Tx,&t0,&Rx,&traveltimes,&m_data,blk_start,blk_end,i]{

                    for ( size_t n=blk_start; n<blk_end; ++n ) {
                        this->raytrace(Tx[n], t0[n], Rx[n], traveltimes[n], m_data[n], i);
                    }
                });

                blk_start = blk_end;
            }
            std::for_each(threads.begin(),threads.end(), std::mem_fn(&std::thread::join));
        }
    }

    template<typename T1, typename T2, typename S>
    void Grid2D<T1,T2,S>::raytrace(const std::vector<std::vector<S>>& Tx,
                                   const std::vector<std::vector<T1>>& t0,
                                   const std::vector<std::vector<S>>& Rx,
                                   std::vector<std::vector<T1>>& traveltimes,
                                   std::vector<std::vector<std::vector<S>>>& r_data,
                                   std::vector<std::vector<std::vector<siv2<T1>>>>& l_data) const {

        if ( Tx.size() == 1 ) {
            this->raytrace(Tx[0], t0[0], Rx[0], traveltimes[0], r_data[0], l_data[0], 0);
        } else if ( nThreads == 1 ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                this->raytrace(Tx[n], t0[n], Rx[n], traveltimes[n], r_data[n], l_data[n], 0);
            }
        } else if ( usePool ) {
            std::vector<std::future<void>> results(Tx.size());
            for ( size_t n=0; n<Tx.size(); ++n ) {
                results[n] = pool.push(std::ref(*this),
                                       std::ref(Tx[n]),
                                       std::ref(t0[n]),
                                       std::ref(Rx[n]),
                                       std::ref(traveltimes[n]),
                                       std::ref(r_data[n]),
                                       std::ref(l_data[n]));
            }
            for ( size_t n=0; n<Tx.size(); ++n ) {
                results[n].get();
            }
        } else {
            std::vector<size_t> blk_size = get_blk_size(Tx.size());

            std::vector<std::thread> threads(blk_size.size());
            size_t blk_start = 0;
            for ( size_t i=0; i<blk_size.size(); ++i ) {

                size_t blk_end = blk_start + blk_size[i];
                threads[i]=std::thread( [this,&Tx,&t0,&Rx,&traveltimes,&r_data,&l_data,blk_start,blk_end,i]{

                    for ( size_t n=blk_start; n<blk_end; ++n ) {
                        this->raytrace(Tx[n], t0[n], Rx[n], traveltimes[n], r_data[n], l_data[n], i);
                    }
                });

                blk_start = blk_end;
            }
            std::for_each(threads.begin(),threads.end(), std::mem_fn(&std::thread::join));
        }
    }

    template<typename T1, typename T2, typename S>
    void Grid2D<T1,T2,S>::raytrace(const std::vector<std::vector<S>>& Tx,
                                   const std::vector<std::vector<T1>>& t0,
                                   const std::vector<std::vector<S>>& Rx,
                                   std::vector<std::vector<T1>>& traveltimes,
                                   std::vector<std::vector<std::vector<S>>>& r_data,
                                   std::vector<std::vector<std::vector<siv<T1>>>>& l_data) const {

        if ( Tx.size() == 1 ) {
            this->raytrace(Tx[0], t0[0], Rx[0], traveltimes[0], r_data[0], l_data[0], 0);
        } else if ( nThreads == 1 ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                this->raytrace(Tx[n], t0[n], Rx[n], traveltimes[n], r_data[n], l_data[n], 0);
            }
        } else if ( usePool ) {
            std::vector<std::future<void>> results(Tx.size());
            for ( size_t n=0; n<Tx.size(); ++n ) {
                results[n] = pool.push(std::ref(*this),
                                       std::ref(Tx[n]),
                                       std::ref(t0[n]),
                                       std::ref(Rx[n]),
                                       std::ref(traveltimes[n]),
                                       std::ref(r_data[n]),
                                       std::ref(l_data[n]));
            }
            for ( size_t n=0; n<Tx.size(); ++n ) {
                results[n].get();
            }
        } else {
            std::vector<size_t> blk_size = get_blk_size(Tx.size());

            std::vector<std::thread> threads(blk_size.size());
            size_t blk_start = 0;
            for ( size_t i=0; i<blk_size.size(); ++i ) {

                size_t blk_end = blk_start + blk_size[i];
                threads[i]=std::thread( [this,&Tx,&t0,&Rx,&traveltimes,&r_data,&l_data,blk_start,blk_end,i]{

                    for ( size_t n=blk_start; n<blk_end; ++n ) {
                        this->raytrace(Tx[n], t0[n], Rx[n], traveltimes[n], r_data[n], l_data[n], i);
                    }
                });

                blk_start = blk_end;
            }
            std::for_each(threads.begin(),threads.end(), std::mem_fn(&std::thread::join));
        }
    }

    template<typename T1, typename T2, typename S>
    void Grid2D<T1,T2,S>::raytrace(const std::vector<std::vector<S>>& Tx,
                                   const std::vector<std::vector<T1>>& t0,
                                   const std::vector<std::vector<S>>& Rx,
                                   std::vector<std::vector<T1>>& traveltimes,
                                   std::vector<std::vector<std::vector<S>>>& r_data,
                                   std::vector<std::vector<std::vector<sijv<T1>>>>& m_data) const {

        if ( Tx.size() == 1 ) {
            this->raytrace(Tx[0], t0[0], Rx[0], traveltimes[0], r_data[0], m_data[0], 0);
        } else if ( nThreads == 1 ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                this->raytrace(Tx[n], t0[n], Rx[n], traveltimes[n], r_data[n], m_data[n], 0);
            }
        } else if ( usePool ) {
            std::vector<std::future<void>> results(Tx.size());
            for ( size_t n=0; n<Tx.size(); ++n ) {
                results[n] = pool.push(std::ref(*this),
                                       std::ref(Tx[n]),
                                       std::ref(t0[n]),
                                       std::ref(Rx[n]),
                                       std::ref(traveltimes[n]),
                                       std::ref(r_data[n]),
                                       std::ref(m_data[n]));
            }
            for ( size_t n=0; n<Tx.size(); ++n ) {
                results[n].get();
            }
        } else {
            std::vector<size_t> blk_size = get_blk_size(Tx.size());

            std::vector<std::thread> threads(blk_size.size());
            size_t blk_start = 0;
            for ( size_t i=0; i<blk_size.size(); ++i ) {

                size_t blk_end = blk_start + blk_size[i];
                threads[i]=std::thread( [this,&Tx,&t0,&Rx,&traveltimes,&r_data,&m_data,blk_start,blk_end,i]{

                    for ( size_t n=blk_start; n<blk_end; ++n ) {
                        this->raytrace(Tx[n], t0[n], Rx[n], traveltimes[n], r_data[n], m_data[n], i);
                    }
                });

                blk_start = blk_end;
            }
            std::for_each(threads.begin(),threads.end(), std::mem_fn(&std::thread::join));
        }
    }

}

#endif

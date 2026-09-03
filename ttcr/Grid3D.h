// Grid3D is the class template of an object containing the 3D grid and the function
// to perform raytracing
//
//  Grid.h
//  ttcr
//
//  Created by Bernard Giroux on 2013-01-10.
//  Copyright (c) 2013 Bernard Giroux. All rights reserved.
//

/**
 * @file Grid3D.h
 * @brief Abstract base class template for 3D traveltime grids and raytracing.
 *
 * @ref ttcr::Grid3D is the common interface shared by every 3D grid/mesh
 * implementation in ttcr (rectilinear, node-based, dynamic-shortest-path,
 * fast-sweeping, etc.). It stores the data common to all of them (thread count,
 * optional origin translation, cell-to-node neighbor lists) and provides the
 * public @c raytrace() entry points, including the multi-source threaded
 * variants that fan work out over a thread pool or raw @c std::thread objects.
 *
 * Concrete grids override the protected virtual hooks (@c raytrace,
 * @c getTraveltime, @c getRaypath, ...) and the accessors describing their
 * geometry; the non-virtual public overloads in this header handle origin
 * translation, output-buffer sizing, and thread scheduling before delegating
 * to those hooks.
 */

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

#ifndef ttcr_Grid3D_h
#define ttcr_Grid3D_h

#include <algorithm>
#include <exception>
#include <functional>
#include <fstream>
#include <future>
#include <iostream>
#include <thread>
#include <vector>

#include "ctpl_stl.h"

#include "ttcr_t.h"

namespace ttcr {

    /**
     * @brief Abstract base for 3D traveltime grids and raytracing.
     * @tparam T1 Floating-point type for coordinates, slowness and traveltimes.
     * @tparam T2 Integer type used for node/cell indices.
     */
    template<typename T1, typename T2>
    class Grid3D {
    public:
        /**
         * @brief Constructs the shared grid state.
         * @param ttrp             If true, receiver traveltimes are recomputed
         *                         by integrating along the raypath rather than
         *                         read directly from the traveltime field.
         * @param ncells           Number of cells (sizes the neighbor table).
         * @param nt               Number of worker threads.
         * @param _translateOrigin If true, coordinates are shifted by @ref origin
         *                         before computation (improves conditioning).
         * @param _usePool         If true, use a persistent thread pool; otherwise
         *                         spawn raw threads per raytrace call.
         */
        Grid3D(const bool ttrp,
               const size_t ncells,
               const size_t nt=1,
               const bool _translateOrigin=false,
               const bool _usePool=1) :
        nThreads(nt), tt_from_rp(ttrp),
        translateOrigin(_translateOrigin),
        origin({0.0, 0.0, 0.0}),
        neighbors(std::vector<std::vector<T2>>(ncells)),
        usePool(_usePool)
        {
            if ( nThreads > 1 && usePool ) {
                pool.resize(static_cast<int>(nThreads));
            }
        }

        /// Virtual destructor.
        virtual ~Grid3D() {}

        /**
         * @name Thread-pool call operators
         * Callable entry points invoked by the thread pool; @p id is the
         * worker/thread number. Each simply forwards to the matching
         * @c raytrace() overload. The trailing output arguments select what
         * is computed in addition to traveltimes: @c r_data (raypaths),
         * @c m_data (model-parameter sensitivities), @c l_data (per-cell ray
         * length sensitivities).
         * @{
         */
        /// Traveltimes only.
        void operator()(int id, const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<sxyz<T1>>& Rx,
                        std::vector<T1>& traveltimes) const {
            this->raytrace(Tx, t0, Rx, traveltimes, id);
        }

        /// Traveltimes and raypaths.
        void operator()(int id, const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<sxyz<T1>>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<sxyz<T1>>>& r_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, r_data, id);
        }

        /// Traveltimes and model-parameter sensitivities.
        void operator()(int id, const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<sxyz<T1>>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<sijv<T1>>>& m_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, m_data, id);
        }

        /// Traveltimes, raypaths and model-parameter sensitivities.
        void operator()(int id, const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<sxyz<T1>>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<sxyz<T1>>>& r_data,
                        std::vector<std::vector<sijv<T1>>>& m_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, r_data, m_data, id);
        }

        /// Traveltimes and per-cell ray-length sensitivities.
        void operator()(int id, const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<sxyz<T1>>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<siv<T1>>>& l_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, l_data, id);
        }

        /// Traveltimes, raypaths and per-cell ray-length sensitivities.
        void operator()(int id, const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<sxyz<T1>>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<sxyz<T1>>>& r_data,
                        std::vector<std::vector<siv<T1>>>& l_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, r_data, l_data, id);
        }

        /// Traveltimes and two-component sensitivities (@ref siv2).
        void operator()(int id, const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<sxyz<T1>>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<siv2<T1>>>& l_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, l_data, id);
        }

        /// Traveltimes and four-component sensitivities (@ref siv4).
        void operator()(int id, const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<sxyz<T1>>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<siv4<T1>>>& l_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, l_data, id);
        }

        /// Traveltimes, raypaths and two-component sensitivities (@ref siv2).
        void operator()(int id, const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<sxyz<T1>>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<sxyz<T1>>>& r_data,
                        std::vector<std::vector<siv2<T1>>>& l_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, r_data, l_data, id);
        }

        /// Traveltimes, raypaths and four-component sensitivities (@ref siv4).
        void operator()(int id, const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<sxyz<T1>>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<sxyz<T1>>>& r_data,
                        std::vector<std::vector<siv4<T1>>>& l_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, r_data, l_data, id);
        }
        /** @} */

        /**
         * @brief Compute traveltimes from sources to receivers.
         * @param Tx          Source locations.
         * @param t0          Excitation time of each source.
         * @param Rx          Receiver locations.
         * @param[out] traveltimes Traveltime to each receiver (resized as needed).
         * @param threadNo    Worker/thread number.
         */
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              const size_t threadNo=0) const;

        /**
         * @brief Compute traveltimes to grouped receivers (one group per source).
         * @param Tx          Source locations.
         * @param t0          Excitation time of each source.
         * @param Rx          Receiver locations, grouped per source.
         * @param[out] traveltimes Pointers to per-group traveltime vectors.
         * @param threadNo    Worker/thread number.
         */
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<std::vector<sxyz<T1>>>& Rx,
                              std::vector<std::vector<T1>*>& traveltimes,
                              const size_t threadNo=0) const;

        /// As raytrace(), additionally returning the raypath to each receiver in @p r_data.
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<sxyz<T1>>>& r_data,
                              const size_t threadNo=0) const;

        /// Grouped-receiver variant additionally returning raypaths in @p r_data.
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<std::vector<sxyz<T1>>>& Rx,
                              std::vector<std::vector<T1>*>& traveltimes,
                              std::vector<std::vector<std::vector<sxyz<T1>>>*>& r_data,
                              const size_t threadNo=0) const;

        /// As raytrace(), additionally returning raypaths (@p r_data) and model sensitivities (@p m_data).
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<sxyz<T1>>>& r_data,
                              std::vector<std::vector<sijv<T1>>>& m_data,
                              const size_t threadNo=0) const;

        /// As raytrace(), additionally returning model-parameter sensitivities in @p m_data.
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<sijv<T1>>>& m_data,
                              const size_t threadNo=0) const;

        /// As raytrace(), additionally returning per-cell ray-length sensitivities in @p l_data.
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<siv<T1>>>& l_data,
                              const size_t threadNo=0) const;

        /// As raytrace(), additionally returning raypaths (@p r_data) and ray-length sensitivities (@p l_data).
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<sxyz<T1>>>& r_data,
                              std::vector<std::vector<siv<T1>>>& l_data,
                              const size_t threadNo=0) const;

        /// @brief Traveltimes and two-component sensitivities (@ref siv2). Subclass must override.
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<siv2<T1>>>& l_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method raytrace should be implemented in subclass");
        }

        /// @brief Traveltimes and four-component sensitivities (@ref siv4). Subclass must override.
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<siv4<T1>>>& l_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method raytrace should be implemented in subclass");
        }

        /// @brief Traveltimes, raypaths and two-component sensitivities (@ref siv2). Subclass must override.
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<sxyz<T1>>>& r_data,
                              std::vector<std::vector<siv2<T1>>>& l_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method raytrace should be implemented in subclass");
        }

        /// @brief Traveltimes, raypaths and four-component sensitivities (@ref siv4). Subclass must override.
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<sxyz<T1>>>& r_data,
                              std::vector<std::vector<siv4<T1>>>& l_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method raytrace should be implemented in subclass");
        }

        /**
         * @name Multi-source threaded raytracing
         * Each source set (outer vector index) is raytraced independently; the
         * work is distributed across @ref nThreads either via the thread pool
         * or raw threads. The output-argument variants mirror the single-source
         * overloads above (@c r_data raypaths, @c m_data / @c l_data
         * sensitivities).
         * @{
         */
        /// Traveltimes only, one source set per outer index.
        void raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes) const;

        /// Traveltimes and raypaths.
        void raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<sxyz<T1>>>>& r_data) const;

        /// Traveltimes and model-parameter sensitivities.
        void raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<sijv<T1>>>>& m_data) const;

        /// Traveltimes, raypaths and model-parameter sensitivities.
        void raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<sxyz<T1>>>>& r_data,
                      std::vector<std::vector<std::vector<sijv<T1>>>>& m_data) const;

        /// Traveltimes and per-cell ray-length sensitivities.
        void raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<siv<T1>>>>& l_data) const;

        /// Traveltimes, raypaths and per-cell ray-length sensitivities.
        void raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<sxyz<T1>>>>& r_data,
                      std::vector<std::vector<std::vector<siv<T1>>>>& l_data) const;

        /// Several sources, two-component sensitivities (@ref siv2).
        void raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<siv2<T1>>>>& l_data) const {
            raytraceMulti(Tx, t0, Rx, traveltimes, l_data);
        }

        /// Several sources, four-component sensitivities (@ref siv4).
        void raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<siv4<T1>>>>& l_data) const {
            raytraceMulti(Tx, t0, Rx, traveltimes, l_data);
        }

        /// Several sources, raypaths and two-component sensitivities (@ref siv2).
        void raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<sxyz<T1>>>>& r_data,
                      std::vector<std::vector<std::vector<siv2<T1>>>>& l_data) const {
            raytraceMulti(Tx, t0, Rx, traveltimes, r_data, l_data);
        }

        /// Several sources, raypaths and four-component sensitivities (@ref siv4).
        void raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<sxyz<T1>>>>& r_data,
                      std::vector<std::vector<std::vector<siv4<T1>>>>& l_data) const {
            raytraceMulti(Tx, t0, Rx, traveltimes, r_data, l_data);
        }
        /** @} */

        /**
         * @brief Compute straight-ray (no refraction) length sensitivities.
         * @param Tx          Source locations.
         * @param Rx          Receiver locations.
         * @param[out] l_data Per-cell ray-length contributions for each ray.
         * @throws std::runtime_error if not overridden by a subclass.
         */
        virtual void getStraightRays(const std::vector<sxyz<T1>>& Tx,
                                     const std::vector<sxyz<T1>>& Rx,
                                     std::vector<std::vector<siv<T1>>>& l_data) const {
            throw std::runtime_error("Method getStraightRays should be implemented in subclass");
        }

        /// Sets the slowness model from a vector (one value per node/cell).
        virtual void setSlowness(const std::vector<T1>& s) {}
        /// Sets the slowness model from a raw array of @p ns values.
        virtual void setSlowness(const T1 *s, const size_t ns) {
            throw std::runtime_error("Method getSlowness should be implemented in subclass");
        }
        /// Retrieves the current slowness model.
        virtual void getSlowness(std::vector<T1>&) const {
            throw std::runtime_error("Method getSlowness should be implemented in subclass");
        }
        /// Sets the anisotropy parameter &chi;, if supported.
        virtual void setChi(const std::vector<T1>& x) {}
        /// Sets the anisotropy parameter &psi;, if supported.
        virtual void setPsi(const std::vector<T1>& x) {}
        /// @brief Set the symmetry-axis tilt angle of the anisotropy. Subclass must override.
        virtual void setTiltAngle(const std::vector<T1>& x) {
            throw std::runtime_error("Method setTiltAngle should be implemented in subclass");
        }
        /// @brief Select the phase to model, qP or qSV. Subclass must override.
        virtual void setPhase(const int p) {
            throw std::runtime_error("Method setPhase should be implemented in subclass");
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

        /// Sets the finite source radius used for source regularization.
        virtual void setSourceRadius(const double) {}

        /// Returns the number of grid nodes.
        virtual size_t getNumberOfNodes() const { return 1; }
        /// Returns the number of grid cells.
        virtual size_t getNumberOfCells() const { return 1; }
        /// Copies the full traveltime field into @p tt for the given thread.
        virtual void getTT(std::vector<T1>& tt, const size_t threadNo=0) const {
            throw std::runtime_error("Method getTT should be implemented in subclass");
        }

        /// Saves the traveltime field to a file (subclass-defined format).
        virtual void saveTT(const std::string &, const int, const size_t nt=0,
                            const int format=1) const {}
        /// Loads a traveltime field from a file (subclass-defined format).
        virtual void loadTT(const std::string &, const int, const size_t nt=0,
                            const int format=1) const {}

        /// Returns the minimum x coordinate of the grid.
        virtual const T1 getXmin() const {
            throw std::runtime_error("Method getXmin should be implemented in subclass");
        }
        /// Returns the maximum x coordinate of the grid.
        virtual const T1 getXmax() const {
            throw std::runtime_error("Method getXmax should be implemented in subclass");
        }
        /// Returns the minimum y coordinate of the grid.
        virtual const T1 getYmin() const {
            throw std::runtime_error("Method getYmin should be implemented in subclass");
        }
        /// Returns the maximum y coordinate of the grid.
        virtual const T1 getYmax() const {
            throw std::runtime_error("Method getYmax should be implemented in subclass");
        }
        /// Returns the minimum z coordinate of the grid.
        virtual const T1 getZmin() const {
            throw std::runtime_error("Method getZmin should be implemented in subclass");
        }
        /// Returns the maximum z coordinate of the grid.
        virtual const T1 getZmax() const {
            throw std::runtime_error("Method getZmax should be implemented in subclass");
        }
        /// Returns the cell size along x.
        virtual const T1 getDx() const {
            throw std::runtime_error("Method getDx should be implemented in subclass");
        }
        /// Returns the cell size along y.
        virtual const T1 getDy() const {
            throw std::runtime_error("Method getDy should be implemented in subclass");
        }
        /// Returns the cell size along z.
        virtual const T1 getDz() const {
            throw std::runtime_error("Method getDz should be implemented in subclass");
        }
        /// Returns the number of cells along x.
        virtual const T2 getNcx() const {
            throw std::runtime_error("Method getNcx should be implemented in subclass");
        }
        /// Returns the number of cells along y.
        virtual const T2 getNcy() const {
            throw std::runtime_error("Method getNcy should be implemented in subclass");
        }
        /// Returns the number of cells along z.
        virtual const T2 getNcz() const {
            throw std::runtime_error("Method getNcz should be implemented in subclass");
        }
        /// Returns the number of secondary nodes per cell edge along x.
        virtual const T2 getNsnx() const {
            throw std::runtime_error("Method getNsnx should be implemented in subclass");
        }
        /// Returns the number of secondary nodes per cell edge along y.
        virtual const T2 getNsny() const {
            throw std::runtime_error("Method getNsny should be implemented in subclass");
        }
        /// Returns the number of secondary nodes per cell edge along z.
        virtual const T2 getNsnz() const {
            throw std::runtime_error("Method getNsnz should be implemented in subclass");
        }

        /// Returns the number of iterations performed (iterative solvers).
        virtual const int get_niter() const { return 0; }
        /// Returns the number of iterations performed within water (iterative solvers).
        virtual const int get_niterw() const { return 0; }

        /// Returns the number of worker threads.
        const size_t getNthreads() const { return nThreads; }

        /// Dumps secondary-node data to a stream (diagnostics; no-op by default).
        virtual void dump_secondary(std::ofstream&) const {}

        /**
         * compute slowness at given point
         *
         * @param pt point to consider
         * @param isTranslated point to consider has been translated (considered only if grid attribute translateOrigin == true)
         * @returns slowness value
         */
        virtual T1 computeSlowness(sxyz<T1> pt, const bool isTranslated=false) const {
            throw std::runtime_error("Method computeSlowness should be implemented in subclass");
        }

        /**
         * gradient of the traveltime field at an arbitrary point
         *
         * The field must already have been computed, i.e. this is only valid
         * after a call to raytrace() for the same thread.  The raypath tracers
         * step down this gradient; @ref computeH walks it back from a receiver
         * to recover the take-off direction at the source.
         *
         * @param[out] g gradient of the traveltime field at @p pt
         * @param pt point to consider
         * @param threadNo thread whose traveltime field to differentiate
         */
        virtual void getTraveltimeGradient(sxyz<T1>& g, const sxyz<T1>& pt,
                                           const size_t threadNo) const {
            throw std::runtime_error("Method getTraveltimeGradient should be implemented in subclass");
        }

        /**
         * traveltimes and the hypocentre-location Jacobian H
         *
         * This calls raytrace() itself: it returns the traveltimes it computed
         * along with H, so the caller needs a single call rather than a
         * raytrace() followed by a separate assembly of the derivatives.
         *
         * H holds the partial derivatives of the arrival time with respect to
         * the hypocentre parameters, one row per receiver.  With @p full,
         * the columns are
         *
         *     [ 1, dT/dx, dT/dy, dT/dz ]
         *
         * the leading 1 being the derivative with respect to the origin time;
         * otherwise the two columns are [ dT/dx, dT/dy ].
         *
         * The spatial derivatives follow from reciprocity: the traveltime
         * gradient at a point is the slowness vector of the ray through it, so
         *
         *     dT/dx_s = -s(x_s) * e
         *
         * with @c e the unit take-off direction at the source.  @c e is
         * obtained by walking down the traveltime field from the receiver
         * until the source is approached, rather than from the stored raypath:
         * the raypath endpoint convention differs between solvers, and its
         * last segment is a noisy estimate of the tangent.
         *
         * @param Tx source positions
         * @param t0 source times
         * @param Rx receiver positions
         * @param[out] traveltimes traveltimes computed by the internal raytrace()
         * @param[out] H Jacobian, Rx.size() rows of 4 (or 2) values
         * @param full build the four-column form rather than the two-column one
         * @param radius_factor take-off is measured where the backward walk
         *        first comes within radius_factor average edge lengths of the
         *        source
         * @param threadNo thread to use
         * @throws runtime_error if the grid cannot supply traveltime gradients
         */
        void computeH(const std::vector<sxyz<T1>>& Tx,
                      const std::vector<T1>& t0,
                      const std::vector<sxyz<T1>>& Rx,
                      std::vector<T1>& traveltimes,
                      std::vector<std::vector<T1>>& H,
                      const bool full=true,
                      const T1 radius_factor=3.0,
                      const size_t threadNo=0) const;

        /// Enables/disables the persistent thread pool, resizing it if needed.
        void setUsePool(const bool up) {
            usePool = up;
            if ( nThreads > 1 && usePool && pool.size() != nThreads ) {
                pool.resize(static_cast<int>(nThreads));
            }
        }

        /// Selects whether receiver traveltimes are obtained by integrating along the raypath.
        void setTraveltimeFromRaypath(const bool ttrp) { tt_from_rp = ttrp; }

        /**
         * check if points are inside grid
         *
         * @param pts points to consider
         * @param translated point to consider has been translated (considered only if grid attribute translateOrigin == true)
         * @throws runtime_error if a point is outside grid
         */
        virtual void checkPts(std::vector<sxyz<T1>> pts, const bool translated=false) const {
            throw std::runtime_error("Method checkPts should be implemented in subclass");
        }

        /**
         * compute terms of matrix of interpolation weights for velocity data points constraint
         *
         * @param pts points to consider
         * @param d_data values to fill matrix D
         * @param translated point to consider has been translated (considered only if grid attribute translateOrigin == true)
         * @throws runtime_error if a point is outside grid
         */
        virtual void computeD(std::vector<sxyz<T1>> pts,
                              std::vector<std::vector<sijv<T1>>> &d_data,
                              const bool translated=false) const {
            throw std::runtime_error("Method computeD should be implemented in subclass");
        }

        /**
         * @brief Compute the smoothing/roughness operator K for the model.
         * @param[out] d_data          Sparse operator terms, per dimension.
         * @param order                Derivative order of the operator.
         * @param taylorSeriesOrder    Order of the Taylor expansion used.
         * @param weighting            Enable distance-based weighting.
         * @param s0inside             Include the central point in the stencil.
         * @param additionnalPoints    Extra stencil points on each side.
         * @throws std::runtime_error if not overridden by a subclass.
         */
        virtual void computeK(std::vector<std::vector<std::vector<siv<T1>>>>& d_data,
                              const int order, const int taylorSeriesOrder,
                              const bool weighting, const bool s0inside,
                              const int additionnalPoints) const {
            throw std::runtime_error("Method computeK should be implemented in subclass");
        }

        /// Returns the average edge length of the grid/mesh.
        virtual const T1 getAverageEdgeLength() const {
            throw std::runtime_error("Method getAverageEdgeLength should be implemented in subclass");
        }

#ifdef VTK
        /// Saves the model as an unstructured-grid VTK file (.vtu). Requires VTK.
        virtual void saveModelVTU(const std::string &, const bool saveSlowness=true,
                                  const bool savePhysicalEntity=false) const {}
        /// Saves the model as a rectilinear-grid VTK file (.vtr). Requires VTK.
        virtual void saveModelVTR(const std::string &,
                                  const bool saveSlowness=true) const {}
        /// Saves the model as a .vtr file using externally supplied data. Requires VTK.
        virtual void saveModelVTR(const std::string &, const double*,
                                  const bool saveSlowness=true) const {}
#endif
        /**
         * @brief Spread several sources over the threads, whatever the container
         *
         * The work is distributed as in the overloads taking @ref siv, which
         * predate the wider containers and keep their own copy of it.
         */
        template<typename SIV>
        void raytraceMulti(const std::vector<std::vector<sxyz<T1>>>& Tx,
                           const std::vector<std::vector<T1>>& t0,
                           const std::vector<std::vector<sxyz<T1>>>& Rx,
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
        void raytraceMulti(const std::vector<std::vector<sxyz<T1>>>& Tx,
                           const std::vector<std::vector<T1>>& t0,
                           const std::vector<std::vector<sxyz<T1>>>& Rx,
                           std::vector<std::vector<T1>>& traveltimes,
                           std::vector<std::vector<std::vector<sxyz<T1>>>>& r_data,
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
        size_t nThreads;         ///< Number of worker threads.
        bool tt_from_rp;         ///< If true, receiver traveltimes come from raypath integration.
        bool translateOrigin;    ///< If true, coordinates are shifted by @ref origin.
        sxyz<T1> origin;         ///< Origin offset applied when @ref translateOrigin is set.
        std::vector<std::vector<T2>> neighbors;  ///< For each cell, the indices of its nodes.

        /**
         * @brief Populate @ref neighbors from the nodes' cell-ownership lists.
         * @tparam N   Node type exposing @c getOwners().
         * @param nodes The grid nodes.
         */
        template<typename N>
        void buildGridNeighbors(const std::vector<N>& nodes) {
            //Index the neighbors nodes of each cell
            for ( T2 n=0; n<nodes.size(); ++n ) {
                for ( size_t n2=0; n2<nodes[n].getOwners().size(); ++n2) {
                    neighbors[ nodes[n].getOwners()[n2] ].push_back(n);
                }
            }
        }

        /**
         * @brief Core solver hook: propagate traveltimes from the sources.
         * @param Tx       Source locations (already origin-translated).
         * @param t0       Source excitation times.
         * @param Rx       Receiver locations (used by some solvers to bound work).
         * @param threadNo Worker/thread number.
         * @note Must be overridden; the public raytrace() overloads call this
         *       after translating coordinates and sizing output buffers.
         */
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method raytrace should be implemented in subclass");
        }

        /// Core solver hook, grouped-receiver form.
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<std::vector<sxyz<T1>>>& Rx,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method raytrace should be implemented in subclass");
        }

        /// Returns the traveltime interpolated at point @p pt from the solved field.
        virtual T1 getTraveltime(const sxyz<T1>& pt,
                                 const size_t threadNo) const {
            throw std::runtime_error("Method getTraveltime should be implemented in subclass");
        }

        /// Returns the receiver traveltime obtained by integrating along the raypath.
        virtual T1 getTraveltimeFromRaypath(const std::vector<sxyz<T1>>& Tx,
                                            const std::vector<T1>& t0,
                                            const sxyz<T1>& Rx,
                                            const size_t threadNo) const {
            throw std::runtime_error("Method getTraveltimeFromRaypath should be implemented in subclass");
        }

        /**
         * @brief Trace the raypath from a receiver back to the source.
         * @param Tx       Source locations.
         * @param t0       Source excitation times.
         * @param Rx       Receiver location.
         * @param[out] r_data Ordered points of the raypath.
         * @param[out] tt  Traveltime along the recovered raypath.
         * @param threadNo Worker/thread number.
         */
        virtual void getRaypath(const std::vector<sxyz<T1>>& Tx,
                                const std::vector<T1>& t0,
                                const sxyz<T1>& Rx,
                                std::vector<sxyz<T1>>& r_data,
                                T1 &tt,
                                const size_t threadNo) const {
            throw std::runtime_error("Method getRaypath should be implemented in subclass");
        }

        /// Raypath variant returning model-parameter sensitivities @p m_data for receiver @p RxNo.
        virtual void getRaypath(const std::vector<sxyz<T1>>& Tx,
                                const std::vector<T1>& t0,
                                const sxyz<T1>& Rx,
                                std::vector<sijv<T1>>& m_data,
                                const size_t RxNo,
                                T1 &tt,
                                const size_t threadNo) const {
            throw std::runtime_error("Method getRaypath should be implemented in subclass");
        }

        /// Raypath variant returning both raypath points @p r_data and model sensitivities @p m_data.
        virtual void getRaypath(const std::vector<sxyz<T1>>& Tx,
                                const std::vector<T1>& t0,
                                const sxyz<T1>& Rx,
                                std::vector<sxyz<T1>>& r_data,
                                std::vector<sijv<T1>>& m_data,
                                const size_t RxNo,
                                T1 &tt,
                                const size_t threadNo) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }

        /// Raypath variant returning per-cell ray-length sensitivities @p l_data.
        virtual void getRaypath(const std::vector<sxyz<T1>>& Tx,
                                const std::vector<T1>& t0,
                                const sxyz<T1> &Rx,
                                std::vector<siv<T1>> &l_data,
                                T1 &tt,
                                const size_t threadNo) const {
            throw std::runtime_error("Method getRaypathshould be implemented in subclass");
        }

        /// Raypath variant returning both raypath points @p r_data and ray-length sensitivities @p l_data.
        virtual void getRaypath(const std::vector<sxyz<T1>>& Tx,
                                const std::vector<T1>& t0,
                                const sxyz<T1> &Rx,
                                std::vector<sxyz<T1>> &r_data,
                                std::vector<siv<T1>> &l_data,
                                T1 &tt,
                                const size_t threadNo) const {
            throw std::runtime_error("Method getRaypath should be implemented in subclass");
        }

    private:
        bool usePool;                     ///< Whether to use the persistent thread pool.
        mutable ctpl::thread_pool pool;   ///< Reusable pool of worker threads.

        /**
         * @brief Partition @p nTx source sets into per-thread block sizes.
         * @param nTx Number of source sets to distribute.
         * @return Block size for each thread (round-robin, near-even split).
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


    template<typename T1, typename T2>
    void Grid3D<T1,T2>::computeH(const std::vector<sxyz<T1>>& Tx,
                                 const std::vector<T1>& t0,
                                 const std::vector<sxyz<T1>>& Rx,
                                 std::vector<T1>& traveltimes,
                                 std::vector<std::vector<T1>>& H,
                                 const bool full,
                                 const T1 radius_factor,
                                 const size_t threadNo) const {

        if ( Tx.empty() ) {
            throw std::runtime_error("computeH: no source");
        }
        if ( t0.size() != Tx.size() ) {
            throw std::runtime_error("computeH: Tx and t0 must have the same size");
        }
        if ( radius_factor <= static_cast<T1>(0.0) ) {
            throw std::runtime_error("computeH: radius_factor must be positive");
        }

        // this also leaves the traveltime field that the walk below descends
        this->raytrace(Tx, t0, Rx, traveltimes, threadNo);

        // raytrace() shifts its arguments internally, so that field lives in
        // translated coordinates; the walk has to work in the same frame.
        std::vector<sxyz<T1>> sTx = Tx;
        std::vector<sxyz<T1>> sRx = Rx;
        if ( translateOrigin ) {
            for ( size_t n=0; n<sTx.size(); ++n ) {
                sTx[n] -= origin;
            }
            for ( size_t n=0; n<sRx.size(); ++n ) {
                sRx[n] -= origin;
            }
        }

        const size_t ncol = full ? 4 : 2;
        const T1 h = getAverageEdgeLength();
        const T1 r_stop = radius_factor * h;
        const T1 step = static_cast<T1>(0.5) * h;
        const T1 small = static_cast<T1>(1.e-30);

        H.resize( Rx.size() );

        for ( size_t n=0; n<Rx.size(); ++n ) {

            H[n].assign(ncol, static_cast<T1>(0.0));
            if ( full ) {
                H[n][0] = static_cast<T1>(1.0);  // d(arrival)/d(t0)
            }

            // source this receiver is closest to, used both as the walk's
            // target and as the point the take-off is measured from
            size_t iTx = 0;
            T1 d_src = sRx[n].getDistance(sTx[0]);
            for ( size_t k=1; k<sTx.size(); ++k ) {
                T1 d = sRx[n].getDistance(sTx[k]);
                if ( d < d_src ) {
                    d_src = d;
                    iTx = k;
                }
            }

            // walk down the traveltime field, receiver towards source, and
            // take the direction where we first come within r_stop.  Closer
            // than that the field is radially degenerate about the source and
            // carries no information about which ray reaches this receiver.
            sxyz<T1> p = sRx[n];
            sxyz<T1> e = sRx[n] - sTx[iTx];      // straight-line fallback
            const size_t max_steps = static_cast<size_t>(20.0 * d_src / step) + 100;

            for ( size_t it=0; it<max_steps; ++it ) {

                T1 d = p.getDistance(sTx[iTx]);
                for ( size_t k=0; k<sTx.size(); ++k ) {
                    T1 dk = p.getDistance(sTx[k]);
                    if ( dk < d ) {
                        d = dk;
                        iTx = k;
                    }
                }
                if ( d <= r_stop ) {
                    e = p - sTx[iTx];
                    break;
                }

                sxyz<T1> g;
                getTraveltimeGradient(g, p, threadNo);
                T1 gn = norm(g);
                if ( gn < small ) {
                    break;                        // stalled: keep the fallback
                }
                p = p - (step/gn) * g;
            }

            T1 en = norm(e);
            if ( en < small ) {
                continue;                         // receiver sits on the source
            }
            e = e / en;

            const T1 s0 = computeSlowness(sTx[iTx], translateOrigin);
            if ( full ) {
                H[n][1] = -s0 * e.x;
                H[n][2] = -s0 * e.y;
                H[n][3] = -s0 * e.z;
            } else {
                H[n][0] = -s0 * e.x;
                H[n][1] = -s0 * e.y;
            }
        }
    }

    template<typename T1, typename T2>
    void Grid3D<T1,T2>::raytrace(const std::vector<sxyz<T1>>& _Tx,
                                 const std::vector<T1>& t0,
                                 const std::vector<sxyz<T1>>& _Rx,
                                 std::vector<T1>& traveltimes,
                                 const size_t threadNo) const {

        std::vector<sxyz<T1>> Tx = _Tx;
        std::vector<sxyz<T1>> Rx = _Rx;
        if ( translateOrigin ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                Tx[n] -= origin;
            }
            for ( size_t n=0; n<Rx.size(); ++n ) {
                Rx[n] -= origin;
            }
        }

        this->raytrace(Tx, t0, Rx, threadNo);

        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }

        if ( this->tt_from_rp ) {
            for (size_t n=0; n<Rx.size(); ++n) {
                traveltimes[n] = this->getTraveltimeFromRaypath(Tx, t0, Rx[n], threadNo);
            }
        } else {
            for (size_t n=0; n<Rx.size(); ++n) {
                traveltimes[n] = this->getTraveltime(Rx[n], threadNo);
            }
        }
    }

    template<typename T1, typename T2>
    void Grid3D<T1,T2>::raytrace(const std::vector<sxyz<T1>>& _Tx,
                                 const std::vector<T1>& t0,
                                 const std::vector<std::vector<sxyz<T1>>>& _Rx,
                                 std::vector<std::vector<T1>*>& traveltimes,
                                 const size_t threadNo) const {
        std::vector<sxyz<T1>> Tx = _Tx;
        std::vector<std::vector<sxyz<T1>>> Rx = _Rx;
        if ( translateOrigin ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                Tx[n] -= origin;
            }
            for ( size_t n=0; n<Rx.size(); ++n ) {
                for ( size_t nn=0; nn<Rx[n].size(); ++nn ) {
                    Rx[n][nn] = _Rx[n][nn] - origin;
                }
            }
        }

        this->raytrace(Tx, t0, Rx, threadNo);

        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }

        if ( this->tt_from_rp ) {
            for (size_t nr=0; nr<Rx.size(); ++nr) {
                traveltimes[nr]->resize( Rx[nr].size() );
                for (size_t n=0; n<Rx[nr].size(); ++n) {
                    (*traveltimes[nr])[n] = this->getTraveltimeFromRaypath(Tx, t0, Rx[nr][n], threadNo);
                }
            }
        } else {
            for (size_t nr=0; nr<Rx.size(); ++nr) {
                traveltimes[nr]->resize( Rx[nr].size() );
                for (size_t n=0; n<Rx[nr].size(); ++n) {
                    (*traveltimes[nr])[n] = this->getTraveltime(Rx[nr][n], threadNo);
                }
            }
        }
    }

    template<typename T1, typename T2>
    void Grid3D<T1,T2>::raytrace(const std::vector<sxyz<T1>>& _Tx,
                                 const std::vector<T1>& t0,
                                 const std::vector<sxyz<T1>>& _Rx,
                                 std::vector<T1>& traveltimes,
                                 std::vector<std::vector<sxyz<T1>>>& r_data,
                                 const size_t threadNo) const {

        std::vector<sxyz<T1>> Tx = _Tx;
        std::vector<sxyz<T1>> Rx = _Rx;
        if ( translateOrigin ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                Tx[n] -= origin;
            }
            for ( size_t n=0; n<Rx.size(); ++n ) {
                Rx[n] -= origin;
            }
        }

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
        if ( translateOrigin ) {
            for (size_t n=0; n<r_data.size(); ++n) {
                for (size_t nn=0; nn<r_data[n].size(); ++nn) {
                    r_data[n][nn] += origin;
                }
            }
        }
    }

    template<typename T1, typename T2>
    void Grid3D<T1,T2>::raytrace(const std::vector<sxyz<T1>>& _Tx,
                                 const std::vector<T1>& t0,
                                 const std::vector<std::vector<sxyz<T1>>>& _Rx,
                                 std::vector<std::vector<T1>*>& traveltimes,
                                 std::vector<std::vector<std::vector<sxyz<T1>>>*>& r_data,
                                 const size_t threadNo) const {
        std::vector<sxyz<T1>> Tx = _Tx;
        std::vector<std::vector<sxyz<T1>>> Rx = _Rx;
        if ( translateOrigin ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                Tx[n] -= origin;
            }
            for ( size_t n=0; n<Rx.size(); ++n ) {
                for ( size_t nn=0; nn<Rx[n].size(); ++nn ) {
                    Rx[n][nn] = _Rx[n][nn] - origin;
                }
            }
        }
        
        this->raytrace(Tx, t0, Rx, threadNo);

        if ( r_data.size() != Rx.size() ) {
            r_data.resize( Rx.size() );
        }
        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }

        for (size_t nr=0; nr<Rx.size(); ++nr) {
            r_data[nr]->resize( Rx[nr].size() );
            for ( size_t ni=0; ni<r_data[nr]->size(); ++ni ) {
                (*r_data[nr])[ni].resize( 0 );
            }
            traveltimes[nr]->resize( Rx[nr].size() );

            for (size_t n=0; n<Rx[nr].size(); ++n) {
                this->getRaypath(Tx, t0, Rx[nr][n], (*r_data[nr])[n],
                                 (*traveltimes[nr])[n], threadNo);
            }
        }
        if ( translateOrigin ) {
            for (size_t n=0; n<r_data.size(); ++n) {
                for (size_t nn=0; nn<r_data[n]->size(); ++nn) {
                    for (size_t nnn=0; nnn<(*r_data[n])[nn].size(); ++nnn) {
                        (*r_data[n])[nn][nnn] += origin;
                    }
                }
            }
        }
    }

    template<typename T1, typename T2>
    void Grid3D<T1,T2>::raytrace(const std::vector<sxyz<T1>>& _Tx,
                                 const std::vector<T1>& t0,
                                 const std::vector<sxyz<T1>>& _Rx,
                                 std::vector<T1>& traveltimes,
                                 std::vector<std::vector<sxyz<T1>>>& r_data,
                                 std::vector<std::vector<sijv<T1>>>& m_data,
                                 const size_t threadNo) const {
        std::vector<sxyz<T1>> Tx = _Tx;
        std::vector<sxyz<T1>> Rx = _Rx;
        if ( translateOrigin ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                Tx[n] -= origin;
            }
            for ( size_t n=0; n<Rx.size(); ++n ) {
                Rx[n] -= origin;
            }
        }

        this->raytrace(Tx, t0, Rx, threadNo);

        if ( r_data.size() != Rx.size() ) {
            r_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<r_data.size(); ++ni ) {
            r_data[ni].resize( 0 );
        }
        if ( m_data.size() != Rx.size() ) {
            m_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<m_data.size(); ++ni ) {
            m_data[ni].resize( 0 );
        }
        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }

        for (size_t n=0; n<Rx.size(); ++n) {
            this->getRaypath(Tx, t0, Rx[n], r_data[n], m_data[n], n, traveltimes[n], threadNo);
        }
        if ( translateOrigin ) {
            for (size_t n=0; n<r_data.size(); ++n) {
                for (size_t nn=0; nn<r_data[n].size(); ++nn) {
                    r_data[n][nn] += origin;
                }
            }
        }
    }

    template<typename T1, typename T2>
    void Grid3D<T1,T2>::raytrace(const std::vector<sxyz<T1>>& _Tx,
                                 const std::vector<T1>& t0,
                                 const std::vector<sxyz<T1>>& _Rx,
                                 std::vector<T1>& traveltimes,
                                 std::vector<std::vector<sxyz<T1>>>& r_data,
                                 std::vector<std::vector<siv<T1>>>& l_data,
                                 const size_t threadNo) const {
        std::vector<sxyz<T1>> Tx = _Tx;
        std::vector<sxyz<T1>> Rx = _Rx;
        if ( translateOrigin ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                Tx[n] -= origin;
            }
            for ( size_t n=0; n<Rx.size(); ++n ) {
                Rx[n] -= origin;
            }
        }

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
        }
        if ( translateOrigin ) {
            for (size_t n=0; n<r_data.size(); ++n) {
                for (size_t nn=0; nn<r_data[n].size(); ++nn) {
                    r_data[n][nn] += origin;
                }
            }
        }
    }

    template<typename T1, typename T2>
    void Grid3D<T1,T2>::raytrace(const std::vector<sxyz<T1>>& _Tx,
                                 const std::vector<T1>& t0,
                                 const std::vector<sxyz<T1>>& _Rx,
                                 std::vector<T1>& traveltimes,
                                 std::vector<std::vector<sijv<T1>>>& m_data,
                                 const size_t threadNo) const {
        std::vector<sxyz<T1>> Tx = _Tx;
        std::vector<sxyz<T1>> Rx = _Rx;
        if ( translateOrigin ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                Tx[n] -= origin;
            }
            for ( size_t n=0; n<Rx.size(); ++n ) {
                Rx[n] -= origin;
            }
        }

        this->raytrace(Tx, t0, Rx, threadNo);

        if ( m_data.size() != Rx.size() ) {
            m_data.resize( Rx.size() );
        }
        for ( size_t ni=0; ni<m_data.size(); ++ni ) {
            m_data[ni].resize( 0 );
        }
        if ( traveltimes.size() != Rx.size() ) {
            traveltimes.resize( Rx.size() );
        }

        for (size_t n=0; n<Rx.size(); ++n) {
            this->getRaypath(Tx, t0, Rx[n], m_data[n], n, traveltimes[n], threadNo);
        }
    }

    template<typename T1, typename T2>
    void Grid3D<T1,T2>::raytrace(const std::vector<sxyz<T1>>& _Tx,
                                 const std::vector<T1>& t0,
                                 const std::vector<sxyz<T1>>& _Rx,
                                 std::vector<T1>& traveltimes,
                                 std::vector<std::vector<siv<T1>>>& l_data,
                                 const size_t threadNo) const {

        std::vector<sxyz<T1>> Tx = _Tx;
        std::vector<sxyz<T1>> Rx = _Rx;
        if ( translateOrigin ) {
            for ( size_t n=0; n<Tx.size(); ++n ) {
                Tx[n] -= origin;
            }
            for ( size_t n=0; n<Rx.size(); ++n ) {
                Rx[n] -= origin;
            }
        }

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
        }
    }

    template<typename T1, typename T2>
    void Grid3D<T1,T2>::raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                                 const std::vector<std::vector<T1>>& t0,
                                 const std::vector<std::vector<sxyz<T1>>>& Rx,
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

    template<typename T1, typename T2>
    void Grid3D<T1,T2>::raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                                 const std::vector<std::vector<T1>>& t0,
                                 const std::vector<std::vector<sxyz<T1>>>& Rx,
                                 std::vector<std::vector<T1>>& traveltimes,
                                 std::vector<std::vector<std::vector<sxyz<T1>>>>& r_data) const {

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

    template<typename T1, typename T2>
    void Grid3D<T1,T2>::raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                                 const std::vector<std::vector<T1>>& t0,
                                 const std::vector<std::vector<sxyz<T1>>>& Rx,
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

    template<typename T1, typename T2>
    void Grid3D<T1,T2>::raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                                 const std::vector<std::vector<T1>>& t0,
                                 const std::vector<std::vector<sxyz<T1>>>& Rx,
                                 std::vector<std::vector<T1>>& traveltimes,
                                 std::vector<std::vector<std::vector<sxyz<T1>>>>& r_data,
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

    template<typename T1, typename T2>
    void Grid3D<T1,T2>::raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                                 const std::vector<std::vector<T1>>& t0,
                                 const std::vector<std::vector<sxyz<T1>>>& Rx,
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

    template<typename T1, typename T2>
    void Grid3D<T1,T2>::raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                                 const std::vector<std::vector<T1>>& t0,
                                 const std::vector<std::vector<sxyz<T1>>>& Rx,
                                 std::vector<std::vector<T1>>& traveltimes,
                                 std::vector<std::vector<std::vector<sxyz<T1>>>>& r_data,
                                 std::vector<std::vector<std::vector<siv<T1>>>>& l_data) const {
        if ( verbose > 2 ) {
            std::cout << "\nIn Grid3D::raytrace\n" << std::endl;
        }
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

}


#endif

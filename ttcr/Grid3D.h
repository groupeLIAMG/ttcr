// Grid3D is the class template of an object containing the 3D grid and the function
// to perform raytracing
//
//  Grid.h
//  ttcr
//
//  Created by Bernard Giroux on 2013-01-10.
//  Copyright (c) 2013 Bernard Giroux. All rights reserved.
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

#ifndef ttcr_Grid3D_h
#define ttcr_Grid3D_h

#include <algorithm>
#include <exception>
#include <functional>
#include <fstream>
#include <future>
#include <thread>
#include <vector>

#include "ctpl_stl.h"

#include "ttcr_t.h"

namespace ttcr {

    template<typename T1, typename T2>
    class Grid3D {
    public:
        Grid3D(const bool ttrp,
               const size_t ncells,
               const size_t nt=1,
               const bool _translateOrigin=false,
               const bool _usePool=1) :
        usePool(_usePool), nThreads(nt), tt_from_rp(ttrp),
        translateOrigin(_translateOrigin),
        origin({0.0, 0.0, 0.0}),
        neighbors(std::vector<std::vector<T2>>(ncells))
        {
            if ( nThreads > 1 && usePool ) {
                pool.resize(static_cast<int>(nThreads));
            }
        }

        virtual ~Grid3D() {}

        // operators: for usage with thread pool
        void operator()(int id, const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<sxyz<T1>>& Rx,
                        std::vector<T1>& traveltimes) const {
            this->raytrace(Tx, t0, Rx, traveltimes, id);
        }

        void operator()(int id, const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<sxyz<T1>>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<sxyz<T1>>>& r_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, r_data, id);
        }

        void operator()(int id, const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<sxyz<T1>>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<sijv<T1>>>& m_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, m_data, id);
        }

        void operator()(int id, const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<sxyz<T1>>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<sxyz<T1>>>& r_data,
                        std::vector<std::vector<sijv<T1>>>& m_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, r_data, m_data, id);
        }

        void operator()(int id, const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<sxyz<T1>>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<siv<T1>>>& l_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, l_data, id);
        }

        void operator()(int id, const std::vector<sxyz<T1>>& Tx,
                        const std::vector<T1>& t0,
                        const std::vector<sxyz<T1>>& Rx,
                        std::vector<T1>& traveltimes,
                        std::vector<std::vector<sxyz<T1>>>& r_data,
                        std::vector<std::vector<siv<T1>>>& l_data) const {
            this->raytrace(Tx, t0, Rx, traveltimes, r_data, l_data, id);
        }

        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              const size_t threadNo=0) const;

        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<std::vector<sxyz<T1>>>& Rx,
                              std::vector<std::vector<T1>*>& traveltimes,
                              const size_t threadNo=0) const;

        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<sxyz<T1>>>& r_data,
                              const size_t threadNo=0) const;

        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<std::vector<sxyz<T1>>>& Rx,
                              std::vector<std::vector<T1>*>& traveltimes,
                              std::vector<std::vector<std::vector<sxyz<T1>>>*>& r_data,
                              const size_t threadNo=0) const;

        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<sxyz<T1>>>& r_data,
                              std::vector<std::vector<sijv<T1>>>& m_data,
                              const size_t threadNo=0) const;

        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<sijv<T1>>>& m_data,
                              const size_t threadNo=0) const;

        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<siv<T1>>>& l_data,
                              const size_t threadNo=0) const;

        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<sxyz<T1>>>& r_data,
                              std::vector<std::vector<siv<T1>>>& l_data,
                              const size_t threadNo=0) const;

        // methods for threaded raytracing
        void raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes) const;

        void raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<sxyz<T1>>>>& r_data) const;

        void raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<sijv<T1>>>>& m_data) const;

        void raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<sxyz<T1>>>>& r_data,
                      std::vector<std::vector<std::vector<sijv<T1>>>>& m_data) const;

        void raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<siv<T1>>>>& l_data) const;

        void raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                      const std::vector<std::vector<T1>>& t0,
                      const std::vector<std::vector<sxyz<T1>>>& Rx,
                      std::vector<std::vector<T1>>& traveltimes,
                      std::vector<std::vector<std::vector<sxyz<T1>>>>& r_data,
                      std::vector<std::vector<std::vector<siv<T1>>>>& l_data) const;

        virtual void getStraightRays(const std::vector<sxyz<T1>>& Tx,
                                     const std::vector<sxyz<T1>>& Rx,
                                     std::vector<std::vector<siv<T1>>>& l_data) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }

        virtual void setSlowness(const std::vector<T1>& s) {}
        virtual void setSlowness(const T1 *s, const size_t ns) {
            throw std::runtime_error("Method should be implemented in subclass");
        }
        virtual void getSlowness(std::vector<T1>&) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }
        virtual void setChi(const std::vector<T1>& x) {}
        virtual void setPsi(const std::vector<T1>& x) {}

        virtual void setSourceRadius(const double) {}

        virtual size_t getNumberOfNodes() const { return 1; }
        virtual size_t getNumberOfCells() const { return 1; }
        virtual void getTT(std::vector<T1>& tt, const size_t threadNo=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }

        virtual void saveTT(const std::string &, const int, const size_t nt=0,
                            const int format=1) const {}
        virtual void loadTT(const std::string &, const int, const size_t nt=0,
                            const int format=1) const {}

        virtual const T1 getXmin() const {
            throw std::runtime_error("Method computeSlowness should be implemented in subclass");
        }
        virtual const T1 getXmax() const {
            throw std::runtime_error("Method computeSlowness should be implemented in subclass");
        }
        virtual const T1 getYmin() const {
            throw std::runtime_error("Method computeSlowness should be implemented in subclass");
        }
        virtual const T1 getYmax() const {
            throw std::runtime_error("Method computeSlowness should be implemented in subclass");
        }
        virtual const T1 getZmin() const {
            throw std::runtime_error("Method computeSlowness should be implemented in subclass");
        }
        virtual const T1 getZmax() const {
            throw std::runtime_error("Method computeSlowness should be implemented in subclass");
        }
        virtual const T1 getDx() const {
            throw std::runtime_error("Method computeSlowness should be implemented in subclass");
        }
        virtual const T1 getDy() const {
            throw std::runtime_error("Method computeSlowness should be implemented in subclass");
        }
        virtual const T1 getDz() const {
            throw std::runtime_error("Method computeSlowness should be implemented in subclass");
        }
        virtual const T2 getNcx() const {
            throw std::runtime_error("Method computeSlowness should be implemented in subclass");
        }
        virtual const T2 getNcy() const {
            throw std::runtime_error("Method computeSlowness should be implemented in subclass");
        }
        virtual const T2 getNcz() const {
            throw std::runtime_error("Method computeSlowness should be implemented in subclass");
        }
        virtual const T2 getNsnx() const {
            throw std::runtime_error("Method computeSlowness should be implemented in subclass");
        }
        virtual const T2 getNsny() const {
            throw std::runtime_error("Method computeSlowness should be implemented in subclass");
        }
        virtual const T2 getNsnz() const {
            throw std::runtime_error("Method computeSlowness should be implemented in subclass");
        }

        virtual const int get_niter() const { return 0; }
        virtual const int get_niterw() const { return 0; }

        const size_t getNthreads() const { return nThreads; }

        virtual void dump_secondary(std::ofstream&) const {}

        virtual T1 computeSlowness(const sxyz<T1>&) const {
            throw std::runtime_error("Method computeSlowness should be implemented in subclass");
        }

        void setUsePool(const bool up) {
            usePool = up;
            if ( nThreads > 1 && usePool && pool.size() != nThreads ) {
                pool.resize(static_cast<int>(nThreads));
            }
        }

        void setTraveltimeFromRaypath(const bool ttrp) { tt_from_rp = ttrp; }

        virtual void checkPts(const std::vector<sxyz<T1>>&) const {
            throw std::runtime_error("Method checkPts should be implemented in subclass");
        }

        virtual void computeD(const std::vector<sxyz<T1>> &pts,
                              std::vector<std::vector<sijv<T1>>> &d_data) const {
            throw std::runtime_error("Method computeD should be implemented in subclass");
        }

        virtual void computeK(std::vector<std::vector<std::vector<siv<T1>>>>& d_data,
                              const int order, const int taylorSeriesOrder,
                              const bool weighting, const bool s0inside,
                              const int additionnalPoints) const {
            throw std::runtime_error("Method computeK should be implemented in subclass");
        }

        virtual const T1 getAverageEdgeLength() const {
            throw std::runtime_error("Method computeSlowness should be implemented in subclass");
        }

#ifdef VTK
        virtual void saveModelVTU(const std::string &, const bool saveSlowness=true,
                                  const bool savePhysicalEntity=false) const {}
        virtual void saveModelVTR(const std::string &,
                                  const bool saveSlowness=true) const {}
        virtual void saveModelVTR(const std::string &, const double*,
                                  const bool saveSlowness=true) const {}
#endif
    protected:
        size_t nThreads;         // number of threads
        bool tt_from_rp;
        bool translateOrigin;
        sxyz<T1> origin;
        std::vector<std::vector<T2>> neighbors;  // nodes common to a cell

        template<typename N>
        void buildGridNeighbors(const std::vector<N>& nodes) {
            //Index the neighbors nodes of each cell
            for ( T2 n=0; n<nodes.size(); ++n ) {
                for ( size_t n2=0; n2<nodes[n].getOwners().size(); ++n2) {
                    neighbors[ nodes[n].getOwners()[n2] ].push_back(n);
                }
            }
        }

        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }

        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<std::vector<sxyz<T1>>>& Rx,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }

        virtual T1 getTraveltime(const sxyz<T1>& pt,
                                 const size_t threadNo) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }

        virtual T1 getTraveltimeFromRaypath(const std::vector<sxyz<T1>>& Tx,
                                            const std::vector<T1>& t0,
                                            const sxyz<T1>& Rx,
                                            const size_t threadNo) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }

        virtual void getRaypath(const std::vector<sxyz<T1>>& Tx,
                                const std::vector<T1>& t0,
                                const sxyz<T1>& Rx,
                                std::vector<sxyz<T1>>& r_data,
                                T1 &tt,
                                const size_t threadNo) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }

        virtual void getRaypath(const std::vector<sxyz<T1>>& Tx,
                                const std::vector<T1>& t0,
                                const sxyz<T1>& Rx,
                                std::vector<sijv<T1>>& m_data,
                                const size_t RxNo,
                                T1 &tt,
                                const size_t threadNo) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }

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

        virtual void getRaypath(const std::vector<sxyz<T1>>& Tx,
                                const std::vector<T1>& t0,
                                const sxyz<T1> &Rx,
                                std::vector<siv<T1>> &l_data,
                                T1 &tt,
                                const size_t threadNo) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }

        virtual void getRaypath(const std::vector<sxyz<T1>>& Tx,
                                const std::vector<T1>& t0,
                                const sxyz<T1> &Rx,
                                std::vector<sxyz<T1>> &r_data,
                                std::vector<siv<T1>> &l_data,
                                T1 &tt,
                                const size_t threadNo) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }

    private:
        bool usePool;
        mutable ctpl::thread_pool pool;

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

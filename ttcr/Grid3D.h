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
#include <thread>

#include "ttcr_t.h"

namespace ttcr {
    
    template<typename T1, typename T2>
    class Grid3D {
    public:
        
        virtual ~Grid3D() {}
        
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }
        
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                              std::vector<std::vector<T1>*>& traveltimes,
                              const size_t=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }
        
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<sxyz<T1>>>& r_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }
        
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<const std::vector<sxyz<T1>>*>& Rx,
                              std::vector<std::vector<T1>*>& traveltimes,
                              std::vector<std::vector<std::vector<sxyz<T1>>>*>& r_data,
                              const size_t=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }
        
        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<sxyz<T1>>>& r_data,
                              std::vector<std::vector<sijv<T1>>>& m_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }

        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<sijv<T1>>>& m_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }

        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<siv<T1>>>& l_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }

        virtual void raytrace(const std::vector<sxyz<T1>>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<sxyz<T1>>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<sxyz<T1>>>& r_data,
                              std::vector<std::vector<siv<T1>>>& l_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }

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

        virtual void getRaypath(const std::vector<sxyz<T1>>& Tx,
                                const std::vector<T1>& t0,
                                const sxyz<T1>& Rx,
                                std::vector<sxyz<T1>>& r_data,
                                T1 &tt,
                                const size_t threadNo=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }
        
        virtual void setSlowness(const std::vector<T1>& s) {}
        virtual void setChi(const std::vector<T1>& x) {}
        virtual void setPsi(const std::vector<T1>& x) {}
        
        virtual void setSourceRadius(const double) {}
        
        virtual size_t getNumberOfNodes() const { return 1; }
        virtual size_t getNumberOfCells() const { return 1; }
        
        virtual void saveTT(const std::string &, const int, const size_t nt=0,
                            const int format=1) const {}
        virtual void loadTT(const std::string &, const int, const size_t nt=0,
                            const int format=1) const {}

        virtual const T1 getXmin() const { return 1; }
        virtual const T1 getXmax() const { return 1; }
        virtual const T1 getYmin() const { return 1; }
        virtual const T1 getYmax() const { return 1; }
        virtual const T1 getZmin() const { return 1; }
        virtual const T1 getZmax() const { return 1; }
        
        virtual const int get_niter() const { return 0; }
        virtual const int get_niterw() const { return 0; }
        
        virtual const size_t getNthreads() const { return 1; }
        
        virtual void dump_secondary(std::ofstream&) const {}

        virtual T1 computeSlowness(const sxyz<T1>&) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }
        
#ifdef VTK
        virtual void saveModelVTU(const std::string &, const bool saveSlowness=true,
                                  const bool savePhysicalEntity=false) const {}
        virtual void saveModelVTR(const std::string &,
                                  const bool saveSlowness=true) const {}
        virtual void saveModelVTR(const std::string &, const double*,
                                  const bool saveSlowness=true) const {}
#endif
    private:
        const std::vector<size_t> get_blk_size(const size_t nTx) const {
            size_t n_blk = this->getNthreads() > nTx ? this->getNthreads() : nTx;
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
    void Grid3D<T1,T2>::raytrace(const std::vector<std::vector<sxyz<T1>>>& Tx,
                                 const std::vector<std::vector<T1>>& t0,
                                 const std::vector<std::vector<sxyz<T1>>>& Rx,
                                 std::vector<std::vector<T1>>& traveltimes) const {

        if ( Tx.size() == 0 ) {
            this->raytrace(Tx[0], t0[0], Rx[0], traveltimes[0], 0);
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

        if ( Tx.size() == 0 ) {
            this->raytrace(Tx[0], t0[0], Rx[0], traveltimes[0], r_data[0], 0);
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

        if ( Tx.size() == 0 ) {
            this->raytrace(Tx[0], t0[0], Rx[0], traveltimes[0], m_data[0], 0);
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

        if ( Tx.size() == 0 ) {
            this->raytrace(Tx[0], t0[0], Rx[0], traveltimes[0], r_data[0], m_data[0], 0);
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

        if ( Tx.size() == 0 ) {
            this->raytrace(Tx[0], t0[0], Rx[0], traveltimes[0], l_data[0], 0);
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

        if ( Tx.size() == 0 ) {
            this->raytrace(Tx[0], t0[0], Rx[0], traveltimes[0], r_data[0], l_data[0], 0);
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

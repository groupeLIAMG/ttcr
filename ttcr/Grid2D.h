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


#ifndef ttcr_Grid2D_h
#define ttcr_Grid2D_h

#include <array>
#include <fstream>
#include <functional>
#include <thread>
#include <vector>

#include "ttcr_t.h"

namespace ttcr {
    
    template<typename T1 = double, typename T2 = uint32_t, typename S = sxz<double>>
    class Grid2D {
    public:
        Grid2D(const size_t ncells, const size_t nt=1) :
            nThreads(nt), neighbors(std::vector<std::vector<T2>>(ncells)) {}

        virtual ~Grid2D() {}
        
        const size_t getNthreads() const { return nThreads; }
        
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              std::vector<T1>& traveltimes,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }
        
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<const std::vector<S>*>& Rx,
                              std::vector<std::vector<T1>*>& traveltimes,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }
        
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<S>>& r_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }
        
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<const std::vector<S>*>& Rx,
                              std::vector<std::vector<T1>*>& traveltimes,
                              std::vector<std::vector<std::vector<S>>*>& r_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }
        
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<S>>& r_data,
                              std::vector<std::vector<siv2<T1>>>& l_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }
        
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<siv2<T1>>>& l_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }
        
        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<S>>& r_data,
                              T1& v0,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }

        virtual void raytrace(const std::vector<S>& Tx,
                              const std::vector<T1>& t0,
                              const std::vector<S>& Rx,
                              std::vector<T1>& traveltimes,
                              std::vector<std::vector<S>>& r_data,
                              T1& v0,
                              std::vector<std::vector<sijv<T1>>>& m_data,
                              const size_t threadNo=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }

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
                      std::vector<std::vector<std::vector<S>>>& r_data,
                      std::vector<std::vector<std::vector<siv2<T1>>>>& l_data) const;

        virtual void setSlowness(const std::vector<T1>& s) {}
        virtual void setXi(const std::vector<T1>& x) {
                throw std::runtime_error("Method should be implemented in subclass");
            }
        virtual void setTiltAngle(const std::vector<T1>& x) {
                throw std::runtime_error("Method should be implemented in subclass");
            }
        virtual void setVp0(const std::vector<T1>& s) {
                throw std::runtime_error("Method should be implemented in subclass");
            }
        virtual void setVs0(const std::vector<T1>& s) {
                throw std::runtime_error("Method should be implemented in subclass");
            }
        virtual void setDelta(const std::vector<T1>& s) {
                throw std::runtime_error("Method should be implemented in subclass");
            }
        virtual void setEpsilon(const std::vector<T1>& s) {
                throw std::runtime_error("Method should be implemented in subclass");
            }
        virtual void setGamma(const std::vector<T1>& s) {
                throw std::runtime_error("Method should be implemented in subclass");
            }
        
        virtual size_t getNumberOfNodes(const bool primary=false) const { return 1; }
        virtual size_t getNumberOfCells() const { return 1; }
        virtual void getTT(std::vector<T1>& tt, const size_t threadNo=0) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }

        virtual void saveTT(const std::string &, const int, const size_t nt=0,
                            const int format=1) const {}
        
        virtual void saveTTgrad(const std::string &, const size_t nt=0,
                                const bool vtkFormat=0) const {}
        
        virtual const T1 getXmin() const { return 1; }
        virtual const T1 getXmax() const { return 1; }
        virtual const T1 getYmin() const { return 1; }
        virtual const T1 getYmax() const { return 1; }
        virtual const T1 getZmin() const { return 1; }
        virtual const T1 getZmax() const { return 1; }
        virtual const T1 getDx() const { return 1; }
        virtual const T1 getDz() const { return 1; }
        virtual const T2 getNcx() const { return 1; }
        virtual const T2 getNcz() const { return 1; }
        
        
        virtual const int get_niter() const { return 0; }
        virtual const int get_niterw() const { return 0; }
        
        virtual int projectPts(std::vector<S>&) const { return 1; }
        
        virtual void interpolateAtNodes(std::vector<T1> &) const {}
        virtual void interpolateAtNodes(const std::vector<T1> &,
                                        std::vector<T1> &) const {}
        
        virtual void dump_secondary(std::ofstream&) const {};
        
        virtual void getNodes(std::vector<S>& nodes) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }
        virtual void getTriangles(std::vector<std::array<T2, 3>>&) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }
        // keep next method until array is supported by cython
        virtual void getTriangles(std::vector<std::vector<T2>>&) const {
            throw std::runtime_error("Method should be implemented in subclass");
        }

#ifdef VTK
        virtual void saveModelVTU(const std::string &, const bool saveSlowness=true,
                                  const bool savePhysicalEntity=false) const {}
        virtual void saveModelVTR(const std::string &, const double*,
                                  const bool saveSlowness=true) const {}
#endif
    protected:
        size_t nThreads;
        
        std::vector<std::vector<T2>> neighbors;  // nodes common to a cell
        
        template<typename N>
        void buildGridNeighbors(std::vector<N>& nodes) {
            for ( T2 n=0; n<nodes.size(); ++n ) {
                for ( size_t n2=0; n2<nodes[n].getOwners().size(); ++n2) {
                    neighbors[ nodes[n].getOwners()[n2] ].push_back(n);
                }
            }
        }

    private:
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
    void Grid2D<T1,T2,S>::raytrace(const std::vector<std::vector<S>>& Tx,
                                   const std::vector<std::vector<T1>>& t0,
                                   const std::vector<std::vector<S>>& Rx,
                                   std::vector<std::vector<T1>>& traveltimes) const {
        
        if ( Tx.size() == 1 ) {
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

    template<typename T1, typename T2, typename S>
    void Grid2D<T1,T2,S>::raytrace(const std::vector<std::vector<S>>& Tx,
                                   const std::vector<std::vector<T1>>& t0,
                                   const std::vector<std::vector<S>>& Rx,
                                   std::vector<std::vector<T1>>& traveltimes,
                                   std::vector<std::vector<std::vector<S>>>& r_data) const {

        if ( Tx.size() == 1 ) {
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

    template<typename T1, typename T2, typename S>
    void Grid2D<T1,T2,S>::raytrace(const std::vector<std::vector<S>>& Tx,
                                   const std::vector<std::vector<T1>>& t0,
                                   const std::vector<std::vector<S>>& Rx,
                                   std::vector<std::vector<T1>>& traveltimes,
                                   std::vector<std::vector<std::vector<siv2<T1>>>>& l_data) const {
        
        if ( Tx.size() == 1 ) {
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

    template<typename T1, typename T2, typename S>
    void Grid2D<T1,T2,S>::raytrace(const std::vector<std::vector<S>>& Tx,
                                   const std::vector<std::vector<T1>>& t0,
                                   const std::vector<std::vector<S>>& Rx,
                                   std::vector<std::vector<T1>>& traveltimes,
                                   std::vector<std::vector<std::vector<S>>>& r_data,
                                   std::vector<std::vector<std::vector<siv2<T1>>>>& l_data) const {
        
        if ( Tx.size() == 1 ) {
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

//
//  Grid3Dttcr.cpp
//  ttcr
//
//  Created by Bernard Giroux on 16-10-11.
//  Copyright © 2016 Bernard Giroux. All rights reserved.
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

#include <thread>

#include "Grid3Dttcr.h"

using namespace std;


namespace ttcr {

    Grid3Dttcr::Grid3Dttcr(std::string& t,
                           const uint32_t nx, const uint32_t ny, const uint32_t nz,
                           const double ddx,
                           const double minx, const double miny, const double minz,
                           const double eps, const int maxit, const bool w,
                           const size_t nt=1) : type(t) {
        if ( type.compare("node") == 0 )
          grid_instance = new gridi(nx, ny, nz, ddx, minx, miny, minz, eps, maxit, w, nt);
        else if ( type.compare("cell") == 0 )
          grid_instance = new gridc(nx, ny, nz, ddx, minx, miny, minz, eps, maxit, w, nt);
        else
          // error: type not defined
          throw bad_cast();
      }

    void Grid3Dttcr::setSlowness(const std::vector<double>& slowness) {
        if ( grid_instance->setSlowness(slowness) == 1 ) {
            throw out_of_range("Number of slowness values not consistent with grid.");
        }
    }

    int Grid3Dttcr::raytrace(const std::vector<sxyz<double>>& Tx,
                             const std::vector<double>& tTx,
                             const std::vector<sxyz<double>>& Rx,
                             double* traveltimes) const {
        /*
         Looking for redundants Tx pts
         */

        size_t nTx = Tx.size();
        size_t nRx = Rx.size();
        vector<vector<sxyz<double>>> vTx;
        vector<vector<double>> t0;
        vector<vector<size_t>> iTx;
        vTx.push_back( vector<sxyz<double> >(1, Tx[0]) );
        t0.push_back( vector<double>(1, tTx[0]) );
        iTx.push_back( vector<size_t>(1, 0) );  // indices of Rx corresponding to current Tx
        for ( size_t ntx=1; ntx<nTx; ++ntx ) {
            bool found = false;

            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                if ( vTx[nv][0]==Tx[ntx] ) {
                    found = true;
                    iTx[nv].push_back( ntx ) ;
                    break;
                }
            }
            if ( !found ) {
                vTx.push_back( vector<sxyz<double>>(1, Tx[ntx]) );
                t0.push_back( vector<double>(1, tTx[ntx]) );
                iTx.push_back( vector<size_t>(1, ntx) );
            }
        }

        /*
         Looping over all non redundant Tx
         */

        vector<sxyz<double>> vRx;
        vector<vector<double>> tt( vTx.size() );

        if ( grid_instance->getNthreads() == 1 || vTx.size() <= grid_instance->getNthreads() ) {
            for ( size_t nv=0; nv<vTx.size(); ++nv ) {

                vRx.resize( 0 );
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ] );
                }
                if ( grid_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv]) == 1 ) {
                    throw runtime_error("Problem while raytracing.");
                }
            }
        } else {
            size_t num_threads = grid_instance->getNthreads();
            size_t blk_size = vTx.size()/num_threads;
            if ( blk_size == 0 ) blk_size++;

            vector<thread> threads(num_threads-1);
            size_t blk_start = 0;
            for ( size_t i=0; i<num_threads-1; ++i ) {

                size_t blk_end = blk_start + blk_size;
                grid *grid_ref = grid_instance;
                threads[i]=thread( [&grid_ref,&vTx,&tt,&t0,&Rx,&iTx,&nRx,
                                    blk_start,blk_end,i]{

                    for ( size_t nv=blk_start; nv<blk_end; ++nv ) {

                        vector<sxyz<double>> vRx;
                        for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                            vRx.push_back( Rx[ iTx[nv][ni] ] );
                        }
                        if ( grid_ref->raytrace(vTx[nv], t0[nv], vRx, tt[nv], i+1) == 1 ) {
                            throw runtime_error("Problem while raytracing.");
                        }
                    }
                });

                blk_start = blk_end;
            }

            for ( size_t nv=blk_start; nv<vTx.size(); ++nv ) {
                vector<sxyz<double>> vRx;
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ] );
                }
                if ( grid_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], 0) == 1 ) {
                    throw runtime_error("Problem while raytracing.");
                }
            }

            std::for_each(threads.begin(),threads.end(),
                          std::mem_fn(&std::thread::join));
        }

        for ( size_t nv=0; nv<vTx.size(); ++nv ) {
            for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                traveltimes[ iTx[nv][ni] ] = tt[nv][ni];
            }
        }
        return 0;
    }

    int Grid3Dttcr::raytrace(const std::vector<sxyz<double>>& Tx,
                             const std::vector<double>& tTx,
                             const std::vector<sxyz<double>>& Rx,
                             double* traveltimes,
                             PyObject* rays,
                             double* V0) const {
        // rays must be a pointer to a tuple object of size nRx

        /*
         Looking for redundants Tx pts
         */

        size_t nTx = Tx.size();
        size_t nRx = Rx.size();
        vector<vector<sxyz<double>>> vTx;
        vector<vector<double>> t0;
        vector<vector<size_t>> iTx;
        vTx.push_back( vector<sxyz<double> >(1, Tx[0]) );
        t0.push_back( vector<double>(1, tTx[0]) );
        iTx.push_back( vector<size_t>(1, 0) );  // indices of Rx corresponding to current Tx
        for ( size_t ntx=1; ntx<nTx; ++ntx ) {
            bool found = false;

            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                if ( vTx[nv][0]==Tx[ntx] ) {
                    found = true;
                    iTx[nv].push_back( ntx ) ;
                    break;
                }
            }
            if ( !found ) {
                vTx.push_back( vector<sxyz<double>>(1, Tx[ntx]) );
                t0.push_back( vector<double>(1, tTx[ntx]) );
                iTx.push_back( vector<size_t>(1, ntx) );
            }
        }

        /*
         Looping over all non redundant Tx
         */

        vector<sxyz<double>> vRx;
        vector<vector<double>> tt( vTx.size() );
        vector<vector<vector<sxyz<double>>>> r_data( vTx.size() );
        vector<double> v0( vTx.size() );

        if ( grid_instance->getNthreads() == 1 || vTx.size() <= grid_instance->getNthreads() ) {
            for ( size_t nv=0; nv<vTx.size(); ++nv ) {

                vRx.resize( 0 );
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ] );
                }

                if ( grid_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv]) == 1 ) {
                    throw runtime_error("Problem while raytracing.");
                }
            }
        } else {
            size_t num_threads = grid_instance->getNthreads();
            size_t blk_size = vTx.size()/num_threads;
            if ( blk_size == 0 ) blk_size++;

            vector<thread> threads(num_threads-1);
            size_t blk_start = 0;
            for ( size_t i=0; i<num_threads-1; ++i ) {

                size_t blk_end = blk_start + blk_size;
                grid *grid_ref = grid_instance;
                threads[i]=thread( [&grid_ref,&vTx,&tt,&t0,&Rx,&iTx,&nRx,
                                    &r_data,&v0,blk_start,blk_end,i]{

                    for ( size_t nv=blk_start; nv<blk_end; ++nv ) {

                        vector<sxyz<double>> vRx;
                        for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                            vRx.push_back( Rx[ iTx[nv][ni] ] );
                        }
                        if ( grid_ref->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], i+1) == 1 ) {
                            throw runtime_error("Problem while raytracing.");
                        }
                    }
                });

                blk_start = blk_end;
            }

            for ( size_t nv=blk_start; nv<vTx.size(); ++nv ) {
                vector<sxyz<double>> vRx;
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ] );
                }
                if ( grid_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], 0) == 1 ) {
                    throw runtime_error("Problem while raytracing.");
                }
            }

            std::for_each(threads.begin(),threads.end(),
                          std::mem_fn(&std::thread::join));
        }

        for ( size_t nv=0; nv<vTx.size(); ++nv ) {
            for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                traveltimes[ iTx[nv][ni] ] = tt[nv][ni];
                V0[ iTx[nv][ni] ] = v0[nv];
            }
        }

        // rays
        import_array();  // to use PyArray_SimpleNewFromData

        for ( size_t nv=0; nv<vTx.size(); ++nv ) {
            for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                size_t npts = r_data[nv][ni].size();
                npy_intp dims[] = {static_cast<npy_intp>(npts), 3};
                double* ray_p = new double[3*npts];
                PyObject* ray = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, ray_p);

                for ( size_t np=0; np<npts; ++np ) {
                    ray_p[3*np] = r_data[nv][ni][np].x;
                    ray_p[3*np+1] = r_data[nv][ni][np].y;
                    ray_p[3*np+2] = r_data[nv][ni][np].z;
                }

                PyTuple_SetItem(rays, iTx[nv][ni], ray);
            }
        }
        return 0;
    }

    int Grid3Dttcr::raytrace(const std::vector<sxyz<double>>& Tx,
                             const std::vector<double>& tTx,
                             const std::vector<sxyz<double>>& Rx,
                             double* traveltimes,
                             PyObject* rays,
                             double* V0,
                             PyObject* M) const {
        // rays must be a pointer to a tuple object of size nRx

        /*
         Looking for redundants Tx pts
         */

        size_t nTx = Tx.size();
        size_t nRx = Rx.size();
        vector<vector<sxyz<double>>> vTx;
        vector<vector<double>> t0;
        vector<vector<size_t>> iTx;
        vTx.push_back( vector<sxyz<double> >(1, Tx[0]) );
        t0.push_back( vector<double>(1, tTx[0]) );
        iTx.push_back( vector<size_t>(1, 0) );  // indices of Rx corresponding to current Tx
        for ( size_t ntx=1; ntx<nTx; ++ntx ) {
            bool found = false;

            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                if ( vTx[nv][0]==Tx[ntx] ) {
                    found = true;
                    iTx[nv].push_back( ntx ) ;
                    break;
                }
            }
            if ( !found ) {
                vTx.push_back( vector<sxyz<double>>(1, Tx[ntx]) );
                t0.push_back( vector<double>(1, tTx[ntx]) );
                iTx.push_back( vector<size_t>(1, ntx) );
            }
        }

        /*
         Looping over all non redundant Tx
         */

        vector<sxyz<double>> vRx;
        vector<vector<double>> tt( vTx.size() );
        vector<vector<vector<sxyz<double>>>> r_data( vTx.size() );
        vector<double> v0( vTx.size() );
        vector<vector<vector<sijv<double>>>> m_data( vTx.size() );

        if ( grid_instance->getNthreads() == 1 || vTx.size() <= grid_instance->getNthreads() ) {
            for ( size_t nv=0; nv<vTx.size(); ++nv ) {

                vRx.resize( 0 );
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ] );
                }

                if ( grid_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], m_data[nv]) == 1 ) {
                    throw runtime_error("Problem while raytracing.");
                }
            }
        } else {
            size_t num_threads = grid_instance->getNthreads();
            size_t blk_size = vTx.size()/num_threads;
            if ( blk_size == 0 ) blk_size++;

            vector<thread> threads(num_threads-1);
            size_t blk_start = 0;
            for ( size_t i=0; i<num_threads-1; ++i ) {

                size_t blk_end = blk_start + blk_size;
                grid *grid_ref = grid_instance;
                threads[i]=thread( [&grid_ref,&vTx,&tt,&t0,&Rx,&iTx,&nRx,
                                    &r_data,&v0,&m_data,blk_start,blk_end,i]{

                    for ( size_t nv=blk_start; nv<blk_end; ++nv ) {

                        vector<sxyz<double>> vRx;
                        for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                            vRx.push_back( Rx[ iTx[nv][ni] ] );
                        }
                        if ( grid_ref->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], m_data[nv], i+1) == 1 ) {
                            throw runtime_error("Problem while raytracing.");
                        }
                    }
                });

                blk_start = blk_end;
            }

            for ( size_t nv=blk_start; nv<vTx.size(); ++nv ) {
                vector<sxyz<double>> vRx;
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ] );
                }
                if ( grid_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], m_data[nv], 0) == 1 ) {
                    throw runtime_error("Problem while raytracing.");
                }
            }

            std::for_each(threads.begin(),threads.end(),
                          std::mem_fn(&std::thread::join));
        }

        for ( size_t nv=0; nv<vTx.size(); ++nv ) {
            for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                traveltimes[ iTx[nv][ni] ] = tt[nv][ni];
                V0[ iTx[nv][ni] ] = v0[nv];
            }
        }

        // rays
        import_array();  // to use PyArray_SimpleNewFromData

        for ( size_t nv=0; nv<vTx.size(); ++nv ) {
            for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                size_t npts = r_data[nv][ni].size();
                npy_intp dims[] = {static_cast<npy_intp>(npts), 3};
                double* ray_p = new double[3*npts];
                PyObject* ray = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, ray_p);

                for ( size_t np=0; np<npts; ++np ) {
                    ray_p[3*np] = r_data[nv][ni][np].x;
                    ray_p[3*np+1] = r_data[nv][ni][np].y;
                    ray_p[3*np+2] = r_data[nv][ni][np].z;
                }

                PyTuple_SetItem(rays, iTx[nv][ni], ray);
            }
        }

        // data for matrix M
        // M is a tuple of tuples
        // first element of tuple contains data, size is nnz
        // second element contains column indices for rows of the matrix, size is nnz
        // third element contains pointers for indices, size is nrow+1 (nTx+1)

        size_t nnodes = grid_instance->getNumberOfNodes();

        for ( size_t nv=0; nv<vTx.size(); ++nv ) {

            PyObject* tuple = PyTuple_New(3);

            size_t nRcv = m_data[nv].size();
            size_t nnz = 0;
            for ( size_t ni=0; ni<m_data[nv].size(); ++ni ) {
                nnz += m_data[nv][ni].size();
            }


            npy_intp dims[] = {static_cast<npy_intp>(nnz)};
            double* data_p = new double[nnz];
            PyObject* data = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, data_p);

            int64_t* indices_p = new int64_t[nnz];
            PyObject* indices = PyArray_SimpleNewFromData(1, dims, NPY_INT64, indices_p);

            dims[0] = nRcv+1;
            int64_t* indptr_p = new int64_t[nRcv+1];
            PyObject* indptr = PyArray_SimpleNewFromData(1, dims, NPY_INT64, indptr_p);

            size_t k = 0;
            for ( size_t ni=0; ni<nRcv; ++ni ) {
                indptr_p[ni] = k;
                for ( size_t j=0; j<nnodes; ++j ) {
                    for ( size_t n=0; n<m_data[nv][ni].size(); ++n) {
                        if ( m_data[nv][ni][n].j == j && m_data[nv][ni][n].i == ni ) {
                            indices_p[k] = j;
                            data_p[k] = m_data[nv][ni][n].v;
                            k++;
                        }
                    }
                }
            }
            indptr_p[nRcv] = k;

            PyTuple_SetItem(tuple, 0, data);
            PyTuple_SetItem(tuple, 1, indices);
            PyTuple_SetItem(tuple, 2, indptr);

            PyTuple_SetItem(M, nv, tuple);
        }
        return 0;
    }
    
    
    int Grid3Dttcr::raytrace(const std::vector<sxyz<double>>& Tx,
                             const std::vector<double>& tTx,
                             const std::vector<sxyz<double>>& Rx,
                             double* traveltimes,
                             PyObject* rays,
                             PyObject* L) const {
        
        // rays must be a pointer to a tuple object of size nRx
        
        /*
         Looking for redundants Tx pts
         */
        
        size_t nTx = Tx.size();
        size_t nRx = Rx.size();
        vector<vector<sxyz<double>>> vTx;
        vector<vector<double>> t0;
        vector<vector<size_t>> iTx;
        vTx.push_back( vector<sxyz<double> >(1, Tx[0]) );
        t0.push_back( vector<double>(1, tTx[0]) );
        iTx.push_back( vector<size_t>(1, 0) );  // indices of Rx corresponding to current Tx
        for ( size_t ntx=1; ntx<nTx; ++ntx ) {
            bool found = false;
            
            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                if ( vTx[nv][0]==Tx[ntx] ) {
                    found = true;
                    iTx[nv].push_back( ntx ) ;
                    break;
                }
            }
            if ( !found ) {
                vTx.push_back( vector<sxyz<double>>(1, Tx[ntx]) );
                t0.push_back( vector<double>(1, tTx[ntx]) );
                iTx.push_back( vector<size_t>(1, ntx) );
            }
        }
        
        /*
         Looping over all non redundant Tx
         */
        
        vector<sxyz<double>> vRx;
        vector<vector<double>> tt( vTx.size() );
        vector<vector<vector<sxyz<double>>>> r_data( vTx.size() );
        vector<vector<siv<double>>> L_data(nTx);
        vector<vector<vector<siv<double>>>> l_data( vTx.size() );
        
        if ( grid_instance->getNthreads() == 1 || vTx.size() <= grid_instance->getNthreads() ) {
            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                
                vRx.resize( 0 );
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ] );
                }
                
                if ( grid_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], l_data[nv]) == 1 ) {
                    throw runtime_error("Problem while raytracing.");
                }
            }
        } else {
            size_t num_threads = grid_instance->getNthreads();
            size_t blk_size = vTx.size()/num_threads;
            if ( blk_size == 0 ) blk_size++;
            
            vector<thread> threads(num_threads-1);
            size_t blk_start = 0;
            for ( size_t i=0; i<num_threads-1; ++i ) {
                
                size_t blk_end = blk_start + blk_size;
                grid *grid_ref = grid_instance;
                threads[i]=thread( [&grid_ref,&vTx,&tt,&t0,&Rx,&iTx,&nRx,
                                    &r_data,&l_data,blk_start,blk_end,i]{
                    
                    for ( size_t nv=blk_start; nv<blk_end; ++nv ) {
                        
                        vector<sxyz<double>> vRx;
                        for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                            vRx.push_back( Rx[ iTx[nv][ni] ] );
                        }
                        if ( grid_ref->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], l_data[nv], i+1) == 1 ) {
                            throw runtime_error("Problem while raytracing.");
                        }
                    }
                });
                
                blk_start = blk_end;
            }
            
            for ( size_t nv=blk_start; nv<vTx.size(); ++nv ) {
                vector<sxyz<double>> vRx;
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ] );
                }
                if ( grid_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], l_data[nv], 0) == 1 ) {
                    throw runtime_error("Problem while raytracing.");
                }
            }
            
            std::for_each(threads.begin(),threads.end(),
                          std::mem_fn(&std::thread::join));
        }
        
        for ( size_t nv=0; nv<vTx.size(); ++nv ) {
            for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                traveltimes[ iTx[nv][ni] ] = tt[nv][ni];
            }
        }
        
        // rays
        import_array();  // to use PyArray_SimpleNewFromData
        
        for ( size_t nv=0; nv<vTx.size(); ++nv ) {
            for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                size_t npts = r_data[nv][ni].size();
                npy_intp dims[] = {static_cast<npy_intp>(npts), 3};
                double* ray_p = new double[3*npts];
                PyObject* ray = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, ray_p);
                PyArray_ENABLEFLAGS((PyArrayObject*)ray, NPY_ARRAY_OWNDATA);
                
                for ( size_t np=0; np<npts; ++np ) {
                    ray_p[3*np] = r_data[nv][ni][np].x;
                    ray_p[3*np+1] = r_data[nv][ni][np].y;
                    ray_p[3*np+2] = r_data[nv][ni][np].z;
                }
                
                PyTuple_SetItem(rays, iTx[nv][ni], ray);
            }
        }
        
        // L
        // first element of tuple contains data, size is nnz
        // second element contains column indices for rows of the matrix, size is nnz
        // third element contains pointers for indices, size is nrow+1 (nTx+1)
        
        for ( size_t nv=0; nv<vTx.size(); ++nv ) {
            for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                L_data[ iTx[nv][ni] ] = l_data[nv][ni];
            }
        }
        
        size_t nnz = 0;
        for ( size_t n=0; n<L_data.size(); ++n ) {
            nnz += L_data[n].size();
        }
        size_t ncell = grid_instance->getNumberOfCells();
        
        npy_intp dims[] = {static_cast<npy_intp>(nnz)};
        double* data_p = new double[nnz];
        PyObject* data = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, data_p);
        PyArray_ENABLEFLAGS((PyArrayObject*)data, NPY_ARRAY_OWNDATA);
        
        int64_t* indices_p = new int64_t[nnz];
        PyObject* indices = PyArray_SimpleNewFromData(1, dims, NPY_INT64, indices_p);
        PyArray_ENABLEFLAGS((PyArrayObject*)indices, NPY_ARRAY_OWNDATA);
        
        dims[0] = nTx+1;
        int64_t* indptr_p = new int64_t[nTx+1];
        PyObject* indptr = PyArray_SimpleNewFromData(1, dims, NPY_INT64, indptr_p);
        PyArray_ENABLEFLAGS((PyArrayObject*)indptr, NPY_ARRAY_OWNDATA);
        
        size_t k = 0;
        for ( size_t i=0; i<nTx; ++i ) {
            indptr_p[i] = k;
            for ( size_t j=0; j<ncell; ++j ) {
                for ( size_t n=0; n<L_data[i].size(); ++n ) {
                    if ( L_data[i][n].i == j ) {
                        indices_p[k] = j;
                        data_p[k] = L_data[i][n].v;
                        k++;
                    }
                }
            }
        }
        indptr_p[nTx] = k;
        
        PyTuple_SetItem(L, 0, data);
        PyTuple_SetItem(L, 1, indices);
        PyTuple_SetItem(L, 2, indptr);
        
        return 0;
    }
    
    int Grid3Dttcr::raytrace(const std::vector<sxyz<double>>& Tx,
                             const std::vector<double>& tTx,
                             const std::vector<sxyz<double>>& Rx,
                             double* traveltimes,
                             PyObject* L) const {
        
        /*
         Looking for redundants Tx pts
         */
        
        size_t nTx = Tx.size();
        size_t nRx = Rx.size();
        vector<vector<sxyz<double>>> vTx;
        vector<vector<double>> t0;
        vector<vector<size_t>> iTx;
        vTx.push_back( vector<sxyz<double> >(1, Tx[0]) );
        t0.push_back( vector<double>(1, tTx[0]) );
        iTx.push_back( vector<size_t>(1, 0) );  // indices of Rx corresponding to current Tx
        for ( size_t ntx=1; ntx<nTx; ++ntx ) {
            bool found = false;
            
            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                if ( vTx[nv][0]==Tx[ntx] ) {
                    found = true;
                    iTx[nv].push_back( ntx ) ;
                    break;
                }
            }
            if ( !found ) {
                vTx.push_back( vector<sxyz<double>>(1, Tx[ntx]) );
                t0.push_back( vector<double>(1, tTx[ntx]) );
                iTx.push_back( vector<size_t>(1, ntx) );
            }
        }
        
        /*
         Looping over all non redundant Tx
         */
        
        vector<sxyz<double>> vRx;
        vector<vector<double>> tt( vTx.size() );
        vector<vector<siv<double>>> L_data(nTx);
        vector<vector<vector<siv<double>>>> l_data( vTx.size() );
        
        if ( grid_instance->getNthreads() == 1 || vTx.size()<=grid_instance->getNthreads() ) {
            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                
                vRx.resize( 0 );
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ] );
                }
                
                if ( grid_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], l_data[nv]) == 1 ) {
                    throw runtime_error("Problem while raytracing.");
                }
            }
        } else {
            size_t num_threads = grid_instance->getNthreads() < vTx.size() ? grid_instance->getNthreads() : vTx.size();
            size_t blk_size = vTx.size()/num_threads;
            if ( blk_size == 0 ) blk_size++;
            
            vector<thread> threads(num_threads-1);
            size_t blk_start = 0;
            for ( size_t i=0; i<num_threads-1; ++i ) {
                
                size_t blk_end = blk_start + blk_size;
                grid *grid_ref = grid_instance;
                threads[i]=thread( [&grid_ref,&vTx,&tt,&t0,&Rx,&iTx,&nRx,
                                    &l_data,blk_start,blk_end,i]{
                    
                    for ( size_t nv=blk_start; nv<blk_end; ++nv ) {
                        
                        vector<sxyz<double>> vRx;
                        for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                            vRx.push_back( Rx[ iTx[nv][ni] ] );
                        }
                        if ( grid_ref->raytrace(vTx[nv], t0[nv], vRx, tt[nv], l_data[nv], i+1) == 1 ) {
                            throw runtime_error("Problem while raytracing.");
                        }
                    }
                });
                
                blk_start = blk_end;
            }
            
            for ( size_t nv=blk_start; nv<vTx.size(); ++nv ) {
                vector<sxyz<double>> vRx;
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ] );
                }
                if ( grid_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], l_data[nv], 0) == 1 ) {
                    throw runtime_error("Problem while raytracing.");
                }
            }
            
            std::for_each(threads.begin(),threads.end(),
                          std::mem_fn(&std::thread::join));
        }
        
        for ( size_t nv=0; nv<vTx.size(); ++nv ) {
            for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                traveltimes[ iTx[nv][ni] ] = tt[nv][ni];
            }
        }
        
        import_array();  // to use PyArray_SimpleNewFromData
        
        // L
        // first element of tuple contains data, size is nnz
        // second element contains column indices for rows of the matrix, size is nnz
        // third element contains pointers for indices, size is nrow+1 (nTx+1)
        
        for ( size_t nv=0; nv<vTx.size(); ++nv ) {
            for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                L_data[ iTx[nv][ni] ] = l_data[nv][ni];
            }
        }
        
        size_t nnz = 0;
        for ( size_t n=0; n<L_data.size(); ++n ) {
            nnz += L_data[n].size();
        }
        size_t ncell = grid_instance->getNumberOfCells();
        
        npy_intp dims[] = {static_cast<npy_intp>(nnz)};
        double* data_p = new double[nnz];
        PyObject* data = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, data_p);
        PyArray_ENABLEFLAGS((PyArrayObject*)data, NPY_ARRAY_OWNDATA);
        
        int64_t* indices_p = new int64_t[nnz];
        PyObject* indices = PyArray_SimpleNewFromData(1, dims, NPY_INT64, indices_p);
        PyArray_ENABLEFLAGS((PyArrayObject*)indices, NPY_ARRAY_OWNDATA);
        
        dims[0] = nTx+1;
        int64_t* indptr_p = new int64_t[nTx+1];
        PyObject* indptr = PyArray_SimpleNewFromData(1, dims, NPY_INT64, indptr_p);
        PyArray_ENABLEFLAGS((PyArrayObject*)indptr, NPY_ARRAY_OWNDATA);
        
        size_t k = 0;
        for ( size_t i=0; i<nTx; ++i ) {
            indptr_p[i] = k;
            for ( size_t j=0; j<ncell; ++j ) {
                for ( size_t n=0; n<L_data[i].size(); ++n ) {
                    if ( L_data[i][n].i == j ) {
                        indices_p[k] = j;
                        data_p[k] = L_data[i][n].v;
                        k++;
                    }
                }
            }
        }
        indptr_p[nTx] = k;
        
        PyTuple_SetItem(L, 0, data);
        PyTuple_SetItem(L, 1, indices);
        PyTuple_SetItem(L, 2, indptr);
        
        
        return 0;
    }


}

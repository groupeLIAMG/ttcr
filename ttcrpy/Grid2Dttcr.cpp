//
//  Grid2Dttcr.cpp
//  ttcr
//
//  Created by Bernard Giroux on 16-07-28.
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <functional>
#include <thread>

#include "Grid2Dttcr.h"

using namespace std;

#define import_array_throw() {if (_import_array() < 0) {PyErr_Print(); PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import"); throw NUMPY_IMPORT_ARRAY_RETVAL; } }

namespace ttcr {

    Grid2Dttcr::Grid2Dttcr(std::string& _type,
                           uint32_t nx, uint32_t nz,
                           double dx, double dz,
                           double xmin, double zmin,
                           uint32_t nsnx, uint32_t nsnz,
                           size_t nthreads) : type(_type){

        if ( type.compare("iso")==0 ) {
            grid_instance = new gridiso(nx, nz,
                                        dx, dz,
                                        xmin, zmin,
                                        nsnx, nsnz,
                                        nthreads);
        } else if ( type.compare("elliptical")==0 ) {
            grid_instance = new gridaniso(nx, nz,
                                          dx, dz,
                                          xmin, zmin,
                                          nsnx, nsnz,
                                          nthreads);
        } else if ( type.compare("tilted")==0 ) {
            grid_instance = new gridtilted(nx, nz,
                                           dx, dz,
                                           xmin, zmin,
                                           nsnx, nsnz,
                                           nthreads);
        } else {
            // error: type not defined
            throw bad_cast();
        }
    }

    void Grid2Dttcr::setSlowness(const std::vector<double>& slowness) {
        try {
            grid_instance->setSlowness(slowness);
        } catch (length_error& e) {
            throw;
        }
    }

    void Grid2Dttcr::setXi(const std::vector<double>& xi) {
        try {
            grid_instance->setXi(xi);
        } catch (std::exception& e) {
            throw;
        }
    }

    void Grid2Dttcr::setTheta(const std::vector<double>& theta) {
        try {
            grid_instance->setTiltAngle(theta);
        } catch (std::exception& e) {
            throw;
        }
    }

    void Grid2Dttcr::raytrace(const std::vector<sxz<double>>& Tx,
                              const std::vector<double>& tTx,
                              const std::vector<sxz<double>>& Rx,
                              double* traveltimes,
                              PyObject* rays,
                              PyObject* L) const {

        // rays must be a pointer to a tuple object of size nRx

        /*
         Looking for redundants Tx pts
         */

        size_t nTx = Tx.size();
        vector<vector<sxz<double>>> vTx;
        vector<vector<double>> t0;
        vector<vector<size_t>> iTx;
        vTx.push_back( vector<sxz<double> >(1, Tx[0]) );
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
                vTx.push_back( vector<sxz<double>>(1, Tx[ntx]) );
                t0.push_back( vector<double>(1, tTx[ntx]) );
                iTx.push_back( vector<size_t>(1, ntx) );
            }
        }

        /*
         Looping over all non redundant Tx
         */

        vector<sxz<double>> vRx;
        vector<vector<double>> tt( vTx.size() );
        vector<vector<vector<sxz<double>>>> r_data( vTx.size() );
        vector<vector<siv2<double>>> L_data(nTx);
        vector<vector<vector<siv2<double>>>> l_data( vTx.size() );

        if ( grid_instance->getNthreads() == 1 || vTx.size() <= grid_instance->getNthreads() ) {
            for ( size_t nv=0; nv<vTx.size(); ++nv ) {

                vRx.resize( 0 );
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ] );
                }

                try {
                    grid_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], l_data[nv]);
                } catch (std::exception& e) {
                    throw;
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
                threads[i]=thread( [&grid_ref,&vTx,&tt,&t0,&Rx,&iTx,
                                    &r_data,&l_data,blk_start,blk_end,i]{

                    for ( size_t nv=blk_start; nv<blk_end; ++nv ) {

                        vector<sxz<double>> vRx;
                        for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                            vRx.push_back( Rx[ iTx[nv][ni] ] );
                        }
                        try {
                            grid_ref->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], l_data[nv], i+1);
                        } catch (std::exception& e) {
                            throw;
                        }
                    }
                });

                blk_start = blk_end;
            }

            for ( size_t nv=blk_start; nv<vTx.size(); ++nv ) {
                vector<sxz<double>> vRx;
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ] );
                }
                try {
                    grid_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], l_data[nv], 0);
                } catch (std::exception& e) {
                    throw;
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
        import_array_throw();  // to use PyArray_SimpleNewFromData

        for ( size_t nv=0; nv<vTx.size(); ++nv ) {
            for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                size_t npts = r_data[nv][ni].size();
                npy_intp dims[] = {static_cast<npy_intp>(npts), 2};
                double* ray_p = new double[2*npts];
                PyObject* ray = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, ray_p);
                PyArray_ENABLEFLAGS((PyArrayObject*)ray, NPY_ARRAY_OWNDATA);

                for ( size_t np=0; np<npts; ++np ) {
                    ray_p[2*np] = r_data[nv][ni][np].x;
                    ray_p[2*np+1] = r_data[nv][ni][np].z;
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

        if ( type.compare("iso")==0 ) {
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

        } else {
            nnz *= 2;

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
                for ( size_t j=0; j<2*ncell; ++j ) {
                    for ( size_t n=0; n<L_data[i].size(); ++n ) {
                        if ( L_data[i][n].i == j ) {
                            indices_p[k] = j;
                            data_p[k] = L_data[i][n].v;
                            k++;
                        }	else if ( L_data[i][n].i+ncell == j ) {
                            indices_p[k] = j;
                            data_p[k] = L_data[i][n].v2;
                            k++;
                        }
                    }
                }
            }
            indptr_p[nTx] = k;

            PyTuple_SetItem(L, 0, data);
            PyTuple_SetItem(L, 1, indices);
            PyTuple_SetItem(L, 2, indptr);
        }
    }

    void Grid2Dttcr::raytrace(const std::vector<sxz<double>>& Tx,
                              const std::vector<double>& tTx,
                              const std::vector<sxz<double>>& Rx,
                              double* traveltimes) const {

        /*
         Looking for redundants Tx pts
         */

        size_t nTx = Tx.size();
        vector<vector<sxz<double>>> vTx;
        vector<vector<double>> t0;
        vector<vector<size_t>> iTx;
        vTx.push_back( vector<sxz<double> >(1, Tx[0]) );
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
                vTx.push_back( vector<sxz<double>>(1, Tx[ntx]) );
                t0.push_back( vector<double>(1, tTx[ntx]) );
                iTx.push_back( vector<size_t>(1, ntx) );
            }
        }

        /*
         Looping over all non redundant Tx
         */

        vector<sxz<double>> vRx;
        vector<vector<double>> tt( vTx.size() );

        if ( grid_instance->getNthreads() == 1 || vTx.size() <= grid_instance->getNthreads() ) {
            for ( size_t nv=0; nv<vTx.size(); ++nv ) {

                vRx.resize( 0 );
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ] );
                }

                try {
                    grid_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv]);
                } catch (std::exception& e) {
                    throw;
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
                threads[i]=thread( [&grid_ref,&vTx,&tt,&t0,&Rx,&iTx,
                                    blk_start,blk_end,i]{

                    for ( size_t nv=blk_start; nv<blk_end; ++nv ) {

                        vector<sxz<double>> vRx;
                        for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                            vRx.push_back( Rx[ iTx[nv][ni] ] );
                        }
                        try {
                            grid_ref->raytrace(vTx[nv], t0[nv], vRx, tt[nv], i+1);
                        } catch (std::exception& e) {
                            throw;
                        }
                    }
                });

                blk_start = blk_end;
            }

            for ( size_t nv=blk_start; nv<vTx.size(); ++nv ) {
                vector<sxz<double>> vRx;
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ] );
                }
                try {
                    grid_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], 0);
                } catch (std::exception& e) {
                    throw;
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
    }
    
    void Grid2Dttcr::raytrace(const std::vector<sxz<double>>& Tx,
                              const std::vector<double>& tTx,
                              const std::vector<sxz<double>>& Rx,
                              double* traveltimes,
                              PyObject* L) const {
        
        /*
         Looking for redundants Tx pts
         */

        size_t nTx = Tx.size();
        vector<vector<sxz<double>>> vTx;
        vector<vector<double>> t0;
        vector<vector<size_t>> iTx;
        vTx.push_back( vector<sxz<double> >(1, Tx[0]) );
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
                vTx.push_back( vector<sxz<double>>(1, Tx[ntx]) );
                t0.push_back( vector<double>(1, tTx[ntx]) );
                iTx.push_back( vector<size_t>(1, ntx) );
            }
        }

        /*
         Looping over all non redundant Tx
         */
        
        vector<sxz<double>> vRx;
        vector<vector<double>> tt( vTx.size() );
        vector<vector<siv2<double>>> L_data(nTx);
        vector<vector<vector<siv2<double>>>> l_data( vTx.size() );
        
        if ( grid_instance->getNthreads() == 1 || vTx.size() <= grid_instance->getNthreads() ) {
            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                
                vRx.resize( 0 );
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ] );
                }
                
                try {
                    grid_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], l_data[nv]);
                } catch (std::exception& e) {
                    throw;
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
                threads[i]=thread( [&grid_ref,&vTx,&tt,&t0,&Rx,&iTx,
                                    &l_data,blk_start,blk_end,i]{
                    
                    for ( size_t nv=blk_start; nv<blk_end; ++nv ) {
                        
                        vector<sxz<double>> vRx;
                        for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                            vRx.push_back( Rx[ iTx[nv][ni] ] );
                        }
                        try {
                            grid_ref->raytrace(vTx[nv], t0[nv], vRx, tt[nv], l_data[nv], i+1);
                        } catch (std::exception& e) {
                            throw;
                        }
                    }
                });
                
                blk_start = blk_end;
            }
            
            for ( size_t nv=blk_start; nv<vTx.size(); ++nv ) {
                vector<sxz<double>> vRx;
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ] );
                }
                try {
                    grid_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], l_data[nv], 0);
                } catch (std::exception& e) {
                    throw;
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
        import_array_throw();  // to use PyArray_SimpleNewFromData

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
        
        if ( type.compare("iso")==0 ) {
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
            
        } else {
            nnz *= 2;
            
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
                for ( size_t j=0; j<2*ncell; ++j ) {
                    for ( size_t n=0; n<L_data[i].size(); ++n ) {
                        if ( L_data[i][n].i == j ) {
                            indices_p[k] = j;
                            data_p[k] = L_data[i][n].v;
                            k++;
                        }    else if ( L_data[i][n].i+ncell == j ) {
                            indices_p[k] = j;
                            data_p[k] = L_data[i][n].v2;
                            k++;
                        }
                    }
                }
            }
            indptr_p[nTx] = k;
            
            PyTuple_SetItem(L, 0, data);
            PyTuple_SetItem(L, 1, indices);
            PyTuple_SetItem(L, 2, indptr);
        }
    }

    int Grid2Dttcr::Lsr2d(const double* Tx,
                          const double* Rx,
                          const size_t nTx,
                          const double* grx,
                          const size_t n_grx,
                          const double* grz,
                          const size_t n_grz,
                          PyObject* L) {

        const double  small=1.e-10;

        size_t nCells = (n_grx-1)*(n_grz-1);
        size_t nLmax = nTx * n_grx * n_grz/2;
        double percent_sp = (nLmax*1.0)/(nTx*nCells*1.0);

        double* data_p = (double*)malloc( nLmax*sizeof(double) );
        int64_t* indices_p = (int64_t*)malloc( nLmax*sizeof(int64_t) );
        int64_t* indptr_p = (int64_t*)malloc( (nTx+1)*sizeof(int64_t) );

        size_t k = 0;
        size_t ix, iz;
        for ( size_t n=0; n<nTx; ++n ) {
            indptr_p[n] = k;

            double xs = Tx[2*n];
            double zs = Tx[2*n+1];
            double xr = Rx[2*n];
            double zr = Rx[2*n+1];

            if ( xs>xr ) {  /* on va de s à r, on veut x croissant */
                double dtmp = xs;
                xs = xr;
                xr = dtmp;
                dtmp = zs;
                zs = zr;
                zr = dtmp;
            }

            /* points de depart */
            double x = xs;
            double z = zs;

            if ( fabs(zs-zr)<small ) {  /* rai horizontal */

                for ( ix=0; ix<n_grx-1; ++ix ) if ( x < grx[ix+1] ) break;
                for ( iz=0; iz<n_grz-1; ++iz ) if ( z < grz[iz+1] ) break;

                while ( x < xr ) {
                    int64_t iCell = ix*(n_grz-1) + iz;

                    double dlx = ( grx[ix+1]<xr ? grx[ix+1] : xr ) - x;

                    indices_p[k] = iCell;
                    data_p[k] = dlx;
                    k++;

                    if (k>=nLmax){
                        size_t oldnzmax = nLmax;
                        percent_sp += 0.1;
                        nLmax = (size_t)ceil((double)nTx*(double)nCells*percent_sp);

                        /* make sure nzmax increases at least by 1 */
                        if (oldnzmax == nLmax) nLmax++;

                        data_p = (double*)realloc( data_p, nLmax*sizeof(double) );
                        indices_p = (int64_t*)realloc( indices_p, nLmax*sizeof(int64_t) );
                    }

                    ix++;
                    x = grx[ix];
                }
            }
            else if ( fabs(xs-xr)<small ) { /* rai vertical */
                if ( zs > zr ) {  /* on va de s à r, on veut z croissant */
                    double dtmp = zs;
                    zs = zr;
                    zr = dtmp;
                }
                z = zs;

                for ( ix=0; ix<n_grx-1; ++ix ) if ( x < grx[ix+1] ) break;
                for ( iz=0; iz<n_grz-1; ++iz ) if ( z < grz[iz+1] ) break;

                while ( z < zr ) {
                    int64_t iCell = ix*(n_grz-1) + iz;

                    double dlz = ( grz[iz+1]<zr ? grz[iz+1] : zr ) - z;

                    indices_p[k] = iCell;
                    data_p[k] = dlz;
                    k++;

                    if (k>=nLmax){
                        size_t oldnzmax = nLmax;
                        percent_sp += 0.1;
                        nLmax = (size_t)ceil((double)nTx*(double)nCells*percent_sp);

                        /* make sure nzmax increases at least by 1 */
                        if (oldnzmax == nLmax) nLmax++;

                        data_p = (double*)realloc( data_p, nLmax*sizeof(double) );
                        indices_p = (int64_t*)realloc( indices_p, nLmax*sizeof(int64_t) );
                    }

                    iz++;
                    z = grz[iz];
                }
            }
            else { /* rai oblique */
                /* pente du rai */
                double m = (zr-zs)/(xr-xs);
                double b = zr - m*xr;
                bool up = m>0;

                for ( ix=0; ix<n_grx-1; ++ix ) if ( x < grx[ix+1] ) break;
                for ( iz=0; iz<n_grz-1; ++iz ) if ( z < grz[iz+1] ) break;

                while ( x < xr ) {

                    double zi = m*grx[ix+1] + b;

                    if ( up ) {
                        while ( z < zi && z < zr ) {
                            int64_t iCell = ix*(n_grz-1) + iz;

                            double ze = grz[iz+1]<zi ? grz[iz+1] : zi;
                            ze = ze<zr ? ze : zr;
                            double xe = (ze-b)/m;
                            double dlx = xe - x;
                            double dlz = ze - z;
                            double dl = sqrt( dlx*dlx + dlz*dlz );

                            indices_p[k] = iCell;
                            data_p[k] = dl;
                            k++;

                            if (k>=nLmax){
                                size_t oldnzmax = nLmax;
                                percent_sp += 0.1;
                                nLmax = (size_t)ceil((double)nTx*(double)nCells*percent_sp);

                                /* make sure nzmax increases at least by 1 */
                                if (oldnzmax == nLmax) nLmax++;

                                data_p = (double*)realloc( data_p, nLmax*sizeof(double) );
                                indices_p = (int64_t*)realloc( indices_p, nLmax*sizeof(int64_t) );
                            }

                            x = xe;
                            z = ze;
                            if ( fabs(z-grz[iz+1])<small ) iz++;
                        }
                    } else {
                        while ( z > zi && z > zr ) {
                            int64_t iCell = ix*(n_grz-1) + iz;

                            double ze = grz[iz]>zi ? grz[iz] : zi;
                            ze = ze>zr ? ze : zr;
                            double xe = (ze-b)/m;
                            double dlx = xe - x;
                            double dlz = ze - z;
                            double dl = sqrt( dlx*dlx + dlz*dlz );

                            indices_p[k] = iCell;
                            data_p[k] = dl;
                            k++;

                            if (k>=nLmax){
                                size_t oldnzmax = nLmax;
                                percent_sp += 0.1;
                                nLmax = (size_t)ceil((double)nTx*(double)nCells*percent_sp);

                                /* make sure nzmax increases at least by 1 */
                                if (oldnzmax == nLmax) nLmax++;

                                data_p = (double*)realloc( data_p, nLmax*sizeof(double) );
                                indices_p = (int64_t*)realloc( indices_p, nLmax*sizeof(int64_t) );
                            }

                            x = xe;
                            z = ze;
                            if ( fabs(z-grz[iz])<small ) iz--;
                        }
                    }

                    ix++;
                    x = grx[ix];
                }
            }
        }

        indptr_p[nTx] = k;
        size_t nnz = k;

        data_p = (double*)realloc( data_p, nnz*sizeof(double) );
        indices_p = (int64_t*)realloc( indices_p, nnz*sizeof(int64_t) );

        import_array_throw();  // to use PyArray_SimpleNewFromData

        npy_intp dims[] = {static_cast<npy_intp>(nnz)};
        PyObject* data = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, data_p);
        PyArray_ENABLEFLAGS((PyArrayObject*)data, NPY_ARRAY_OWNDATA);
        PyObject* indices = PyArray_SimpleNewFromData(1, dims, NPY_INT64, indices_p);
        PyArray_ENABLEFLAGS((PyArrayObject*)indices, NPY_ARRAY_OWNDATA);
        dims[0] = nTx+1;
        PyObject* indptr = PyArray_SimpleNewFromData(1, dims, NPY_INT64, indptr_p);
        PyArray_ENABLEFLAGS((PyArrayObject*)indptr, NPY_ARRAY_OWNDATA);
        PyTuple_SetItem(L, 0, data);
        PyTuple_SetItem(L, 1, indices);
        PyTuple_SetItem(L, 2, indptr);

        return 0;
    }


    int Grid2Dttcr::Lsr2da(const double* Tx,
                           const double* Rx,
                           const size_t nTx,
                           const double* grx,
                           const size_t n_grx,
                           const double* grz,
                           const size_t n_grz,
                           PyObject* L) {

        const double  small=1.e-10;

        size_t nCells = (n_grx-1)*(n_grz-1);
        size_t nLmax = nTx * n_grx * n_grz/2;
        double percent_sp = (nLmax*1.0)/(nTx*nCells*1.0);

        double* data_p = (double*)malloc( nLmax*sizeof(double) );
        int64_t* indices_p = (int64_t*)malloc( nLmax*sizeof(int64_t) );
        int64_t* indptr_p = (int64_t*)malloc( (nTx+1)*sizeof(int64_t) );

        size_t k = 0;
        size_t ix, iz;
        for ( size_t n=0; n<nTx; ++n ) {
            indptr_p[n] = k;

            double xs = Tx[2*n];
            double zs = Tx[2*n+1];
            double xr = Rx[2*n];
            double zr = Rx[2*n+1];

            if ( xs>xr ) {  /* on va de s à r, on veut x croissant */
                double dtmp = xs;
                xs = xr;
                xr = dtmp;
                dtmp = zs;
                zs = zr;
                zr = dtmp;
            }

            /* points de depart */
            double x = xs;
            double z = zs;

            if ( fabs(zs-zr)<small ) {  /* rai horizontal */

                for ( ix=0; ix<n_grx-1; ++ix ) if ( x < grx[ix+1] ) break;
                for ( iz=0; iz<n_grz-1; ++iz ) if ( z < grz[iz+1] ) break;

                while ( x < xr ) {
                    int64_t iCell = ix*(n_grz-1) + iz;

                    double dlx = ( grx[ix+1]<xr ? grx[ix+1] : xr ) - x;

                    indices_p[k] = iCell;
                    data_p[k] = dlx;
                    k++;

                    if (k>=nLmax){
                        size_t oldnzmax = nLmax;
                        percent_sp += 0.1;
                        nLmax = (size_t)ceil((double)nTx*(double)nCells*percent_sp);

                        /* make sure nzmax increases at least by 1 */
                        if (oldnzmax == nLmax) nLmax++;

                        data_p = (double*)realloc( data_p, nLmax*sizeof(double) );
                        indices_p = (int64_t*)realloc( indices_p, nLmax*sizeof(int64_t) );
                    }

                    ix++;
                    x = grx[ix];
                }
            }
            else if ( fabs(xs-xr)<small ) { /* rai vertical */
                if ( zs > zr ) {  /* on va de s à r, on veut z croissant */
                    double dtmp = zs;
                    zs = zr;
                    zr = dtmp;
                }
                z = zs;

                for ( ix=0; ix<n_grx-1; ++ix ) if ( x < grx[ix+1] ) break;
                for ( iz=0; iz<n_grz-1; ++iz ) if ( z < grz[iz+1] ) break;

                while ( z < zr ) {
                    int64_t iCell = ix*(n_grz-1) + iz;

                    double dlz = ( grz[iz+1]<zr ? grz[iz+1] : zr ) - z;

                    indices_p[k] = iCell+nCells;
                    data_p[k] = dlz;
                    k++;

                    if (k>=nLmax){
                        size_t oldnzmax = nLmax;
                        percent_sp += 0.1;
                        nLmax = (size_t)ceil((double)nTx*(double)nCells*percent_sp);

                        /* make sure nzmax increases at least by 1 */
                        if (oldnzmax == nLmax) nLmax++;

                        data_p = (double*)realloc( data_p, nLmax*sizeof(double) );
                        indices_p = (int64_t*)realloc( indices_p, nLmax*sizeof(int64_t) );
                    }

                    iz++;
                    z = grz[iz];
                }
            }
            else { /* rai oblique */
                /* pente du rai */
                double m = (zr-zs)/(xr-xs);
                double b = zr - m*xr;
                bool up = m>0;

                for ( ix=0; ix<n_grx-1; ++ix ) if ( x < grx[ix+1] ) break;
                for ( iz=0; iz<n_grz-1; ++iz ) if ( z < grz[iz+1] ) break;

                while ( x < xr ) {

                    double zi = m*grx[ix+1] + b;

                    if ( up ) {
                        while ( z < zi && z < zr ) {
                            int64_t iCell = ix*(n_grz-1) + iz;

                            double ze = grz[iz+1]<zi ? grz[iz+1] : zi;
                            ze = ze<zr ? ze : zr;
                            double xe = (ze-b)/m;
                            double dlx = xe - x;
                            double dlz = ze - z;

                            indices_p[k] = iCell;
                            data_p[k] = dlx;
                            k++;
                            indices_p[k] = iCell+nCells;
                            data_p[k] = dlz;
                            k++;

                            if (k>=nLmax){
                                size_t oldnzmax = nLmax;
                                percent_sp += 0.1;
                                nLmax = (size_t)ceil((double)nTx*(double)nCells*percent_sp);

                                /* make sure nzmax increases at least by 2 */
                                if (oldnzmax == nLmax) nLmax+=2;

                                data_p = (double*)realloc( data_p, nLmax*sizeof(double) );
                                indices_p = (int64_t*)realloc( indices_p, nLmax*sizeof(int64_t) );
                            }

                            x = xe;
                            z = ze;
                            if ( fabs(z-grz[iz+1])<small ) iz++;
                        }
                    } else {
                        while ( z > zi && z > zr ) {
                            int64_t iCell = ix*(n_grz-1) + iz;

                            double ze = grz[iz]>zi ? grz[iz] : zi;
                            ze = ze>zr ? ze : zr;
                            double xe = (ze-b)/m;
                            double dlx = xe - x;
                            double dlz = ze - z;

                            indices_p[k] = iCell;
                            data_p[k] = dlx;
                            k++;
                            indices_p[k] = iCell+nCells;
                            data_p[k] = dlz;
                            k++;

                            if (k>=nLmax){
                                size_t oldnzmax = nLmax;
                                percent_sp += 0.1;
                                nLmax = (size_t)ceil((double)nTx*(double)nCells*percent_sp);

                                /* make sure nzmax increases at least by 2 */
                                if (oldnzmax == nLmax) nLmax+=2;

                                data_p = (double*)realloc( data_p, nLmax*sizeof(double) );
                                indices_p = (int64_t*)realloc( indices_p, nLmax*sizeof(int64_t) );
                            }

                            x = xe;
                            z = ze;
                            if ( fabs(z-grz[iz])<small ) iz--;
                        }
                    }

                    ix++;
                    x = grx[ix];
                }
            }
        }
        indptr_p[nTx] = k;
        size_t nnz = k;

        data_p = (double*)realloc( data_p, nnz*sizeof(double) );
        indices_p = (int64_t*)realloc( indices_p, nnz*sizeof(int64_t) );

        import_array_throw();  // to use PyArray_SimpleNewFromData

        npy_intp dims[] = {static_cast<npy_intp>(nnz)};
        PyObject* data = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, data_p);
        PyArray_ENABLEFLAGS((PyArrayObject*)data, NPY_ARRAY_OWNDATA);
        PyObject* indices = PyArray_SimpleNewFromData(1, dims, NPY_INT64, indices_p);
        PyArray_ENABLEFLAGS((PyArrayObject*)indices, NPY_ARRAY_OWNDATA);
        dims[0] = nTx+1;
        PyObject* indptr = PyArray_SimpleNewFromData(1, dims, NPY_INT64, indptr_p);
        PyArray_ENABLEFLAGS((PyArrayObject*)indptr, NPY_ARRAY_OWNDATA);
        PyTuple_SetItem(L, 0, data);
        PyTuple_SetItem(L, 1, indices);
        PyTuple_SetItem(L, 2, indptr);

        return 0;
    }
}

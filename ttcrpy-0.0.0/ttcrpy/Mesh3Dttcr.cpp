//
//  Mesh3Dttcr.cpp
//  ttcr
//
//  Created by Bernard Giroux on 16-10-11.
//  Copyright © 2016 Bernard Giroux. All rights reserved.
//

#include <thread>

#include "Mesh3Dttcr.h"

using namespace std;

#define import_array_throw() {if (_import_array() < 0) {PyErr_Print(); PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import"); throw NUMPY_IMPORT_ARRAY_RETVAL; } }

namespace ttcr {
    
    Mesh3Dttcr::Mesh3Dttcr(const std::vector<sxyz<double>>& no,
                           const std::vector<tetrahedronElem<uint32_t>>& tet,
                           const double eps, const int maxit, const bool rp=false,
                           const size_t nt=1) {
        // find mesh "corners"
        double xmin = no[0].x;
        double xmax = no[0].x;
        double ymin = no[0].y;
        double ymax = no[0].y;
        double zmin = no[0].z;
        double zmax = no[0].z;
        for ( size_t n=1; n<no.size(); ++n ) {
            xmin = xmin < no[n].x ? xmin : no[n].x;
            xmax = xmax > no[n].x ? xmax : no[n].x;
            ymin = ymin < no[n].y ? ymin : no[n].y;
            ymax = ymax > no[n].y ? ymax : no[n].y;
            zmin = zmin < no[n].z ? zmin : no[n].z;
            zmax = zmax > no[n].z ? zmax : no[n].z;
        }
        
        // use corners as ref pts
        std::vector<sxyz<double>> refPts;
        refPts.push_back( {xmin, ymin, zmin} );
        refPts.push_back( {xmin, ymin, zmax} );
        refPts.push_back( {xmin, ymax, zmin} );
        refPts.push_back( {xmin, ymax, zmax} );
        refPts.push_back( {xmax, ymin, zmin} );
        refPts.push_back( {xmax, ymin, zmax} );
        refPts.push_back( {xmax, ymax, zmin} );
        refPts.push_back( {xmax, ymax, zmax} );
        mesh_instance = new mesh(no, tet, eps, maxit, refPts, 2, rp, false, true, 1.e-5, nt);
    }
    
    void Mesh3Dttcr::setSlowness(const std::vector<double>& slowness) {
        try {
            mesh_instance->setSlowness(slowness);
        } catch (length_error& e) {
            throw;
        }
    }

    void Mesh3Dttcr::raytrace(const std::vector<sxyz<double>>& Tx,
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
        
        if ( mesh_instance->getNthreads() == 1 || vTx.size() <= mesh_instance->getNthreads() ) {
            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                
                vRx.resize( 0 );
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ] );
                }
                try {
                    mesh_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv]);
                } catch (std::exception& e) {
                    throw;
                }
            }
        } else {
            size_t num_threads = mesh_instance->getNthreads();
            size_t blk_size = vTx.size()/num_threads;
            if ( blk_size == 0 ) blk_size++;
            
            vector<thread> threads(num_threads-1);
            size_t blk_start = 0;
            for ( size_t i=0; i<num_threads-1; ++i ) {
                
                size_t blk_end = blk_start + blk_size;
                mesh *mesh_ref = mesh_instance;
                threads[i]=thread( [&mesh_ref,&vTx,&tt,&t0,&Rx,&iTx,&nRx,
                                    blk_start,blk_end,i]{
                    
                    for ( size_t nv=blk_start; nv<blk_end; ++nv ) {
                        
                        vector<sxyz<double>> vRx;
                        for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                            vRx.push_back( Rx[ iTx[nv][ni] ] );
                        }
                        try {
                            mesh_ref->raytrace(vTx[nv], t0[nv], vRx, tt[nv], i+1);
                        } catch (std::exception& e) {
                            throw;
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
                try {
                    mesh_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], 0);
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
    
    void Mesh3Dttcr::raytrace(const std::vector<sxyz<double>>& Tx,
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
        
        if ( mesh_instance->getNthreads() == 1 || vTx.size() <= mesh_instance->getNthreads() ) {
            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                
                vRx.resize( 0 );
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ] );
                }
                
                try {
                    mesh_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv]);
                } catch (std::exception& e) {
                    throw;
                }
            }
        } else {
            size_t num_threads = mesh_instance->getNthreads();
            size_t blk_size = vTx.size()/num_threads;
            if ( blk_size == 0 ) blk_size++;
            
            vector<thread> threads(num_threads-1);
            size_t blk_start = 0;
            for ( size_t i=0; i<num_threads-1; ++i ) {
                
                size_t blk_end = blk_start + blk_size;
                mesh *mesh_ref = mesh_instance;
                threads[i]=thread( [&mesh_ref,&vTx,&tt,&t0,&Rx,&iTx,&nRx,
                                    &r_data,&v0,blk_start,blk_end,i]{
                    
                    for ( size_t nv=blk_start; nv<blk_end; ++nv ) {
                        
                        vector<sxyz<double>> vRx;
                        for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                            vRx.push_back( Rx[ iTx[nv][ni] ] );
                        }
                        try {
                            mesh_ref->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], i+1);
                        } catch (std::exception& e) {
                            throw;
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
                try {
                    mesh_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], 0);
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
                V0[ iTx[nv][ni] ] = v0[nv];
            }
        }

        // rays
        import_array_throw();  // to use PyArray_SimpleNewFromData
        
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
    }
    
//    void Mesh3Dttcr::raytrace(const std::vector<sxyz<double>>& Tx,
//                              const std::vector<double>& tTx,
//                              const std::vector<sxyz<double>>& Rx,
//                              double* traveltimes,
//                              PyObject* rays,
//                              double* V0,
//                              PyObject* M) const {
//        // rays must be a pointer to a tuple object of size nRx
//        
//        /*
//         Looking for redundants Tx pts
//         */
//        
//        size_t nTx = Tx.size();
//        size_t nRx = Rx.size();
//        vector<vector<sxyz<double>>> vTx;
//        vector<vector<double>> t0;
//        vector<vector<size_t>> iTx;
//        vTx.push_back( vector<sxyz<double> >(1, Tx[0]) );
//        t0.push_back( vector<double>(1, tTx[0]) );
//        iTx.push_back( vector<size_t>(1, 0) );  // indices of Rx corresponding to current Tx
//        for ( size_t ntx=1; ntx<nTx; ++ntx ) {
//            bool found = false;
//            
//            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
//                if ( vTx[nv][0]==Tx[ntx] ) {
//                    found = true;
//                    iTx[nv].push_back( ntx ) ;
//                    break;
//                }
//            }
//            if ( !found ) {
//                vTx.push_back( vector<sxyz<double>>(1, Tx[ntx]) );
//                t0.push_back( vector<double>(1, tTx[ntx]) );
//                iTx.push_back( vector<size_t>(1, ntx) );
//            }
//        }
//        
//        /*
//         Looping over all non redundant Tx
//         */
//        
//        vector<sxyz<double>> vRx;
//        vector<vector<double>> tt( vTx.size() );
//        vector<vector<vector<sxyz<double>>>> r_data( vTx.size() );
//        vector<double> v0( vTx.size() );
//        vector<vector<vector<sijv<double>>>> m_data( vTx.size() );
//        
//        if ( mesh_instance->getNthreads() == 1 || vTx.size() <= mesh_instance->getNthreads() ) {
//            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
//                
//                vRx.resize( 0 );
//                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
//                    vRx.push_back( Rx[ iTx[nv][ni] ] );
//                }
//                
//                try {
//                    mesh_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], m_data[nv]);
//                } catch (std::exception& e) {
//                    throw;
//                }
//            }
//        } else {
//            size_t num_threads = mesh_instance->getNthreads();
//            size_t blk_size = vTx.size()/num_threads;
//            if ( blk_size == 0 ) blk_size++;
//            
//            vector<thread> threads(num_threads-1);
//            size_t blk_start = 0;
//            for ( size_t i=0; i<num_threads-1; ++i ) {
//                
//                size_t blk_end = blk_start + blk_size;
//                mesh *mesh_ref = mesh_instance;
//                threads[i]=thread( [&mesh_ref,&vTx,&tt,&t0,&Rx,&iTx,&nRx,
//                                    &r_data,&v0,&m_data,blk_start,blk_end,i]{
//                    
//                    for ( size_t nv=blk_start; nv<blk_end; ++nv ) {
//                        
//                        vector<sxyz<double>> vRx;
//                        for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
//                            vRx.push_back( Rx[ iTx[nv][ni] ] );
//                        }
//                        try {
//                            mesh_ref->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], m_data[nv], i+1);
//                        } catch (std::exception& e) {
//                            throw;
//                        }
//                    }
//                });
//                
//                blk_start = blk_end;
//            }
//            
//            for ( size_t nv=blk_start; nv<vTx.size(); ++nv ) {
//                vector<sxyz<double>> vRx;
//                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
//                    vRx.push_back( Rx[ iTx[nv][ni] ] );
//                }
//                try {
//                    mesh_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], m_data[nv], 0);
//                } catch (std::exception& e) {
//                    throw;
//                }
//            }
//            
//            std::for_each(threads.begin(),threads.end(),
//                          std::mem_fn(&std::thread::join));
//        }
//
//        for ( size_t nv=0; nv<vTx.size(); ++nv ) {
//            for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
//                traveltimes[ iTx[nv][ni] ] = tt[nv][ni];
//                V0[ iTx[nv][ni] ] = v0[nv];
//            }
//        }
//        
//        // rays
//        import_array_throw();  // to use PyArray_SimpleNewFromData
//        
//        for ( size_t nv=0; nv<vTx.size(); ++nv ) {
//            for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
//                size_t npts = r_data[nv][ni].size();
//                npy_intp dims[] = {static_cast<npy_intp>(npts), 3};
//                double* ray_p = new double[3*npts];
//                PyObject* ray = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, ray_p);
//                
//                for ( size_t np=0; np<npts; ++np ) {
//                    ray_p[3*np] = r_data[nv][ni][np].x;
//                    ray_p[3*np+1] = r_data[nv][ni][np].y;
//                    ray_p[3*np+2] = r_data[nv][ni][np].z;
//                }
//                
//                PyTuple_SetItem(rays, iTx[nv][ni], ray);
//            }
//        }
//
//        // data for matrix M
//        // M is a tuple of tuples
//        // first element of tuple contains data, size is nnz
//        // second element contains column indices for rows of the matrix, size is nnz
//        // third element contains pointers for indices, size is nrow+1 (nTx+1)
//        
//        size_t nnodes = mesh_instance->getNumberOfNodes();
//        
//        for ( size_t nv=0; nv<vTx.size(); ++nv ) {
//            
//            PyObject* tuple = PyTuple_New(3);
//            
//            size_t nRcv = m_data[nv].size();
//            size_t nnz = 0;
//            for ( size_t ni=0; ni<m_data[nv].size(); ++ni ) {
//                nnz += m_data[nv][ni].size();
//            }
//            
//            
//            npy_intp dims[] = {static_cast<npy_intp>(nnz)};
//            double* data_p = new double[nnz];
//            PyObject* data = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, data_p);
//            
//            int64_t* indices_p = new int64_t[nnz];
//            PyObject* indices = PyArray_SimpleNewFromData(1, dims, NPY_INT64, indices_p);
//            
//            dims[0] = nRcv+1;
//            int64_t* indptr_p = new int64_t[nRcv+1];
//            PyObject* indptr = PyArray_SimpleNewFromData(1, dims, NPY_INT64, indptr_p);
//            
//            size_t k = 0;
//            for ( size_t ni=0; ni<nRcv; ++ni ) {
//                indptr_p[ni] = k;
//                for ( size_t j=0; j<nnodes; ++j ) {
//                    for ( size_t n=0; n<m_data[nv][ni].size(); ++n) {
//                        if ( m_data[nv][ni][n].j == j && m_data[nv][ni][n].i == ni ) {
//                            indices_p[k] = j;
//                            data_p[k] = m_data[nv][ni][n].v;
//                            k++;
//                        }
//                    }
//                }
//            }
//            indptr_p[nRcv] = k;
//            
//            PyTuple_SetItem(tuple, 0, data);
//            PyTuple_SetItem(tuple, 1, indices);
//            PyTuple_SetItem(tuple, 2, indptr);
//            
//            PyTuple_SetItem(M, nv, tuple);
//        }
//    }

}

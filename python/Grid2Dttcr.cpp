//
//  Grid2Dttcr.cpp
//  ttcr
//
//  Created by Bernard Giroux on 16-07-28.
//  Copyright © 2016 Bernard Giroux. All rights reserved.
//

#include <thread>

#include "Grid2Dttcr.h"

using namespace std;


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
        if ( grid_instance->setSlowness(slowness) == 1 ) {
            throw out_of_range("Slowness values must be defined for each grid cell.");
        }
    }
    
    void Grid2Dttcr::setXi(const std::vector<double>& xi) {
        if ( grid_instance->setXi(xi) == 1 ) {
            throw out_of_range("Xi values must be defined for each grid cell.");
        }
    }
    
    void Grid2Dttcr::setTheta(const std::vector<double>& theta) {
        if ( grid_instance->setTiltAngle(theta) == 1 ) {
            throw out_of_range("Theta values must be defined for each grid cell.");
        }
    }
    
    int Grid2Dttcr::raytrace(const std::vector<sxz<double>>& Tx,
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
        size_t nRx = Rx.size();
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
        
        if ( grid_instance->getNthreads() == 1 ) {
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
            size_t num_threads = grid_instance->getNthreads() < vTx.size() ? grid_instance->getNthreads() : vTx.size();
            size_t const blk_size = (vTx.size()%num_threads ? 1 : 0) + vTx.size()/num_threads;
            
            vector<thread> threads(num_threads-1);
            size_t blk_start = 0;
            for ( size_t i=0; i<num_threads-1; ++i ) {
                
                size_t blk_end = blk_start + blk_size;
                grid *grid_ref = grid_instance;
                threads[i]=thread( [&grid_ref,&vTx,&tt,&t0,&Rx,&iTx,&nRx,
                                    &r_data,&l_data,blk_start,blk_end,i]{
                    
                    for ( size_t nv=blk_start; nv<blk_end; ++nv ) {
                        
                        vector<sxz<double>> vRx;
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
                vector<sxz<double>> vRx;
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
                npy_intp dims[] = {static_cast<npy_intp>(npts), 2};
                double* ray_p = new double[2*npts];
                PyObject* ray = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, ray_p);
                
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
            
            int64_t* indices_p = new int64_t[nnz];
            PyObject* indices = PyArray_SimpleNewFromData(1, dims, NPY_INT64, indices_p);
            
            dims[0] = nTx+1;
            int64_t* indptr_p = new int64_t[nTx+1];
            PyObject* indptr = PyArray_SimpleNewFromData(1, dims, NPY_INT64, indptr_p);
            
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
            
            int64_t* indices_p = new int64_t[nnz];
            PyObject* indices = PyArray_SimpleNewFromData(1, dims, NPY_INT64, indices_p);
            
            dims[0] = nTx+1;
            int64_t* indptr_p = new int64_t[nTx+1];
            PyObject* indptr = PyArray_SimpleNewFromData(1, dims, NPY_INT64, indptr_p);
            
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
        
        return 0;
    }
}
//
//  Mesh3Dttcr.cpp
//  ttcr
//
//  Created by Bernard Giroux on 16-10-11.
//  Copyright © 2016 Bernard Giroux. All rights reserved.
//

#include <thread>

#include "Mesh3Dttcr.h"
#include "fstream"
#include  <iostream>
using namespace std;

namespace ttcr {
    
    Mesh3Dttcr::Mesh3Dttcr(const std::vector<sxyz<double>>& no,
                           const std::vector<tetrahedronElem<uint32_t>>& tet,
                           const int ns, const size_t nt=1, const int verbose=0,
                           const size_t ns2=0, const double r=0.0  ) {
        std::cout<<"Hello"<<std::endl;
        MinCorner=no[0];
        for (size_t nn=1;nn<no.size();++nn){
            if(no[nn].x<=MinCorner.x && no[nn].y <= MinCorner.y && no[nn].z <=MinCorner.z)
                MinCorner=no[nn];
        }
        std::vector<sxyz<double>> nodes (no.size());
        for (size_t nn=0;nn<no.size();++nn){
            nodes[nn]=(no[nn]-MinCorner);
        }
        mesh_instance = new mesh(nodes, tet, ns,nt,verbose);
        nst=ns2;
        Radius=r;

    }
    bool Mesh3Dttcr::CheckPoint(const std::vector <sxyz<double>>& Points) const{
        std::vector<sxyz<double>>  rltve_Pnt(Points.size());
        for (size_t np=0;np<Points.size();++np)
            rltve_Pnt[np]=Points[np]-MinCorner;
            //sxyz<double>{Points[np].x-MinCorner.x,Points[np].y-MinCorner.y,Points[np].z-MinCorner.z};
        return (!mesh_instance->checkPts(rltve_Pnt));
    }
    void Mesh3Dttcr::setSlowness(const std::vector<double>& slowness) {
        if ( mesh_instance->setSlowness(slowness) == 1 ) {
            throw out_of_range("Slowness values must be defined for each mesh node.");
        }
    }

    int Mesh3Dttcr::ComputeD(const std::vector<sxyz<double>>& Pnts,
                 PyObject* D)const{
        
        std::vector<std::vector<sijv<double>>> d_data(Pnts.size());
        mesh_instance->computeD(Pnts,d_data);
        D = PyTuple_New(3);
        size_t nnz = 0;
        for ( size_t np=0; np<d_data.size(); ++np ) {
            nnz += d_data[np].size();
        }
        npy_intp dims[] = {static_cast<npy_intp>(nnz)};
        double* data_p = new double[nnz];
        PyObject* data = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, data_p);
        
        int64_t* indices_p = new int64_t[nnz];
        PyObject* indices = PyArray_SimpleNewFromData(1, dims, NPY_INT64, indices_p);
        
        dims[0] = d_data.size()+1;
        int64_t* indptr_p = new int64_t[d_data.size()+1];
        PyObject* indptr = PyArray_SimpleNewFromData(1, dims, NPY_INT64, indptr_p);
        
        //size_t nnodes(mesh_instance->getNumberOfNodes());
        size_t k = 0;
        for ( size_t np=0; np<d_data.size(); ++np ) {
            indptr_p[np] = k;
                for ( size_t n=0; n<d_data[np].size(); ++n) {
                    if(d_data[np][n].i==np){
                    indices_p[k] =  d_data[np][n].j;
                        data_p[k] = d_data[np][n].v;
                        k++;
                    }
            }
        }
        indptr_p[d_data.size()+1] = k;
        
        PyTuple_SetItem(D, 0, data);
        PyTuple_SetItem(D, 1, indices);
        PyTuple_SetItem(D, 2, indptr);
        return 0;
    }


    int Mesh3Dttcr::raytrace(const std::vector<sxyz<double>>& Tx,
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
        vTx.push_back( vector<sxyz<double> >(1, Tx[0]-MinCorner) );
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
                vTx.push_back( vector<sxyz<double>>(1, Tx[ntx]-MinCorner) );
                t0.push_back( vector<double>(1, tTx[ntx]) );
                iTx.push_back( vector<size_t>(1, ntx) );
            }
        }
        
        /*
         Looping over all non redundant Tx
         */
        
        vector<sxyz<double>> vRx;
        vector<vector<double>> tt( vTx.size() );
        
        if ( mesh_instance->getNthreads() == 1 || mesh_instance->getNthreads()>= vTx.size() ) {
            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                
                vRx.resize( 0 );
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ]-MinCorner);
                }
                mesh_instance->AddNewSecondaryNodes(vTx[nv][0],nst, Radius);
                if (mesh_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv]) == 1 )
                        throw runtime_error("Problem while raytracing.");
                mesh_instance->DelateNodes();
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
                double Radius_ref=Radius;
                size_t nst_ref=nst;
                sxyz<double> MinCroner_ref=MinCorner;
        
                threads[i]=thread( [&mesh_ref,&vTx,&tt,&t0,&Rx,&iTx,&nRx,& nst_ref,&Radius_ref,
                                    & MinCroner_ref,blk_start,blk_end,i]{
                    
                    for ( size_t nv=blk_start; nv<blk_end; ++nv ) {
                        
                        vector<sxyz<double>> vRx;
                        for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                            vRx.push_back( Rx[ iTx[nv][ni] ]-MinCroner_ref);
                        }
                        mesh_ref->AddNewSecondaryNodes(vTx[nv][0],nst_ref, Radius_ref);
                        if ( mesh_ref->raytrace(vTx[nv], t0[nv], vRx, tt[nv], i+1) == 1 )
                                throw runtime_error("Problem while raytracing.");
                         mesh_ref->DelateNodes();
                    }
                    
                });
                
                blk_start = blk_end;
            }
            
            for ( size_t nv=blk_start; nv<vTx.size(); ++nv ) {
                vector<sxyz<double>> vRx;
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ] -MinCorner);
                }
                mesh_instance->AddNewSecondaryNodes(vTx[nv][0],nst, Radius);
                if ( mesh_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], 0) == 1 )
                        throw runtime_error("Problem while raytracing.");
                mesh_instance->DelateNodes();
                
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
    
    int Mesh3Dttcr::raytrace(const std::vector<sxyz<double>>& Tx,
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
        vTx.push_back( vector<sxyz<double> >(1, Tx[0]-MinCorner) );
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
                vTx.push_back( vector<sxyz<double>>(1, Tx[ntx]-MinCorner) );
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
        
        if ( mesh_instance->getNthreads() == 1 || mesh_instance->getNthreads()>= vTx.size() ) {
            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                
                vRx.resize( 0 );
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ]-MinCorner);
                }
                mesh_instance->AddNewSecondaryNodes(vTx[nv][0],nst, Radius);
                if ( mesh_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv]) == 1 )
                        throw runtime_error("Problem while raytracing.");
                mesh_instance->DelateNodes();
    
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
                size_t nst_ref=nst;
                double Radius_ref=Radius;
                sxyz<double> MinCorner_ref=MinCorner;
                threads[i]=thread( [&mesh_ref,&vTx,&tt,&t0,&Rx,&iTx,&nRx,& nst_ref,&Radius_ref,& MinCorner_ref,
                                    &r_data,&v0,blk_start,blk_end,i]{
                    
                    for ( size_t nv=blk_start; nv<blk_end; ++nv ) {
                        
                        vector<sxyz<double>> vRx;
                        for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                            vRx.push_back( Rx[ iTx[nv][ni]] - MinCorner_ref);
                        }
                        mesh_ref->AddNewSecondaryNodes(vTx[nv][0],nst_ref, Radius_ref);
                        if ( mesh_ref->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], i+1) == 1 )
                                throw runtime_error("Problem while raytracing.");
                        mesh_ref->DelateNodes();
                        
                    }
                });
                
                blk_start = blk_end;
            }
            
            for ( size_t nv=blk_start; nv<vTx.size(); ++nv ) {
                vector<sxyz<double>> vRx;
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ]-MinCorner);
                }
                mesh_instance->AddNewSecondaryNodes(vTx[nv][0],nst, Radius);
                if ( mesh_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], 0) == 1 )
                        throw runtime_error("Problem while raytracing.");
                mesh_instance->DelateNodes();
                
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
                    ray_p[3*np] = r_data[nv][ni][np].x+MinCorner.x;
                    ray_p[3*np+1] = r_data[nv][ni][np].y+MinCorner.y;
                    ray_p[3*np+2] =r_data[nv][ni][np].z+MinCorner.z;
                }
                
                PyTuple_SetItem(rays, iTx[nv][ni], ray);
            }
        }
        
        return 0;
    }
    
    int Mesh3Dttcr::raytrace(const std::vector<sxyz<double>>& Tx,
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
        vTx.push_back( vector<sxyz<double> >(1, Tx[0]-MinCorner) );
        t0.push_back( vector<double>(1, tTx[0]) );
        iTx.push_back( vector<size_t>(1, 0) );  // indices of Rx corresponding to current Tx
        for ( size_t ntx=1; ntx<nTx; ++ntx ) {
            bool found = false;
            
            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                if ( vTx[nv][0].getDistance({Tx[ntx]-MinCorner})<1.e-10 ) {
                    found = true;
                    iTx[nv].push_back( ntx ) ;
                    break;
                }
            }
            if ( !found ) {
                vTx.push_back( vector<sxyz<double>>(1, Tx[ntx]-MinCorner) );
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
        
        if ( mesh_instance->getNthreads() == 1 || mesh_instance->getNthreads()>= vTx.size() ) {
            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                
                vRx.resize( 0 );
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ]-MinCorner);
                }
                mesh_instance->AddNewSecondaryNodes(vTx[nv][0],nst, Radius);
                if ( mesh_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], m_data[nv],0) == 1 )
                        throw runtime_error("Problem while raytracing.");
                mesh_instance->DelateNodes();
                
//                std::ofstream outfile;
                
//                outfile.open("Dubeg.dat", std::ios_base::app);
//                for (size_t i=0; i<m_data[nv].size();++i){
//                    outfile <<"-----------"<<"\n";
//                    for (size_t j=0;j<m_data[nv][i].size();++j)
//                        outfile << m_data[nv][i][j].v<<"\n";
//                        }
//                outfile.close();
                
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
                double Radius_ref=Radius;
                size_t nst_ref=nst;
                sxyz<double> MinCorner_ref=MinCorner;
                threads[i]=thread( [&mesh_ref,&vTx,&tt,&t0,&Rx,&iTx,&nRx,&Radius_ref,& nst_ref,& MinCorner_ref,
                                    &r_data,&v0,&m_data,blk_start,blk_end,i]{
                    
                    for ( size_t nv=blk_start; nv<blk_end; ++nv ) {
                        
                        vector<sxyz<double>> vRx;
                        for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                            vRx.push_back( Rx[ iTx[nv][ni] ]-MinCorner_ref );
                            vTx[nv][ni]=vTx[nv][ni];
                        }
                        mesh_ref->AddNewSecondaryNodes(vTx[nv][0],nst_ref, Radius_ref);
                            if ( mesh_ref->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], m_data[nv], i+1) == 1 )
                                throw runtime_error("Problem while raytracing.");
                        mesh_ref->DelateNodes();
                        
                    }
                });
                
                blk_start = blk_end;
            }
            
            for ( size_t nv=blk_start; nv<vTx.size(); ++nv ) {
                vector<sxyz<double>> vRx;
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ]-MinCorner);
                }
                mesh_instance->AddNewSecondaryNodes(vTx[nv][0],nst, Radius);
                if ( mesh_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], m_data[nv], 0) == 1 )
                        throw runtime_error("Problem while raytracing.");
                mesh_instance->DelateNodes();
                
            }
            
            std::for_each(threads.begin(),threads.end(),
                          std::mem_fn(&std::thread::join));
        }
//        std::ofstream outfile;
//        outfile.open("Debug.dat", std::ios_base::app);
//        outfile << m_data<<"\n";
//        outfile.close();
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
                    ray_p[3*np] = r_data[nv][ni][np].x+MinCorner.x;
                    ray_p[3*np+1] = r_data[nv][ni][np].y+MinCorner.y;
                    ray_p[3*np+2] = r_data[nv][ni][np].z+MinCorner.z;
                }
                
                PyTuple_SetItem(rays, iTx[nv][ni], ray);
            }
        }
        
        // data for matrix M
        // M is a tuple of tuples
        // first element of tuple contains data, size is nnz
        // second element contains column indices for rows of the matrix, size is nnz
        // third element contains pointers for indices, size is nrow+1 (nTx+1)
        
        size_t nnodes = mesh_instance->getNumberOfNodes();
        for ( size_t nv=0; nv<vTx.size(); ++nv ) {
            
            PyObject* tuple = PyTuple_New(3);
            
            size_t nRcv = m_data[nv].size();
            size_t nnz = nRcv;
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
                            data_p[k] =m_data[nv][ni][n].v;
                            k++;
                        }
                    }
                }
                indices_p[k] = nnodes+ni;
                data_p[k] = 1.0;
                k++;
            }
            indptr_p[nRcv] = k;
            
            PyTuple_SetItem(tuple, 0, data);
            PyTuple_SetItem(tuple, 1, indices);
            PyTuple_SetItem(tuple, 2, indptr);
            
            PyTuple_SetItem(M, nv, tuple);
        }
        
        return 0;
    }

}

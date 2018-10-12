//
//  main.cpp
//  testM
//
//  Created by Bernard Giroux on 17-03-17.
//  Copyright © 2017 Bernard Giroux. All rights reserved.
//

#include <iostream>
#include <vector>
#include <thread>
#include <string>
#include <Eigen/Dense>

#include "MSHReader.h"
#include "Grid3Duifs.h"
#include "Grid3Duifm.h"
#include "Grid3Duisp.h"
#include "Rcv.h"
#include "Src.h"
#include "utils.h"
#include "Node.h"
using namespace ttcr;



int main(int argc, const char * argv[]) {

    std::string method="SPM";
    std::vector<Src<double >> src;
    size_t N=0;
//    for(size_t i=0;i<15;i++){
//        src.push_back( Src<double >("Src"+to_string(i+1+N)+".dat") );
//        src[i].init();
//    }
    src.push_back( Src<double >("Src1.dat") );
    src[0].init();
    
    //bool SecondNodes=false;
    std::vector<sxyz<double >> Ttx;
    std::vector<sxyz<double >> s;
    for (size_t i=0;i<src.size();i++){
        s=src[i].get_coord();
        Ttx.push_back(s.back());
    }
    std::vector<Rcv<double >> R;
    R.push_back( Rcv<double >("rcv.dat") );
    R[0].init(src.size());
    std::vector<sxyz<double >> Rrx;
    for(size_t i=0;i<R.size();i++){
        for(size_t j=0;j<R[i].get_coord().size();j++){
            Rrx.push_back(R[i].get_coord()[j]);
        }
    }
    size_t nTtx=Ttx.size();
    size_t nRrx=Rrx.size();
    std::vector<sxyz<double >> Tx(nTtx*nRrx);
    std::vector<sxyz<double >> Rx(nTtx*nRrx);
    for (unsigned i=0; i<nTtx;i++){
        for(size_t j=0; j<nRrx;j++){
            Tx[i*nRrx+j]=Ttx[i];
            Rx[i*nRrx+j]=Rrx[j];
        }
    }
    MSHReader reader( "Model1.msh" );
    std::vector<sxyz<double >> node(reader.getNumberOfNodes());
    std::vector<tetrahedronElem<uint32_t>> tetrahedra(reader.getNumberOfTetra());
    std::vector<double > slowness(reader.getNumberOfNodes());

    reader.readNodes3D(node);
    std::vector <sxyz<double>> nodes;
    for ( size_t n=0; n<node.size(); ++n ) {
        nodes.push_back(sxyz<double>(node[n].x,node[n].y,node[n].z));
    }
    double ic[] = {1.e20, 1.e20,1.e20};
    for ( size_t n=1; n<nodes.size(); ++n ) {
        if ( nodes[n].x<ic[0]) ic[0] = nodes[n].x;
        if ( nodes[n].y<ic[1]) ic[1] = nodes[n].y;
        if ( nodes[n].z<ic[2]) ic[2] = nodes[n].z;
    }
    for ( size_t n=1; n<nodes.size(); ++n ) {
        nodes[n].x-=ic[0];
        nodes[n].y-=ic[1];
        nodes[n].z-=ic[2];
    }
    sxyz<double> PP(nodes[0]);
    reader.readTetrahedronElements(tetrahedra);

    std::ifstream fin( "Model1.slo" );
    if ( !fin ) {
        std::cerr << "Error: cannot open file model1 " << std::endl;
        exit ( -1);
    }
    std::vector<double > tmp1;
    double dtmp;
    fin >> dtmp;
    while ( fin ) {
        tmp1.push_back( dtmp );
        fin >> dtmp;
    }
    std::vector<double> tmp(tmp1.begin(),tmp1.end());
    fin.close();

    ///////////////////////////////////////////////////////////////////
    size_t start_s=clock();
    if (method=="SPM"){
        Grid3Duisp<double, uint32_t> * mesh_instancesp;
        unsigned ns=7;
        mesh_instancesp=new Grid3Duisp<double, uint32_t>(nodes,tetrahedra,ns,1);
         mesh_instancesp->setSlowness(tmp);
        //    std::vector<double> Slow(nodes.size(),1.0/3.0);
        //    mesh_instancesp->setSlowness(Slow);
        //mesh_instancesp->setSourceRadius(10.0*1.e-3);
        /////////

        size_t nTx = Tx.size();      
        size_t nRx=Rx.size();
        std::vector<std::vector<sxyz<double >>> vTx;
        std::vector<std::vector<double >> t0;
        std::vector<std::vector<size_t>> iTx;
        vTx.push_back( std::vector<sxyz<double >>(1,{Tx[0].x-ic[0],Tx[0].y-ic[1],Tx[0].z-ic[2]}) );
        t0.push_back(std::vector<double >(1,0) );
        iTx.push_back(std::vector<size_t>(1, 0) );  // indices of Rx corresponding to current Tx
        for ( size_t ntx=1; ntx<nTx; ++ntx ) {
            bool found = false;
            
            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                if ( vTx[nv][0].getDistance( {Tx[ntx].x-ic[0],Tx[ntx].y-ic[1],Tx[ntx].z-ic[2]})<0.0000000001) {
                    found = true;
                    iTx[nv].push_back( ntx ) ;
                    break;
                }
            }
            if ( !found ) {
                vTx.push_back( std::vector<sxyz<double >>(1, {Tx[ntx].x-ic[0],Tx[ntx].y-ic[1],Tx[ntx].z-ic[2]}) );
                t0.push_back( std::vector<double >(1,0) );
                iTx.push_back(std::vector<size_t>(1, ntx) );
            }
        }
        /*
         Looping over all non redundant Tx
         */
        //    std::vector<std::vector<sijv<double>>> d_data;
        //    mesh_instance->computeD(Rx, d_data);
        
        std::vector<sxyz<double >> vRx;
        std::vector<std::vector<double >> tt( vTx.size() );
        std::vector<std::vector<std::vector<sxyz<double >>>> r_data( vTx.size() );
        std::vector<double > v0( vTx.size() );
        std::vector<std::vector<std::vector<sijv<double >>>> m_data( vTx.size() );
        
        if ( mesh_instancesp->getNthreads() == 1 ) {
            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                vRx.resize( 0 );
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( {Rx[ iTx[nv][ni] ].x-ic[0],Rx[ iTx[nv][ni] ].y-ic[1],
                        Rx[ iTx[nv][ni] ].z-ic[2]});
                }
                
              mesh_instancesp->AddNewSecondaryNodes(vTx[nv][0],1,0.150,0);
                //mesh_instancesp->AddNewSecondaryNodes(vTx[nv][0], 2, 0.035,0);
                
                if ( mesh_instancesp->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv],v0[nv],m_data[nv],0) == 1 )
                    throw std::runtime_error("Problem while raytracing.");
                size_t stop_s=clock();
                std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << endl;
                mesh_instancesp->DelateNodes();
                for (size_t ii=0;ii<r_data[nv].size();ii++){
                    std::cout<<" size r_data["<<ii<<"] :"<< r_data[nv][ii].size()<<std::endl;
                }
                std::cout<<"\n"<<std::endl;
                for (size_t ii=0;ii<m_data[nv].size();ii++){
                    //for (size_t jj=0;jj<m_data[nv][ii].size();++jj)
                    std::cout<<" size m_data["<<ii<<"] :"<< m_data[nv][ii].size()<<std::endl;
                }
                std::cout<<"\n"<<std::endl;
                
            }
        } else {
            size_t num_threads = mesh_instancesp->getNthreads();
            size_t blk_size = vTx.size()/num_threads;
            if ( blk_size == 0 ) blk_size++;
            
            std::vector<std::thread> threads(num_threads-1);
            size_t blk_start = 0;
            for ( size_t i=0; i<num_threads-1; ++i ) {
                
                size_t blk_end = blk_start + blk_size;
                Grid3Duisp<double , uint32_t> *mesh_ref = mesh_instancesp;
                threads[i]=std::thread( [&mesh_ref,&vTx,&tt,&t0,&Rx,&iTx,&nRx,
                                         &r_data,&v0,&m_data,blk_start,blk_end,i]{
                    
                    for ( size_t nv=blk_start; nv<blk_end; ++nv ) {
                        
                        std::vector<sxyz<double >> vRx;
                        for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                            vRx.push_back( Rx[ iTx[nv][ni] ] );
                        }
                        if ( mesh_ref->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv],m_data[nv],i+1) == 1 ) {
                            throw std::runtime_error("Problem while raytracing.");
                        }
                    }
                });
                
                blk_start = blk_end;
            }
            
            for ( size_t nv=blk_start; nv<vTx.size(); ++nv ) {
                std::vector<sxyz<double >> vRx;
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    vRx.push_back( Rx[ iTx[nv][ni] ] );
                }
                if ( mesh_instancesp->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv],v0[nv],m_data[nv], 0) == 1 ) {
                    throw std::runtime_error("Problem while raytracing.");
                }
            }
            
            std::for_each(threads.begin(),threads.end(),
                          std::mem_fn(&std::thread::join));
        }
        
        for (size_t i=0;i<r_data.size();++i){
            std::string filename="Src"+to_string(N+1+i)+".vtp";
            saveRayPaths(filename, r_data[i]);
        }
        size_t stop_s=clock();
        std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << endl;
        mesh_instancesp->DelateNodes();
        for ( size_t nv=0; nv<vTx.size(); ++nv ) {
            for (size_t ii=0;ii<r_data[nv].size();ii++){
                std::cout<<" size r_data["<<ii<<"] :"<< r_data[nv][ii].size()<<std::endl;
            }
        }
        std::cout<<"\n"<<std::endl;
        cout<<"\n"<<endl;
        mesh_instancesp->saveTT("time",0,0,false);
    }
    if (method=="FMM"){///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
            Grid3Duifm<double, uint32_t> * mesh_instancefm;
            mesh_instancefm=new Grid3Duifm<double, uint32_t>(nodes,tetrahedra,true,1);
            mesh_instancefm->setSlowness(tmp);
            //    std::vector<double> Slow(nodes.size(),1.0/3.0);
            //    mesh_instancesp->setSlowness(Slow);
            //mesh_instancefm->setSourceRadius(30.0*1.e-3);
            /////////
            
            size_t nTx = Tx.size();
            size_t nRx=Rx.size();
            std::vector<std::vector<sxyz<double >>> vTx;
            std::vector<std::vector<double >> t0;
            std::vector<std::vector<size_t>> iTx;
            vTx.push_back( std::vector<sxyz<double >>(1,Tx[0]) );
            t0.push_back(std::vector<double >(1,0) );
            iTx.push_back(std::vector<size_t>(1, 0) );  // indices of Rx corresponding to current Tx
            for ( size_t ntx=1; ntx<nTx; ++ntx ) {
                bool found = false;
                
                for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                    if ( vTx[nv][0].getDistance(Tx[ntx])<0.00000001) {
                        found = true;
                        iTx[nv].push_back( ntx ) ;
                        break;
                    }
                }
                if ( !found ) {
                    vTx.push_back( std::vector<sxyz<double >>(1, Tx[ntx]) );
                    t0.push_back( std::vector<double >(1,0) );
                    iTx.push_back(std::vector<size_t>(1, ntx) );
                }
            }
            /*
             Looping over all non redundant Tx
             */
            //    std::vector<std::vector<sijv<double>>> d_data;
            //    mesh_instance->computeD(Rx, d_data);
            
            std::vector<sxyz<double >> vRx;
            std::vector<std::vector<double >> tt( vTx.size() );
            std::vector<std::vector<std::vector<sxyz<double >>>> r_data( vTx.size() );
            std::vector<double > v0( vTx.size() );
            std::vector<std::vector<std::vector<sijv<double >>>> m_data( vTx.size() );
            
            if ( mesh_instancefm->getNthreads() == 1 || mesh_instancefm->getNthreads()<= vTx.size() ) {
                for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                    vRx.resize( 0 );
                    for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                        vRx.push_back( Rx[ iTx[nv][ni] ]);
                    }
                    
                    
                    if ( mesh_instancefm->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv],0) == 1 )
                        throw std::runtime_error("Problem while raytracing.");
                    size_t stop_s=clock();
                    std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << endl;
                    for (size_t ii=0;ii<r_data[nv].size();ii++){
                        std::cout<<" size r_data["<<ii<<"] :"<< r_data[nv][ii].size()<<std::endl;
                    }
                    std::cout<<"\n"<<std::endl;
                    for (size_t ii=0;ii<m_data[nv].size();ii++){
                        //for (size_t jj=0;jj<m_data[nv][ii].size();++jj)
                        std::cout<<" size m_data["<<ii<<"] :"<< m_data[nv][ii].size()<<std::endl;
                    }
                    std::cout<<"\n"<<std::endl;
                    
                }
            } else {
                size_t num_threads = mesh_instancefm->getNthreads();
                size_t blk_size = vTx.size()/num_threads;
                if ( blk_size == 0 ) blk_size++;
                
                std::vector<std::thread> threads(num_threads-1);
                size_t blk_start = 0;
                for ( size_t i=0; i<num_threads-1; ++i ) {
                    
                    size_t blk_end = blk_start + blk_size;
                    Grid3Duifm<double , uint32_t> *mesh_ref = mesh_instancefm;
                    threads[i]=std::thread( [&mesh_ref,&vTx,&tt,&t0,&Rx,&iTx,&nRx,
                                             &r_data,&v0,&m_data,blk_start,blk_end,i]{
                        
                        for ( size_t nv=blk_start; nv<blk_end; ++nv ) {
                            
                            std::vector<sxyz<double >> vRx;
                            for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                                vRx.push_back( Rx[ iTx[nv][ni] ] );
                            }
                            if ( mesh_ref->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], i+1) == 1 ) {
                                throw std::runtime_error("Problem while raytracing.");
                            }
                        }
                    });
                    
                    blk_start = blk_end;
                }
                
                for ( size_t nv=blk_start; nv<vTx.size(); ++nv ) {
                    std::vector<sxyz<double >> vRx;
                    for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                        vRx.push_back( Rx[ iTx[nv][ni] ] );
                    }
                    if ( mesh_instancefm->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], 0) == 1 ) {
                        throw std::runtime_error("Problem while raytracing.");
                    }
                }
                
                std::for_each(threads.begin(),threads.end(),
                              std::mem_fn(&std::thread::join));
            }
            
            for (size_t i=0;i<r_data.size();++i){
                std::string filename="Src"+to_string(N+1+i)+".vtp";
                saveRayPaths(filename, r_data[i]);
            }
            
            cout<<"\n"<<endl;
            mesh_instancefm->saveTT("time",0,0,false);
        }
///////////////////////////////////////////
    std::cout << "Hello, World!\n";
    return 0;
}

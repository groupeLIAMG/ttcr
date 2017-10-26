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
#include "Rcv.h"
#include "Src.h"
#include "utils.h"
#include "Node.h"
using namespace ttcr;



int main(int argc, const char * argv[]) {


    std::vector<Src<double >> src;
    size_t N=0;
//    for(size_t i=0;i<15;i++){
//
//        src.push_back( Src<double >("Src"+to_string(i+1+N)+".dat") );
//        src[i].init();
//    }
    src.push_back( Src<double >("Src1.dat") );
    src[0].init();
    
    bool SecondNodes =false;
    std::vector<sxyz<double >> Ttx;
    std::vector<sxyz<double >> s;
    for (size_t i=0;i<src.size();i++){
        s=src[i].get_coord();
        Ttx.push_back(s.back());
    }
    std::vector<Rcv<double >> R;
    R.push_back( Rcv<double >("rcv1.dat") );
    //Rx.init(1);
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
    ///////////////////////////////////////////////////////////////////
    // find mesh "corners"
    double xmin = nodes[0].x;
    double xmax = nodes[0].x;
    double ymin = nodes[0].y;
    double ymax = nodes[0].y;
    double zmin = nodes[0].z;
    double zmax = nodes[0].z;
    for ( size_t n=1; n<nodes.size(); ++n ) {
        xmin = xmin < nodes[n].x ? xmin : nodes[n].x;
        xmax = xmax > nodes[n].x ? xmax : nodes[n].x;
        ymin = ymin < nodes[n].y ? ymin : nodes[n].y;
        ymax = ymax > nodes[n].y ? ymax : nodes[n].y;
        zmin = zmin < nodes[n].z ? zmin : nodes[n].z;
        zmax = zmax > nodes[n].z ? zmax : nodes[n].z;
    }
    // use corners as ref pts
    std::vector<sxyz<double >> refPts;
    refPts.push_back( {xmin, ymin, zmin} );
    refPts.push_back( {xmin, ymin, zmax} );
    refPts.push_back( {xmin, ymax, zmin} );
    refPts.push_back( {xmin, ymax, zmax} );
    refPts.push_back( {xmax, ymin, zmin} );
    refPts.push_back( {xmax, ymin, zmax} );
    refPts.push_back( {xmax, ymax, zmin} );
    refPts.push_back( {xmax, ymax, zmax} );
    double eps(1.e-20);
    std::vector<sxyz<double>> no;
    for (size_t n=0;n<nodes.size();++n){
        no.push_back({nodes[n].x-xmin,nodes[n].y-ymin,nodes[n].z-zmin});
    }
    std::vector<sxyz<double >> refPts2;
    for (size_t rr=0;rr<refPts.size();++rr)
        refPts2.push_back(refPts[rr]-refPts[0]);
    Grid3Duifs<double , uint32_t>* mesh_instance;
    if (SecondNodes==true){
       mesh_instance = new Grid3Duifs<double , uint32_t>(no, tetrahedra, eps,20,true, 1);
    }else{
        std::vector<sxyz<double >> refPts2;
        for (size_t rr=0;rr<refPts.size();++rr)
            refPts2.push_back(refPts[rr]-refPts[0]);
       mesh_instance = new Grid3Duifs<double , uint32_t>(no, tetrahedra, eps,20,refPts2,2,true, 1);
    }
    mesh_instance->setSlowness(tmp);
    //mesh_instance->setSourceRadius(20.0*1.e-3);
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
    
    if ( mesh_instance->getNthreads() == 1 || mesh_instance->getNthreads()<= vTx.size() ) {
        for ( size_t nv=0; nv<vTx.size(); ++nv ) {
            vRx.resize( 0 );
            for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                vRx.push_back( Rx[ iTx[nv][ni] ]-refPts[0]);
            }
            vTx[nv][0]=vTx[nv][0]-refPts[0];
            if (SecondNodes){
                // add secondary nodes
                std::vector<tetrahedronElem<uint32_t>> DelatedCellsTx;
                mesh_instance->addNodes(vTx[nv][0],DelatedCellsTx,mesh_instance->getNthreads());
                mesh_instance->addNodes(vTx[nv][0],DelatedCellsTx,mesh_instance->getNthreads());
                std::vector<sxyz<double >> refPts2;
                for (size_t rr=0;rr<refPts.size();++rr)
                    refPts2.push_back(refPts[rr]-refPts[0]);
                mesh_instance->initOrdering(refPts2,2);
//                if ( mesh_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], m_data[nv], 0) == 1 )
//                    throw std::runtime_error("Problem while raytracing.");
                if ( mesh_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], 0) == 1 )
                    throw std::runtime_error("Problem while raytracing.");
                mesh_instance->delateCells(DelatedCellsTx);
            }else{
                if ( mesh_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], m_data[nv], 0) == 1 )
                    throw std::runtime_error("Problem while raytracing.");
//                if ( mesh_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], 0) == 1 )
//                    throw std::runtime_error("Problem while raytracing.");
            }
            
            for (size_t ii=0;ii<r_data[nv].size();ii++){
                std::cout<<" size r_data["<<ii<<"] :"<< r_data[nv][ii].size()<<std::endl;
                if (r_data[nv][ii].size()>498)
                    std::cout<<"Infinit loop"<<std::endl;
            }
            std::cout<<"\n"<<std::endl;
            for (size_t ii=0;ii<m_data[nv].size();ii++){
                //for (size_t jj=0;jj<m_data[nv][ii].size();++jj)
                std::cout<<" size m_data["<<ii<<"] :"<< m_data[nv][ii].size()<<std::endl;
            }
             std::cout<<"\n"<<std::endl;

        }
    } else {
        size_t num_threads = mesh_instance->getNthreads();
        size_t blk_size = vTx.size()/num_threads;
        if ( blk_size == 0 ) blk_size++;

        std::vector<std::thread> threads(num_threads-1);
        size_t blk_start = 0;
        for ( size_t i=0; i<num_threads-1; ++i ) {

            size_t blk_end = blk_start + blk_size;
            Grid3Duifs<double , uint32_t> *mesh_ref = mesh_instance;
            threads[i]=std::thread( [&mesh_ref,&vTx,&tt,&t0,&Rx,&iTx,&nRx,
                                     &r_data,&v0,&m_data,blk_start,blk_end,i]{

                for ( size_t nv=blk_start; nv<blk_end; ++nv ) {

                    std::vector<sxyz<double >> vRx;
                    for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                        vRx.push_back( Rx[ iTx[nv][ni] ] );
                    }
                    if ( mesh_ref->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], m_data[nv], i+1) == 1 ) {
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
            if ( mesh_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], m_data[nv], 0) == 1 ) {
                throw std::runtime_error("Problem while raytracing.");
            }
        }

        std::for_each(threads.begin(),threads.end(),
                      std::mem_fn(&std::thread::join));
    }
    double traveltimes [nRx];
    double V0 [nRx];
    for ( size_t nv=0; nv<vTx.size(); ++nv ) {
        for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
            traveltimes[ iTx[nv][ni] ] = tt[nv][ni];
            V0[ iTx[nv][ni] ] = v0[nv];
        }
    }
    for (size_t nv=0;nv<r_data.size();++nv){
        for(size_t i=0;i<r_data[nv].size();++i){
            for(size_t P=0;P<r_data[nv][i].size();++P)
            r_data[nv][i][P]=r_data[nv][i][P]+refPts[0];
        }
    }


        for (size_t i=0;i<r_data.size();++i){
            std::string filename="Src"+to_string(N+1+i)+".vtp";
            saveRayPaths(filename, r_data[i]);
        }

    
    std::cout << "Hello, World!\n";
    return 0;
}

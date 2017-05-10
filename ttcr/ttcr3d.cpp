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
    


    std::vector<Src<long double >> src;
    src.push_back( Src<long double >("Src12.dat") );
    src[0].init();
//    src.push_back( Src<long double >("Src2.dat") );
//    src[1].init();
//    src.push_back( Src<long double >("Src3.dat") );
//    src[2].init();
//    src.push_back( Src<long double >("Src4.dat") );
//    src[3].init();
//    src.push_back( Src<long double >("Src5.dat") );
//    src[4].init();
//    src.push_back( Src<long double >("Src6.dat") );
//    src[5].init();
//    src.push_back( Src<long double >("Src7.dat") );
//    src[6].init();
//    src.push_back( Src<long double >("Src8.dat") );
//    src[7].init();
//    src.push_back( Src<long double >("Src9.dat") );
//    src[8].init();
//    src.push_back( Src<long double >("Src10.dat") );
//    src[9].init();
    std::vector<sxyz<long double >> Ttx;
    std::vector<sxyz<long double >> s;
    for (size_t i=0;i<src.size();i++){
        s=src[i].get_coord();
        Ttx.push_back(s.back());
    }
    std::vector<Rcv<long double >> R;
    R.push_back( Rcv<long double >("rcv.dat") );
    //Rx.init(1);
    R[0].init(src.size());
    std::vector<sxyz<long double >> Rrx;
    for(size_t i=0;i<R.size();i++){
        for(size_t j=0;j<R[i].get_coord().size();j++){
            Rrx.push_back(R[i].get_coord()[j]);
        }
    }
    size_t nTtx=Ttx.size();
    size_t nRrx=Rrx.size();
    std::vector<sxyz<long double >> Tx(nTtx*nRrx);
    std::vector<sxyz<long double >> Rx(nTtx*nRrx);
    for (unsigned i=0; i<nTtx;i++){
        for(size_t j=0; j<nRrx;j++){
            Tx[i*nRrx+j]=Ttx[i];
            Rx[i*nRrx+j]=Rrx[j];
        }
    }
    MSHReader reader( "Model1.msh" );
    std::vector<sxyz<double >> node(reader.getNumberOfNodes());
    std::vector<tetrahedronElem<uint32_t>> tetrahedra(reader.getNumberOfTetra());
    std::vector<long double > slowness(reader.getNumberOfNodes());

    reader.readNodes3D(node);
    std::vector <sxyz<long double>> nodes;
    for ( size_t n=0; n<node.size(); ++n ) {
        nodes.push_back(sxyz<long double>(node[n].x,node[n].y,node[n].z));
    }
    sxyz<long double> PP(nodes[0]);
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
    std::vector<long double> tmp(tmp1.begin(),tmp1.end());
    fin.close();

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
    std::vector<sxyz<long double >> refPts;
    refPts.push_back( {xmin, ymin, zmin} );
    refPts.push_back( {xmin, ymin, zmax} );
    refPts.push_back( {xmin, ymax, zmin} );
    refPts.push_back( {xmin, ymax, zmax} );
    refPts.push_back( {xmax, ymin, zmin} );
    refPts.push_back( {xmax, ymin, zmax} );
    refPts.push_back( {xmax, ymax, zmin} );
    refPts.push_back( {xmax, ymax, zmax} );
    double eps(1.e-15);
    for (size_t n=0;n<nodes.size();++n){
        nodes[n]=nodes[n]-refPts[0];
    }
    for (size_t n=0;n<Rx.size();++n){
        Rx[n]=Rx[n]-refPts[0];
    }
    for (size_t n=0;n<Tx.size();++n){
        Tx[n]=Tx[n]-refPts[0];
    }
    Grid3Duifs<long double , uint32_t>*mesh_instance = new Grid3Duifs<long double , uint32_t>(nodes, tetrahedra, eps,20,true, 1);
    mesh_instance->setSlowness(tmp);
    mesh_instance-> setSourceRadius(0.0);

//
//   std::vector<std::vector<tetrahedronElem<uint32_t>>> DelatedCellsRx(Rrx.size());
//    for(size_t r=0;r<Rrx.size();++r){
//        mesh_instance->addNodes(Rrx[r]-refPts[0],DelatedCellsRx[r],mesh_instance->getNthreads());
////        mesh_instance->addNodes(Rrx[r]-refPts[0],DelatedCellsRx[r],mesh_instance->getNthreads());
////        mesh_instance->addNodes(Rrx[r]-refPts[0],DelatedCellsRx[r],mesh_instance->getNthreads());
//    }
    //////////////

    size_t nTx = Tx.size();
    size_t nRx=Rx.size();
    std::vector<std::vector<sxyz<long double >>> vTx;
    std::vector<std::vector<long double >> t0;
    std::vector<std::vector<size_t>> iTx;
    vTx.push_back( std::vector<sxyz<long double >>(1,Tx[0]) );
    t0.push_back(std::vector<long double >(1,0) );
    iTx.push_back(std::vector<size_t>(1, 0) );  // indices of Rx corresponding to current Tx
    for ( size_t ntx=1; ntx<nTx; ++ntx ) {
        bool found = false;

        for ( size_t nv=0; nv<vTx.size(); ++nv ) {
            if ( vTx[nv][0].getDistance(Tx[ntx])<0.000001) {
                found = true;
                iTx[nv].push_back( ntx ) ;
                break;
            }
        }
        if ( !found ) {
            vTx.push_back( std::vector<sxyz<long double >>(1, Tx[ntx]) );
            t0.push_back( std::vector<long double >(1,0) );
            iTx.push_back(std::vector<size_t>(1, ntx) );
        }
    }
    /*
     Looping over all non redundant Tx
     */
//    std::vector<std::vector<sijv<long double>>> d_data;
//    mesh_instance->computeD(Rx, d_data);

    std::vector<sxyz<long double >> vRx;
    std::vector<std::vector<long double >> tt( vTx.size() );
    std::vector<std::vector<std::vector<sxyz<long double >>>> r_data( vTx.size() );
    std::vector<long double > v0( vTx.size() );
    std::vector<std::vector<std::vector<sijv<long double >>>> m_data( vTx.size() );
  ///////////////////////////// test triangle
//    Node3Di<long double, uint32_t> vertexA(mesh_instance->getNthreads());
//    vertexA.setXYZindex(nodes[1455], 1);
//    Node3Di<long double, uint32_t> vertexB(mesh_instance->getNthreads());
//    vertexB.setXYZindex(nodes[5010], 2);
//    Node3Di<long double, uint32_t> vertexC(mesh_instance->getNthreads());
//    vertexC.setXYZindex(nodes[11128], 3);
//    bool v=mesh_instance->testInTriangle(&vertexA, &vertexB, &vertexC, nodes[3489]);
    
//    std::array<uint32_t, 3> faces={1455,5010,11128};
//    std::set<uint32_t> Cells;
//    mesh_instance->findAdjacentCell3(faces, 28734, Cells);
    
    if ( mesh_instance->getNthreads() == 1 || mesh_instance->getNthreads()<= vTx.size() ) {
        for ( size_t nv=0; nv<vTx.size(); ++nv ) {
            vRx.resize( 0 );
            for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                vRx.push_back( Rx[ iTx[nv][ni] ] );
            }
           //  add new nodes arround  sources ans receivers
            std::vector<tetrahedronElem<uint32_t>> DelatedCellsTx;
            mesh_instance->addNodes(vTx[nv][0],DelatedCellsTx,mesh_instance->getNthreads());
            //mesh_instance->addNodes(vTx[nv][0],DelatedCellsTx,mesh_instance->getNthreads());
            mesh_instance->initOrdering(refPts,2);
            if ( mesh_instance->raytrace(vTx[nv], t0[nv], vRx, tt[nv], r_data[nv], v0[nv], m_data[nv]) == 1 ) {
                throw std::runtime_error("Problem while raytracing.");
            }
//            mesh_instance->delateCells(DelatedCellsTx);
            //            cout<< "size r_data"<< r_data[nv].size()<<endl;
            //            for (size_t ii=0;ii<r_data[nv].size();ii++){
            //                cout<<" size r_data ["<< ii<<"]= "<<r_data[nv][ii].size()<<endl;
            //            }
            cout<<m_data[nv].size()<<endl;
            for (size_t ii=0;ii<r_data[nv].size();ii++){
                cout<<" size r_data["<<ii<<"] :"<< r_data[nv][ii].size()<<endl;
                //                 cout<<"t0 ("<<nv<<","<<ii<<") = "<< t0[nv][ii]<<endl;
            }
            //            for (size_t i=0;i<tt[nv].size();i++){
            //                //               cout<<" size r_data["<<ii<<"] :"<< t0[nv].size()<<endl;
            //                cout<<"tt ("<<nv<<","<<i<<") = "<< tt[nv][i]<<endl;
            //            }
        }
    } else {
        size_t num_threads = mesh_instance->getNthreads();
        size_t blk_size = vTx.size()/num_threads;
        if ( blk_size == 0 ) blk_size++;

        std::vector<std::thread> threads(num_threads-1);
        size_t blk_start = 0;
        for ( size_t i=0; i<num_threads-1; ++i ) {

            size_t blk_end = blk_start + blk_size;
            Grid3Duifs<long double , uint32_t> *mesh_ref = mesh_instance;
            threads[i]=std::thread( [&mesh_ref,&vTx,&tt,&t0,&Rx,&iTx,&nRx,
                                     &r_data,&v0,&m_data,blk_start,blk_end,i]{

                for ( size_t nv=blk_start; nv<blk_end; ++nv ) {

                    std::vector<sxyz<long double >> vRx;
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
            std::vector<sxyz<long double >> vRx;
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
        size_t half_size (r_data[i].size()/2);
       std::vector<std::vector<sxyz<long double >>> r_data1(r_data[i].begin(),r_data[i].begin()+half_size);
       std::vector<std::vector<sxyz<long double >>> r_data2(r_data[i].begin()+half_size,r_data[i].end());
        std::string filename="Sourcea.vtp";//+to_string(i+11)+".vtp";
        saveRayPaths(filename, r_data1);
        filename="Sourceb.vtp";//+to_string(i+11)+".vtp";
        saveRayPaths(filename, r_data2);
    }
//
//        for (size_t i=0;i<r_data.size();++i){
//            std::string filename="Source"+to_string(i+1)+".vtp";
//            saveRayPaths(filename, r_data[i]);
//        }
    // insert code here...
    std::cout << "Hello, World!\n";
    return 0;
}

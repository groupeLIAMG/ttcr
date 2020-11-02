//
//  grid3duifm_mex.cpp
//  ttcr
//
//  Created by Bernard Giroux on 2015-03-04.
//  Copyright (c) 2015 Bernard Giroux. All rights reserved.
//

#include <exception>
#include <thread>

#include "mex.h"
#include "class_handle.hpp"

#include "Grid3Dunfs.h"

using namespace std;
using namespace ttcr;

typedef Grid3D<double,uint32_t> g3d;
typedef Grid3Dunfs<double,uint32_t> grid;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Get the command string
    char cmd[64];
    if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
        mexErrMsgTxt("First input should be a command string less than 64 characters long.");

    //  ---------------------------------------------------------------------------
    // New
    if (!strcmp("new", cmd)) {
        // Check parameters
        if (nlhs != 1) {
            mexErrMsgTxt("New: One output expected.");
        }
        if (nrhs != 3 || nrhs != 4 ) {
            mexErrMsgTxt("New: 2 or 3 input arguments needed.");
        }
        // Return a handle to a new C++ instance

        //
        // Nodes
        //
        if (!(mxIsDouble(prhs[1]))) {
            mexErrMsgTxt("Node coordinates must be double precision.");
        }
        double *xyz = static_cast<double*>( mxGetPr(prhs[1]) );
        mwSize number_of_dims = mxGetNumberOfDimensions(prhs[1]);
        if ( number_of_dims != 2 ) {
            mexErrMsgTxt("Node coordinates must be a matrix (nNodes by 3).");
        }
        const mwSize *dim_array = mxGetDimensions(prhs[1]);
        size_t nXYZ = static_cast<size_t>( dim_array[0] );
        if ( dim_array[1] != 3 ) {
            mexErrMsgTxt("Node coordinates must be a matrix (nNodes by 3).");
        }
        vector<sxyz<double>> nodes(nXYZ);
        for ( size_t n=0; n<nXYZ; ++n ) {
            nodes[n].x = xyz[n];
            nodes[n].y = xyz[n+nXYZ];
            nodes[n].z = xyz[n+2*nXYZ];
        }

        //
        // Tetrahedra
        //
        if (!(mxIsDouble(prhs[2]))) {
            mexErrMsgTxt("Tetrahedra node indices must be double precision.");
        }
        double *ind = static_cast<double*>( mxGetPr(prhs[2]) );
        number_of_dims = mxGetNumberOfDimensions(prhs[2]);
        if ( number_of_dims != 2 ) {
            mexErrMsgTxt("Tetrahedra node indices coordiates must be a matrix (nTetrahedra by 4).");
        }
        dim_array = mxGetDimensions(prhs[2]);
        size_t nTet = static_cast<size_t>( dim_array[0] );
        if ( dim_array[1] != 4 ) {
            mexErrMsgTxt("Tetrahedra node indices coordiates must be a matrix (nTetrahedra by 4).");
        }
        vector<tetrahedronElem<uint32_t>> tetrahedra(nTet);
        for ( size_t n=0; n<nTet; ++n ) {
            tetrahedra[n].i[0] = static_cast<uint32_t>( ind[n]-1 );
            tetrahedra[n].i[1] = static_cast<uint32_t>( ind[n+nTet]-1 );
            tetrahedra[n].i[2] = static_cast<uint32_t>( ind[n+2*nTet]-1 );
            tetrahedra[n].i[3] = static_cast<uint32_t>( ind[n+3*nTet]-1 );
        }

		// ------------------------------------------------------
		// number of threads
		// ------------------------------------------------------
		size_t nthreads = 1;
		if ( nrhs == 4 ) {
			size_t mrows = mxGetM(prhs[3]);
			size_t ncols = mxGetN(prhs[3]);
			if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
			   !(mrows==1 && ncols==1) ) {
				mexErrMsgTxt("Input must be a noncomplex scalar double.");
			}

			double *dtmp = mxGetPr( prhs[3] );
			nthreads = round( *dtmp );
		}

		int order=2; // use order 2 for metric
        double xmin = nodes[0].x;
        double xmax = nodes[0].x;
        double ymin = nodes[0].y;
        double ymax = nodes[0].y;
        double zmin = nodes[0].z;
        double zmax = nodes[0].z;
        for (size_t n=1; n<nodes.size(); ++n) {
            xmin = (xmin<nodes[n].x) ? xmin : nodes[n].x;
            xmax = (xmax>nodes[n].x) ? xmax : nodes[n].x;
            ymin = (ymin<nodes[n].y) ? ymin : nodes[n].y;
            ymax = (ymax>nodes[n].y) ? ymax : nodes[n].y;
            zmin = (zmin<nodes[n].z) ? zmin : nodes[n].z;
            zmax = (zmax>nodes[n].z) ? zmax : nodes[n].z;
        }

        std::vector<sxyz<double>> ptsRef;
        ptsRef.push_back( {xmin, ymin, zmin} );
        ptsRef.push_back( {xmin, ymin, zmax} );
        ptsRef.push_back( {xmin, ymax, zmin} );
        ptsRef.push_back( {xmin, ymax, zmax} );
        ptsRef.push_back( {xmax, ymin, zmin} );
        ptsRef.push_back( {xmax, ymin, zmax} );
        ptsRef.push_back( {xmax, ymax, zmin} );
        ptsRef.push_back( {xmax, ymax, zmax} );

        int rp_method = 1;
        bool invert_vel = false;
        bool rp_from_tt = true;
        double min_dist = 1.e-8;

        plhs[0] = convertPtr2Mat<g3d>(new grid(nodes, tetrahedra, 1.e-15, 20,
                                                ptsRef, order, rp_method,
                                                invert_vel, rp_from_tt,
                                                min_dist, nthreads));
        return;
    }

    // Check there is a second input, which should be the class instance handle
    if (nrhs < 2)
        mexErrMsgTxt("Second input should be a class instance handle.");

    // ---------------------------------------------------------------------------
    // Delete
    //
    if (!strcmp("delete", cmd)) {
        // Destroy the C++ object
        destroyObject<g3d>(prhs[1]);
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
    }

    // Get the class instance pointer from the second input
    g3d *grid_instance = convertMat2Ptr<g3d>(prhs[1]);



    // Call the various class methods
    // ---------------------------------------------------------------------------
    // setSlowness
    //
    if (!strcmp("setSlowness", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs != 3)
            mexErrMsgTxt("setSlowness: Unexpected arguments.");
        // Call the method

        if (!(mxIsDouble(prhs[2]))) {
            mexErrMsgTxt("Slowness must be double precision.");
        }
        double *slowness = static_cast<double*>( mxGetPr(prhs[2]) );
        mwSize number_of_dims = mxGetNumberOfDimensions(prhs[2]);
        if ( number_of_dims != 2 ) {
            mexErrMsgTxt("Slowness must be a vector (nSlowness by 1).");
        }
        const mwSize *dim_array = mxGetDimensions(prhs[2]);
        size_t nSlowness = static_cast<size_t>( dim_array[0] );
        if ( dim_array[1] != 1 ) {
            mexErrMsgTxt("Slowness must be a vector (nSlowness by 1).");
        }

        try {
            grid_instance->setSlowness(slowness, nSlowness);
        } catch (std::exception& e) {
            mexErrMsgTxt("Slowness values must be defined for each grid node.");
        }

        return;
    }

    //  ---------------------------------------------------------------------------
    // raytrace
    if (!strcmp("raytrace", cmd)) {
        // Check parameters

        if ( nrhs != 5 && nrhs != 6 ) {
            mexErrMsgTxt("raytrace: Unexpected arguments.");
        }
        if (nlhs > 3) {
            mexErrMsgTxt("raytrace has a maximum of three output argument.");
        }
        //
        // Slowness
        //
        if (!(mxIsDouble(prhs[2]))) {
            mexErrMsgTxt("Slowness must be double precision.");
        }
        double *slowness = static_cast<double*>( mxGetPr(prhs[2]) );
        mwSize number_of_dims = mxGetNumberOfDimensions(prhs[2]);
        if ( number_of_dims != 2 ) {
            mexErrMsgTxt("Slowness must be a vector (nSlowness by 1).");
        }
        const mwSize *dim_array = mxGetDimensions(prhs[2]);
        size_t nSlowness = static_cast<size_t>( dim_array[0] );
        if ( dim_array[1] != 1 ) {
            mexErrMsgTxt("Slowness must be a vector (nSlowness by 1).");
        }

        try {
            grid_instance->setSlowness(slowness, nSlowness);
        } catch (std::exception& e) {
            mexErrMsgTxt("Slowness values must be defined for each grid node.");
        }

        //
        // Tx
        //
        if (!(mxIsDouble(prhs[3]))) {
            mexErrMsgTxt("Tx must be double precision.");
        }
        number_of_dims = mxGetNumberOfDimensions(prhs[3]);
        if ( number_of_dims != 2 ){
            mexErrMsgTxt("Tx must be a rank 2 matrix.");
        }
        dim_array = mxGetDimensions(prhs[3]);
        size_t nTx = static_cast<size_t>( dim_array[0] );
        if ( dim_array[1] != 3 ) {
            mexErrMsgTxt("Tx: matrix nTx by 3.");
        }
        double *Tx = static_cast<double*>( mxGetPr(prhs[3]) );

        //
        // Rx
        //
        if (!(mxIsDouble(prhs[4]))) {
            mexErrMsgTxt("Rx must be double precision.");
        }
        number_of_dims = mxGetNumberOfDimensions(prhs[4]);
        if ( number_of_dims != 2 ){
            mexErrMsgTxt("Rx must be a rank 2 matrix.");
        }
        dim_array = mxGetDimensions(prhs[4]);
        size_t nRx = static_cast<size_t>( dim_array[0] );
        if ( dim_array[1] != 3 ) {
            mexErrMsgTxt("Rx: matrix nRx by 3.");
        }
        double *Rx = static_cast<double*>( mxGetPr(prhs[4]) );

        if ( nTx != nRx ) {
            mexErrMsgTxt("nTx should be equal to nRx.");
        }

        //
        // t0
        //
        double *tTx;
        if ( nrhs == 6 ) {
            if (!(mxIsDouble(prhs[5]))) {
                mexErrMsgTxt("t0 must be double precision.");
            }
            number_of_dims = mxGetNumberOfDimensions(prhs[5]);
            if ( number_of_dims != 2 ){
                mexErrMsgTxt("t0 must be a rank 2 matrix.");
            }
            dim_array = mxGetDimensions(prhs[5]);
            size_t nT0 = static_cast<size_t>( dim_array[0] );
            if ( dim_array[1] != 1 || nT0 != nTx ) {
                mexErrMsgTxt("t0: matrix nTx by 1.");
            }
            tTx = static_cast<double*>( mxGetPr(prhs[5]) );
        } else {
            tTx = new double [nTx];
            for ( size_t n=0; n<nTx; ++n ) tTx[n] = 0.0;
        }

        /* ------------------------------------------------------
         Output variable
         ------------------------------------------------------ */

        plhs[0] = mxCreateDoubleMatrix(nRx, 1, mxREAL);
        double *t_arr = mxGetPr(plhs[0]);

        /*
         Looking for redundants Tx pts
         */
        vector<vector<sxyz<double>>> vTx;
        vector<vector<double>> t0;
        vector<vector<size_t>> iTx;
        sxyz<double> sxyz_tmp;
        sxyz_tmp.x = Tx[0];
        sxyz_tmp.y = Tx[nTx];
        sxyz_tmp.z = Tx[2*nTx];
        vTx.push_back( vector<sxyz<double> >(1, sxyz_tmp) );
        t0.push_back( vector<double>(1, tTx[0]) );
        iTx.push_back( vector<size_t>(1, 0) );  // indices of Rx corresponding to current Tx
        for ( size_t ntx=1; ntx<nTx; ++ntx ) {
            sxyz_tmp.x = Tx[ntx];
            sxyz_tmp.y = Tx[ntx+nTx];
            sxyz_tmp.z = Tx[ntx+2*nTx];
            bool found = false;

            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                if ( vTx[nv][0].x==sxyz_tmp.x &&
                    vTx[nv][0].y==sxyz_tmp.y &&
                    vTx[nv][0].z==sxyz_tmp.z ) {
                    found = true;
                    iTx[nv].push_back( ntx ) ;
                    break;
                }
            }
            if ( !found ) {
                vTx.push_back( vector<sxyz<double>>(1, sxyz_tmp) );
                t0.push_back( vector<double>(1, tTx[ntx]) );
                iTx.push_back( vector<size_t>(1, ntx) );
            }
        }

        if ( nrhs == 5 ) {
            delete [] tTx;
        }

        /* ------------------------------------------------------
         Optional output variables
         ------------------------------------------------------ */

        mxArray **Rays;
        mxArray **M;
        if ( nlhs >= 2 ) {
            // 2rd arg: rays.
            plhs[1] = mxCreateCellMatrix(nRx, 1);
            Rays = (mxArray **) mxCalloc(nRx, sizeof(mxArray *));
        }
        if ( nlhs >= 3 ) {
            plhs[2] = mxCreateCellMatrix(vTx.size(), 1);
            M = (mxArray **) mxCalloc(vTx.size(), sizeof(mxArray *));
        }

        /*
         Looping over all non redundant Tx
         */

        vector<vector<sxyz<double>>> vRx( vTx.size() );
        vector<vector<double>> tt( vTx.size() );
        vector<vector<vector<sxyz<double>>>> r_data( vTx.size() );
        vector<vector<vector<sijv<double>>>> m_data( vTx.size() );

        for ( size_t nv=0; nv<vTx.size(); ++nv ) {

            vRx[nv].resize( 0 );
            for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                sxyz_tmp.x = Rx[ iTx[nv][ni] ];
                sxyz_tmp.y = Rx[ iTx[nv][ni]+nRx ];
                sxyz_tmp.z = Rx[ iTx[nv][ni]+2*nRx ];
                vRx[nv].push_back( sxyz_tmp );
            }
        }

        if ( nlhs == 3 ) {
            try {
                grid_instance->raytrace(vTx, t0, vRx, tt, r_data, m_data);
            } catch (...) {
                mexErrMsgTxt("Problem while raytracing.");
            }
        } else if ( nlhs == 2 ) {
            try {
                grid_instance->raytrace(vTx, t0, vRx, tt, r_data);
            } catch (...) {
                mexErrMsgTxt("Problem while raytracing.");
            }
        }
        else {
            try {grid_instance->raytrace(vTx, t0, vRx, tt);
            } catch (...) {
                mexErrMsgTxt("Problem while raytracing.");
            }
        }

        for ( size_t nv=0; nv<vTx.size(); ++nv ) {
            for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                t_arr[ iTx[nv][ni] ] = tt[nv][ni];
            }
        }

        if ( nlhs >= 2 ) {
            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                for ( size_t ni=0; ni<iTx[nv].size(); ++ni ) {
                    size_t npts = r_data[nv][ni].size();
                    Rays[ iTx[nv][ni] ] = mxCreateDoubleMatrix(npts, 3, mxREAL);
                    double *rays_p = (double*) mxGetData(Rays[ iTx[nv][ni] ]);
                    for ( size_t np=0; np<npts; ++np ) {
                        rays_p[np] = r_data[nv][ni][np].x;
                        rays_p[np+npts] = r_data[nv][ni][np].y;
                        rays_p[np+2*npts] = r_data[nv][ni][np].z;
                    }
                    mxSetCell( plhs[1], iTx[nv][ni], Rays[ iTx[nv][ni] ] );
                }
            }
        }
        if ( nlhs >= 3 ) {
            // TODO : remove elements for static corrections
            // for this to work, Tx & Rx data should be ordered so that redundant Tx should be contiguous
            for ( size_t nv=0; nv<vTx.size(); ++nv ) {
                size_t nRcv = m_data[nv].size();
                size_t nMmax = nRcv;
                for ( size_t ni=0; ni<m_data[nv].size(); ++ni ) {
                    nMmax += m_data[nv][ni].size();
                }
                M[ nv ] = mxCreateSparse(nRcv, nSlowness+nRcv, nMmax, mxREAL);
                double *Mval = mxGetPr( M[ nv ] );
                mwIndex *irM  = mxGetIr( M[ nv ] );
                mwIndex *jcM  = mxGetJc( M[ nv ] );
                size_t k = 0;
                for ( size_t j=0; j<nSlowness; ++j ) {
                    jcM[j] = k;
                    for ( size_t ni=0; ni<m_data[nv].size(); ++ni ) {
                        for ( size_t n=0; n<m_data[nv][ni].size(); ++n) {
                            if ( m_data[nv][ni][n].j == j && m_data[nv][ni][n].i == ni ) {
                                irM[k] = ni;
                                Mval[k] = m_data[nv][ni][n].v;
                                k++;
                            }
                        }
                    }
                }
                for ( size_t j=0; j<nRcv; ++j ) {  // derivative of t w/r to static correction
                    jcM[nSlowness+j] = k;
                    irM[k] = j;
                    Mval[k] = 1.0;
                    k++;
                }
                jcM[nSlowness+nRcv] = k;
                mxSetCell( plhs[2], nv, M[ nv ] );
            }
        }
        return;
    }

    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}

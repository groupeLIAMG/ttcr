//
//  Mesh3Dttcr.h
//  ttcr
//
//  Created by Bernard Giroux on 16-10-11.
//  Copyright © 2016 Bernard Giroux. All rights reserved.
//

#ifndef Mesh3Dttcr_h
#define Mesh3Dttcr_h

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "Python.h"
#include "numpy/ndarrayobject.h"

#include <vector>

#include "Grid3Dunfs.h"

namespace ttcr {
    typedef Grid3Dunfs<double,uint32_t> mesh;

    class Mesh3Dttcr {
    public:
        Mesh3Dttcr(const std::vector<sxyz<double>>&,
                   const std::vector<tetrahedronElem<uint32_t>>&,
                   const double, const int, const bool,
                   const size_t);
        
        ~Mesh3Dttcr() {
            delete mesh_instance;
        }
        
        void setSlowness(const std::vector<double>& slowness);
        
        void raytrace(const std::vector<sxyz<double>>& Tx,
                     const std::vector<double>& tTx,
                     const std::vector<sxyz<double>>& Rx,
                     double* traveltimes) const;

        void raytrace(const std::vector<sxyz<double>>& Tx,
                     const std::vector<double>& tTx,
                     const std::vector<sxyz<double>>& Rx,
                     double* traveltimes,
                     PyObject* rays,
                     double* v0) const;

//        void raytrace(const std::vector<sxyz<double>>& Tx,
//                     const std::vector<double>& tTx,
//                     const std::vector<sxyz<double>>& Rx,
//                     double* traveltimes,
//                     PyObject* rays,
//                     double* v0,
//                     PyObject* M) const;

    private:
        mesh *mesh_instance;
        
        Mesh3Dttcr() {}
    };
}

#endif /* Mesh3Dttcr_h */

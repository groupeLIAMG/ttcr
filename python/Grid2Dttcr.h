//
//  Grid2Dttcr.hpp
//  ttcr
//
//  Created by Bernard Giroux on 16-07-28.
//  Copyright © 2016 Bernard Giroux. All rights reserved.
//

#ifndef Grid2Dttcr_h
#define Grid2Dttcr_h



extern "C" {
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "Python.h"
#include "numpy/ndarrayobject.h"
}

#include <string>
#include <vector>

#include "Cell.h"
#include "Grid2Drcsp.h"


namespace ttcr {
    
    typedef Grid2D<double,uint32_t,sxz<double>> grid;
    typedef Grid2Drcsp<double,uint32_t,Cell<double,Node2Dcsp<double,uint32_t>,sxz<double>>> gridiso;
    typedef Grid2Drcsp<double,uint32_t,CellElliptical<double,Node2Dcsp<double,uint32_t>,sxz<double>>> gridaniso;
    typedef Grid2Drcsp<double,uint32_t,CellTiltedElliptical<double,Node2Dcsp<double,uint32_t>,sxz<double>>> gridtilted;
    
    class Grid2Dttcr {
    public:
        Grid2Dttcr(std::string&, uint32_t, uint32_t, double, double, double, double, uint32_t, uint32_t, size_t);
        ~Grid2Dttcr() {
            delete grid_instance;
        }
        std::string getType() const { return type; }
        
        void setSlowness(const std::vector<double>& slowness);
        void setXi(const std::vector<double>& xi);
        void setTheta(const std::vector<double>& theta);
        
        int raytrace(const std::vector<sxz<double>>& Tx,
                     const std::vector<double>& tTx,
                     const std::vector<sxz<double>>& Rx,
                     double* traveltimes,
                     PyObject* rays,
                     PyObject* L) const;
        
    private:
        const std::string type;
        grid *grid_instance;
        
        Grid2Dttcr() {}
    };
    
}

#endif /* Grid2Dttcr_h */

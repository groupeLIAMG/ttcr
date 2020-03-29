//
//  Grid2Dttcr.h
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef Grid2Dttcr_h
#define Grid2Dttcr_h

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "Python.h"
#include "numpy/ndarrayobject.h"

#include <string>
#include <vector>

#include "Cell.h"
#include "Grid2Drcsp.h"

namespace ttcr {

    typedef Grid2D<double,uint32_t,sxz<double>> grid;
    typedef Grid2Drcsp<double,uint32_t,sxz<double>,Cell<double,Node2Dcsp<double,uint32_t>,sxz<double>>> gridiso;
    typedef Grid2Drcsp<double,uint32_t,sxz<double>,CellElliptical<double,Node2Dcsp<double,uint32_t>,sxz<double>>> gridaniso;
    typedef Grid2Drcsp<double,uint32_t,sxz<double>,CellTiltedElliptical<double,Node2Dcsp<double,uint32_t>,sxz<double>>> gridtilted;

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

        void raytrace(const std::vector<sxz<double>>& Tx,
                      const std::vector<double>& tTx,
                      const std::vector<sxz<double>>& Rx,
                      double* traveltimes,
                      PyObject* rays,
                      PyObject* L) const;

        void raytrace(const std::vector<sxz<double>>& Tx,
                      const std::vector<double>& tTx,
                      const std::vector<sxz<double>>& Rx,
                      double* traveltimes,
                      PyObject* L) const;

        void raytrace(const std::vector<sxz<double>>& Tx,
                      const std::vector<double>& tTx,
                      const std::vector<sxz<double>>& Rx,
                      double* traveltimes) const;

        static int Lsr2d(const double* Tx,
                          const double* Rx,
                          const size_t nTx,
                          const double* grx,
                          const size_t n_grx,
                          const double* grz,
                          const size_t n_grz,
                          PyObject* L);

		static int Lsr2da(const double* Tx,
						   const double* Rx,
						   const size_t nTx,
						   const double* grx,
						   const size_t n_grx,
						   const double* grz,
						   const size_t n_grz,
						   PyObject* L);


    private:
        const std::string type;
        grid *grid_instance;

        Grid2Dttcr() {}
    };

}

#endif /* Grid2Dttcr_h */

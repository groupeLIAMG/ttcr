//
//  Grid3Dttcr.h
//  ttcr
//
//  Created by Bernard Giroux on 16-10-11.
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

#ifndef Grid3Dttcr_h
#define Grid3Dttcr_h

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "Python.h"
#include "numpy/ndarrayobject.h"

#include <string>
#include <vector>

#include "Cell.h"
#include "Grid3Drifs.h"
#include "Grid3Drcfs.h"

namespace ttcr {

    typedef Grid3D<double,uint32_t> grid;
    typedef Grid3Drifs<double,uint32_t> gridi;
    typedef Grid3Drcfs<double,uint32_t> gridc;

    class Grid3Dttcr {
    public:
        Grid3Dttcr(const std::string&, const uint32_t, const uint32_t, const uint32_t, const double,
                   const double, const double, const double,
                   const double, const int, const bool,
                   const size_t);

        ~Grid3Dttcr() {
            delete grid_instance;
        }
        std::string getType() const { return type; }

        void setSlowness(const std::vector<double>& slowness);

        int raytrace(const std::vector<sxyz<double>>& Tx,
                     const std::vector<double>& tTx,
                     const std::vector<sxyz<double>>& Rx,
                     double* traveltimes) const;

        int raytrace(const std::vector<sxyz<double>>& Tx,
                     const std::vector<double>& tTx,
                     const std::vector<sxyz<double>>& Rx,
                     double* traveltimes,
                     PyObject* rays,
                     double* v0) const;

        int raytrace(const std::vector<sxyz<double>>& Tx,
                     const std::vector<double>& tTx,
                     const std::vector<sxyz<double>>& Rx,
                     double* traveltimes,
                     PyObject* rays,
                     double* v0,
                     PyObject* M) const;

       int raytrace(const std::vector<sxyz<double>>& Tx,
                    const std::vector<double>& tTx,
                    const std::vector<sxyz<double>>& Rx,
                    double* traveltimes,
                    PyObject* rays,
                    PyObject* L) const;

       int raytrace(const std::vector<sxyz<double>>& Tx,
                    const std::vector<double>& tTx,
                    const std::vector<sxyz<double>>& Rx,
                    double* traveltimes,
                    PyObject* L) const;

    private:
        const std::string type;
        grid *grid_instance;

        Grid3Dttcr() {}
    };


}

#endif /* Grid3Dttcr_h */

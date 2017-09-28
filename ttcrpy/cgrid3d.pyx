# -*- coding: utf-8 -*-

"""
Copyright 2017 Bernard Giroux
email: bernard.giroux@ete.inrs.ca

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport uint32_t
from libcpp cimport bool

import numpy as np
cimport numpy as np

from scipy.sparse import csr_matrix

cdef extern from "ttcr_t.h" namespace "ttcr":
    cdef cppclass sxyz[T]:
        sxyz(T, T, T) except +

cdef extern from "Grid3Dttcr.h" namespace "ttcr":
    cdef cppclass Grid3Dttcr:
        Grid3Dttcr(string&, uint32_t, uint32_t, uint32_t, double, double, double, double, double, int, bool, size_t) except +
        string getType()
        void setSlowness(const vector[double]&) except +
        int raytrace(const vector[sxyz[double]]&,
                     const vector[double]&,
                     const vector[sxyz[double]]&,
                     double*)
        int raytrace(const vector[sxyz[double]]&,
                     const vector[double]&,
                     const vector[sxyz[double]]&,
                     double*, object, double*)
        int raytrace(const vector[sxyz[double]]& Tx,
                     const vector[double]& tTx,
                     const vector[sxyz[double]]& Rx,
                     double*, object, double*, object)
        int raytrace(const vector[sxyz[double]]& Tx,
                     const vector[double]& tTx,
                     const vector[sxyz[double]]& Rx,
                     double*, object, object)
        int raytrace(const vector[sxyz[double]]& Tx,
                     const vector[double]& tTx,
                     const vector[sxyz[double]]& Rx,
                     double*, object)


cdef class Grid3Dcpp:
    cdef uint32_t nx
    cdef uint32_t ny
    cdef uint32_t nz
    cdef Grid3Dttcr* grid
    def __cinit__(self, gridType, uint32_t nx, uint32_t ny, uint32_t nz, double dx,
                  double xmin, double ymin, double zmin,
                  double eps, int maxit, bool weno, size_t nthreads):
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.grid = new Grid3Dttcr(gridType, nx, ny, nz, dx, xmin, ymin, zmin, eps, maxit, weno, nthreads)

    def __dealloc__(self):
        del self.grid

    def getType(self):
        return self.grid.getType()

    def raytrace(self, slowness, Tx, Rx, t0, nout):

        # assing model data
        cdef vector[double] slown
        if self.grid.getType() == b'node':
            nx = self.nx+1
            ny = self.ny+1
            nz = self.nz+1
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        # slowness is in 'C' order and we must pass it in 'F' order
                        slown.push_back(slowness[(i*ny + j)*nz + k])
        elif self.grid.getType() == b'cell':
            for k in range(self.nz):
                for j in range(self.ny):
                    for i in range(self.nx):
                        # slowness is in 'C' order and we must pass it in 'F' order
                        slown.push_back(slowness[(i*self.ny + j)*self.nz + k])
        self.grid.setSlowness(slown)

        # create C++ input variables
        cdef vector[sxyz[double]] cTx
        for t in Tx:
            cTx.push_back(sxyz[double](t[0], t[1], t[2]))

        cdef vector[sxyz[double]] cRx
        for r in Rx:
            cRx.push_back(sxyz[double](r[0], r[1], r[2]))

        cdef vector[double] ct0
        for tmp in t0:
            ct0.push_back(tmp)

        # instantiate output variables
        cdef np.ndarray tt = np.empty([Rx.shape[0],], dtype=np.double)  # tt should be the right size
        cdef np.ndarray v0 = np.empty([Rx.shape[0],], dtype=np.double)

        if nout == 1:
            if self.grid.raytrace(cTx, ct0, cRx, <double*> np.PyArray_DATA(tt)) != 0:
                raise RuntimeError()
            return tt

        elif nout == 2:
            if self.grid.getType() == b'node':
                raise RuntimeError('raytracing with 2 output arguments not defined for grids with velocity defined at nodes')

            Ldata = ([0.0], [0.0], [0.0])

            if self.grid.raytrace(cTx, ct0, cRx, <double*> np.PyArray_DATA(tt), Ldata) != 0:
                raise RuntimeError()

            M = Rx.shape[0]
            N = len(slowness)
            L = csr_matrix(Ldata, shape=(M,N))

            return tt,L

        elif nout == 3:
            rays = tuple([ [0.0] for i in range(Rx.shape[0]) ])
            if self.grid.getType() == b'node':

                if self.grid.raytrace(cTx, ct0, cRx, <double*> np.PyArray_DATA(tt), rays, <double*> np.PyArray_DATA(v0)) != 0:
                    raise RuntimeError()

                return tt, rays, v0
            elif self.grid.getType() == b'cell':
                Ldata = ([0.0], [0.0], [0.0])

                if self.grid.raytrace(cTx, ct0, cRx, <double*> np.PyArray_DATA(tt), rays, Ldata) != 0:
                    raise RuntimeError()

                M = Rx.shape[0]
                N = len(slowness)
                L = csr_matrix(Ldata, shape=(M,N))

                return tt, L, rays

        elif nout == 4:
            if self.grid.getType() == b'cell':
                raise RuntimeError('raytracing with 4 output arguments not defined for grids with constant velocity cells')

            rays = tuple([ [0.0] for i in range(Rx.shape[0]) ])

            nTx = np.unique(Tx, axis=0).shape[0]

            M = tuple([ ([0.0],[0.0],[0.0]) for i in range(nTx) ])

            if self.grid.raytrace(cTx, ct0, cRx, <double*> np.PyArray_DATA(tt), rays, <double*> np.PyArray_DATA(v0), M) != 0:
                raise RuntimeError()

            return tt, rays, v0, M

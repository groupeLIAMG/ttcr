# -*- coding: utf-8 -*-

"""
    Copyright 2016 Bernard Giroux
    email: bernard.giroux@ete.inrs.ca
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it /will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport uint32_t

import numpy as np
cimport numpy as np

from scipy.sparse import csr_matrix

from utils import nargout

cdef extern from "ttcr_t.h" namespace "ttcr":
    cdef cppclass sxz[T]:
        sxz(T, T) except +


cdef extern from "Grid2Dttcr.h" namespace "ttcr":
    cdef cppclass Grid2Dttcr:
        Grid2Dttcr(string&, uint32_t, uint32_t, double, double, double, double, uint32_t, uint32_t, size_t) except +
        string getType()
        void setSlowness(const vector[double]&) except +
        void setXi(const vector[double]&) except +
        void setTheta(const vector[double]&) except +
        int raytrace(vector[sxz[double]]&,vector[double]&,vector[sxz[double]]&,double*,object,object)
        int raytrace(vector[sxz[double]]&,vector[double]&,vector[sxz[double]]&,double*,object)
        @staticmethod
        void Lsr2d(double*,double*,size_t,double*,size_t,double*,size_t,object)
        @staticmethod
        void Lsr2da(double*,double*,size_t,double*,size_t,double*,size_t,object)



cdef class Grid2Dcpp:
    cdef Grid2Dttcr* grid
    def __cinit__(self, gridType, uint32_t nx, uint32_t nz, double dx, double dz,double xmin, double zmin,uint32_t nsnx, uint32_t nsnz,size_t nthreads):
        self.grid = new Grid2Dttcr(gridType, nx, nz, dx, dz, xmin, zmin, nsnx, nsnz, nthreads)

    def __dealloc__(self):
        del self.grid

    def getType(self):
        return self.grid.getType()

    def raytrace(self, slowness, xi, theta, Tx, Rx, t0):

        nout = nargout()
        # check if types are consistent with input data
        if len(xi) != 0:
            if len(theta) != 0:
                if self.grid.getType() != b'tilted':
                    raise TypeError('Grid should handle raytracing in tilted elliptically anisotropic media')
            else:
                if self.grid.getType() != b'elliptical':
                    raise TypeError('Grid should handle raytracing in elliptically anisotropic media')
        else:
            if self.grid.getType() != b'iso':
                raise TypeError('Grid should handle raytracing in isotropic media')

        # assing model data
        cdef vector[double] slown# = slowness.tolist()
        for tmp in slowness:
            slown.push_back(tmp)
        self.grid.setSlowness(slown)
        cdef vector[double] x
        if len(xi) != 0:
            for tmp in xi:
                x.push_back(tmp)
            self.grid.setXi(x)
        cdef vector[double] t
        if len(theta) != 0:
            for tmp in theta:
                t.push_back(tmp)
            self.grid.setTheta(t)

        # create C++ input variables
        cdef vector[sxz[double]] cTx
        for t in Tx:
            cTx.push_back(sxz[double](t[0], t[2]))

        cdef vector[sxz[double]] cRx
        for r in Rx:
            cRx.push_back(sxz[double](r[0], r[2]))

        cdef vector[double] ct0# = t0.tolist()
        for tmp in t0:
            ct0.push_back(tmp)

        # instantiate output variables
        cdef np.ndarray tt = np.empty([Rx.shape[0],], dtype=np.double)  # tt should be the right size

        if nout==2:
            Ldata = ([0.0], [0.0], [0.0])

            if self.grid.raytrace(cTx, ct0, cRx, <double*> np.PyArray_DATA(tt), Ldata) != 0:
                raise RuntimeError()

            M = Rx.shape[0]
            N = len(slowness)
            if self.grid.getType() != b'iso':
                N = 2*N

            L = csr_matrix(Ldata, shape=(M,N))

            return tt,L
            
        elif nout==3:
            rays = tuple([ [0.0] for i in range(Rx.shape[0]) ])
            Ldata = ([0.0], [0.0], [0.0])

            if self.grid.raytrace(cTx, ct0, cRx, <double*> np.PyArray_DATA(tt), rays, Ldata) != 0:
                raise RuntimeError()

            M = Rx.shape[0]
            N = len(slowness)
            if self.grid.getType() != b'iso':
                N = 2*N

            L = csr_matrix(Ldata, shape=(M,N))

            return tt,L,rays


    @staticmethod
    def Lsr2d(Tx, Rx, grx, grz):

        cdef size_t nTx = Tx.shape[0]
        cdef size_t n_grx = grx.shape[0]
        cdef size_t n_grz = grz.shape[0]

        Ldata = ([0.0], [0.0], [0.0])

        Grid2Dttcr.Lsr2d(<double*> np.PyArray_DATA(Tx), <double*> np.PyArray_DATA(Rx), nTx, <double*> np.PyArray_DATA(grx), n_grx, <double*> np.PyArray_DATA(grz), n_grz, Ldata)

        M = nTx
        N = (n_grx-1)*(n_grz-1)
        L = csr_matrix(Ldata, shape=(M,N))

        return L


    @staticmethod
    def Lsr2da(Tx, Rx, grx, grz):

        cdef size_t nTx = Tx.shape[0]
        cdef size_t n_grx = grx.shape[0]
        cdef size_t n_grz = grz.shape[0]

        Ldata = ([0.0], [0.0], [0.0])

        Grid2Dttcr.Lsr2da(<double*> np.PyArray_DATA(Tx), <double*> np.PyArray_DATA(Rx), nTx, <double*> np.PyArray_DATA(grx), n_grx, <double*> np.PyArray_DATA(grz), n_grz, Ldata)

        M = nTx
        N = 2*(n_grx-1)*(n_grz-1)
        L = csr_matrix(Ldata, shape=(M,N))

        return L

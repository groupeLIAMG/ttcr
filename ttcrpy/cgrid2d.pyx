# -*- coding: utf-8 -*-

"""
    2D grid for raytracing based on the shortest path method

    Important: the raytracing codes are based on a column-major order
                for the slowness vector (Z is the "fast" axis).
                To visualize the slowness model with Z axis vertical and X horizontal,
                the vector should be reshaped as
                slowness.reshape(nx,nz).T

    This code is part of ttcr ( https://github.com/groupeLIAMG/ttcr )
"""

"""
    Copyright 2016 Bernard Giroux
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

import numpy as np
cimport numpy as np

from scipy.sparse import csr_matrix

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
        void raytrace(vector[sxz[double]]&,vector[double]&,vector[sxz[double]]&,double*,object,object) except +
        void raytrace(vector[sxz[double]]&,vector[double]&,vector[sxz[double]]&,double*,object) except +
        void raytrace(vector[sxz[double]]&,vector[double]&,vector[sxz[double]]&,double*) except +
        @staticmethod
        int Lsr2d(double*,double*,size_t,double*,size_t,double*,size_t,object)
        @staticmethod
        int Lsr2da(double*,double*,size_t,double*,size_t,double*,size_t,object)


cdef extern from "Cell.h" namespace "ttcr":
    cdef cppclass Cell[T, NODE, S]:
        Cell(size_t)
    cdef cppclass CellElliptical[T, NODE, S]:
        CellElliptical(size_t)
    cdef cppclass CellTiltedElliptical[T, NODE, S]:
        CellTiltedElliptical(size_t)


cdef extern from "Node2Dcsp.h" namespace "ttcr":
    cdef cppclass Node2Dcsp[T1, T2]:
        Node2Dcsp(size_t)


cdef extern from "Grid2Drcsp.h" namespace "ttcr":
    cdef cppclass Grid2Drcsp[T1, T2, S, CELL]:
        Grid2Drcsp(T2, T2, T1, T1, T1, T1, T2, T2, size_t)
        size_t getNthreads()
        void setSlowness(vector[T1]&) except +


cdef class Grid2Dcpp:
    """
    Grid2Dcpp(type, nx, nz, dx, dz, xmin, zmin, nsnx, nsnz, nthreads)

    Parameters
    ----------
    type : type of media
             'iso'          isotorpic
             'elliptical'   elliptically anisotropic
             'tilted'       tilted elliptically anisotropic
    nx : number of cells in x
    nz : number of cells in z
    dx : cell size in x
    dz : cell size in z
    xmin : origin inx
    zmin : origin in z
    nsnx : number of secondary nodes in x
    nsnz : number of secondary nodes in z
    nthreads : number of threads
    """
    cdef Grid2Dttcr* grid
    def __cinit__(self, gridType, uint32_t nx, uint32_t nz, double dx, double dz,
                  double xmin, double zmin, uint32_t nsnx, uint32_t nsnz,
                  size_t nthreads):
        self.grid = new Grid2Dttcr(gridType, nx, nz, dx, dz, xmin, zmin, nsnx, nsnz, nthreads)

    def __dealloc__(self):
        del self.grid

    def getType(self):
        return self.grid.getType()

    def raytrace(self, slowness, xi, theta, Tx, Rx, t0, nout):
        """
        raytrace(slowness, xi, theta, Tx, Rx, t0, nout) -> tt,L,rays

        Parameters
        ----------
        slowness : vector of slowness at grid cells
        xi : vector of anisotropy ratio at grid cells
        theta : vector of anisotropy angle at grid cells
        Tx : source coordinates, nTx by 2
               1st column contains X coordinates,
               2nd contains Z coordinates
        Rx : receiver coordinates, nTx by 2
               1st column contains X coordinates,
               2nd contains Z coordinates
        t0 : origin time
        nout : number of output variables (1, 2 or 3)

        Returns
        -------
        tt : travel times
        L : (if nout > 1) data kernel matrix (tt = L*slowness)
        rays : (if nout == 3) coordinates of rays
        """

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
            cTx.push_back(sxz[double](t[0], t[1]))

        cdef vector[sxz[double]] cRx
        for r in Rx:
            cRx.push_back(sxz[double](r[0], r[1]))

        cdef vector[double] ct0# = t0.tolist()
        for tmp in t0:
            ct0.push_back(tmp)

        # instantiate output variables
        cdef np.ndarray tt = np.empty([Rx.shape[0],], dtype=np.double)  # tt should be the right size

        if nout==1:
            self.grid.raytrace(cTx, ct0, cRx, <double*> np.PyArray_DATA(tt))
            return tt

        elif nout==2:
            Ldata = ([0.0], [0.0], [0.0])
            self.grid.raytrace(cTx, ct0, cRx, <double*> np.PyArray_DATA(tt), Ldata)

            M = Rx.shape[0]
            N = len(slowness)
            if self.grid.getType() != b'iso':
                N = 2*N

            L = csr_matrix(Ldata, shape=(M,N))

            return tt,L

        elif nout==3:
            rays = tuple([ [0.0] for i in range(Rx.shape[0]) ])
            Ldata = ([0.0], [0.0], [0.0])

            self.grid.raytrace(cTx, ct0, cRx, <double*> np.PyArray_DATA(tt), rays, Ldata)

            M = Rx.shape[0]
            N = len(slowness)
            if self.grid.getType() != b'iso':
                N = 2*N

            L = csr_matrix(Ldata, shape=(M,N))

            return tt,L,rays


    @staticmethod
    def Lsr2d(Tx, Rx, grx, grz):
        """
        Lsr2d(Tx, Rx, grx, grz) -> L

        Raytracing with straight rays in 2D

        Parameters
        ----------
        Tx : source coordinates, nTx by 2
              1st column contains X coordinates,
              2nd contains Z coordinates
        Rx : receiver coordinates, nTx by 2
              1st column contains X coordinates,
              2nd contains Z coordinates
        grx : grid node coordinates along x
        grz : grid node coordinates along z

        Returns
        -------
        L : data kernel matrix (tt = L*slowness)
        """

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
        """
        Lsr2da(Tx, Rx, grx, grz) -> L

        Raytracing with straight rays in 2D elliptically anisotropic media

        Parameters
        ----------
        Tx : source coordinates, nTx by 2
               1st column contains X coordinates,
               2nd contains Z coordinates
        Rx : receiver coordinates, nTx by 2
               1st column contains X coordinates,
               2nd contains Z coordinates
        grx : grid node coordinates along x
        grz : grid node coordinates along z

        Returns
        -------
        L : data kernel matrix (tt = L*slowness)
        """

        cdef size_t nTx = Tx.shape[0]
        cdef size_t n_grx = grx.shape[0]
        cdef size_t n_grz = grz.shape[0]

        Ldata = ([0.0], [0.0], [0.0])

        Grid2Dttcr.Lsr2da(<double*> np.PyArray_DATA(Tx), <double*> np.PyArray_DATA(Rx), nTx, <double*> np.PyArray_DATA(grx), n_grx, <double*> np.PyArray_DATA(grz), n_grz, Ldata)

        M = nTx
        N = 2*(n_grx-1)*(n_grz-1)
        L = csr_matrix(Ldata, shape=(M,N))

        return L


cdef class Grid2Diso:
    cdef uint32_t nx
    cdef uint32_t nz
    cdef Grid2Drcsp[double, uint32_t, sxz[double], Cell[double, Node2Dcsp[double, uint32_t], sxz[double]]]* grid

    def __cinit__(self, uint32_t nx, uint32_t nz, double dx, double dz,
                  double xmin, double zmin, uint32_t nsnx, uint32_t nsnz,
                  size_t nthreads):
        self.grid = new Grid2Drcsp[double, uint32_t, sxz[double], Cell[double, Node2Dcsp[double, uint32_t], sxz[double]]](nx, nz, dx, dz, xmin, zmin, nsnx, nsnz, nthreads)

    def __dealloc__(self):
        del self.grid

    def get_type(self):
        return 'iso'

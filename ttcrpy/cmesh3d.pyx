# -*- coding: utf-8 -*-

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
from libcpp cimport bool

import numpy as np
cimport numpy as np

from scipy.sparse import csr_matrix

cdef extern from "ttcr_t.h" namespace "ttcr":
    cdef cppclass sxyz[T]:
        sxyz(T, T, T) except +

cdef extern from "ttcr_t.h" namespace "ttcr":
    cdef cppclass tetrahedronElem[T]:
        tetrahedronElem(T, T, T, T) except +


cdef extern from "Mesh3Dttcr.h" namespace "ttcr":
    cdef cppclass Mesh3Dttcr:
        Mesh3Dttcr(vector[sxyz[double]]&, vector[tetrahedronElem[uint32_t]]&,
                   double, int, bool rp, size_t) except +
        void setSlowness(const vector[double]&) except +
        void raytrace(const vector[sxyz[double]]&,
                     const vector[double]&,
                     const vector[sxyz[double]]&,
                     double*) except +
        void raytrace(const vector[sxyz[double]]&,
                     const vector[double]&,
                     const vector[sxyz[double]]&,
                     double*, object, double*) except +
#        void raytrace(const vector[sxyz[double]]& Tx,
#                     const vector[double]& tTx,
#                     const vector[sxyz[double]]& Rx,
#                     double*, object, double*, object) except +

cdef class Mesh3Dcpp:
    cdef Mesh3Dttcr* mesh
    def __cinit__(self, nodes, tetra, const double eps, const int maxit, const bool rp, const size_t nt):
        cdef vector[sxyz[double]] nod
        for no in nodes:
            nod.push_back(sxyz[double](no[0], no[1], no[2]))
        cdef vector[tetrahedronElem[uint32_t]] tet
        for t in tetra:
            tet.push_back(tetrahedronElem[uint32_t](t[0], t[1], t[2], t[3]))
        self.mesh = new Mesh3Dttcr(nod, tet, eps, maxit, rp, nt)

    def __dealloc__(self):
        del self.mesh

    def raytrace(self, slowness, Tx, Rx, t0):
        nout = 4

        # assing model data
        cdef vector[double] slown
        for tmp in slowness:
            slown.push_back(tmp)
        self.mesh.setSlowness(slown)

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
            self.mesh.raytrace(cTx, ct0, cRx, <double*> np.PyArray_DATA(tt))
            return tt

        elif nout == 3:
            rays = tuple([ [0.0] for i in range(Rx.shape[0]) ])

            self.mesh.raytrace(cTx, ct0, cRx, <double*> np.PyArray_DATA(tt), rays, <double*> np.PyArray_DATA(v0))
            return tt, rays, v0

#        elif nout == 4:
#            rays = tuple([ [0.0] for i in range(Rx.shape[0]) ])
#
#            nTx = np.unique(Tx, axis=0).shape[0]
#
#            M = tuple([ ([0.0],[0.0],[0.0]) for i in range(nTx) ])
#
#            self.mesh.raytrace(cTx, ct0, cRx, <double*> np.PyArray_DATA(tt), rays, <double*> np.PyArray_DATA(v0), M)
#            return tt, rays, v0, M

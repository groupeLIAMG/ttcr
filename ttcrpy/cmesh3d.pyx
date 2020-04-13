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
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.string cimport memcpy
import numpy as np
cimport numpy as np
import scipy.sparse as sp


#from scipy.sparse import csr_matrix


cdef extern from "ttcr_t.h" namespace "ttcr":
    cdef cppclass sxyz[T]:
        sxyz(T, T, T) except +

cdef extern from "ttcr_t.h" namespace "ttcr":
    cdef cppclass tetrahedronElem[T]:
        tetrahedronElem(T, T, T, T) except +


cdef extern from "Mesh3Dttcr.h" namespace "ttcr":
    cdef cppclass Mesh3Dttcr:
        Mesh3Dttcr(vector[sxyz[double]]&, vector[tetrahedronElem[uint32_t]]&,
                   const int, const size_t, const int ,  const size_t ,
                   const double) except +
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
                     const vector[uint32_t]& ,
                     const bool &,const bool &,double*, object, double*, object)
        object ComputeD(const vector[sxyz[double]]& Pts)
        int ComputeK(const int &, const string &,const int &,const size_t &,const bool &, object)
        bool CheckPoint(const vector[sxyz[double]] & Point)

cdef class Mesh3Dcpp:
    cdef Mesh3Dttcr* mesh
    def __cinit__(self, nodes, tetra, const int ns,
                  const size_t nt, const  int verbose, const size_t nst,
                  const double R_ratio):
        cdef vector[sxyz[double]] nod
        for no in nodes:
            nod.push_back(sxyz[double](no[0], no[1], no[2]))
        cdef vector[tetrahedronElem[uint32_t]] tet
        for t in tetra:
            tet.push_back(tetrahedronElem[uint32_t](t[0], t[1], t[2], t[3]))
        self.mesh = new Mesh3Dttcr(nod, tet, ns, nt,verbose,nst,R_ratio)

    def __dealloc__(self):
        del self.mesh
    def ComputeD(self,Points):
        cdef vector[sxyz[double]] Pnts
        for P in Points:
            Pnts.push_back(sxyz[double](P[0],P[1],P[2]))
        D=self.mesh.ComputeD(Pnts)
        return D
    def BuildK(self, order=1, method='4D' , expansion=2, weighting=False, minpnts=10):
        if order !=1 and order !=2:
            raise RuntimeError('invalid order')
        if expansion !=1 and expansion !=2 and expansion!=12:
            raise RuntimeError('invalid expansion order')
        cdef char* methodc
        if method=='3d' or method=='3D':
            methodc='3D'
        elif method=='4d' or method=='4D':
            methodc='4D'
        elif method=='ABM':
            methodc='ABM'
            print('Warning: the ABM supports only a simple expansion')
            expansion=1
        else:
            raise RuntimeError('invalid method')
        if order==2 and expansion!=2:
            if expansion ==12:
                expansion=2
                order=1
            elif expansion ==1:
                expansion=1
                order=1
            method_object = <bytes> methodc
            K = tuple([ ([0.0],[0.0],[0.0]) for i in range(3) ])
            if self.mesh.ComputeK(order,method_object,expansion,minpnts,weighting,K)!=0:
                raise RuntimeError()
            kx=sp.csr_matrix((K[0][0],K[0][1], K[0][2]))
            ky=sp.csr_matrix((K[1][0],K[1][1], K[1][2]))
            kz=sp.csr_matrix((K[2][0],K[2][1], K[2][2]))
            kx2=kx.dot(kx)
            ky2=ky.dot(ky)
            kz2=kz.dot(kz)
            return kx2,ky2,kz2
        method_object = <bytes> methodc
        K = tuple([ ([0.0],[0.0],[0.0]) for i in range(3) ])
        if self.mesh.ComputeK(order,method_object,expansion,minpnts,weighting,K)!=0:
            raise RuntimeError()
        kx=sp.csr_matrix((K[0][0],K[0][1], K[0][2]))
        ky=sp.csr_matrix((K[1][0],K[1][1], K[1][2]))
        kz=sp.csr_matrix((K[2][0],K[2][1], K[2][2]))
        return kx,ky,kz
    def CheckPoint(self,Point):
        cdef vector [sxyz[double]] Pnts
        Pnts.push_back(sxyz[double](Point[0],Point[1],Point[2]))
        return self.mesh.CheckPoint(Pnts)
    def raytrace1(self, slowness, Tx, Rx, t0,Rxindices,st_corr=True,slow=False):

        # assing model data
        cdef vector[double] slown
        for tmp in slowness:
            slown.push_back(tmp)
        self.mesh.setSlowness(slown)

        # create C++ input variables
        cdef vector[sxyz[double]] cTx
        for t in Tx:
            cTx.push_back(sxyz[double](t[0], t[1], t[2]))

        cdef vector[uint32_t] Rxidx
        for ir in Rxindices:
            Rxidx.push_back(ir)
        cdef vector[sxyz[double]] cRx
        for r in Rx:
            cRx.push_back(sxyz[double](r[0], r[1], r[2]))

        cdef vector[double] ct0
        for tmp in t0:
            ct0.push_back(tmp)
        # instantiate output variables
        cdef np.ndarray tt = np.empty([Rx.shape[0],], dtype=np.double)  # tt should be the right size
        cdef np.ndarray v0 = np.empty([Rx.shape[0],], dtype=np.double)
        rays = tuple([ [0.0] for i in range(Rx.shape[0]) ])

        nTx = np.unique(Tx.view([('',Tx.dtype)]*Tx.shape[1])).view(Tx.dtype).reshape(-1,Tx.shape[1]).shape[0]   # get number of unique Tx

        M = tuple([ ([0.0],[0.0],[0.0]) for i in range(nTx) ])

#        if self.mesh.raytrace(cTx, ct0, cRx, <double*> np.PyArray_DATA(tt), rays, <double*> np.PyArray_DATA(v0), M) != 0:
#            raise RuntimeError()
        if self.mesh.raytrace(cTx, ct0, cRx,Rxidx, st_corr,slow,<double*> np.PyArray_DATA(tt), rays, <double*> np.PyArray_DATA(v0), M) != 0:
            raise RuntimeError()

        return tt, rays, v0,M

    def raytrace2(self, Tx, Rx, t0):

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
        rays = tuple([ [0.0] for i in range(Rx.shape[0]) ])

        if self.mesh.raytrace(cTx, ct0, cRx, <double*> np.PyArray_DATA(tt), rays, <double*> np.PyArray_DATA(v0)) != 0:
            raise RuntimeError()

        return (tt,rays,v0)
    def SetSlowness(self,slowness):
        # assing model data
        cdef vector[double] slown
        for tmp in slowness:
            slown.push_back(tmp)
        self.mesh.setSlowness(slown)
    def raytrace3(self,slowness, Tx, Rx, t0):

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
        rays = tuple([ [0.0] for i in range(Rx.shape[0]) ])

        if self.mesh.raytrace(cTx, ct0, cRx, <double*> np.PyArray_DATA(tt), rays, <double*> np.PyArray_DATA(v0)) != 0:
            raise RuntimeError()

        return (tt,rays,v0)

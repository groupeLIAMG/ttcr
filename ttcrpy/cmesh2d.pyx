# -*- coding: utf-8 -*-

"""
Code for raytracing on 2D triangular meshes based on the shortest path method


This code is part of ttcr ( https://github.com/groupeLIAMG/ttcr )
"""

"""
Copyright 2018 Bernard Giroux
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

import vtk
from vtk.util import numpy_support

cdef extern from "ttcr_t.h" namespace "ttcr":
    cdef cppclass sxz[T]:
        sxz(T, T) except +
        T x
        T z

    cdef cppclass siv[T]:
        siv(size_t, T) except +
        size_t i
        T v

    cdef cppclass triangleElem[T]:
        triangleElem(T, T, T) except +
        T i[3]
        T physical_entity


cdef extern from "Node2Dcsp.h" namespace "ttcr":
    cdef cppclass Node2Dcsp[T1,T2]:
        Node2Dcsp(size_t) except +


cdef extern from "Grid2Ducsp.h" namespace "ttcr":
    cdef cppclass Grid2Ducsp[T1,T2,NODE,S]:
        Grid2Ducsp(vector[S], vector[triangleElem[T2]], T2, size_t) except +
        size_t getNthreads()
        size_t getNumberOfCells()
        size_t getNumberOfNodes(bool)
        void setSlowness(T1 *, size_t) except +
        void raytrace(vector[S]&, vector[T1]&, vector[S]&, vector[T1]&, size_t) except +
        void raytrace(vector[S]&, vector[T1]&, vector[S]&, vector[T1]&, vector[vector[S]]&, vector[vector[siv[T1]]]&, size_t) except +
        void calculateArea(vector[T1] &) except +
        void interpolateAtNodes(vector[T1] &) except +
        void getTT(vector[T1]&, size_t)
        void getNodes(vector[S]&)
        void getTriangles(vector[vector[T2]]&)


cdef class Mesh2D:
    """
        Wrapper class for the C++ code

        Constructor

            Mesh2D(nodes, triangles, nsecondary, nthreads)

            Parameters
            ----------
            nodes : coordinates of nodes (ndarray of size nnodes x 2, double)
                   1st column contains X coordinates,
                   2nd contains Z coordinates
            triangles : indices of nodes forming triangles (ndarray of size ntriange x 3, int)
            nsecondary : number of secondary nodes in triangle edges (int)
            nthreads : number of threads/processes for raytracing (int)
    """
    cdef Grid2Ducsp[double, uint32_t, Node2Dcsp[double, uint32_t], sxz[double]] *mesh
    def __cinit__(self, nodes, triangles, uint32_t nsecondary, size_t nthreads):
        cdef vector[sxz[double]] no
        for n in nodes:
            no.push_back(sxz[double](n[0], n[1]))
        cdef vector[triangleElem[uint32_t]] tri
        for t in triangles:
            tri.push_back(triangleElem[uint32_t](t[0], t[1], t[2]))
        self.mesh = new Grid2Ducsp[double, uint32_t, Node2Dcsp[double, uint32_t], sxz[double]](no, tri, nsecondary, nthreads)

    def __dealloc__(self):
        del self.mesh

    def get_grid_traveltimes(self, thread_no=0):
        """
        Obtain traveltimes computed at primary grid nodes

        Parameters
        ----------
        thread_no : int
            thread used to computed traveltimes (default is 0)

        Returns
        -------
        tt: np ndarray
        """
        if thread_no >= self.mesh.getNthreads():
            raise ValueError('Thread number is larger than number of threads')
        cdef vector[double] tmp
        cdef int n
        self.mesh.getTT(tmp, thread_no)
        tt = np.empty((tmp.size(),))
        for n in range(tmp.size()):
            tt[n] = tmp[n]
        return tt

    def set_slowness(self, slowness):
        """
        Assign slowness at cells

        set_slowness(slowness)

        Parameters
        ----------
        slowness : ndarray of size ntriangle.  Indices of slowness values
                   must correspond to triangle indices given in constructor
        """
        if not slowness.flags['C_CONTIGUOUS']:
            slowness = np.ascontiguousarray(slowness)

        # ::1 tells cython that data is contiguous
        cdef double[::1] slowness_memview = slowness

        self.mesh.setSlowness(&slowness_memview[0], slowness_memview.shape[0])

    def raytrace(self, slowness, Tx, Rx, t0=0.0, nout=1, thread_no=0):
        """
        raytrace(slowness, Tx, Rx, t0, nout, thread_no) -> tt,L,rays

        Compute traveltimes for a single source.  The source can comprise many
        points, each having its own otigin time (t0).

        Parameters
        ----------
        slowness : vector of slowness at cells
        Tx : source coordinates, ndarray of size nTx by 2
            1st column contains X coordinates,
            2nd contains Z coordinates
        Rx : receiver coordinates, ndarray of size nRx by 2
            1st column contains X coordinates,
            2nd contains Z coordinates
        t0 : origin times, scalar or ndarray of size nTx
        nout : number of output variables (1, 2 or 3)
        thread_no : thread/process number on which computation should be run

        Returns
        -------
        tt : travel times
        L : (if nout > 1) data kernel matrix (scipy csr) (tt = L*slowness)
        rays : (if nout == 3) coordinates of rays
        """
        cdef double[::1] slowness_memview
        if slowness is not None:
            if not slowness.flags['C_CONTIGUOUS']:
                slowness = np.ascontiguousarray(slowness)

            slowness_memview = slowness
            self.mesh.setSlowness(&slowness_memview[0], slowness_memview.shape[0])

        cdef vector[sxz[double]] vTx
        cdef vector[sxz[double]] vRx
        cdef vector[double] vt0
        cdef vector[double] vtt

        for t in Tx:
            vTx.push_back(sxz[double](t[0], t[1]))
        for r in Rx:
            vRx.push_back(sxz[double](r[0], r[1]))

        if np.isscalar(t0):
            for n in range(Tx.shape[0]):
                vt0.push_back(t0)
        else:
            for t in t0:
                vt0.push_back(t)

        vtt.resize(Rx.shape[0])

        cdef vector[vector[sxz[double]]] r_data
        cdef vector[vector[siv[double]]] l_data
        cdef double v0 = 0.0

        if nout == 1:
            self.mesh.raytrace(vTx, vt0, vRx, vtt, thread_no)

            tt = np.empty((Rx.shape[0],))
            for n in range(Rx.shape[0]):
                tt[n] = vtt[n]

            return tt

        elif nout==2:
            self.mesh.raytrace(vTx, vt0, vRx, vtt, r_data, l_data, thread_no)

            tt = np.empty((Rx.shape[0],))
            for n in range(Rx.shape[0]):
                tt[n] = vtt[n]

            indptr = np.empty((Rx.shape[0]+1,), dtype=np.int64)
            nel = 0
            for n in range(l_data.size()):
                nel += l_data[n].size()
            indices = np.empty((nel,), dtype=np.int64)
            val = np.empty((nel,))

            k = 0
            M = Rx.shape[0]
            N = self.mesh.getNumberOfCells()
            for i in range(M):
                indptr[i] = k
                for j in range(N):
                    for n in range(l_data[i].size()):
                        if l_data[i][n].i == j:
                            indices[k] = j
                            val[k] = l_data[i][n].v
                            k += 1

            indptr[M] = k
            L = csr_matrix((val, indices, indptr), shape=(M,N))

            return tt, L

        elif nout == 3:
            self.mesh.raytrace(vTx, vt0, vRx, vtt, r_data, l_data, thread_no)

            rays = [ [0.0] for i in range(Rx.shape[0]) ]
            tt = np.empty((Rx.shape[0],))
            for n in range(Rx.shape[0]):
                tt[n] = vtt[n]
                rays[n] = np.empty((r_data[n].size(), 2))
                for nn in range(r_data[n].size()):
                    rays[n][nn, 0] = r_data[n][nn].x
                    rays[n][nn, 1] = r_data[n][nn].z

            indptr = np.empty((Rx.shape[0]+1,), dtype=np.int64)
            nel = 0
            for n in range(l_data.size()):
                nel += l_data[n].size()
            indices = np.empty((nel,), dtype=np.int64)
            val = np.empty((nel,))

            k = 0
            M = Rx.shape[0]
            N = self.mesh.getNumberOfCells()
            for i in range(M):
                indptr[i] = k
                for j in range(N):
                    for n in range(l_data[i].size()):
                        if l_data[i][n].i == j:
                            indices[k] = j
                            val[k] = l_data[i][n].v
                            k += 1

            indptr[M] = k
            L = csr_matrix((val, indices, indptr), shape=(M,N))

            return tt, L, rays

    def to_vtk(self, fields, filename):
        """
        Save grid variables to VTK format

        Parameters
        ----------
        fields: dict
            dict of variables to save to file
            variables should be numpy ndarrays of size equal to either the
            number of nodes of the number of cells of the grid
        filename: str
            Name of file without extension for saving (extension vtr will be
            added)

        Notes
        -----
        VTK files can be visualized with Paraview (https://www.paraview.org)
        """
        cdef vector[sxz[double]] nodes
        cdef vector[vector[uint32_t]] triangles
        self.mesh.getNodes(nodes)
        self.mesh.getTriangles(triangles)

        cdef int n
        ugrid = vtk.vtkUnstructuredGrid()
        pts = vtk.vtkPoints()
        tri = vtk.vtkTriangle()
        for n in range(nodes.size()):
            pts.InsertNextPoint(nodes[n].x, 0.0, nodes[n].z)
        ugrid.SetPoints(pts)
        for n in range(triangles.size()):
            tri.GetPointIds().SetId(0, triangles[n][0])
            tri.GetPointIds().SetId(1, triangles[n][1])
            tri.GetPointIds().SetId(2, triangles[n][2])
            ugrid.InsertNextCell(tri.GetCellType(), tri.GetPointIds())

        for fn in fields:
            data = fields[fn]

            scalar = vtk.vtkDoubleArray()
            scalar.SetName(fn)
            scalar.SetNumberOfComponents(1)
            scalar.SetNumberOfTuples(data.size)
            if data.size == self.mesh.getNumberOfNodes(1):
                for n in range(nodes.size()):
                    scalar.SetTuple1(n, data[n])
                ugrid.GetPointData().AddArray(scalar)
            elif data.size == self.mesh.getNumberOfCells():
                for n in range(triangles.size()):
                    scalar.SetTuple1(n, data[n])
                ugrid.GetCellData().AddArray(scalar)
            else:
                raise ValueError('Field {0:s} has incorrect size'.format(fn))

            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetFileName(filename+'.vtu')
            writer.SetInputData(ugrid)
            writer.SetDataModeToBinary()
            writer.Update()

#    def get_Dx_Dz(self):
#        """
#            Compute spatial derivatives on mesh
#        """
#        cdef vector[double] area
#        self.mesh.calculateArea(area)
#
#        # TODO: complete

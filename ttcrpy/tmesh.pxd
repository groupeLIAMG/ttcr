
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport uint32_t, int64_t
from libc.math cimport sqrt
from libcpp cimport bool

from ttcrpy.common cimport sxz, sxyz, siv, siv2, sijv, Node3Dc, Node3Dcsp, \
Node3Dn, Node3Dnsp, Cell, Node2Dc, Node2Dcsp, Node2Dn, Node2Dnsp


cdef extern from "ttcr_t.h" namespace "ttcr" nogil:
    cdef cppclass triangleElem[T]:
        triangleElem(T, T, T) except +
        T i[3]
        T physical_entity

cdef extern from "ttcr_t.h" namespace "ttcr" nogil:
    cdef cppclass tetrahedronElem[T]:
        tetrahedronElem(T, T, T, T) except +
        T i[4]
        T physical_entity


cdef extern from "Grid3D.h" namespace "ttcr" nogil:
    cdef cppclass Grid3D[T1,T2]:
        void setUsePool(bool up)
        void setTraveltimeFromRaypath(bool)
        void checkPts(vector[sxyz[T1]]&) except +
        void computeD(vector[sxyz[T1]] &, vector[vector[sijv[T1]]] &)
        void computeK(vector[vector[vector[siv[T1]]]]&, int, int, bool, bool, int) except +
        size_t getNthreads()
        void setSlowness(vector[T1]&) except +
        T1 computeSlowness(sxyz[T1]&) except +
        void getTT(vector[T1]& tt, size_t threadNo) except +
        void raytrace(vector[sxyz[T1]]& Tx,
                      vector[T1]& t0,
                      vector[sxyz[T1]]& Rx,
                      vector[T1]& tt,
                      size_t thread_no) except +
        void raytrace(vector[sxyz[T1]]& Tx,
                      vector[T1]& t0,
                      vector[sxyz[T1]]& Rx,
                      vector[T1]& tt,
                      vector[vector[sxyz[T1]]]& r_data,
                      size_t thread_no) except +
        void raytrace(vector[sxyz[T1]]& Tx,
                      vector[T1]& t0,
                      vector[sxyz[T1]]& Rx,
                      vector[T1]& tt,
                      vector[vector[siv[T1]]]& l_data,
                      size_t thread_no) except +
        void raytrace(vector[sxyz[T1]]& Tx,
                      vector[T1]& t0,
                      vector[sxyz[T1]]& Rx,
                      vector[T1]& tt,
                      vector[vector[sxyz[T1]]]& r_data,
                      vector[vector[siv[T1]]]& l_data,
                      size_t thread_no) except +
        void raytrace(vector[sxyz[T1]]& Tx,
                      vector[T1]& t0,
                      vector[sxyz[T1]]& Rx,
                      vector[T1]& tt,
                      vector[vector[sijv[T1]]]& m_data,
                      size_t thread_no) except +
        void raytrace(vector[sxyz[T1]]& Tx,
                      vector[T1]& t0,
                      vector[sxyz[T1]]& Rx,
                      vector[T1]& tt,
                      vector[vector[sxyz[T1]]]& r_data,
                      vector[vector[sijv[T1]]]& m_data,
                      size_t thread_no) except +
        void raytrace(vector[vector[sxyz[T1]]]& Tx,
                      vector[vector[T1]]& t0,
                      vector[vector[sxyz[T1]]]& Rx,
                      vector[vector[T1]]& traveltimes) except +
        void raytrace(vector[vector[sxyz[T1]]]& Tx,
                      vector[vector[T1]]& t0,
                      vector[vector[sxyz[T1]]]& Rx,
                      vector[vector[T1]]& traveltimes,
                      vector[vector[vector[sxyz[T1]]]]& r_data) except +
        void raytrace(vector[vector[sxyz[T1]]]& Tx,
                      vector[vector[T1]]& t0,
                      vector[vector[sxyz[T1]]]& Rx,
                      vector[vector[T1]]& traveltimes,
                      vector[vector[vector[sijv[T1]]]]& m_data) except +
        void raytrace(vector[vector[sxyz[T1]]]& Tx,
                      vector[vector[T1]]& t0,
                      vector[vector[sxyz[T1]]]& Rx,
                      vector[vector[T1]]& traveltimes,
                      vector[vector[vector[sxyz[T1]]]]& r_data,
                      vector[vector[vector[sijv[T1]]]]& m_data) except +
        void raytrace(vector[vector[sxyz[T1]]]& Tx,
                      vector[vector[T1]]& t0,
                      vector[vector[sxyz[T1]]]& Rx,
                      vector[vector[T1]]& traveltimes,
                      vector[vector[vector[siv[T1]]]]& l_data) except +
        void raytrace(vector[vector[sxyz[T1]]]& Tx,
                      vector[vector[T1]]& t0,
                      vector[vector[sxyz[T1]]]& Rx,
                      vector[vector[T1]]& traveltimes,
                      vector[vector[vector[sxyz[T1]]]]& r_data,
                      vector[vector[vector[siv[T1]]]]& l_data) except +
        void getStraightRays(vector[sxyz[T1]]& Tx,
                             vector[sxyz[T1]]& Rx,
                             vector[vector[siv[T1]]]& l_data) except +

cdef extern from "Grid3Duc.h" namespace "ttcr" nogil:
    cdef cppclass Grid3Duc[T1,T2,N](Grid3D[T1,T2]):
        pass

cdef extern from "Grid3Dun.h" namespace "ttcr" nogil:
    cdef cppclass Grid3Dun[T1,T2,N](Grid3D[T1,T2]):
        pass

cdef extern from "Grid3Ducfs.h" namespace "ttcr" nogil:
    cdef cppclass Grid3Ducfs[T1, T2](Grid3Duc[T1,T2,Node3Dc[T1,T2]]):
        Grid3Ducfs(vector[sxyz[T1]], vector[tetrahedronElem[T2]], T1, int,
                   vector[sxyz[T1]]& ref_pts, int order, bool, bool, T1,
                   size_t, bool) except +

cdef extern from "Grid3Ducsp.h" namespace "ttcr" nogil:
    cdef cppclass Grid3Ducsp[T1, T2](Grid3Duc[T1,T2,Node3Dcsp[T1,T2]]):
        Grid3Ducsp(vector[sxyz[T1]], vector[tetrahedronElem[T2]],
                   int, bool, T1, size_t, bool) except +

cdef extern from "Grid3Ducdsp.h" namespace "ttcr" nogil:
    cdef cppclass Grid3Ducdsp[T1, T2](Grid3Duc[T1,T2,Node3Dc[T1,T2]]):
        Grid3Ducdsp(vector[sxyz[T1]], vector[tetrahedronElem[T2]],
                    int, int, T1, int, bool, T1, T1, bool, size_t, bool) except +

cdef extern from "Grid3Dunfs.h" namespace "ttcr" nogil:
    cdef cppclass Grid3Dunfs[T1, T2](Grid3Dun[T1,T2,Node3Dn[T1,T2]]):
        Grid3Dunfs(vector[sxyz[T1]], vector[tetrahedronElem[T2]], T1, int,
                   vector[sxyz[T1]]& ref_pts, int order, int, bool, bool, T1,
                   size_t, bool) except +

cdef extern from "Grid3Dunsp.h" namespace "ttcr" nogil:
    cdef cppclass Grid3Dunsp[T1, T2](Grid3Dun[T1,T2,Node3Dnsp[T1,T2]]):
        Grid3Dunsp(vector[sxyz[T1]], vector[tetrahedronElem[T2]], int, bool,
                   bool, T1, size_t, bool) except +

cdef extern from "Grid3Dundsp.h" namespace "ttcr" nogil:
    cdef cppclass Grid3Dundsp[T1, T2](Grid3Dun[T1,T2,Node3Dn[T1,T2]]):
        Grid3Dundsp(vector[sxyz[T1]], vector[tetrahedronElem[T2]], int, int, T1,
                    bool, int, bool, T1, T1, bool, size_t, bool) except +


cdef extern from "Grid2D.h" namespace "ttcr" nogil:
    cdef cppclass Grid2D[T1,T2,S]:
        void setUsePool(bool up)
        void setTraveltimeFromRaypath(bool)
        size_t getNthreads()
        void setSlowness(vector[T1]&) except +
        void getSlowness(vector[T1]&) except +
        void getTT(vector[T1]& tt, size_t threadNo) except +
        void raytrace(vector[S]& Tx,
                      vector[T1]& t0,
                      vector[S]& Rx,
                      vector[T1]& traveltimes,
                      size_t threadNo) except +
        void raytrace(vector[S]& Tx,
                      vector[T1]& t0,
                      vector[S]& Rx,
                      vector[T1]& traveltimes,
                      vector[vector[sxz[T1]]]& r_data,
                      size_t threadNo) except +
        void raytrace(vector[S]& Tx,
                      vector[T1]& t0,
                      vector[S]& Rx,
                      vector[T1]& traveltimes,
                      vector[vector[siv2[T1]]]& l_data,
                      size_t threadNo) except +
        void raytrace(vector[S]& Tx,
                      vector[T1]& t0,
                      vector[S]& Rx,
                      vector[T1]& traveltimes,
                      vector[vector[S]]& r_data,
                      vector[vector[siv2[T1]]]& l_data,
                      size_t threadNo) except +
        void raytrace(vector[vector[S]]& Tx,
                      vector[vector[T1]]& t0,
                      vector[vector[S]]& Rx,
                      vector[vector[T1]]& traveltimes) except +
        void raytrace(vector[vector[S]]& Tx,
                      vector[vector[T1]]& t0,
                      vector[vector[S]]& Rx,
                      vector[vector[T1]]& traveltimes,
                      vector[vector[vector[S]]]& r_data) except +
        void raytrace(vector[vector[S]]& Tx,
                      vector[vector[T1]]& t0,
                      vector[vector[S]]& Rx,
                      vector[vector[T1]]& traveltimes,
                      vector[vector[vector[siv2[T1]]]]& l_data) except +
        void raytrace(vector[vector[S]]& Tx,
                      vector[vector[T1]]& t0,
                      vector[vector[S]]& Rx,
                      vector[vector[T1]]& traveltimes,
                      vector[vector[vector[S]]]& r_data,
                      vector[vector[vector[siv2[T1]]]]& l_data) except +

cdef extern from "Grid2Duc.h" namespace "ttcr" nogil:
    cdef cppclass Grid2Duc[T1,T2,NODE,S](Grid2D[T1,T2,S]):
        pass

cdef extern from "Grid2Dun.h" namespace "ttcr" nogil:
    cdef cppclass Grid2Dun[T1,T2,NODE,S](Grid2D[T1,T2,S]):
        pass

cdef extern from "Grid2Ducsp.h" namespace "ttcr" nogil:
    cdef cppclass Grid2Ducsp[T1,T2,NODE,S](Grid2Duc[T1,T2,NODE,S]):
        Grid2Ducsp(vector[S]&, vector[triangleElem[T2]]&, T2, bool, size_t) except +

cdef extern from "Grid2Ducdsp.h" namespace "ttcr" nogil:
    cdef cppclass Grid2Ducdsp[T1,T2,S](Grid2Duc[T1,T2,Node2Dc[T1,T2],S]):
        Grid2Ducdsp(vector[S]&, vector[triangleElem[T2]]&, T2, int, T1, bool, bool, size_t) except +

cdef extern from "Grid2Ducfs.h" namespace "ttcr" nogil:
    cdef cppclass Grid2Ducfs[T1,T2,NODE,S](Grid2Duc[T1,T2,NODE,S]):
        Grid2Ducfs(vector[S]&, vector[triangleElem[T2]]&, T1, int,
                   vector[S]& ref_pts, int order, bool, size_t, bool) except +

cdef extern from "Grid2Dunsp.h" namespace "ttcr" nogil:
    cdef cppclass Grid2Dunsp[T1,T2,NODE,S](Grid2Dun[T1,T2,NODE,S]):
        Grid2Dunsp(vector[S]&, vector[triangleElem[T2]]&, T2, bool, size_t) except +

cdef extern from "Grid2Dundsp.h" namespace "ttcr" nogil:
    cdef cppclass Grid2Dundsp[T1,T2,S](Grid2Dun[T1,T2,Node2Dn[T1,T2],S]):
        Grid2Dundsp(vector[S]&, vector[triangleElem[T2]]&, int, int, T1, bool, bool, size_t) except +

cdef extern from "Grid2Dunfs.h" namespace "ttcr" nogil:
    cdef cppclass Grid2Dunfs[T1,T2,NODE,S](Grid2Dun[T1,T2,NODE,S]):
        Grid2Dunfs(vector[S]&, vector[triangleElem[T2]]&, T1, int,
                   vector[S]& ref_pts, int order, bool, size_t, bool) except +

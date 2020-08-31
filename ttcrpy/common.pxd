

cdef extern from "ttcr_t.h" namespace "ttcr" nogil:
    cdef cppclass sxyz[T]:
        sxyz(T, T, T) except +
        T x
        T y
        T z

    cdef cppclass sxz[T]:
        sxz(T, T) except +
        T x
        T z

    cdef cppclass siv[T]:
        siv(size_t, T) except +
        size_t i
        T v

    cdef cppclass siv2[T]:
        siv2(size_t, T) except +
        size_t i
        T v
        T v2

    cdef cppclass sijv[T]:
        sijv(size_t, T) except +
        size_t i
        size_t j
        T v

cdef extern from "Node3Dn.h" namespace "ttcr" nogil:
    cdef cppclass Node3Dn[T1,T2]:
        pass

cdef extern from "Node3Dnsp.h" namespace "ttcr" nogil:
    cdef cppclass Node3Dnsp[T1,T2]:
        pass

cdef extern from "Node3Dc.h" namespace "ttcr" nogil:
    cdef cppclass Node3Dc[T1,T2]:
        pass

cdef extern from "Node3Dcsp.h" namespace "ttcr" nogil:
    cdef cppclass Node3Dcsp[T1,T2]:
        pass

cdef extern from "Cell.h" namespace "ttcr" nogil:
    cdef cppclass Cell[T,NODE,S]:
        pass

cdef extern from "Node2Dn.h" namespace "ttcr" nogil:
    cdef cppclass Node2Dn[T1,T2]:
        pass

cdef extern from "Node2Dc.h" namespace "ttcr" nogil:
    cdef cppclass Node2Dc[T1,T2]:
        pass

cdef extern from "Node2Dnsp.h" namespace "ttcr" nogil:
    cdef cppclass Node2Dnsp[T1,T2]:
        pass

cdef extern from "Node2Dcsp.h" namespace "ttcr" nogil:
    cdef cppclass Node2Dcsp[T1,T2]:
        pass

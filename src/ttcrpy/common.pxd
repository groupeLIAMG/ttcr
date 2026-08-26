

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

    cdef cppclass siv4[T]:
        siv4() except +
        size_t i
        T v
        T v2
        T v3
        T v4

    cdef cppclass siv5[T]:
        siv5() except +
        size_t i
        T v
        T v2
        T v3
        T v4
        T v5

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
    cdef cppclass CellElliptical[T,NODE,S]:
        pass
    cdef cppclass CellTiltedElliptical[T,NODE,S]:
        pass
    cdef cppclass CellVTI_PSV[T,NODE,S]:
        pass
    cdef cppclass CellVTI_SH[T,NODE,S]:
        pass
    cdef cppclass CellTTI_PSV[T,NODE,S]:
        pass
    cdef cppclass CellTTI_SH[T,NODE,S]:
        pass
    cdef cppclass CellWeaklyAnelliptical[T,NODE,S]:
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


cdef inline int l_nparams(char iso):
    """Number of medium parameters of an anisotropy model, and so the number of
    blocks of columns the matrix of sensitivities holds.

    The codes are those the grid and mesh classes store in their `iso` member:
    e elliptical, h VTI SH, t tilted elliptical, H TTI SH, w weakly
    anelliptical, p VTI qP/qSV, P TTI qP/qSV, anything else isotropic.

    E is the 3D ellipsoid of CellElliptical3D, which takes a vertical slowness
    and two ratios where the 2D ellipse takes a horizontal slowness and one.
    The upper case marks it as the three-dimensional model, not a tilted one.
    """
    if iso == b'e' or iso == b'h':
        return 2
    elif iso == b't' or iso == b'H' or iso == b'w' or iso == b'E':
        return 3
    elif iso == b'p':
        return 4
    elif iso == b'P':
        return 5
    return 1

# -*- coding: utf-8 -*-
"""
Raytracing on rectilinear grids

This module contains two classes to perform traveltime computation and
raytracing on rectilinear grids:
    - `Grid2d` for 2D media
    - `Grid3d` for 3D media

Three algorithms are implemented
    - the Shortest-Path Method
    - the Fast-Sweeping Method
    - the Dynamic Shortest-Path Method

Slowness model can be defined in two ways:
    1) slowness constant within the voxels of the grid (the default)
    2) slowness defined at nodes of the grid

This code is part of ttcr ( https://github.com/groupeLIAMG/ttcr )
"""

# distutils: language = c++

import numpy as np
cimport numpy as np
import scipy.sparse as sp

import vtk
from vtk.util import numpy_support

from ttcrpy.rgrid cimport Grid3D, Grid3Drcfs, Grid3Drcsp, Grid3Drcdsp, \
    Grid3Drnfs, Grid3Drnsp, Grid3Drndsp, Grid2D, Grid2Drc, Grid2Drn, \
    Grid2Drcsp, Grid2Drcfs, Grid2Drnsp, Grid2Drnfs

cdef extern from "verbose.h" namespace "ttcr" nogil:
    void setVerbose(int)

def set_verbose(v):
    """
    Set verbosity level for C++ code

    Parameters
    ----------
    v: int
        verbosity level
    """
    setVerbose(v)

cdef class Grid3d:
    """
    class to perform raytracing with 3D rectilinear grids

    Attributes
    ----------
    x: np.ndarray
        node coordinates along x
    y: np.ndarray
        node coordinates along y
    z: np.ndarray
        node coordinates along z
    dx: float
        node separation along x
    dy: float
        node separation along y
    dz: float
        node separation along z
    shape: (int, int, int)
        number of parameters along each dimension
    nparams: int
        total number of parameters for grid
    n_threads: int
        number of threads for raytracing

    Constructor:

    Grid3d(x, y, z, n_threads=1, cell_slowness=1, method='FSM', tt_from_rp=1, interp_vel=0, eps=1.e-15, maxit=20, weno=1, nsnx=5, nsny=5, nsnz=5, n_secondary=2, n_tertiary=2, radius_factor_tertiary=3.0, translate_grid=False) -> Grid3d

        Parameters
        ----------
        x : np.ndarray
            node coordinates along x
        y : np.ndarray
            node coordinates along y
        z : np.ndarray
            node coordinates along z
        n_threads : int
            number of threads for raytracing (default is 1)
        cell_slowness : bool
            slowness defined for cells (True) or nodes (False) (default is 1)
        method : string
            raytracing method (default is FSM)
                - 'FSM' : fast marching method
                - 'SPM' : shortest path method
                - 'DSPM' : dynamic shortest path
        tt_from_rp : bool
            compute traveltimes from raypaths (FSM or DSPM only) (default is 1)
        interp_vel : bool
            interpolate velocity instead of slowness at nodes (for
            cell_slowness == False or FSM) (defauls is False)
        eps : double
            convergence criterion (FSM) (default is 1e-15)
        maxit : int
            max number of sweeping iterations (FSM) (default is 20)
        weno : bool
            use 3rd order weighted essentially non-oscillatory operator (FSM)
            (default is True)
        nsnx : int
            number of secondary nodes in x (SPM) (default is 5)
        nsny : int
            number of secondary nodes in y (SPM) (default is 5)
        nsnz : int
            number of secondary nodes in z (SPM) (default is 5)
        n_secondary : int
            number of secondary nodes (DSPM) (default is 2)
        n_tertiary : int
            number of tertiary nodes (DSPM) (default is 2)
        radius_factor_tertiary : double
            multiplication factor used to compute radius of sphere around source
            that includes tertiary nodes (DSPM).  The radius is the average edge
            length multiplied by this factor (default is 3)
        translate_grid : bool
            Translate the grid such that origin is (0, 0, 0) to perform
            computations, which may increase accuracy when large values, e.g.
            UTM coordinates, are used.  When raytracing, src and rcv should be
            given in the original system, and output raypath coordinates are
            also given in the original system (default if False)
    """
    cdef vector[double] _x
    cdef vector[double] _y
    cdef vector[double] _z
    cdef double _dx
    cdef double _dy
    cdef double _dz
    cdef bool cell_slowness
    cdef size_t _n_threads
    cdef char method
    cdef bool tt_from_rp
    cdef bool interp_vel
    cdef bool translate_grid
    cdef double eps
    cdef int maxit
    cdef bool weno
    cdef uint32_t nsnx
    cdef uint32_t nsny
    cdef uint32_t nsnz
    cdef uint32_t n_secondary
    cdef uint32_t n_tertiary
    cdef double radius_factor_tertiary
    cdef Grid3D[double, uint32_t]* grid

    def __cinit__(self, np.ndarray[np.double_t, ndim=1] x,
                  np.ndarray[np.double_t, ndim=1] y,
                  np.ndarray[np.double_t, ndim=1] z,
                  size_t n_threads=1,
                  bool cell_slowness=1, str method='FSM',
                  bool tt_from_rp=1, bool interp_vel=0,
                  double eps=1.e-15, int maxit=20, bool weno=1,
                  uint32_t nsnx=5, uint32_t nsny=5, uint32_t nsnz=5,
                  uint32_t n_secondary=2, uint32_t n_tertiary=2,
                  double radius_factor_tertiary=3.0,
                  bool translate_grid=0):

        cdef uint32_t nx = x.size-1
        cdef uint32_t ny = y.size-1
        cdef uint32_t nz = z.size-1
        self._dx = x[1] - x[0]
        self._dy = y[1] - y[0]
        self._dz = z[1] - z[0]
        cdef double xmin = x[0]
        cdef double ymin = y[0]
        cdef double zmin = z[0]
        self.cell_slowness = cell_slowness
        self._n_threads = n_threads
        self.tt_from_rp = tt_from_rp
        self.interp_vel = interp_vel
        self.eps = eps
        self.maxit = maxit
        self.weno = weno
        self.nsnx = nsnx
        self.nsny = nsny
        self.nsnz = nsnz
        self.n_secondary = n_secondary
        self.n_tertiary = n_tertiary
        self.radius_factor_tertiary = radius_factor_tertiary
        self.translate_grid = translate_grid

        cdef use_edge_length = True

        if method == 'FSM':
            if np.abs(self._dx - self._dy)>0.000001 or np.abs(self._dx - self._dz)>0.000001:
                raise ValueError('FSM: Grid cells must be cubic')

        for val in x:
            self._x.push_back(val)
        for val in y:
            self._y.push_back(val)
        for val in z:
            self._z.push_back(val)
        self._x.shrink_to_fit()
        self._y.shrink_to_fit()
        self._z.shrink_to_fit()

        if cell_slowness:
            if method == 'FSM':
                self.method = b'f'
                self.grid = new Grid3Drcfs[double,uint32_t](nx, ny, nz, self._dx,
                                                            xmin, ymin, zmin,
                                                            eps, maxit, weno,
                                                            tt_from_rp, interp_vel,
                                                            n_threads,
                                                            translate_grid)
            elif method == 'SPM':
                self.method = b's'
                self.grid = new Grid3Drcsp[double,uint32_t,Cell[double,Node3Dcsp[double,uint32_t],sxyz[double]]](nx, ny, nz,
                                                            self._dx, self._dy, self._dz,
                                                            xmin, ymin, zmin,
                                                            nsnx, nsny, nsnz,
                                                            tt_from_rp,
                                                            n_threads,
                                                            translate_grid)
            elif method == 'DSPM':
                self.method = b'd'
                self.grid = new Grid3Drcdsp[double,uint32_t,Cell[double,Node3Dc[double,uint32_t],sxyz[double]]](nx, ny, nz,
                                                           self._dx, self._dy, self._dz,
                                                           xmin, ymin, zmin,
                                                           n_secondary, tt_from_rp,
                                                           n_tertiary, radius_factor_tertiary,
                                                           use_edge_length, n_threads,
                                                           translate_grid)
            else:
                raise ValueError('Method {0:s} undefined'.format(method))
        else:
            if method == 'FSM':
                self.method = b'f'
                self.grid = new Grid3Drnfs[double,uint32_t](nx, ny, nz, self._dx,
                                                            xmin, ymin, zmin,
                                                            eps, maxit, weno,
                                                            tt_from_rp, interp_vel,
                                                            n_threads,
                                                            translate_grid)
            elif method == 'SPM':
                self.method = b's'
                self.grid = new Grid3Drnsp[double,uint32_t](nx, ny, nz,
                                                            self._dx, self._dy, self._dz,
                                                            xmin, ymin, zmin,
                                                            nsnx, nsny, nsnz,
                                                            tt_from_rp, interp_vel,
                                                            n_threads,
                                                            translate_grid)
            elif method == 'DSPM':
                self.method = b'd'
                self.grid = new Grid3Drndsp[double,uint32_t](nx, ny, nz,
                                                             self._dx, self._dy, self._dz,
                                                             xmin, ymin, zmin,
                                                             n_secondary, tt_from_rp,
                                                             n_tertiary, radius_factor_tertiary,
                                                             interp_vel, use_edge_length,
                                                             n_threads,
                                                             translate_grid)
            else:
                raise ValueError('Method {0:s} undefined'.format(method))

    def __dealloc__(self):
        del self.grid

    def __reduce__(self):
        if self.method == b'f':
            method = 'FSM'
        elif self.method == b's':
            method = 'SPM'
        elif self.method == b'd':
            method = 'DSPM'

        constructor_params = (self.n_threads, self.cell_slowness, method,
                              self.tt_from_rp, self.interp_vel, self.eps,
                              self.maxit, self.weno, self.nsnx, self.nsny,
                              self.nsnz, self.n_secondary, self.n_tertiary,
                              self.radius_factor_tertiary, self.translate_grid)
        return (_rebuild3d, (self.x, self.y, self.z, constructor_params))

    @property
    def x(self):
        """np.ndarray: node coordinates along x"""
        tmp = np.empty((self._x.size(),))
        cdef int n
        for n in range(self._x.size()):
            tmp[n] = self._x[n]
        return tmp

    @property
    def y(self):
        """np.ndarray: node coordinates along y"""
        tmp = np.empty((self._y.size(),))
        cdef int n
        for n in range(self._y.size()):
            tmp[n] = self._y[n]
        return tmp

    @property
    def z(self):
        """np.ndarray: node coordinates along z"""
        tmp = np.empty((self._z.size(),))
        cdef int n
        for n in range(self._z.size()):
            tmp[n] = self._z[n]
        return tmp

    @property
    def dx(self):
        """float: node separation along x"""
        return self._dx

    @property
    def dy(self):
        """float: node separation along y"""
        return self._dy

    @property
    def dz(self):
        """float: node separation along z"""
        return self._dz

    @property
    def shape(self):
        """:obj:`list` of :obj:`int`: number of parameters along each dimension"""
        if self.cell_slowness:
            return (self._x.size()-1, self._y.size()-1, self._z.size()-1)
        else:
            return (self._x.size(), self._y.size(), self._z.size())

    @property
    def n_threads(self):
        """int: number of threads for raytracing"""
        return self._n_threads

    @property
    def nparams(self):
        """int: total number of parameters for grid"""
        if self.cell_slowness:
            return (self._x.size()-1) * (self._y.size()-1) * (self._z.size()-1)
        else:
            return self._x.size() * self._y.size() * self._z.size()

    def set_use_thread_pool(self, use_thread_pool):
        """
        set_use_thread_pool(use_thread_pool)

        Set option to use thread pool instead of parallel loop

        Parameters
        ----------
        use_thread_pool : bool
            option value
        """
        self.grid.setUsePool(use_thread_pool)

    def set_traveltime_from_raypath(self, traveltime_from_raypath):
        """
        set_traveltime_from_raypath(ttrp)

        Set option to compute traveltime using raypath

        Parameters
        ----------
        ttrp : bool
            option value
        """
        self.grid.setTraveltimeFromRaypath(traveltime_from_raypath)

    def get_number_of_nodes(self):
        """
        Returns
        -------
        int:
            number of nodes in grid
        """
        return self._x.size() * self._y.size() * self._z.size()

    def get_number_of_cells(self):
        """
        Returns
        -------
        int:
            number of cells in grid
        """
        return (self._x.size()-1) * (self._y.size()-1) * (self._z.size()-1)

    def get_grid_traveltimes(self, thread_no=0):
        """
        get_grid_traveltimes(thread_no=0)

        Obtain traveltimes computed at primary grid nodes

        Parameters
        ----------
        thread_no : int
            thread used to computed traveltimes (default is 0)

        Returns
        -------
        tt: np ndarray, shape (nx, ny, nz)
            traveltimes
        """
        if thread_no >= self._n_threads:
            raise ValueError('Thread number is larger than number of threads')
        cdef vector[double] tmp
        cdef int n
        self.grid.getTT(tmp, thread_no)
        tt = np.empty((tmp.size(),), order='F')
        for n in range(tmp.size()):
            tt[n] = tmp[n]
        shape = (self._x.size(), self._y.size(), self._z.size())
        return tt.reshape(shape)

    def ind(self, int i, int j, int k):
        """
        ind(i, j, k)

        Return node index

        Parameters
        ----------
        i : int
            index of node along x
        j : int
            index of node along y
        k : int
            index of node along z

        Returns
        -------
        int:
            node index for a "flattened" grid
        """
        return (i*self._y.size() + j)*self._z.size() + k

    def indc(self, int i, int j, int k):
        """
        return cell index

        Parameters
        ----------
        i : int
            index of cell along x
        j : int
            index of cell along y
        k : int
            index of cell along z

        Returns
        -------
        int:
            cell index for a "flattened" grid
        """
        return (i*(self._y.size()-1) + j)*(self._z.size()-1) + k

    cdef uint32_t _f2c_ind(self, uint32_t ind):
        """Convert cell index from 'F' order to 'C' order"""
        cdef uint32_t i, j, k
        k = <uint32_t>(ind / ((self._x.size()-1) * (self._y.size()-1)))
        j = <uint32_t>((ind - k * (self._x.size()-1) * (self._y.size()-1)) / (self._x.size()-1))
        i = <uint32_t>(ind - ((k * (self._y.size()-1) + j) * (self._x.size()-1)))
        return (i * (self._y.size()-1) + j) * (self._z.size()-1) + k

    def is_outside(self, np.ndarray[np.double_t, ndim=2] pts):
        """
        is_outside(pts)

        Check if points are outside grid

        Parameters
        ----------
        pts : np ndarray, shape(npts, 3)
            coordinates of points to check

        Returns
        -------
        bool:
            True if at least one point outside grid
        """
        return ( np.min(pts[:,0]) < self._x.front() or np.max(pts[:,0]) > self._x.back() or
                np.min(pts[:,1]) < self._y.front() or np.max(pts[:,1]) > self._y.back() or
                np.min(pts[:,2]) < self._z.front() or np.max(pts[:,2]) > self._z.back() )

    def get_slowness(self):
        """Returns slowness of grid

        Returns
        -------
        slowness : np ndarray, shape (nx, ny, nz)
        """
        cdef int i
        cdef vector[double] slown
        self.grid.getSlowness(slown)
        nx, ny, nz = self.shape
        slowness = np.ndarray((nx*ny*nz,), order='F')
        for i in range(slown.size()):
            slowness[i] = slown[i]
        return slowness.reshape((nx, ny, nz))

    def set_slowness(self, slowness):
        """
        set_slowness(slowness)

        Assign slowness to grid

        Parameters
        ----------
        slowness : np ndarray, shape (nx, ny, nz)
            slowness may also have been flattened (with default 'C' order)
        """
        if self.cell_slowness:
            nx = self._x.size()-1
            ny = self._y.size()-1
            nz = self._z.size()-1
        else:
            nx = self._x.size()
            ny = self._y.size()
            nz = self._z.size()
        if slowness.size != nx*ny*nz:
            raise ValueError('Slowness vector has wrong size')

        cdef vector[double] slown
        cdef int i
        if slowness.ndim == 3:
            if slowness.shape != (nx, ny, nz):
                raise ValueError('Slowness has wrong shape')
            tmp = slowness.flatten('F')
            for i in range(nx*ny*nz):
                slown.push_back(tmp[i])
        elif slowness.ndim == 1:
            # slowness is in 'C' order and we must pass it in 'F' order
            tmp = slowness.reshape((nx, ny, nz)).flatten('F')
            for i in range(nx*ny*nz):
                slown.push_back(tmp[i])
        else:
            raise ValueError('Slowness must be 1D or 3D ndarray')
        self.grid.setSlowness(slown)

    def set_velocity(self, velocity):
        """
        set_velocity(slowness)

        Assign velocity to grid

        Parameters
        ----------
        velocity : np ndarray, shape (nx, ny, nz)
            velocity may also have been flattened (with default 'C' order)
        """
        if self.cell_slowness:
            nx = self._x.size()-1
            ny = self._y.size()-1
            nz = self._z.size()-1
        else:
            nx = self._x.size()
            ny = self._y.size()
            nz = self._z.size()
        if velocity.size != nx*ny*nz:
            raise ValueError('velocity vector has wrong size')

        cdef vector[double] slown
        cdef int i
        if velocity.ndim == 3:
            if velocity.shape != (nx, ny, nz):
                raise ValueError('velocity has wrong shape')
            tmp = velocity.flatten('F')
            for i in range(nx*ny*nz):
                slown.push_back(1./tmp[i])
        elif velocity.ndim == 1:
            # slowness is in 'C' order and we must pass it in 'F' order
            tmp = velocity.reshape((nx, ny, nz)).flatten('F')
            for i in range(nx*ny*nz):
                slown.push_back(1./tmp[i])
        else:
            raise ValueError('velocity must be 1D or 3D ndarray')
        self.grid.setSlowness(slown)

    def compute_D(self, coord):
        """
        compute_D(coord)

        Return matrix of interpolation weights for velocity data points
        constraint

        Parameters
        ----------
        coord : np.ndarray, shape (npts, 3)
            coordinates of data points

        Returns
        -------
        D : scipy csr_matrix, shape (npts, nparams)
            Matrix of interpolation weights
        """
        if self.is_outside(coord):
            raise ValueError('Velocity data point outside grid')

        # D is npts x nparams

        cdef Py_ssize_t n, i, j, k, i1, i2, j1, j2, k1, k2, ii

        if self.cell_slowness:
            # for each point in coord, we have 1 values in D
            ivec = np.arange(coord.shape[0], dtype=np.int64)
            jvec = np.zeros(ivec.shape, dtype=np.int64)
            vec = np.ones(ivec.shape)
            for n in np.arange(coord.shape[0]):
                i = int((coord[n,0]-self.x[0])/self._dx )
                j = int((coord[n,1]-self.y[0])/self._dx )
                k = int((coord[n,2]-self.z[0])/self._dx )

                jvec[n] = self.indc(i,j,k)

            return sp.csr_matrix((vec, (ivec,jvec)),
                                 shape=(coord.shape[0], self.get_number_of_cells()))
        else:
            # for each point in coord, we have 8 values in D
            ivec = np.kron(np.arange(coord.shape[0], dtype=np.int64),np.ones(8, dtype=np.int64))
            jvec = np.zeros(ivec.shape, dtype=np.int64)
            vec = np.zeros(ivec.shape)
            for n in np.arange(coord.shape[0]):
                i1 = int(1.e-6 + (coord[n,0]-self._x[0])/self._dx )
                i2 = i1 + 1
                j1 = int(1.e-6 + (coord[n,1]-self._y[0])/self._dy )
                j2 = j1 + 1
                k1 = int(1.e-6 + (coord[n,2]-self._z[0])/self._dz )
                k2 = k1 + 1

                ii = 0
                for i in (i1, i2):
                    for j in (j1, j2):
                        for k in (k1, k2):
                            jvec[n*8+ii] = self.ind(i,j,k)
                            vec[n*8+ii] = ((1. - np.abs(coord[n,0]-self._x[i])/self._dx) *
                                           (1. - np.abs(coord[n,1]-self._y[j])/self._dy) *
                                           (1. - np.abs(coord[n,2]-self._z[k])/self._dz))
                            ii += 1

            return sp.csr_matrix((vec, (ivec,jvec)),
                                 shape=(coord.shape[0], self.get_number_of_nodes()))

    def compute_K(self):
        """
        Compute smoothing matrices (2nd order derivative)

        Returns
        -------
        Kx, Ky, Kz : :obj:`tuple` of :obj:`csr_matrix`
            matrices for derivatives along x, y, & z
        """
        # central operator f"(x) = (f(x+h)-2f(x)+f(x-h))/h^2
        # forward operator f"(x) = (f(x+2h)-2f(x+h)+f(x))/h^2
        # backward operator f"(x) = (f(x)-2f(x-h)+f(x-2h))/h^2

        cdef Py_ssize_t i, j, k

        nx, ny, nz = self.shape
        # Kx
        iK = np.kron(np.arange(nx*ny*nz, dtype=np.int64), np.ones(3,dtype=np.int64))
        val = np.tile(np.array([1., -2., 1.]), nx*ny*nz) / (self._dx*self._dx)
        # i=0 -> forward op
        jK = np.vstack((np.arange(ny*nz,dtype=np.int64),
                        np.arange(ny*nz, 2*ny*nz,dtype=np.int64),
                        np.arange(2*ny*nz, 3*ny*nz,dtype=np.int64))).T.flatten()

        for i in np.arange(1,nx-1):
            jK = np.hstack((jK,
                            np.vstack((np.arange((i-1)*ny*nz,i*ny*nz, dtype=np.int64),
                                       np.arange(i*ny*nz, (i+1)*ny*nz,dtype=np.int64),
                                       np.arange((i+1)*ny*nz, (i+2)*ny*nz,dtype=np.int64))).T.flatten()))
        # i=nx-1 -> backward op
        jK = np.hstack((jK,
                        np.vstack((np.arange((nx-3)*ny*nz, (nx-2)*ny*nz, dtype=np.int64),
                                   np.arange((nx-2)*ny*nz, (nx-1)*ny*nz,dtype=np.int64),
                                   np.arange((nx-1)*ny*nz, nx*ny*nz,dtype=np.int64))).T.flatten()))

        Kx = sp.csr_matrix((val,(iK,jK)))

        # j=0 -> forward op

        jK = np.vstack((np.arange(nz,dtype=np.int64),
                        np.arange(nz,2*nz,dtype=np.int64),
                        np.arange(2*nz,3*nz,dtype=np.int64))).T.flatten()
        val = np.tile(np.array([1., -2., 1.]), nx*ny*nz) / (self._dy*self._dy)
        for j in np.arange(1,ny-1):
            jK = np.hstack((jK,
                            np.vstack((np.arange((j-1)*nz,j*nz,dtype=np.int64),
                                       np.arange(j*nz,(j+1)*nz,dtype=np.int64),
                                       np.arange((j+1)*nz,(j+2)*nz,dtype=np.int64))).T.flatten()))
        # j=ny-1 -> backward op
        jK = np.hstack((jK,
                        np.vstack((np.arange((ny-3)*nz,(ny-2)*nz,dtype=np.int64),
                                   np.arange((ny-2)*nz,(ny-1)*nz,dtype=np.int64),
                                   np.arange((ny-1)*nz,ny*nz,dtype=np.int64))).T.flatten()))
        tmp = jK.copy()
        for i in np.arange(1,nx):
            jK = np.hstack((jK, i*ny*nz+tmp))

        Ky = sp.csr_matrix((val,(iK,jK)))

        # k=0
        val = np.tile(np.array([1., -2., 1.]), nx*ny*nz) / (self._dz*self._dz)
        jK = np.arange(3,dtype=np.int64)
        for k in np.arange(1,nz-1):
            jK = np.hstack((jK,(k-1)+np.arange(3,dtype=np.int64)))
        # k=nz-1
        jK = np.hstack((jK, (nz-3)+np.arange(3,dtype=np.int64)))

        tmp = jK.copy()
        for j in np.arange(1,ny):
            jK = np.hstack((jK, j*nz+tmp))

        tmp = jK.copy()
        for i in np.arange(1,nx):
            jK = np.hstack((jK, i*ny*nz+tmp))

        Kz = sp.csr_matrix((val,(iK,jK)))

        return Kx, Ky, Kz

    def get_s0(self, hypo, slowness=None):
        """
        get_s0(hypo, slowness=None)

        Return slowness at source points

        Parameters
        ----------
        hypo : np.ndarray with 5 columns
            hypo holds source information, i.e.
                - 1st column is event ID number
                - 2nd column is origin time
                - 3rd column is source easting
                - 4th column is source northing
                - 5th column is source elevation

        slowness : np ndarray, shape (nx, ny, nz) (optional)
            slowness at grid nodes or cells (depending on cell_slowness)
            slowness may also have been flattened (with default 'C' order)

        Returns
        -------
        s0 : np.ndarray
            slowness at source points
        """

        cdef Py_ssize_t n, nn

        if hypo.shape[1] != 5:
            raise ValueError('hypo should be npts x 5')
        src = hypo[:,2:5]
        evID = hypo[:,0]
        eid = np.sort(np.unique(evID))
        nTx = len(eid)

        if slowness is not None:
            self.set_slowness(slowness)

        i0 = np.empty((nTx,), dtype=np.int64)
        for n in range(nTx):
            for nn in range(evID.size):
                if eid[n] == evID[nn]:
                    i0[n] = nn
                    break

        cdef vector[vector[sxyz[double]]] vTx
        vTx.resize(nTx)

        iTx = []
        i0 = 0
        for n in range(nTx):
            for nn in range(evID.size):
                if eid[n] == evID[nn]:
                    i0 = nn
                    break
            vTx[n].push_back(sxyz[double](src[i0,0], src[i0,1], src[i0,2]))

        for i in eid:
            ii = evID == i
            iTx.append(np.nonzero(ii)[0])

        s0 = np.zeros((src.shape[0],))

        for n in range(nTx):
            s = 0.0
            for nn in range(vTx[n].size()):
                s += self.grid.computeSlowness(vTx[n][nn])
            s0[iTx[n]] = s/vTx[n].size()
        return s0

    def raytrace(self, source, rcv, slowness=None, thread_no=None,
                 aggregate_src=False, compute_L=False, compute_M=False,
                 return_rays=False):
        """
        raytrace(source, rcv, slowness=None, thread_no=None,
                 aggregate_src=False, compute_L=False, compute_M=False,
                 return_rays=False) -> tt, rays, M, L

        Perform raytracing

        Parameters
        ----------
        source : 2D np.ndarray with 3, 4 or 5 columns
            see notes below
        rcv : 2D np.ndarray with 3 columns
            Columns correspond to x, y and z coordinates
        slowness : np ndarray, shape (nx, ny, nz) (None by default)
            slowness at grid nodes or cells (depending on cell_slowness)
            slowness may also have been flattened (with default 'C' order)
            if None, slowness must have been assigned previously
        thread_no : int (None by default)
            Perform calculations in thread number "thread_no"
            if None, attempt to run in parallel if warranted by number of
            sources and value of n_threads in constructor
        aggregate_src : bool (False by default)
            if True, all source coordinates belong to a single event
        compute_L : bool (False by default)
            Compute matrices of partial derivative of travel time w/r to slowness
        compute_M : bool (False by default)
            Compute matrices of partial derivative of travel time w/r to velocity
            Note : compute_M and compute_L are mutually exclusive
        return_rays : bool (False by default)
            Return raypaths

        Returns
        -------
        tt : np.ndarray
            travel times for the appropriate source-rcv  (see Notes below)
        rays : :obj:`list` of :obj:`np.ndarray`
            Coordinates of segments forming raypaths (if return_rays is True)
        M : :obj:`list` of :obj:`csr_matrix`
            matrices of partial derivative of travel time w/r to velocity.
            the number of matrices is equal to the number of sources
        L : scipy csr_matrix
            Matrix of partial derivative of travel time w/r to slowness.
            if input argument source has 5 columns, L is a list of matrices and
            the number of matrices is equal to the number of sources
            otherwise, L is a single csr_matrix

        Notes
        -----
        If source has 3 columns:
            - Columns correspond to x, y and z coordinates
            - Origin time (t0) is 0 for all points
        If source has 4 columns:
            - 1st column corresponds to origin times
            - 2nd, 3rd & 4th columns correspond to x, y and z coordinates
        If source has 5 columns:
            - 1st column corresponds to event ID
            - 2nd column corresponds to origin times
            - 3rd, 4th & 5th columns correspond to x, y and z coordinates

        For the latter case (5 columns), source and rcv should contain the same
        number of rows, each row corresponding to a source-receiver pair.
        For the 2 other cases, source and rcv can contain the same number of
        rows, each row corresponding to a source-receiver pair, or the number
        of rows may differ if aggregate_src is True or if all rows in source
        are identical.
        """

        # check input data consistency

        if source.ndim != 2 or rcv.ndim != 2:
            raise ValueError('source and rcv should be 2D arrays')

        if self.method == b'd' and aggregate_src:
            raise ValueError('Cannot aggregate source with DSPM raytracing')

        if compute_L and compute_M:
            raise ValueError('compute_L and compute_M are mutually exclusive')

        if self.cell_slowness and compute_M:
            raise NotImplementedError('compute_M not defined for grids with slowness defined for cells')

        if compute_L and not self.cell_slowness:
            raise NotImplementedError('compute_L defined only for grids with slowness defined for cells')

        evID = None
        if source.shape[1] == 5:
            src = source[:,2:5]
            t0 = source[:,1]
            evID = source[:,0]
            eid = np.sort(np.unique(evID))
            nTx = len(eid)
        elif source.shape[1] == 3:
            src = source
            _, ind = np.unique(source, axis=0, return_index=True)
            Tx = source[np.sort(ind), :]     # this to keep the original order
            t0 = np.zeros((Tx.shape[0], 1))
            nTx = Tx.shape[0]
        elif source.shape[1] == 4:
            src = source[:,1:4]
            _, ind = np.unique(source, axis=0, return_index=True)
            tmp = source[np.sort(ind), :]    # this to keep the original order
            nTx = tmp.shape[0]
            Tx = tmp[:,1:4]
            t0 = tmp[:,0]
        else:
            raise ValueError('source should be either nsrc x 3, 4 or 5')

        if src.shape[1] != 3 or rcv.shape[1] != 3:
            raise ValueError('src and rcv should be ndata x 3')

        if self.is_outside(src):
            raise ValueError('Source point outside grid')

        if self.is_outside(rcv):
            raise ValueError('Receiver outside grid')

        if slowness is not None:
            self.set_slowness(slowness)

        cdef vector[vector[sxyz[double]]] vTx
        cdef vector[vector[sxyz[double]]] vRx
        cdef vector[vector[double]] vt0
        cdef vector[vector[double]] vtt

        cdef vector[vector[vector[sxyz[double]]]] r_data
        cdef vector[vector[vector[siv[double]]]] l_data
        cdef vector[vector[vector[sijv[double]]]] m_data
        cdef size_t thread_nb

        cdef int i, j, k, n, nn, MM, NN

        vTx.resize(nTx)
        vRx.resize(nTx)
        vt0.resize(nTx)
        vtt.resize(nTx)
        if compute_L:
            l_data.resize(nTx)
        if compute_M:
            m_data.resize(nTx)
        if return_rays:
            r_data.resize(nTx)

        iRx = []
        if evID is None:
            if nTx == 1:
                vTx[0].push_back(sxyz[double](src[0,0], src[0,1], src[0,2]))
                for r in rcv:
                    vRx[0].push_back(sxyz[double](r[0], r[1], r[2]))
                vt0[0].push_back(t0[0])
                vtt[0].resize(rcv.shape[0])
                iRx.append(np.arange(rcv.shape[0]))
            elif aggregate_src:
                for t in Tx:
                    vTx[0].push_back(sxyz[double](t[0], t[1], t[2]))
                for t in t0:
                    vt0[0].push_back(t)
                for r in rcv:
                    vRx[0].push_back(sxyz[double](r[0], r[1], r[2]))
                vtt[0].resize(rcv.shape[0])
                nTx = 1
                iRx.append(np.arange(rcv.shape[0]))
            else:
                if src.shape != rcv.shape:
                    raise ValueError('src and rcv should be of equal size')

                for n in range(nTx):
                    ind = np.sum(Tx[n,:] == src, axis=1) == 3
                    iRx.append(np.nonzero(ind)[0])
                    vTx[n].push_back(sxyz[double](Tx[n,0], Tx[n,1], Tx[n,2]))
                    vt0[n].push_back(t0[n])
                    for r in rcv[ind,:]:
                        vRx[n].push_back(sxyz[double](r[0], r[1], r[2]))
                    vtt[n].resize(vRx[n].size())
        else:
            if src.shape != rcv.shape:
                raise ValueError('src and rcv should be of equal size')

            i0 = 0
            for n in range(nTx):
                for nn in range(evID.size):
                    if eid[n] == evID[nn]:
                        i0 = nn
                        break
                vTx[n].push_back(sxyz[double](src[i0,0], src[i0,1], src[i0,2]))
                vt0[n].push_back(t0[i0])

            for i in eid:
                ii = evID == i
                iRx.append(np.nonzero(ii)[0])

            for n in range(nTx):
                for r in rcv[iRx[n],:]:
                    vRx[n].push_back(sxyz[double](r[0], r[1], r[2]))
                vtt[n].resize(vRx[n].size())

        # for n in range(nTx):
        #     print(n, '\nTx')
        #     for nn in range(vTx[n].size()):
        #         print('  ',vTx[n][nn].x, vTx[n][nn].y, vTx[n][nn].z)
        #     print('Rx')
        #     for nn in range(vRx[n].size()):
        #         print('  ',vRx[n][nn].x, vRx[n][nn].y, vRx[n][nn].z)
        #     print('t0')
        #     for nn in range(vt0[n].size()):
        #         print('  ',vt0[n][nn])
        #     print('tt')
        #     print('  size of tt = ',vtt[n].size())

        tt = np.zeros((rcv.shape[0],))
        if self._n_threads == 1:
            if compute_L==False and compute_M==False and return_rays==False:
                for n in range(nTx):
                    self.grid.raytrace(vTx[n], vt0[n], vRx[n], vtt[n], 0)
            elif compute_M and return_rays:
                for n in range(nTx):
                    self.grid.raytrace(vTx[n], vt0[n], vRx[n], vtt[n], r_data[n], m_data[n], 0)
            elif compute_L and return_rays:
                for n in range(nTx):
                    self.grid.raytrace(vTx[n], vt0[n], vRx[n], vtt[n], r_data[n], l_data[n], 0)
            elif compute_L:
                for n in range(nTx):
                    self.grid.raytrace(vTx[n], vt0[n], vRx[n], vtt[n], l_data[n], 0)
            elif compute_M:
                for n in range(nTx):
                    self.grid.raytrace(vTx[n], vt0[n], vRx[n], vtt[n], m_data[n], 0)
            else:
                for n in range(nTx):
                    self.grid.raytrace(vTx[n], vt0[n], vRx[n], vtt[n], r_data[n], 0)

        elif thread_no is not None:
            # we should be here for just one event
            assert nTx == 1
            # normally we should not need to compute M or L
            assert compute_L is False
            assert compute_M is False
            thread_nb = thread_no

            if return_rays:
                self.grid.raytrace(vTx[0], vt0[0], vRx[0], vtt[0], r_data[0], thread_nb)
                for nt in range(vtt[0].size()):
                    tt[nt] = vtt[0][nt]
                rays = []
                for n2 in range(vRx.size()):
                    r = np.empty((r_data[0][n2].size(), 3))
                    for nn in range(r_data[0][n2].size()):
                        r[nn, 0] = r_data[0][n2][nn].x
                        r[nn, 1] = r_data[0][n2][nn].y
                        r[nn, 2] = r_data[0][n2][nn].z
                    rays.append(r)
                return tt, rays

            else:
                self.grid.raytrace(vTx[0], vt0[0], vRx[0], vtt[0], thread_nb)
                for nt in range(vtt[0].size()):
                    tt[nt] = vtt[0][nt]
                return tt

        else:
            if compute_L==False and compute_M==False and return_rays==False:
                self.grid.raytrace(vTx, vt0, vRx, vtt)
            elif compute_M and return_rays:
                self.grid.raytrace(vTx, vt0, vRx, vtt, r_data, m_data)
            elif compute_L and return_rays:
                self.grid.raytrace(vTx, vt0, vRx, vtt, r_data, l_data)
            elif compute_L:
                self.grid.raytrace(vTx, vt0, vRx, vtt, l_data)
            elif compute_M:
                self.grid.raytrace(vTx, vt0, vRx, vtt, m_data)
            else:
                self.grid.raytrace(vTx, vt0, vRx, vtt, r_data)

        for n in range(nTx):
            for nt in range(vtt[n].size()):
                tt[iRx[n][nt]] = vtt[n][nt]

        if return_rays:
            rays = [ [0.0] for n in range(rcv.shape[0])]
            for n in range(nTx):
                r = [ [0.0] for i in range(vRx[n].size())]
                for n2 in range(vRx[n].size()):
                    r[n2] = np.empty((r_data[n][n2].size(), 3))
                    for nn in range(r_data[n][n2].size()):
                        r[n2][nn, 0] = r_data[n][n2][nn].x
                        r[n2][nn, 1] = r_data[n][n2][nn].y
                        r[n2][nn, 2] = r_data[n][n2][nn].z
                for nt in range(vtt[n].size()):
                    rays[iRx[n][nt]] = r[nt]

        if compute_L:

            # we build an array of matrices, for each event
            L = []
            for n in range(nTx):
                nnz = 0
                for ni in range(l_data[n].size()):
                    nnz += l_data[n][ni].size()
                indptr = np.empty((vRx[n].size()+1,), dtype=np.int64)
                indices = np.empty((nnz,), dtype=np.int64)
                val = np.empty((nnz,))

                k = 0
                MM = vRx[n].size()
                NN = self.get_number_of_cells()
                for i in range(MM):
                    indptr[i] = k
                    for j in range(NN):
                        for nn in range(l_data[n][i].size()):
                            if self._f2c_ind(l_data[n][i][nn].i) == j:
                                indices[k] = j
                                val[k] = l_data[n][i][nn].v
                                k += 1

                indptr[MM] = k
                L.append( sp.csr_matrix((val, indices, indptr), shape=(MM,NN)) )

            if evID is None:
                # we want a single matrix
                tmp = sp.vstack(L)
                itmp = []
                for n in range(nTx):
                    for nt in range(vtt[n].size()):
                        itmp.append(iRx[n][nt])
                L = tmp[itmp,:]

        if compute_M:
            # we return array of matrices, one for each event
            M = []
            for n in range(nTx):
                nnz = 0
                for ni in range(m_data[n].size()):
                    nnz += m_data[n][ni].size()
                indptr = np.empty((vRx[n].size()+1,), dtype=np.int64)
                indices = np.empty((nnz,), dtype=np.int64)
                val = np.empty((nnz,))

                k = 0
                MM = vRx[n].size()
                NN = self.get_number_of_nodes()
                for i in range(MM):
                    indptr[i] = k
                    for j in range(NN):
                        for nn in range(m_data[n][i].size()):
                            if m_data[n][i][nn].i == i and m_data[n][i][nn].j == j:
                                indices[k] = j
                                val[k] = m_data[n][i][nn].v
                                k += 1

                indptr[MM] = k
                M.append( sp.csr_matrix((val, indices, indptr), shape=(MM,NN)) )

        if compute_L==False and compute_M==False and return_rays==False:
            return tt
        elif compute_M and return_rays:
            return tt, rays, M
        elif compute_L and return_rays:
            return tt, rays, L
        elif compute_L:
            return tt, L
        elif compute_M:
            return tt, M
        else:
            return tt, rays

    def to_vtk(self, fields, filename):
        """
        to_vtk(fields, filename)

        Save grid variables and/or raypaths to VTK format

        Parameters
        ----------
        fields: dict
            dict of variables to save to file. Variables should be np.ndarray of
            size equal to either the number of nodes of the number of cells of
            the grid, or a list of raypath coordinates.
        filename: str
            Name of file without extension for saving (extension vtr will be
            added).  Raypaths are saved in separate files, and filename will
            be appended by the dict key and have a vtp extension.

        Notes
        -----
        VTK files can be visualized with Paraview (https://www.paraview.org)
        """
        xCoords = numpy_support.numpy_to_vtk(self.x)
        yCoords = numpy_support.numpy_to_vtk(self.y)
        zCoords = numpy_support.numpy_to_vtk(self.z)

        rgrid = vtk.vtkRectilinearGrid()
        rgrid.SetDimensions(self._x.size(), self._y.size(), self._z.size())
        rgrid.SetXCoordinates(xCoords)
        rgrid.SetYCoordinates(yCoords)
        rgrid.SetZCoordinates(zCoords)

        save_grid = False
        for fn in fields:
            data = fields[fn]

            if isinstance(data, list):
                self._save_raypaths(data, filename+'_'+fn+'.vtp')
            else:
                save_grid = True
                scalar = vtk.vtkDoubleArray()
                scalar.SetName(fn)
                scalar.SetNumberOfComponents(1)
                scalar.SetNumberOfTuples(data.size)
                if data.size == self.get_number_of_nodes():
                    shape = (self._x.size(), self._y.size(), self._z.size())
                    if data.ndim == 3:
                        if data.shape != shape:
                            raise ValueError('Field {0:s} has incorrect shape'.format(fn))
                        # VTK stores in 'F' order
                        tmp = data.flatten(order='F')
                    elif data.ndim == 1 or (data.ndim == 2 and data.shape[0] == data.size):
                        # 'C' order assumed, reshape back and flatten to 'F' order
                        tmp = data.reshape(shape).flatten(order='F')
                    else:
                        raise ValueError('Field {0:s} has incorrect ndim ({1:d})'.format(fn, data.ndim))
                    for n in range(data.size):
                        scalar.SetTuple1(n, tmp[n])
                    rgrid.GetPointData().AddArray(scalar)
                elif data.size == self.get_number_of_cells():
                    shape = (self._x.size()-1, self._y.size()-1, self._z.size()-1)
                    if data.ndim == 3:
                        if data.shape != shape:
                            raise ValueError('Field {0:s} has incorrect shape'.format(fn))
                        # VTK stores in 'F' order
                        tmp = data.flatten(order='F')
                    elif data.ndim == 1 or (data.ndim == 2 and data.shape[0] == data.size):
                        # 'C' order assumed, reshape back and flatten to 'F' order
                        tmp = data.reshape(shape).flatten(order='F')
                    else:
                        raise ValueError('Field {0:s} has incorrect ndim ({1:d})'.format(fn, data.ndim))
                    for n in range(data.size):
                        scalar.SetTuple1(n, tmp[n])
                    rgrid.GetCellData().AddArray(scalar)
                else:
                    raise ValueError('Field {0:s} has incorrect size'.format(fn))

        if save_grid:
            writer = vtk.vtkXMLRectilinearGridWriter()
            writer.SetFileName(filename+'.vtr')
            writer.SetInputData(rgrid)
            writer.SetDataModeToBinary()
            writer.Update()

    def _save_raypaths(self, rays, filename):
        polydata = vtk.vtkPolyData()
        cellarray = vtk.vtkCellArray()
        pts = vtk.vtkPoints()
        npts = 0
        for n in range(len(rays)):
            npts += rays[n].shape[0]
        pts.SetNumberOfPoints(npts)
        npts = 0
        for n in range(len(rays)):
            for p in range(rays[n].shape[0]):
                pts.InsertPoint(npts, rays[n][p, 0], rays[n][p, 1], rays[n][p, 2])
                npts += 1
        polydata.SetPoints(pts)
        npts = 0
        for n in range(len(rays)):
            line = vtk.vtkPolyLine()
            line.GetPointIds().SetNumberOfIds(rays[n].shape[0])
            for p in range(rays[n].shape[0]):
                line.GetPointIds().SetId(p, npts)
                npts += 1
            cellarray.InsertNextCell(line)
        polydata.SetLines(cellarray)
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(filename)
        writer.SetInputData(polydata)
        writer.SetDataModeToBinary()
        writer.Update()


    @staticmethod
    def builder(filename, size_t n_threads=1, str method='FSM',
                bool tt_from_rp=1, bool interp_vel=0,
                double eps=1.e-15, int maxit=20, bool weno=1,
                uint32_t nsnx=5, uint32_t nsny=5, uint32_t nsnz=5,
                uint32_t n_secondary=2, uint32_t n_tertiary=2,
                double radius_factor_tertiary=3.0,
                bool translate_grid=0):
        """
        builder(filename, n_threads=1, method='FSM', tt_from_rp=1, interp_vel=0, eps=1.e-15, maxit=20, weno=1, nsnx=5, nsny=5, nsnz=5, n_secondary=2, n_tertiary=2, radius_factor_tertiary=3.0, translate_grid=0)

        Build instance of Grid3d from VTK file

        Parameters
        ----------
        filename : str
            Name of file holding a vtkRectilinearGrid.
            The grid must have point or cell attribute named either
            'Slowness', 'slowness', 'Velocity', 'velocity', or
            'P-wave velocity'

        Other parameters are defined in Constructor

        Returns
        -------
        grid: :obj:`Grid3d`
            grid instance
        """

        reader = vtk.vtkXMLRectilinearGridReader()
        reader.SetFileName(filename)
        reader.Update()

        data = reader.GetOutput()
        x = numpy_support.vtk_to_numpy(data.GetXCoordinates())
        y = numpy_support.vtk_to_numpy(data.GetYCoordinates())
        z = numpy_support.vtk_to_numpy(data.GetZCoordinates())

        names = ('Slowness', 'slowness', 'Velocity', 'velocity',
                 'P-wave velocity')
        for name in names:
            if data.GetPointData().HasArray(name):
                cell_slowness = 0
                data = numpy_support.vtk_to_numpy(data.GetPointData().GetArray(name))
                dim = (x.size, y.size, z.size)
                break
            if data.GetCellData().HasArray(name):
                cell_slowness = 1
                data = numpy_support.vtk_to_numpy(data.GetCellData().GetArray(name))
                dim = (x.size-1, y.size-1, z.size-1)
                break
        else:
            raise ValueError('File should contain slowness or velocity data')

        if 'lowness' in name:
            slowness = data.reshape(dim, order='F').flatten()
        else:
            slowness = 1.0 / data.reshape(dim, order='F').flatten()

        g = Grid3d(x, y, z, n_threads, cell_slowness, method, tt_from_rp,
                   interp_vel, eps, maxit, weno, nsnx, nsny, nsnz,
                   n_secondary, n_tertiary, radius_factor_tertiary,
                   translate_grid)
        g.set_slowness(slowness)
        return g

    @staticmethod
    def data_kernel_straight_rays(np.ndarray[np.double_t, ndim=2] Tx,
                                  np.ndarray[np.double_t, ndim=2] Rx,
                                  np.ndarray[np.double_t, ndim=1] grx,
                                  np.ndarray[np.double_t, ndim=1] gry,
                                  np.ndarray[np.double_t, ndim=1] grz,
                                  centers=False):
        """
        data_kernel_straight_rays(Tx, Rx, grx, gry, grz, centers) -> L, (xc, yc, zc)

        Raytracing with straight rays in 3D

        Parameters
        ----------
        Tx : np.ndarray
            source coordinates, nTx by 3
                - 1st column contains X coordinates,
                - 2nd contains Y coordinates
                - 3rd contains Z coordinates
        Rx : np.ndarray
            receiver coordinates, nTx by 3
                - 1st column contains X coordinates,
                - 2nd contains Y coordinates
                - 3rd contains Z coordinates
        grx : np.ndarray
            grid node coordinates along x
        gry : np.ndarray
            grid node coordinates along y
        grz : np.ndarray
            grid node coordinates along z
        centers : bool
            return coordinates of center of cells (False by default)

        Returns
        -------
        L : scipy csr_matrix
            data kernel matrix (tt = L*slowness)
        (xc, yc, zc) : :obj:`tuple` of `np.ndarray`
            vectors of coordinates of center of cells

        Note
        ----
        Tx and Rx should contain the same number of rows, each row corresponding
        to a source-receiver pair

        """

        cdef size_t nTx = Tx.shape[0]
        cdef size_t n_grx = grx.shape[0]
        cdef size_t n_gry = gry.shape[0]
        cdef size_t n_grz = grz.shape[0]

        cdef double small = 1.e-10

        data_p = []
        indices_p = []
        indptr_p = []

        cdef size_t k = 0
        cdef size_t ix = 0
        cdef size_t iy = 0
        cdef size_t iz = 0

        cdef double x, y, z, x1, y1, z1, x2, y2, z2, dtmp, d, l, m, n
        cdef double m_y, b_y, m_z, b_z, dlx, dly, dlz, dl, xe, ye, ze
        cdef int sy, sz
        cdef bool up_y, up_z
        cdef Py_ssize_t nt

        for nt in range(nTx):
            indptr_p.append(k)

            x1 = Tx[nt, 0]
            y1 = Tx[nt, 1]
            z1 = Tx[nt, 2]
            x2 = Rx[nt, 0]
            y2 = Rx[nt, 1]
            z2 = Rx[nt, 2]

            if x1 > x2:  # on veut x croissant
                dtmp = x1
                x1 = x2
                x2 = dtmp
                dtmp = y1
                y1 = y2
                y2 = dtmp
                dtmp = z1
                z1 = z2
                z2 = dtmp

            d = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) )
            # cosinus directeurs
            l = (x2-x1)/d
            m = (y2-y1)/d
            n = (z2-z1)/d

            sy = (m > 0) - (m < 0)  # sign of m
            sz = (n > 0) - (n < 0)

            x = x1
            y = y1
            z = z1

            for ix in range(n_grx-1):
                if x < grx[ix+1] and x >= grx[ix]:
                    break
            for iy in range(n_gry-1):
                if y < gry[iy+1] and y >= gry[iy]:
                    break
            for iz in range(n_grz-1):
                if z < grz[iz+1] and z >= grz[iz]:
                    break

            if abs(l) < small:
                if abs(m) < small:
                    # X & Y constants
                    if z1 > z2:
                        dtmp = z1
                        z1 = z2
                        z2 = dtmp
                        z = z1
                        for iz in range(n_grz-1):
                            if z < grz[iz+1] and z >= grz[iz]:
                                break
                    while z < z2:
                        iCell = (ix*(n_gry-1)+iy)*(n_grz-1) + iz

                        dlz = ( grz[iz+1] if grz[iz+1] < z2 else z2 ) - z

                        indices_p.append(iCell)
                        data_p.append(dlz)
                        k += 1

                        iz += 1
                        z = grz[iz]

                elif abs(n) < small:
                    # X & Z constants
                    if y1 > y2:
                        dtmp = y1
                        y1 = y2
                        y2 = dtmp
                        y = y1
                        for iy in range(n_gry-1):
                            if y < gry[iy+1] and y >= gry[iy]:
                                break
                    while y < y2:
                        iCell = (ix*(n_gry-1)+iy)*(n_grz-1) + iz

                        dly = ( gry[iy+1] if gry[iy+1] < y2 else y2 ) - y

                        indices_p.append(iCell)
                        data_p.append(dly)
                        k += 1

                        iy += 1
                        y = gry[iy]

                else:
                    # seul X constant
                    if y1 > y2:
                        dtmp = y1
                        y1 = y2
                        y2 = dtmp
                        dtmp = z1
                        z1 = z2
                        z2 = dtmp
                        y = y1
                        for iy in range(n_gry-1):
                            if y < gry[iy+1] and y >= gry[iy]:
                                break
                        z = z1
                        for iz in range(n_grz-1):
                            if z < grz[iz+1] and z >= grz[iz]:
                                break

                    m_z = (z2-z1) / (y2-y1)
                    b_z = z2 - m_z*y2
                    up_z = m_z > 0

                    while y < y2:

                        zi = m_z*gry[iy+1] + b_z
                        if up_z:
                            while z < zi and z < z2:
                                iCell = (ix*(n_gry-1)+iy)*(n_grz-1) + iz

                                ze = grz[iz+1] if grz[iz+1]<zi else zi
                                ze = ze if ze < z2 else z2
                                ye = (ze-b_z) / m_z
                                dly = ye - y
                                dlz = ze - z
                                dl = sqrt( dly*dly + dlz*dlz )

                                indices_p.append(iCell)
                                data_p.append(dl)
                                k += 1

                                y = ye
                                z = ze

                                if abs(z-grz[iz+1]) < small:
                                    iz += 1

                        else: # down
                            while z > zi and z > z2:
                                iCell = (ix*(n_gry-1)+iy)*(n_grz-1) + iz

                                ze = grz[iz] if grz[iz] > zi else zi
                                ze = ze if ze > z2 else z2
                                ye = (ze-b_z)/m_z
                                dly = ye - y
                                dlz = ze - z
                                dl = sqrt( dly*dly + dlz*dlz )

                                indices_p.append(iCell)
                                data_p.append(dl)
                                k += 1

                                y = ye
                                z = ze
                                if abs(z-grz[iz]) < small:
                                    iz -= 1

                        iy += 1
                        y = gry[iy]
            else:
                if abs(m) < small and abs(n) < small:
                    # Y & Z constants
                    while x < x2:
                        iCell = (ix*(n_gry-1)+iy)*(n_grz-1) + iz

                        dlx = ( grx[ix+1] if grx[ix+1] < x2 else x2 ) - x

                        indices_p.append(iCell)
                        data_p.append(dlx)
                        k += 1

                        ix += 1
                        x = grx[ix]

                elif abs(m) < small:
                    # seul Y constant
                    m_z = (z2-z1) / (x2-x1)
                    b_z = z2 - m_z*x2
                    up_z = m_z>0

                    while x < x2:

                        zi = m_z*grx[ix+1] + b_z

                        if up_z:
                            while z < zi and z < z2:
                                iCell = (ix*(n_gry-1)+iy)*(n_grz-1) + iz

                                ze = grz[iz+1] if grz[iz+1] < zi else zi
                                ze = ze if ze < z2 else z2
                                xe = (ze-b_z) / m_z
                                dlx = xe - x
                                dlz = ze - z
                                dl = sqrt( dlx*dlx + dlz*dlz )

                                indices_p.append(iCell)
                                data_p.append(dl)
                                k += 1

                                x = xe
                                z = ze
                                if abs(z-grz[iz+1]) < small:
                                    iz += 1
                        else: # down
                            while z > zi and z > z2:
                                iCell = (ix*(n_gry-1)+iy)*(n_grz-1) + iz

                                ze = grz[iz] if grz[iz] > zi else zi
                                ze = ze if ze > z2 else z2
                                xe = (ze-b_z) / m_z
                                dlx = xe - x
                                dlz = ze - z
                                dl = sqrt( dlx*dlx + dlz*dlz )

                                indices_p.append(iCell)
                                data_p.append(dl)
                                k += 1

                                x = xe
                                z = ze
                                if abs(z-grz[iz]) < small:
                                    iz -= 1

                        ix += 1
                        x = grx[ix]
                elif abs(n) < small:
                    # seul Z constant
                    m_y = (y2-y1) / (x2-x1)
                    b_y = y2 - m_y*x2
                    up_y = m_y > 0

                    while x < x2:

                        yi = m_y*grx[ix+1] + b_y

                        if up_y:
                            while y < yi and y < y2:
                                iCell = (ix*(n_gry-1)+iy)*(n_grz-1) + iz

                                ye = gry[iy+1] if gry[iy+1] < yi else yi
                                ye = ye if ye < y2 else y2
                                xe = (ye-b_y) / m_y
                                dlx = xe - x
                                dly = ye - y
                                dl = sqrt( dlx*dlx + dly*dly )

                                indices_p.append(iCell)
                                data_p.append(dl)
                                k += 1

                                x = xe
                                y = ye
                                if abs(y-gry[iy+1]) < small:
                                    iy += 1
                        else:  # down
                            while y > yi and y > y2:
                                iCell = (ix*(n_gry-1)+iy)*(n_grz-1) + iz

                                ye = gry[iy] if gry[iy] > yi else yi
                                ye = ye if ye > y2 else y2
                                xe = (ye-b_y)/m_y
                                dlx = xe - x
                                dly = ye - y
                                dl = sqrt( dlx*dlx + dly*dly )

                                indices_p.append(iCell)
                                data_p.append(dl)
                                k += 1

                                x = xe
                                y = ye
                                if abs(y-gry[iy]) < small:
                                    iy -= 1

                        ix += 1
                        x = grx[ix]
                else:
                    while x < x2:

                        m_y = (y2-y1) / (x2-x1)
                        b_y = y2 - m_y*x2
                        up_y = m_y > 0

                        m_z = (z2-z1) / (x2-x1)
                        b_z = z2 - m_z*x2
                        up_z = m_z > 0

                        yi = m_y*grx[ix+1] + b_y
                        zi = m_z*grx[ix+1] + b_z

                        while (sy*(yi-y)) > 0 and (sy*(y2-y)) > 0 and \
							  (sz*(zi-z)) > 0 and (sz*(z2-z)) > 0:

                            if up_y:
                                ye = gry[iy+1] if gry[iy+1] < yi else yi
                                ye = ye if ye<y2 else y2
                            else:
                                ye = gry[iy] if gry[iy] > yi else yi
                                ye = ye if ye>y2 else y2
                            if up_z:
                                ze = grz[iz+1] if grz[iz+1] < zi else zi
                                ze = ze if ze < z2 else z2
                            else:
                                ze = grz[iz] if grz[iz] > zi else zi
                                ze = ze if ze > z2 else z2

                            if (ze-b_z)/m_z < (ye-b_y)/m_y:  # we cross z before y
                                iCell = (ix*(n_gry-1)+iy)*(n_grz-1) + iz

                                xe = (ze-b_z) / m_z
                                ye = m_y*xe + b_y
                                dlx = xe - x
                                dly = ye - y
                                dlz = ze - z
                                dl = sqrt( dlx*dlx + dly*dly + dlz*dlz )

                                indices_p.append(iCell)
                                data_p.append(dl)
                                k += 1

                                x = xe
                                y = ye
                                z = ze
                                if up_z:
                                    if abs(z-grz[iz+1]) < small:
                                        iz += 1
                                else:
                                    if abs(z-grz[iz]) < small:
                                        iz -= 1

                            else:  # we cross y before z
                                iCell = (ix*(n_gry-1)+iy)*(n_grz-1) + iz

                                xe = (ye-b_y)/m_y
                                ze = m_z*xe + b_z
                                dlx = xe - x
                                dly = ye - y
                                dlz = ze - z

                                dl = sqrt( dlx*dlx + dly*dly + dlz*dlz )

                                indices_p.append(iCell)
                                data_p.append(dl)
                                k += 1

                                x = xe
                                y = ye
                                z = ze
                                if up_y:
                                    if abs(y-gry[iy+1]) < small:
                                        iy += 1
                                else:
                                    if abs(y-gry[iy]) < small:
                                        iy -= 1

                        ix += 1
                        x = grx[ix]

        indptr_p.append(k)
        L = sp.csr_matrix((data_p, indices_p, indptr_p),
                          shape=(nTx, (n_grx-1)*(n_gry-1)*(n_grz-1)))

        if centers:
            xc = (grx[1:]+grx[:-1])/2
            yc = (gry[1:]+gry[:-1])/2
            zc = (grz[1:]+grz[:-1])/2
            return L, (xc, yc, zc)
        else:
            return L


cdef class Grid2d:
    """
    class to perform raytracing with 2D rectilinear grids

    Attributes
    ----------
    x: np.ndarray
        node coordinates along x
    z: np.ndarray
        node coordinates along z
    dx: float
        node separation along x
    dz: float
        node separation along z
    shape: (int, int)
        number of parameters along each dimension
    nparams: int
        total number of parameters for grid
    n_threads: int
        number of threads for raytracing

    Constructor:

    Grid2d(x, z, n_threads=1, cell_slowness=1, method='SPM', aniso='iso', eps=1.e-15, maxit=20, weno=1, rotated_template=0, nsnx=10, nsnz=10, n_secondary=3, n_tertiary=3, radius_factor_tertiary=3.0, tt_from_rp=0) -> Grid2d

    Parameters
    ----------
    x : np.ndarray
        node coordinates along x
    z : np.ndarray
        node coordinates along z
    n_threads : int
        number of threads for raytracing (default is 1)
    cell_slowness : bool
        slowness defined for cells (True) or nodes (False) (default is 1)
    method : string
        raytracing method (default is SPM)
            - 'FSM' : fast marching method
            - 'SPM' : shortest path method
            - 'DSPM' : dynamic shortest path method
    aniso : string
        type of anisotropy (implemented only for the SPM method)
            - 'iso' : isotropic medium
            - 'elliptical' : elliptical anisotropy
            - 'tilted_elliptical' : tilted elliptical anisotropy
            - 'vti_psv' : vertical transverse isotropy, P and SV waves
            - 'vti_sh' : vertical transverse isotropy, SH waves
    eps : double
        convergence criterion (FSM) (default is 1e-15)
    maxit : int
        max number of sweeping iterations (FSM) (default is 20)
    weno : bool
        use 3rd order weighted essentially non-oscillatory operator (FSM)
        (default is True)
    rotated_template : bool
        use rotated templates (FSM)
    nsnx : int
        number of secondary nodes in x (SPM) (default is 10)
    nsnz : int
        number of secondary nodes in z (SPM) (default is 10)
    n_secondary : int
        number of secondary nodes (DSPM) (default is 3)
    n_tertiary : int
        number of tertiary nodes (DSPM) (default is 3)
    radius_factor_tertiary : double
            multiplication factor used to compute radius of sphere around source
            that includes tertiary nodes (DSPM).  The radius is the average edge
            length multiplied by this factor (default is 3)
    tt_from_rp : bool
        compute traveltime using raypaths (available for FSM and DSPM only)
        (default is False)
    """
    cdef vector[double] _x
    cdef vector[double] _z
    cdef double _dx
    cdef double _dz
    cdef bool cell_slowness
    cdef size_t _n_threads
    cdef char method
    cdef char iso
    cdef double eps
    cdef int maxit
    cdef bool weno
    cdef bool rotated_template
    cdef bool tt_from_rp
    cdef uint32_t nsnx
    cdef uint32_t nsnz
    cdef uint32_t n_secondary
    cdef uint32_t n_tertiary
    cdef double radius_factor_tertiary

    cdef Grid2D[double,uint32_t,sxz[double]]* grid

    def __cinit__(self, np.ndarray[np.double_t, ndim=1] x,
                  np.ndarray[np.double_t, ndim=1] z, size_t n_threads=1,
                  bool cell_slowness=1, str method='SPM', str aniso='iso',
                  double eps=1.e-15, int maxit=20, bool weno=1,
                  bool rotated_template=0, uint32_t nsnx=10, uint32_t nsnz=10,
                  uint32_t n_secondary=3, uint32_t n_tertiary=3,
                  double radius_factor_tertiary=3.0, bool tt_from_rp=0):

        cdef uint32_t nx = x.size-1
        cdef uint32_t nz = z.size-1
        self._dx = x[1] - x[0]
        self._dz = z[1] - z[0]
        cdef double xmin = x[0]
        cdef double zmin = z[0]
        self.cell_slowness = cell_slowness
        self._n_threads = n_threads
        self.eps = eps
        self.maxit = maxit
        self.weno = weno
        self.rotated_template = rotated_template
        self.nsnx = nsnx
        self.nsnz = nsnz
        self.iso = b'i'
        self.n_secondary = n_secondary
        self.n_tertiary = n_tertiary
        self.radius_factor_tertiary = radius_factor_tertiary
        self.tt_from_rp = tt_from_rp

        cdef use_edge_length = True

        for val in x:
            self._x.push_back(val)
        for val in z:
            self._z.push_back(val)
        self._x.shrink_to_fit()
        self._z.shrink_to_fit()

        if cell_slowness:
            if method == 'SPM':
                self.method = b's'
                if aniso == 'iso':
                    self.grid = new Grid2Drcsp[double,uint32_t,sxz[double],cell2d](
                                    nx, nz, self._dx, self._dz,
                                    xmin, zmin, nsnx, nsnz, tt_from_rp, n_threads)
                elif aniso == 'elliptical':
                    self.iso = b'e'
                    self.grid = new Grid2Drcsp[double,uint32_t,sxz[double],cell2d_e](
                                    nx, nz, self._dx, self._dz,
                                    xmin, zmin, nsnx, nsnz, tt_from_rp, n_threads)
                elif aniso == 'tilted_elliptical':
                    self.iso = b't'
                    self.grid = new Grid2Drcsp[double,uint32_t,sxz[double],cell2d_te](
                                    nx, nz, self._dx, self._dz,
                                    xmin, zmin, nsnx, nsnz, tt_from_rp, n_threads)
                elif aniso == 'vti_psv':
                    self.iso = b'p'
                    self.grid = new Grid2Drcsp[double,uint32_t,sxz[double],cell2d_p](
                                    nx, nz, self._dx, self._dz,
                                    xmin, zmin, nsnx, nsnz, tt_from_rp, n_threads)
                elif aniso == 'vti_sh':
                    self.iso = b'h'
                    self.grid = new Grid2Drcsp[double,uint32_t,sxz[double],cell2d_h](
                                    nx, nz, self._dx, self._dz,
                                    xmin, zmin, nsnx, nsnz, tt_from_rp, n_threads)
                else:
                    raise ValueError('Anisotropy model not implemented')
            elif method == 'FSM':
                self.method = b'f'
                self.grid = new Grid2Drcfs[double,uint32_t,sxz[double]](nx, nz,
                                self._dx, self._dz, xmin, zmin, eps,
                                maxit, weno, rotated_template, tt_from_rp, n_threads)
            elif method == 'DSPM':
                self.method = b'd'
                self.grid = new Grid2Drcdsp[double,uint32_t,sxz[double],cell2d](
                                    nx, nz, self._dx, self._dz,
                                    xmin, zmin, n_secondary, n_tertiary,
                                    radius_factor_tertiary, tt_from_rp,
                                    use_edge_length, n_threads)
            else:
                raise ValueError('Method {0:s} undefined'.format(method))
        else:
            if method == 'SPM':
                self.method = b's'
                self.grid = new Grid2Drnsp[double,uint32_t,sxz[double]](nx, nz,
                                self._dx, self._dz, xmin, zmin,
                                nsnx, nsnz, tt_from_rp, n_threads)
            elif method == 'FSM':
                self.method = b'f'
                self.grid = new Grid2Drnfs[double,uint32_t,sxz[double]](nx, nz,
                                self._dx, self._dz, xmin, zmin, eps,
                                maxit, weno, rotated_template, tt_from_rp, n_threads)
            elif method == 'DSPM':
                self.method = b'd'
                self.grid = new Grid2Drndsp[double,uint32_t,sxz[double]](nx, nz,
                                self._dx, self._dz, xmin, zmin,
                                n_secondary, n_tertiary,
                                radius_factor_tertiary, tt_from_rp,
                                use_edge_length, n_threads)
            else:
                raise ValueError('Method {0:s} undefined'.format(method))

    def __dealloc__(self):
        del self.grid

    def __reduce__(self):
        if self.method == b'f':
            method = 'FSM'
        elif self.method == b's':
            method = 'SPM'
        elif self.method == b'd':
            method = 'DSPM'

        if self.iso == b'i':
            aniso = 'iso'
        elif self.iso == b'e':
            aniso = 'elliptical'
        elif self.iso == b't':
            aniso = 'tilted_elliptical'
        elif self.iso == b'p':
            aniso = 'vti_psv'
        elif self.iso == b'h':
            aniso = 'vti_sh'

        constructor_params = (self.n_threads, self.cell_slowness, method,
                              aniso, self.eps, self.maxit, self.weno,
                              self.rotated_template, self.nsnx, self.nsnz,
                              self.n_secondary, self.n_tertiary,
                              self.radius_factor_tertiary, self.tt_from_rp)
        return (_rebuild2d, (self.x, self.z, constructor_params))

    @property
    def x(self):
        """np.ndarray: node coordinates along x"""
        tmp = np.empty((self._x.size(),))
        cdef int n
        for n in range(self._x.size()):
            tmp[n] = self._x[n]
        return tmp

    @property
    def z(self):
        """np.ndarray: node coordinates along z"""
        tmp = np.empty((self._z.size(),))
        cdef int n
        for n in range(self._z.size()):
            tmp[n] = self._z[n]
        return tmp

    @property
    def dx(self):
        """float: node separation along x"""
        return self._dx

    @property
    def dz(self):
        """float: node separation along x"""
        return self._dz

    @property
    def shape(self):
        """:obj:`list` of :obj:`int`: number of parameters along each dimension"""
        if self.cell_slowness:
            return (self._x.size()-1, self._z.size()-1)
        else:
            return (self._x.size(), self._z.size())

    @property
    def n_threads(self):
        """int: number of threads for raytracing"""
        return self._n_threads

    @property
    def nparams(self):
        """int: total number of parameters for grid"""
        if self.cell_slowness:
            return (self._x.size()-1) * (self._z.size()-1)
        else:
            return self._x.size() * self._z.size()

    def set_use_thread_pool(self, use_thread_pool):
        """
        set_use_thread_pool(use_thread_pool)

        Set option to use thread pool instead of parallel loop

        Parameters
        ----------
        use_thread_pool : bool
            option value
        """
        self.grid.setUsePool(use_thread_pool)

    def set_traveltime_from_raypath(self, traveltime_from_raypath):
        """
        set_traveltime_from_raypath(ttrp)

        Set option to compute traveltime using raypath

        Parameters
        ----------
        ttrp : bool
            option value
        """
        self.grid.setTraveltimeFromRaypath(traveltime_from_raypath)

    def get_number_of_nodes(self):
        """
        Returns
        -------
        int:
            number of nodes in grid
        """
        return self._x.size() * self._z.size()

    def get_number_of_cells(self):
        """
        Returns
        -------
        int:
            number of cells in grid
        """
        return (self._x.size()-1) * (self._z.size()-1)

    def get_grid_traveltimes(self, thread_no=0):
        """
        get_grid_traveltimes(thread_no=0)

        Obtain traveltimes computed at primary grid nodes

        Parameters
        ----------
        thread_no : int
            thread used to computed traveltimes (default is 0)

        Returns
        -------
        tt: np ndarray, shape (nx, nz)
        """
        if thread_no >= self._n_threads:
            raise ValueError('Thread number is larger than number of threads')
        cdef vector[double] tmp
        cdef int n
        self.grid.getTT(tmp, thread_no)
        tt = np.empty((tmp.size(),))
        for n in range(tmp.size()):
            tt[n] = tmp[n]
        shape = (self._x.size(), self._z.size())
        return tt.reshape(shape)

    def is_outside(self, np.ndarray[np.double_t, ndim=2] pts):
        """
        is_outside(pts)

        Check if points are outside grid

        Parameters
        ----------
        pts : np ndarray, shape(npts, 3)
            coordinates of points to check

        Returns
        -------
        bool:
            True if at least one point outside grid
        """
        return ( np.min(pts[:,0]) < self._x.front() or np.max(pts[:,0]) > self._x.back() or
                np.min(pts[:,1]) < self._z.front() or np.max(pts[:,1]) > self._z.back() )

    def get_slowness(self):
        """Returns slowness of grid

        Returns
        -------
        slowness : np ndarray, shape (nx, nz)
        """
        cdef int i
        cdef vector[double] slown
        self.grid.getSlowness(slown)
        nx, nz = self.shape
        slowness = np.ndarray((nx*nz,))
        for i in range(slown.size()):
            slowness[i] = slown[i]
        return slowness.reshape((nx, nz))

    def set_slowness(self, slowness):
        """
        set_slowness(slowness)

        Assign slowness to grid

        Parameters
        ----------
        slowness : np ndarray, shape (nx, nz)
            slowness may also have been flattened (with default 'C' order)
        """
        if self.cell_slowness:
            nx = self._x.size()-1
            nz = self._z.size()-1
        else:
            nx = self._x.size()
            nz = self._z.size()
        if slowness.size != nx*nz:
            raise ValueError('Slowness vector has wrong size')

        cdef vector[double] data
        cdef int i
        if slowness.ndim == 2:
            if slowness.shape != (nx, nz):
                raise ValueError('Slowness has wrong shape')
            tmp = slowness.flatten()
            for i in range(nx*nz):
                data.push_back(tmp[i])
        elif slowness.ndim == 1:
            for i in range(nx*nz):
                data.push_back(slowness[i])
        else:
            raise ValueError('Slowness must be 1D or 3D ndarray')
        self.grid.setSlowness(data)

    def set_velocity(self, velocity):
        """
        set_velocity(velocity)

        Assign velocity to grid

        Parameters
        ----------
        velocity : np ndarray, shape (nx, nz)
            velocity may also have been flattened (with default 'C' order)
        """
        if self.cell_slowness:
            nx = self._x.size()-1
            nz = self._z.size()-1
        else:
            nx = self._x.size()
            nz = self._z.size()
        if velocity.size != nx*nz:
            raise ValueError('velocity vector has wrong size')

        cdef vector[double] data
        cdef int i
        if velocity.ndim == 2:
            if velocity.shape != (nx, nz):
                raise ValueError('velocity has wrong shape')
            tmp = velocity.flatten()
            for i in range(nx*nz):
                data.push_back(1./tmp[i])
        elif velocity.ndim == 1:
            for i in range(nx*nz):
                data.push_back(1./velocity[i])
        else:
            raise ValueError('velocity must be 1D or 3D ndarray')
        self.grid.setSlowness(data)

    def set_xi(self, xi):
        """
        set_xi(xi)

        Assign elliptical anisotropy ratio to grid

        Parameters
        ----------
        xi : np ndarray, shape (nx, nz)
            xi may also have been flattened (with default 'C' order)
        """
        if self.cell_slowness:
            nx = self._x.size()-1
            nz = self._z.size()-1
        else:
            nx = self._x.size()
            nz = self._z.size()
        if xi.size != nx*nz:
            raise ValueError('xi vector has wrong size')

        cdef vector[double] data
        cdef int i
        if xi.ndim == 2:
            if xi.shape != (nx, nz):
                raise ValueError('xi has wrong shape')
            tmp = xi.flatten()
            for i in range(nx*nz):
                data.push_back(tmp[i])
        elif xi.ndim == 1:
            tmp = xi.reshape((nx, nz)).flatten()
            for i in range(nx*nz):
                data.push_back(tmp[i])
        else:
            raise ValueError('xi must be 1D or 3D ndarray')
        self.grid.setXi(data)

    def set_tilt_angle(self, theta):
        """
        set_tilt_angle(theta)

        Assign tilted elliptical anisotropy angle to grid

        Parameters
        ----------
        theta : np ndarray, shape (nx, nz)
            theta may also have been flattened (with default 'C' order)
        """
        if self.cell_slowness:
            nx = self._x.size()-1
            nz = self._z.size()-1
        else:
            nx = self._x.size()
            nz = self._z.size()
        if theta.size != nx*nz:
            raise ValueError('theta vector has wrong size')

        cdef vector[double] data
        cdef int i
        if theta.ndim == 2:
            if theta.shape != (nx, nz):
                raise ValueError('theta has wrong shape')
            tmp = theta.flatten()
            for i in range(nx*nz):
                data.push_back(tmp[i])
        elif theta.ndim == 1:
            tmp = theta.reshape((nx, nz)).flatten()
            for i in range(nx*nz):
                data.push_back(tmp[i])
        else:
            raise ValueError('theta must be 1D or 3D ndarray')
        self.grid.setTiltAngle(data)

    def set_Vp0(self, v):
        """
        set_Vp0(v)

        Assign vertical Vp to grid (VTI medium)

        Parameters
        ----------
        v : np ndarray, shape (nx, nz)
            v may also have been flattened (with default 'C' order)
        """
        if self.cell_slowness:
            nx = self._x.size()-1
            nz = self._z.size()-1
        else:
            nx = self._x.size()
            nz = self._z.size()
        if v.size != nx*nz:
            raise ValueError('v vector has wrong size')

        cdef vector[double] data
        cdef int i
        if v.ndim == 2:
            if v.shape != (nx, nz):
                raise ValueError('v has wrong shape')
            tmp = v.flatten()
            for i in range(nx*nz):
                data.push_back(tmp[i])
        elif v.ndim == 1:
            tmp = v.reshape((nx, nz)).flatten()
            for i in range(nx*nz):
                data.push_back(tmp[i])
        else:
            raise ValueError('v must be 1D or 3D ndarray')
        self.grid.setVp0(data)

    def set_Vs0(self, v):
        """
        set_Vs0(v)

        Assign vertical Vs to grid (VTI medium)

        Parameters
        ----------
        v : np ndarray, shape (nx, nz)
            v may also have been flattened (with default 'C' order)
        """
        if self.cell_slowness:
            nx = self._x.size()-1
            nz = self._z.size()-1
        else:
            nx = self._x.size()
            nz = self._z.size()
        if v.size != nx*nz:
            raise ValueError('v vector has wrong size')

        cdef vector[double] data
        cdef int i
        if v.ndim == 2:
            if v.shape != (nx, nz):
                raise ValueError('v has wrong shape')
            tmp = v.flatten()
            for i in range(nx*nz):
                data.push_back(tmp[i])
        elif v.ndim == 1:
            tmp = v.reshape((nx, nz)).flatten()
            for i in range(nx*nz):
                data.push_back(tmp[i])
        else:
            raise ValueError('v must be 1D or 3D ndarray')
        self.grid.setVs0(data)

    def set_delta(self, v):
        """
        set_delta(d)

        Assign Thomsen delta parameter to grid (VTI medium, P-SV waves)

        Parameters
        ----------
        d : np ndarray, shape (nx, nz)
            d may also have been flattened (with default 'C' order)
        """
        if self.cell_slowness:
            nx = self._x.size()-1
            nz = self._z.size()-1
        else:
            nx = self._x.size()
            nz = self._z.size()
        if v.size != nx*nz:
            raise ValueError('v vector has wrong size')

        cdef vector[double] data
        cdef int i
        if v.ndim == 2:
            if v.shape != (nx, nz):
                raise ValueError('v has wrong shape')
            tmp = v.flatten()
            for i in range(nx*nz):
                data.push_back(tmp[i])
        elif v.ndim == 1:
            tmp = v.reshape((nx, nz)).flatten()
            for i in range(nx*nz):
                data.push_back(tmp[i])
        else:
            raise ValueError('v must be 1D or 3D ndarray')
        self.grid.setDelta(data)

    def set_epsilon(self, v):
        """
        set_epsilon(e)

        Assign Thomsen epsilon parameter to grid (VTI medium, P-SV waves)

        Parameters
        ----------
        e : np ndarray, shape (nx, nz)
            e may also have been flattened (with default 'C' order)
        """
        if self.cell_slowness:
            nx = self._x.size()-1
            nz = self._z.size()-1
        else:
            nx = self._x.size()
            nz = self._z.size()
        if v.size != nx*nz:
            raise ValueError('v vector has wrong size')

        cdef vector[double] data
        cdef int i
        if v.ndim == 2:
            if v.shape != (nx, nz):
                raise ValueError('v has wrong shape')
            tmp = v.flatten()
            for i in range(nx*nz):
                data.push_back(tmp[i])
        elif v.ndim == 1:
            tmp = v.reshape((nx, nz)).flatten()
            for i in range(nx*nz):
                data.push_back(tmp[i])
        else:
            raise ValueError('v must be 1D or 3D ndarray')
        self.grid.setEpsilon(data)

    def set_gamma(self, v):
        """
        set_gamma(g)

        Assign Thomsen gamma parameter to grid (VTI medium, SH waves)

        Parameters
        ----------
        g : np ndarray, shape (nx, nz)
            g may also have been flattened (with default 'C' order)
        """
        if self.cell_slowness:
            nx = self._x.size()-1
            nz = self._z.size()-1
        else:
            nx = self._x.size()
            nz = self._z.size()
        if v.size != nx*nz:
            raise ValueError('v vector has wrong size')

        cdef vector[double] data
        cdef int i
        if v.ndim == 2:
            if v.shape != (nx, nz):
                raise ValueError('v has wrong shape')
            tmp = v.flatten()
            for i in range(nx*nz):
                data.push_back(tmp[i])
        elif v.ndim == 1:
            tmp = v.reshape((nx, nz)).flatten()
            for i in range(nx*nz):
                data.push_back(tmp[i])
        else:
            raise ValueError('v must be 1D or 3D ndarray')
        self.grid.setGamma(data)

    def compute_K(self, order=1):
        """
        Compute smoothing matrices

        Parameters
        ----------
        order : int
            order of smoothing operator, accept 1 or 2 (1 by default)

        Returns
        -------
        Kx, Kz : :obj:`tuple` of :obj:`csr_matrix`
            matrices for derivatives along x & z
        """

        nx, nz = self.shape
        if order == 1:
            # forward operator is (u_{i+1} - u_i)/dx
            # centered operator is (u_{i+1} - u_{i-1})/(2dx)
            # backward operator is (u_i - u_{i-1})/dx

            idx = 1 / self.dx
            idz = 1 / self.dz

            i = np.kron(np.arange(nx * nz), np.ones((2, ), dtype=np.int32))
            j = np.zeros((nz * nx * 2, ), dtype=np.int32)
            v = np.zeros((nz * nx * 2, ))

            jj = np.vstack((np.arange(nz), nz + np.arange(nz))).T
            jj = jj.flatten()
            j[:2 * nz] = jj
            vd = idx * np.tile(np.array([-1, 1]), (nz, ))
            v[:2 * nz] = vd

            jj = np.vstack((-nz + np.arange(nz), nz + np.arange(nz))).T
            jj = jj.flatten()
            for n in range(1, nx - 1):
                j[n * 2 * nz:(n + 1) * 2 * nz] = n * nz + jj
                v[n * 2 * nz:(n + 1) * 2 * nz] = 0.5 * vd

            jj = np.vstack((-nz + np.arange(nz), np.arange(nz))).T
            jj = jj.flatten()
            j[(nx - 1) * 2 * nz:nx * 2 * nz] = (nx - 1) * nz + jj
            v[(nx - 1) * 2 * nz:nx * 2 * nz] = vd

            Kx = sp.csr_matrix((v, (i, j)))

            jj = np.vstack((np.hstack((0, np.arange(nz - 1))),
                            np.hstack((np.arange(1, nz), nz - 1)))).T
            jj = jj.flatten()
            vd = idz * np.hstack((np.array([-1, 1]),
                                  np.tile(np.array([-0.5, 0.5]), (nz - 2,)),
                                  np.array([-1, 1])))

            for n in range(nx):
                j[n * 2 * nz:(n + 1) * 2 * nz] = n * nz + jj
                v[n * 2 * nz:(n + 1) * 2 * nz] = vd

            Kz = sp.csr_matrix((v, (i, j)))
        elif order == 2:
            # forward operator is (u_i - 2u_{i+1} + u_{i+2})/dx^2
            # centered operator is (u_{i-1} - 2u_i + u_{i+1})/dx^2
            # backward operator is (u_{i-2} - 2u_{i-1} + u_i)/dx^2

            idx2 = 1 / (self.dx * self.dx)
            idz2 = 1 / (self.dz * self.dz)

            i = np.kron(np.arange(nx * nz), np.ones((3, ), dtype=np.int32))
            j = np.zeros((nz * nx * 3, ), dtype=np.int32)
            v = np.zeros((nz * nx * 3, ))

            jj = np.vstack((np.arange(nz), nz + np.arange(nz),
                            2 * nz + np.arange(nz))).T
            jj = jj.flatten()
            j[:3 * nz] = jj
            vd = idx2 * np.tile(np.array([1.0, -2.0, 1.0]), (nz, ))
            v[:3 * nz] = vd

            for n in range(1, nx - 1):
                j[n * 3 * nz:(n + 1) * 3 * nz] = (n - 1) * nz + jj
                v[n * 3 * nz:(n + 1) * 3 * nz] = vd

            j[(nx - 1) * 3 * nz:nx * 3 * nz] = (nx - 3) * nz + jj
            v[(nx - 1) * 3 * nz:nx * 3 * nz] = vd

            Kx = sp.csr_matrix((v, (i, j)))

            jj = np.vstack((np.hstack((0, np.arange(nz - 2), nz - 3)),
                            np.hstack((1, np.arange(1, nz - 1), nz - 2)),
                            np.hstack((2, np.arange(2, nz), nz - 1)))).T
            jj = jj.flatten()
            vd = vd * idz2 / idx2

            for n in range(nx):
                j[n * 3 * nz:(n + 1) * 3 * nz] = n * nz + jj
                v[n * 3 * nz:(n + 1) * 3 * nz] = vd

            Kz = sp.csr_matrix((v, (i, j)))

        else:
            raise ValueError('order value not valid (1 or 2 accepted)')

        return Kx, Kz

    def raytrace(self, source, rcv, slowness=None, xi=None, theta=None,
                 Vp0=None, Vs0=None, delta=None, epsilon=None, gamma=None,
                 thread_no=None, aggregate_src=False, compute_L=False,
                 return_rays=False):
        """
        raytrace(source, rcv, slowness=None, xi=None, theta=None, Vp0=None, Vs0=None, delta=None, epsilon=None, gamma=None, thread_no=None, aggregate_src=False, compute_L=False, return_rays=False) -> tt, rays, L

        Perform raytracing

        Parameters
        ----------
        source : 2D np.ndarray with 2 or 3 columns
            see notes below
        rcv : 2D np.ndarray with 2 columns
            Columns correspond to x, y and z coordinates
        slowness : np ndarray, shape (nx, nz) (None by default)
            slowness at grid nodes or cells (depending on cell_slowness)
            slowness may also have been flattened (with default 'C' order)
            if None, slowness must have been assigned previously
        xi : np ndarray, shape (nx, nz) (None by default)
            xi at grid cells (only for SPM & cell_slowness=True)
            xi may also have been flattened (with default 'C' order)
            if None, xi must have been assigned previously
        theta : np ndarray, shape (nx, nz) (None by default)
            theta at grid cells (only for SPM & cell_slowness=True)
            theta may also have been flattened (with default 'C' order)
            if None, theta must have been assigned previously
        Vp0 : np ndarray, shape (nx, nz) (None by default)
            Vp0 at grid cells (only for SPM & cell_slowness=True)
            Vp0 may also have been flattened (with default 'C' order)
            if None, Vp0 must have been assigned previously
        Vs0 : np ndarray, shape (nx, nz) (None by default)
            Vs0 at grid cells (only for SPM & cell_slowness=True)
            Vs0 may also have been flattened (with default 'C' order)
            if None, Vs0 must have been assigned previously
        delta : np ndarray, shape (nx, nz) (None by default)
            delta at grid cells (only for SPM & cell_slowness=True)
            delta may also have been flattened (with default 'C' order)
            if None, delta must have been assigned previously
        epsilon : np ndarray, shape (nx, nz) (None by default)
            epsilon at grid cells (only for SPM & cell_slowness=True)
            epsilon may also have been flattened (with default 'C' order)
            if None, epsilon must have been assigned previously
        gamma : np ndarray, shape (nx, nz) (None by default)
            gamma at grid cells (only for SPM & cell_slowness=True)
            gamma may also have been flattened (with default 'C' order)
            if None, gamma must have been assigned previously
        thread_no : int (None by default)
            Perform calculations in thread number "thread_no"
            if None, attempt to run in parallel if warranted by number of
            sources and value of n_threads in constructor
        aggregate_src : bool (False by default)
            if True, all source coordinates belong to a single event
        compute_L : bool (False by default)
            Compute matrices of partial derivative of travel time w/r to slowness
        return_rays : bool (False by default)
            Return raypaths

        Returns
        -------
        tt : np.ndarray
            travel times for the appropriate source-rcv  (see Notes below)
        rays : :obj:`list` of :obj:`np.ndarray`
            Coordinates of segments forming raypaths (if return_rays is True)
        L : scipy csr_matrix
            Matrix of partial derivative of travel time w/r to slowness

        Notes
        -----
        If source has 2 columns:
            - Columns correspond to x and z coordinates
            - Origin time (t0) is 0 for all points
        If source has 3 columns:
            - 1st column corresponds to origin times
            - 2nd & 3rd columns correspond to x and z coordinates

        source and rcv can contain the same number of rows, each row
        corresponding to a source-receiver pair, or the number of rows may
        differ if aggregate_src is True or if all rows in source are identical.
        """
        # check input data consistency

        if source.ndim != 2 or rcv.ndim != 2:
            raise ValueError('source and rcv should be 2D arrays')

        if self.iso != b'i' and self.method != b's' and not self.cell_slowness:
            raise NotImplementedError('Anisotropic raytracing implemented only for SPM')

        if compute_L and not self.cell_slowness:
            raise NotImplementedError('compute_L defined only for grids with slowness defined for cells')

        if compute_L and (self.iso == b'p' or self.iso == b'h'):
            raise NotImplementedError('compute_L not implemented for VTI media')

        if source.shape[1] == 2:
            src = source
            _, ind = np.unique(source, axis=0, return_index=True)
            Tx = source[np.sort(ind), :]     # this to keep the original order
            t0 = np.zeros((Tx.shape[0], 1))
            nTx = Tx.shape[0]
        elif source.shape[1] == 3:
            src = source[:,1:3]
            _, ind = np.unique(source, axis=0, return_index=True)
            tmp = source[np.sort(ind), :]    # this to keep the original order
            nTx = tmp.shape[0]
            Tx = tmp[:,1:3]
            t0 = tmp[:,0]
        else:
            raise ValueError('source should be either nsrc x 2 or 3')

        if src.shape[1] != 2 or rcv.shape[1] != 2:
            raise ValueError('src and rcv should be ndata x 2')

        if self.is_outside(src):
            raise ValueError('Source point outside grid')

        if self.is_outside(rcv):
            raise ValueError('Receiver outside grid')

        if slowness is not None:
            self.set_slowness(slowness)
        if xi is not None:
            self.set_xi(xi)
        if theta is not None:
            self.set_tilt_angle(theta)
        if Vp0 is not None:
            self.set_Vp0(Vp0)
        if Vs0 is not None:
            self.set_Vs0(Vs0)
        if delta is not None:
            self.set_delta(delta)
        if epsilon is not None:
            self.set_epsilon(epsilon)
        if gamma is not None:
            self.set_gamma(gamma)

        cdef vector[vector[sxz[double]]] vTx
        cdef vector[vector[sxz[double]]] vRx
        cdef vector[vector[double]] vt0
        cdef vector[vector[double]] vtt

        cdef vector[vector[vector[sxz[double]]]] r_data
        cdef vector[vector[vector[siv2[double]]]] l_data
        cdef size_t thread_nb

        cdef int i, j, k, n, nn, MM, NN

        vTx.resize(nTx)
        vRx.resize(nTx)
        vt0.resize(nTx)
        vtt.resize(nTx)
        if compute_L:
            l_data.resize(nTx)
        if return_rays:
            r_data.resize(nTx)

        iRx = []
        if nTx == 1:
            vTx[0].push_back(sxz[double](src[0,0], src[0,1]))
            for r in rcv:
                vRx[0].push_back(sxz[double](r[0], r[1]))
            vt0[0].push_back(t0[0])
            vtt[0].resize(rcv.shape[0])
            iRx.append(np.arange(rcv.shape[0]))
        elif aggregate_src:
            for t in Tx:
                vTx[0].push_back(sxz[double](t[0], t[1]))
            for t in t0:
                vt0[0].push_back(t)
            for r in rcv:
                vRx[0].push_back(sxz[double](r[0], r[1]))
            vtt[0].resize(rcv.shape[0])
            nTx = 1
            iRx.append(np.arange(rcv.shape[0]))
        else:
            if src.shape != rcv.shape:
                raise ValueError('src and rcv should be of equal size')

            for n in range(nTx):
                ind = np.sum(Tx[n,:] == src, axis=1) == 2
                iRx.append(np.nonzero(ind)[0])
                vTx[n].push_back(sxz[double](Tx[n,0], Tx[n,1]))
                vt0[n].push_back(t0[n])
                for r in rcv[ind,:]:
                    vRx[n].push_back(sxz[double](r[0], r[1]))
                vtt[n].resize(vRx[n].size())

        tt = np.zeros((rcv.shape[0],))
        if self._n_threads == 1:
            if compute_L==False and return_rays==False:
                for n in range(nTx):
                    self.grid.raytrace(vTx[n], vt0[n], vRx[n], vtt[n], 0)
            elif compute_L and return_rays:
                for n in range(nTx):
                    self.grid.raytrace(vTx[n], vt0[n], vRx[n], vtt[n], r_data[n], l_data[n], 0)
            elif compute_L:
                for n in range(nTx):
                    self.grid.raytrace(vTx[n], vt0[n], vRx[n], vtt[n], l_data[n], 0)
            else:
                for n in range(nTx):
                    self.grid.raytrace(vTx[n], vt0[n], vRx[n], vtt[n], r_data[n], 0)

        elif thread_no is not None:
            # we should be here for just one event
            assert nTx == 1
            # normally we should not need to compute L
            assert compute_L is False
            thread_nb = thread_no

            if return_rays:
                self.grid.raytrace(vTx[0], vt0[0], vRx[0], vtt[0], r_data[0], thread_nb)
                for nt in range(vtt[0].size()):
                    tt[nt] = vtt[0][nt]
                rays = []
                for n2 in range(vRx.size()):
                    r = np.empty((r_data[0][n2].size(), 2))
                    for nn in range(r_data[0][n2].size()):
                        r[nn, 0] = r_data[0][n2][nn].x
                        r[nn, 1] = r_data[0][n2][nn].z
                    rays.append(r)
                return tt, rays

            else:
                self.grid.raytrace(vTx[0], vt0[0], vRx[0], vtt[0], thread_nb)
                for nt in range(vtt[0].size()):
                    tt[nt] = vtt[0][nt]
                return tt

        else:
            if compute_L==False and return_rays==False:
                self.grid.raytrace(vTx, vt0, vRx, vtt)
            elif compute_L and return_rays:
                self.grid.raytrace(vTx, vt0, vRx, vtt, r_data, l_data)
            elif compute_L:
                self.grid.raytrace(vTx, vt0, vRx, vtt, l_data)
            else:
                self.grid.raytrace(vTx, vt0, vRx, vtt, r_data)

        for n in range(nTx):
            for nt in range(vtt[n].size()):
                tt[iRx[n][nt]] = vtt[n][nt]

        if return_rays:
            rays = [ [0.0] for n in range(rcv.shape[0])]
            for n in range(nTx):
                r = [ [0.0] for i in range(vRx[n].size())]
                for n2 in range(vRx[n].size()):
                    r[n2] = np.empty((r_data[n][n2].size(), 2))
                    for nn in range(r_data[n][n2].size()):
                        r[n2][nn, 0] = r_data[n][n2][nn].x
                        r[n2][nn, 1] = r_data[n][n2][nn].z
                for nt in range(vtt[n].size()):
                    rays[iRx[n][nt]] = r[nt]

        if compute_L:
            # first build an array of matrices, for each event
            L = []
            ncells = self.get_number_of_cells()
            for n in range(nTx):
                nnz = 0
                for ni in range(l_data[n].size()):
                    nnz += l_data[n][ni].size()

                if self.iso == b'i':
                    indptr = np.empty((vRx[n].size()+1,), dtype=np.int64)
                    indices = np.empty((nnz,), dtype=np.int64)
                    val = np.empty((nnz,))

                    k = 0
                    MM = vRx[n].size()
                    NN = ncells
                    for i in range(MM):
                        indptr[i] = k
                        for j in range(NN):
                            for nn in range(l_data[n][i].size()):
                                if l_data[n][i][nn].i == j:
                                    indices[k] = j
                                    val[k] = l_data[n][i][nn].v
                                    k += 1

                    indptr[MM] = k
                    L.append( sp.csr_matrix((val, indices, indptr), shape=(MM,NN)) )
                else:
                    nnz *= 2
                    indptr = np.empty((vRx[n].size()+1,), dtype=np.int64)
                    indices = np.empty((nnz,), dtype=np.int64)
                    val = np.empty((nnz,))

                    k = 0
                    MM = vRx[n].size()
                    NN = 2*ncells
                    for i in range(MM):
                        indptr[i] = k
                        for j in range(NN):
                            for nn in range(l_data[n][i].size()):
                                if l_data[n][i][nn].i == j:
                                    indices[k] = j
                                    val[k] = l_data[n][i][nn].v
                                    k += 1
                                elif l_data[n][i][nn].i+ncells == j:
                                    indices[k] = j
                                    val[k] = l_data[n][i][nn].v2
                                    k += 1

                    indptr[MM] = k
                    L.append( sp.csr_matrix((val, indices, indptr), shape=(MM,NN)) )
            # we want a single matrix
            tmp = sp.vstack(L)
            itmp = []
            for n in range(nTx):
                for nt in range(vtt[n].size()):
                    itmp.append(iRx[n][nt])
            L = tmp[itmp,:]

        if compute_L==False and return_rays==False:
            return tt
        elif compute_L and return_rays:
            return tt, rays, L
        elif compute_L:
            return tt, L
        else:
            return tt, rays

    def to_vtk(self, fields, filename):
        """
        to_vtk(fields, filename)

        Save grid variables and/or raypaths to VTK format

        Parameters
        ----------
        fields: dict
            dict of variables to save to file. Variables should be np.ndarray of
            size equal to either the number of nodes of the number of cells of
            the grid, or a list of raypath coordinates.
        filename: str
            Name of file without extension for saving (extension vtr will be
            added).  Raypaths are saved in separate files, and filename will
            be appended by the dict key and have a vtp extension.

        Notes
        -----
        VTK files can be visualized with Paraview (https://www.paraview.org)
        """
        xCoords = numpy_support.numpy_to_vtk(self.x)
        yCoords = numpy_support.numpy_to_vtk(np.array([0.0]))
        zCoords = numpy_support.numpy_to_vtk(self.z)

        rgrid = vtk.vtkRectilinearGrid()
        rgrid.SetDimensions(self._x.size(), 1, self._z.size())
        rgrid.SetXCoordinates(xCoords)
        rgrid.SetYCoordinates(yCoords)
        rgrid.SetZCoordinates(zCoords)

        save_grid = False
        for fn in fields:
            data = fields[fn]

            if isinstance(data, list):
                self._save_raypaths(data, filename+'_'+fn+'.vtp')
            else:
                save_grid = True
                scalar = vtk.vtkDoubleArray()
                scalar.SetName(fn)
                scalar.SetNumberOfComponents(1)
                scalar.SetNumberOfTuples(data.size)
                if data.size == self.get_number_of_nodes():
                    shape = (self._x.size(), self._z.size())
                    if data.ndim == 2:
                        if data.shape != shape:
                            raise ValueError('Field {0:s} has incorrect shape'.format(fn))
                        # VTK stores in 'F' order
                        tmp = data.flatten(order='F')
                    elif data.ndim == 1 or (data.ndim == 2 and data.shape[0] == data.size):
                        # 'C' order assumed, reshape back and flatten to 'F' order
                        tmp = data.reshape(shape).flatten(order='F')
                    else:
                        raise ValueError('Field {0:s} has incorrect ndim'.format(fn))
                    for n in range(data.size):
                        scalar.SetTuple1(n, tmp[n])
                    rgrid.GetPointData().AddArray(scalar)
                elif data.size == self.get_number_of_cells():
                    shape = (self._x.size()-1, self._z.size()-1)
                    if data.ndim == 2:
                        if data.shape != shape:
                            raise ValueError('Field {0:s} has incorrect shape'.format(fn))
                        # VTK stores in 'F' order
                        tmp = data.flatten(order='F')
                    elif data.ndim == 1 or (data.ndim == 2 and data.shape[0] == data.size):
                        # 'C' order assumed, reshape back and flatten to 'F' order
                        tmp = data.reshape(shape).flatten(order='F')
                    else:
                        raise ValueError('Field {0:s} has incorrect ndim'.format(fn))
                    for n in range(data.size):
                        scalar.SetTuple1(n, tmp[n])
                    rgrid.GetCellData().AddArray(scalar)
                else:
                    raise ValueError('Field {0:s} has incorrect size'.format(fn))

        if save_grid:
            writer = vtk.vtkXMLRectilinearGridWriter()
            writer.SetFileName(filename+'.vtr')
            writer.SetInputData(rgrid)
            writer.SetDataModeToBinary()
            writer.Update()

    def _save_raypaths(self, rays, filename):
        polydata = vtk.vtkPolyData()
        cellarray = vtk.vtkCellArray()
        pts = vtk.vtkPoints()
        npts = 0
        for n in range(len(rays)):
            npts += rays[n].shape[0]
        pts.SetNumberOfPoints(npts)
        npts = 0
        for n in range(len(rays)):
            for p in range(rays[n].shape[0]):
                pts.InsertPoint(npts, rays[n][p, 0], 0.0, rays[n][p, 1])
                npts += 1
        polydata.SetPoints(pts)
        npts = 0
        for n in range(len(rays)):
            line = vtk.vtkPolyLine()
            line.GetPointIds().SetNumberOfIds(rays[n].shape[0])
            for p in range(rays[n].shape[0]):
                line.GetPointIds().SetId(p, npts)
                npts += 1
            cellarray.InsertNextCell(line)
        polydata.SetLines(cellarray)
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(filename)
        writer.SetInputData(polydata)
        writer.SetDataModeToBinary()
        writer.Update()


    @staticmethod
    def data_kernel_straight_rays(np.ndarray[np.double_t, ndim=2] Tx,
                                  np.ndarray[np.double_t, ndim=2] Rx,
                                  np.ndarray[np.double_t, ndim=1] grx,
                                  np.ndarray[np.double_t, ndim=1] grz,
                                  aniso=False):
        """
        data_kernel_straight_rays(Tx, Rx, grx, grz, aniso=False) -> L

        Raytracing with straight rays in 2D

        Parameters
        ----------
        Tx : np.ndarray
            source coordinates, nTx by 2
                - 1st column contains X coordinates,
                - 2nd contains Z coordinates
        Rx : np.ndarray
            receiver coordinates, nTx by 2
                - 1st column contains X coordinates,
                - 2nd contains Z coordinates
        grx : np.ndarray
            grid node coordinates along x
        grz : np.ndarray
            grid node coordinates along z
        aniso : bool
            compute L for elliptically anisotropic medium (True) or isotropic
            medium (False)

        Returns
        -------
        L : scipy csr_matrix
            data kernel matrix (tt = L*slowness)

        Note
        ----
        Tx and Rx should contain the same number of rows, each row corresponding
        to a source-receiver pair
        """

        cdef size_t nTx = Tx.shape[0]
        cdef size_t n_grx = grx.shape[0]
        cdef size_t n_grz = grz.shape[0]
        cdef size_t nCells = (n_grx-1)*(n_grz-1)

        cdef double small = 1.e-10

        data_p = []
        indices_p = []
        indptr_p = []

        cdef size_t k = 0
        cdef size_t ix = 0
        cdef size_t iz = 0

        cdef double x, z, xs, xr, zs, zr, dtmp, dlx, dlz, m, b, ze, xe
        cdef int64_t iCell
        cdef bool up
        cdef Py_ssize_t n

        for n in range(nTx):
            indptr_p.append(k)

            xs = Tx[n, 0]
            zs = Tx[n, 1]
            xr = Rx[n, 0]
            zr = Rx[n, 1]

            if xs > xr:  # on va de s à r, on veut x croissant
                dtmp = xs
                xs = xr
                xr = dtmp
                dtmp = zs
                zs = zr
                zr = dtmp

            # point de départ
            x = xs
            z = zs

            if abs(zs-zr)<small: # rai horizontal

                for ix in range(n_grx-1):
                    if x < grx[ix+1]:
                        break
                for iz in range(n_grz-1):
                    if z < grz[iz+1]:
                        break

                while x < xr:
                    iCell = ix*(n_grz-1) + iz

                    dlx = ( grx[ix+1] if grx[ix+1]<xr else xr ) - x

                    indices_p.append(iCell)
                    data_p.append(dlx)
                    k += 1

                    ix += 1
                    x = grx[ix]

            elif abs(xs-xr) < small:   # rai vertical
                if zs > zr:
                    dtmp = zs
                    zs = zr
                    zr = dtmp
                z = zs

                for ix in range(n_grx-1):
                    if x < grx[ix+1]:
                        break
                for iz in range(n_grz-1):
                    if z < grz[iz+1]:
                        break

                while z < zr:
                    iCell = ix*(n_grz-1) + iz

                    dlz = ( grz[iz+1] if grz[iz+1]<zr else zr ) - z

                    if not aniso:
                        indices_p.append(iCell)
                    else:
                        indices_p.append(iCell + nCells)
                    data_p.append(dlz)
                    k += 1

                    iz += 1
                    z = grz[iz]
            else:   # rai oblique
                m = (zr-zs)/(xr-xs)
                b = zr - m*xr
                up = m>0

                for ix in range(n_grx-1):
                    if x < grx[ix+1]:
                        break
                for iz in range(n_grz-1):
                    if z < grz[iz+1]:
                        break

                while x < xr:

                    zi = m*grx[ix+1] + b

                    if up:
                        while z < zi and z < zr:
                            iCell = ix*(n_grz-1) + iz

                            ze = grz[iz+1] if grz[iz+1]<zi else zi
                            ze = ze if ze<zr else zr
                            xe = (ze-b)/m
                            dlx = xe - x
                            dlz = ze - z

                            if not aniso:
                                dl = sqrt( dlx*dlx + dlz*dlz )
                                indices_p.append(iCell)
                                data_p.append(dl)
                                k += 1
                            else:
                                indices_p.append(iCell)
                                data_p.append(dlx)
                                k += 1
                                indices_p.append(iCell+nCells)
                                data_p.append(dlz)
                                k += 1

                            x = xe
                            z = ze
                            if abs(z-grz[iz+1])<small:
                                iz += 1
                    else:
                        while z > zi and z > zr:
                            iCell = ix*(n_grz-1) + iz

                            ze = grz[iz] if grz[iz]>zi else zi
                            ze = ze if ze>zr else zr
                            xe = (ze-b)/m
                            dlx = xe - x
                            dlz = ze - z

                            if not aniso:
                                dl = sqrt( dlx*dlx + dlz*dlz )
                                indices_p.append(iCell)
                                data_p.append(dl)
                                k += 1
                            else:
                                indices_p.append(iCell)
                                data_p.append(dlx)
                                k += 1
                                indices_p.append(iCell+nCells)
                                data_p.append(dlz)
                                k += 1

                            x = xe
                            z = ze
                            if abs(z-grz[iz])<small:
                                iz -= 1

                    ix += 1
                    x = grx[ix]

        indptr_p.append(k)
        if not aniso:
            L = sp.csr_matrix((data_p, indices_p, indptr_p),
                              shape=(nTx, nCells))
        else:
            L = sp.csr_matrix((data_p, indices_p, indptr_p),
                              shape=(nTx, 2*nCells))
        return L


def _rebuild3d(x, y, z, constructor_params):
    (n_threads, cell_slowness, method, tt_from_rp, interp_vel, eps, maxit,
     weno, nsnx, nsny, nsnz, n_secondary,
     n_tertiary, radius_factor_tertiary, translate_grid) = constructor_params
    g = Grid3d(x, y, z, n_threads, cell_slowness, method, tt_from_rp,
               interp_vel, eps, maxit, weno, nsnx, nsny, nsnz, n_secondary,
               n_tertiary, radius_factor_tertiary, translate_grid)
    return g

def _rebuild2d(x, z, constructor_params):
    (n_threads, cell_slowness, method, aniso, eps, maxit, weno,
     rotated_template, nsnx, nsnz, n_secondary, n_tertiary,
     radius_factor_tertiary, tt_from_rp) = constructor_params
    g = Grid2d(x, z, n_threads, cell_slowness, method, aniso, eps, maxit, weno,
               rotated_template, nsnx, nsnz, n_secondary, n_tertiary,
               radius_factor_tertiary, tt_from_rp)
    return g

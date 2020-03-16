# -*- coding: utf-8 -*-
"""
    Raytracing on rectilinear grids

    This code is part of ttcr ( https://github.com/groupeLIAMG/ttcr )
"""

# distutils: language = c++

import numpy as np
cimport numpy as np
import scipy.sparse as sp

import vtk
from vtk.util import numpy_support

from ttcrpy.rgrid cimport Grid3D, Grid3Drcfs, Grid3Drcsp, Grid3Drcdsp, \
    Grid3Drnfs, Grid3Drnsp, Grid3Drndsp, Grid2D, Grid2Drcsp, Grid2Drcfs, \
    Grid2Drnsp, Grid2Drnfs, Grid2Drc, Grid2Drn


cdef class Grid3d:
    """
    class to perform raytracing with 3D rectilinear grids

    Slowness model can be defined in two ways:
        1) slowness constant within the voxels of the grid (the default)
        2) slowness defined at nodes of the grid

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
    nthreads: int
        number of threads for raytracing

    Constructor:

    Grid3d(x, y, z, nthreads=1, cell_slowness=1, method='FSM', tt_from_rp=1,
           interp_vel=0, eps=1.e-15, maxit=20, weno=1, nsnx=5, nsny=5, nsnz=5,
           n_secondary=2, n_tertiary=2, radius_tertiary=1.0) -> Grid3d

       Parameters
       ----------
       x : numpy ndarray
             node coordinates along x
       y : numpy ndarray
             node coordinates along y
       z : numpy ndarray
             node coordinates along z
       nthreads : int
             number of threads for raytracing (default is 1)
       cell_slowness : bool
             slowness defined for cells (True) or nodes (False) (default is 1)
       method : string
             raytracing method (default is FSM)
               'FSM' : fast marching method
               'SPM' : shortest path method
               'DSPM' : dynamic shortest path
       tt_from_rp : bool
             compute traveltimes from raypaths (SPM od DSPM) (default is 1)
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
       radius_tertiary : double
             radius of sphere around source that includes tertiary nodes (DSPM)
              (default is 1)
    """
    cdef vector[double] _x
    cdef vector[double] _y
    cdef vector[double] _z
    cdef double _dx
    cdef double _dy
    cdef double _dz
    cdef bool cell_slowness
    cdef size_t _nthreads
    cdef char method
    cdef Grid3D[double, uint32_t]* grid

    def __cinit__(self, np.ndarray[np.double_t, ndim=1] x,
                  np.ndarray[np.double_t, ndim=1] y,
                  np.ndarray[np.double_t, ndim=1] z,
                  size_t nthreads=1,
                  bool cell_slowness=1, str method='FSM',
                  bool tt_from_rp=1, bool interp_vel=0,
                  double eps=1.e-15, int maxit=20, bool weno=1,
                  uint32_t nsnx=5, uint32_t nsny=5, uint32_t nsnz=5,
                  uint32_t n_secondary=2, uint32_t n_tertiary=2,
                  double radius_tertiary=1.0):

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
        self._nthreads = nthreads

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
                                                            nthreads)
            elif method == 'SPM':
                self.method = b's'
                self.grid = new Grid3Drcsp[double,uint32_t,Cell[double,Node3Dcsp[double,uint32_t],sxyz[double]]](nx, ny, nz,
                                                            self._dx, self._dy, self._dz,
                                                            xmin, ymin, zmin,
                                                            nsnx, nsny, nsnz,
                                                            tt_from_rp,
                                                            nthreads)
            elif method == 'DSPM':
                self.method = b'd'
                self.grid = new Grid3Drcdsp[double,uint32_t,Cell[double,Node3Dc[double,uint32_t],sxyz[double]]](nx, ny, nz,
                                                           self._dx, self._dy, self._dz,
                                                           xmin, ymin, zmin,
                                                           n_secondary, tt_from_rp,
                                                           n_tertiary, radius_tertiary,
                                                           nthreads)
            else:
                raise ValueError('Method {0:s} undefined'.format(method))
        else:
            if method == 'FSM':
                self.method = b'f'
                self.grid = new Grid3Drnfs[double,uint32_t](nx, ny, nz, self._dx,
                                                            xmin, ymin, zmin,
                                                            eps, maxit, weno,
                                                            tt_from_rp, interp_vel,
                                                            nthreads)
            elif method == 'SPM':
                self.method = b's'
                self.grid = new Grid3Drnsp[double,uint32_t](nx, ny, nz,
                                                            self._dx, self._dy, self._dz,
                                                            xmin, ymin, zmin,
                                                            nsnx, nsny, nsnz,
                                                            tt_from_rp, interp_vel,
                                                            nthreads)
            elif method == 'DSPM':
                self.method = b'd'
                self.grid = new Grid3Drndsp[double,uint32_t](nx, ny, nz,
                                                             self._dx, self._dy, self._dz,
                                                             xmin, ymin, zmin,
                                                             n_secondary, tt_from_rp,
                                                             n_tertiary, radius_tertiary,
                                                             interp_vel, nthreads)
            else:
                raise ValueError('Method {0:s} undefined'.format(method))

    def __dealloc__(self):
        del self.grid

    @property
    def x(self):
        tmp = np.empty((self._x.size(),))
        cdef int n
        for n in range(self._x.size()):
            tmp[n] = self._x[n]
        return tmp

    @property
    def y(self):
        tmp = np.empty((self._y.size(),))
        cdef int n
        for n in range(self._y.size()):
            tmp[n] = self._y[n]
        return tmp

    @property
    def z(self):
        tmp = np.empty((self._z.size(),))
        cdef int n
        for n in range(self._z.size()):
            tmp[n] = self._z[n]
        return tmp

    @property
    def dx(self):
        return self._dx

    @property
    def dy(self):
        return self._dy

    @property
    def dz(self):
        return self._dz

    @property
    def shape(self):
        if self.cell_slowness:
            return (self._x.size()-1, self._y.size()-1, self._z.size()-1)
        else:
            return (self._x.size(), self._y.size(), self._z.size())

    @property
    def nthreads(self):
        return self._nthreads

    @property
    def nparams(self):
        if self.cell_slowness:
            return (self._x.size()-1) * (self._y.size()-1) * (self._z.size()-1)
        else:
            return self._x.size() * self._y.size() * self._z.size()

    def get_number_of_nodes(self):
        """
        Return number of nodes in grid
        """
        return self._x.size() * self._y.size() * self._z.size()

    def get_number_of_cells(self):
        """
        Return number of cells in grid
        """
        return (self._x.size()-1) * (self._y.size()-1) * (self._z.size()-1)

    def get_grid_traveltimes(self, thread_no=0):
        """
        Obtain traveltimes computed at primary grid nodes

        Parameters
        ----------
        thread_no : int
            thread used to computed traveltimes (default is 0)

        Returns
        -------
        tt: np ndarray, shape (nx, ny, nz)
        """
        if thread_no >= self._nthreads:
            raise ValueError('Thread number is larger than number of threads')
        cdef vector[double] tmp
        cdef int n
        self.grid.getTT(tmp, thread_no)
        tt = np.empty((tmp.size(),), order='F')
        for n in range(tmp.size()):
            tt[n] = tmp[n]
        shape = (self._x.size(), self._y.size(), self._z.size())
        return tt.reshape(shape)

    def ind(self, i, j, k):
        """
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
        Node index for a "flattened" grid
        """
        return (i*self._y.size() + j)*self._z.size() + k

    def indc(self, i, j, k):
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
        Cell index for a "flattened" grid
        """
        return (i*(self._y.size()-1) + j)*(self._z.size()-1) + k

    def is_outside(self, np.ndarray[np.double_t, ndim=2] pts):
        """
        Return True if at least one point outside grid
        """
        return ( np.min(pts[:,0]) < self._x.front() or np.max(pts[:,0]) > self._x.back() or
                np.min(pts[:,1]) < self._y.front() or np.max(pts[:,1]) > self._y.back() or
                np.min(pts[:,2]) < self._z.front() or np.max(pts[:,2]) > self._z.back() )

    def set_slowness(self, slowness):
        """
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

    def compute_D(self, coord):
        """
        Return matrix of interpolation weights for velocity data points
        constraint

        Parameters
        ----------
        coord : numpy ndarray
            coordinates of data points (npts x 3)

        Returns
        -------
        D : scipy csr_matrix
            Matrix of interpolation weights, of size npts x nparams
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
        Kx, Ky, Kz : tuple of sparse (csr) matrices for derivatives
            along x, y, & z
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
        Return slowness at source points

        Parameters
        ----------
        hypo : numpy ndarray with 5 columns
            1st column is event ID number
            2nd column is origin time
            3rd column is source easting
            4th column is source northing
            5th column is source elevation

        slowness : np ndarray, shape (nx, ny, nz) (optional)
            slowness at grid nodes or cells (depending on cell_slowness)
            slowness may also have been flattened (with default 'C' order)

        Returns
        -------
        s0 : numpy ndarray
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
        Perform raytracing

        Parameters
        ----------
        source : 2D numpy ndarray with 3, 4 or 5 columns
            see notes below
        rcv : 2D numpy ndarray with 3 columns
            Columns correspond to x, y and z coordinates
        slowness : np ndarray, shape (nx, ny, nz) (None by default)
            slowness at grid nodes or cells (depending on cell_slowness)
            slowness may also have been flattened (with default 'C' order)
                if None, slowness must have been assigned previously
        thread_no : int (None by default)
            Perform calculations in thread number "thread_no"
                if None, attempt to run in parallel if warranted by number of
                sources and value of nthreads in constructor
        aggregate_src : bool (False by default)
            it True, all source coordinates belong to a single event
        compute_L : bool (False by default)
            Compute matrices of partial derivative of travel time w/r to slowness
        compute_M : bool (False by default)
            Compute matrices of partial derivative of travel time w/r to velocity
                Note : compute_M and compute_L are mutually exclusive
        return_rays : bool (False by default)
            Return raypaths

        Returns
        -------
        tt : numpy ndarray
            travel times for the appropriate source-rcv  (see Notes below)
        rays : list of numpy ndarray (if return_rays is True)
            Coordinates of segments forming raypaths
        M : list of scipy csr_matrix
            list of matrices of partial derivative of travel time w/r to velocity
                the number of matrices is equal to the number of sources
        L : scipy csr_matrix
            Matrix of partial derivative of travel time w/r to slowness
            if input argument source has 5 columns, L is a list of matrices and
            the number of matrices is equal to the number of sources
            otherwise, L is a single csr_matrix

        Notes
        -----
        If source has 3 columns:
                Columns correspond to x, y and z coordinates
                Origin time (t0) is 0 for all points
        If source has 4 columns:
                1st column corresponds to origin times
                2nd, 3rd & 4th columns correspond to x, y and z coordinates
        If source has 5 columns:
                1st column corresponds to event ID
                2nd column corresponds to origin times
                3rd, 4th & 5th columns correspond to x, y and z coordinates

        For the latter case (5 columns), source and rcv should contain the same
            number of rows, each row corresponding to a source-receiver pair
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
            Tx = np.unique(source, axis=0)
            t0 = np.zeros((Tx.shape[0], 1))
            nTx = Tx.shape[0]
        elif source.shape[1] == 4:
            src = source[:,1:4]
            tmp = np.unique(source, axis=0)
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
        if nTx < self._nthreads or self._nthreads == 1:
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
                            if l_data[n][i][nn].i == j:
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
        xCoords = numpy_support.numpy_to_vtk(self.x)
        yCoords = numpy_support.numpy_to_vtk(self.y)
        zCoords = numpy_support.numpy_to_vtk(self.z)

        rgrid = vtk.vtkRectilinearGrid()
        rgrid.SetDimensions(self._x.size(), self._y.size(), self._z.size())
        rgrid.SetXCoordinates(xCoords)
        rgrid.SetYCoordinates(yCoords)
        rgrid.SetZCoordinates(zCoords)

        for fn in fields:
            data = fields[fn]

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
                elif data.ndim == 1:
                    # 'C' order assumed, reshape back and flatten to 'F' order
                    tmp = data.reshape(shape).flatten(order='F')
                else:
                    raise ValueError('Field {0:s} has incorrect ndim'.format(fn))
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
                elif data.ndim == 1:
                    # 'C' order assumed, reshape back and flatten to 'F' order
                    tmp = data.reshape(shape).flatten(order='F')
                else:
                    raise ValueError('Field {0:s} has incorrect ndim'.format(fn))
                for n in range(data.size):
                    scalar.SetTuple1(n, tmp[n])
                rgrid.GetCellData().AddArray(scalar)
            else:
                raise ValueError('Field {0:s} has incorrect size'.format(fn))

        writer = vtk.vtkXMLRectilinearGridWriter()
        writer.SetFileName(filename+'.vtr')
        writer.SetInputData(rgrid)
        writer.SetDataModeToBinary()
        writer.Update()

    @staticmethod
    def builder(filename, nthreads=1, method='FSM', tt_from_rp=1, interp_vel=0,
                eps=1.e-15, maxit=20, weno=1, nsnx=5, nsny=5, nsnz=5,
                n_secondary=2, n_tertiary=2, radius_tertiary=1.0):
        """
        Build instance of Grid3d from VTK file

        Parameters
        ----------
        filename : str
            Name of file holding rectilinear grid
            The grid must have point or cell attribute named either
            'Slowness', 'slowness', 'Velocity', 'velocity', or
            'P-wave velocity'

        Returns
        -------
        Instance of Grid3d
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

        g = Grid3d(x, y, z, nthreads, cell_slowness, method, tt_from_rp,
                   interp_vel, eps, maxit, weno, nsnx, nsny, nsnz,
                   n_secondary, n_tertiary, radius_tertiary)
        g.set_slowness(slowness)
        return g


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
    nthreads: int
        number of threads for raytracing
    """
    cdef vector[double] _x
    cdef vector[double] _z
    cdef double _dx
    cdef double _dz
    cdef bool cell_slowness
    cdef size_t _nthreads
    cdef char method

#    cdef Grid2D[double,uint32_t,sxz[double]]* grid
    cdef Grid2Drcsp[double,uint32_t, cell2d]* grid

    def __cinit__(self, np.ndarray[np.double_t, ndim=1] x,
                  np.ndarray[np.double_t, ndim=1] z,
                  size_t nthreads=1,
                  bool cell_slowness=1, str method='SPM',
                  double eps=1.e-15, int maxit=20, bool weno=1,
                  bool rotated_template=0, uint32_t nsnx=10, uint32_t nsnz=10):

        cdef uint32_t nx = x.size-1
        cdef uint32_t nz = z.size-1
        self._dx = x[1] - x[0]
        self._dz = z[1] - z[0]
        cdef double xmin = x[0]
        cdef double zmin = z[0]
        self.cell_slowness = cell_slowness
        self._nthreads = nthreads

        for val in x:
            self._x.push_back(val)
        for val in z:
            self._z.push_back(val)
        self._x.shrink_to_fit()
        self._z.shrink_to_fit()

        if cell_slowness:
            if method == 'SPM':
                self.grid = new Grid2Drcsp[double,uint32_t,cell2d](
                                nx, nz, self._dx, self._dz,
                                xmin, zmin, nsnx, nsnz, nthreads)
            elif method == 'FSM':
                raise ValueError('Implementation issue at this time')
                # self.grid = new Grid2Drcfs[double,uint32_t](nx, nz,
                #                                             self._dx, self._dz,
                #                                             xmin, zmin,
                #                                             eps, maxit, weno,
                #                                             rotated_template,
                #                                             nthreads)
            else:
                raise ValueError('Method {0:s} undefined'.format(method))
        else:
            #
            raise ValueError('Implementation issue at this time')
            # if method == 'SPM':
            #     self.grid = new Grid2Drnsp[double,uint32_t](nx, nz,
            #                                                 self._dx, self._dz,
            #                                                 xmin, zmin,
            #                                                 nsnx, nsnz,
            #                                                 nthreads)

    def __dealloc__(self):
        del self.grid

    @property
    def x(self):
        tmp = np.empty((self._x.size(),))
        cdef int n
        for n in range(self._x.size()):
            tmp[n] = self._x[n]
        return tmp

    @property
    def z(self):
        tmp = np.empty((self._z.size(),))
        cdef int n
        for n in range(self._z.size()):
            tmp[n] = self._z[n]
        return tmp

    @property
    def dx(self):
        return self._dx

    @property
    def dz(self):
        return self._dz

    @property
    def shape(self):
        if self.cell_slowness:
            return (self._x.size()-1, self._z.size()-1)
        else:
            return (self._x.size(), self._z.size())

    @property
    def nthreads(self):
        return self._nthreads

    @property
    def nparams(self):
        if self.cell_slowness:
            return (self._x.size()-1) * (self._z.size()-1)
        else:
            return self._x.size() * self._z.size()

    def get_number_of_nodes(self):
        """
        Return number of nodes in grid
        """
        return self._x.size() * self._z.size()

    def get_number_of_cells(self):
        """
        Return number of cells in grid
        """
        return (self._x.size()-1) * (self._z.size()-1)






    @staticmethod
    def data_kernel_straight_rays(np.ndarray[np.double_t, ndim=2] Tx,
                                  np.ndarray[np.double_t, ndim=2] Rx,
                                  np.ndarray[np.double_t, ndim=1] grx,
                                  np.ndarray[np.double_t, ndim=1] grz):
        """
        data_kernel_straight_rays(Tx, Rx, grx, grz) -> L

        Raytracing with straight rays in 2D

        Parameters
        ----------
        Tx : numpy ndarray
              source coordinates, nTx by 2
                1st column contains X coordinates,
                2nd contains Z coordinates
        Rx : numpy ndarray
              receiver coordinates, nTx by 2
                1st column contains X coordinates,
                2nd contains Z coordinates
        grx : numpy ndarray
               grid node coordinates along x
        grz : numpy ndarray
               grid node coordinates along z

        Returns
        -------
        L : scipy csr_matrix
               data kernel matrix (tt = L*slowness)
        """


        cdef size_t nTx = Tx.shape[0]
        cdef size_t n_grx = grx.shape[0]
        cdef size_t n_grz = grz.shape[0]

        cdef double small = 1.e-10

        cdef size_t nCells = (n_grx-1)*(n_grz-1)
        cdef size_t nLmax = 1 + nTx * n_grx * int(n_grz/2)
        cdef double percent_sp = (nLmax*1.0)/(nTx*nCells*1.0)

        data_p = []
        indices_p = []
        indptr_p = []

        cdef size_t k = 0
        cdef size_t ix, iz

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

                    indices_p.append(iCell)
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
                            dl = sqrt( dlx*dlx + dlz*dlz )

                            indices_p.append(iCell)
                            data_p.append(dl)
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
                            dl = sqrt( dlx*dlx + dlz*dlz )

                            indices_p.append(iCell)
                            data_p.append(dl)
                            k += 1

                            x = xe
                            z = ze
                            if abs(z-grz[iz])<small:
                                iz -= 1

                    ix += 1
                    x = grx[ix]

        indptr_p.append(k)
        L = sp.csr_matrix((data_p, indices_p, indptr_p), shape=(nTx, (n_grx-1)*(n_grz-1)))
        return L

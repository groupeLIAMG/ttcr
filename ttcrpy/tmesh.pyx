# -*- coding: utf-8 -*-
"""
Raytracing on unstructured triangular and tetrahedral meshes

This module contains two classes to perform traveltime computation and
raytracing on unstructured meshes:
    - `Mesh2d` for 2D media
    - `Mesh3d` for 3D media

    Three algorithms are implemented
        - the Shortest-Path Method
        - the Fast-Sweeping Method
        - the Dynamic Shortest-Path Method

    Slowness model can be defined in two ways:
        1) slowness constant within the voxels of the mesh (the default)
        2) slowness defined at nodes of the mesh

    This code is part of ttcr ( https://github.com/groupeLIAMG/ttcr )
"""

# distutils: language = c++

import numpy as np
cimport numpy as np
import scipy.sparse as sp

import vtk
from vtk.util import numpy_support

from ttcrpy.tmesh cimport Grid3D, Grid3Ducfs, Grid3Ducsp, Grid3Ducdsp, \
    Grid3Dunfs, Grid3Dunsp, Grid3Dundsp, Grid2D, Grid2Duc, Grid2Dun, \
    Grid2Ducsp, Grid2Ducfs, Grid2Dunsp, Grid2Dunfs

cdef extern from "verbose.h" namespace "ttcr" nogil:
    void setVerbose(int)

cdef extern from "utils_cython.h":
    int build_matrix_siv[T](size_t, size_t, vector[vector[siv[T]]]&, object)

def set_verbose(v):
    """Set verbosity level for C++ code

    Parameters
    ----------
    v: int
        verbosity level
    """
    setVerbose(v)


cdef class Mesh3d:
    """class to perform raytracing with tetrahedral meshes

    Attributes
    ----------
    nparams: int
        total number of parameters for grid
    n_threads: int
        number of threads for raytracing

    Constructor:

    Mesh3d(nodes, tetra, n_threads, cell_slowness, method, gradient_method, tt_from_rp, process_vel, eps, maxit, min_dist, n_secondary, n_tertiary, radius_factor_tertiary, translate_grid=False) -> Mesh3d

        Parameters
        ----------
        nodes : np.ndarray, shape (nnodes, 3)
            node coordinates
        tetra : np.ndarray of int, shape (ntetra, 4)
            indices of nodes forming the tetrahedra
        n_threads : int
            number of threads for raytracing (default is 1)
        cell_slowness : bool
            slowness defined for cells (True) or nodes (False) (default is 1)
        method : string
            raytracing method (default is FSM)
                - 'FSM' : fast marching method
                - 'SPM' : shortest path method
                - 'DSPM' : dynamic shortest path
        gradient_method : int
            method to compute traveltime gradient (default is 1)
                - 0 : least-squares first-order
                - 1 : least-squares second-order
                - 2 : Averaging-Based method
        tt_from_rp : bool
            compute traveltimes from raypaths (FSM or DSPM only) (default is 1)
        process_vel : bool
            process velocity instead of slowness at nodes when interpolating and
            computing matrix of partial derivative of traveltime w/r to model
            parameters (interpolation: for cell_slowness == False or FSM)
            (defauls is False)
        eps : double
            convergence criterion (FSM) (default is 1e-15)
        maxit : int
            max number of sweeping iterations (FSM) (default is 20)
        min_dist : double
            tolerance for backward raytracing (default is 1e-5)
        n_secondary : int
            number of secondary nodes (SPM & DSPM) (default is 2)
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
    cdef bool cell_slowness
    cdef size_t _n_threads
    cdef char method
    cdef bool tt_from_rp
    cdef bool process_vel
    cdef bool translate_grid
    cdef double eps
    cdef int maxit
    cdef int gradient_method
    cdef double min_dist
    cdef uint32_t n_secondary
    cdef uint32_t n_tertiary
    cdef double radius_factor_tertiary
    cdef vector[sxyz[double]] no
    cdef vector[tetrahedronElem[uint32_t]] tet
    cdef Grid3D[double, uint32_t]* grid

    def __cinit__(self, np.ndarray[np.double_t, ndim=2] nodes,
                  np.ndarray[np.int64_t, ndim=2] tetra,
                  size_t n_threads=1, bool cell_slowness=1,
                  str method='FSM', int gradient_method=1,
                  bool tt_from_rp=1, bool process_vel=0,
                  double eps=1.e-15, int maxit=20, double min_dist=1.e-5,
                  uint32_t n_secondary=2, uint32_t n_tertiary=2,
                  double radius_factor_tertiary=3.0,
                  bool translate_grid=0):

        self.cell_slowness = cell_slowness
        self._n_threads = n_threads
        self.tt_from_rp = tt_from_rp
        self.process_vel = process_vel
        self.eps = eps
        self.maxit = maxit
        self.gradient_method = gradient_method
        self.min_dist = min_dist
        self.n_secondary = n_secondary
        self.n_tertiary = n_tertiary
        self.radius_factor_tertiary = radius_factor_tertiary
        self.translate_grid = translate_grid

        cdef double source_radius = 0.0
        cdef vector[sxyz[double]] pts_ref
        cdef int n
        cdef use_edge_length = True

        if not nodes.flags['C_CONTIGUOUS']:
            nodes = np.ascontiguousarray(nodes)
        if not tetra.flags['C_CONTIGUOUS']:
            tetra = np.ascontiguousarray(tetra)

        for n in range(nodes.shape[0]):
            self.no.push_back(sxyz[double](nodes[n, 0],
                                           nodes[n, 1],
                                           nodes[n, 2]))
        for n in range(tetra.shape[0]):
            self.tet.push_back(tetrahedronElem[uint32_t](tetra[n, 0],
                                                         tetra[n, 1],
                                                         tetra[n, 2],
                                                         tetra[n, 3]))

        if method == 'FSM':
            xmin = np.min(nodes[:, 0])
            xmax = np.max(nodes[:, 0])
            ymin = np.min(nodes[:, 1])
            ymax = np.max(nodes[:, 1])
            zmin = np.min(nodes[:, 2])
            zmax = np.max(nodes[:, 2])
            pts_ref.push_back(sxyz[double](xmin, ymin, zmin))
            pts_ref.push_back(sxyz[double](xmin, ymin, zmax))
            pts_ref.push_back(sxyz[double](xmin, ymax, zmin))
            pts_ref.push_back(sxyz[double](xmin, ymax, zmax))
            pts_ref.push_back(sxyz[double](xmax, ymin, zmin))
            pts_ref.push_back(sxyz[double](xmax, ymin, zmax))
            pts_ref.push_back(sxyz[double](xmax, ymax, zmin))
            pts_ref.push_back(sxyz[double](xmax, ymax, zmax))

        if cell_slowness:
            if method == 'FSM':
                self.method = b'f'
                self.grid = new Grid3Ducfs[double,uint32_t](self.no, self.tet,
                                                            eps, maxit,
                                                            pts_ref, 2,
                                                            gradient_method,
                                                            tt_from_rp,
                                                            min_dist, n_threads,
                                                            translate_grid)
            elif method == 'SPM':
                self.method = b's'
                self.grid = new Grid3Ducsp[double,uint32_t](self.no, self.tet,
                                                            n_secondary,
                                                            tt_from_rp,
                                                            min_dist, n_threads,
                                                            translate_grid)
            elif method == 'DSPM':
                self.method = b'd'
                self.grid = new Grid3Ducdsp[double,uint32_t](self.no, self.tet,
                                                             n_secondary,
                                                             n_tertiary,
                                                             source_radius,
                                                             gradient_method,
                                                             tt_from_rp,
                                                             min_dist,
                                                             radius_factor_tertiary,
                                                             use_edge_length,
                                                             n_threads,
                                                             translate_grid)

            else:
                raise ValueError('Method {0:s} undefined'.format(method))
        else:
            if method == 'FSM':
                self.method = b'f'
                self.grid = new Grid3Dunfs[double,uint32_t](self.no, self.tet,
                                                            eps, maxit,
                                                            pts_ref, 2,
                                                            gradient_method,
                                                            process_vel,
                                                            tt_from_rp,
                                                            min_dist, n_threads,
                                                            translate_grid)
            elif method == 'SPM':
                self.method = b's'
                self.grid = new Grid3Dunsp[double,uint32_t](self.no, self.tet,
                                                            n_secondary,
                                                            process_vel,
                                                            tt_from_rp,
                                                            min_dist, n_threads,
                                                            translate_grid)
            elif method == 'DSPM':
                self.method = b'd'
                self.grid = new Grid3Dundsp[double,uint32_t](self.no, self.tet,
                                                             n_secondary,
                                                             n_tertiary,
                                                             source_radius,
                                                             process_vel,
                                                             gradient_method,
                                                             tt_from_rp,
                                                             min_dist,
                                                             radius_factor_tertiary,
                                                             use_edge_length,
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

        nodes = np.ndarray((self.no.size(), 3))
        tetra = np.ndarray((self.tet.size(), 4), dtype=int)
        cdef int n
        cdef int nn
        for n in range(nodes.shape[0]):
            nodes[n, 0] = self.no[n].x
            nodes[n, 1] = self.no[n].y
            nodes[n, 2] = self.no[n].z
        for n in range(tetra.shape[0]):
            for nn in range(4):
                tetra[n, nn] = self.tet[n].i[nn]

        constructor_params = (nodes, tetra, method, self.cell_slowness,
                              self._n_threads, self.tt_from_rp, self.process_vel,
                              self.eps, self.maxit, self.gradient_method,
                              self.min_dist, self.n_secondary, self.n_tertiary,
                              self.radius_factor_tertiary, self.translate_grid)
        return (_rebuild3d, constructor_params)

    @property
    def n_threads(self):
        """int: number of threads for raytracing"""
        return self._n_threads

    @property
    def nparams(self):
        """int: total number of parameters for mesh"""
        if self.cell_slowness:
            return self.tet.size()
        else:
            return self.no.size()

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
        cdef vector[sxyz[double]] vpts
        for n in range(pts.shape[0]):
            vpts.push_back(sxyz[double](pts[n, 0], pts[n, 1], pts[n, 2]))
        try:
            self.grid.checkPts(vpts)  # throws an except if 1 pt is outside
        except RuntimeError:
            return True

        return False

    def get_number_of_nodes(self):
        """
        Returns
        -------
        int:
            number of nodes in grid
        """
        return self.no.size()

    def get_number_of_cells(self):
        """
        Returns
        -------
        int:
            number of cells in grid
        """
        return self.tet.size()

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
        tt: np ndarray, shape (nnodes,)
            traveltimes
        """
        if thread_no >= self._n_threads:
            raise ValueError('Thread number is larger than number of threads')
        cdef vector[double] tmp
        cdef int n
        self.grid.getTT(tmp, thread_no)
        tt = np.empty((tmp.size(),))
        for n in range(tmp.size()):
            tt[n] = tmp[n]
        return tt

    def set_slowness(self, slowness):
        """
        set_slowness(slowness)

        Assign slowness to grid

        Parameters
        ----------
        slowness : np ndarray, shape (nparams, )
        """
        if slowness.size != self.nparams:
            raise ValueError('Slowness vector has wrong size')

        if not slowness.flags['C_CONTIGUOUS']:
            slowness = np.ascontiguousarray(slowness)

        cdef vector[double] slown
        cdef int i
        for i in range(slowness.size):
            slown.push_back(slowness[i])
        self.grid.setSlowness(slown)

    def set_velocity(self, velocity):
        """
        set_velocity(velocity)

        Assign velocity to grid

        Parameters
        ----------
        velocity : np ndarray, shape (nparams, )
        """
        if velocity.size != self.nparams:
            raise ValueError('velocity vector has wrong size')

        if not velocity.flags['C_CONTIGUOUS']:
            velocity = np.ascontiguousarray(velocity)

        cdef vector[double] slown
        cdef int i
        for i in range(velocity.size):
            slown.push_back(1./velocity[i])
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
        cdef vector[sxyz[double]] vpts
        for n in range(coord.shape[0]):
            vpts.push_back(sxyz[double](coord[n, 0], coord[n, 1], coord[n, 2]))

        cdef vector[vector[sijv[double]]] d_data

        self.grid.computeD(vpts, d_data)

        nnz = 0
        for ni in range(d_data.size()):
            nnz += d_data[ni].size()

        MM = vpts.size()
        NN = self.nparams
        indptr = np.empty((MM+1,), dtype=np.int64)
        indices = np.empty((nnz,), dtype=np.int64)
        val = np.empty((nnz,))

        k = 0
        for i in range(MM):
            indptr[i] = k
            for j in range(NN):
                for nn in range(d_data[i].size()):
                    if d_data[i][nn].i == i and d_data[i][nn].j == j:
                        indices[k] = j
                        val[k] = d_data[i][nn].v
                        k += 1

        indptr[MM] = k
        return sp.csr_matrix((val, indices, indptr), shape=(MM,NN))

    def compute_K(self, order=2, taylor_order=2, weighting=True, squared=True,
                  s0inside=False, additional_points=0):
        """
        Compute smoothing matrices (spatial derivative)

        Parameters
        ----------
        order : int
            order of derivative (1 or 2, 2 by default)
        taylor_order : int
            order of taylors series expansion (1 or 2, 2 by default)
        weighting : bool
            apply inverse distance weighting (True by default)
        squared : bool
            Second derivative evaluated by taking the square of first
            derivative.  Applied only if order == 2 (True by default)
        s0inside : bool
            (experimental) ignore slowness value at local node (value is a
            filtered estimate) (False by default)
        additional_points : int
            use additional points to compute derivatives (minimum sometimes
            yield noisy results when rays are close to domain limits)
            (0 by default)

        Returns
        -------
        Kx, Ky, Kz : :obj:`tuple` of :obj:`csr_matrix`
            matrices for derivatives along x, y, & z
        """

        cdef int o = order
        cdef int to = taylor_order
        cdef bool w = weighting
        cdef bool s = s0inside
        cdef int add_pts = additional_points
        cdef vector[vector[vector[siv[double]]]] k_data

        if order == 2 and squared:
            o = 1
        self.grid.computeK(k_data, o, to, w, s, add_pts)

        K = []
        cdef size_t MM = self.nparams
        cdef size_t NN = self.nparams
        for nk in range(3):
            m_tuple = ([0.0], [0.0], [0.0])
            build_matrix_siv(MM, NN, k_data[nk], m_tuple)
            K.append( sp.csr_matrix(m_tuple, shape=(MM,NN)) )
            if order == 2 and squared:
                K[-1] = K[-1] * K[-1]

        return tuple(K)

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
        slowness : np ndarray, shape (nparams, ) (optional)
            slowness at grid nodes or cells (depending on cell_slowness)

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
                 aggregate_src=False, compute_L=False, return_rays=False):
        """
        raytrace(source, rcv, slowness=None, thread_no=None, aggregate_src=False, compute_L=False, return_rays=False) -> tt, rays, L

        Perform raytracing

        Parameters
        ----------
        source : 2D np.ndarray with 3, 4 or 5 columns
            see notes below
        rcv : 2D np.ndarray with 3 columns
            Columns correspond to x, y and z coordinates
        slowness : np ndarray, (None by default)
            slowness at grid nodes or cells (depending on cell_slowness)
            if None, slowness must have been assigned previously
        thread_no : int (None by default)
            Perform calculations in thread number "thread_no"
            if None, attempt to run in parallel if warranted by number of
            sources and value of n_threads in constructor
        aggregate_src : bool (False by default)
            if True, all source coordinates belong to a single event
        compute_L : bool (False by default)
            Compute matrices of partial derivative of travel time w/r to
            slowness (or velocity if process_vel == True in constructor)
        return_rays : bool (False by default)
            Return raypaths

        Returns
        -------
        tt : np.ndarray
            travel times for the appropriate source-rcv  (see Notes below)
        rays : :obj:`list` of :obj:`np.ndarray`
            Coordinates of segments forming raypaths (if return_rays is True)
        L :  :obj:`list` of :obj:`csr_matrix`  or  scipy csr_matrix
            Matrix of partial derivative of travel time w/r to slowness.
            if input argument source has 5 columns or if slowness is defined at
            nodes, L is a list of matrices and the number of matrices is equal
            to the number of sources otherwise, L is a single csr_matrix

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

        cdef int i, n, nn, n2, nt, nnz, num_row, index

        vTx.resize(nTx)
        vRx.resize(nTx)
        vt0.resize(nTx)
        vtt.resize(nTx)
        if compute_L and self.cell_slowness:
            raise NotImplementedError('compute_L not implemented for mesh with slowness defined in cells')
        elif compute_L and not self.cell_slowness:
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

        tt = np.zeros((rcv.shape[0],))
        if self._n_threads == 1:
            if compute_L==False and return_rays==False:
                for n in range(nTx):
                    self.grid.raytrace(vTx[n], vt0[n], vRx[n], vtt[n], 0)
            elif compute_L and return_rays:
                if self.cell_slowness:
                    # TODO: implement this in C++ codebase!
                    for n in range(nTx):
                        self.grid.raytrace(vTx[n], vt0[n], vRx[n], vtt[n], r_data[n], l_data[n], 0)
                else:
                    for n in range(nTx):
                        self.grid.raytrace(vTx[n], vt0[n], vRx[n], vtt[n], r_data[n], m_data[n], 0)
            elif compute_L:
                if self.cell_slowness:
                    for n in range(nTx):
                        self.grid.raytrace(vTx[n], vt0[n], vRx[n], vtt[n], l_data[n], 0)
                else:
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
            if compute_L==False and return_rays==False:
                self.grid.raytrace(vTx, vt0, vRx, vtt)
            elif compute_L and return_rays:
                if self.cell_slowness:
                    self.grid.raytrace(vTx, vt0, vRx, vtt, r_data, l_data)
                else:
                    self.grid.raytrace(vTx, vt0, vRx, vtt, r_data, m_data)
            elif compute_L:
                if self.cell_slowness:
                    self.grid.raytrace(vTx, vt0, vRx, vtt, l_data)
                else:
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

        if compute_L and not self.cell_slowness:
            # we return array of matrices, one for each event
            L = []
            for n in range(nTx):
                nnz = 0
                num_row = vRx[n].size()+1
                for ni in range(m_data[n].size()):
                    nnz += m_data[n][ni].size()
                    if m_data[n][ni].size() == 0:
                         num_row -= 1
                indptr = np.empty((num_row,), dtype=np.int64)
                indices = np.empty((nnz,), dtype=np.int64)
                val = np.empty((nnz,))

                k = 0
                NN = self.get_number_of_nodes()
                index = 0
                for i in range(m_data[n].size()):
                    if m_data[n][i].size() == 0:
                         continue
                    indptr[index] = k
                    index +=1
                    for nn in range(m_data[n][i].size()):
                        indices[k] = m_data[n][i][nn].j
                        val[k] = m_data[n][i][nn].v
                        k += 1

                indptr[index] = k
                L.append(sp.csr_matrix((val, indices, indptr),
                         shape=(indptr.size - 1, NN)))

        if compute_L==False and return_rays==False:
            return tt
        elif compute_L and return_rays:
            return tt, rays, L
        elif compute_L:
            return tt, L
        else:
            return tt, rays

    def data_kernel_straight_rays(self, np.ndarray[np.double_t, ndim=2] Tx,
                                  np.ndarray[np.double_t, ndim=2] Rx):
        """
        data_kernel_straight_rays(Tx, Rx) -> L

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
        Returns
        -------
        L : scipy csr_matrix
            data kernel matrix (tt = L*slowness)

        Note
        ----
        Tx and Rx should contain the same number of rows, each row corresponding
        to a source-receiver pair
        """
        if Tx.ndim != 2 or Rx.ndim != 2:
            raise ValueError('Tx and Rx should be 2D arrays')
        if Tx.shape[1] != 3 or Rx.shape[1] != 3:
            raise ValueError('Tx and Rx should be ndata x 3')
        if Tx.shape[0] != Rx.shape[0]:
            raise ValueError('Tx and Rx should be of equal size')

        cdef long int n, ni, nnz, k
        cdef vector[sxyz[double]] vTx
        cdef vector[sxyz[double]] vRx
        cdef vector[vector[siv[double]]] l_data

        vTx.resize(Tx.shape[0])
        vRx.resize(Tx.shape[0])
        l_data.resize(Tx.shape[0])

        for n in range(Tx.shape[0]):
            vTx.push_back(sxyz[double](Tx[n, 0], Tx[n, 1], Tx[n, 2]))
            vRx.push_back(sxyz[double](Rx[n, 0], Rx[n, 1], Rx[n, 2]))

        self.grid.getStraightRays(vTx, vRx, l_data)

        nnz = 0
        for ni in range(l_data.size()):
            nnz += l_data[ni].size()
        indptr = np.empty((vRx.size()+1,), dtype=np.int64)
        indices = np.empty((nnz,), dtype=np.int64)
        val = np.empty((nnz,))

        k = 0
        MM = vRx.size()
        NN = self.get_number_of_cells()
        for i in range(MM):
            indptr[i] = k
            for j in range(NN):
                for nn in range(l_data[i].size()):
                    if l_data[i][nn].i == j:
                        indices[k] = j
                        val[k] = l_data[i][nn].v
                        k += 1
        indptr[MM] = k
        return sp.csr_matrix((val, indices, indptr), shape=(MM,NN))

    def to_vtk(self, fields, filename):
        """
        to_vtk(fields, filename)

        Save mesh variables and/or raypaths to VTK format

        Parameters
        ----------
        fields: dict
            dict of variables to save to file. Variables should be np.ndarray of
            size equal to either the number of nodes of the number of cells of
            the mesh, or a list of raypath coordinates.
        filename: str
            Name of file without extension for saving (extension vtu will be
            added).  Raypaths are saved in separate files, and filename will
            be appended by the dict key and have a vtp extension.

        Notes
        -----
        VTK files can be visualized with Paraview (https://www.paraview.org)
        """
        cdef int n, nn
        ugrid = vtk.vtkUnstructuredGrid()
        tPts = vtk.vtkPoints()
        tPts.SetNumberOfPoints(self.no.size())
        for n in range(self.no.size()):
            tPts.InsertPoint(n, self.no[n].x, self.no[n].y, self.no[n].z)
        ugrid.SetPoints(tPts)
        tet = vtk.vtkTetra()
        for n in range(self.tet.size()):
            for nn in range(4):
                tet.GetPointIds().SetId(nn, self.tet[n].i[nn])
            ugrid.InsertNextCell(tet.GetCellType(), tet.GetPointIds())

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
                    for n in range(data.size):
                        scalar.SetTuple1(n, data[n])
                    ugrid.GetPointData().AddArray(scalar)
                elif data.size == self.get_number_of_cells():
                    for n in range(data.size):
                        scalar.SetTuple1(n, data[n])
                    ugrid.GetCellData().AddArray(scalar)
                else:
                    raise ValueError('Field {0:s} has incorrect size'.format(fn))

        if save_grid:
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetFileName(filename+'.vtu')
            writer.SetInputData(ugrid)
            writer.SetDataModeToBinary()
            writer.Update()

    def  _save_raypaths(self, rays, filename):
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
    def builder(filename, size_t n_threads=1, bool cell_slowness=1,
                str method='FSM', int gradient_method=1,
                bool tt_from_rp=1, bool process_vel=0,
                double eps=1.e-15, int maxit=20, double min_dist=1.e-5,
                uint32_t n_secondary=2, uint32_t n_tertiary=2,
                double radius_factor_tertiary=3.0, bool translate_grid=0):
        """
        builder(filename, n_threads, cell_slowness, method, gradient_method, tt_from_rp, process_vel, eps, maxit, min_dist, n_secondary, n_tertiary, radius_factor_tertiary, translate_grid=0)

        Build instance of Mesh3d from VTK file

        Parameters
        ----------
        filename : str
            Name of file holding a vtkUnstructuredGrid.
            The grid must have point or cell attribute named either
            'Slowness', 'slowness', 'Velocity', 'velocity', or
            'P-wave velocity'.  All cells must be of type vtkTetra

        Other parameters are defined in Constructor

        Returns
        -------
        mesh: :obj:`Mesh3d`
            mesh instance
        """

        cdef int n, nn
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(filename)
        reader.Update()

        data = reader.GetOutput()
        nodes = numpy_support.vtk_to_numpy(data.GetPoints().GetData())
        for n in range(data.GetNumberOfCells()):
            if data.GetCellType(n) != vtk.VTK_TETRA:
                raise ValueError('{0:s} should only contain tetrahedra')
        tet = np.ndarray((data.GetNumberOfCells(), 4))
        for n in range(data.GetNumberOfCells()):
            for nn in range(4):
                tet[n, nn] = data.GetCell(n).GetPointIds().GetId(nn)

        names = ('Slowness', 'slowness', 'Velocity', 'velocity',
                 'P-wave velocity')
        for name in names:
            if data.GetPointData().HasArray(name):
                cell_slowness = 0
                data = numpy_support.vtk_to_numpy(data.GetPointData().GetArray(name))
                break
            if data.GetCellData().HasArray(name):
                cell_slowness = 1
                data = numpy_support.vtk_to_numpy(data.GetCellData().GetArray(name))
                break
        else:
            raise ValueError('File should contain slowness or velocity data')

        if 'lowness' in name:
            slowness = data
        else:
            slowness = 1.0 / data

        m = Mesh3d(nodes, tet, n_threads, cell_slowness, method, gradient_method,
                   tt_from_rp, process_vel, eps, maxit, min_dist, n_secondary,
                   n_tertiary, radius_factor_tertiary, translate_grid)
        m.set_slowness(slowness)
        return m


cdef class Mesh2d:
    """class to perform raytracing with triangular meshes

    Attributes
    ----------
    nparams: int
        total number of parameters for grid
    n_threads: int
        number of threads for raytracing

    Constructor:

    Mesh2d(nodes, triangles, n_threads, cell_slowness, method, eps, maxit, process_obtuse, n_secondary, n_tertiary, radius_factor_tertiary, tt_from_rp) -> Mesh2d

        Parameters
        ----------
        nodes : np.ndarray, shape (nnodes, 2)
            node coordinates
        triangles : np.ndarray of int, shape (ntriangles, 3)
            indices of nodes forming the triangles
        n_threads : int
            number of threads for raytracing (default is 1)
        cell_slowness : bool
            slowness defined for cells (True) or nodes (False) (default is 1)
        method : string
            raytracing method (default is FSM)
                - 'FSM' : fast marching method
                - 'SPM' : shortest path method
                - 'DSPM' : dynamic shortest path
        eps : double
            convergence criterion (FSM) (default is 1e-15)
        maxit : int
            max number of sweeping iterations (FSM) (default is 20)
        process_obtuse : bool
            use method of Qian et al (2007) to improve accuracy for triangles
            with obtuse angle (default is True)
        n_secondary : int
            number of secondary nodes (SPM) (default is 5)
        n_tertiary : int
            number of tertiary nodes (DSPM) (default is 2)
        radius_factor_tertiary : double
            multiplication factor used to compute radius of sphere around source
            that includes tertiary nodes (DSPM).  The radius is the average edge
            length multiplied by this factor (default is 2)
        tt_from_rp : bool
            compute traveltimes using raypaths (default is False)
    """
    cdef bool cell_slowness
    cdef bool process_obtuse
    cdef bool tt_from_rp
    cdef size_t _n_threads
    cdef double eps
    cdef int maxit
    cdef char method
    cdef uint32_t n_secondary
    cdef uint32_t n_tertiary
    cdef double radius_factor_tertiary
    cdef vector[sxz[double]] no
    cdef vector[triangleElem[uint32_t]] tri
    cdef Grid2D[double, uint32_t,sxz[double]]* grid

    def __cinit__(self, np.ndarray[np.double_t, ndim=2] nodes,
                  np.ndarray[np.int64_t, ndim=2] triangles,
                  size_t n_threads=1, bool cell_slowness=1,
                  str method='FSM', double eps=1.e-15, int maxit=20,
                  bool process_obtuse=1, uint32_t n_secondary=5,
                  uint32_t n_tertiary=2, double radius_factor_tertiary=2.0,
                  bool tt_from_rp=0):

        self.cell_slowness = cell_slowness
        self._n_threads = n_threads
        self.eps = eps
        self.maxit = maxit
        self.process_obtuse = process_obtuse
        self.n_secondary = n_secondary
        self.n_tertiary = n_tertiary
        self.radius_factor_tertiary = radius_factor_tertiary
        self.tt_from_rp = tt_from_rp

        cdef vector[sxz[double]] pts_ref
        cdef int n
        cdef use_edge_length = True

        if not nodes.flags['C_CONTIGUOUS']:
            nodes = np.ascontiguousarray(nodes)
        if not triangles.flags['C_CONTIGUOUS']:
            triangles = np.ascontiguousarray(triangles)

        for n in range(nodes.shape[0]):
            self.no.push_back(sxz[double](nodes[n, 0],
                                          nodes[n, 1]))
        for n in range(triangles.shape[0]):
            self.tri.push_back(triangleElem[uint32_t](triangles[n, 0],
                                                      triangles[n, 1],
                                                      triangles[n, 2]))

        if method == 'FSM':
            xmin = np.min(nodes[:, 0])
            xmax = np.max(nodes[:, 0])
            zmin = np.min(nodes[:, 1])
            zmax = np.max(nodes[:, 1])
            pts_ref.push_back(sxz[double](xmin, zmin))
            pts_ref.push_back(sxz[double](xmin, zmax))
            pts_ref.push_back(sxz[double](xmax, zmin))
            pts_ref.push_back(sxz[double](xmax, zmax))

        if cell_slowness:
            if method == 'FSM':
                self.method = b'f'
                self.grid = new Grid2Ducfs[double,uint32_t,Node2Dc[double,uint32_t],sxz[double]](self.no,
                                                                                                 self.tri,
                                                                                                 eps,
                                                                                                 maxit,
                                                                                                 pts_ref, 2,
                                                                                                 tt_from_rp,
                                                                                                 n_threads,
                                                                                                 process_obtuse)
            elif method == 'SPM':
                self.method = b's'
                self.grid = new Grid2Ducsp[double,uint32_t,Node2Dcsp[double,uint32_t],sxz[double]](self.no,
                                                                                                   self.tri,
                                                                                                   n_secondary,
                                                                                                   tt_from_rp,
                                                                                                   n_threads)
            elif method == 'DSPM':
                self.method = b'd'
                self.grid = new Grid2Ducdsp[double,uint32_t,sxz[double]](self.no,
                                                                         self.tri,
                                                                         n_secondary,
                                                                         n_tertiary,
                                                                         radius_factor_tertiary,
                                                                         tt_from_rp,
                                                                         use_edge_length,
                                                                         n_threads)
            else:
                raise ValueError('Method {0:s} undefined'.format(method))
        else:
            if method == 'FSM':
                self.method = b'f'
                self.grid = new Grid2Dunfs[double,uint32_t,Node2Dn[double,uint32_t],sxz[double]](self.no,
                                                                                                 self.tri,
                                                                                                 eps,
                                                                                                 maxit,
                                                                                                 pts_ref, 2,
                                                                                                 tt_from_rp,
                                                                                                 n_threads,
                                                                                                 process_obtuse)
            elif method == 'SPM':
                self.method = b's'
                self.grid = new Grid2Dunsp[double,uint32_t,Node2Dnsp[double,uint32_t],sxz[double]](self.no,
                                                                                                   self.tri,
                                                                                                   n_secondary,
                                                                                                   tt_from_rp,
                                                                                                   n_threads)
            elif method == 'DSPM':
                self.method = b'd'
                self.grid = new Grid2Dundsp[double,uint32_t,sxz[double]](self.no,
                                                                         self.tri,
                                                                         n_secondary,
                                                                         n_tertiary,
                                                                         radius_factor_tertiary,
                                                                         tt_from_rp,
                                                                         use_edge_length,
                                                                         n_threads)
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

        nodes = np.ndarray((self.no.size(), 2))
        triangles = np.ndarray((self.tri.size(), 3), dtype=int)
        cdef int n
        cdef int nn
        for n in range(nodes.shape[0]):
            nodes[n, 0] = self.no[n].x
            nodes[n, 1] = self.no[n].z
        for n in range(triangles.shape[0]):
            for nn in range(3):
                triangles[n, nn] = self.tri[n].i[nn]

        constructor_params = (nodes, triangles,
                              method, self.cell_slowness,
                              self._n_threads,
                              self.eps, self.maxit, self.process_obtuse,
                              self.n_secondary, self.n_tertiary,
                              self.radius_factor_tertiary, self.tt_from_rp)
        return (_rebuild2d, constructor_params)

    @property
    def n_threads(self):
        """int: number of threads for raytracing"""
        return self._n_threads

    @property
    def nparams(self):
        """int: total number of parameters for mesh"""
        if self.cell_slowness:
            return self.tri.size()
        else:
            return self.no.size()

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
        return self.no.size()

    def get_number_of_cells(self):
        """
        Returns
        -------
        int:
            number of cells in grid
        """
        return self.tri.size()

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
        tt: np ndarray, shape (nnodes,)
            traveltimes
        """
        if thread_no >= self._n_threads:
            raise ValueError('Thread number is larger than number of threads')
        cdef vector[double] tmp
        cdef int n
        self.grid.getTT(tmp, thread_no)
        tt = np.empty((tmp.size(),))
        for n in range(tmp.size()):
            tt[n] = tmp[n]
        return tt

    def set_slowness(self, slowness):
        """
        set_slowness(slowness)

        Assign slowness to grid

        Parameters
        ----------
        slowness : np ndarray, shape (nparams, )
        """
        if slowness.size != self.nparams:
            raise ValueError('Slowness vector has wrong size')

        if not slowness.flags['C_CONTIGUOUS']:
            slowness = np.ascontiguousarray(slowness)

        cdef vector[double] slown
        cdef int i
        for i in range(slowness.size):
            slown.push_back(slowness[i])
        self.grid.setSlowness(slown)

    def set_velocity(self, velocity):
        """
        set_velocity(velocity)

        Assign velocity to grid

        Parameters
        ----------
        velocity : np ndarray, shape (nparams, )
        """
        if velocity.size != self.nparams:
            raise ValueError('velocity vector has wrong size')

        if not velocity.flags['C_CONTIGUOUS']:
            velocity = np.ascontiguousarray(velocity)

        cdef vector[double] slown
        cdef int i
        for i in range(velocity.size):
            slown.push_back(1./velocity[i])
        self.grid.setSlowness(slown)

    def raytrace(self, source, rcv, slowness=None, thread_no=None,
                 aggregate_src=False, return_rays=False):
        """
        raytrace(source, rcv, slowness=None, thread_no=None, aggregate_src=False, return_rays=False) -> tt, rays

        Perform raytracing

        Parameters
        ----------
        source : 2D np.ndarray with 2 or 3 columns
            see notes below
        rcv : 2D np.ndarray with 2 columns
            Columns correspond to x and z coordinates
        slowness : np ndarray, (None by default)
            slowness at grid nodes or cells (depending on cell_slowness)
            if None, slowness must have been assigned previously
        thread_no : int (None by default)
            Perform calculations in thread number "thread_no"
            if None, attempt to run in parallel if warranted by number of
            sources and value of n_threads in constructor
        aggregate_src : bool (False by default)
            if True, all source coordinates belong to a single event
        return_rays : bool (False by default)
            Return raypaths

        Returns
        -------
        tt : np.ndarray
            travel times for the appropriate source-rcv  (see Notes below)
        rays : :obj:`list` of :obj:`np.ndarray`
            Coordinates of segments forming raypaths (if return_rays is True)

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

        if slowness is not None:
            self.set_slowness(slowness)

        cdef vector[vector[sxz[double]]] vTx
        cdef vector[vector[sxz[double]]] vRx
        cdef vector[vector[double]] vt0
        cdef vector[vector[double]] vtt

        cdef vector[vector[vector[sxz[double]]]] r_data
        cdef size_t thread_nb

        cdef int i, n, n2, nt

        vTx.resize(nTx)
        vRx.resize(nTx)
        vt0.resize(nTx)
        vtt.resize(nTx)
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
            if return_rays==False:
                for n in range(nTx):
                    self.grid.raytrace(vTx[n], vt0[n], vRx[n], vtt[n], 0)
            else:
                for n in range(nTx):
                    self.grid.raytrace(vTx[n], vt0[n], vRx[n], vtt[n], r_data[n], 0)

        elif thread_no is not None:
            # we should be here for just one event
            assert nTx == 1
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
            if return_rays==False:
                self.grid.raytrace(vTx, vt0, vRx, vtt)
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

        if return_rays==False:
            return tt
        else:
            return tt, rays

    def to_vtk(self, fields, filename):
        """
        to_vtk(fields, filename)

        Save mesh variables and/or raypaths to VTK format

        Parameters
        ----------
        fields: dict
            dict of variables to save to file. Variables should be np.ndarray of
            size equal to either the number of nodes of the number of cells of
            the mesh, or a list of raypath coordinates.
        filename: str
            Name of file without extension for saving (extension vtu will be
            added).  Raypaths are saved in separate files, and filename will
            be appended by the dict key and have a vtp extension.

        Notes
        -----
        VTK files can be visualized with Paraview (https://www.paraview.org)
        """
        cdef int n, nn
        ugrid = vtk.vtkUnstructuredGrid()
        tPts = vtk.vtkPoints()
        tPts.SetNumberOfPoints(self.no.size())
        for n in range(self.no.size()):
            tPts.InsertPoint(n, self.no[n].x, 0.0, self.no[n].z)
        ugrid.SetPoints(tPts)
        tri = vtk.vtkTriangle()
        for n in range(self.tri.size()):
            for nn in range(3):
                tri.GetPointIds().SetId(nn, self.tri[n].i[nn])
            ugrid.InsertNextCell(tri.GetCellType(), tri.GetPointIds())

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
                    for n in range(data.size):
                        scalar.SetTuple1(n, data[n])
                    ugrid.GetPointData().AddArray(scalar)
                elif data.size == self.get_number_of_cells():
                    for n in range(data.size):
                        scalar.SetTuple1(n, data[n])
                    ugrid.GetCellData().AddArray(scalar)
                else:
                    raise ValueError('Field {0:s} has incorrect size'.format(fn))

        if save_grid:
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetFileName(filename+'.vtu')
            writer.SetInputData(ugrid)
            writer.SetDataModeToBinary()
            writer.Update()

    def  _save_raypaths(self, rays, filename):
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
                pts.InsertPoint(npts, rays[n][p, 0], 0.0, rays[n][p, 2])
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
    def builder(filename, size_t n_threads=1, bool cell_slowness=1,
                str method='FSM',double eps=1.e-15, int maxit=20,
                bool process_obtuse=1, uint32_t n_secondary=5,
                uint32_t n_tertiary=2, double radius_factor_tertiary=3.0,
                bool tt_from_rp=0):
        """
        builder(filename, n_threads, cell_slowness, method, eps, maxit, process_obtuse, n_secondary, n_tertiary, radius_factor_tertiary, tt_from_rp)

        Build instance of Mesh2d from VTK file

        Parameters
        ----------
        filename : str
            Name of file holding a vtkUnstructuredGrid.
            The grid must have point or cell attribute named either
            'Slowness', 'slowness', 'Velocity', 'velocity', or
            'P-wave velocity'.  All cells must be of type vtkTriangle

        Other parameters are defined in Constructor

        Returns
        -------
        mesh: :obj:`Mesh2d`
            mesh instance
        """

        cdef int n, nn
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(filename)
        reader.Update()

        data = reader.GetOutput()
        nodes = numpy_support.vtk_to_numpy(data.GetPoints().GetData())
        nod = np.ascontiguousarray(nodes[:,[0, 2]])
        for n in range(data.GetNumberOfCells()):
            if data.GetCellType(n) != vtk.VTK_TRIANGLE:
                raise ValueError('{0:s} should only contain triangles')
        tri = np.ndarray((data.GetNumberOfCells(), 3))
        for n in range(data.GetNumberOfCells()):
            for nn in range(3):
                tri[n, nn] = data.GetCell(n).GetPointIds().GetId(nn)

        names = ('Slowness', 'slowness', 'Velocity', 'velocity',
                 'P-wave velocity')
        for name in names:
            if data.GetPointData().HasArray(name):
                cell_slowness = 0
                data = numpy_support.vtk_to_numpy(data.GetPointData().GetArray(name))
                break
            if data.GetCellData().HasArray(name):
                cell_slowness = 1
                data = numpy_support.vtk_to_numpy(data.GetCellData().GetArray(name))
                break
        else:
            raise ValueError('File should contain slowness or velocity data')

        if 'lowness' in name:
            slowness = data
        else:
            slowness = 1.0 / data

        m = Mesh2d(nod, tri, n_threads, cell_slowness, method, eps, maxit,
                   process_obtuse, n_secondary, n_tertiary, radius_factor_tertiary, tt_from_rp)
        m.set_slowness(slowness)
        return m


def _rebuild3d(constructor_params):
    (nodes, tetra, method, cell_slowness, n_threads, tt_from_rp, process_vel, eps,
     maxit, gradient_method, min_dist, n_secondary, n_tertiary,
     radius_factor_tertiary, translate_grid) = constructor_params

    g = Mesh3d(nodes, tetra, n_threads, cell_slowness, method, gradient_method,
               tt_from_rp, process_vel, eps, maxit, min_dist, n_secondary,
               n_tertiary, radius_factor_tertiary, translate_grid)
    return g

def _rebuild2d(constructor_params):
    (nodes, triangles, method, cell_slowness, n_threads, eps, maxit,
     process_obtuse, n_secondary, n_tertiary, radius_factor_tertiary,
     tt_from_rp) = constructor_params

    g = Mesh2d(nodes, triangles, n_threads, cell_slowness, method, eps, maxit,
        process_obtuse, n_secondary, n_tertiary, radius_factor_tertiary, tt_from_rp)
    return g

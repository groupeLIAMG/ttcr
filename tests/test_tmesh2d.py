# -*- coding: utf-8 -*-
"""Tests for verifying python wrappers, module tmesh in 2D"""

import unittest

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

import ttcrpy.tmesh as tm


def get_tt(filename):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    data = reader.GetOutput()
    names = ('Travel Time', 'Travel time', 'travel time')
    for name in names:
        if data.GetPointData().HasArray(name):
            break
    tt = vtk_to_numpy(data.GetPointData().GetArray(name))
    return tt


class TestMesh2dc(unittest.TestCase):

    def setUp(self):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName('./files/layers_fine2d.vtu')
        reader.Update()

        self.nodes = np.empty((reader.GetOutput().GetNumberOfPoints(), 2 ))
        for n in range(reader.GetOutput().GetNumberOfPoints()):
            x = reader.GetOutput().GetPoint(n)
            self.nodes[n, 0] = x[0]
            self.nodes[n, 1] = x[2]

        self.tri = np.empty((reader.GetOutput().GetNumberOfCells(), 3 ),
                            dtype=int)
        ind = vtk.vtkIdList()
        for n in range(reader.GetOutput().GetNumberOfCells()):
            reader.GetOutput().GetCellPoints(n, ind)
            self.tri[n, 0] = ind.GetId(0)
            self.tri[n, 1] = ind.GetId(1)
            self.tri[n, 2] = ind.GetId(2)

        data = reader.GetOutput()
        self.slowness = vtk_to_numpy(data.GetCellData().GetArray('Slowness'))

        self.src = np.loadtxt('./files/src2d.dat',skiprows=1)
        # we roll because file has x z t0 and we want t0 x z
        self.src = np.roll(self.src, 1).reshape((1, 3))
        self.rcv = np.loadtxt('./files/rcv2d.dat',skiprows=1)
        self.src_in = np.loadtxt('./files/src2d_in.dat',skiprows=1)
        self.src_in = np.roll(self.src_in, 1).reshape((1, 3))
        self.rcv_in = np.loadtxt('./files/rcv2d_in.dat',skiprows=1)
        
        # tm.set_verbose(True)

    def test_Mesh2Dfs(self):
        g = tm.Mesh2d(self.nodes, self.tri, method='FSM')
        tt = g.raytrace(self.src, self.rcv, slowness=self.slowness)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid2Ducfs_tt_grid.vtu')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'FSM accuracy failed (slowness in cells)')

    def test_Mesh2Dsp(self):
        g = tm.Mesh2d(self.nodes, self.tri, method='SPM', n_secondary=10)
        tt = g.raytrace(self.src, self.rcv, slowness=self.slowness)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid2Ducsp_tt_grid.vtu')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'SPM accuracy failed (slowness in cells)')

    def test_Mesh2Ddsp(self):
        g = tm.Mesh2d(self.nodes, self.tri, method='DSPM', n_secondary=3,
                      n_tertiary=3, radius_factor_tertiary=3.0)
        tt = g.raytrace(self.src, self.rcv, slowness=self.slowness)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid2Ducdsp_tt_grid.vtu')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'DSPM accuracy failed (slowness in cells)')

    def test_Mesh2Dfs_L(self):
        g = tm.Mesh2d(self.nodes, self.tri, method='FSM')
        tt, L = g.raytrace(self.src_in, self.rcv_in, slowness=self.slowness, compute_L=True)
        tt2 = L @ self.slowness
        self.assertLess(np.sum(np.abs(tt-tt2))/tt.size, 0.01,
                        'FSM_L accuracy failed (slowness in cells)')

    def test_Mesh2Dsp_L(self):
        g = tm.Mesh2d(self.nodes, self.tri, method='SPM', n_secondary=10)
        tt, L = g.raytrace(self.src_in, self.rcv_in, slowness=self.slowness, compute_L=True)
        tt2 = L @ self.slowness
        self.assertLess(np.sum(np.abs(tt-tt2))/tt.size, 0.01,
                        'SPM_L accuracy failed (slowness in cells)')

    def test_Mesh2Ddsp_L(self):
        g = tm.Mesh2d(self.nodes, self.tri, method='DSPM', n_secondary=3,
                      n_tertiary=3, radius_factor_tertiary=3.0)
        tt, L = g.raytrace(self.src_in, self.rcv_in, slowness=self.slowness, compute_L=True)
        tt2 = L @ self.slowness
        self.assertLess(np.sum(np.abs(tt-tt2))/tt.size, 0.01,
                        'DSPM_L accuracy failed (slowness in cells)')

class TestMesh2dn(unittest.TestCase):

    def setUp(self):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName('./files/gradient_fine2d.vtu')
        reader.Update()

        self.nodes = np.empty((reader.GetOutput().GetNumberOfPoints(), 2 ))
        for n in range(reader.GetOutput().GetNumberOfPoints()):
            x = reader.GetOutput().GetPoint(n)
            self.nodes[n, 0] = x[0]
            self.nodes[n, 1] = x[2]

        self.tri = np.empty((reader.GetOutput().GetNumberOfCells(), 3 ),
                            dtype=int)
        ind = vtk.vtkIdList()
        for n in range(reader.GetOutput().GetNumberOfCells()):
            reader.GetOutput().GetCellPoints(n, ind)
            self.tri[n, 0] = ind.GetId(0)
            self.tri[n, 1] = ind.GetId(1)
            self.tri[n, 2] = ind.GetId(2)

        data = reader.GetOutput()
        self.slowness = vtk_to_numpy(data.GetPointData().GetArray('Slowness'))

        self.src = np.loadtxt('./files/src2d.dat',skiprows=1)
        self.src = np.roll(self.src, 1).reshape((1, 3))
        self.rcv = np.loadtxt('./files/rcv2d.dat',skiprows=1)

    def test_Mesh2Dfs(self):
        g = tm.Mesh2d(self.nodes, self.tri, method='FSM', cell_slowness=0)
        tt = g.raytrace(self.src, self.rcv, slowness=self.slowness)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid2Dunfs_tt_grid.vtu')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'FSM accuracy failed (slowness at nodes)')

    def test_Mesh2Dsp(self):
        g = tm.Mesh2d(self.nodes, self.tri, method='SPM', n_secondary=10,
                      cell_slowness=0)
        tt = g.raytrace(self.src, self.rcv, slowness=self.slowness)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid2Dunsp_tt_grid.vtu')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'SPM accuracy failed (slowness at nodes)')

    def test_Mesh2Ddsp(self):
        g = tm.Mesh2d(self.nodes, self.tri, method='DSPM', n_secondary=3,
                      n_tertiary=3, radius_factor_tertiary=3.0, cell_slowness=0)
        tt = g.raytrace(self.src, self.rcv, slowness=self.slowness)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid2Dundsp_tt_grid.vtu')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'DSPM accuracy failed (slowness at nodes)')

class TestAniso(unittest.TestCase):

    def setUp(self):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName('./files/elliptical_fine2d.vtu')
        reader.Update()

        self.nodes = np.empty((reader.GetOutput().GetNumberOfPoints(), 2 ))
        for n in range(reader.GetOutput().GetNumberOfPoints()):
            x = reader.GetOutput().GetPoint(n)
            self.nodes[n, 0] = x[0]
            self.nodes[n, 1] = x[2]

        self.tri = np.empty((reader.GetOutput().GetNumberOfCells(), 3 ),
                            dtype=int)
        ind = vtk.vtkIdList()
        for n in range(reader.GetOutput().GetNumberOfCells()):
            reader.GetOutput().GetCellPoints(n, ind)
            self.tri[n, 0] = ind.GetId(0)
            self.tri[n, 1] = ind.GetId(1)
            self.tri[n, 2] = ind.GetId(2)

        data = reader.GetOutput()
        self.slowness = vtk_to_numpy(data.GetCellData().GetArray('Slowness'))
        self.xi = vtk_to_numpy(data.GetCellData().GetArray('xi'))

        self.src = np.loadtxt('./files/src2d.dat',skiprows=1)
        # we roll because file has x z t0 and we want t0 x z
        self.src = np.roll(self.src, 1).reshape((1, 3))
        self.rcv = np.loadtxt('./files/rcv2daniso.dat',skiprows=1)

    def test_Mesh2Dsp(self):
        g = tm.Mesh2d(self.nodes, self.tri, method='SPM', n_secondary=10, aniso='elliptical')
        g.set_slowness(self.slowness)
        g.set_xi(self.xi)
        tt = g.raytrace(self.src, self.rcv)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid2Ducsp_tt_grid_elliptical.vtu')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'SPM accuracy failed (elliptical anisotropy)')

class TestSensitivity(unittest.TestCase):
    """The matrix returned by compute_L on a triangular mesh.

    L holds one block of ncells columns per medium parameter.  Three things are
    asked of it: that it has the width the model calls for, that each block
    matches a finite difference of the traveltimes with respect to that
    parameter, and that it satisfies the homogeneity identities, which hold
    exactly and need no finite difference.

    That last one is the check the mesh could not pass while it took the length
    of every segment itself instead of asking the cells: the matrix was one of
    path lengths, so L times the slowness reproduced the traveltime of an
    isotropic medium rather than the one the call returned.
    """

    # aniso : (setter, base value, step, tolerance) per parameter
    MEDIA = {
        'iso': [('set_slowness', 0.5, 1.e-5, 5.e-3)],
        'elliptical': [('set_slowness', 0.5, 1.e-5, 5.e-3),
                       ('set_xi', 1.1, 1.e-5, 5.e-3)],
        'vti_sh': [('set_Vs0', 1.8, 1.e-5, 5.e-3),
                   ('set_gamma', 0.15, 1.e-5, 5.e-3)],
        'tilted_elliptical': [('set_slowness', 0.5, 1.e-5, 5.e-3),
                              ('set_xi', 1.1, 1.e-5, 5.e-3),
                              ('set_tilt_angle', 0.3, 1.e-5, 5.e-3)],
        'tti_sh': [('set_Vs0', 1.8, 1.e-5, 5.e-3),
                   ('set_gamma', 0.15, 1.e-5, 5.e-3),
                   ('set_tilt_angle', 0.3, 1.e-5, 5.e-3)],
        'weakly_anelliptical': [('set_slowness', 0.5, 1.e-5, 5.e-3),
                                ('set_s2', 0.05, 1.e-5, 5.e-3),
                                ('set_s4', 0.01, 1.e-5, 5.e-3)],
        # the group velocity of the coupled cells is tabulated every tenth of a
        # degree, so a finite difference has to take a much larger step, which
        # in turn costs it some accuracy
        'vti_psv': [('set_Vp0', 3.094, 0.02, 5.e-2),
                    ('set_Vs0', 1.51, 0.02, 5.e-2),
                    ('set_epsilon', 0.256, 0.02, 5.e-2),
                    ('set_delta', -0.0505, 0.02, 5.e-2)],
        'tti_psv': [('set_Vp0', 3.094, 0.02, 5.e-2),
                    ('set_Vs0', 1.51, 0.02, 5.e-2),
                    ('set_epsilon', 0.256, 0.02, 5.e-2),
                    ('set_delta', -0.0505, 0.02, 5.e-2),
                    ('set_tilt_angle', 0.3, 0.02, 5.e-2)],
    }

    # the parameter the traveltime is homogeneous in, the block holding its
    # derivative, and the degree: +1 for a slowness, -1 for a velocity
    HOMOGENEOUS = {
        'iso': (0.5, 0, 1.), 'elliptical': (0.5, 0, 1.),
        'tilted_elliptical': (0.5, 0, 1.), 'weakly_anelliptical': (0.5, 0, 1.),
        'vti_sh': (1.8, 0, -1.), 'tti_sh': (1.8, 0, -1.),
    }

    def setUp(self):
        # a triangulated square, built here so that the medium is uniform and
        # the identities above are exact
        n = 9
        xs = np.linspace(0., 1., n)
        X, Z = np.meshgrid(xs, xs, indexing='ij')
        self.nodes = np.column_stack([X.ravel(), Z.ravel()])
        tri = []
        for i in range(n-1):
            for j in range(n-1):
                a, b = i*n+j, (i+1)*n+j
                c, d = i*n+j+1, (i+1)*n+j+1
                tri.append([a, b, c])
                tri.append([b, d, c])
        self.tri = np.array(tri, dtype=np.int64)
        self.ncells = self.tri.shape[0]
        self.src = np.array([[0.15, 0.15]])
        self.rcv = np.array([[0.85, 0.8], [0.85, 0.25], [0.3, 0.85]])

    def _mesh(self, aniso, bump=None):
        """bump = (setter, value) replacing the base value of one parameter"""
        m = tm.Mesh2d(self.nodes, self.tri, method='SPM', aniso=aniso,
                      n_secondary=10, cell_slowness=1)
        for setter, base, _, _ in self.MEDIA[aniso]:
            value = base
            if bump is not None and bump[0] == setter:
                value = bump[1]
            getattr(m, setter)(np.full(self.ncells, value))
        return m

    def test_L_shape(self):
        for aniso, params in self.MEDIA.items():
            with self.subTest(aniso=aniso):
                _, L = self._mesh(aniso).raytrace(self.src, self.rcv,
                                                  compute_L=True)
                self.assertEqual(L.shape,
                                 (self.rcv.shape[0], len(params)*self.ncells))

    def test_L_against_finite_differences(self):
        for aniso, params in self.MEDIA.items():
            _, L = self._mesh(aniso).raytrace(self.src, self.rcv,
                                              compute_L=True)
            L = L.toarray()
            for blk, (setter, base, h, tol) in enumerate(params):
                with self.subTest(aniso=aniso, param=setter):
                    up = self._mesh(aniso, (setter, base+h))
                    dn = self._mesh(aniso, (setter, base-h))
                    fd = (up.raytrace(self.src, self.rcv) -
                          dn.raytrace(self.src, self.rcv))/(2*h)
                    ana = L[:, blk*self.ncells:(blk+1)*self.ncells].sum(axis=1)
                    den = np.maximum(np.abs(fd), 1.e-6)
                    self.assertLess(np.max(np.abs(ana-fd)/den), tol)

    def test_L_homogeneity(self):
        # the traveltime is homogeneous of degree one in a slowness and of
        # degree minus one in a velocity
        for aniso, (value, blk, sign) in self.HOMOGENEOUS.items():
            with self.subTest(aniso=aniso):
                tt, L = self._mesh(aniso).raytrace(self.src, self.rcv,
                                                   compute_L=True)
                L = L.toarray()
                got = L[:, blk*self.ncells:(blk+1)*self.ncells] @ \
                    np.full(self.ncells, value)
                self.assertLess(np.max(np.abs(got - sign*tt)/tt), 1.e-9)

        # the group velocity of a coupled cell scales with both of its vertical
        # velocities together
        for aniso in ('vti_psv', 'tti_psv'):
            with self.subTest(aniso=aniso):
                tt, L = self._mesh(aniso).raytrace(self.src, self.rcv,
                                                   compute_L=True)
                L = L.toarray()
                got = (L[:, :self.ncells] @ np.full(self.ncells, 3.094) +
                       L[:, self.ncells:2*self.ncells] @
                       np.full(self.ncells, 1.51))
                self.assertLess(np.max(np.abs(got + tt)/tt), 1.e-9)

    def test_L_does_not_change_traveltimes(self):
        # asking for the matrix must not change the traveltimes
        for aniso in self.MEDIA:
            with self.subTest(aniso=aniso):
                tt_plain = self._mesh(aniso).raytrace(self.src, self.rcv)
                tt_L, _ = self._mesh(aniso).raytrace(self.src, self.rcv,
                                                     compute_L=True)
                self.assertLess(np.max(np.abs(tt_L-tt_plain)), 1.e-9)


class TestPhaseAndPickle(unittest.TestCase):
    """Choosing the wave modelled, and rebuilding a mesh from a pickle.

    A mesh could not be pickled at all before: __reduce__ handed the whole
    parameter tuple to _rebuild2d as positional arguments rather than as one,
    and _rebuild2d did not take the anisotropy of the mesh from it.
    """

    def setUp(self):
        n = 7
        xs = np.linspace(0., 1., n)
        X, Z = np.meshgrid(xs, xs, indexing='ij')
        self.nodes = np.column_stack([X.ravel(), Z.ravel()])
        tri = []
        for i in range(n-1):
            for j in range(n-1):
                a, b = i*n+j, (i+1)*n+j
                c, d = i*n+j+1, (i+1)*n+j+1
                tri.append([a, b, c])
                tri.append([b, d, c])
        self.tri = np.array(tri, dtype=np.int64)
        self.ncells = self.tri.shape[0]
        self.src = np.array([[0.15, 0.15]])
        self.rcv = np.array([[0.85, 0.8], [0.85, 0.25]])

    def _mesh(self, aniso='vti_psv', phase=None):
        m = tm.Mesh2d(self.nodes, self.tri, method='SPM', aniso=aniso,
                      n_secondary=8)
        if aniso in ('vti_psv', 'tti_psv'):
            m.set_Vp0(np.full(self.ncells, 3.094))
            m.set_Vs0(np.full(self.ncells, 1.51))
            m.set_epsilon(np.full(self.ncells, 0.256))
            m.set_delta(np.full(self.ncells, -0.0505))
            if aniso == 'tti_psv':
                m.set_tilt_angle(np.full(self.ncells, 0.3))
        else:
            m.set_slowness(np.full(self.ncells, 0.5))
            m.set_xi(np.full(self.ncells, 1.1))
        if phase is not None:
            m.set_phase(phase)
        return m

    def test_qSV_is_slower_than_qP(self):
        tp = self._mesh('vti_psv', 'qP').raytrace(self.src, self.rcv)
        ts = self._mesh('vti_psv', 'qSV').raytrace(self.src, self.rcv)
        self.assertTrue(np.all(ts > tp))

    def test_sensitivity_follows_the_phase(self):
        for phase in ('qP', 'qSV'):
            with self.subTest(phase=phase):
                tt, L = self._mesh('vti_psv', phase).raytrace(
                    self.src, self.rcv, compute_L=True)
                L = L.toarray()
                got = (L[:, :self.ncells] @ np.full(self.ncells, 3.094) +
                       L[:, self.ncells:2*self.ncells] @
                       np.full(self.ncells, 1.51))
                self.assertLess(np.max(np.abs(got + tt)/tt), 1.e-9)

    def test_pickle_keeps_the_anisotropy(self):
        import pickle
        m = self._mesh('elliptical')
        tt = m.raytrace(self.src, self.rcv)
        m2 = pickle.loads(pickle.dumps(m))
        # a pickle carries the geometry and the options, not the medium
        m2.set_slowness(np.full(self.ncells, 0.5))
        m2.set_xi(np.full(self.ncells, 1.1))
        self.assertTrue(np.allclose(tt, m2.raytrace(self.src, self.rcv)))

    def test_pickle_keeps_the_phase(self):
        import pickle
        for phase in ('qP', 'qSV'):
            with self.subTest(phase=phase):
                m = self._mesh('vti_psv', phase)
                tt = m.raytrace(self.src, self.rcv)
                m2 = pickle.loads(pickle.dumps(m))
                m2.set_Vp0(np.full(self.ncells, 3.094))
                m2.set_Vs0(np.full(self.ncells, 1.51))
                m2.set_epsilon(np.full(self.ncells, 0.256))
                m2.set_delta(np.full(self.ncells, -0.0505))
                self.assertTrue(np.allclose(
                    tt, m2.raytrace(self.src, self.rcv)))


class TestWeakly(unittest.TestCase):

    def setUp(self):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName('./files/weakly_an_fine2d.vtu')
        reader.Update()

        self.nodes = np.empty((reader.GetOutput().GetNumberOfPoints(), 2 ))
        for n in range(reader.GetOutput().GetNumberOfPoints()):
            x = reader.GetOutput().GetPoint(n)
            self.nodes[n, 0] = x[0]
            self.nodes[n, 1] = x[2]

        self.tri = np.empty((reader.GetOutput().GetNumberOfCells(), 3 ),
                            dtype=int)
        ind = vtk.vtkIdList()
        for n in range(reader.GetOutput().GetNumberOfCells()):
            reader.GetOutput().GetCellPoints(n, ind)
            self.tri[n, 0] = ind.GetId(0)
            self.tri[n, 1] = ind.GetId(1)
            self.tri[n, 2] = ind.GetId(2)

        data = reader.GetOutput()
        self.slowness = 1/vtk_to_numpy(data.GetCellData().GetArray('Velocity'))
        self.s2 = vtk_to_numpy(data.GetCellData().GetArray('s2'))
        self.s4 = vtk_to_numpy(data.GetCellData().GetArray('s4'))

        self.src = np.loadtxt('./files/src2d.dat',skiprows=1)
        # we roll because file has x z t0 and we want t0 x z
        self.src = np.roll(self.src, 1).reshape((1, 3))
        self.rcv = np.loadtxt('./files/rcv2daniso.dat',skiprows=1)

    def test_Mesh2Dsp(self):
        g = tm.Mesh2d(self.nodes, self.tri, method='SPM', n_secondary=10, aniso='weakly_anelliptical')
        g.set_slowness(self.slowness)
        g.set_s2(self.s2)
        g.set_s4(self.s4)
        tt = g.raytrace(self.src, self.rcv)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid2Ducsp_tt_grid_weakly.vtu')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'SPM accuracy failed (weakly anelliptical)')

#
# class Data_kernel(unittest.TestCase):
#
#     def test_2d(self):
#
#         V = np.ones((11, 13))
#         V[:, 7:] = 2
#         slowness = 1. / V.flatten()
#
#         grx = np.arange(12.)
#         grz = np.arange(14.)
#
#         z = 0.5 + np.arange(13.)
#         Tx = np.vstack((0.5+np.zeros((13,)),
#                         z)).T
#         Rx = np.vstack((10.5+np.zeros((13,)),
#                         z)).T
#         nTx = Tx.shape[0]
#         nRx = Rx.shape[0]
#         Tx = np.kron(Tx, np.ones((nRx,1)))
#         Rx = np.kron(np.ones((nTx,1)), Rx)
#
#         L = tm.Mesh2d.data_kernel_straight_rays(Tx, Rx, grx, grz)
#         tt = L.dot(slowness)
#
#         tt2 = np.zeros(tt.shape)
#         d = np.sqrt(np.sum((Tx-Rx)**2, axis=1))
#
#         ind = np.logical_and(Tx[:,1]>7, Rx[:,1]>7)
#         tt2[ind] = d[ind]/2
#
#         ind2 = np.logical_and(Tx[:,1]<7, Rx[:,1]<7)
#         tt2[ind2] = d[ind2]
#
#         ind3 = np.logical_and(np.logical_not(ind), np.logical_not(ind2))
#
#         f = (7-Tx[ind3,1]) / (Rx[ind3,1]-Tx[ind3,1])
#         ind = (Rx[ind3,1]-Tx[ind3,1]) < 0
#         f[ind] = 1-f[ind]
#         tt2[ind3] = d[ind3]*f + d[ind3]*(1-f)/2
#
#         self.assertAlmostEqual(np.sum(np.abs(tt-tt2)), 0.0 )

if __name__ == '__main__':

    unittest.main()

# -*- coding: utf-8 -*-
"""Tests for verifying python wrappers, module rgrid in 3D"""

import unittest
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from scipy.io import mmread

vtk.vtkObject.GlobalWarningDisplayOff()

import ttcrpy.rgrid as rg


def get_tt(filename):
    reader = vtk.vtkXMLRectilinearGridReader()
    reader.SetFileName(filename)
    reader.Update()
    data = reader.GetOutput()
    names = ('Travel Time', 'Travel time', 'travel time')
    for name in names:
        if data.GetPointData().HasArray(name):
            break
    dim = data.GetDimensions()
    tt = vtk_to_numpy(data.GetPointData().GetArray(name)).reshape(dim, order='F')
    return tt.flatten()


class TestGrid3dc(unittest.TestCase):

    def setUp(self):
        reader = vtk.vtkXMLRectilinearGridReader()
        reader.SetFileName('./files/layers_medium.vtr')
        reader.Update()

        data = reader.GetOutput()
        self.x = vtk_to_numpy(data.GetXCoordinates())
        self.y = vtk_to_numpy(data.GetYCoordinates())
        self.z = vtk_to_numpy(data.GetZCoordinates())

        self.slowness = vtk_to_numpy(data.GetCellData().GetArray('Slowness'))
        dim = (self.x.size-1, self.y.size-1, self.z.size-1)
        self.slowness = self.slowness.reshape(dim, order='F').flatten()
        self.src = np.loadtxt('./files/src.dat',skiprows=1)
        self.src = self.src.reshape((1, 4))
        self.rcv = np.loadtxt('./files/rcv.dat',skiprows=1)

    def test_Grid3Dfs(self):
        g = rg.Grid3d(self.x, self.y, self.z, method='FSM', weno=1,
                      tt_from_rp=False)
        tt = g.raytrace(self.src, self.rcv, self.slowness)
        dim = (self.x.size, self.y.size, self.z.size)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid3Drcfs_tt_grid.vtr')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'FSM accuracy failed (slowness in cells)')

    def test_Grid3Dsp(self):
        g = rg.Grid3d(self.x, self.y, self.z, method='SPM', tt_from_rp=False,
                      nsnx=5, nsny=5, nsnz=5)
        tt = g.raytrace(self.src, self.rcv, self.slowness)
        dim = (self.x.size, self.y.size, self.z.size)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid3Drcsp_tt_grid.vtr')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.1,
                        'SPM accuracy failed (slowness in cells)')

    def test_Grid3Ddsp(self):
        g = rg.Grid3d(self.x, self.y, self.z, method='DSPM', tt_from_rp=False,
                      n_secondary=2, n_tertiary=3, radius_factor_tertiary=3.0)
        tt = g.raytrace(self.src, self.rcv, self.slowness)
        dim = (self.x.size, self.y.size, self.z.size)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid3Drcdsp_tt_grid.vtr')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.1,
                        'DSPM accuracy failed (slowness in cells)')


class TestGrid3dc_L(unittest.TestCase):

    def setUp(self):
        reader = vtk.vtkXMLRectilinearGridReader()
        reader.SetFileName('./files/layers_medium.vtr')
        reader.Update()

        data = reader.GetOutput()
        self.x = vtk_to_numpy(data.GetXCoordinates())
        self.y = vtk_to_numpy(data.GetYCoordinates())
        self.z = vtk_to_numpy(data.GetZCoordinates())

        self.slowness = vtk_to_numpy(data.GetCellData().GetArray('Slowness'))
        dim = (self.x.size-1, self.y.size-1, self.z.size-1)
        self.slowness = self.slowness.reshape(dim, order='F').flatten()
        self.src = np.loadtxt('./files/src3d_in.dat',skiprows=1)
        self.src = np.roll(self.src, 1).reshape((1, 4))
        self.rcv = np.loadtxt('./files/rcv3d_in.dat',skiprows=1)

    # def test_Grid3Dfs(self):
    #     g = rg.Grid3d(self.x, self.y, self.z, method='FSM', weno=1,
    #                   tt_from_rp=False)
    #     _, L = g.raytrace(self.src, self.rcv, self.slowness, compute_L=True)
    #     L2 = mmread('./files/Grid3Drcfs_L')
    #     s2 = np.loadtxt('./files/Grid3Drcfs_slo')
    #     tt = L @ self.slowness
    #     tt2 = L2 @ s2
    #     err = np.sum(np.abs(tt - tt2)) / tt.size
    #     self.assertLess(err, 0.0001, 'FSM accuracy failed (slowness in cells)')

    def test_Grid3Dsp(self):
        g = rg.Grid3d(self.x, self.y, self.z, method='SPM', tt_from_rp=False,
                      nsnx=5, nsny=5, nsnz=5)
        _, L = g.raytrace(self.src, self.rcv, self.slowness, compute_L=True)
        L2 = mmread('./files/Grid3Drcsp_L')
        s2 = np.loadtxt('./files/Grid3Drcsp_slo')
        tt = L @ self.slowness
        tt2 = L2 @ s2
        err = np.sum(np.abs(tt-tt2)) / tt.size
        self.assertLess(err, 0.0001, 'SPM accuracy failed (slowness in cells)')

    def test_Grid3Ddsp(self):
        g = rg.Grid3d(self.x, self.y, self.z, method='DSPM', tt_from_rp=False,
                      n_secondary=2, n_tertiary=3, radius_factor_tertiary=3.0)
        _, L = g.raytrace(self.src, self.rcv, self.slowness, compute_L=True)
        L2 = mmread('./files/Grid3Drcdsp_L')
        s2 = np.loadtxt('./files/Grid3Drcdsp_slo')
        tt = L @ self.slowness
        tt2 = L2 @ s2
        err = np.sum(np.abs(tt - tt2)) / tt.size
        self.assertLess(err, 0.0001, 'DSPM accuracy failed (slowness in cells)')


class TestGrid3dn(unittest.TestCase):

    def setUp(self):
        reader = vtk.vtkXMLRectilinearGridReader()
        reader.SetFileName('./files/gradient_medium.vtr')
        reader.Update()

        data = reader.GetOutput()
        self.x = vtk_to_numpy(data.GetXCoordinates())
        self.y = vtk_to_numpy(data.GetYCoordinates())
        self.z = vtk_to_numpy(data.GetZCoordinates())

        self.slowness = vtk_to_numpy(data.GetPointData().GetArray('Slowness'))
        dim = (self.x.size, self.y.size, self.z.size)
        self.slowness = self.slowness.reshape(dim, order='F').flatten()
        self.src = np.loadtxt('./files/src.dat',skiprows=1)
        self.src = self.src.reshape((1, 4))
        self.rcv = np.loadtxt('./files/rcv.dat',skiprows=1)

    def test_Grid3Dfs(self):
        g = rg.Grid3d(self.x, self.y, self.z, method='FSM', tt_from_rp=False,
                      cell_slowness=0, weno=1)
        tt = g.raytrace(self.src, self.rcv, self.slowness)
        dim = (self.x.size, self.y.size, self.z.size)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid3Drnfs_tt_grid.vtr')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'FSM accuracy failed (slowness at nodes)')

    def test_Grid3Dsp(self):
        g = rg.Grid3d(self.x, self.y, self.z, method='SPM', tt_from_rp=False,
                      nsnx=5, nsny=5, nsnz=5, cell_slowness=0)
        tt = g.raytrace(self.src, self.rcv, self.slowness)
        dim = (self.x.size, self.y.size, self.z.size)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid3Drnsp_tt_grid.vtr')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.1,
                        'SPM accuracy failed (slowness at nodes)')

    def test_Grid3Ddsp(self):
        g = rg.Grid3d(self.x, self.y, self.z, method='DSPM', tt_from_rp=False,
                      n_secondary=2, n_tertiary=3, radius_factor_tertiary=3.0,
                      cell_slowness=0)
        tt = g.raytrace(self.src, self.rcv, self.slowness)
        dim = (self.x.size, self.y.size, self.z.size)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid3Drndsp_tt_grid.vtr')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.1,
                        'DSPM accuracy failed (slowness at nodes)')


class Data_kernel(unittest.TestCase):

    def test_3d(self):

        V = np.ones((11, 12, 13))
        V[:, :, 7:] = 2
        slowness = 1. / V.flatten()

        grx = np.arange(12.)
        gry = np.arange(13.)
        grz = np.arange(14.)

        z = 0.5 + np.arange(13.)
        Tx = np.vstack((0.5+np.zeros((13,)),
                        0.5+np.zeros((13,)),
                        z)).T
        Rx = np.vstack((10.5+np.zeros((13,)),
                        11.5+np.zeros((13,)),
                        z)).T
        nTx = Tx.shape[0]
        nRx = Rx.shape[0]
        Tx = np.kron(Tx, np.ones((nRx,1)))
        Rx = np.kron(np.ones((nTx,1)), Rx)

        L = rg.Grid3d.data_kernel_straight_rays(Tx, Rx, grx, gry, grz)
        tt = L.dot(slowness)

        tt2 = np.zeros(tt.shape)
        d = np.sqrt(np.sum((Tx-Rx)**2, axis=1))

        ind = np.logical_and(Tx[:,2]>7, Rx[:,2]>7)
        tt2[ind] = d[ind]/2

        ind2 = np.logical_and(Tx[:,2]<7, Rx[:,2]<7)
        tt2[ind2] = d[ind2]

        ind3 = np.logical_and(np.logical_not(ind), np.logical_not(ind2))

        f = (7-Tx[ind3,2]) / (Rx[ind3,2]-Tx[ind3,2])
        ind = (Rx[ind3,2]-Tx[ind3,2]) < 0
        f[ind] = 1-f[ind]
        tt2[ind3] = d[ind3]*f + d[ind3]*(1-f)/2

        self.assertAlmostEqual(np.sum(np.abs(tt-tt2)), 0.0 )


if __name__ == '__main__':

    unittest.main()


class TestSensitivity3d(unittest.TestCase):
    """The matrix returned by compute_L, for every 3D anisotropy model.

    L holds one block of ncells columns per medium parameter.  Two things are
    asked of it: that each block matches a finite difference of the traveltimes
    with respect to that parameter, and that it satisfies the homogeneity
    identities, which hold exactly and need no finite difference.
    """

    # aniso : (setter, base value, step, tolerance) per parameter, in the
    # order the blocks of columns appear
    MEDIA = {
        'iso': [('set_slowness', 0.5, 1.e-5, 5.e-3)],
        'elliptical': [('set_slowness', 0.5, 1.e-5, 5.e-3),
                       ('set_chi', 1.1, 1.e-5, 5.e-3),
                       ('set_psi', 0.9, 1.e-5, 5.e-3)],
        'vti_sh': [('set_Vs0', 1.8, 1.e-5, 5.e-3),
                   ('set_gamma', 0.15, 1.e-5, 5.e-3)],
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
    }

    def setUp(self):
        n, h = 13, 0.1
        self.x = np.arange(n)*h
        self.y = np.arange(n)*h
        self.z = np.arange(n)*h
        self.ncells = (n-1)**3
        self.src = np.array([[0.3, 0.3, 0.3]])
        self.rcv = np.array([[0.9, 0.8, 0.9], [1.0, 0.4, 0.7]])

    def _grid(self, aniso, bump=None):
        """bump = (setter, value) replacing the base value of one parameter"""
        g = rg.Grid3d(self.x, self.y, self.z, method='SPM', aniso=aniso,
                      nsnx=2, nsny=2, nsnz=2, tt_from_rp=False)
        for setter, base, _, _ in self.MEDIA[aniso]:
            value = base
            if bump is not None and bump[0] == setter:
                value = bump[1]
            getattr(g, setter)(np.full(self.ncells, value))
        return g

    def test_L_shape(self):
        for aniso, params in self.MEDIA.items():
            with self.subTest(aniso=aniso):
                _, L = self._grid(aniso).raytrace(self.src, self.rcv,
                                                  compute_L=True)
                self.assertEqual(L.shape,
                                 (self.rcv.shape[0], len(params)*self.ncells))

    def test_L_against_finite_differences(self):
        for aniso, params in self.MEDIA.items():
            _, L = self._grid(aniso).raytrace(self.src, self.rcv,
                                              compute_L=True)
            L = L.toarray()
            for blk, (setter, base, h, tol) in enumerate(params):
                with self.subTest(aniso=aniso, param=setter):
                    up = self._grid(aniso, (setter, base+h))
                    dn = self._grid(aniso, (setter, base-h))
                    fd = (up.raytrace(self.src, self.rcv) -
                          dn.raytrace(self.src, self.rcv))/(2*h)
                    ana = L[:, blk*self.ncells:(blk+1)*self.ncells].sum(axis=1)
                    den = np.maximum(np.abs(fd), 1.e-6)
                    self.assertLess(np.max(np.abs(ana-fd)/den), tol)

    def test_L_homogeneity(self):
        # the traveltime is homogeneous of degree one in a slowness and of
        # degree minus one in a velocity, and the group velocity of a coupled
        # cell scales with both of its vertical velocities together.  This is
        # also the check that L is not a matrix of path lengths: it ties L back
        # to the traveltimes returned by the very same call.
        for aniso, value, sign in (('iso', 0.5, 1.),
                                   ('elliptical', 0.5, 1.),
                                   ('weakly_anelliptical', 0.5, 1.),
                                   ('vti_sh', 1.8, -1.)):
            with self.subTest(aniso=aniso):
                tt, L = self._grid(aniso).raytrace(self.src, self.rcv,
                                                   compute_L=True)
                L = L.toarray()
                got = L[:, :self.ncells] @ np.full(self.ncells, value)
                self.assertLess(np.max(np.abs(got - sign*tt)/tt), 1.e-9)

        with self.subTest(aniso='vti_psv'):
            tt, L = self._grid('vti_psv').raytrace(self.src, self.rcv,
                                                   compute_L=True)
            L = L.toarray()
            got = (L[:, :self.ncells] @ np.full(self.ncells, 3.094) +
                   L[:, self.ncells:2*self.ncells] @
                   np.full(self.ncells, 1.51))
            self.assertLess(np.max(np.abs(got + tt)/tt), 1.e-9)

    def test_L_heterogeneous(self):
        """The homogeneity identities on a medium that varies cell to cell.

        With every cell holding the same value, a column written into the wrong
        cell of its block would go unnoticed, since the row sum is unchanged.
        Making the model vary ties each column to its own cell, which is what
        checks the F-to-C reindexing the 3D grids apply when they assemble L.
        """
        rng = np.random.default_rng(12345)

        def model(aniso, scale):
            return {setter: base*(1. + scale*rng.random(self.ncells))
                    for setter, base, _, _ in self.MEDIA[aniso]}

        def grid(aniso, vals):
            g = rg.Grid3d(self.x, self.y, self.z, method='SPM', aniso=aniso,
                          nsnx=2, nsny=2, nsnz=2, tt_from_rp=False)
            for setter, v in vals.items():
                getattr(g, setter)(v)
            return g

        for aniso, sign in (('iso', 1.), ('elliptical', 1.),
                            ('weakly_anelliptical', 1.), ('vti_sh', -1.)):
            with self.subTest(aniso=aniso):
                vals = model(aniso, 0.3)
                tt, L = grid(aniso, vals).raytrace(self.src, self.rcv,
                                                   compute_L=True)
                first = self.MEDIA[aniso][0][0]
                got = L.toarray()[:, :self.ncells] @ vals[first]
                self.assertLess(np.max(np.abs(got - sign*tt)/tt), 1.e-9)

        with self.subTest(aniso='vti_psv'):
            vals = model('vti_psv', 0.2)
            tt, L = grid('vti_psv', vals).raytrace(self.src, self.rcv,
                                                   compute_L=True)
            L = L.toarray()
            got = (L[:, :self.ncells] @ vals['set_Vp0'] +
                   L[:, self.ncells:2*self.ncells] @ vals['set_Vs0'])
            self.assertLess(np.max(np.abs(got + tt)/tt), 1.e-9)

    def test_azimuthal_symmetry(self):
        """A transversely isotropic medium is symmetric about the vertical
        axis, and a cubic grid maps onto itself under a quarter turn, so
        rotating the shot and the receivers by 90 degrees about the centre must
        reproduce the traveltimes exactly.  This is the test of the reduction
        the 3D coupled cells make, from the segment to its polar angle.
        """
        c = 0.5*(self.x[0] + self.x[-1])

        def turn(p):
            q = p.copy()
            q[:, 0] = c + (p[:, 1] - c)
            q[:, 1] = c - (p[:, 0] - c)
            return q

        for aniso in ('vti_psv', 'vti_sh', 'weakly_anelliptical'):
            with self.subTest(aniso=aniso):
                g = self._grid(aniso)
                tt = g.raytrace(self.src, self.rcv)
                tt_turned = g.raytrace(turn(self.src), turn(self.rcv))
                np.testing.assert_allclose(tt, tt_turned, rtol=1.e-12)

    def test_phase_and_pickle(self):
        import pickle
        g = self._grid('vti_psv')
        tt_qp = g.raytrace(self.src, self.rcv)
        g.set_phase('qSV')
        tt_qsv = g.raytrace(self.src, self.rcv)
        # the shear wave is the slower one here
        self.assertTrue(np.all(tt_qsv > tt_qp))

        g2 = pickle.loads(pickle.dumps(g))
        for setter, base, _, _ in self.MEDIA['vti_psv']:
            getattr(g2, setter)(np.full(self.ncells, base))
        np.testing.assert_allclose(g2.raytrace(self.src, self.rcv), tt_qsv,
                                   rtol=1.e-12)

    def test_refused_where_not_implemented(self):
        with self.assertRaises(ValueError):
            rg.Grid3d(self.x, self.y, self.z, method='FSM', aniso='vti_psv')
        with self.assertRaises(ValueError):
            rg.Grid3d(self.x, self.y, self.z, method='SPM', aniso='vti_psv',
                      cell_slowness=0)
        with self.assertRaises(ValueError):
            rg.Grid3d(self.x, self.y, self.z, method='SPM', aniso='tti_psv')

    def test_float_grid(self):
        """The single-precision class takes the same models."""
        g = rg.Grid3d(self.x, self.y, self.z, method='SPM', aniso='vti_psv',
                      nsnx=2, nsny=2, nsnz=2, tt_from_rp=False,
                      dtype=np.float32)
        for setter, base, _, _ in self.MEDIA['vti_psv']:
            getattr(g, setter)(np.full(self.ncells, base, dtype=np.float32))
        tt, L = g.raytrace(self.src.astype(np.float32),
                           self.rcv.astype(np.float32), compute_L=True)
        self.assertEqual(L.shape, (self.rcv.shape[0], 4*self.ncells))
        L = L.toarray()
        got = (L[:, :self.ncells] @ np.full(self.ncells, 3.094) +
               L[:, self.ncells:2*self.ncells] @ np.full(self.ncells, 1.51))
        self.assertLess(np.max(np.abs(got + tt)/tt), 1.e-4)

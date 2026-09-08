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


class TestComputeH(unittest.TestCase):
    """The Jacobian returned by compute_H.

    H holds the derivatives of the arrival time with respect to the hypocentre
    parameters.  The spatial columns are -s(x_s) times the unit take-off
    direction, so they are checked against a model whose take-off is known
    exactly: a constant velocity gradient, for which

        T = arccosh[1 + a^2 r^2 / (2 v(z_s) v(z))] / a

    and dT/dx_s follows by differentiating that closed form.
    """

    V0 = 2.0
    A = 1.5
    h = 0.05
    n = 41

    def setUp(self):
        self.x = np.arange(self.n) * self.h
        self.y = self.x.copy()
        self.z = self.x.copy()
        X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')
        self.V = self.V0 + self.A * Z
        self.Vc = self.V0 + self.A * (Z[:-1, :-1, :-1] + 0.5 * self.h)
        self.xs = np.array([1.0, 1.0, 1.5])
        self.rcv = np.array([[0.30, 0.55, 0.10],
                             [1.70, 0.40, 0.15],
                             [0.60, 1.75, 0.05],
                             [1.80, 1.80, 0.20]])
        self.src = np.repeat(self.xs.reshape(1, 3), self.rcv.shape[0], axis=0)

    def _T_analytic(self, p, s):
        v1 = self.V0 + self.A * s[2]
        v2 = self.V0 + self.A * p[..., 2]
        r2 = np.sum((p - s) ** 2, axis=-1)
        return np.arccosh(1.0 + self.A * self.A * r2 / (2 * v1 * v2)) / self.A

    def _dT_analytic(self):
        """dT/dx_s by central differences of the closed form."""
        d = 1.e-6
        out = np.zeros((self.rcv.shape[0], 3))
        for k in range(self.rcv.shape[0]):
            for j in range(3):
                a = self.xs.copy()
                b = self.xs.copy()
                a[j] += d
                b[j] -= d
                out[k, j] = (self._T_analytic(self.rcv[k], a) -
                             self._T_analytic(self.rcv[k], b)) / (2 * d)
        return out

    def _grid(self, cell_slowness, **kwargs):
        g = rg.Grid3d(self.x, self.y, self.z, n_threads=1,
                      cell_slowness=cell_slowness, **kwargs)
        g.set_slowness((1.0 / (self.Vc if cell_slowness else self.V)).ravel())
        return g

    def _configs(self):
        return (('node FSM', 0, dict(method='FSM')),
                ('node SPM', 0, dict(method='SPM', nsnx=3, nsny=3, nsnz=3)),
                ('node DSPM', 0, dict(method='DSPM', n_secondary=3,
                                      n_tertiary=3)),
                ('cell FSM', 1, dict(method='FSM')),
                ('cell SPM', 1, dict(method='SPM', nsnx=3, nsny=3, nsnz=3)))

    def test_against_analytic(self):
        """Direction and magnitude of the spatial columns."""
        dT = self._dT_analytic()
        for label, cs, kwargs in self._configs():
            with self.subTest(config=label):
                tt, H = self._grid(cs, **kwargs).compute_H(self.src, self.rcv)
                for k in range(self.rcv.shape[0]):
                    hk = H[k, 1:]
                    ang = np.degrees(np.arccos(np.clip(
                        hk @ dT[k] / (np.linalg.norm(hk) *
                                      np.linalg.norm(dT[k])), -1, 1)))
                    self.assertLess(ang, 2.5,
                                    '%s rcv %d: take-off off by %.2f deg'
                                    % (label, k, ang))
                    rel = abs(np.linalg.norm(hk) - np.linalg.norm(dT[k])) / \
                        np.linalg.norm(dT[k])
                    self.assertLess(rel, 0.02,
                                    '%s rcv %d: |H| off by %.1f%%'
                                    % (label, k, 100 * rel))

    def test_against_finite_difference(self):
        """H against a finite difference of the solver itself.

        Independent of the analytic solution, so it catches sign and
        column-order errors that the analytic test could absorb.
        """
        g = self._grid(0, method='FSM')
        delta = 4 * self.h
        _, H = g.compute_H(self.src, self.rcv)
        for j in range(3):
            sp = self.xs.copy()
            sm = self.xs.copy()
            sp[j] += delta
            sm[j] -= delta
            tp = g.raytrace(np.repeat(sp.reshape(1, 3), 4, axis=0), self.rcv)
            tm = g.raytrace(np.repeat(sm.reshape(1, 3), 4, axis=0), self.rcv)
            fd = (tp - tm) / (2 * delta)
            for k in range(self.rcv.shape[0]):
                self.assertAlmostEqual(H[k, 1 + j], fd[k], places=2,
                                       msg='column %d, rcv %d' % (1 + j, k))

    def test_columns(self):
        """Shape and column convention of both forms."""
        g = self._grid(0, method='FSM')
        tt, H = g.compute_H(self.src, self.rcv)
        self.assertEqual(H.shape, (self.rcv.shape[0], 4))
        np.testing.assert_allclose(H[:, 0], 1.0)

        tt2, H2 = g.compute_H(self.src, self.rcv, full=False)
        self.assertEqual(H2.shape, (self.rcv.shape[0], 2))
        np.testing.assert_allclose(H2, H[:, 1:3])
        np.testing.assert_allclose(tt2, tt)

    def test_traveltimes_match_raytrace(self):
        """compute_H raytraces itself; those traveltimes must be the usual ones."""
        for label, cs, kwargs in self._configs():
            with self.subTest(config=label):
                g = self._grid(cs, **kwargs)
                tt, H = g.compute_H(self.src, self.rcv)
                np.testing.assert_allclose(tt, g.raytrace(self.src, self.rcv))

    def test_method_independence(self):
        """H must not depend on which solver produced the field.

        Regression test for the raypath endpoint convention, which differs
        between SPM and FSM/DSPM: building H from rays[1] instead of the source
        end fails here by tens of degrees.
        """
        ref = None
        for label, cs, kwargs in self._configs():
            _, H = self._grid(cs, **kwargs).compute_H(self.src, self.rcv)
            if ref is None:
                ref = H
                continue
            for k in range(self.rcv.shape[0]):
                ang = np.degrees(np.arccos(np.clip(
                    H[k, 1:] @ ref[k, 1:] / (np.linalg.norm(H[k, 1:]) *
                                             np.linalg.norm(ref[k, 1:])),
                    -1, 1)))
                self.assertLess(ang, 3.0,
                                '%s rcv %d differs from reference by %.2f deg'
                                % (label, k, ang))

    def test_translated_grid(self):
        """translate_grid must not change H.

        raytrace() shifts coordinates by the grid origin internally, so the
        traveltime field lives in translated space; compute_H has to walk it in
        the same frame.  Checked at UTM-like offsets, where getting this wrong
        moves the evaluation point by hundreds of kilometres.
        """
        off = np.array([500000.0, 4800000.0, 0.0])
        xo = self.x + off[0]
        yo = self.y + off[1]
        zo = self.z + off[2]
        src = self.src + off
        rcv = self.rcv + off

        out = {}
        for tg in (False, True):
            g = rg.Grid3d(xo, yo, zo, n_threads=1, method='FSM',
                          cell_slowness=0, translate_grid=tg)
            g.set_slowness((1.0 / self.V).ravel())
            out[tg] = g.compute_H(src, rcv)

        np.testing.assert_allclose(out[True][0], out[False][0])
        np.testing.assert_allclose(out[True][1], out[False][1], atol=1e-6)


class TestRaypathOrder(unittest.TestCase):
    """Raypath coordinates run from source to receiver, whatever the solver.

    The shortest-path solvers have always returned that order (they reorder
    explicitly, "the order should be from Tx to Rx").  The steepest-descent
    raypath builders used by FSM and DSPM walked from the receiver down to the
    source and returned the result as built, so the convention depended on the
    method.  Code indexing one end silently got the other.
    """

    h = 0.05

    def setUp(self):
        self.x = np.arange(21) * self.h
        self.y = self.x.copy()
        self.z = self.x.copy()
        X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')
        self.V = 2.0 + 3.0 * Z
        self.Vc = 2.0 + 3.0 * (Z[:-1, :-1, :-1] + 0.5 * self.h)
        self.src = np.array([[0.5, 0.5, 0.8]])
        self.rcv = np.array([[0.9, 0.5, 0.0]])

    def test_source_first(self):
        for method, kwargs in (('FSM', {}),
                               ('SPM', dict(nsnx=3, nsny=3, nsnz=3)),
                               ('DSPM', dict(n_secondary=3, n_tertiary=3))):
            for cell_slowness in (0, 1):
                with self.subTest(method=method, cell_slowness=cell_slowness):
                    g = rg.Grid3d(self.x, self.y, self.z, n_threads=1,
                                  method=method, cell_slowness=cell_slowness,
                                  **kwargs)
                    g.set_slowness((1.0 / (self.Vc if cell_slowness
                                           else self.V)).ravel())
                    r = np.asarray(g.raytrace(self.src, self.rcv,
                                              return_rays=True)[1][0])
                    self.assertGreater(r.shape[0], 2)
                    # the path runs from the source to the receiver, and both
                    # endpoints are present exactly
                    np.testing.assert_allclose(r[0], self.src[0], atol=1e-9)
                    np.testing.assert_allclose(r[-1], self.rcv[0], atol=1e-9)


class TestSourceOffNode(unittest.TestCase):
    """A source that does not coincide with a grid node.

    The shortest-path solvers add an extra node for such a source and number it
    nodes.size()+txNodes.size(); the walk back from the receiver recognises it
    by that number being >= nodes.size(), which is what ends the walk.
    Grid3Drcsp numbered it one too low, so the walk read the source as an
    ordinary grid node, never reached the end of the parent chain and cycled,
    growing l_data and the raypath on every turn until the process ran out of
    memory.

    A source landing exactly on a node takes the other branch of initQueue and
    creates no extra node, which is why TestRaypathOrder -- whose source sits
    on a node -- stayed green throughout.  Every coordinate below is chosen to
    miss the primary nodes and the secondary ones between them.

    The two checks are deliberately named so that unittest, which runs the
    methods of a class in alphabetical order, reaches compute_L first: should
    the cycle return, l_data stays bounded (accumulate() merges entries by cell
    index) and the run hangs, whereas the raypath grows without limit and takes
    the machine down with it -- measured at 4 GB within three seconds.  Keep
    test_compute_L_off_node sorting ahead of test_raypath_off_node.
    """

    h = 0.05
    src = np.array([[0.517, 0.483, 0.762]])
    rcv = np.array([[0.913, 0.526, 0.038],
                    [0.237, 0.688, 0.114]])

    def setUp(self):
        self.x = np.arange(21) * self.h
        self.y = self.x.copy()
        self.z = self.x.copy()
        _, _, Z = np.meshgrid(self.x, self.y, self.z, indexing='ij')
        self.V = 2.0 + 3.0 * Z
        self.Vc = 2.0 + 3.0 * (Z[:-1, :-1, :-1] + 0.5 * self.h)

    def _grid(self, method, cell_slowness, **kwargs):
        g = rg.Grid3d(self.x, self.y, self.z, n_threads=1, method=method,
                      cell_slowness=cell_slowness, **kwargs)
        s = (1.0 / (self.Vc if cell_slowness else self.V)).ravel()
        g.set_slowness(s)
        return g, s

    def test_compute_L_off_node(self):
        # L holds the length the ray spends in each cell, so L @ s is the
        # traveltime; a truncated or looping walk cannot satisfy that.
        for method, kwargs in (('SPM', dict(nsnx=3, nsny=3, nsnz=3)),
                               ('DSPM', dict(n_secondary=3, n_tertiary=3))):
            with self.subTest(method=method):
                g, s = self._grid(method, 1, **kwargs)
                src = np.repeat(self.src, self.rcv.shape[0], axis=0)
                tt, L = g.raytrace(src, self.rcv, compute_L=True)
                self.assertEqual(L.shape, (self.rcv.shape[0], s.size))
                self.assertGreater(L.nnz, 0)
                np.testing.assert_allclose(L @ s, tt, rtol=1e-9)

    def test_raypath_off_node(self):
        for method, kwargs in (('FSM', {}),
                               ('SPM', dict(nsnx=3, nsny=3, nsnz=3)),
                               ('DSPM', dict(n_secondary=3, n_tertiary=3))):
            for cell_slowness in (0, 1):
                with self.subTest(method=method, cell_slowness=cell_slowness):
                    g, _ = self._grid(method, cell_slowness, **kwargs)
                    src = np.repeat(self.src, self.rcv.shape[0], axis=0)
                    _, rays = g.raytrace(src, self.rcv, return_rays=True)
                    self.assertEqual(len(rays), self.rcv.shape[0])
                    for n, ray in enumerate(rays):
                        r = np.asarray(ray)
                        self.assertGreater(r.shape[0], 2)
                        # a finite walk that starts at the source it was given
                        # and stops at the receiver
                        np.testing.assert_allclose(r[0], self.src[0], atol=1e-9)
                        np.testing.assert_allclose(r[-1], self.rcv[n], atol=1e-9)


class TestUniformSpacing3d(unittest.TestCase):
    """Node spacing must be constant along each axis, but may differ between
    axes.  Every method takes the spacing as ``x[1] - x[0]`` and assumes it
    holds throughout, so an unevenly spaced axis would be solved on a grid
    other than the one asked for.
    """

    methods = ('FSM', 'SPM', 'DSPM')
    dtypes = (np.float64, np.float32)

    def _uniform(self, dtype):
        """dx, dy and dz all different, each axis evenly spaced."""
        return (np.arange(0.0, 11.0, 1.00, dtype=dtype),
                np.arange(0.0, 22.0, 2.00, dtype=dtype),
                np.arange(0.0, 5.50, 0.50, dtype=dtype))

    def test_unequal_but_uniform_accepted(self):
        """dx != dy != dz is fine for every method, in either precision."""
        for dtype in self.dtypes:
            x, y, z = self._uniform(dtype)
            for method in self.methods:
                with self.subTest(method=method, dtype=np.dtype(dtype).name):
                    g = rg.Grid3d(x, y, z, method=method, n_threads=1,
                                  dtype=dtype)
                    # each axis keeps its own spacing
                    self.assertAlmostEqual(g.dx, 1.0, places=5)
                    self.assertAlmostEqual(g.dy, 2.0, places=5)
                    self.assertAlmostEqual(g.dz, 0.5, places=5)

    def test_non_uniform_rejected(self):
        """Each axis is checked, and the message names the offending one."""
        for dtype in self.dtypes:
            x, y, z = self._uniform(dtype)
            uneven = np.array([0.0, 1.0, 2.0, 4.0, 5.0, 6.0], dtype=dtype)
            for axis, args in (('x', (uneven, y, z)),
                               ('y', (x, uneven, z)),
                               ('z', (x, y, uneven))):
                for method in self.methods:
                    with self.subTest(axis=axis, method=method,
                                      dtype=np.dtype(dtype).name):
                        with self.assertRaises(ValueError) as cm:
                            rg.Grid3d(*args, method=method, n_threads=1,
                                      dtype=dtype)
                        self.assertIn('uniformly spaced', str(cm.exception))
                        self.assertIn('along ' + axis, str(cm.exception))

    def test_small_departure_from_uniform_rejected(self):
        """A single displaced node is caught, not just an obvious gap."""
        x, y, z = self._uniform(np.float64)
        nudged = np.arange(0.0, 11.0)
        nudged[5] += 0.05
        with self.assertRaises(ValueError):
            rg.Grid3d(nudged, y, z, method='FSM', n_threads=1)

    def test_float32_uniform_not_rejected(self):
        """The tolerance must follow the dtype, not be a fixed absolute value.

        A float32 axis cannot be exactly uniform: with a step of 0.05 its
        node-to-node deviation reaches 3e-6, so a fixed tolerance of 1e-6
        would reject a perfectly even axis.
        """
        fine = np.arange(0.0, 50.0, 0.05, dtype=np.float32)
        coarse = np.arange(0.0, 11.0, 1.0, dtype=np.float32)
        rg.Grid3d(fine, coarse, coarse, method='FSM', n_threads=1,
                  dtype=np.float32)

    def test_single_cell_axis_accepted(self):
        """Two nodes are one interval, uniform by construction."""
        two = np.array([0.0, 1.0])
        rest = np.arange(0.0, 11.0)
        rg.Grid3d(two, rest, rest, method='FSM', n_threads=1)

    def test_single_precision_values_in_double_array(self):
        """Coordinates that came from a single-precision source.

        Reading a float32 VTK file gives a float64 array whose values still
        carry only float32 precision, so its steps wobble at the 1e-6 level.
        Judging those by float64 epsilon would reject an axis that is as
        uniform as its source allows -- one of this project's own fixtures,
        Grid2Drcfs_OpenCL_tt_grid.vtr, has steps from 0.1999988 to 0.2000008.
        The tolerance therefore follows the precision of the values, not the
        dtype of the array holding them.
        """
        x = (np.arange(0.0, 20.2, 0.2, dtype=np.float32)).astype(np.float64)
        self.assertEqual(x.dtype, np.float64)
        self.assertGreater(np.ptp(np.diff(x)), 1e-7)   # genuinely wobbly
        rest = np.arange(0.0, 11.0)
        rg.Grid3d(x, rest, rest, method='FSM', n_threads=1)

    def test_all_bundled_grid_fixtures_accepted(self):
        """No .vtr shipped with the tests may trip the check."""
        import glob
        import vtk
        from vtk.util.numpy_support import vtk_to_numpy
        checked = 0
        for fname in sorted(glob.glob('files/*.vtr')):
            reader = vtk.vtkXMLRectilinearGridReader()
            reader.SetFileName(fname)
            reader.Update()
            data = reader.GetOutput()
            if data is None or data.GetXCoordinates() is None:
                continue
            for name, arr in (('x', data.GetXCoordinates()),
                              ('y', data.GetYCoordinates()),
                              ('z', data.GetZCoordinates())):
                a = vtk_to_numpy(arr)
                if a.size < 3:
                    continue
                with self.subTest(file=fname, axis=name):
                    rg._check_uniform_spacing(a, name)
                checked += 1
        self.assertGreater(checked, 20)

    def test_traveltimes_correct_on_unequal_spacing(self):
        """Guards the spacing being honoured, not merely accepted.

        Constant slowness, receivers along each axis, so the exact answer is
        s*r.  A grid that quietly used dx for all three axes would miss on the
        y and z receivers.
        """
        s = 0.4
        x = np.arange(0.0, 31.0, 1.0)
        y = np.arange(0.0, 62.0, 2.0)
        z = np.arange(0.0, 15.5, 0.5)
        src = np.array([[15.0, 30.0, 7.5]])
        rcv = np.array([[15.0 + 4.0, 30.0, 7.5],
                        [15.0, 30.0 + 8.0, 7.5],
                        [15.0, 30.0, 7.5 + 2.0]])
        exact = s * np.linalg.norm(rcv - src[0], axis=1)
        for weno in (0, 1):
            with self.subTest(weno=weno):
                g = rg.Grid3d(x, y, z, method='FSM', cell_slowness=0,
                              tt_from_rp=0, weno=weno, n_threads=1)
                tt = g.raytrace(src, rcv,
                                slowness=np.full(x.size * y.size * z.size, s))
                np.testing.assert_allclose(tt, exact, rtol=1e-6)

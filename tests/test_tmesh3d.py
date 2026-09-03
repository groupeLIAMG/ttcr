# -*- coding: utf-8 -*-
"""Tests for verifying python wrappers, module tmesh in 3D"""

import unittest

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from scipy.io import mmread

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


class TestMesh3Dc(unittest.TestCase):

    def setUp(self):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName('./files/layers_medium.vtu')
        reader.Update()

        self.nodes = np.empty((reader.GetOutput().GetNumberOfPoints(), 3 ))
        for n in range(reader.GetOutput().GetNumberOfPoints()):
            x = reader.GetOutput().GetPoint(n)
            self.nodes[n, 0] = x[0]
            self.nodes[n, 1] = x[1]
            self.nodes[n, 2] = x[2]

        self.tet = np.empty((reader.GetOutput().GetNumberOfCells(), 4 ), dtype=int)
        ind = vtk.vtkIdList()
        for n in range(reader.GetOutput().GetNumberOfCells()):
            reader.GetOutput().GetCellPoints(n, ind)
            self.tet[n, 0] = ind.GetId(0)
            self.tet[n, 1] = ind.GetId(1)
            self.tet[n, 2] = ind.GetId(2)
            self.tet[n, 3] = ind.GetId(3)

        data = reader.GetOutput()
        self.slowness = vtk_to_numpy(data.GetCellData().GetArray('Slowness'))

        self.src = np.loadtxt('./files/src.dat',skiprows=1)
        self.src = np.roll(self.src, 1).reshape((1, 4))
        self.rcv = np.loadtxt('./files/rcv.dat',skiprows=1)

    def test_Mesh3Dfs(self):
        g = tm.Mesh3d(self.nodes, self.tet, method='FSM', tt_from_rp=0)
        tt = g.raytrace(self.src, self.rcv, slowness=self.slowness)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid3Ducfs_tt_grid.vtu')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'FSM accuracy failed (slowness in cells)')

    def test_Mesh3Dsp(self):
        g = tm.Mesh3d(self.nodes, self.tet, method='SPM',
                      n_secondary=5, tt_from_rp=0)
        tt = g.raytrace(self.src, self.rcv, slowness=self.slowness)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid3Ducsp_tt_grid.vtu')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'SPM accuracy failed (slowness in cells)')

    def test_Mesh3Ddsp(self):
        g = tm.Mesh3d(self.nodes, self.tet, method='DSPM', n_secondary=2,
                      n_tertiary=3, radius_factor_tertiary=3.0, tt_from_rp=0)
        tt = g.raytrace(self.src, self.rcv, slowness=self.slowness)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid3Ducdsp_tt_grid.vtu')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'DSPM accuracy failed (slowness in cells)')


class TestMesh3Dc_L(unittest.TestCase):

    def setUp(self):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName('./files/layers_medium.vtu')
        reader.Update()

        self.nodes = np.empty((reader.GetOutput().GetNumberOfPoints(), 3 ))
        for n in range(reader.GetOutput().GetNumberOfPoints()):
            x = reader.GetOutput().GetPoint(n)
            self.nodes[n, 0] = x[0]
            self.nodes[n, 1] = x[1]
            self.nodes[n, 2] = x[2]

        self.tet = np.empty((reader.GetOutput().GetNumberOfCells(), 4 ), dtype=int)
        ind = vtk.vtkIdList()
        for n in range(reader.GetOutput().GetNumberOfCells()):
            reader.GetOutput().GetCellPoints(n, ind)
            self.tet[n, 0] = ind.GetId(0)
            self.tet[n, 1] = ind.GetId(1)
            self.tet[n, 2] = ind.GetId(2)
            self.tet[n, 3] = ind.GetId(3)

        data = reader.GetOutput()
        self.slowness = vtk_to_numpy(data.GetCellData().GetArray('Slowness'))

        self.src = np.loadtxt('./files/src3d_in.dat', skiprows=1)
        self.src = np.roll(self.src, 1).reshape((1, 4))
        self.rcv = np.loadtxt('./files/rcv3d_in.dat', skiprows=1)

    def test_Mesh3Dfs(self):
        g = tm.Mesh3d(self.nodes, self.tet, method='FSM', tt_from_rp=0)
        _, L = g.raytrace(self.src, self.rcv, slowness=self.slowness, compute_L=True)
        L2 = mmread('./files/Grid3Ducfs_L')
        dL = L - L2
        err = np.sum(np.abs(dL))
        self.assertLess(err, 0.0001, 'FSM (L) accuracy failed (slowness in cells)')

    def test_Mesh3Dsp(self):
        g = tm.Mesh3d(self.nodes, self.tet, method='SPM', n_secondary=5, tt_from_rp=0)
        _, L = g.raytrace(self.src, self.rcv, slowness=self.slowness, compute_L=True)
        L2 = mmread('./files/Grid3Ducsp_L')
        dL = L - L2
        err = np.sum(np.abs(dL))
        self.assertLess(err, 0.0001, 'SPM (L) accuracy failed (slowness in cells)')

    def test_Mesh3Ddsp(self):
        g = tm.Mesh3d(self.nodes, self.tet, method='DSPM', n_secondary=2,
                      n_tertiary=3, radius_factor_tertiary=3.0, tt_from_rp=0)
        _, L = g.raytrace(self.src, self.rcv, slowness=self.slowness, compute_L=True)
        L2 = mmread('./files/Grid3Ducdsp_L')
        dL = L - L2
        err = np.sum(np.abs(dL))
        self.assertLess(err, 0.0001, 'DSPM (L) accuracy failed (slowness in cells)')


class TestMesh3Dn(unittest.TestCase):

    def setUp(self):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName('./files/gradient_medium.vtu')
        reader.Update()

        self.nodes = np.empty((reader.GetOutput().GetNumberOfPoints(), 3 ))
        for n in range(reader.GetOutput().GetNumberOfPoints()):
            x = reader.GetOutput().GetPoint(n)
            self.nodes[n, 0] = x[0]
            self.nodes[n, 1] = x[1]
            self.nodes[n, 2] = x[2]

        self.tet = np.empty((reader.GetOutput().GetNumberOfCells(), 4 ), dtype=int)
        ind = vtk.vtkIdList()
        for n in range(reader.GetOutput().GetNumberOfCells()):
            reader.GetOutput().GetCellPoints(n, ind)
            self.tet[n, 0] = ind.GetId(0)
            self.tet[n, 1] = ind.GetId(1)
            self.tet[n, 2] = ind.GetId(2)
            self.tet[n, 3] = ind.GetId(3)

        data = reader.GetOutput()
        self.slowness = vtk_to_numpy(data.GetPointData().GetArray('Slowness'))

        self.src = np.loadtxt('./files/src.dat',skiprows=1)
        self.src = np.roll(self.src, 1).reshape((1, 4))
        self.rcv = np.loadtxt('./files/rcv.dat',skiprows=1)

    def test_Mesh3Dfs(self):
        g = tm.Mesh3d(self.nodes, self.tet, method='FSM', cell_slowness=0,
                      tt_from_rp=0)
        tt = g.raytrace(self.src, self.rcv, slowness=self.slowness)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid3Dunfs_tt_grid.vtu')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'FSM accuracy failed (slowness at nodes)')

    def test_Mesh3Dsp(self):
        g = tm.Mesh3d(self.nodes, self.tet, method='SPM', n_secondary=5,
                      cell_slowness=0, tt_from_rp=0)
        tt = g.raytrace(self.src, self.rcv, slowness=self.slowness)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid3Dunsp_tt_grid.vtu')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'SPM accuracy failed (slowness at nodes)')

    def test_Mesh3Ddsp(self):
        g = tm.Mesh3d(self.nodes, self.tet, method='DSPM', n_secondary=2,
                      n_tertiary=3, radius_factor_tertiary=3.0, cell_slowness=0,
                      tt_from_rp=0)
        tt = g.raytrace(self.src, self.rcv, slowness=self.slowness)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid3Dundsp_tt_grid.vtu')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'DSPM accuracy failed (slowness at nodes)')

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


def _unit_cube_mesh(n=5, h=0.25):
    """A conforming tetrahedral mesh of a cube, n nodes per side.

    Each cube is split the Freudenthal way, into the six tetrahedra sharing the
    (0,0,0)-(1,1,1) diagonal, which keeps the mesh conforming across faces.
    """
    idx = lambda i, j, k: (i*n + j)*n + k
    nodes = np.empty(((n)**3, 3))
    for i in range(n):
        for j in range(n):
            for k in range(n):
                nodes[idx(i, j, k)] = (i*h, j*h, k*h)

    steps = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
    tet = []
    for i in range(n-1):
        for j in range(n-1):
            for k in range(n-1):
                base = np.array([i, j, k])
                for perm in steps:
                    c = base.copy()
                    verts = [idx(*c)]
                    for ax in perm:
                        c = c.copy()
                        c[ax] += 1
                        verts.append(idx(*c))
                    # keep a positive volume, as the solver expects
                    p = nodes[verts]
                    if np.dot(np.cross(p[1]-p[0], p[2]-p[0]), p[3]-p[0]) < 0:
                        verts[1], verts[2] = verts[2], verts[1]
                    tet.append(verts)
    return nodes, np.array(tet, dtype=int)


class TestSensitivity3d(unittest.TestCase):
    """The matrix returned by compute_L, for every 3D anisotropy model.

    L holds one block of ncells columns per medium parameter.  Two things are
    asked of it: that each block matches a finite difference of the traveltimes
    with respect to that parameter, and that it satisfies the homogeneity
    identities, which hold exactly and need no finite difference.
    """

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
        self.nodes, self.tet = _unit_cube_mesh(5, 0.25)
        self.ncells = self.tet.shape[0]
        self.src = np.array([[0.15, 0.15, 0.15]])
        self.rcv = np.array([[0.85, 0.7, 0.85], [0.9, 0.3, 0.6]])

    def _mesh(self, aniso, bump=None, values=None):
        g = tm.Mesh3d(self.nodes, self.tet, method='SPM', aniso=aniso,
                      n_secondary=2, tt_from_rp=False)
        for setter, base, _, _ in self.MEDIA[aniso]:
            if values is not None:
                v = values[setter]
            else:
                v = base
                if bump is not None and bump[0] == setter:
                    v = bump[1]
                v = np.full(self.ncells, v)
            getattr(g, setter)(v)
        return g

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
        """L times the medium must reproduce the traveltimes of the same call.

        The traveltime is homogeneous of degree one in a slowness and of degree
        minus one in a velocity, and the group velocity of a coupled cell
        scales with both of its vertical velocities together.  This is also the
        check that L is not a matrix of path lengths.
        """
        for aniso, value, sign in (('iso', 0.5, 1.),
                                   ('elliptical', 0.5, 1.),
                                   ('weakly_anelliptical', 0.5, 1.),
                                   ('vti_sh', 1.8, -1.)):
            with self.subTest(aniso=aniso):
                tt, L = self._mesh(aniso).raytrace(self.src, self.rcv,
                                                   compute_L=True)
                L = L.toarray()
                got = L[:, :self.ncells] @ np.full(self.ncells, value)
                self.assertLess(np.max(np.abs(got - sign*tt)/tt), 1.e-9)

        with self.subTest(aniso='vti_psv'):
            tt, L = self._mesh('vti_psv').raytrace(self.src, self.rcv,
                                                   compute_L=True)
            L = L.toarray()
            got = (L[:, :self.ncells] @ np.full(self.ncells, 3.094) +
                   L[:, self.ncells:2*self.ncells] @ np.full(self.ncells, 1.51))
            self.assertLess(np.max(np.abs(got + tt)/tt), 1.e-9)

    def test_L_heterogeneous(self):
        """The same identities on a medium that varies tetrahedron to tetrahedron.

        With every cell holding the same value, a column written into the wrong
        cell of its block would go unnoticed, since the row sum is unchanged.
        Making the model vary ties each column to its own cell.
        """
        rng = np.random.default_rng(20240826)

        def model(aniso, scale):
            return {s: b*(1. + scale*rng.random(self.ncells))
                    for s, b, _, _ in self.MEDIA[aniso]}

        for aniso, sign in (('iso', 1.), ('elliptical', 1.),
                            ('weakly_anelliptical', 1.), ('vti_sh', -1.)):
            with self.subTest(aniso=aniso):
                vals = model(aniso, 0.3)
                tt, L = self._mesh(aniso, values=vals).raytrace(
                    self.src, self.rcv, compute_L=True)
                first = self.MEDIA[aniso][0][0]
                got = L.toarray()[:, :self.ncells] @ vals[first]
                self.assertLess(np.max(np.abs(got - sign*tt)/tt), 1.e-9)

        with self.subTest(aniso='vti_psv'):
            vals = model('vti_psv', 0.2)
            tt, L = self._mesh('vti_psv', values=vals).raytrace(
                self.src, self.rcv, compute_L=True)
            L = L.toarray()
            got = (L[:, :self.ncells] @ vals['set_Vp0'] +
                   L[:, self.ncells:2*self.ncells] @ vals['set_Vs0'])
            self.assertLess(np.max(np.abs(got + tt)/tt), 1.e-9)

    def test_vertical_and_horizontal_are_the_axial_velocities(self):
        """Along the symmetry axis the group velocity is Vp0 exactly, and along
        the horizontal it is Vp0*sqrt(1+2*epsilon).  A straight path between two
        nodes of the mesh must reproduce those, up to what the shortest-path
        discretisation costs."""
        Vp0, eps = 3.094, 0.256
        g = self._mesh('vti_psv')
        src = np.array([[0.5, 0.5, 0.0]])
        rcv = np.array([[0.5, 0.5, 1.0],       # along the symmetry axis
                        [1.0, 0.5, 0.0]])      # perpendicular to it
        tt = g.raytrace(src, rcv)
        self.assertAlmostEqual(tt[0], 1.0/Vp0, delta=1.e-3*1.0/Vp0)
        v_h = Vp0*np.sqrt(1. + 2.*eps)
        self.assertAlmostEqual(tt[1], 0.5/v_h, delta=1.e-3*0.5/v_h)

    def test_phase_and_pickle(self):
        import pickle
        g = self._mesh('vti_psv')
        tt_qp = g.raytrace(self.src, self.rcv)
        g.set_phase('qSV')
        tt_qsv = g.raytrace(self.src, self.rcv)
        self.assertTrue(np.all(tt_qsv > tt_qp))

        g2 = pickle.loads(pickle.dumps(g))
        for setter, base, _, _ in self.MEDIA['vti_psv']:
            getattr(g2, setter)(np.full(self.ncells, base))
        np.testing.assert_allclose(g2.raytrace(self.src, self.rcv), tt_qsv,
                                   rtol=1.e-12)

    def test_refused_where_not_implemented(self):
        with self.assertRaises(ValueError):
            tm.Mesh3d(self.nodes, self.tet, method='FSM', aniso='vti_psv')
        with self.assertRaises(ValueError):
            tm.Mesh3d(self.nodes, self.tet, method='SPM', aniso='vti_psv',
                      cell_slowness=0)
        with self.assertRaises(ValueError):
            tm.Mesh3d(self.nodes, self.tet, method='SPM', aniso='tti_psv')


class TestRaypathOrder(unittest.TestCase):
    """Raypath coordinates run from source to receiver, whatever the solver.

    Grid3Dun already reordered explicitly; Grid3Duc returned the raypath as
    the steepest descent built it, from the receiver down.  Both are checked
    here, i.e. node and cell slowness.
    """

    def setUp(self):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName('./files/layers_medium.vtu')
        reader.Update()
        data = reader.GetOutput()
        self.nodes = np.array([data.GetPoint(n)
                               for n in range(data.GetNumberOfPoints())])
        self.tet = np.empty((data.GetNumberOfCells(), 4), dtype=int)
        ind = vtk.vtkIdList()
        for n in range(data.GetNumberOfCells()):
            data.GetCellPoints(n, ind)
            self.tet[n] = [ind.GetId(i) for i in range(4)]
        self.slow_c = vtk_to_numpy(data.GetCellData().GetArray('Slowness'))
        self.slow_n = 1.0 / (2.0 + 0.15 * self.nodes[:, 2])
        self.src = np.array([[1.0, 1.0, 1.0]])
        self.rcv = np.array([[8.0, 8.0, 8.0]])

    def test_source_first(self):
        for method, kwargs in (('FSM', {}),
                               ('SPM', dict(n_secondary=3)),
                               ('DSPM', dict(n_secondary=2, n_tertiary=2))):
            for cell_slowness in (0, 1):
                with self.subTest(method=method, cell_slowness=cell_slowness):
                    g = tm.Mesh3d(self.nodes, self.tet, n_threads=1,
                                  method=method, cell_slowness=cell_slowness,
                                  **kwargs)
                    g.set_slowness(self.slow_c if cell_slowness else self.slow_n)
                    r = np.asarray(g.raytrace(self.src, self.rcv,
                                              return_rays=True)[1][0])
                    self.assertGreater(r.shape[0], 2)
                    np.testing.assert_allclose(r[-1], self.rcv[0], atol=1e-6)
                    self.assertLess(np.linalg.norm(r[0] - self.src[0]),
                                    np.linalg.norm(r[0] - self.rcv[0]))

# -*- coding: utf-8 -*-
"""Tests for verifying python wrappers, module tmesh in 3D"""

import itertools
import unittest

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from scipy.io import mmread
from scipy.spatial import Delaunay

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


def tet_volumes(nodes, tet):
    """Volumes of tetrahedra."""
    edges = nodes[tet[:, 1:]] - nodes[tet[:, :1]]
    return np.abs(np.linalg.det(edges)) / 6.


def without_flat_cells(nodes, tet):
    """Drop the flat cells that Delaunay leaves on the boundary of a box."""
    return tet[tet_volumes(nodes, tet) > 1.e-9].astype(np.int64)


def cube_mesh(n):
    """Unit cubes split in 6 tetrahedra along their main diagonal, none flat."""
    x = np.arange(float(n))
    nodes = np.array(np.meshgrid(x, x, x, indexing='ij')).reshape(3, -1).T
    index = lambda c: c[0]*n*n + c[1]*n + c[2]
    tet = []
    for corner in itertools.product(range(n-1), repeat=3):
        for axes in itertools.permutations(range(3)):
            c = list(corner)
            cell = [index(c)]
            for ax in axes:
                c[ax] += 1
                cell.append(index(c))
            tet.append(cell)
    return nodes, np.array(tet, dtype=np.int64)


def jittered_mesh(n):
    """Irregular interior nodes of a box of unit cells, without flat cells."""
    x = np.arange(float(n))
    nodes = np.array(np.meshgrid(x, x, x, indexing='ij')).reshape(3, -1).T
    inner = np.all((nodes > 0.) & (nodes < n - 1.), axis=1)
    nodes[inner] += np.random.default_rng(1).uniform(-.2, .2, (inner.sum(), 3))
    return nodes, without_flat_cells(nodes, Delaunay(nodes).simplices)


METHODS = (('FSM', {}),
           ('SPM', {'n_secondary': 3}),
           ('DSPM', {'n_secondary': 2, 'n_tertiary': 2}))


class TestFlatTetrahedra(unittest.TestCase):
    """Solvers return nan, crash or never return on flat cells: refuse them."""

    def lattice(self):
        x = np.arange(4.)
        return np.array(np.meshgrid(x, x, x, indexing='ij')).reshape(3, -1).T

    def test_regular_lattice_refused(self):
        nodes = self.lattice()
        tet = Delaunay(nodes).simplices.astype(np.int64)
        self.assertLess(tet_volumes(nodes, tet).min(), 1.e-9)
        for cell_slowness in (True, False):
            with self.subTest(cell_slowness=cell_slowness):
                with self.assertRaisesRegex(ValueError, 'no volume'):
                    tm.Mesh3d(nodes, tet, cell_slowness=cell_slowness)

    def test_flat_boundary_cells_refused(self):
        # interior nodes are irregular, but the faces of the box are coplanar
        nodes = jittered_mesh(5)[0]
        with self.assertRaisesRegex(ValueError, 'no volume'):
            tm.Mesh3d(nodes, Delaunay(nodes).simplices.astype(np.int64))

    def test_scale_does_not_matter(self):
        nodes, tet = cube_mesh(3)
        for scale in (1.e-4, 1.e4):
            with self.subTest(scale=scale):
                tm.Mesh3d(nodes * scale, tet)
                flat = tet.copy()
                flat[0, 3] = flat[0, 2]
                with self.assertRaisesRegex(ValueError, 'no volume'):
                    tm.Mesh3d(nodes * scale, flat)


class TestPointsOnMeshFeatures(unittest.TestCase):
    """Raypaths starting or ending exactly on a node, a face or an edge."""

    def setUp(self):
        self.nodes, self.tet = cube_mesh(9)
        # for V = a + b z, the traveltime is known exactly
        self.a, self.b = 2., .5

    def exact(self, src, rcv):
        v1 = self.a + self.b * src[2]
        v2 = self.a + self.b * rcv[:, 2]
        d2 = np.sum((rcv - src)**2, axis=1)
        return np.arccosh(1. + self.b**2 * d2 / (2. * v1 * v2)) / self.b

    def mesh(self, nodes, tet, cell_slowness, method, kwargs):
        g = tm.Mesh3d(nodes, tet, cell_slowness=cell_slowness, method=method,
                      tt_from_rp=True, **kwargs)
        z = nodes[tet].mean(axis=1)[:, 2] if cell_slowness else nodes[:, 2]
        g.set_slowness(1. / (self.a + self.b * z))
        return g

    def test_source_and_receivers_on_nodes(self):
        # such a source had no cell, and a receiver on a node was taken to be
        # in it: the raypath stopped at once, with a wrong slowness.  The source
        # stays away from cell 0, where this could go unnoticed.
        src = np.array([4., 4., 4.])
        rcv = np.array([[7., 7., 7.], [7., 1., 7.], [1., 6., 1.], [4., 4., 7.]])
        exact = self.exact(src, rcv)
        for cell_slowness, (method, kwargs) in itertools.product((True, False),
                                                                 METHODS):
            with self.subTest(cell_slowness=cell_slowness, method=method):
                g = self.mesh(self.nodes, self.tet, cell_slowness, method, kwargs)
                tt = g.raytrace(src[np.newaxis], rcv)
                np.testing.assert_allclose(tt, exact, rtol=0.1)

    def test_source_on_node_in_a_uniform_model(self):
        src = np.array([[4., 4., 4.]])
        rcv = np.array([[7., 7., 7.], [4., 4., 7.], [8., 4., 4.]])
        for cell_slowness, (method, kwargs) in itertools.product((True, False),
                                                                 METHODS):
            with self.subTest(cell_slowness=cell_slowness, method=method):
                g = tm.Mesh3d(self.nodes, self.tet, cell_slowness=cell_slowness,
                              method=method, tt_from_rp=True, **kwargs)
                g.set_slowness(np.ones(g.nparams))
                tt, rays = g.raytrace(src, rcv, return_rays=True)
                length = [np.linalg.norm(np.diff(r, axis=0), axis=1).sum()
                          for r in rays]
                np.testing.assert_allclose(tt, length, rtol=1.e-6)
                np.testing.assert_allclose(tt, np.linalg.norm(rcv - src, axis=1),
                                           rtol=0.07)

    def test_raypath_along_an_inner_edge(self):
        # from the surface, the raypath runs exactly along edges inside the mesh:
        # looking for an external face to project on then crashed
        src = np.array([[4., 4., 4.]])
        rcv = np.array([[4., 4., 8.]])
        for cell_slowness, (method, kwargs) in itertools.product((True, False),
                                                                 METHODS):
            with self.subTest(cell_slowness=cell_slowness, method=method):
                g = tm.Mesh3d(self.nodes, self.tet, cell_slowness=cell_slowness,
                              method=method, tt_from_rp=True, **kwargs)
                g.set_slowness(np.ones(g.nparams))
                tt = g.raytrace(src, rcv)
                np.testing.assert_allclose(tt, 4., rtol=0.1)

    def test_receivers_on_the_surface(self):
        # with slowness in cells, raypath segments along the surface took no time
        nodes, tet = jittered_mesh(7)
        src = np.array([[3., 3., 3.]])
        rcv = nodes[np.any((nodes == 0.) | (nodes == 6.), axis=1)][::5]
        for cell_slowness, (method, kwargs) in itertools.product((True, False),
                                                                 METHODS):
            with self.subTest(cell_slowness=cell_slowness, method=method):
                g = tm.Mesh3d(nodes, tet, cell_slowness=cell_slowness,
                              method=method, tt_from_rp=True, **kwargs)
                g.set_slowness(np.ones(g.nparams))
                tt, rays = g.raytrace(src, rcv, return_rays=True)
                length = np.array([np.linalg.norm(np.diff(r, axis=0), axis=1).sum()
                                   for r in rays])
                # a raypath may give up, with a warning and no points
                traced = length > 0.
                self.assertGreaterEqual(traced.sum(), rcv.shape[0] - 1)
                np.testing.assert_allclose(tt[traced], length[traced], rtol=1.e-6)
                if cell_slowness and method != 'SPM':
                    # the same segments were missing from L
                    tt, L = g.raytrace(src, rcv, compute_L=True)
                    np.testing.assert_allclose((L @ np.ones(g.nparams))[traced],
                                               tt[traced], rtol=1.e-6)

    def test_receiver_on_a_face(self):
        # the cell first picked may be the wrong one of the two sharing the face:
        # with slowness at nodes, moving to the other one then crashed
        src = np.array([1., 1., 1.])
        rcv = self.nodes[self.tet[2000:2040, :3]].mean(axis=1)
        exact = self.exact(src, rcv)
        for cell_slowness in (True, False):
            with self.subTest(cell_slowness=cell_slowness):
                g = self.mesh(self.nodes, self.tet, cell_slowness, 'FSM', {})
                tt = g.raytrace(src[np.newaxis], rcv)
                np.testing.assert_allclose(tt, exact, rtol=0.06)

    def test_last_segment_on_an_interface(self):
        # the last segment runs along an edge shared by cells of slowness 10
        # and 1: it was timed in whichever of them came first, 10 or 1 times
        # the distance depending on the order of the cells
        src = np.array([[4., 4., 4.]])
        rcv = np.array([[4., 5., 5.]])
        orders = (np.arange(self.tet.shape[0]),
                  np.random.default_rng(0).permutation(self.tet.shape[0]))
        for order, (method, kwargs) in itertools.product(orders, METHODS[::2]):
            with self.subTest(method=method):
                tet = self.tet[order]
                slowness = np.where(self.nodes[tet].mean(axis=1)[:, 0] < 4., 10., 1.)
                g = tm.Mesh3d(self.nodes, tet, cell_slowness=True, method=method,
                              tt_from_rp=True, **kwargs)
                g.set_slowness(slowness)
                tt, L = g.raytrace(src, rcv, compute_L=True)
                # the walk along the interface still picks up some slowness 10
                self.assertLess(tt[0], 2. * np.sqrt(2.))
                np.testing.assert_allclose(L @ slowness, tt, rtol=1.e-6)

    def test_nodal_kernel_from_a_face(self):
        # starting on a face, the first segment was weighted on node 0, far away,
        # instead of the nodes of that face
        src = np.array([[6.2, 6.3, 6.1]])
        rcv = self.nodes[self.tet[3000, :3]].mean(axis=0)[np.newaxis]
        g = tm.Mesh3d(self.nodes, self.tet, cell_slowness=False, method='FSM',
                      tt_from_rp=True)
        g.set_slowness(np.ones(g.nparams))
        tt, L = g.raytrace(src, rcv, compute_L=True)
        L = L[0].toarray()
        np.testing.assert_allclose(L.sum(axis=1), tt, rtol=1.e-6)
        far = np.linalg.norm(self.nodes - rcv, axis=1) > 2.
        far &= np.linalg.norm(self.nodes - src, axis=1) > 2.
        self.assertEqual(np.abs(L[0, far]).max(), 0.)

    def test_receiver_on_a_secondary_node(self):
        # it was taken for a primary node, and part of the raypath was left out
        nodes, tet = jittered_mesh(9)
        src = nodes[tet[1000]].mean(axis=0)[np.newaxis]
        # with 2 secondary nodes per edge, one lies at the centroid of each face
        rcv = nodes[tet[2000, :3]].mean(axis=0)[np.newaxis]
        for cell_slowness in (True, False):
            with self.subTest(cell_slowness=cell_slowness):
                g = tm.Mesh3d(nodes, tet, cell_slowness=cell_slowness, method='DSPM',
                              tt_from_rp=True, n_secondary=2, n_tertiary=2)
                g.set_slowness(np.ones(g.nparams))
                tt, rays = g.raytrace(src, rcv, return_rays=True)
                length = np.linalg.norm(np.diff(rays[0], axis=0), axis=1).sum()
                self.assertAlmostEqual(tt[0], length, places=6)
                self.assertGreaterEqual(tt[0], np.linalg.norm(rcv - src) - 1.e-9)


class TestComputeKSquared(unittest.TestCase):
    """compute_K(order=2, squared=True) is the matrix square of order 1.

    The squaring used to read K * K.  That was a matrix product while ttcrpy
    returned sparse matrices; on the sparse arrays it returns now, * is
    elementwise, and the second-derivative operator came back as the
    elementwise square without any error.  Nothing exercised this path, so
    139 tests passed across the change.
    """

    def setUp(self):
        # a small mesh: compute_K takes minutes on layers_medium.vtu.  The
        # interior nodes are jittered, as compute_K rejects a regular lattice
        # as poorly conditioned.
        rng = np.random.default_rng(0)
        x = np.linspace(0., 10., 6)
        h = x[1] - x[0]
        nodes = np.array(np.meshgrid(x, x, x, indexing='ij')).reshape(3, -1).T
        inner = np.all((nodes > 0.) & (nodes < 10.), axis=1)
        nodes[inner] += rng.uniform(-0.3*h, 0.3*h, (inner.sum(), 3))
        self.nodes = nodes
        self.tet = without_flat_cells(nodes, Delaunay(nodes).simplices)

    def test_squared_is_the_matrix_square(self):
        # compute_K exists for slowness at nodes only
        g = tm.Mesh3d(self.nodes, self.tet, cell_slowness=False,
                      method='FSM', tt_from_rp=0)
        first = g.compute_K(order=1)
        squared = g.compute_K(order=2, squared=True)
        for name, a, b in zip('xyz', first, squared):
            with self.subTest(axis=name):
                scale = abs(b).max()
                self.assertLessEqual(abs(a @ a - b).max(), 1e-12 * scale)
                # and the check can tell the two apart
                self.assertGreater(abs(a.multiply(a) - b).max(), 1e-6 * scale)


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


class TestComputeH(unittest.TestCase):
    """The Jacobian returned by compute_H, on a tetrahedral mesh.

    Same identity as the rectilinear case: the spatial columns are
    -s(x_s) times the unit take-off direction.  Checked against a constant
    velocity gradient imposed on the mesh, for which the take-off is known in
    closed form.
    """

    V0 = 2.0
    A = 0.15

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
        self.slow = 1.0 / (self.V0 + self.A * self.nodes[:, 2])
        self.xs = np.array([1.0, 1.0, 1.0])
        # geometries chosen to avoid a pre-existing crash in the mesh solver
        # for certain source/receiver pairs, unrelated to compute_H
        self.rcv = np.array([[8.0, 8.0, 8.0], [6.0, 9.0, 7.0],
                             [7.0, 7.0, 11.0]])
        self.src = np.repeat(self.xs.reshape(1, 3), self.rcv.shape[0], axis=0)

    def _dT_analytic(self):
        d = 1.e-5
        out = np.zeros((self.rcv.shape[0], 3))

        def T(p, s):
            v1 = self.V0 + self.A * s[2]
            v2 = self.V0 + self.A * p[2]
            r2 = np.sum((p - s) ** 2)
            return np.arccosh(1.0 + self.A * self.A * r2 / (2 * v1 * v2)) / self.A

        for k in range(self.rcv.shape[0]):
            for j in range(3):
                a = self.xs.copy()
                b = self.xs.copy()
                a[j] += d
                b[j] -= d
                out[k, j] = (T(self.rcv[k], a) - T(self.rcv[k], b)) / (2 * d)
        return out

    def test_against_analytic(self):
        dT = self._dT_analytic()
        g = tm.Mesh3d(self.nodes, self.tet, n_threads=1, cell_slowness=0,
                      method='SPM', n_secondary=3)
        g.set_slowness(self.slow)
        tt, H = g.compute_H(self.src, self.rcv)

        self.assertEqual(H.shape, (self.rcv.shape[0], 4))
        np.testing.assert_allclose(H[:, 0], 1.0)
        for k in range(self.rcv.shape[0]):
            hk = H[k, 1:]
            ang = np.degrees(np.arccos(np.clip(
                hk @ dT[k] / (np.linalg.norm(hk) * np.linalg.norm(dT[k])),
                -1, 1)))
            self.assertLess(ang, 3.0,
                            'rcv %d: take-off off by %.2f deg' % (k, ang))
            # |grad T| = s is the eikonal equation itself
            rel = abs(np.linalg.norm(hk) - np.linalg.norm(dT[k])) / \
                np.linalg.norm(dT[k])
            self.assertLess(rel, 0.02, 'rcv %d: |H| off by %.1f%%' % (k, 100 * rel))

    def test_columns_and_traveltimes(self):
        g = tm.Mesh3d(self.nodes, self.tet, n_threads=1, cell_slowness=0,
                      method='SPM', n_secondary=3)
        g.set_slowness(self.slow)
        tt, H = g.compute_H(self.src, self.rcv)
        tt2, H2 = g.compute_H(self.src, self.rcv, full=False)
        self.assertEqual(H2.shape, (self.rcv.shape[0], 2))
        np.testing.assert_allclose(H2, H[:, 1:3])
        np.testing.assert_allclose(tt, g.raytrace(self.src, self.rcv))

    def test_thread_no_range(self):
        """A thread number the mesh does not have is refused, not used."""
        for n_threads in (1, 2):
            g = tm.Mesh3d(self.nodes, self.tet, n_threads=n_threads,
                          cell_slowness=0, method='SPM', n_secondary=3)
            g.set_slowness(self.slow)
            tt, H = g.compute_H(self.src, self.rcv)
            for thread_no in (-1, n_threads, 1000):
                with self.subTest(n_threads=n_threads, thread_no=thread_no):
                    with self.assertRaises(ValueError):
                        g.compute_H(self.src, self.rcv, slowness=2 * self.slow,
                                    thread_no=thread_no)
            # the last valid thread gives the same answer, so the calls
            # refused above did not change the slowness either
            tt2, H2 = g.compute_H(self.src, self.rcv, thread_no=n_threads - 1)
            np.testing.assert_allclose(tt2, tt)
            np.testing.assert_allclose(H2, H)

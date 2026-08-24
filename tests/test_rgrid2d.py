# -*- coding: utf-8 -*-
"""Tests for verifying python wrappers, module rgrid in 2D"""

import unittest

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

import ttcrpy.rgrid as rg


def get_tt(filename):
    reader = vtk.vtkXMLRectilinearGridReader()
    reader.SetFileName(filename)
    reader.Update()
    data = reader.GetOutput()
    x = vtk_to_numpy(data.GetXCoordinates())
    z = vtk_to_numpy(data.GetZCoordinates())
    names = ('Travel Time', 'Travel time', 'travel time')
    for name in names:
        if data.GetPointData().HasArray(name):
            break
    tt = vtk_to_numpy(data.GetPointData().GetArray(name))
    dim = (x.size, z.size)
    return tt.reshape(dim, order='F').flatten()


class TestGrid2dc(unittest.TestCase):

    def setUp(self):
        reader = vtk.vtkXMLRectilinearGridReader()
        reader.SetFileName('./files/layers_fine2d.vtr')
        reader.Update()

        data = reader.GetOutput()
        self.x = vtk_to_numpy(data.GetXCoordinates())
        self.z = vtk_to_numpy(data.GetZCoordinates())

        self.slowness = vtk_to_numpy(data.GetCellData().GetArray('Slowness'))
        dim = (self.x.size-1, self.z.size-1)
        self.slowness = self.slowness.reshape(dim, order='F').flatten()

        self.src = np.loadtxt('./files/src2d.dat', skiprows=1)
        # we roll because file has x z t0 and we want t0 x z
        self.src = np.roll(self.src, 1).reshape((1, 3))
        self.rcv = np.loadtxt('./files/rcv2d.dat', skiprows=1)
        self.src_in = np.loadtxt('./files/src2d_in.dat',skiprows=1)
        self.src_in = np.roll(self.src_in, 1).reshape((1, 3))
        self.rcv_in = np.loadtxt('./files/rcv2d_in.dat',skiprows=1)

    def test_Grid2Dfs(self):
        g = rg.Grid2d(self.x, self.z, method='FSM')
        tt = g.raytrace(self.src, self.rcv, slowness=self.slowness)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid2Drcfs_tt_grid.vtr')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'FSM accuracy failed (slowness in cells)')

    def test_Grid2Dsp(self):
        g = rg.Grid2d(self.x, self.z, method='SPM', nsnx=10, nsnz=10)
        tt = g.raytrace(self.src, self.rcv, slowness=self.slowness)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid2Drcsp_tt_grid.vtr')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'SPM accuracy failed (slowness in cells)')

    def test_Grid2Ddsp(self):
        g = rg.Grid2d(self.x, self.z, method='DSPM', n_secondary=3,
                      n_tertiary=3, radius_factor_tertiary=3.0)
        tt = g.raytrace(self.src, self.rcv, slowness=self.slowness)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid2Drcdsp_tt_grid.vtr')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'DSPM accuracy failed (slowness in cells)')

    def test_Grid2Dfs_L(self):
        g = rg.Grid2d(self.x, self.z, method='FSM')
        tt, L = g.raytrace(self.src_in, self.rcv_in, slowness=self.slowness, compute_L=True)
        tt2 = L @ self.slowness
        self.assertLess(np.sum(np.abs(tt-tt2))/tt.size, 0.01,
                        'FSM_L accuracy failed (slowness in cells)')

    def test_Grid2Dsp_L(self):
        g = rg.Grid2d(self.x, self.z, method='SPM', nsnx=10, nsnz=10)
        tt, L = g.raytrace(self.src_in, self.rcv_in, slowness=self.slowness, compute_L=True)
        tt2 = L @ self.slowness
        self.assertLess(np.sum(np.abs(tt-tt2))/tt.size, 0.01,
                        'SPM_L accuracy failed (slowness in cells)')

    def test_Grid2Ddsp_L(self):
        g = rg.Grid2d(self.x, self.z, method='DSPM', n_secondary=3,
                      n_tertiary=3, radius_factor_tertiary=3.0)
        tt, L = g.raytrace(self.src_in, self.rcv_in, slowness=self.slowness, compute_L=True)
        tt2 = L @ self.slowness
        self.assertLess(np.sum(np.abs(tt-tt2))/tt.size, 0.01,
                        'DSPM_L accuracy failed (slowness in cells)')


class TestGrid2dn(unittest.TestCase):

    def setUp(self):
        reader = vtk.vtkXMLRectilinearGridReader()
        reader.SetFileName('./files/gradient_fine2d.vtr')
        reader.Update()

        data = reader.GetOutput()
        self.x = vtk_to_numpy(data.GetXCoordinates())
        self.z = vtk_to_numpy(data.GetZCoordinates())

        self.slowness = vtk_to_numpy(data.GetPointData().GetArray('Slowness'))
        dim = (self.x.size, self.z.size)
        self.slowness = self.slowness.reshape(dim, order='F').flatten()

        self.src = np.loadtxt('./files/src2d.dat', skiprows=1)
        self.src = np.roll(self.src, 1).reshape((1, 3))
        self.rcv = np.loadtxt('./files/rcv2d.dat', skiprows=1)

    def test_Grid2Dfs(self):
        g = rg.Grid2d(self.x, self.z, method='FSM', cell_slowness=0)
        tt = g.raytrace(self.src, self.rcv, slowness=self.slowness)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid2Drnfs_tt_grid.vtr')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'FSM accuracy failed (slowness at nodes)')

    def test_Grid2Dsp(self):
        g = rg.Grid2d(self.x, self.z, method='SPM', nsnx=10, nsnz=10,
                      cell_slowness=0)
        tt = g.raytrace(self.src, self.rcv, slowness=self.slowness)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid2Drnsp_tt_grid.vtr')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'SPM accuracy failed (slowness at nodes)')

    def test_Grid2Ddsp(self):
        g = rg.Grid2d(self.x, self.z, method='DSPM', n_secondary=3,
                      n_tertiary=3, radius_factor_tertiary=3.0, cell_slowness=0)
        tt = g.raytrace(self.src, self.rcv, slowness=self.slowness)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid2Drndsp_tt_grid.vtr')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'DSPM accuracy failed (slowness at nodes)')

class TestAniso(unittest.TestCase):

    def setUp(self):
        reader = vtk.vtkXMLRectilinearGridReader()
        reader.SetFileName('./files/elliptical_fine2d.vtr')
        reader.Update()

        data = reader.GetOutput()
        self.x = vtk_to_numpy(data.GetXCoordinates())
        self.z = vtk_to_numpy(data.GetZCoordinates())

        self.slowness = vtk_to_numpy(data.GetCellData().GetArray('Slowness'))
        dim = (self.x.size-1, self.z.size-1)
        self.slowness = self.slowness.reshape(dim, order='F').flatten()

        self.xi = vtk_to_numpy(data.GetCellData().GetArray('xi'))
        self.xi = self.xi.reshape(dim, order='F').flatten()

        self.src = np.loadtxt('./files/src2d.dat', skiprows=1)
        self.src = np.roll(self.src, 1).reshape((1, 3))
        self.rcv = np.loadtxt('./files/rcv2daniso.dat', skiprows=1)

    def test_Grid2Dsp(self):
        g = rg.Grid2d(self.x, self.z, method='SPM', nsnx=10, nsnz=10, aniso='elliptical')
        g.set_slowness(self.slowness)
        g.set_xi(self.xi)
        tt = g.raytrace(self.src, self.rcv)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid2Drcsp_tt_grid_elliptical.vtr')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'SPM accuracy failed (elliptical anisotropy)')

class TestWeakly(unittest.TestCase):

    def setUp(self):
        reader = vtk.vtkXMLRectilinearGridReader()
        reader.SetFileName('./files/weakly_an_fine2d.vtr')
        reader.Update()

        data = reader.GetOutput()
        self.x = vtk_to_numpy(data.GetXCoordinates())
        self.z = vtk_to_numpy(data.GetZCoordinates())

        v = vtk_to_numpy(data.GetCellData().GetArray('Velocity'))
        dim = (self.x.size-1, self.z.size-1)
        self.slowness = 1/v.reshape(dim, order='F').flatten()

        self.s2 = vtk_to_numpy(data.GetCellData().GetArray('s2'))
        self.s2 = self.s2.reshape(dim, order='F').flatten()
        self.s4 = vtk_to_numpy(data.GetCellData().GetArray('s4'))
        self.s4 = self.s4.reshape(dim, order='F').flatten()

        self.src = np.loadtxt('./files/src2d.dat', skiprows=1)
        self.src = np.roll(self.src, 1).reshape((1, 3))
        self.rcv = np.loadtxt('./files/rcv2daniso.dat', skiprows=1)

    def test_Grid2Dsp(self):
        g = rg.Grid2d(self.x, self.z, method='SPM', nsnx=10, nsnz=10, aniso='weakly_anelliptical')
        g.set_slowness(self.slowness)
        g.set_s2(self.s2)
        g.set_s4(self.s4)
        tt = g.raytrace(self.src, self.rcv)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten()
        tt_ref = get_tt('./files/Grid2Drcsp_tt_grid_weakly.vtr')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'SPM accuracy failed (weakly anelliptical)')
#

class TestComputeD(unittest.TestCase):

    def test_cell_slowness(self):
        reader = vtk.vtkXMLRectilinearGridReader()

        reader.SetFileName('./files/layers_fine2d.vtr')
        reader.Update()

        data = reader.GetOutput()
        x = vtk_to_numpy(data.GetXCoordinates())
        z = vtk_to_numpy(data.GetZCoordinates())

        slowness = vtk_to_numpy(data.GetCellData().GetArray('Slowness'))
        dim = (x.size - 1, z.size - 1)
        slowness = slowness.reshape(dim, order='F').flatten()

        g = rg.Grid2d(x, z, method='FSM')

        dx = x[1] - x[0]
        xc = x[:-1] + 0.5*(dx)
        zc = z[:-1] + 0.5*(z[1] - z[0])

        rng = np.random.default_rng()
        tmp = np.meshgrid(xc, zc, indexing='ij')
        coord = np.c_[tmp[0].ravel(), tmp[1].ravel()]
        coord += rng.uniform(low=-dx/3, high=dx/3, size=coord.shape)

        pts = vtk.vtkPoints()
        for n in range(coord.shape[0]):
            pts.InsertNextPoint(coord[n, 0], coord[n, 1], 0.0)

        ppts = vtk.vtkPolyData()
        ppts.SetPoints(pts)

        delaunay = vtk.vtkDelaunay2D()
        delaunay.SetInputData(ppts)
        delaunay.Update()

        rotation = vtk.vtkTransform()
        rotation.RotateX(90)
        rotation.Scale(1.0, 0.0, 1.0)

        transformFilter = vtk.vtkTransformPolyDataFilter()
        transformFilter.SetInputConnection(delaunay.GetOutputPort())
        transformFilter.SetTransform(rotation)
        transformFilter.Update()

        ppts = transformFilter.GetOutput()

        pts = vtk_to_numpy(ppts.GetPoints().GetData())
        # we must recreate coord, order of pts was changed by triangulation
        coord = np.array(np.c_[pts[:, 0], pts[:, 2]], dtype=np.float64)

        D = g.compute_D(coord)
        s1 = D @ slowness

        pf = vtk.vtkProbeFilter()
        pf.SetInputData(ppts)
        pf.SetSourceData(data)
        pf.Update()

        temp = pf.GetOutput().GetPointData().GetArray('Slowness')

        s2 = vtk_to_numpy(temp)

        self.assertAlmostEqual(np.sum(np.abs(s1 - s2)), 0.0)

    def test_node_slowness(self):
        reader = vtk.vtkXMLRectilinearGridReader()

        reader.SetFileName('./files/gradient_fine2d.vtr')
        reader.Update()

        data = reader.GetOutput()
        x = vtk_to_numpy(data.GetXCoordinates())
        z = vtk_to_numpy(data.GetZCoordinates())

        slowness = vtk_to_numpy(data.GetPointData().GetArray('Slowness'))
        dim = (x.size, z.size)
        slowness = slowness.reshape(dim, order='F').flatten()

        g = rg.Grid2d(x, z, method='FSM', cell_slowness=0)

        rng = np.random.default_rng()
        coord = np.c_[rng.uniform(0.001, 20.0, 100), rng.uniform(0.001, 20.0, 100)]
        D = g.compute_D(coord)
        s1 = D @ slowness

        pts = vtk.vtkPoints()
        for n in range(coord.shape[0]):
            pts.InsertNextPoint(coord[n, 0], 0.0, coord[n, 1])

        ppts = vtk.vtkPolyData()
        ppts.SetPoints(pts)

        interpolator = vtk.vtkPointInterpolator()
        interpolator.SetInputData(ppts)
        interpolator.SetSourceData(data)
        interpolator.Update()

        s2 = vtk_to_numpy(interpolator.GetOutput().GetPointData().GetArray('Slowness'))

        self.assertLess(np.sum(np.abs(s1-s2))/coord.shape[0], 0.01, "compute_D accuracy failed slowness at nodes")

class Data_kernel(unittest.TestCase):

    def test_2d(self):

        V = np.ones((11, 13))
        V[:, 7:] = 2
        slowness = 1. / V.flatten()

        grx = np.arange(12.)
        grz = np.arange(14.)

        z = 0.5 + np.arange(13.)
        Tx = np.vstack((0.5+np.zeros((13,)), z)).T
        Rx = np.vstack((10.5+np.zeros((13,)), z)).T
        nTx = Tx.shape[0]
        nRx = Rx.shape[0]
        Tx = np.kron(Tx, np.ones((nRx, 1)))
        Rx = np.kron(np.ones((nTx, 1)), Rx)

        L = rg.Grid2d.data_kernel_straight_rays(Tx, Rx, grx, grz)
        tt = L.dot(slowness)

        tt2 = np.zeros(tt.shape)
        d = np.sqrt(np.sum((Tx-Rx)**2, axis=1))

        ind = np.logical_and(Tx[:, 1] > 7, Rx[:, 1] > 7)
        tt2[ind] = d[ind]/2

        ind2 = np.logical_and(Tx[:, 1] < 7, Rx[:, 1] < 7)
        tt2[ind2] = d[ind2]

        ind3 = np.logical_and(np.logical_not(ind), np.logical_not(ind2))

        f = (7-Tx[ind3, 1]) / (Rx[ind3, 1]-Tx[ind3, 1])
        ind = (Rx[ind3, 1]-Tx[ind3, 1]) < 0
        f[ind] = 1-f[ind]
        tt2[ind3] = d[ind3]*f + d[ind3]*(1-f)/2

        self.assertAlmostEqual(np.sum(np.abs(tt-tt2)), 0.0)


class TestVTI(unittest.TestCase):
    """Traveltimes in homogeneous VTI media, qP/qSV and SH phases.

    The traveltime along a ray segment is its length divided by the *group*
    (energy) velocity taken in the direction of the segment.  Dividing by the
    phase velocity instead underestimates it by several percent away from the
    symmetry axes, phase and group velocities coinciding only along those axes.
    """

    vp, vs, eps, dlt, gam = 3.094, 1.51, 0.256, -0.0505, 0.15

    @classmethod
    def _phase(cls, theta, sign):
        """Phase velocity of the qP (sign=+1) or qSV (sign=-1) wave, and its
        derivative with respect to the phase angle."""
        f = 1. - cls.vs**2 / cls.vp**2
        s = np.sin(theta)**2
        a = 1. + 2.*cls.eps*s/f
        r = a*a - 8.*(cls.eps-cls.dlt)*s*(1.-s)/f
        g = 1. + cls.eps*s - f/2. + sign*(f/2.)*np.sqrt(np.maximum(r, 0.))
        drds = 2.*a*(2.*cls.eps/f) - 8.*(cls.eps-cls.dlt)*(1.-2.*s)/f
        dgds = cls.eps + sign*(f/4.)*drds/np.sqrt(np.maximum(r, 1.e-300))
        v = cls.vp*np.sqrt(g)
        return v, cls.vp*dgds*np.sin(2.*theta)/(2.*np.sqrt(g))

    @classmethod
    def _group_vel(cls, psi, sign):
        """Group velocity of the first arrival, for a ray angle psi."""
        theta = np.linspace(-np.pi, np.pi, 400001)
        v, dv = cls._phase(theta, sign)
        vg = np.sqrt(v*v + dv*dv)
        ray = theta + np.arctan2(dv, v)
        d = ((ray - psi + np.pi) % (2.*np.pi)) - np.pi
        return max(vg[k] for k in np.where(np.diff(np.sign(d)) != 0)[0])

    def setUp(self):
        n, h = 101, 0.02
        self.x = np.arange(n)*h
        self.z = np.arange(n)*h
        self.ncells = (n-1)*(n-1)
        self.radius = 0.7
        self.angles = np.radians(np.arange(5., 90., 5.))
        self.src = np.array([[1., 1.]])
        self.rcv = np.column_stack([1. + self.radius*np.sin(self.angles),
                                    1. + self.radius*np.cos(self.angles)])

    def test_vti_psv_homogeneous(self):
        # the python wrapper does not call setPhase(), so the qP phase is used
        g = rg.Grid2d(self.x, self.z, method='SPM', nsnx=15, nsnz=15,
                      aniso='vti_psv')
        g.set_Vp0(np.full(self.ncells, self.vp))
        g.set_Vs0(np.full(self.ncells, self.vs))
        g.set_epsilon(np.full(self.ncells, self.eps))
        g.set_delta(np.full(self.ncells, self.dlt))
        tt = g.raytrace(self.src, self.rcv)
        ref = np.array([self.radius/self._group_vel(a, 1.) for a in self.angles])
        self.assertLess(np.max(np.abs(tt-ref)/ref), 0.005)

    def test_vti_sh_homogeneous(self):
        vs0 = 1.8
        g = rg.Grid2d(self.x, self.z, method='SPM', nsnx=15, nsnz=15,
                      aniso='vti_sh')
        g.set_Vs0(np.full(self.ncells, vs0))
        g.set_gamma(np.full(self.ncells, self.gam))
        tt = g.raytrace(self.src, self.rcv)
        # the SH phase velocity is elliptical, so the traveltime is closed form
        lx = self.radius*np.sin(self.angles)
        lz = self.radius*np.cos(self.angles)
        ref = np.sqrt(lx*lx/(1. + 2.*self.gam) + lz*lz)/vs0
        self.assertLess(np.max(np.abs(tt-ref)/ref), 0.005)

    def test_vti_sh_matches_elliptical(self):
        # the same medium described through two different cell classes
        vs0 = 1.8
        xi = np.sqrt(1. + 2.*self.gam)     # set_xi squares its argument
        gh = rg.Grid2d(self.x, self.z, method='SPM', nsnx=10, nsnz=10,
                       aniso='vti_sh')
        gh.set_Vs0(np.full(self.ncells, vs0))
        gh.set_gamma(np.full(self.ncells, self.gam))
        ge = rg.Grid2d(self.x, self.z, method='SPM', nsnx=10, nsnz=10,
                       aniso='elliptical')
        ge.set_slowness(np.full(self.ncells, 1./(vs0*xi)))
        ge.set_xi(np.full(self.ncells, xi))
        self.assertLess(np.max(np.abs(gh.raytrace(self.src, self.rcv) -
                                      ge.raytrace(self.src, self.rcv))), 1.e-12)


class TestTTI(unittest.TestCase):
    """Traveltimes in homogeneous TTI media.

    Tilting a transversely isotropic medium rotates its group-velocity surface
    rigidly, so a TTI cell with a zero tilt angle must reproduce its VTI
    counterpart, and a tilted one must reproduce it with the ray angle shifted
    by the tilt.  The tilt convention is that of CellTiltedElliptical: the
    symmetry axis lies at -theta from the vertical.
    """

    vp, vs, eps, dlt, gam = 3.094, 1.51, 0.256, -0.0505, 0.15
    vs0 = 1.8

    def setUp(self):
        n, h = 101, 0.02
        self.x = np.arange(n)*h
        self.z = np.arange(n)*h
        self.ncells = (n-1)*(n-1)
        self.radius = 0.7
        self.angles = np.radians(np.arange(-85., 90., 5.))
        self.src = np.array([[1., 1.]])
        self.rcv = np.column_stack([1. + self.radius*np.sin(self.angles),
                                    1. + self.radius*np.cos(self.angles)])

    def _psv(self, aniso, tilt=None):
        g = rg.Grid2d(self.x, self.z, method='SPM', nsnx=15, nsnz=15,
                      aniso=aniso)
        g.set_Vp0(np.full(self.ncells, self.vp))
        g.set_Vs0(np.full(self.ncells, self.vs))
        g.set_epsilon(np.full(self.ncells, self.eps))
        g.set_delta(np.full(self.ncells, self.dlt))
        if tilt is not None:
            g.set_tilt_angle(np.full(self.ncells, tilt))
        return g

    def _sh(self, aniso, tilt=None):
        g = rg.Grid2d(self.x, self.z, method='SPM', nsnx=15, nsnz=15,
                      aniso=aniso)
        g.set_Vs0(np.full(self.ncells, self.vs0))
        g.set_gamma(np.full(self.ncells, self.gam))
        if tilt is not None:
            g.set_tilt_angle(np.full(self.ncells, tilt))
        return g

    def test_tti_psv_zero_tilt(self):
        tt_tti = self._psv('tti_psv', 0.0).raytrace(self.src, self.rcv)
        tt_vti = self._psv('vti_psv').raytrace(self.src, self.rcv)
        self.assertLess(np.max(np.abs(tt_tti-tt_vti)), 1.e-12)

    def test_tti_sh_zero_tilt(self):
        tt_tti = self._sh('tti_sh', 0.0).raytrace(self.src, self.rcv)
        tt_vti = self._sh('vti_sh').raytrace(self.src, self.rcv)
        self.assertLess(np.max(np.abs(tt_tti-tt_vti)), 1.e-12)

    def test_tti_sh_matches_tilted_elliptical(self):
        # the SH phase velocity is elliptical, so a tilted SH medium is a
        # tilted elliptical one described through a different cell class
        tilt = np.radians(30.)
        xi = np.sqrt(1. + 2.*self.gam)     # set_xi squares its argument
        gh = self._sh('tti_sh', tilt)
        ge = rg.Grid2d(self.x, self.z, method='SPM', nsnx=15, nsnz=15,
                       aniso='tilted_elliptical')
        ge.set_slowness(np.full(self.ncells, 1./(self.vs0*xi)))
        ge.set_xi(np.full(self.ncells, xi))
        ge.set_tilt_angle(np.full(self.ncells, tilt))
        self.assertLess(np.max(np.abs(gh.raytrace(self.src, self.rcv) -
                                      ge.raytrace(self.src, self.rcv))), 1.e-12)

    def test_tti_psv_rotation(self):
        # a tilted medium probed at psi must match the untilted one at psi+tilt
        tilt = np.radians(30.)
        g_tti = self._psv('tti_psv', tilt)
        g_vti = self._psv('vti_psv')
        rcv_rot = np.column_stack([
            1. + self.radius*np.sin(self.angles+tilt),
            1. + self.radius*np.cos(self.angles+tilt)])
        tt_tti = g_tti.raytrace(self.src, self.rcv)
        tt_vti = g_vti.raytrace(self.src, rcv_rot)
        self.assertLess(np.max(np.abs(tt_tti-tt_vti)/tt_vti), 0.002)

    def test_tti_psv_symmetry_axis(self):
        # qP is slowest along the symmetry axis, which lies at -tilt
        tilt = np.radians(35.)
        g = self._psv('tti_psv', tilt)
        rcv = np.array([[1.-self.radius*np.sin(tilt), 1.+self.radius*np.cos(tilt)],
                        [1.+self.radius*np.cos(tilt), 1.+self.radius*np.sin(tilt)]])
        tt = g.raytrace(self.src, rcv)
        v_axis, v_normal = self.radius/tt[0], self.radius/tt[1]
        self.assertLess(abs(v_axis/self.vp - 1.), 0.005)
        self.assertLess(abs(v_normal/(self.vp*np.sqrt(1.+2.*self.eps)) - 1.), 0.005)

    def test_tti_pickle(self):
        import pickle
        g = self._psv('tti_psv', np.radians(20.))
        self.assertEqual(pickle.loads(pickle.dumps(g)).__class__, g.__class__)


if __name__ == '__main__':

    unittest.main()

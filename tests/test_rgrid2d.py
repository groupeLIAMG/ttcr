# -*- coding: utf-8 -*-
"""Tests for verifying python wrappers, module rgrid in 2D"""

import contextlib
import os
import tempfile
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


class TestComputeL(unittest.TestCase):
    """Sensitivity matrix returned by compute_L.

    Two invariants are checked.  First, compute_L must not change the
    traveltimes: every raytrace() mode has to return the same values.  Second,
    for a homogeneous medium the rays are straight, so the row sum of the
    matrix is the length of the ray.

    What the anisotropic media report, and how wide a matrix they give, is
    checked per medium by TestSensitivity.
    """

    def setUp(self):
        n, h = 41, 0.05
        self.x = np.arange(n)*h
        self.z = np.arange(n)*h
        self.ncells = (n-1)*(n-1)
        self.src = np.array([[0.2, 0.2]])
        self.rcv = np.array([[1.7, 1.6], [1.8, 0.5], [0.5, 1.8]])
        self.dx = np.abs(self.rcv[:, 0] - self.src[0, 0])
        self.dz = np.abs(self.rcv[:, 1] - self.src[0, 1])

    def _grid(self, aniso):
        g = rg.Grid2d(self.x, self.z, method='SPM', nsnx=10, nsnz=10,
                      aniso=aniso)
        if aniso == 'iso':
            g.set_slowness(np.full(self.ncells, 0.5))
        elif aniso in ('elliptical', 'tilted_elliptical'):
            g.set_slowness(np.full(self.ncells, 0.5))
            g.set_xi(np.full(self.ncells, 1.1))
            if aniso == 'tilted_elliptical':
                g.set_tilt_angle(np.full(self.ncells, np.radians(20.)))
        elif aniso == 'weakly_anelliptical':
            g.set_slowness(np.full(self.ncells, 0.5))
            g.set_s2(np.full(self.ncells, 0.05))
            g.set_s4(np.full(self.ncells, 0.01))
        return g

    def test_compute_L_does_not_change_traveltimes(self):
        # regression: without its own l_data override, Grid2Drcsp fell back to
        # the gradient-following raypath of Grid2D and returned traveltimes
        # several percent higher than the other raytrace() modes
        for aniso in ('iso', 'elliptical', 'tilted_elliptical',
                      'weakly_anelliptical'):
            with self.subTest(aniso=aniso):
                tt_plain = self._grid(aniso).raytrace(self.src, self.rcv)
                tt_L, _ = self._grid(aniso).raytrace(self.src, self.rcv,
                                                     compute_L=True)
                out = self._grid(aniso).raytrace(self.src, self.rcv,
                                                 compute_L=True,
                                                 return_rays=True)
                self.assertLess(np.max(np.abs(tt_L-tt_plain)), 1.e-9)
                self.assertLess(np.max(np.abs(out[0]-tt_plain)), 1.e-9)

    def test_compute_L_iso(self):
        g = self._grid('iso')
        tt, L = g.raytrace(self.src, self.rcv, compute_L=True)
        self.assertEqual(L.shape, (self.rcv.shape[0], self.ncells))
        # straight rays in a homogeneous medium: row sum is the ray length
        length = np.sqrt(self.dx**2 + self.dz**2)
        self.assertLess(np.max(np.abs(np.asarray(L.sum(axis=1)).ravel()
                                      - length)/length), 0.01)
        # and L must reproduce the traveltimes
        self.assertLess(np.max(np.abs(L @ np.full(self.ncells, 0.5) - tt)),
                        1.e-9)

    def test_compute_L_receiver_on_source(self):
        # a receiver coinciding with the source has no path to walk; without a
        # guard the parent walk dereferenced nodes[T2::max()] and crashed
        rcv = np.vstack([self.src, self.rcv])
        for aniso in ('iso', 'elliptical', 'tilted_elliptical',
                      'weakly_anelliptical'):
            with self.subTest(aniso=aniso):
                tt, L = self._grid(aniso).raytrace(self.src, rcv,
                                                   compute_L=True)
                self.assertEqual(tt[0], 0.0)
                self.assertEqual(L[0].nnz, 0)      # empty row, zero path
                self.assertTrue(np.all(tt[1:] > 0.))

    def test_compute_L_elliptical(self):
        # the elliptical cells report derivatives of the traveltime with
        # respect to the slowness and the anisotropy ratio
        sx, xi = 0.5, 1.1
        g = self._grid('elliptical')
        tt, L = g.raytrace(self.src, self.rcv, compute_L=True)
        self.assertEqual(L.shape, (self.rcv.shape[0], 2*self.ncells))
        L = L.toarray()
        # dt is homogeneous of degree one in the slowness, so Euler's identity
        # makes the first block reproduce the traveltimes exactly
        self.assertLess(np.max(np.abs(L[:, :self.ncells] @
                                      np.full(self.ncells, sx) - tt)), 1.e-9)
        # and both blocks must match a finite difference of the traveltimes
        for blk, val, name in ((0, sx, 'slowness'), (1, xi, 'xi')):
            h = 1.e-4
            up, dn = [self._grid('elliptical') for _ in range(2)]
            if blk == 0:
                up.set_slowness(np.full(self.ncells, sx+h))
                dn.set_slowness(np.full(self.ncells, sx-h))
            else:
                up.set_xi(np.full(self.ncells, xi+h))
                dn.set_xi(np.full(self.ncells, xi-h))
            fd = (up.raytrace(self.src, self.rcv) -
                  dn.raytrace(self.src, self.rcv))/(2*h)
            ana = L[:, blk*self.ncells:(blk+1)*self.ncells].sum(axis=1)
            self.assertLess(np.max(np.abs(ana-fd)/np.abs(fd)), 5.e-3,
                            'd(tt)/d(%s) disagrees with a finite difference'
                            % name)



class TestSensitivity(unittest.TestCase):
    """The matrix returned by compute_L, for every anisotropy model.

    L holds one block of ncells columns per medium parameter.  Two things are
    asked of it: that each block matches a finite difference of the traveltimes
    with respect to that parameter, and that it satisfies the homogeneity
    identities, which hold exactly and need no finite difference.
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

    def setUp(self):
        n, h = 41, 0.05
        self.x = np.arange(n)*h
        self.z = np.arange(n)*h
        self.ncells = (n-1)*(n-1)
        self.src = np.array([[0.3, 0.3]])
        self.rcv = np.array([[1.7, 1.6], [1.8, 0.6], [0.6, 1.8]])

    def _grid(self, aniso, bump=None):
        """bump = (setter, value) replacing the base value of one parameter"""
        g = rg.Grid2d(self.x, self.z, method='SPM', nsnx=10, nsnz=10,
                      aniso=aniso)
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
        # cell scales with both of its vertical velocities together
        for aniso, blk, value, sign in (
                ('iso', 0, 0.5, 1.), ('elliptical', 0, 0.5, 1.),
                ('tilted_elliptical', 0, 0.5, 1.),
                ('weakly_anelliptical', 0, 0.5, 1.),
                ('vti_sh', 0, 1.8, -1.), ('tti_sh', 0, 1.8, -1.)):
            with self.subTest(aniso=aniso):
                tt, L = self._grid(aniso).raytrace(self.src, self.rcv,
                                                   compute_L=True)
                L = L.toarray()
                got = L[:, blk*self.ncells:(blk+1)*self.ncells] @ \
                    np.full(self.ncells, value)
                self.assertLess(np.max(np.abs(got - sign*tt)/tt), 1.e-9)

        for aniso in ('vti_psv', 'tti_psv'):
            with self.subTest(aniso=aniso):
                tt, L = self._grid(aniso).raytrace(self.src, self.rcv,
                                                   compute_L=True)
                L = L.toarray()
                got = (L[:, :self.ncells] @ np.full(self.ncells, 3.094) +
                       L[:, self.ncells:2*self.ncells] @
                       np.full(self.ncells, 1.51))
                self.assertLess(np.max(np.abs(got + tt)/tt), 1.e-9)


class TestPhase(unittest.TestCase):
    """Choosing the wave modelled in a transversely isotropic medium.

    The coupled cells describe both a quasi-compressional and a quasi-shear
    wave; the qP one is modelled until set_phase says otherwise.
    """

    def setUp(self):
        n, h = 41, 0.05
        self.x = np.arange(n)*h
        self.z = np.arange(n)*h
        self.ncells = (n-1)*(n-1)
        self.src = np.array([[0.3, 0.3]])
        self.rcv = np.array([[1.7, 1.6], [1.8, 0.6]])

    def _grid(self, aniso='vti_psv', phase=None):
        g = rg.Grid2d(self.x, self.z, method='SPM', nsnx=10, nsnz=10,
                      aniso=aniso)
        g.set_Vp0(np.full(self.ncells, 3.094))
        g.set_Vs0(np.full(self.ncells, 1.51))
        g.set_epsilon(np.full(self.ncells, 0.256))
        g.set_delta(np.full(self.ncells, -0.0505))
        if aniso == 'tti_psv':
            g.set_tilt_angle(np.full(self.ncells, 0.3))
        if phase is not None:
            g.set_phase(phase)
        return g

    def test_qSV_is_slower_than_qP(self):
        for aniso in ('vti_psv', 'tti_psv'):
            with self.subTest(aniso=aniso):
                tp = self._grid(aniso, 'qP').raytrace(self.src, self.rcv)
                ts = self._grid(aniso, 'qSV').raytrace(self.src, self.rcv)
                self.assertTrue(np.all(ts > tp))

    def test_qP_is_the_default(self):
        for aniso in ('vti_psv', 'tti_psv'):
            with self.subTest(aniso=aniso):
                self.assertTrue(np.allclose(
                    self._grid(aniso).raytrace(self.src, self.rcv),
                    self._grid(aniso, 'qP').raytrace(self.src, self.rcv)))

    def test_integers_of_the_cpp_setter(self):
        # setPhase takes 1 for qP and anything else for qSV
        self.assertTrue(np.allclose(
            self._grid('vti_psv', 1).raytrace(self.src, self.rcv),
            self._grid('vti_psv', 'qP').raytrace(self.src, self.rcv)))
        self.assertTrue(np.allclose(
            self._grid('vti_psv', 2).raytrace(self.src, self.rcv),
            self._grid('vti_psv', 'qSV').raytrace(self.src, self.rcv)))

    def test_sensitivity_follows_the_phase(self):
        # the identity must hold whichever wave is modelled
        for phase in ('qP', 'qSV'):
            with self.subTest(phase=phase):
                tt, L = self._grid('vti_psv', phase).raytrace(
                    self.src, self.rcv, compute_L=True)
                L = L.toarray()
                got = (L[:, :self.ncells] @ np.full(self.ncells, 3.094) +
                       L[:, self.ncells:2*self.ncells] @
                       np.full(self.ncells, 1.51))
                self.assertLess(np.max(np.abs(got + tt)/tt), 1.e-9)

    def test_a_medium_without_a_phase_refuses(self):
        g = rg.Grid2d(self.x, self.z, method='SPM', nsnx=10, nsnz=10,
                      aniso='elliptical')
        g.set_slowness(np.full(self.ncells, 0.5))
        g.set_xi(np.full(self.ncells, 1.1))
        self.assertRaises(RuntimeError, g.set_phase, 'qSV')
        self.assertRaises(ValueError, self._grid('vti_psv').set_phase, 'qX')

    def test_phase_survives_a_pickle(self):
        import pickle
        for aniso in ('vti_psv', 'tti_psv'):
            for phase in ('qP', 'qSV'):
                with self.subTest(aniso=aniso, phase=phase):
                    g = self._grid(aniso, phase)
                    tt = g.raytrace(self.src, self.rcv)
                    # a pickle carries the geometry and the options, not the
                    # medium, which has to be given to the rebuilt grid again
                    g2 = pickle.loads(pickle.dumps(g))
                    g2.set_Vp0(np.full(self.ncells, 3.094))
                    g2.set_Vs0(np.full(self.ncells, 1.51))
                    g2.set_epsilon(np.full(self.ncells, 0.256))
                    g2.set_delta(np.full(self.ncells, -0.0505))
                    if aniso == 'tti_psv':
                        g2.set_tilt_angle(np.full(self.ncells, 0.3))
                    self.assertTrue(np.allclose(
                        tt, g2.raytrace(self.src, self.rcv)))


class TestFSMconvergence(unittest.TestCase):
    """Convergence of the fast sweeping method with the WENO3 stencil.

    The medium is a constant velocity gradient, v(z) = v0 + a z, whose
    traveltime from a point source is known in closed form, so the comparison
    is against truth rather than against another discretisation.

    Regression test for issue #92, where refining the grid made the answer
    *worse*: the sweeps stopped before converging and said nothing about it.
    The 401 x 401 grid is what caught it -- the sequence was still monotone at
    201 x 201 -- so do not drop that resolution to save time.
    """

    v0, a = 2.0, 0.8
    src = np.array([[0.5, 0.2]])
    sizes = (51, 101, 201, 401)
    # the shipped nitermax; a solve that reaches it has been capped, not
    # converged, whatever the traveltimes look like
    default_maxit = 200

    @classmethod
    def _analytic(cls, x, z, vscale):
        """Traveltime from `src` through a constant velocity gradient."""
        v1 = cls.v0 + cls.a*cls.src[0, 1]
        v2 = cls.v0 + cls.a*z
        r2 = (x - cls.src[0, 0])**2 + (z - cls.src[0, 1])**2
        return np.arccosh(1. + cls.a*cls.a*r2/(2.*v1*v2))/cls.a/vscale

    @classmethod
    def _solve(cls, n, vscale=1., **kwargs):
        """Raytrace to every node of an n x n grid.

        `vscale` multiplies the velocities, which divides every traveltime by
        the same factor: the same physics written in different units.

        Returns the grid, so the iteration counts can be inspected, and the
        mean relative error.  The error is measured away from the source,
        since close in the field is set by the initialisation rather than by
        the sweeps under test.
        """
        x = np.linspace(0., 2., n)
        xx, zz = np.meshgrid(x, x, indexing='ij')
        te = cls._analytic(xx, zz, vscale)
        g = rg.Grid2d(x, x, method='FSM', cell_slowness=0, weno=1,
                      n_threads=1, **kwargs)
        g.set_slowness((1./((cls.v0 + cls.a*zz)*vscale)).ravel())
        tt = g.raytrace(cls.src,
                        np.c_[xx.ravel(), zz.ravel()]).reshape(xx.shape)
        far = te > 0.05/vscale
        return g, np.mean(np.abs(tt[far] - te[far])/te[far])

    @classmethod
    def setUpClass(cls):
        # solve each resolution once; three of the tests read these
        cls.solved = {n: cls._solve(n) for n in cls.sizes}

    def test_error_decreases_under_refinement(self):
        errors = [self.solved[n][1] for n in self.sizes]
        for (nc, ec), (nf, ef) in zip(zip(self.sizes, errors),
                                      zip(self.sizes[1:], errors[1:])):
            with self.subTest(coarse=nc, fine=nf):
                self.assertLess(ef, ec,
                                'error grew when the grid was refined')
                # first order in the limit; well clear of the 0.5 floor at
                # every step, and negative before the fix
                self.assertGreater(np.log2(ec/ef), 0.5,
                                   'convergence slower than half an order')

    def test_sweeps_converge_at_the_default_settings(self):
        for n in self.sizes:
            with self.subTest(n=n):
                g = self.solved[n][0]
                self.assertLess(g.get_niterw(), self.default_maxit,
                                'WENO pass hit the iteration cap')
                self.assertLess(g.get_niter(), self.default_maxit,
                                'first-order pass hit the iteration cap')

    @staticmethod
    @contextlib.contextmanager
    def _captured_stderr():
        """Collect what the solver writes to stderr, as a one-element list.

        The non-convergence warning comes from C++ and goes straight to file
        descriptor 2, so redirecting sys.stderr does not see it.
        """
        collected = []
        with tempfile.TemporaryFile(mode='w+') as tmp:
            saved = os.dup(2)
            try:
                os.dup2(tmp.fileno(), 2)
                yield collected
            finally:
                os.dup2(saved, 2)
                os.close(saved)
                tmp.seek(0)
                collected.append(tmp.read())

    def test_a_capped_solve_reports_the_cap(self):
        # a capped solve is not a converged solve, and it used to be silent:
        # the warning and get_niterw() are the only two signals a caller has
        maxit = 5
        with self._captured_stderr() as stderr:
            g, _ = self._solve(101, maxit=maxit)
        self.assertEqual(g.get_niterw(), maxit)
        # the first-order pass converges in a couple of sweeps either way
        self.assertLess(g.get_niter(), maxit)
        self.assertIn('WENO3', stderr[0])
        self.assertIn('nitermax', stderr[0])

    def test_a_converged_solve_is_silent(self):
        with self._captured_stderr() as stderr:
            self._solve(101)
        self.assertEqual(stderr[0], '', 'converged solve warned anyway')

    def test_tolerance_does_not_depend_on_the_model_units(self):
        # epsilon is relative to the traveltime range, so writing the same
        # model in different units must not change where the sweeps stop.
        # Against an absolute tolerance the x1000 case stopped after a single
        # WENO sweep with 8.6 times the error.
        reference = self.solved[201][1]
        for vscale in (1.e3, 1.e-3):
            with self.subTest(vscale=vscale):
                g, err = self._solve(201, vscale=vscale)
                self.assertGreater(g.get_niterw(), 2)
                # exact agreement is out of reach while the WENO smoothness
                # ratio is regularised with the machine epsilon, which a
                # model scaled far enough down runs into; 25% is wide of that
                # and still an order of magnitude tighter than the old failure
                self.assertLess(abs(err - reference)/reference, 0.25,
                                'accuracy changed with the units of the model')

    def test_weno_pass_is_carried_past_the_rising_change(self):
        # The WENO3 change per sweep rises before it falls.  With a tolerance
        # loose enough to sit above the first sweep's change, a loop that only
        # tested the tolerance stopped after one sweep; the guard has to carry
        # it over the peak, which here takes four.  Two would mean only the
        # minimum-iteration floor is left.
        g, _ = self._solve(101, eps=1.e-3)
        self.assertGreater(g.get_niterw(), 2,
                           'WENO pass stopped while the change was rising')


if __name__ == '__main__':

    unittest.main()


class TestRaypathOrder(unittest.TestCase):
    """Raypath coordinates run from source to receiver, whatever the solver.

    See the companion test in test_rgrid3d.py; the steepest-descent builders
    used by FSM and DSPM returned the reverse of the shortest-path order.
    """

    h = 0.05

    def setUp(self):
        self.x = np.arange(41) * self.h
        self.z = self.x.copy()
        X, Z = np.meshgrid(self.x, self.z, indexing='ij')
        self.V = 2.0 + 3.0 * Z
        self.Vc = 2.0 + 3.0 * (Z[:-1, :-1] + 0.5 * self.h)
        self.src = np.array([[0.5, 0.8]])
        self.rcv = np.array([[1.6, 0.05]])

    def test_source_first(self):
        for method, kwargs in (('FSM', {}),
                               ('SPM', dict(nsnx=3, nsnz=3)),
                               ('DSPM', dict(n_secondary=3, n_tertiary=3))):
            for cell_slowness in (0, 1):
                with self.subTest(method=method, cell_slowness=cell_slowness):
                    g = rg.Grid2d(self.x, self.z, n_threads=1, method=method,
                                  cell_slowness=cell_slowness, **kwargs)
                    g.set_slowness((1.0 / (self.Vc if cell_slowness
                                           else self.V)).ravel())
                    r = np.asarray(g.raytrace(self.src, self.rcv,
                                              return_rays=True)[1][0])
                    self.assertGreater(r.shape[0], 2)
                    # both endpoints present exactly; Tx used to be dropped
                    # from the cell-slowness FSM path
                    np.testing.assert_allclose(r[0], self.src[0], atol=1e-9)
                    np.testing.assert_allclose(r[-1], self.rcv[0], atol=1e-9)

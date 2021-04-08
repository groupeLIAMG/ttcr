# -*- coding: utf-8 -*-
"""Tests for verifying python wrappers, module rgrid in 3D"""

import unittest
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

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
    tt = vtk_to_numpy(data.GetPointData().GetArray(name))
    return tt


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

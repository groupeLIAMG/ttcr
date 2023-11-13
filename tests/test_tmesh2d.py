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
        self.r2 = vtk_to_numpy(data.GetCellData().GetArray('r2'))
        self.r4 = vtk_to_numpy(data.GetCellData().GetArray('r4'))

        self.src = np.loadtxt('./files/src2d.dat',skiprows=1)
        # we roll because file has x z t0 and we want t0 x z
        self.src = np.roll(self.src, 1).reshape((1, 3))
        self.rcv = np.loadtxt('./files/rcv2daniso.dat',skiprows=1)

    def test_Mesh2Dsp(self):
        g = tm.Mesh2d(self.nodes, self.tri, method='SPM', n_secondary=10, aniso='weakly_anelliptical')
        g.set_slowness(self.slowness)
        g.set_r2(self.r2)
        g.set_r4(self.r4)
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

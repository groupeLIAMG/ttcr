# -*- coding: utf-8 -*-

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
        reader.SetFileName('layers_coarse.vtr')
        reader.Update()

        data = reader.GetOutput()
        self.x = vtk_to_numpy(data.GetXCoordinates())
        self.y = vtk_to_numpy(data.GetYCoordinates())
        self.z = vtk_to_numpy(data.GetZCoordinates())

        self.slowness = vtk_to_numpy(data.GetCellData().GetArray('Slowness'))
        dim = (self.x.size-1, self.y.size-1, self.z.size-1)
        self.slowness = self.slowness.reshape(dim, order='F').flatten()
        self.src = np.loadtxt('src.dat',skiprows=1)
        self.src = self.src.reshape((1, 4))
        self.rcv = np.loadtxt('rcv.dat',skiprows=1)

    def test_Grid3Dfs(self):
        g = rg.Grid3d(self.x, self.y, self.z, method='FSM', tt_from_rp=False)
        tt = g.raytrace(self.src, self.rcv, self.slowness)
        dim = (self.x.size, self.y.size, self.z.size)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten(order='F')
        tt_ref = get_tt('fsm_d_p_lc_src_all_tt.vtr')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'FSM accuracy failed (slowness in cells)')

    def test_Grid3Dsp(self):
        g = rg.Grid3d(self.x, self.y, self.z, method='SPM', tt_from_rp=False,
                      nsnx=5, nsny=5, nsnz=5)
        tt = g.raytrace(self.src, self.rcv, self.slowness)
        dim = (self.x.size, self.y.size, self.z.size)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten(order='F')
        tt_ref = get_tt('spm_d_p_lc_05_src_all_tt.vtr')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.1,
                        'SPM accuracy failed (slowness in cells)')

    def test_Grid3Ddsp(self):
        g = rg.Grid3d(self.x, self.y, self.z, method='DSPM', tt_from_rp=False,
                      n_secondary=2, n_tertiary=2, radius_tertiary=2)
        tt = g.raytrace(self.src, self.rcv, self.slowness)
        dim = (self.x.size, self.y.size, self.z.size)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten(order='F')
        tt_ref = get_tt('dspm_d_p_lc_2_2_2_src_all_tt.vtr')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.1,
                        'SPM accuracy failed (slowness in cells)')

class TestGrid3dn(unittest.TestCase):

    def setUp(self):
        reader = vtk.vtkXMLRectilinearGridReader()
        reader.SetFileName('gradient_coarse.vtr')
        reader.Update()

        data = reader.GetOutput()
        self.x = vtk_to_numpy(data.GetXCoordinates())
        self.y = vtk_to_numpy(data.GetYCoordinates())
        self.z = vtk_to_numpy(data.GetZCoordinates())

        self.slowness = vtk_to_numpy(data.GetPointData().GetArray('Slowness'))
        dim = (self.x.size, self.y.size, self.z.size)
        self.slowness = self.slowness.reshape(dim, order='F').flatten()
        self.src = np.loadtxt('src.dat',skiprows=1)
        self.src = self.src.reshape((1, 4))
        self.rcv = np.loadtxt('rcv.dat',skiprows=1)

    def test_Grid3Dfs(self):
        g = rg.Grid3d(self.x, self.y, self.z, method='FSM', tt_from_rp=False,
                      cell_slowness=0)
        tt = g.raytrace(self.src, self.rcv, self.slowness)
        dim = (self.x.size, self.y.size, self.z.size)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten(order='F')
        tt_ref = get_tt('fsm_d_p_gc_src_all_tt.vtr')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.01,
                        'FSM accuracy failed (slowness at nodes)')

    def test_Grid3Dsp(self):
        g = rg.Grid3d(self.x, self.y, self.z, method='SPM', tt_from_rp=False,
                      nsnx=5, nsny=5, nsnz=5, cell_slowness=0)
        tt = g.raytrace(self.src, self.rcv, self.slowness)
        dim = (self.x.size, self.y.size, self.z.size)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten(order='F')
        tt_ref = get_tt('spm_d_p_gc_05_src_all_tt.vtr')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.1,
                        'SPM accuracy failed (slowness at nodes)')

    def test_Grid3Ddsp(self):
        g = rg.Grid3d(self.x, self.y, self.z, method='DSPM', tt_from_rp=False,
                      n_secondary=2, n_tertiary=2, radius_tertiary=2,
                      cell_slowness=0)
        tt = g.raytrace(self.src, self.rcv, self.slowness)
        dim = (self.x.size, self.y.size, self.z.size)
        tt = g.get_grid_traveltimes()
        tt = tt.flatten(order='F')
        tt_ref = get_tt('dspm_d_p_gc_2_2_2_src_all_tt.vtr')
        self.assertLess(np.sum(np.abs(tt-tt_ref))/tt.size, 0.1,
                        'SPM accuracy failed (slowness at nodes)')

if __name__ == '__main__':

    unittest.main()

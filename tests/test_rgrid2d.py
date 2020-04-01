# -*- coding: utf-8 -*-

import unittest

import numpy as np
import scipy.sparse as sp

import ttcrpy.rgrid as rg


# class TestGrid2dc(unittest.TestCase):
#
#     def setUp(self):
#         self.x = np.arange(0., 10.1, 0.1)
#         self.z = np.arange(0., 10.1, 0.1)
#
#
#
#     def test_Grid2Dfs(self):
        

class Data_kernel(unittest.TestCase):

    def test_2d(self):

        V = np.ones((11, 13))
        V[:, 7:] = 2
        slowness = 1. / V.flatten()

        grx = np.arange(12.)
        grz = np.arange(14.)

        z = 0.5 + np.arange(13.)
        Tx = np.vstack((0.5+np.zeros((13,)),
                        z)).T
        Rx = np.vstack((10.5+np.zeros((13,)),
                        z)).T
        nTx = Tx.shape[0]
        nRx = Rx.shape[0]
        Tx = np.kron(Tx, np.ones((nRx,1)))
        Rx = np.kron(np.ones((nTx,1)), Rx)

        L = rg.Grid2d.data_kernel_straight_rays(Tx, Rx, grx, grz)
        tt = L.dot(slowness)

        tt2 = np.zeros(tt.shape)
        d = np.sqrt(np.sum((Tx-Rx)**2, axis=1))

        ind = np.logical_and(Tx[:,1]>7, Rx[:,1]>7)
        tt2[ind] = d[ind]/2

        ind2 = np.logical_and(Tx[:,1]<7, Rx[:,1]<7)
        tt2[ind2] = d[ind2]

        ind3 = np.logical_and(np.logical_not(ind), np.logical_not(ind2))

        f = (7-Tx[ind3,1]) / (Rx[ind3,1]-Tx[ind3,1])
        ind = (Rx[ind3,1]-Tx[ind3,1]) < 0
        f[ind] = 1-f[ind]
        tt2[ind3] = d[ind3]*f + d[ind3]*(1-f)/2

        self.assertAlmostEqual(np.sum(np.abs(tt-tt2)), 0.0 )

if __name__ == '__main__':

    unittest.main()

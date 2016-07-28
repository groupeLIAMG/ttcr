# -*- coding: utf-8 -*-

"""
Created on Tue Jun 21 20:55:29 2016

@author: giroux

Copyright 2016 Bernard Giroux
email: bernard.giroux@ete.inrs.ca

This file is part of BhTomoPy.

BhTomoPy is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it /will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np

import cgrid2d

class Grid:
    """
    Superclass for 2D and 3D grids
    """
    
    def __init__(self):
        
        self.grx = np.array([])
        self.gry = np.array([])
        self.grz = np.array([])
        self.cont = None
        self.Tx = np.array([])
        self.Rx = np.array([])
        self.TxCosDir = np.array([])
        self.RxCosDir = np.array([])
        self.bord = np.array([])
        self.Tx_Z_water = 0.0
        self.Rx_Z_water = 0.0
        self.inGrid = np.array([])
    
    def getNumberOfCells(self):
        """
        Returns the number of cells of the grid
        """
        nc = (self.grx.size-1)*(self.grz.size-1)
        if self.gry.size>1:
            nc *= self.gry.size-1
        return nc
    
    def getNcell(self):
        """
        Returns a tuple with the number of cells in each dimension
        """
        if self.gry.size>1:
            return (self.grx.size-1, self.gry.size-1, self.grz.size-1)
        else:
            return (self.grx.size-1, self.grz.size-1)


@staticmethod
    def lsplane(X):
        """
        Least-squares plane (orthogonal distance regression) to a cloud of points
        
        translation of the matlab function lsplane.m by I M Smith
        """
        
        m = X.shape[0]
        if m<3:
            raise ValueError('At least 3 data points required')
        
        x0 = np.mean(X, axis=0)
        
        A = np.vstack((X[:,0]-x0[0], X[:,1]-x0[1], X[:,2]-x0[2])).T
        (U,S,Vh) = np.linalg.svd(A)
        V = Vh.T
        
        i = np.argmin(S)
        a = V[:, i]
        
        # s = np.amin(S)
        # d = U[:,i]*s
        # normd = np.linalg.norm(d)
        # return (x0, a, d, normd)
        
        return (x0, a)


class Grid2D(Grid):
    """
    Class for 2D grids
    """
    
    def __init__(self, grx, grz, nthreads=1):
        Grid.__init__(self)
        self.grx = grx
        self.grz = grz
        self.nthreads = nthreads
        self.nsnx = 10
        self.nsnz = 10
        self.cgrid = None
    
    def raytrace(self, slowness, Tx, Rx, t0=(), xi=(), theta=()):
        """
        Compute traveltimes, raypaths and build ray projection matrix
        
        Usage:
            tt,rays,L = grid.raytrace(slowness,Tx,Rx,t0,xi,theta)
        
        Input:
            slowness: vector of slowness values at grid cells (ncell x 1)
            Tx: coordinates of sources points (ndata x 3)
            Rx: coordinates of receivers      (ndata x 3)
            t0 (optional): initial time at sources points (ndata x 1)
            xi (optional): anisotropy ratio vector ( ncell x 1 )
                values are ratio of slowness in Z over slowness in X
            theta (optional): angle of rotation of the ellipse of anisotropy ( ncell x 1 ),
                counter-clockwise from horizontal, units in radian
        Output:
            tt: vector of traveltimes, ndata by 1
            rays: tuple containing the matrices of coordinates of the ray
                  paths, ndata by 1.  Each matrix is nPts by 2
            L: ray projection matrix, ndata by ncell (ndata x 2*ncell for anisotropic media)
        """
        
        # check input data consistency
        
        if Tx.ndim != 2 or Rx.ndim != 2:
            raise ValueError('Tx and Rx should be 2D arrays')
        
        if Tx.shape[1] !=3 or Rx.shape[1] != 3:
            raise ValueError('Tx and Rx should be ndata x 3')
        
        if Tx.shape != Rx.shape:
            raise ValueError('Tx and Rx should be of equal size')
        
        if len(slowness) != self.getNumberOfCells():
            raise ValueError('Length of slowness vector should equal number of cells')
        
        if len(xi) != 0 and len(xi) != len(slowness):
            raise ValueError('Length of xi should equal length of slowness')
        
        if len(theta) != 0 and len(theta) != len(slowness):
            raise ValueError('Length of theta should equal length of slowness')
        
        if len(t0) == 0:
            t0 = np.zeros([Tx.shape[0],])
        elif len(t0) != Tx.shape[0]:
            raise ValueError('Length of t0 should equal number of Tx')
    
        if self.cgrid == None:
            nx = len(self.grx)-1
            nz = len(self.grz)-1
            dx = self.grx[1]-self.grx[0]
            dz = self.grz[1]-self.grz[0]
            
            typeG = b'iso'
            if len(xi)!=0:
                if len(theta)!=0:
                    typeG = b'tilted'
                else:
                    typeG = b'elliptical'
            self.cgrid = cgrid2d.Grid2Dcpp(typeG,nx,nz,dx,dz,self.grx[0],self.grz[0],self.nsnx,self.nsnz,self.nthreads)
        
                return self.cgrid.raytrace(slowness,xi,theta,Tx,Rx,t0)



if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    
    
    grx = np.linspace(0,10,num=21)
    grz = np.linspace(0,15,num=31)
    
    grid = Grid2D(grx,grz)
    
    slowness = np.ones((grid.getNumberOfCells(),))
    
    Tx = np.array([[0.2, 0.0, 0.2],[0.2, 0.0, 0.2],[0.2, 0.0, 0.2],[0.2, 0.0, 1.2],[0.2, 0.0, 1.2],[0.2, 0.0, 1.2]])
    Rx = np.array([[9.8, 0.0, 0.2],[9.8, 0.0, 3.2],[9.8, 0.0, 6.2],[9.8, 0.0, 0.2],[9.8, 0.0, 3.2],[9.8, 0.0, 6.2]])
    t0 = np.zeros([6,])
    
    tt1,rays1,L1 = grid.raytrace(slowness, Tx, Rx, t0)
    
    tt2,rays2,L2 = grid.raytrace(slowness, Tx, Rx)
    
    d = np.sqrt(np.sum((Tx-Rx)**2,axis=1))
    
    print(d-tt1)
    print(d-tt2)   
    
    print(L1)
    
    plt.matshow(L1.toarray())
    
    plt.show()
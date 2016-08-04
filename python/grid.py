# -*- coding: utf-8 -*-

"""
Created on Tue Jun 21 20:55:29 2016

@author: giroux

Copyright 2016 Bernard Giroux
email: bernard.giroux@ete.inrs.ca

This program is free software: you can redistribute it and/or modify
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

import math
import numpy as np

import cgrid2d

from utils import nargout


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
            
    def dx(self):
        return self.grx[1]-self.grx[0]
        
    def dy(self):
        if len(self.gry)==0:
            return 0
        else:
            return self.gry[1]-self.gry[0]
        
    def dz(self):
        return self.grz[1]-self.grz[0]
            


    @staticmethod
    def lsplane(X):
        """
        Least-squares plane (orthogonal distance regression) to a cloud of points
        Usages:
            x0,a = lsplane(x)
            x0,a,d,normd = lsplane(x)
            
        Input:
            x: cloud of points (n x 3)
            
        Output:
            x0: point on plane (3,)
            a: direction cosine of normal to plane (1x3)
            d: residuals
            normd: norm of residuals
            
        translation of the matlab function lsplane.m by I M Smith
        """
        
        nout = nargout()
        
        m = X.shape[0]
        if m<3:
            raise ValueError('At least 3 data points required')
        
        x0 = np.mean(X, axis=0).T
        
        A = np.vstack([X[:,0]-x0[0], X[:,1]-x0[1], X[:,2]-x0[2]]).T
        (U,S,V) = np.linalg.svd(A)
        
        i = np.argmin(S)
        a = V[i,:].T
        
        if nout==4:
            s = np.amin(S)
            d = U[:,i]*s
            normd = np.linalg.norm(d)
            return (x0, a, d, normd)
        elif nout==2:
            return (x0, a)


    @staticmethod
    def boreholes_order(bh):
        """
        Compute order of many boreholes so that they best align
        
        Usage:
            boreholes_order(bh) -> order
            
        Assumption: boreholes are more or less vertical
        """
        nd = len(bh)
        x = [ b.X for b in bh ]
        y = [ b.Y for b in bh ]
        
        dx = max(x)-min(x)
        dy = max(y)-min(y)
        
        # find along which dimension we should work
        if dx>dy:
            x = np.array(x)
            y = np.array(y)
        else:
            tmp = x
            x = np.array(y)
            y = np.array(tmp)
            
        order = np.arange(nd)
        ix = np.argsort(x)
        x = x[ix]
        y = y[ix]
        order = order[ix]
        
        # check if xmin is repeated
        ind = np.equal(x[0], x)
        if ind.sum()>1:
            iy = np.argsort(y[ind])
            y[ind] = y[iy]
            order[ind] = order[iy]
        
        # sort sequentially according to distance from previous bh
        for n in np.arange(nd-2):
            dist = np.sqrt( (x[n]-x[n+1:])**2 + (y[n]-y[n+1:])**2)
            ind = np.argsort(dist)
            tmp = x[n+1:]
            x[n+1:] = tmp[ind]
            tmp = y[n+1:]
            y[n+1:] = tmp[ind]
            tmp = order[n+1:]
            order[n+1:] = tmp[ind]

        return order
        
        
    @staticmethod
    def proj_plane(data, x0, a):
        """
        p_data = proj_plane(data, x0, a)
        
        Input:
            data : coordinates to project on plane (nx3)
            x0   : pt on the plane (3,)
            a    : direction cosine of normal to plane (3,)

        Output:
            p_data : projected coordinates (nx3)
        """
        p_data = data.copy()
        for n in np.arange(data.shape[0]):
            r = x0 - data[n,:]
            p = np.dot(a,r)
            p_data[n,:] = data[n,:] + p*a
        
        return p_data
    
    
    @staticmethod
    def proj_planes(data, planes):
        """
        p_data,no_plane = proj_planes(data, planes)
        
        Input:
            data : coordinates to project on plane (nx3)
            planes : list of plane objects with the following attributes
                x0   : pt on the plane (3,)
                a    : direction cosine of normal to plane (3,)
        
        Output:
            p_data : projected coordinates (nx3)
            no_plane : index of plane on which data is projected (n,)
        """
        p_data = data.copy()
        no_plane = np.zeros((data.shape[0],),dtype=int)
        p = np.zeros((len(planes),))
        for n in np.arange(data.shape[0]):
            for nn in np.arange(len(planes)):
                r = planes[nn].x0 - data[n,:]
                p[nn] = np.dot(planes[nn].a,r)
                
            no = np.argmin(np.abs(p))  # find plane closest to current point
            p_data[n,:] = data[n,:] + p[no]*planes[no].a  # project on closest plane
            no_plane[n] = no

        return p_data,no_plane
        
        
    @staticmethod
    def transl_rotat(data, origin, az, dip):
        """
        Translate with respect to origin and rotate w/r to azimuth and dip
        
        m_data = transl_rotat(data, origin, az, dip)
        
        """
        
        # translation w/r to origin        
        m_data = np.vstack([data[:,0]-origin[0], data[:,1]-origin[1], data[:,2]-origin[2]]).T
        # rotation w/r to azimuth
        if np.abs(az)>(math.pi/720.0):  # if azimuth larger that 1/4 of a degree
            rot = np.array([[np.cos(az), -np.sin(az)],[np.sin(az), np.cos(az)]])
            for n in np.arange(m_data.shape[0]):
                m_data[n,:2] = np.dot(m_data[n,:2], rot.T)
                # np.dot(m_data[n,:2], rot.T) is equal to np.dot(rot, m_data[n,:2].T).T
                
        # rotation w/r to dip
        if np.abs(dip)>(math.pi/720.0):  # if azimuth larger that 1/4 of a degree
            rot = np.array([[np.cos(dip), -np.sin(dip)],[np.sin(dip), np.cos(dip)]])
            for n in np.arange(m_data.shape[0]):
                m_data[n,1:] = np.dot(m_data[n,1:], rot.T)
        
        return m_data
        
        
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
        
        Usages:
            tt,L,rays = grid.raytrace(slowness,Tx,Rx,t0,xi,theta)
            tt,L = grid.raytrace(slowness,Tx,Rx,t0,xi,theta)  {Note: rays can be omitted from output}
        
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
            L: ray projection matrix, ndata by ncell (ndata x 2*ncell for anisotropic media)
            rays: tuple containing the matrices of coordinates of the ray
                  paths, ndata by 1.  Each matrix is nPts by 2
        """
        
        nout = nargout()
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
        
        if nout==2:
            tt,L = self.cgrid.raytrace(slowness,xi,theta,Tx,Rx,t0)
            return tt,L
        elif nout==3:
            tt,L,rays = self.cgrid.raytrace(slowness,xi,theta,Tx,Rx,t0)
            return tt,L,rays
            

    def getForwardStraightRays(self, ind=None, dx=None, dy=None, dz=None, aniso=False):
        """
        Build ray projection matrix for straight rays
        
        Input:
            ind: indices of Tx-Rx pairs for which matrix is built
            dx: grid cell size along X (default is size of grid instance)
            dz: grid cell size along Z (default is size of grid instance)
            aniso: if true build matrix for anisotropic slowness
            
        Output:
            L: ray projection matrix, ndata by ncell (ndata x 2*ncell for anisotropic media)
        """
        if ind is None:
            ind = np.ones((self.Tx.shape[0],),dtype=bool)

        small = 0.00001            
        if dx is None:
            grx = self.grx
        else:
            grx = np.arange(self.grx[0],self.grx[-1]+small, dx)  
        
        if dz is None:
            grz = self.grz
        else:
            grz = np.arange(self.grz[0],self.grz[-1]+small, dz)

        if aniso==False:
            return cgrid2d.Grid2Dcpp.Lsr2d(self.Tx[np.ix_(ind,[0,2])], self.Rx[np.ix_(ind,[0,2])], grx, grz)
        else:
            return cgrid2d.Grid2Dcpp.Lsr2da(self.Tx[np.ix_(ind,[0,2])], self.Rx[np.ix_(ind,[0,2])], grx, grz)
            
            

if __name__ == '__main__':
    
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    
    testRaytrace = False
    testStatic = True
    
    if testRaytrace:
        grx = np.linspace(0,10,num=21)
        grz = np.linspace(0,15,num=31)
        
        grid = Grid2D(grx,grz)
        
        nc = grid.getNumberOfCells()
        slowness = np.ones((nc,))
        
        Tx = np.array([[0.2, 0.0, 0.2],[0.2, 0.0, 0.2],[0.2, 0.0, 0.2],[0.2, 0.0, 1.2],[0.2, 0.0, 1.2],[0.2, 0.0, 1.2]])
        Rx = np.array([[9.8, 0.0, 0.2],[9.8, 0.0, 3.2],[9.8, 0.0, 6.2],[9.8, 0.0, 0.2],[9.8, 0.0, 3.2],[9.8, 0.0, 6.2]])
        t0 = np.zeros([6,])
        
        tt1,L1,rays1 = grid.raytrace(slowness, Tx, Rx, t0)
        tt1b,L1b = grid.raytrace(slowness, Tx, Rx, t0)    
        tt2,L2,rays2 = grid.raytrace(slowness, Tx, Rx)
        
        d = np.sqrt(np.sum((Tx-Rx)**2,axis=1))
        
        print(d-tt1)
        print(d-tt1b)
        
        print(d-tt2)   
        
        dL = L1-L1b
        print( dL==0 )
        
        plt.matshow(L1.toarray())
        
        plt.show()
        
        
        grid.Tx = Tx
        grid.Rx = Rx
        
        Lsr = grid.getForwardStraightRays()
        ttsr = Lsr*slowness
        
        print(d-ttsr)
    
    
        Lsr2 = grid.getForwardStraightRays(aniso=True)
        slowness2 = np.ones((2*nc,))
        ttsr2 = np.sqrt( (Lsr2[:,:nc]*slowness2[:nc])**2 + (Lsr2[:,nc:]*slowness2[nc:])**2)
        
        print(d-ttsr2)
    
        
        plt.matshow(Lsr.toarray())
        
        plt.show()
        
    if testStatic:
        
        class Borehole:
            X = 0
            Y = 0
        
        bh1 = Borehole()
        bh2 = Borehole()
        bh3 = Borehole()
        bh4 = Borehole()
        
        bh3.X = 4.0
        bh3.Y = 1.0
        bh2.X = 7.6
        bh2.Y = 0.67
        bh1.Y = 1.0
        
        bh=[bh1, bh2, bh3, bh4]
        
        order = Grid.boreholes_order(bh)
        
        print(order)
        
        
        import numpy.random
        x = numpy.random.rand(100,3)
        x[:,0] *= 30
        x[:,1] *= 25
        x[:,2] -= 0.5
        
        x0,a = Grid.lsplane(x)
        
        xp = Grid.proj_plane(x,x0,a)
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        ax.scatter(xs=x[:,0],  ys=x[:,1],  zs=x[:,2],  c='r')
        ax.scatter(xs=xp[:,0], ys=xp[:,1], zs=xp[:,2], c='b')        
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        
        ax.view_init(0,45)

        plt.show()

        class Plane:
            x0 = 0
            a = 0
        
        plane1 = Plane()
        plane1.a = a
        plane1.x0 = x0
        
        x = numpy.random.rand(100,3)
        x[:,0] *= 30
        x[:,1] *= 25
        x[:,2] -= 0.25
        
        x0,a = Grid.lsplane(x)
        plane2 = Plane()
        plane2.a = a
        plane2.x0 = x0
        
        planes = [plane1, plane2]
        
        xp,np_plane = Grid.proj_planes(x, planes)
        
        
        origin = np.array([np.min(xp[:,0]), np.min(xp[:,1]), 0])
        az = 30*3.14159/180
        dip = -45*3.14159/180
        
        
        xd = Grid.transl_rotat(xp, origin, az, dip)
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        ax.scatter(xs=xd[:,0], ys=xd[:,1], zs=xd[:,2],  c='r')
        ax.scatter(xs=xp[:,0], ys=xp[:,1], zs=xp[:,2], c='b')        
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        
        ax.view_init(0,45)

        plt.show()
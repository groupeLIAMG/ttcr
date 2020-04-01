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
from scipy.sparse import csr_matrix
import h5py

from cutils import cgrid2d # @UnresolvedImport

from utils import nargout
import covar

class Grid:
    """
    Superclass for 2D and 3D grids

    Grids are regular, i.e. constant step size, although not restricted to
        square or cubic cells (dx, dy and dz may differ from each other)
    """

    def __init__(self):

        self.grx = np.array([])  # coordinates of grid nodes along X
        self.gry = np.array([])  # coordinates of grid nodes along Y
        self.grz = np.array([])  # coordinates of grid nodes along Z
        self.cont = None
        self.Tx = np.array([])
        self.Rx = np.array([])
        self.TxCosDir = np.array([])
        self.RxCosDir = np.array([])
        self.border = np.array([])
        self.Tx_Z_water = np.nan
        self.Rx_Z_water = np.nan
        self.in_vect = np.array([])

    def getNumberOfCells(self):
        """
        Returns the number of cells of the grid
        """
        nc = (self.grx.size-1)*(self.grz.size-1)
        if self.gry.size > 1:
            nc *= self.gry.size-1
        return nc

    def getNcell(self):
        """
        Returns a tuple with the number of cells in each dimension
        """
        if self.gry.size > 1:
            return (self.grx.size-1, self.gry.size-1, self.grz.size-1)
        else:
            return (self.grx.size-1, self.grz.size-1)

    @property
    def dx(self):
        return self.grx[1] - self.grx[0]

    @property
    def dy(self):
        if len(self.gry) == 0:
            return 0
        else:
            return self.gry[1] - self.gry[0]

    @property
    def dz(self):
        return self.grz[1] - self.grz[0]



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
        if m < 3:
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

    Important: the raytracing codes are based on a column-major order
            for the slowness vector (Z is the "fast" axis).
            To visualize the slowness model with Z axis vertical and X horizontal,
            the vector should be reshaped as
            slowness.reshape(nx,nz).T

    """

    def __init__(self, grx=None, grz=None, nthreads=1):
        Grid.__init__(self)
        if grx is not None:
            self.grx = grx
        if grz is not None:
            self.grz = grz
        self.nthreads = nthreads
        self.nsnx = 10
        self.nsnz = 10
        self.cgrid = None
        self.border = np.array([1, 1, 1, 1])
        self.flip = 0
        self.borehole_x0 = 1
        self.x0 = np.array([])
        self.type = None

    def __reduce__(self):
        # cgrid excluded volontarily, it will be set to None after unpickling
        # (cgrid is instantiated when needed in method raytrace)
        # this is done to avoid writing code to pickle cython class Grid2Dcpp
        return (Grid2D.rebuild,(self.grx, self.grz, self.cont, self.Tx, self.Rx,
                                self.TxCosDir, self.RxCosDir, self.border,
                                self.Tx_Z_water, self.Rx_Z_water, self.in_vect,
                                self.nthreads, self.nsnx, self.nsnz, self.flip,
                                self.borehole_x0, self.x0, self.type))

    @staticmethod
    def rebuild(grx, grz, cont, Tx, Rx, TxCosDir, RxCosDir, border, Tx_Z_water,
                Rx_Z_water, in_vect, nthreads, nsnx, nsnz, flip, borehole_x0, x0, _type):

        g = Grid2D(grx, grz, nthreads)

        g.cont = cont
        g.Tx = Tx
        g.Rx = Rx
        g.TxCosDir = TxCosDir
        g.RxCosDir = RxCosDir
        g.border = border
        g.Tx_Z_water = Tx_Z_water
        g.Rx_Z_water = Rx_Z_water
        g.in_vect = in_vect
        g.nsnx = nsnx
        g.nsnz = nsnz
        g.flip = flip
        g.borehole_x0 = borehole_x0
        g.x0 = x0
        g.type = _type

        return g


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
            self.cgrid = cgrid2d.Grid2Dcpp(typeG,nx,nz,dx,dz,self.grx[0],self.grz[0], #  @UndefinedVariable
                                           self.nsnx,self.nsnz,self.nthreads)

        if nout==2:
            tt,L = self.cgrid.raytrace(slowness,xi,theta,Tx,Rx,t0,nout)
            return tt,L
        elif nout==3:
            tt,L,rays = self.cgrid.raytrace(slowness,xi,theta,Tx,Rx,t0,nout)
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
            return cgrid2d.Grid2Dcpp.Lsr2d(self.Tx[np.ix_(ind,[0,2])], self.Rx[np.ix_(ind,[0,2])], grx, grz) #  @UndefinedVariable
        else:
            return cgrid2d.Grid2Dcpp.Lsr2da(self.Tx[np.ix_(ind,[0,2])], self.Rx[np.ix_(ind,[0,2])], grx, grz) #  @UndefinedVariable


    def getCellCenter(self, dx=None, dz=None):
        """
        Returns a nCell x 2 array containing the coordinates of the center of the cells
        """
        if dx is None:
            dx = self.grx[1] - self.grx[0]
        if dz is None:
            dz = self.grz[1] - self.grz[0]

        xmin = self.grx[0] + dx/2.0
        zmin = self.grz[0] + dz/2.0
        xmax = self.grx[-1] - dx/3.0  # divide by 3 to avoid truncation error
        zmax = self.grz[-1] = dz/3.0
        nx = np.ceil((xmax-xmin)/dx)
        nz = np.ceil((zmax-zmin)/dz)

        c = np.vstack( [xmin+np.kron(np.ones((nz,)), np.arange(nx)*dx),
                        zmin+np.kron(np.arange(nz), np.ones((nx,))*dz)] ).T
        return c


    def checkCenter(self, x, y, z):
        """
        Verify if given coordinates correspond to the center of the cells
        """
        if len(y)>0:
            return False    # Y should be empty for 2D grids

        dx = self.grx[1] - self.grx[0]
        xmin = self.grx[0] + dx/2.0
        xmax = self.grx[-1] - dx/3.0  # divide by 3 to avoid truncation error
        xx = np.arange(xmin, xmax, dx)
        if xx.size() != x.size():
            return False

        if np.any(np.abs(x-xx)>1000.0*np.finfo(float).eps):
            return False

        dz = self.grz[1] - self.grz[0]
        zmin = self.grz[0] + dz/2.0
        zmax = self.grz[-1] - dz/3.0  # divide by 3 to avoid truncation error
        zz = np.arange(zmin, zmax, dz)
        if zz.size() != z.size():
            return False

        if np.any(np.abs(z-zz)>1000.0*np.finfo(float).eps):
            return False

        return True

    def derivative(self, order, normalize = False):
        """
        Compute spatial derivative operators for grid _cells_

        For 1st order:
            forward operator is (u_{i+1} - u_i)/dx
            centered operator is (u_{i+1} - u_{i-1})/(2dx)
            backward operator is (u_i - u_{i-1})/dx

        For 2nd order:
            forward operator is (u_i - 2u_{i+1} + u_{i+2})/dx^2
            centered operator is (u_{i-1} - 2u_i + u_{i+1})/dx^2
            backward operator is (u_{i-2} - 2u_{i-1} + u_i)/dx^2
        """
        dx = 1
        dz = 1
        if normalize:
            dx = self.dx()
            dz = self.dz()

        nx = len(self.grx) - 1
        nz = len(self.grz) - 1

        if order == 1:

            # forward operator is (u_{i+1} - u_i)/dx
            # centered operator is (u_{i+1} - u_{i-1})/(2dx)
            # backward operator is (u_i - u_{i-1})/dx

            idx = 1/dx
            idz = 1/dz

            i = np.kron(np.arange(nx*nz),np.ones((2,)))
            j = np.zeros((nz * nx * 2,))
            v = np.zeros((nz * nx * 2,))


            jj = np.vstack((np.arange(nz),nz+np.arange(nz))).T
            jj = jj.flatten()
            j[:2*nz] = jj
            vd = idx * np.tile(np.array([-1,1]),(nz,))
            v[:2*nz] = vd

            jj = np.vstack((-nz+np.arange(nz),nz+np.arange(nz))).T
            jj = jj.flatten()
            for n in range(1,nx-1):
                j[n*2*nz:(n+1)*2*nz] = n*nz + jj
                v[n*2*nz:(n+1)*2*nz] = 0.5*vd

            jj = np.vstack((-nz+np.arange(nz),np.arange(nz))).T
            jj = jj.flatten()
            j[(nx-1)*2*nz:nx*2*nz] = (nx-1)*nz + jj
            v[(nx-1)*2*nz:nx*2*nz] = vd

            Dx = csr_matrix((v,(i,j)))


            jj = np.vstack((np.hstack((0,np.arange(nz-1))),
                            np.hstack((np.arange(1,nz), nz-1)))).T
            jj = jj.flatten()
            vd = idz*np.hstack((np.array([-1,1]),
                                np.tile(np.array([-0.5,0.5]),(nz-2,)),np.array([-1,1])))

            for n in range(nx):
                j[n*2*nz:(n+1)*2*nz] = n*nz + jj
                v[n*2*nz:(n+1)*2*nz] = vd

            Dz = csr_matrix((v,(i,j)))
        else:  # 2nd order

            # forward operator is (u_i - 2u_{i+1} + u_{i+2})/dx^2
            # centered operator is (u_{i-1} - 2u_i + u_{i+1})/dx^2
            # backward operator is (u_{i-2} - 2u_{i-1} + u_i)/dx^2

            idx2 = 1/(dx*dx)
            idz2 = 1/(dz*dz)

            i = np.kron(np.arange(nx*nz),np.ones((3,)))
            j = np.zeros((nz * nx * 3,))
            v = np.zeros((nz * nx * 3,))

            jj = np.vstack((np.arange(nz),nz+np.arange(nz),2*nz+np.arange(nz))).T
            jj = jj.flatten()
            j[:3*nz] = jj
            vd = idx2*np.tile(np.array([1.0,-2.0,1.0]),(nz,))
            v[:3*nz] = vd

            for n in range(1,nx-1):
                j[n*3*nz:(n+1)*3*nz] = (n-1)*nz + jj
                v[n*3*nz:(n+1)*3*nz] = vd

            j[(nx-1)*3*nz:nx*3*nz] = (nx-3)*nz + jj
            v[(nx-1)*3*nz:nx*3*nz] = vd

            Dx = csr_matrix((v,(i,j)))


            jj = np.vstack((np.hstack((0,np.arange(nz-2),nz-3)),
                            np.hstack((1,np.arange(1,nz-1), nz-2)),
                            np.hstack((2,np.arange(2,nz), nz-1)))).T
            jj = jj.flatten()
            vd = vd*idz2/idx2

            for n in range(nx):
                j[n*3*nz:(n+1)*3*nz] = n*nz + jj
                v[n*3*nz:(n+1)*3*nz] = vd

            Dz = csr_matrix((v,(i,j)))

        # no derivative along y, empty Dy matrix
        Dy = csr_matrix(Dx.shape)

        return Dx,Dy,Dz

    def preFFTMA(self, cm):
        """
        Compute matrix G for FFT-MA simulations

        INPUT
            cm: list of covariance models

        OUTPUT
            G: covariance matrix in spectral domain
        """
        small = 1.0e-6
        Nx = 2*self.grx.size
        Nz = 2*self.grz.size

        Nx2 = Nx/2
        Nz2 = Nz/2

        x = self.dx * np.hstack((np.arange(Nx2), np.arange(-Nx2+1,1)))
        z = self.dz * np.hstack((np.arange(Nz2), np.arange(-Nz2+1,1)))

        x = np.kron(x,np.ones((Nz,)))
        z = np.kron(z,np.ones((1,Nx)).T).flatten()

        d = 0
        for c in cm:
            d = d + c.compute(np.vstack((x,z)).T, np.zeros((1,2)))
        K = d.reshape(Nx,Nz)

        mk = True
        while mk:
            mk = False
            if np.min(K[0,:])>small:
                # Enlarge grid to make sure that covariance falls to zero
                Nz = 2*Nz
                mk = True

            if np.min(K[:,0])>small:
                Nx = 2*Nx
                mk = True

            if mk:
                Nx2 = Nx/2
                Nz2 = Nz/2

                x = self.dx * np.hstack((np.arange(Nx2), np.arange(-Nx2+1,1)))
                z = self.dz * np.hstack((np.arange(Nz2), np.arange(-Nz2+1,1)))

                x = np.kron(x,np.ones((Nz,)))
                z = np.kron(z,np.ones((1,Nx)).T).flatten()

                d = 0
                for c in cm:
                    d = d + c.compute(np.vstack((x,z)).T, np.zeros((1,2)))
                K = d.reshape(Nx,Nz)

        return np.sqrt(np.fft.fft2(K))

    def FFTMA(self, G):
        """
        Perform FFT-MA simulation using pre-computed spectral matrix

        INPUT
            G: covariance matrix in spectral domain as return by preFFTMA

        OUTPUT
            Z: simulated field of size nx x nz
        """
        Nx,Nz = G.shape
        U = np.random.randn(G.shape[0], G.shape[1])
        Z = np.real(np.fft.ifft2(G*U))

        return Z[int(round((Nx+2)/2)):(int(round((Nx+2)/2))+self.grx.size-1),
                 int(round((Nz+2)/2)):(int(round((Nz+2)/2))+self.grz.size-1)]

    def toXdmf(self, field, fieldname, filename):
        """
        Save a field in xdmf format

        INPUT
            field: data array of size equal to the number of cells in the grid
            fieldname: name to be assinged to the data (string)
            filename: name of xdmf file (string)
        """
        nx = self.grx.size-1
        nz = self.grz.size-1
        ox = self.grx[0] + self.dx/2
        oz = self.grz[0] + self.dz/2

        f = open(filename,'w')

        f.write('<?xml version="1.0" ?>\n')
        f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
        f.write('<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">\n')
        f.write(' <Domain>\n')
        f.write('   <Grid Name="Structured Grid" GridType="Uniform">\n')
        f.write('     <Topology TopologyType="2DCORECTMesh" NumberOfElements="'+
        repr(nz+1)+' '+repr(nx+1)+'"/>\n')
        f.write('     <Geometry GeometryType="ORIGIN_DXDY">\n')
        f.write('       <DataItem Dimensions="2 " NumberType="Float" Precision="4" Format="XML">\n')
        f.write('          '+repr(oz)+' '+repr(ox)+'\n')
        f.write('       </DataItem>\n')
        f.write('       <DataItem Dimensions="2 " NumberType="Float" Precision="4" Format="XML">\n')
        f.write('        '+repr(self.dz)+' '+repr(self.dx)+'\n')
        f.write('       </DataItem>\n')
        f.write('     </Geometry>\n')
        f.write('     <Attribute Name="'+fieldname+'" AttributeType="Scalar" Center="Cell">\n')
        f.write('       <DataItem Dimensions="'+repr(nz)+' '+repr(nx)+
        '" NumberType="Float" Precision="4" Format="HDF">'+filename+'.h5:/'+fieldname+'</DataItem>\n')
        f.write('     </Attribute>\n')
        f.write('   </Grid>\n')
        f.write(' </Domain>\n')
        f.write('</Xdmf>\n')

        f.close()

        h5f = h5py.File(filename+'.h5', 'w')
        h5f.create_dataset(fieldname, data=field.reshape((nx,nz)).T.astype(np.float32))
        h5f.close()

if __name__ == '__main__':

    from mpl_toolkits.mplot3d import Axes3D # @UnresolvedImport
    import matplotlib.pyplot as plt

    testRaytrace = False
    testRaytrace2 = False
    testStatic = False
    testDeriv = False
    testFFTMA = True
    testPickle = False

    if testRaytrace:
        grx = np.linspace(0,10,num=21)
        grz = np.linspace(0,15,num=31)

        grid = Grid2D(grx,grz)

        print(grid.dx)

        nc = grid.getNumberOfCells()
        slowness = np.ones((nc,))

        Tx = np.array([[0.2, 0.0, 0.2],
                       [0.2, 0.0, 0.2],
                       [0.2, 0.0, 0.2],
                       [0.2, 0.0, 1.2],
                       [0.2, 0.0, 1.2],
                       [0.2, 0.0, 1.2]])
        Rx = np.array([[9.8, 0.0, 0.2],
                       [9.8, 0.0, 3.2],
                       [9.8, 0.0, 6.2],
                       [9.8, 0.0, 0.2],
                       [9.8, 0.0, 3.2],
                       [9.8, 0.0, 6.2]])
        t0 = np.zeros([6,])

        tt1,L1,rays1 = grid.raytrace(slowness, Tx, Rx, t0)
        tt1b,L1b = grid.raytrace(slowness, Tx, Rx, t0)
        tt2,L2,rays2 = grid.raytrace(slowness, Tx, Rx)

        d = np.sqrt(np.sum((Tx-Rx)**2,axis=1))

        print(d-tt1)
        print(d-tt1b)

        print(d-tt2)
        tt2b = L2*slowness
        print(d-tt2b)

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

    if testRaytrace2:
        grx = np.linspace(0,10,num=21)
        grz = np.linspace(0,15,num=31)

        grid = Grid2D(grx,grz)

        nc = grid.getNumberOfCells()
        slowness = np.ones((nc,))
        slowness = slowness.reshape((20,30))
        slowness[:,:10] = 0.5

        plt.matshow(slowness.T)
        plt.show()

        s = slowness.flatten()
        print(s[8:13])

        z = np.arange(1,15)

        Tx = np.vstack((np.ones(z.size), np.zeros(z.size),z)).T
        Rx = np.vstack((9+np.ones(z.size), np.zeros(z.size),z)).T

        grid.Tx = Tx
        grid.Rx = Rx

        Lsr = grid.getForwardStraightRays()
        ttsr = Lsr*s

        print(ttsr)

        tt1,L1,rays1 = grid.raytrace(s, Tx, Rx)



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

    if testDeriv:
        grx = np.linspace(0,3,num=8)
        grz = np.linspace(0,6,num=9)
        nx=len(grx)-1
        nz=len(grz)-1

        grid = Grid2D(grx,grz)

        Dx,Dy,Dz = grid.derivative(1)

        plt.matshow(Dx.toarray())
        plt.title('Dx')
        plt.show()

        nc = grid.getNumberOfCells()
        s = np.ones((nc,))
        s[0] = 2
        s[27] = 2
        s[-1] = 2

        dsx = Dx*s
        dsz = Dz*s

        plt.matshow(s.reshape((nx,nz)).T)
        plt.title('s')
        plt.colorbar()
        plt.show()

        plt.matshow(dsx.reshape((nx,nz)).T)
        plt.title('dsx')
        plt.colorbar()
        plt.show()

        plt.matshow(Dz.toarray())
        plt.title('Dz')
        plt.show()

        plt.matshow(dsz.reshape((nx,nz)).T)
        plt.title('dsz')
        plt.colorbar()
        plt.show()



        Dx,Dy,Dz = grid.derivative(2)

        plt.matshow(Dx.toarray())
        plt.title('Dx O2')
        plt.show()

        dsx = Dx*s
        dsz = Dz*s

        plt.matshow(dsx.reshape((nx,nz)).T)
        plt.title('dsx O2')
        plt.colorbar()
        plt.show()

        plt.matshow(Dz.toarray())
        plt.title('Dz O2')
        plt.show()

        plt.matshow(dsz.reshape((nx,nz)).T)
        plt.title('dsz O2')
        plt.colorbar()
        plt.show()

    if testFFTMA:
        grx = np.linspace(0,3,num=50)
        grz = np.linspace(0,6,num=100)
        nx=len(grx)-1
        nz=len(grz)-1

        grid = Grid2D(grx,grz)

        cm = [covar.CovarianceExponential(np.array([5.0,2.0]), np.array([0]), 2.5)]

        G = grid.preFFTMA(cm)

        Z = grid.FFTMA(G)

        plt.matshow(Z.T)
        plt.show()

        grid.toXdmf(Z,'fftma','test.xmf')

    if testPickle:
        import pickle

        grx = np.linspace(0,10,num=21)
        grz = np.linspace(0,15,num=31)

        grid = Grid2D(grx,grz)

        z = np.arange(1,15)

        Tx = np.vstack((np.ones(z.size), np.zeros(z.size),z)).T
        Rx = np.vstack((9+np.ones(z.size), np.zeros(z.size),z)).T

        grid.Tx = Tx
        grid.Rx = Rx

        nc = grid.getNumberOfCells()
        slowness = np.ones((nc,))
        slowness = slowness.reshape((20,30))
        slowness[:,:10] = 0.5

        s = slowness.flatten()


        Lsr = grid.getForwardStraightRays()
        ttsr = Lsr*s

        with open('/tmp/data.pickle', 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(grid, f, pickle.HIGHEST_PROTOCOL)


        del grid

        with open('/tmp/data.pickle', 'rb') as f:
            # The protocol version used is detected automatically, so we do not
            # have to specify it.
            grid = pickle.load(f)



        Lsr1 = grid.getForwardStraightRays()
        ttsr1 = Lsr1*s

        tt1,L1,rays1 = grid.raytrace(s, Tx, Rx)


        with open('/tmp/data.pickle', 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(grid, f, pickle.HIGHEST_PROTOCOL)


        del grid

        with open('/tmp/data.pickle', 'rb') as f:
            # The protocol version used is detected automatically, so we do not
            # have to specify it.
            grid = pickle.load(f)

        Lsr2 = grid.getForwardStraightRays()
        ttsr2 = Lsr2*s

        tt2,L2,rays2 = grid.raytrace(s, Tx, Rx)

        print(ttsr)
        print(ttsr1)
        print(ttsr2)
        print(tt1)
        print(tt2)

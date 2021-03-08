#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 08:44:16 2019

@author: giroux
"""

import numpy as np
import scipy.optimize as spo

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



def fun(xi, h, v, x):
    f = np.zeros(h.shape)
    for n in np.arange(h.size-1):
        f[n] = xi[n] * v[n+1] * np.sqrt(h[n+1]**2 + xi[n+1]**2) - \
            xi[n+1] * v[n] * np.sqrt(h[n]**2 + xi[n]**2)
    f[-1] = np.sum(xi)-x
    return f

def solve_nlayers(h, V, x):
    '''
    Calcul le temps dans un milieu 1D, Ã  un offset x au mur de la couche la
      plus basse

      input
        h: Épaisseurs des couches
        V: vitesse d'intervalle des couches
        x: deport (offset) sur la surface la plus profonde

      output
        t : temps
        rl: longueur du rai
        xi: distances horizontales pour chaque couche
    '''
    # starting value is the average
    xi0 = np.ones(h.shape) * x/h.size

    xi = spo.fsolve(fun, xi0, (h, V, x))

    l = np.sqrt( xi**2 + h**2 )
    t = np.sum(l/V)
    rl = np.sum(l)
    return t, rl, xi



def clean_nan(x):
    i = np.logical_not(np.isnan(x))
    x = x[i]
    x = np.reshape(x, (int(x.size/3 +0.000000001), 3))
    return x


def ind(i,j,k):
    return (i*(N+1)+j)*(N+1) + k


# %%

N = 20
x = np.arange(N+1)
y = np.arange(N+1)
z = np.arange(N+1)

# Vitesse a z=0
a = 1      # (km/s)
Vmax = 3   # vitesse a 20 km  (Rapport 1:3, tel que Podvin & Lecomte, 1991, fig 13)
b = (Vmax-a)/z[-1]

# Source a (0,0,0)
x0 = 0
y0 = 0
z0 = 0

zc = 0.5+z
V = a+b*zc     #   Vitesse des couches


#
# epaisseur des couches
#
h = np.diff(z)


#
# dist horizontale
#

xx, yy = np.meshgrid(x,y)
r = np.sqrt( (xx-x0)**2 + (yy-y0)**2 )
ir = np.arctan2( (yy-y0), (xx-x0) )

#
# angles critiques aux interfaces
#

ic = np.arcsin(V[:-1]/V[1:])
i = np.zeros((N, N))
for nc in np.arange(N):
    i[nc,nc] = ic[nc]

for nc in np.arange(1, N):
    for na in np.arange(nc-1,-1,-1):
        i[na,nc] = np.arcsin( V[na]*np.sin(i[na+1,nc])/V[na+1] )

tt = np.zeros((x.size, y.size, z.size))

rl = tt.copy()



#  en surface, surf no 1

t1 = r/V[0]
rl1 = r.copy()


rp = [np.nan+np.ones((25,3)) for x in range((N+1)**3)]


for ix in np.arange(N+1):
    for iy in np.arange(N+1):
        rp[ind(ix,iy,0)][0,:] =  np.array([[x0, y0, z0]])
        rp[ind(ix,iy,0)][1,:] =  np.array([[x[ix],y[iy],z[0]]])


# ondes refractees subsequentes
#
#  t = x/v_n + sum_{i=1}^{n-1}\frac{2h_i\cos \theta_i}{V_i}
#
# \theta_i non critique sauf \theta_{n-1}

for n in np.arange(1, N+1):
    t = r/V[n] + np.sum(2*h[:n] * np.cos(i[:n,n-1]) / V[:n])
    imin = t<t1
    t1[imin] = t[imin]
    rl1[imin] = np.sum(2*h[:n] * np.cos(i[:n,n-1])) + V[n]*(t[imin]-np.sum(2*h[:n] * np.cos(i[:n,n-1])/V[:n] ))

    xc = h[:n] * np.tan(i[:n,n-1])
    dx = np.empty((N+1,N+1,len(xc)))
    dy = np.empty((N+1,N+1,len(xc)))
    for nn in np.arange(len(xc)):
        dx[:,:,nn] = xc[nn] * np.sin(ir)
        dy[:,:,nn] = xc[nn] * np.cos(ir)

    for ix in np.arange(N+1):
        for iy in np.arange(N+1):
            if imin[ix,iy]==1:
                rp[ind(ix,iy,0)][0,:] =  np.array([[x0, y0, z0]])
                for nn in np.arange(len(xc)):
                    rp[ind(ix,iy,0)][1+nn,:] = np.array([x0+np.sum(dx[ix,iy,:nn+1]),
                      y0+np.sum(dy[ix,iy,:nn+1]),z[1+nn]])
                    rp[ind(ix,iy,0)][2*n-nn,:] = np.array([x[ix]-np.sum(dx[ix,iy,:nn+1]),
                      y[iy]-np.sum(dy[ix,iy,:nn+1]),z[1+nn]])

                rp[ind(ix,iy,0)][2*(n+1)-1,:] = np.array([x[ix],y[iy],z[0]])

tt[:,:,0] = t1
rl[:,:,0] = rl1



# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ns = 0
# for ix in np.arange(1,N+1):
#     plt.plot(rp[ind(ix,ix,ns)][:,0], rp[ind(ix,ix,ns)][:,1], rp[ind(ix,ix,ns)][:,2])
# ax.invert_zaxis()
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
# ax.view_init(0, -45)
# plt.show()





# sous la 1re couche, surf no 2

d = np.sqrt( (xx-x0)**2 + (yy-y0)**2 + (z[1]-z0)**2 )
t1 = d/V[0]
rl1 = d.copy()

for ix in np.arange(N+1):
    for iy in np.arange(N+1):
        rp[ind(ix,iy,1)][0,:] = np.array([x0,y0,z0])
        rp[ind(ix,iy,1)][1,:] = np.array([x[ix],y[iy],z[1]])


#
# refraction @ l'interface consideree
#
ns=1

# distances horizontales parcourues par les rais obliques
xc = np.tan(i[:ns,ns-1]) * h[:ns]

#
l = np.sqrt( xc**2 + h[:ns]**2 )
dt = l / V[:ns]

x2 = np.abs(r-np.sum(xc))
t = np.sum(dt) + x2/V[ns]
imin = t<t1
t1[imin] = t[imin];
rl1[imin] = np.sum(l) + V[ns]*(t[imin] - np.sum(dt))


dx = xc * np.sin(ir)
dy = xc * np.cos(ir)
for ix in np.arange(N+1):
    for iy in np.arange(N+1):
        if ( imin[ix,iy]==1 ):
            rp[ind(ix,iy,ns)][0,:] = np.array([x0,y0,z0])
            rp[ind(ix,iy,ns)][1,:] = np.array([x0+dx[ix,iy],y0+dy[ix,iy],z[ns]])
            rp[ind(ix,iy,ns)][2,:] = np.array([x[ix],y[iy],z[ns]])


for n in np.arange(ns+1, N+1):

    xc = np.tan(i[:n,n-1]) * h[:n]

    #
    l = np.sqrt( xc**2 + h[:n]**2 )
    dt = l / V[:n]

    x2 = np.abs(r - np.sum(xc)-np.sum(xc[ns:n]))
    t = np.sum(dt) + np.sum(dt[ns:n]) + x2/V[n]
    imin = t<t1
    t1[imin] = t[imin]
    rl1[imin] = np.sum(l)+np.sum(l[ns:n])+V[n]*(t[imin]-np.sum(dt)-np.sum(dt[ns:n]))

    dx = np.empty((N+1,N+1,len(xc)))
    dy = np.empty((N+1,N+1,len(xc)))
    for nn in np.arange(len(xc)):
        dx[:,:,nn] = xc[nn] * np.sin(ir)
        dy[:,:,nn] = xc[nn] * np.cos(ir)

    for ix in np.arange(N+1):
        for iy in np.arange(N+1):
            if imin[ix,iy]==1:
                rp[ind(ix,iy,ns)] = np.nan+np.ones((25,3))
                rp[ind(ix,iy,ns)][0,:] = np.array([x0,y0,z0])
                for nn in np.arange(len(xc)):
                    rp[ind(ix,iy,ns)][1+nn,:] = np.array([x0+np.sum(dx[ix,iy,:nn+1]),
                      y0+np.sum(dy[ix,iy,:nn+1]),z[1+nn]])

                for nn in np.arange(ns-1,len(xc)):
                    rp[ind(ix,iy,ns)][2*n-nn,:] = np.array([x[ix]-np.sum(dx[ix,iy,ns:nn+1]),
                      y[iy]-np.sum(dy[ix,iy,ns:nn+1]),z[1+nn]])

#                rp[ind(ix,iy,ns)][2*(n+1)-1,:] = np.array([x[ix],y[iy],z[ns]])

tt[:,:,ns] = t1
rl[:,:,ns] = rl1





# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# for ix in np.arange(1, N+1):
#     plt.plot(rp[ind(ix,ix,ns)][:,0], rp[ind(ix,ix,ns)][:,1], rp[ind(ix,ix,ns)][:,2])
#     print(rp[ind(ix,ix,ns)])
# ax.invert_zaxis()
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
# ax.set_title(str(ns+1))
# ax.autoscale_view()
# ax.view_init(0, -45)
# plt.show()





for ns in np.arange(2, N+1):
    print(ns)
    # distance horizontale parcourue par les rais obliques, avant d'etre refracte
    xc = np.tan(i[:ns,ns-1])*h[:ns]

    for n1 in np.arange(r.shape[0]):
        for n2 in np.arange(r.shape[1]):

            t1[n1,n2], rl1[n1,n2], xi = solve_nlayers(h[:ns], V[:ns], r[n1,n2])

            dx = np.empty((N+1,N+1,len(xi)))
            dy = np.empty((N+1,N+1,len(xi)))
            for nn in np.arange(len(xi)):
                dx[:,:,nn] = xi[nn] * np.sin(ir)
                dy[:,:,nn] = xi[nn] * np.cos(ir)

            rp[ind(n1,n2,ns)][0,:] = np.array([x0,y0,z0])
            for nn in np.arange(len(xi)):
                rp[ind(n1,n2,ns)][1+nn,:] = np.array([x0+np.sum(dx[n1,n2,:nn+1]),
                  y0+np.sum(dy[n1,n2,:nn+1]),z[1+nn]])

#            rp[ind(n1,n2,ns)][ns+1,:] = np.array([x[n1],y[n2],z[ns]])

            if not np.isreal( t1[n1,n2] ):
                t1[n1,n2] = np.nan
                rl1[n1,n2] = 0



    #
    # temps de parcours le long des rais obliques
    l = np.sqrt( xc**2 + h[:ns]**2 )
    dt = l / V[:ns]



    # distance parcourue par l'onde refractee
    x2 = np.abs(r - np.sum(xc))
    t = np.sum(dt) + x2 / V[ns]
    imin = t<t1
    t1[imin] = t[imin]
    rl1[imin] = np.sum(l) + x2[imin]


    dx = np.empty((N+1,N+1,len(xc)))
    dy = np.empty((N+1,N+1,len(xc)))
    for nn in np.arange(len(xc)):
        dx[:,:,nn] = xc[nn] * np.sin(ir)
        dy[:,:,nn] = xc[nn] * np.cos(ir)

    for ix in np.arange(N+1):
        for iy in np.arange(N+1):
            if ( imin[ix,iy]==1 ):
                rp[ind(ix,iy,ns)][0,:] = np.array([x0,y0,z0])

                for nn in np.arange(len(xc)):
                    rp[ind(ix,iy,ns)][1+nn,:] = np.array([x0+np.sum(dx[ix,iy,:nn+1]),
                      y0+np.sum(dy[ix,iy,:nn+1]),z[1+nn]])

                rp[ind(ix,iy,ns)][ns+1,:] = np.array([x[ix],y[iy],z[ns]])



    for n in np.arange(ns+1, N+1):

        xc = np.tan(i[:n,n-1]) * h[:n]

        #

        l = np.sqrt( xc**2 + h[:n]**2 )
        dt = l / V[:n]

        x2 = np.abs(r - np.sum(xc) - np.sum(xc[ns:n]))
        t = np.sum(dt) + np.sum(dt[ns:n]) + x2 / V[n]
        imin = t<t1
        t1[imin] = t[imin]
        rl1[imin] = np.sum(l) + np.sum(l[ns:n]) + x2[imin]

        dx = np.empty((N+1,N+1,len(xc)))
        dy = np.empty((N+1,N+1,len(xc)))
        for nn in np.arange(len(xc)):
            dx[:,:,nn] = xc[nn] * np.sin(ir)
            dy[:,:,nn] = xc[nn] * np.cos(ir)

        for ix in np.arange(N+1):
            for iy in np.arange(N+1):
                if ( imin[ix,iy]==1 ):
                    rp[ind(ix,iy,ns)] = np.nan+np.ones((25,3))
                    rp[ind(ix,iy,ns)][0,:] = np.array([x0,y0,z0])
                    for nn in np.arange(len(xc)):
                        rp[ind(ix,iy,ns)][1+nn,:] = np.array([x0+np.sum(dx[ix,iy,:nn+1]),
                          y0+np.sum(dy[ix,iy,:nn+1]),z[1+nn]])

                    for nn in np.arange(ns-1,len(xc)):
                        rp[ind(ix,iy,ns)][2*n-nn,:] = np.array([x[ix]-np.sum(dx[ix,iy,ns:nn+1]),
                          y[iy]-np.sum(dy[ix,iy,ns:nn+1]),z[1+nn]])

#                    print(clean_nan(rp[ind(ix,iy,ns)]))

    tt[:,:,ns] = t1
    rl[:,:,ns] = rl1

# %%

for ns in range(N+1):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for ix in np.arange(1,N+1):
        plt.plot(rp[ind(ix,ix,ns)][:,0], rp[ind(ix,ix,ns)][:,1], rp[ind(ix,ix,ns)][:,2])
    ax.invert_zaxis()
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(str(ns+1))
    ax.view_init(0, -45)
    plt.show()

# %%
for ns in range(N+1):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for ix in np.arange(1,N+1):
        plt.plot(rp[ind(ix,0,ns)][:,0], rp[ind(ix,0,ns)][:,1], rp[ind(ix,0,ns)][:,2])
    ax.invert_zaxis()
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(str(ns+1))
    ax.set_ylim(0, 10)
    ax.view_init(0, -90)
    plt.show()

# # %%
# for ns in range(N+1):
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     for ix in np.arange(1,N+1):
#         plt.plot(rp[ind(0,ix,ns)][:,0], rp[ind(0,ix,ns)][:,1], rp[ind(0,ix,ns)][:,2])
#     ax.invert_zaxis()
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_zlabel('Z')
#     ax.set_title(str(ns+1))
#     ax.view_init(0, 0)
#     plt.show()

# %%
import vtk

xCoords = vtk.vtkFloatArray()
for i in x:
    xCoords.InsertNextValue(i)

yCoords = vtk.vtkFloatArray()
for i in y:
    yCoords.InsertNextValue(i)

zCoords = vtk.vtkFloatArray()
for i in z:
    zCoords.InsertNextValue(i)

rgrid = vtk.vtkRectilinearGrid()
rgrid.SetDimensions(len(x), len(y), len(z))
rgrid.SetXCoordinates(xCoords)
rgrid.SetYCoordinates(yCoords)
rgrid.SetZCoordinates(zCoords)

ttData = vtk.vtkDoubleArray()
ttData.SetName('Travel Time')
ttData.SetNumberOfComponents(1)
ttData.SetNumberOfTuples(tt.size)
for ix in range(len(x)):
    for iy in range(len(y)):
        for iz in range(len(z)):
            ii = rgrid.FindPoint(x[ix], y[iy], z[iz])
            ttData.SetTuple1(ii, tt[ix,iy,iz])
rgrid.GetPointData().SetScalars(ttData)

writer = vtk.vtkXMLRectilinearGridWriter()
writer.SetFileName('sol_analytique_couches_tt.vtr')
writer.SetInputData(rgrid)
writer.SetDataModeToBinary()
writer.Update()


yCoords = vtk.vtkFloatArray()
yCoords.InsertNextValue(0)

rgrid = vtk.vtkRectilinearGrid()
rgrid.SetDimensions(len(x), 1, len(z))
rgrid.SetXCoordinates(xCoords)
rgrid.SetYCoordinates(yCoords)
rgrid.SetZCoordinates(zCoords)

ttData = vtk.vtkDoubleArray()
ttData.SetName('Travel Time')
ttData.SetNumberOfComponents(1)
ttData.SetNumberOfTuples(len(x)*len(z))
for ix in range(len(x)):
    for iz in range(len(z)):
        ii = rgrid.FindPoint(x[ix], 0, z[iz])
        ttData.SetTuple1(ii, tt[ix,0,iz])
rgrid.GetPointData().SetScalars(ttData)

writer = vtk.vtkXMLRectilinearGridWriter()
writer.SetFileName('sol_analytique_couches2d_tt.vtr')
writer.SetInputData(rgrid)
writer.SetDataModeToBinary()
writer.Update()

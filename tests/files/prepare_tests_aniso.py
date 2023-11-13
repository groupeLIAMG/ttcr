#! /usr/bin/env python3
#  -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import vtk

import traveltime as wa_tt

# %%

x = np.arange(0.0, 6000.01, 60.0)
z = np.arange(0.0, 3000.01, 60.0)


a = 2271.0
b = 0.88
chi = 0.3

with open('rcv2daniso.dat', 'w') as f:
    print('{0:d}'.format(x.size*z.size), file=f)
    for ix in x:
        for iz in z:
            print('{0:g} {1:g}'.format(ix, iz), file=f)


# %%

x[0] = 0.001

xx, zz = np.meshgrid(x, z)

p = 2*xx / np.sqrt((xx*xx + (1+2*chi) * zz*zz) *
                  ((2*a + b*zz)**2 * (1+2*chi) + b*b * xx*xx))

# p[0, 0] = 0.0

t = 1/b * (np.arctanh(p * b * xx - np.sqrt(1 - (1+2*chi) * p*p * a*a)) +
           np.arctanh(np.sqrt(1 - (1+2*chi) * p*p * a*a)))

# t[0, 0] = 0
# t[:, 0] = z / vz


# %%

h = plt.pcolor(xx, zz, t)
plt.contour(xx, zz, t, colors='k', linewidths=0.5)
plt.gca().invert_yaxis()
plt.gca().axis('equal')
plt.colorbar(h)
plt.tight_layout()
plt.show()

# %%

x = np.arange(30.0, 6000.0, 60.0)
z = np.arange(30.0, 3000.0, 60.0)

xx, zz = np.meshgrid(x, z)

vz = a + b*zz
vx = np.sqrt( chi * 2 * vz**2 + vz**2)


sx = 1/vx
sz = 1/vz

xi = sz / sx

plt.imshow(xi)
plt.colorbar()
plt.show()

plt.imshow(sx)
plt.colorbar()
plt.show()



xCoords = vtk.vtkDoubleArray()
yCoords = vtk.vtkDoubleArray()
zCoords = vtk.vtkDoubleArray()

for i in np.arange(0.0, 6000.1, 60.0):
    xCoords.InsertNextValue(i)
for i in np.arange(0.0, 3000.1, 60.0):
    zCoords.InsertNextValue(i)
yCoords.InsertNextValue(0.0)

rgrid = vtk.vtkRectilinearGrid()
rgrid.SetDimensions(xCoords.GetNumberOfValues(),
                    yCoords.GetNumberOfValues(),
                    zCoords.GetNumberOfValues())
rgrid.SetXCoordinates(xCoords)
rgrid.SetYCoordinates(yCoords)
rgrid.SetZCoordinates(zCoords)

slown = vtk.vtkDoubleArray()
for i in sx.flatten():
    slown.InsertNextValue(i)
slown.SetName("Slowness")
rgrid.GetCellData().AddArray(slown)

xi_vtk = vtk.vtkDoubleArray()
for i in xi.flatten():
    xi_vtk.InsertNextValue(i)
xi_vtk.SetName("xi")
rgrid.GetCellData().AddArray(xi_vtk)

writer = vtk.vtkXMLRectilinearGridWriter()
writer.SetInputData(rgrid)
writer.SetFileName('elliptical_fine2d.vtr')
writer.SetDataModeToBinary()
writer.Write()

# %%
rgrid = vtk.vtkRectilinearGrid()
rgrid.SetDimensions(xCoords.GetNumberOfValues(),
                    yCoords.GetNumberOfValues(),
                    zCoords.GetNumberOfValues())
rgrid.SetXCoordinates(xCoords)
rgrid.SetYCoordinates(yCoords)
rgrid.SetZCoordinates(zCoords)

tt = vtk.vtkDoubleArray()
for i in t.flatten():
    tt.InsertNextValue(i)
tt.SetName("Travel Time")
rgrid.GetPointData().AddArray(tt)

writer = vtk.vtkXMLRectilinearGridWriter()
writer.SetInputData(rgrid)
writer.SetFileName('sol_analytique_elliptical_2d_tt.vtr')
writer.SetDataModeToBinary()
writer.Write()

# %%

x = np.arange(0.0, 6000.01, 60.0)
z = np.arange(0.0, 3000.01, 60.0)

xx, zz = np.meshgrid(x, z)

xx = xx.flatten()
zz = zz.flatten()
t = t.flatten()

npx = x.size
ncx = x.size - 1
ncz = z.size - 1

tri = []
slo = []
xi = []
for nz in range(ncz):
    for nx in range(ncx):
        tri.append([nz*npx + nx, nz*npx + nx+1, (nz+1)*npx + nx])
        
        zc = (zz[tri[-1][0]] + zz[tri[-1][1]] + zz[tri[-1][2]]) / 3.
        vz = a + b*zc
        vx = np.sqrt( chi * 2 * vz**2 + vz**2)

        sx = 1/vx
        sz = 1/vz

        slo.append(sx)
        xi.append( sz / sx )
        
        tri.append([nz*npx + nx+1, (nz+1)*npx + nx, (nz+1)*npx + nx + 1])

        zc = (zz[tri[-1][0]] + zz[tri[-1][1]] + zz[tri[-1][2]]) / 3.
        vz = a + b*zc
        vx = np.sqrt( chi * 2 * vz**2 + vz**2)

        sx = 1/vx
        sz = 1/vz

        slo.append(sx)
        xi.append( sz / sx )
        
# %%

ugrid = vtk.vtkUnstructuredGrid()
pts = vtk.vtkPoints()

for n in range(len(xx)):
    pts.InsertNextPoint(xx[n], 0.0, zz[n])

ugrid.SetPoints(pts)

triangle = vtk.vtkTriangle()
d_s = vtk.vtkDoubleArray()
d_x = vtk.vtkDoubleArray()
d_s.SetName("Slowness")
d_x.SetName("xi")
for n in range(len(tri)):
    triangle.GetPointIds().SetId(0, tri[n][0])
    triangle.GetPointIds().SetId(1, tri[n][1])
    triangle.GetPointIds().SetId(2, tri[n][2])
    ugrid.InsertNextCell( triangle.GetCellType(), triangle.GetPointIds() )
    d_s.InsertNextValue(slo[n])
    d_x.InsertNextValue(xi[n])

ugrid.GetCellData().AddArray(d_s)
ugrid.GetCellData().AddArray(d_x)

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetInputData(ugrid)
writer.SetFileName('elliptical_fine2d.vtu')
writer.SetDataModeToBinary()
writer.Write()

# %%

ugrid = vtk.vtkUnstructuredGrid()
pts = vtk.vtkPoints()

for n in range(len(xx)):
    pts.InsertNextPoint(xx[n], 0.0, zz[n])

ugrid.SetPoints(pts)

triangle = vtk.vtkTriangle()
for n in range(len(tri)):
    triangle.GetPointIds().SetId(0, tri[n][0])
    triangle.GetPointIds().SetId(1, tri[n][1])
    triangle.GetPointIds().SetId(2, tri[n][2])
    ugrid.InsertNextCell( triangle.GetCellType(), triangle.GetPointIds() )

d_t = vtk.vtkDoubleArray()
d_t.SetName("Travel Time")
for n in range(len(t)):
    d_t.InsertNextValue(t[n])
    
ugrid.GetPointData().AddArray(d_t)

# writer = vtk.vtkXMLUnstructuredGridWriter()
# writer.SetInputData(ugrid)
# writer.SetFileName('tests/files/sol_analytique_elliptical_2d_tt.vtu')
# writer.SetDataModeToBinary()
# writer.Write()

# %%

VP0 = 1875.
r2 = 0.1
r4 = 0.15602005

rrr = np.array([1., np.NAN, r2, np.NAN, r4])   # generic anisotropy parameters (here, for Dog Creek Shale)

scs = wa_tt.precompute(v00=VP0, rrr=rrr)                                 # precompute energy velocity and construct an interpolation polynomial

x = np.arange(0.0, 6000.01, 60.0)
z = np.arange(0.0, 3000.01, 60.0)

ncx = x.size - 1
ncz = z.size - 1

xx, zz = np.meshgrid(x, z)

xx = xx.flatten()
zz = zz.flatten()

t = wa_tt.traveltime(xxx=xx, zzz=zz, scs=scs)

t[0] = 0.0


# %%
xCoords = vtk.vtkDoubleArray()
yCoords = vtk.vtkDoubleArray()
zCoords = vtk.vtkDoubleArray()

for i in np.arange(0.0, 6000.1, 60.0):
    xCoords.InsertNextValue(i)
for i in np.arange(0.0, 3000.1, 60.0):
    zCoords.InsertNextValue(i)
yCoords.InsertNextValue(0.0)

rgrid = vtk.vtkRectilinearGrid()
rgrid.SetDimensions(xCoords.GetNumberOfValues(),
                    yCoords.GetNumberOfValues(),
                    zCoords.GetNumberOfValues())
rgrid.SetXCoordinates(xCoords)
rgrid.SetYCoordinates(yCoords)
rgrid.SetZCoordinates(zCoords)

v = vtk.vtkDoubleArray()
r2_v = vtk.vtkDoubleArray()
r4_v = vtk.vtkDoubleArray()

for i in range(ncx*ncz):
    v.InsertNextValue(VP0)
    r2_v.InsertNextValue(r2)
    r4_v.InsertNextValue(r4)
v.SetName("Velocity")
r2_v.SetName("r2")
r4_v.SetName("r4")
rgrid.GetCellData().AddArray(v)
rgrid.GetCellData().AddArray(r2_v)
rgrid.GetCellData().AddArray(r4_v)

writer = vtk.vtkXMLRectilinearGridWriter()
writer.SetInputData(rgrid)
writer.SetFileName('weakly_an_fine2d.vtr')
writer.SetDataModeToBinary()
writer.Write()


rgrid = vtk.vtkRectilinearGrid()
rgrid.SetDimensions(xCoords.GetNumberOfValues(),
                    yCoords.GetNumberOfValues(),
                    zCoords.GetNumberOfValues())
rgrid.SetXCoordinates(xCoords)
rgrid.SetYCoordinates(yCoords)
rgrid.SetZCoordinates(zCoords)

tt = vtk.vtkDoubleArray()
for i in t.flatten():
    tt.InsertNextValue(i)
tt.SetName("Travel Time")
rgrid.GetPointData().AddArray(tt)

writer = vtk.vtkXMLRectilinearGridWriter()
writer.SetInputData(rgrid)
writer.SetFileName('sol_analytique_weakly_an_2d_tt.vtr')
writer.SetDataModeToBinary()
writer.Write()

# %%

npx = x.size

tri = []
for nz in range(ncz):
    for nx in range(ncx):
        tri.append([nz*npx + nx, nz*npx + nx+1, (nz+1)*npx + nx])
        tri.append([nz*npx + nx+1, (nz+1)*npx + nx, (nz+1)*npx + nx + 1])
        
ugrid = vtk.vtkUnstructuredGrid()
pts = vtk.vtkPoints()

for n in range(len(xx)):
    pts.InsertNextPoint(xx[n], 0.0, zz[n])

ugrid.SetPoints(pts)

triangle = vtk.vtkTriangle()
v = vtk.vtkDoubleArray()
r2_v = vtk.vtkDoubleArray()
r4_v = vtk.vtkDoubleArray()
v.SetName("Velocity")
r2_v.SetName("r2")
r4_v.SetName("r4")

for n in range(len(tri)):
    triangle.GetPointIds().SetId(0, tri[n][0])
    triangle.GetPointIds().SetId(1, tri[n][1])
    triangle.GetPointIds().SetId(2, tri[n][2])
    ugrid.InsertNextCell( triangle.GetCellType(), triangle.GetPointIds() )
    v.InsertNextValue(VP0)
    r2_v.InsertNextValue(r2)
    r4_v.InsertNextValue(r4)

ugrid.GetCellData().AddArray(v)
ugrid.GetCellData().AddArray(r2_v)
ugrid.GetCellData().AddArray(r4_v)

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetInputData(ugrid)
writer.SetFileName('weakly_an_fine2d.vtu')
writer.SetDataModeToBinary()
writer.Write()
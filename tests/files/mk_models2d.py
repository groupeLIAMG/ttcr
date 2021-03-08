#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 13:54:30 2021

@author: giroux
"""
import numpy as np

import vtk
from math import floor

import gmsh

N = 20
dx = 1.0


# Rcv

with open('rcv2d.dat', 'w') as f:
    xx = np.arange(0, N+1)
    print('{0:d}'.format(xx.size*xx.size), file=f)
    for x in xx:
        for z in xx:
            print('{0:d} {1:d}'.format(x, z), file=f)

with open('rcv2d_in.dat', 'w') as f:
    xx = np.arange(1, N)
    print('{0:d}'.format(xx.size*xx.size), file=f)
    for x in xx:
        for z in xx:
            print('{0:d} {1:d}'.format(x, z), file=f)

a = 1.0
V20 = 3.0
b = (V20-a)/20.0

#
# Coarse grid
#

frac = 2
N *= frac
dx /= frac


x = np.arange(0.0, 20.0+0.01, dx)
z = np.arange(0.0, 20.0+0.01, dx)


# %% gradient

ugrid = vtk.vtkUnstructuredGrid()
pts = vtk.vtkPoints()
slowness = vtk.vtkDoubleArray()
slowness.SetName('Slowness')

for x in np.arange(0.0, 20.0+0.01, dx):
    for z in np.arange(0.0, 20.0+0.01, dx):
        pts.InsertNextPoint(x, 0.0, z)
        s = 1.0 / (a+b*z)
        slowness.InsertNextValue(s)

ugrid.SetPoints(pts)
ugrid.GetPointData().SetScalars(slowness)

tri = vtk.vtkTriangle()

Np = N + 1
for i in np.arange(N):
    for j in np.arange(N):
        tri.GetPointIds().SetId(0, i*Np+j)
        tri.GetPointIds().SetId(1, (i+1)*Np+j)
        tri.GetPointIds().SetId(2, i*Np+j+1)
        ugrid.InsertNextCell( tri.GetCellType(), tri.GetPointIds() )

        tri.GetPointIds().SetId(0, i*Np+j+1)
        tri.GetPointIds().SetId(1, (i+1)*Np+j)
        tri.GetPointIds().SetId(2, (i+1)*Np+j+1)
        ugrid.InsertNextCell( tri.GetCellType(), tri.GetPointIds() )

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName('gradient_coarse2d.vtu')
writer.SetInputData( ugrid )
writer.SetDataModeToBinary();
writer.Update()

xCoords = vtk.vtkDoubleArray()
for n in range(0, N+1):
    xCoords.InsertNextValue( n*dx )
yCoords = vtk.vtkDoubleArray()
yCoords.InsertNextValue(0.0)
zCoords = xCoords


rgrid = vtk.vtkRectilinearGrid()
rgrid.SetDimensions( N+1, 1, N+1 )
rgrid.SetXCoordinates(xCoords)
rgrid.SetYCoordinates(yCoords)
rgrid.SetZCoordinates(zCoords)

slowness = vtk.vtkDoubleArray()
slowness.SetName("Slowness")
slowness.SetNumberOfComponents(1)
slowness.SetNumberOfTuples( rgrid.GetNumberOfPoints() )

for n in range(0, rgrid.GetNumberOfPoints()):
    x = rgrid.GetPoint(n)
    s = 1.0 / (a+b*x[2])
    slowness.SetTuple1(n, s)

rgrid.GetPointData().SetScalars( slowness );

writer = vtk.vtkXMLRectilinearGridWriter()
writer.SetFileName( "gradient_coarse2d.vtr" )
writer.SetInputData( rgrid )
writer.SetDataModeToBinary()
writer.Update()



# %% layers

ugrid = vtk.vtkUnstructuredGrid()
pts = vtk.vtkPoints()

for x in np.arange(0.0, 20.0+0.01, dx):
    for z in np.arange(0.0, 20.0+0.01, dx):
        pts.InsertNextPoint(x, 0.0, z)

ugrid.SetPoints(pts)

tri = vtk.vtkTriangle()
slowness = vtk.vtkDoubleArray()
slowness.SetName('Slowness')

Np = N + 1
for i in np.arange(N):
    for j in np.arange(N):
        tri.GetPointIds().SetId(0, i*Np+j)
        tri.GetPointIds().SetId(1, (i+1)*Np+j)
        tri.GetPointIds().SetId(2, i*Np+j+1)
        ugrid.InsertNextCell( tri.GetCellType(), tri.GetPointIds() )
        z = floor(j*dx + 0.001) + 0.5
        s = 1.0 / (a+b*z)
        slowness.InsertNextValue(s)

        tri.GetPointIds().SetId(0, i*Np+j+1)
        tri.GetPointIds().SetId(1, (i+1)*Np+j)
        tri.GetPointIds().SetId(2, (i+1)*Np+j+1)
        ugrid.InsertNextCell( tri.GetCellType(), tri.GetPointIds() )
        slowness.InsertNextValue(s)

ugrid.GetCellData().SetScalars(slowness)

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName('layers_coarse2d.vtu')
writer.SetInputData( ugrid )
writer.SetDataModeToBinary();
writer.Update()

xCoords = vtk.vtkDoubleArray()
for n in range(0, N+1):
    xCoords.InsertNextValue( n*dx )
yCoords = vtk.vtkDoubleArray()
yCoords.InsertNextValue(0.0)
zCoords = xCoords


rgrid = vtk.vtkRectilinearGrid()
rgrid.SetDimensions( N+1, 1, N+1 )
rgrid.SetXCoordinates(xCoords)
rgrid.SetYCoordinates(yCoords)
rgrid.SetZCoordinates(zCoords)

slowness = vtk.vtkDoubleArray()
slowness.SetName("Slowness")
slowness.SetNumberOfComponents(1)
slowness.SetNumberOfTuples( rgrid.GetNumberOfCells() )

for n in range(0, rgrid.GetNumberOfCells()):
    bo = rgrid.GetCell(n).GetBounds()
    z = floor(bo[4]) + 0.5
    s = 1.0 / (a+b*z)
    slowness.SetTuple1(n, s)

rgrid.GetCellData().SetScalars( slowness );

writer = vtk.vtkXMLRectilinearGridWriter()
writer.SetFileName( "layers_coarse2d.vtr" )
writer.SetInputData( rgrid )
writer.SetDataModeToBinary()
writer.Update()

# %%
# Finer grid
#


frac = 2.5
N = int(frac*N + 0.000000001)
dx /= frac


x = np.arange(0.0, 20.0+0.01, dx)
z = np.arange(0.0, 20.0+0.01, dx)


# %% gradient

ugrid = vtk.vtkUnstructuredGrid()
pts = vtk.vtkPoints()
slowness = vtk.vtkDoubleArray()
slowness.SetName('Slowness')

for x in np.arange(0.0, 20.0+0.01, dx):
    for z in np.arange(0.0, 20.0+0.01, dx):
        pts.InsertNextPoint(x, 0.0, z)
        s = 1.0 / (a+b*z)
        slowness.InsertNextValue(s)

ugrid.SetPoints(pts)
ugrid.GetPointData().SetScalars(slowness)

tri = vtk.vtkTriangle()

Np = N + 1
for i in np.arange(N):
    for j in np.arange(N):
        tri.GetPointIds().SetId(0, i*Np+j)
        tri.GetPointIds().SetId(1, (i+1)*Np+j)
        tri.GetPointIds().SetId(2, i*Np+j+1)
        ugrid.InsertNextCell( tri.GetCellType(), tri.GetPointIds() )

        tri.GetPointIds().SetId(0, i*Np+j+1)
        tri.GetPointIds().SetId(1, (i+1)*Np+j)
        tri.GetPointIds().SetId(2, (i+1)*Np+j+1)
        ugrid.InsertNextCell( tri.GetCellType(), tri.GetPointIds() )

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName('gradient_fine2d.vtu')
writer.SetInputData( ugrid )
writer.SetDataModeToBinary();
writer.Update()


xCoords = vtk.vtkDoubleArray()
for n in range(0, N+1):
    xCoords.InsertNextValue( n*dx )
yCoords = vtk.vtkDoubleArray()
yCoords.InsertNextValue(0.0)
zCoords = xCoords


rgrid = vtk.vtkRectilinearGrid()
rgrid.SetDimensions( N+1, 1, N+1 )
rgrid.SetXCoordinates(xCoords)
rgrid.SetYCoordinates(yCoords)
rgrid.SetZCoordinates(zCoords)

slowness = vtk.vtkDoubleArray()
slowness.SetName("Slowness")
slowness.SetNumberOfComponents(1)
slowness.SetNumberOfTuples( rgrid.GetNumberOfPoints() )

for n in range(0, rgrid.GetNumberOfPoints()):
    x = rgrid.GetPoint(n)
    s = 1.0 / (a+b*x[2])
    slowness.SetTuple1(n, s)

rgrid.GetPointData().SetScalars( slowness );

writer = vtk.vtkXMLRectilinearGridWriter()
writer.SetFileName( "gradient_fine2d.vtr" )
writer.SetInputData( rgrid )
writer.SetDataModeToBinary()
writer.Update()


# %% layers

ugrid = vtk.vtkUnstructuredGrid()
pts = vtk.vtkPoints()

for x in np.arange(0.0, 20.0+0.01, dx):
    for z in np.arange(0.0, 20.0+0.01, dx):
        pts.InsertNextPoint(x, 0.0, z)

ugrid.SetPoints(pts)

tri = vtk.vtkTriangle()
slowness = vtk.vtkDoubleArray()
slowness.SetName('Slowness')

Np = N + 1
for i in np.arange(N):
    for j in np.arange(N):
        tri.GetPointIds().SetId(0, i*Np+j)
        tri.GetPointIds().SetId(1, (i+1)*Np+j)
        tri.GetPointIds().SetId(2, i*Np+j+1)
        ugrid.InsertNextCell( tri.GetCellType(), tri.GetPointIds() )
        z = floor(j*dx + 0.001) + 0.5
        s = 1.0 / (a+b*z)
        slowness.InsertNextValue(s)

        tri.GetPointIds().SetId(0, i*Np+j+1)
        tri.GetPointIds().SetId(1, (i+1)*Np+j)
        tri.GetPointIds().SetId(2, (i+1)*Np+j+1)
        ugrid.InsertNextCell( tri.GetCellType(), tri.GetPointIds() )
        slowness.InsertNextValue(s)

ugrid.GetCellData().SetScalars(slowness)

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName('layers_fine2d.vtu')
writer.SetInputData( ugrid )
writer.SetDataModeToBinary();
writer.Update()


xCoords = vtk.vtkDoubleArray()
for n in range(0, N+1):
    xCoords.InsertNextValue( n*dx )
yCoords = vtk.vtkDoubleArray()
yCoords.InsertNextValue(0.0)
zCoords = xCoords


rgrid = vtk.vtkRectilinearGrid()
rgrid.SetDimensions( N+1, 1, N+1 )
rgrid.SetXCoordinates(xCoords)
rgrid.SetYCoordinates(yCoords)
rgrid.SetZCoordinates(zCoords)

slowness = vtk.vtkDoubleArray()
slowness.SetName("Slowness")
slowness.SetNumberOfComponents(1)
slowness.SetNumberOfTuples( rgrid.GetNumberOfCells() )

for n in range(0, rgrid.GetNumberOfCells()):
    bo = rgrid.GetCell(n).GetBounds()
    z = floor(bo[4]) + 0.5
    s = 1.0 / (a+b*z)
    slowness.SetTuple1(n, s)

rgrid.GetCellData().SetScalars( slowness );

writer = vtk.vtkXMLRectilinearGridWriter()
writer.SetFileName( "layers_fine2d.vtr" )
writer.SetInputData( rgrid )
writer.SetDataModeToBinary()
writer.Update()


# %%
gmsh.initialize()

# %%
meshSize = 0.43
gmsh.clear()
tag = 1
# for x in range(21):
#     for z in range(21):
#         gmsh.model.geo.addPoint(x, 0, z, meshSize=meshSize, tag=tag)
#         tag += 1

# pts = []
# for k in range(21):
#     pts.append(k)
# for i in range(1, 21):
#     pts.append(i*21+20)
# for k in range(19, -1, -1):
#     pts.append(20*21 + k)
# for i in range(19, -1, -1):
#     pts.append(i*21)

gmsh.model.geo.addPoint(0, 0, 0, meshSize=meshSize, tag=1)
gmsh.model.geo.addPoint(0, 0, 20, meshSize=meshSize, tag=2)
gmsh.model.geo.addPoint(20, 0, 20, meshSize=meshSize, tag=3)
gmsh.model.geo.addPoint(20, 0, 0, meshSize=meshSize, tag=4)

# add a few points "randomly" for have a more irregular mesh
gmsh.model.geo.addPoint(13.43, 0, 3.3312, meshSize=meshSize, tag=5)
gmsh.model.geo.addPoint(3.43, 0, 13.3312, meshSize=meshSize, tag=6)
gmsh.model.geo.addPoint(8.43, 0, 8.3312, meshSize=meshSize, tag=7)
gmsh.model.geo.addPoint(17.43, 0, 7.3312, meshSize=meshSize, tag=8)
gmsh.model.geo.addPoint(1.43, 0, 17.12, meshSize=meshSize, tag=9)
gmsh.model.geo.addPoint(1.43, 0, 1.3312, meshSize=meshSize, tag=10)
gmsh.model.geo.addPoint(9.43, 0, 15.3312, meshSize=meshSize, tag=11)


pts = (0, 1, 2, 3, 0)

lines = []
for n in range(len(pts)-1):
    lines.append(gmsh.model.geo.addLine(pts[n]+1, pts[n+1]+1))

loop = gmsh.model.geo.addCurveLoop(lines)

srf = gmsh.model.geo.addPlaneSurface((loop,))

gmsh.model.geo.synchronize()

gmsh.model.mesh.embed(0, [5, 6, 7, 8, 9, 10, 11], 2, srf)

gmsh.model.geo.synchronize()

gmsh.model.geo.addPhysicalGroup(2, (srf,))

gmsh.model.mesh.generate(2)
gmsh.write('/tmp/mesh.vtk')

# %%

reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName('/tmp/mesh.vtk')
reader.Update()

slowness = vtk.vtkDoubleArray()
slowness.SetName('Slowness')
ugrid = reader.GetOutput()
for n in range(ugrid.GetNumberOfPoints()):
    x, y, z = ugrid.GetPoint(n)
    s = 1.0 / (a+b*z)
    slowness.InsertNextValue(s)

ugrid.GetPointData().SetScalars(slowness)

ugrid2 = vtk.vtkUnstructuredGrid()
ugrid2.SetPoints(ugrid.GetPoints())
ugrid2.GetPointData().SetScalars(slowness)
for n in range(ugrid.GetNumberOfCells()):
    cell = ugrid.GetCell(n)
    if type(cell) is vtk.vtkTriangle:
        ugrid2.InsertNextCell(cell.GetCellType(), cell.GetPointIds())

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName('gradient_gmsh2d.vtu')
writer.SetInputData( ugrid2 )
writer.SetDataModeToBinary();
writer.Update()


# %%
# nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes()
# elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(2)


# elemTags = elemTags[0]
# elemNodeTags = elemNodeTags[0]

# nodeCoords = np.reshape(nodeCoords, (-1, 3))
# elemNodeTags = np.reshape(elemNodeTags, (-1, 3))

# ugrid = vtk.vtkUnstructuredGrid()
# pts = vtk.vtkPoints()
# slowness = vtk.vtkDoubleArray()
# slowness.SetName('Slowness')

# for n in range(nodeCoords.shape[0]):
#     x = nodeCoords[n, 0]
#     z = nodeCoords[n, 2]
#     pts.InsertNextPoint(x, 0.0, z)
#     s = 1.0 / (a+b*z)
#     slowness.InsertNextValue(s)

# ugrid.SetPoints(pts)
# ugrid.GetPointData().SetScalars(slowness)

# tri = vtk.vtkTriangle()

# for n in np.arange(elemNodeTags.shape[0]):
#     tri.GetPointIds().SetId(0, int(elemNodeTags[n, 0]-1))
#     tri.GetPointIds().SetId(1, int(elemNodeTags[n, 1]-1))
#     tri.GetPointIds().SetId(2, int(elemNodeTags[n, 2]-1))
#     ugrid.InsertNextCell( tri.GetCellType(), tri.GetPointIds() )

# writer = vtk.vtkXMLUnstructuredGridWriter()
# writer.SetFileName('gradient_gmsh2d.vtu')
# writer.SetInputData( ugrid )
# writer.SetDataModeToBinary();
# writer.Update()

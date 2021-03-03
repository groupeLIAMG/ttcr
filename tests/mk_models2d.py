#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 13:54:30 2021

@author: giroux
"""
import numpy as np

from vtk import *
from math import floor

N = 20
dx = 1.0


# Rcv

with open('./files/rcv2d.dat', 'w') as f:
    xx = np.arange(0, N+1)
    print('{0:d}'.format(xx.size*xx.size), file=f)
    for x in xx:
        for z in xx:
            print('{0:d} {1:d}'.format(x, z), file=f)

with open('./files/rcv2d_in.dat', 'w') as f:
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

ugrid = vtkUnstructuredGrid()
pts = vtkPoints()
slowness = vtkDoubleArray()
slowness.SetName('Slowness')

for x in np.arange(0.0, 20.0+0.01, dx):
    for z in np.arange(0.0, 20.0+0.01, dx):
        pts.InsertNextPoint(x, 0.0, z)
        s = 1.0 / (a+b*z)
        slowness.InsertNextValue(s)

ugrid.SetPoints(pts)
ugrid.GetPointData().SetScalars(slowness)

tri = vtkTriangle()

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

writer = vtkXMLUnstructuredGridWriter()
writer.SetFileName('./files/gradient_coarse2d.vtu')
writer.SetInputData( ugrid )
writer.SetDataModeToBinary();
writer.Update()

xCoords = vtkDoubleArray()
for n in range(0, N+1):
    xCoords.InsertNextValue( n*dx )
yCoords = vtkDoubleArray()
yCoords.InsertNextValue(0.0)
zCoords = xCoords


rgrid = vtkRectilinearGrid()
rgrid.SetDimensions( N+1, 1, N+1 )
rgrid.SetXCoordinates(xCoords)
rgrid.SetYCoordinates(yCoords)
rgrid.SetZCoordinates(zCoords)

slowness = vtkDoubleArray()
slowness.SetName("Slowness")
slowness.SetNumberOfComponents(1)
slowness.SetNumberOfTuples( rgrid.GetNumberOfPoints() )

for n in range(0, rgrid.GetNumberOfPoints()):
    x = rgrid.GetPoint(n)
    s = 1.0 / (a+b*x[2])
    slowness.SetTuple1(n, s)

rgrid.GetPointData().SetScalars( slowness );

writer = vtkXMLRectilinearGridWriter()
writer.SetFileName( "./files/gradient_coarse2d.vtr" )
writer.SetInputData( rgrid )
writer.SetDataModeToBinary()
writer.Update()



# %% layers

ugrid = vtkUnstructuredGrid()
pts = vtkPoints()

for x in np.arange(0.0, 20.0+0.01, dx):
    for z in np.arange(0.0, 20.0+0.01, dx):
        pts.InsertNextPoint(x, 0.0, z)

ugrid.SetPoints(pts)

tri = vtkTriangle()
slowness = vtkDoubleArray()
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

writer = vtkXMLUnstructuredGridWriter()
writer.SetFileName('./files/layers_coarse2d.vtu')
writer.SetInputData( ugrid )
writer.SetDataModeToBinary();
writer.Update()

xCoords = vtkDoubleArray()
for n in range(0, N+1):
    xCoords.InsertNextValue( n*dx )
yCoords = vtkDoubleArray()
yCoords.InsertNextValue(0.0)
zCoords = xCoords


rgrid = vtkRectilinearGrid()
rgrid.SetDimensions( N+1, 1, N+1 )
rgrid.SetXCoordinates(xCoords)
rgrid.SetYCoordinates(yCoords)
rgrid.SetZCoordinates(zCoords)

slowness = vtkDoubleArray()
slowness.SetName("Slowness")
slowness.SetNumberOfComponents(1)
slowness.SetNumberOfTuples( rgrid.GetNumberOfCells() )

for n in range(0, rgrid.GetNumberOfCells()):
    bo = rgrid.GetCell(n).GetBounds()
    z = floor(bo[4]) + 0.5
    s = 1.0 / (a+b*z)
    slowness.SetTuple1(n, s)

rgrid.GetCellData().SetScalars( slowness );

writer = vtkXMLRectilinearGridWriter()
writer.SetFileName( "./files/layers_coarse2d.vtr" )
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

ugrid = vtkUnstructuredGrid()
pts = vtkPoints()
slowness = vtkDoubleArray()
slowness.SetName('Slowness')

for x in np.arange(0.0, 20.0+0.01, dx):
    for z in np.arange(0.0, 20.0+0.01, dx):
        pts.InsertNextPoint(x, 0.0, z)
        s = 1.0 / (a+b*z)
        slowness.InsertNextValue(s)

ugrid.SetPoints(pts)
ugrid.GetPointData().SetScalars(slowness)

tri = vtkTriangle()

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

writer = vtkXMLUnstructuredGridWriter()
writer.SetFileName('./files/gradient_fine2d.vtu')
writer.SetInputData( ugrid )
writer.SetDataModeToBinary();
writer.Update()


xCoords = vtkDoubleArray()
for n in range(0, N+1):
    xCoords.InsertNextValue( n*dx )
yCoords = vtkDoubleArray()
yCoords.InsertNextValue(0.0)
zCoords = xCoords


rgrid = vtkRectilinearGrid()
rgrid.SetDimensions( N+1, 1, N+1 )
rgrid.SetXCoordinates(xCoords)
rgrid.SetYCoordinates(yCoords)
rgrid.SetZCoordinates(zCoords)

slowness = vtkDoubleArray()
slowness.SetName("Slowness")
slowness.SetNumberOfComponents(1)
slowness.SetNumberOfTuples( rgrid.GetNumberOfPoints() )

for n in range(0, rgrid.GetNumberOfPoints()):
    x = rgrid.GetPoint(n)
    s = 1.0 / (a+b*x[2])
    slowness.SetTuple1(n, s)

rgrid.GetPointData().SetScalars( slowness );

writer = vtkXMLRectilinearGridWriter()
writer.SetFileName( "./files/gradient_fine2d.vtr" )
writer.SetInputData( rgrid )
writer.SetDataModeToBinary()
writer.Update()


# %% layers

ugrid = vtkUnstructuredGrid()
pts = vtkPoints()

for x in np.arange(0.0, 20.0+0.01, dx):
    for z in np.arange(0.0, 20.0+0.01, dx):
        pts.InsertNextPoint(x, 0.0, z)

ugrid.SetPoints(pts)

tri = vtkTriangle()
slowness = vtkDoubleArray()
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

writer = vtkXMLUnstructuredGridWriter()
writer.SetFileName('./files/layers_fine2d.vtu')
writer.SetInputData( ugrid )
writer.SetDataModeToBinary();
writer.Update()


xCoords = vtkDoubleArray()
for n in range(0, N+1):
    xCoords.InsertNextValue( n*dx )
yCoords = vtkDoubleArray()
yCoords.InsertNextValue(0.0)
zCoords = xCoords


rgrid = vtkRectilinearGrid()
rgrid.SetDimensions( N+1, 1, N+1 )
rgrid.SetXCoordinates(xCoords)
rgrid.SetYCoordinates(yCoords)
rgrid.SetZCoordinates(zCoords)

slowness = vtkDoubleArray()
slowness.SetName("Slowness")
slowness.SetNumberOfComponents(1)
slowness.SetNumberOfTuples( rgrid.GetNumberOfCells() )

for n in range(0, rgrid.GetNumberOfCells()):
    bo = rgrid.GetCell(n).GetBounds()
    z = floor(bo[4]) + 0.5
    s = 1.0 / (a+b*z)
    slowness.SetTuple1(n, s)

rgrid.GetCellData().SetScalars( slowness );

writer = vtkXMLRectilinearGridWriter()
writer.SetFileName( "./files/layers_fine2d.vtr" )
writer.SetInputData( rgrid )
writer.SetDataModeToBinary()
writer.Update()

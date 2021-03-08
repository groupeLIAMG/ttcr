#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 08:06:41 2021

@author: giroux
"""

import gmsh
import vtk

import numpy as np
from math import floor

N = 20
dx = 1.0


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

# %% Rcv

with open('rcv3d_in.dat', 'w') as f:
    xx = np.arange(1, 20)
    print('{0:d}'.format(xx.size*xx.size), file=f)
    for x in xx:
        for z in xx:
            print('{0:d} {1:d} {2:d}'.format(x, x, z), file=f)


# %% gradient

xCoords = vtk.vtkDoubleArray()
for n in range(0, N+1):
    xCoords.InsertNextValue( n*dx )
yCoords = vtk.vtkDoubleArray()
yCoords = xCoords
zCoords = xCoords

rgrid = vtk.vtkRectilinearGrid()
rgrid.SetDimensions( N+1, N+1, N+1 )
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
writer.SetFileName( "gradient_medium.vtr" )
writer.SetInputData( rgrid )
writer.SetDataModeToBinary()
writer.Update()


# %%
gmsh.initialize()

# %%
meshSize = 0.43
gmsh.clear()

gmsh.model.occ.addBox(0, 0, 0, 20, 20, 20)

gmsh.model.occ.synchronize()

gmsh.model.mesh.setSize(gmsh.model.getEntities(0), meshSize)

gmsh.model.mesh.generate(3)
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
    if type(cell) is vtk.vtkTetra:
        ugrid2.InsertNextCell(cell.GetCellType(), cell.GetPointIds())

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName('gradient_medium.vtu')
writer.SetInputData( ugrid2 )
writer.SetDataModeToBinary();
writer.Update()


# %%

gmsh.clear()

tags = []
for n in range(20):
    tags.append((3, gmsh.model.occ.addBox(0, 0, n, 20, 20, 1)))
gmsh.model.occ.synchronize()
gmsh.model.occ.fragment(tags, ())
gmsh.model.occ.synchronize()

gmsh.model.mesh.setSize(gmsh.model.getEntities(0), meshSize)

gmsh.model.mesh.generate(3)
gmsh.write('/tmp/mesh.vtk')

# %%
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName('/tmp/mesh.vtk')
reader.Update()

ugrid = reader.GetOutput()
s = []
for n in range(ugrid.GetNumberOfCells()):
    bnds = ugrid.GetCell(n).GetBounds()
    z = floor(bnds[4]) + 0.5
    s.append(1.0 / (a+b*z))

ugrid2 = vtk.vtkUnstructuredGrid()
ugrid2.SetPoints(ugrid.GetPoints())
slowness = vtk.vtkDoubleArray()
slowness.SetName('Slowness')
for n in range(ugrid.GetNumberOfCells()):
    cell = ugrid.GetCell(n)
    if type(cell) is vtk.vtkTetra:
        ugrid2.InsertNextCell(cell.GetCellType(), cell.GetPointIds())
        slowness.InsertNextValue(s[n])
ugrid2.GetCellData().SetScalars(slowness)

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName('layers_medium.vtu')
writer.SetInputData( ugrid2 )
writer.SetDataModeToBinary();
writer.Update()


# %%

rgrid = vtk.vtkRectilinearGrid()
rgrid.SetDimensions( N+1, N+1, N+1 )
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
writer.SetFileName( "layers_medium.vtr" )
writer.SetInputData( rgrid )
writer.SetDataModeToBinary()
writer.Update()

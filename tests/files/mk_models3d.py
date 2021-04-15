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

################################################################################
#
# Medium grid
#
################################################################################

frac = 2
N *= frac
dx /= frac


x = np.arange(0.0, 20.0+0.01, dx)
z = np.arange(0.0, 20.0+0.01, dx)

x_offset = 319640.45
y_offset = 5180561.39

# %% Rcv

with open('rcv3d_in.dat', 'w') as f:
    xx = np.arange(1, 20)
    print('{0:d}'.format(xx.size*xx.size), file=f)
    for x in xx:
        for z in xx:
            print('{0:d} {1:d} {2:d}'.format(x, x, z), file=f)

with open('rcv3d_in2.dat', 'w') as f:
    xx = np.arange(1, 20)
    print('{0:d}'.format(xx.size*xx.size), file=f)
    for x in xx:
        for z in xx:
            print('{0:f} {1:f} {2:d}'.format(x+x_offset, x+y_offset, z), file=f)

with open('src3d_in2.dat', 'w') as f:
    print('1\n{0:f} {1:f} 1.0 0.0'.format(1+x_offset, 1+y_offset), file=f)

################################################################################
#
# %% gradient - Rectilinear

xCoords = vtk.vtkDoubleArray()
yCoords = vtk.vtkDoubleArray()
zCoords = vtk.vtkDoubleArray()
for n in range(0, N+1):
    xCoords.InsertNextValue( n*dx )
yCoords = vtk.vtkDoubleArray()
yCoords.DeepCopy(xCoords)
zCoords.DeepCopy(xCoords)

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

#
# %% shift grid
#
for n in range(0, N+1):
    xCoords.InsertValue(n, n*dx+x_offset)
    yCoords.InsertValue(n, n*dx+y_offset)
rgrid.SetXCoordinates(xCoords)
rgrid.SetYCoordinates(yCoords)

writer = vtk.vtkXMLRectilinearGridWriter()
writer.SetFileName( "gradient_medium2.vtr" )
writer.SetInputData( rgrid )
writer.SetDataModeToBinary()
writer.Update()

for n in range(0, N+1):
    xCoords.InsertValue(n, n*dx)
    yCoords.InsertValue(n, n*dx)
rgrid.SetXCoordinates(xCoords)
rgrid.SetYCoordinates(yCoords)


################################################################################
#
# %% gradient - Unstructured grid
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

#
# %% shift grid
#
tr = vtk.vtkTransform()
tr.Translate(x_offset, y_offset, 0.0)
tf = vtk.vtkTransformFilter()
tf.SetInputData(ugrid2)
tf.SetTransform(tr)
tf.Update()

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName( "gradient_medium2.vtu" )
writer.SetInputData( tf.GetOutput() )
writer.SetDataModeToBinary()
writer.Update()

################################################################################
#
# %%  layers - unstructured grid

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

#
# %% shift grid
#
tr = vtk.vtkTransform()
tr.Translate(x_offset, y_offset, 0.0)
tf = vtk.vtkTransformFilter()
tf.SetInputData(ugrid2)
tf.SetTransform(tr)
tf.Update()

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName( "layers_medium2.vtu" )
writer.SetInputData( tf.GetOutput() )
writer.SetDataModeToBinary()
writer.Update()

################################################################################
#
# %% layers - rectilinear

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

#
# %% shift grid
#
for n in range(0, N+1):
    xCoords.InsertValue(n, n*dx+x_offset)
    yCoords.InsertValue(n, n*dx+y_offset)
rgrid.SetXCoordinates(xCoords)
rgrid.SetYCoordinates(yCoords)

writer = vtk.vtkXMLRectilinearGridWriter()
writer.SetFileName( "layers_medium2.vtr" )
writer.SetInputData( rgrid )
writer.SetDataModeToBinary()
writer.Update()

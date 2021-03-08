#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 13:44:15 2019

@author: giroux
"""

import numpy as np
import vtk


# %%

N = 20
x = np.arange(N+1)
y = np.arange(N+1)
z = np.arange(N+1)

Va = 1.0
V20 = 3.0
b = (V20-Va)/20.0

# src is at 0,0,0

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
ttData.SetNumberOfTuples(x.size * y.size * z.size)

for ix in range(len(x)):
    for iy in range(len(y)):
        for iz in range(len(z)):

            Vb = Va+b*z[iz]
            r = np.sqrt(x[ix]*x[ix] + y[iy]*y[iy] + z[iz]*z[iz])
            tt = np.abs(np.arccosh(1.+(b*b*r*r)/(2.*Va*Vb))/b)

            ii = rgrid.FindPoint(x[ix], y[iy], z[iz])
            ttData.SetTuple1(ii, tt)
rgrid.GetPointData().SetScalars(ttData)

writer = vtk.vtkXMLRectilinearGridWriter()
writer.SetFileName('sol_analytique_gradient_tt.vtr')
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
        Vb = Va+b*z[iz]
        r = np.sqrt(x[ix]*x[ix] + z[iz]*z[iz])
        tt = np.abs(np.arccosh(1.+(b*b*r*r)/(2.*Va*Vb))/b)
        ii = rgrid.FindPoint(x[ix], 0, z[iz])
        ttData.SetTuple1(ii, tt)
rgrid.GetPointData().SetScalars(ttData)

writer = vtk.vtkXMLRectilinearGridWriter()
writer.SetFileName('sol_analytique_gradient2d_tt.vtr')
writer.SetInputData(rgrid)
writer.SetDataModeToBinary()
writer.Update()

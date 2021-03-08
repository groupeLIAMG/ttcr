#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 09:44:16 2021

@author: giroux
"""
import subprocess
import numpy as np

import vtk


N = 21

x = np.arange(N)
z = np.arange(N)

tt = np.zeros((N, N))

for h in z:
    for km in x:
        if h == 0 and km == 0:
            continue
        if h == 0:
            hh = km * 0.04 / 20.0
        else:
            hh = h
        cmd = ['taup_time', '-mod', 'couches', '-h', str(hh), '-km', str(km),
                '-ph', 'p', '--time']

        result = subprocess.run(cmd, stdout=subprocess.PIPE,
                                env={"PATH": "/usr/local/TauP/bin:/usr/bin/java:/usr/bin"})
        try:
            tt[km, h] = float(result.stdout)
        except ValueError:
            tt[km, h] = np.nan
            print(result)


xCoords = vtk.vtkFloatArray()
for i in x:
    xCoords.InsertNextValue(i)

yCoords = vtk.vtkFloatArray()
yCoords.InsertNextValue(0)

zCoords = vtk.vtkFloatArray()
for i in z:
    zCoords.InsertNextValue(i)

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
        ttData.SetTuple1(ii, tt[ix,iz])
rgrid.GetPointData().SetScalars(ttData)

writer = vtk.vtkXMLRectilinearGridWriter()
writer.SetFileName('sol_analytique_couches2d_taup_tt.vtr')
writer.SetInputData(rgrid)
writer.SetDataModeToBinary()
writer.Update()

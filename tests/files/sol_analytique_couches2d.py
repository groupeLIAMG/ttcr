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

a = 1.0
V20 = 3.0
b = (V20-a)/20.0
v = []
for j in np.arange(N+1):
    zc = j+0.5
    v.append(a + b*zc)

with open('couches.tvel', 'w') as f:
    print('/* Modèle en couche */', file=f)
    print('/* Solution avec TauP */', file=f)
    for j in np.arange(N+1):
        print('{0:f} {1:f} {2:f} {3:f}'.format(j, v[j], 0.6*v[j], 2.67), file=f)

cmd = ['taup_create', '-tvel', 'couches.tvel']
subprocess.run(cmd, stdout=subprocess.PIPE,
               env={"PATH": "/usr/local/TauP/bin:/usr/bin/java:/usr/bin"}) 
for h in z:
    for km in x:
        if h == 0 and km == 0:
            continue
        if h == 0:
            hh = km * 0.04 / 20.0
        else:
            hh = h
        cmd = ['taup_time', '-mod', 'couches', '-h', '0', '--stadepth', str(hh), '-km', str(km),
                '-ph', 'p', '--time']

        result = subprocess.run(cmd, stdout=subprocess.PIPE,
                                env={"PATH": "/usr/local/TauP/bin:/usr/bin/java:/usr/bin"})
        print(cmd)
        print(h, km, result)
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

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate constant-velocity rectilinear models (medium and fine resolution)
used by the accuracy assessment of the 3D rectilinear grid methods.

Same 20 x 20 x 20 domain and node counts as the layers/gradient models
produced by mk_models3d.py:
    medium : frac = 2  ->  41^3 nodes, dx = 0.5
    fine   : frac = 8  -> 161^3 nodes, dx = 0.125

The velocity is constant, so the analytic travel time from any source is
simply  t = slowness * euclidean_distance, which the C++ program computes
directly (no reference .vtr file needed).
"""

import vtk

# constant velocity (km/s); slowness = 1/V is written at the grid nodes
V0 = 3.0


def make_model(frac, filename):
    N = 20 * frac
    dx = 1.0 / frac

    coords = vtk.vtkDoubleArray()
    for n in range(0, N + 1):
        coords.InsertNextValue(n * dx)
    ycoords = vtk.vtkDoubleArray()
    ycoords.DeepCopy(coords)
    zcoords = vtk.vtkDoubleArray()
    zcoords.DeepCopy(coords)

    rgrid = vtk.vtkRectilinearGrid()
    rgrid.SetDimensions(N + 1, N + 1, N + 1)
    rgrid.SetXCoordinates(coords)
    rgrid.SetYCoordinates(ycoords)
    rgrid.SetZCoordinates(zcoords)

    slowness = vtk.vtkDoubleArray()
    slowness.SetName("Slowness")
    slowness.SetNumberOfComponents(1)
    slowness.SetNumberOfTuples(rgrid.GetNumberOfPoints())
    for n in range(0, rgrid.GetNumberOfPoints()):
        slowness.SetTuple1(n, 1.0 / V0)
    rgrid.GetPointData().SetScalars(slowness)

    writer = vtk.vtkXMLRectilinearGridWriter()
    writer.SetFileName(filename)
    writer.SetInputData(rgrid)
    writer.SetDataModeToBinary()
    writer.Update()
    print(f"wrote {filename}  ({N + 1}^3 nodes, dx = {dx})")


if __name__ == "__main__":
    make_model(2, "constant_medium.vtr")
    make_model(8, "constant_fine.vtr")

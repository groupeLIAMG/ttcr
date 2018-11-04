# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 09:48:08 2016

@author: giroux
"""

import math
import re
import sys

from numba import jit
import numpy as np
from scipy.sparse import csr_matrix
import h5py

import vtk

import cmesh3d

from utils import nargout


class MSHReader:
    def __init__(self, filename):
        self.filename = filename
        self.valid = self.checkFormat()

    def checkFormat(self):
        try:
            f = open(self.filename, 'r')
            line = f.readline()
            if not line.startswith('$MeshFormat'):
                return False
            tmp = re.split(r' ', f.readline())
            if tmp[0] != '2.2' or tmp[1] != '0':
                return False

        except OSError as err:
            print("OS error: {0}".format(err))
            return False
        except ValueError:
            print("Could not convert data to an integer.")
            return False
        except:
            print("Unexpected error:", sys.exc_info()[0])
            return False

    def is2D(self):
        try:
            nodes = self.readNodes()
        except:
            raise

        ymin = np.min(nodes[:,1])
        ymax = np.max(nodes[:,1])
        return ymin==ymax

    def getNumberOfNodes(self):
        try:
            f = open(self.filename, 'r')
            for line in f:
                if line.startswith('$Nodes'):
                    nnodes = int(f.readline())
                    break
            f.close()
            return nnodes

        except OSError as err:
            print("OS error: {0}".format(err))
        except ValueError:
            print("Could not convert data to an integer.")
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise

    def readNodes(self):
        try:
            f = open(self.filename, 'r')
            for line in f:
                if line.startswith('$Nodes'):
                    nnodes = int(f.readline())
                    nodes = np.zeros((nnodes,3))
                    for n in np.arange(nnodes):
                        tmp = re.split(r' ', f.readline())
                        nodes[n,0] = tmp[1]
                        nodes[n,1] = tmp[2]
                        nodes[n,2] = tmp[3]

                    break
            f.close()
            return nodes

        except OSError as err:
            print("OS error: {0}".format(err))
        except ValueError:
            print("Could not convert data to an integer.")
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise


    def getNumberOfElements(self, eltype=0):
        try:
            f = open(self.filename, 'r')
            for line in f:
                if line.startswith('$Elements'):
                    nelem = int(f.readline())
                    if eltype != 0:
                        nel = 0
                        for n in np.arange(nelem):
                            tmp = re.split(r' ', f.readline())
                            if eltype == int(tmp[1]):
                                nel += 1
                        nelem = nel
                    break
            f.close()
            return nelem

        except OSError as err:
            print("OS error: {0}".format(err))
        except ValueError:
            print("Could not convert data to an integer.")
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise

    def readTriangleElements(self):
        try:
            f = open(self.filename, 'r')
            for line in f:
                if line.startswith('$Elements'):
                    nelem = int(f.readline())
                    triangles = np.ndarray((nelem,3), dtype = np.int64)
                    nt = 0
                    for n in np.arange(nelem):
                        tmp = re.split(r' ', f.readline())
                        if tmp[1] == '2':
                            nTags = int(tmp[2])
                            triangles[nt,0] = tmp[nTags+3]
                            triangles[nt,1] = tmp[nTags+4]
                            triangles[nt,2] = tmp[nTags+5]
                            nt += 1
                    break
            return triangles[:nt,:]-1  # indices start at 0 in python

        except OSError as err:
            print("OS error: {0}".format(err))
        except ValueError:
            print("Could not convert data to an integer.")
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise

    def readTetraherdonElements(self):
        try:
            f = open(self.filename, 'r')
            for line in f:
                if line.startswith('$Elements'):
                    nelem = int(f.readline())
                    tetrahedra = np.ndarray((nelem,4), dtype = np.int64)
                    nt = 0
                    for n in np.arange(nelem):
                        tmp = re.split(r' ', f.readline())
                        if tmp[1] == '4':
                            nTags = int(tmp[2])
                            tetrahedra[nt,0] = tmp[nTags+3]
                            tetrahedra[nt,1] = tmp[nTags+4]
                            tetrahedra[nt,2] = tmp[nTags+5]
                            tetrahedra[nt,3] = tmp[nTags+6]
                            nt += 1
                    break
            return tetrahedra[:nt,:]-1  # indices start at 0 in python

        except OSError as err:
            print("OS error: {0}".format(err))
        except ValueError:
            print("Could not convert data to an integer.")
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise







class Mesh:
    """
    Superclass for 2D and 3D Meshes
    """

    def __init__(self):
        self.nodes = np.ndarray([])
        self.cells = np.ndarray([])   # can be triangles (2D/3D) or tetrahedra (3D)

    def getNumberOfNodes(self):
        return self.nodes.shape[0]

    def getNumberOfCells(self):
        return self.cells.shape[0]





class meshTriangle3D(Mesh):
    def __init__(self, nod, tri):
        self.cells = np.array(tri)
        itr = np.unique(tri)
        for n in np.arange(itr.size):
            i,j = np.nonzero(tri == itr[n])
            self.cells[i,j] = n
        self.nodes = np.array(nod[itr,:])

    def projz(self, pt):
        ind = self.ni(pt)
        if np.size( ind ) > 0:
            return np.array(self.nodes[ind,:])
        else:
            p = np.array(pt)
            p[2] = 0.0

            d = np.sqrt(  (self.nodes[:,0]-p[0])**2 + (self.nodes[:,1]-p[1])**2 )
            iclosest = np.argmin(d)
            # find second closest
            d[iclosest] *= 1e6
            isecond = np.argmin(d)
            itri,jtri = np.nonzero(np.logical_or(self.cells == iclosest, self.cells == isecond))

            for i in itri:
                a = np.array(self.nodes[self.cells[i,0],:])
                b = np.array(self.nodes[self.cells[i,1],:])
                c = np.array(self.nodes[self.cells[i,2],:])
                a[2] = 0.0
                b[2] = 0.0
                c[2] = 0.0
                if meshTriangle3D.insideTriangle(p, a, b, c) == True:
                    a[:] = self.nodes[self.cells[i,0],:]
                    b[:] = self.nodes[self.cells[i,1],:]
                    c[:] = self.nodes[self.cells[i,2],:]

#                    import matplotlib.pyplot as plt
#                    fig = plt.figure()
#                    ax = fig.add_subplot(1,1,1)
#                    ax.plot([a[0], b[0], c[0], a[0]], [a[1], b[1], c[1], a[1]])
#                    ax.hold(True)
#                    ax.plot(p[0], p[1],'rv')
#                    ax.hold(False)
#                    plt.show()

                    p[2] = np.min(self.nodes[:,2]) - 100.0  # make sure pt is below mesh
                    # now project vertically
                    v = np.array([0.0, 0.0, 1.0])  # point upward
                    pa = a-p
                    pb = b-p
                    pc = c-p

                    A = np.vstack((pa, pb, pc)).T
                    x = np.linalg.solve(A, v)
                    G = (x[0]*a + x[1]*b + x[2]*c)/np.sum(x)
                    return G
            else:
                
                import matplotlib.pyplot as plt
                fig = plt.figure()
                ax = fig.add_subplot(1,1,1)
                ax.plot(p[0], p[1],'rv')
                ax.axis('equal')
                ax.hold(True)
                ax.plot(self.nodes[iclosest,0], self.nodes[iclosest,1],'ko')
                for i in itri:
                    a = np.array(self.nodes[self.cells[i,0],:])
                    b = np.array(self.nodes[self.cells[i,1],:])
                    c = np.array(self.nodes[self.cells[i,2],:])
                    a[2] = 0.0
                    b[2] = 0.0
                    c[2] = 0.0
                    ax.plot([a[0], b[0], c[0], a[0]], [a[1], b[1], c[1], a[1]])
                ax.hold(False)
                plt.show()
                
                
                print('Houston! We have a problem')



    def ni(self, pt):
        """
        Check if pt is one of the nodes
        """
        ind = []
        for n in np.arange(self.nodes.shape[0]):
            if np.sqrt( np.sum( (self.nodes[n,:]-pt)**2 ) ) < 0.0001:
                ind.append(n)
        return ind


    def getVtkUnstructuredGrid(self):
        tPts = vtk.vtkPoints()
        ugrid = vtk.vtkUnstructuredGrid()
        tPts.SetNumberOfPoints(self.nodes.shape[0])
        for n in np.arange(self.nodes.shape[0]):
            tPts.InsertPoint(n, self.nodes[n,0], self.nodes[n,1], self.nodes[n,2])
        ugrid.SetPoints(tPts)
        tri = vtk.vtkTriangle()
        for n in np.arange(self.cells.shape[0]):
            tri.GetPointIds().SetId(0, self.cells[n,0])
            tri.GetPointIds().SetId(1, self.cells[n,1])
            tri.GetPointIds().SetId(2, self.cells[n,2])
            ugrid.InsertNextCell( tri.GetCellType(), tri.GetPointIds() )

        return ugrid


    @staticmethod
    @jit
    def insideTriangle(p, a, b, c):
        a -= p
        b -= p
        c -= p
        u = np.cross(b, c)
        v = np.cross(c, a)
        if np.dot(u, v) < 0.0:
            return False
        w = np.cross(a, b)
        if np.dot(u, w) < 0.0:
            return False
        return True


class MeshTetrahedra(Mesh):

    def __init__(self, nthreads=1):
        Mesh.__init__(self)
        self.nthreads = nthreads
        self.cmesh = None

    @jit
    def buildFromMSH(self, filename):
        reader = MSHReader(filename)
        nodes = reader.readNodes()
        tet = reader.readTetraherdonElements()
        self.cells = np.array(tet)
        itet = np.unique(tet)
        for n in np.arange(itet.size):
            i,j = np.nonzero(tet == itet[n])
            self.cells[i,j] = n
        self.nodes = np.array(nodes[itet,:])


    def raytrace(self, slowness, Tx, Rx, t0=()):
        nout = nargout()
        if nout != 1 and nout != 3 and nout != 4:
            raise SyntaxError('MeshTetrahedra.raytrace: 1, 3 or 4 output arguments allowed')
        # check input data consistency

        if Tx.ndim != 2 or Rx.ndim != 2:
            raise ValueError('Tx and Rx should be 2D arrays')

        if Tx.shape[1] != 3 or Rx.shape[1] != 3:
            raise ValueError('Tx and Rx should be ndata x 3')

        if Tx.shape != Rx.shape:
            raise ValueError('Tx and Rx should be of equal size')

        if len(slowness) != self.getNumberOfNodes():
            raise ValueError('Length of slowness vector should equal number of nodes')

        if len(t0) == 0:
            t0 = np.zeros([Tx.shape[0],])
        elif len(t0) != Tx.shape[0]:
            raise ValueError('Length of t0 should equal number of Tx')

        if self.cmesh == None:
            eps = 1.0e-15
            rp_method = 1
            maxit = 20
            self.cmesh = cmesh3d.Mesh3Dcpp(self.nodes, self.cells, eps, maxit, rp_method, self.nthreads)

        if nout == 1:
            tt = self.cmesh.raytrace(slowness, Tx, Rx, t0)
            return tt
        elif nout == 3:
            tt, rays, v0 = self.cmesh.raytrace(slowness, Tx, Rx, t0)
            return tt, rays, v0
        elif nout == 4:
            tt, rays, v0, M = self.cmesh.raytrace(slowness, Tx, Rx, t0)
            return tt, rays, v0, M





if __name__ == '__main__':
    reader = MSHReader('topo.msh')

    print(reader.is2D())
    print(reader.getNumberOfElements())
    print(reader.getNumberOfElements(2))

    nodes = reader.readNodes()
    tri = reader.readTriangleElements()
    print(tri.shape)

    import mayavi.mlab as mm
    mm.figure()
    mm.triangular_mesh(nodes[:,0], nodes[:,1], nodes[:,2], tri, color=(0, 0, 1))
    mm.show()



    mesh = meshTriangle3D(nodes, tri)

    p = np.array([305105.0, 5723733.0, 100.0])

    mesh.projz(p)

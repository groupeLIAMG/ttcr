# -*- coding: utf-8 -*-

import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

# for parallel computations
from multiprocessing import Process, Queue

# ttcrpy must be in your PYTHONPATH
from ttcrpy import mesh, cmesh2d


# Some functions used later

def f(x, y):
    '''
    slowness of Model 2 in Qian et al., 2007
    '''
    return 2.*np.pi*np.sqrt( np.cos(2.*np.pi*x)*np.sin(2.*np.pi*y)*np.cos(2.*np.pi*x)*np.sin(2.*np.pi*y) +
							 np.sin(2.*np.pi*x)*np.cos(2.*np.pi*y)*np.sin(2.*np.pi*x)*np.cos(2.*np.pi*y) )

def _rt_worker(idx, mesh, istart, iend, Tx, Rx, iRx, t0, nout, tt_queue):
    """
    worker for spanning raytracing to different processes
    """
    # slowness must have been set before raytracing
    for n in range(istart, iend):
        t = mesh.raytrace(None,
                          np.atleast_2d(Tx[n, :]),
                          np.atleast_2d(Rx),
                          t0[n],
                          nout=nout,
                          thread_no=idx)
        tt_queue.put((t, iRx[n], n))
    tt_queue.close()



# Read mesh
reader = mesh.MSHReader('model2_qian.msh')

nodes = reader.readNodes()
triangles = reader.readTriangleElements()

# *** Important ***
#  keep only x and z coordinates of nodes

nodes = np.vstack((nodes[:,0], nodes[:,2])).T

# get centroid of triangles

x = []
z = []
for t in triangles:
    x.append((nodes[t[0],0] + nodes[t[1],0] + nodes[t[2],0])/3.0)
    z.append((nodes[t[0],1] + nodes[t[1],1] + nodes[t[2],1])/3.0)

# compute slowness at centroids, these values will be assigned to triangles
slowness = f(np.array(x), np.array(z))

# %%
# create mesh for serial calculations

m1 = cmesh2d.Mesh2D(nodes, triangles, nsecondary=5, nthreads=1)

# Source in Qian et al 2007 (5 points each having its t0)
Tx = np.array([[0.25, 0.25],
               [0.75, 0.75],
               [0.25, 0.75],
               [0.25, 0.25],
               [0.5, 0.5]])
t0 = np.array([1.0, 1.0, -1.0, -1.0, 0.0])

# Compute tt at all nodes
Rx = nodes

# %%
tt = m1.raytrace(slowness, Tx, Rx, t0)

# %% show traveltimes

plt.scatter(nodes[:,0], nodes[:,1], c=tt, cmap='seismic')
plt.colorbar()
plt.show(block=False)


# %% Second example showing how to use the data kernel matrix

slowness[:] = 1.0
m1.set_slowness(slowness)

# Create crosshole-like acquisition geometry
nTx = 17
nRx = 35
Tx = np.vstack((0.1+np.zeros((nTx,)), np.linspace(0.1, 0.9, nTx))).T
Rx = np.vstack((0.9+np.zeros((nRx,)), np.linspace(0.1, 0.9, nRx))).T
t0 = np.zeros((nTx,))

tt = np.zeros((nTx*nRx,))
L = []
# loop over sources
for n in range(Tx.shape[0]):
    ind = n*nRx + np.arange(nRx)
    tt[ind], LL = m1.raytrace(slowness, np.atleast_2d(Tx[n, :]), Rx, t0, nout=2)
    L.append(LL)
# assemble matrices
L = sp.vstack(L)

# Compute times with data kernel matrix
tt2 = L.dot(slowness)

# Check results
plt.figure()
plt.plot(tt, 'o')
plt.plot(tt2)
plt.show(block=False)

# %% Example parallel Computation

nthreads = 4

# We must create a new mesh with the right number of threads
m2 = cmesh2d.Mesh2D(nodes, triangles, nsecondary=10, nthreads=nthreads)

tt = np.zeros((nTx*nRx,))
t0 = np.zeros((nTx,))

# with parallel calculations, we don't know in which order transmitters will be
# processed, we have to build a vector of Rx indices, to store tt at the correct
# position
iRx = []
for n in range(nTx):
    iRx.append(n*nRx + np.arange(nRx))

# size of parallel blocks
blk_size = np.zeros((nthreads,), dtype=np.int64)
nj = nTx
while nj > 0:
    for n in range(nthreads):
        blk_size[n] += 1
        nj -= 1
        if nj == 0:
            break

m2.set_slowness(slowness)

nout = 3
if nout >=2:
    L = [ [] for i in range(nTx) ]
if nout == 3:
    rays = [ [] for n in range(nTx)]


processes = []
blk_start = 0
tt_queue = Queue()  # we need a queue to retrieve the results

# distribute jobs over threads
for n in range(nthreads):
    blk_end = blk_start + blk_size[n]
    p = Process(target=_rt_worker,
                args=(n, m2, blk_start, blk_end, Tx, Rx,
                      iRx, t0, nout, tt_queue),
                daemon=True)
    processes.append(p)
    p.start()
    blk_start += blk_size[n]

if nout == 1:
    for i in range(nTx):
        t, ind, n = tt_queue.get()
        tt[ind] = t
if nout == 2:
    for i in range(nTx):
        t, ind, n = tt_queue.get()
        tt[ind] = t[0]
        L[n] = t[1]
    L = sp.vstack(L)

elif nout == 3:
    for i in range(nTx):
        t, ind, n = tt_queue.get()
        tt[ind] = t[0]
        L[n] = t[1]
        rays[n] = t[2]
    L = sp.vstack(L)


plt.figure()
plt.plot(tt, 'o')
if nout >= 2:
    tt2 = L.dot(slowness)
    plt.plot(tt2)
if nout == 3:
    plt.figure()
    for nt in range(nTx):
        for nr in range(nRx):
            plt.plot(rays[nt][nr][:,0], rays[nt][nr][:,1], 'k')
plt.show()

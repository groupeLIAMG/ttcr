# -*- coding: utf-8 -*-

import numpy as np

import ttcrpy.rgrid as rg

x = np.arange(5.0)
y = np.arange(1.0, 8.0, 1.1)
z = np.arange(-1.0, 4.0, 1.2)
g1 = rg.Grid3d(x, y, z, method='SPM')

slowness = np.ones((g1.nparams,))

src = np.array([[0.5, 1.5, 0.0]])
rcv = np.array([[3.5, 1.5, 0.0],
                [3.5, 1.5, 1.0],
                [3.5, 1.5, 2.0]])

print(g1.nparams, g1.shape, src.ndim, rcv.ndim)

print('\nRaytrace 1')
tt = g1.raytrace(src, rcv, slowness)
print(tt)

src = np.array([[0.5, 1.5, 0.0],
                [0.5, 1.5, 1.0]])

try:
    print('\nRaytrace 2')
    g1.raytrace(src, rcv, slowness)
except ValueError:
    print('Exception caught: OK')

print('\nRaytrace 3')
tt = g1.raytrace(src, rcv, slowness, aggregate_src=True)
print(tt)

src = np.kron(src, np.ones((3,1)))
rcv = np.kron(np.ones((2,1)), rcv)

print('\nRaytrace 4')
tt = g1.raytrace(src, rcv, slowness)
print(tt)

src = np.array([[5.0, 0.5, 1.5, 0.0],
                [6.0, 0.5, 1.5, 1.0]])

src = np.kron(src, np.ones((3,1)))

print('\nRaytrace 5')
tt = g1.raytrace(src, rcv, slowness)
print(tt)

src = src[:-1,:]
rcv = rcv[:-1,:]
print('\nRaytrace 6')
tt = g1.raytrace(src, rcv, slowness)
print(tt)

src = np.array([[99, 5.5, 0.5, 1.5, 0.0],
                [105, 6.5, 0.5, 1.5, 1.0]])

src = np.kron(src, np.ones((3,1)))
src = src[:-1,:]

print('\nRaytrace 7')
tt = g1.raytrace(src, rcv, slowness)
print(tt)

src = np.roll(src, 1, axis=0)
rcv = np.roll(rcv, 1, axis=0)

print('\nRaytrace 8')
tt = g1.raytrace(src, rcv, slowness)
print(tt)

src = np.array([[5.0, 0.5, 1.5, 0.0],
                [5.2, 0.5, 1.5, 0.2],
                [5.4, 0.5, 1.5, 0.4],
                [5.6, 0.5, 1.5, 0.6],
                [5.8, 0.5, 1.5, 0.8],
                [6.0, 0.5, 1.5, 1.0]])
nsrc = src.shape[0]
rcv = np.array([[3.5, 1.5, 0.0],
                [3.5, 1.5, 1.0],
                [3.5, 1.5, 2.0]])
nrcv = rcv.shape[0]

src = np.kron(src, np.ones((nrcv,1)))
rcv = np.kron(np.ones((nsrc,1)), rcv)

print('\nRaytrace 9')
tt = g1.raytrace(src, rcv, slowness)
print(tt)


g1 = rg.Grid3d(x, y, z, method='SPM', nthreads=2)
print('\nRaytrace 10 (parallel)')
tt = g1.raytrace(src, rcv, slowness)
print(tt)

slowness[:5] = 2.
toto = np.zeros((g1.get_number_of_nodes(),))
toto[:7] = 1.

g1.to_vtk({'slowness': slowness, 'toto': toto}, 'testVTK')

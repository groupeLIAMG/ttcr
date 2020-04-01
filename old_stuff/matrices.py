# -*- coding: utf-8 -*-

import platform
import sys
if platform.system() == 'Darwin':
    sys.path.insert(0, "/Users/giroux/GitHub/hypopy")
    sys.path.insert(0, "/Users/giroux/src/ttcr")
    sys.path.insert(0, "/Users/giroux/src/ttcr/ttcrpy")
elif platform.system() == 'Linux':
    sys.path.insert(0, "/home/giroux/jacques/GitHub/hypopy")
    sys.path.insert(0, "/Users/giroux/src/ttcr")

import time
import numpy as np
import matplotlib.pyplot as plt

import scipy.sparse as sp
import ttcrpy.rgrid as rg
import hypo

np.random.seed(1)

xmin = 90.0
xmax = 211.0
ymin = 80.0
ymax = 211.0
zmin = 0.0
zmax = 101.0

dx = 5.0   # grid cell size, we use cubic cells here

x = np.arange(xmin, xmax, dx)
y = np.arange(ymin, ymax, dx)
z = np.arange(zmin, zmax, dx)

# np.savetxt('x.dat', x)
# np.savetxt('y.dat', y)
# np.savetxt('z.dat', z)

# %%

weno = True

nthreads = 8   # do calculations in parallel
gh = hypo.Grid3D(x, y, z, nthreads=nthreads)
gr = rg.Grid3d(x, y, z, nthreads=nthreads, method='FSM', cell_slowness=False, weno=weno)


rcv = np.array([[112., 115., 13.],
                [111., 116., 40.],
                [111., 113., 90.],
                [151., 117., 17.],
                [180., 115., 16.],
                [113., 145., 11.],
                [160., 150., 17.],
                [185., 149., 15.],
                [117., 184., 11.],
                [155., 192.,  9.],
                [198., 198., 10.],
                [198., 196., 40.],
                [198., 193., 90.]])
ircv = np.arange(rcv.shape[0]).reshape(-1,1)   # vector of rcv indices
nsta = rcv.shape[0]

nev = 20
src = np.vstack((np.arange(nev),                                            # event ID
                 np.linspace(0., 50., nev) + np.random.randn(nev),          # origin time
                 160. + 10.*np.random.randn(nev),                           # x
                 140. + 10.*np.random.randn(nev),                           # y
                  60. + 10.*np.random.randn(nev))).T                        # z

# np.savetxt('src.dat', src)
# np.savetxt('rcv.dat', rcv)

# %%
h_true = src.copy()

def Vz(z):
    return 4000.0 + 10.0*(z-50.0)
Vp = np.kron(Vz(z), np.ones((gr.shape[0], gr.shape[1], 1)))

slowness = 1./Vp.flatten()
# np.savetxt('slowness.dat', slowness)

src = np.kron(src,np.ones((nsta,1)))   # use kron to replicate src-rcv pairs correctly
rcv_data = np.kron(np.ones((nev,1)), rcv)
ircv_data = np.kron(np.ones((nev,1)), ircv)


Kh = sp.hstack(gh.computeK())
Kr = sp.hstack(gr.compute_K())

d = Kh - Kr
print('K', d.count_nonzero())


Vpts = np.array([[Vz(1.), 100.0, 100.0, 1.],
                 [Vz(1.), 100.0, 200.0, 1.],
                 [Vz(1.), 200.0, 100.0, 1.],
                 [Vz(1.), 200.0, 200.0, 1.],
                 [Vz(11.), 112.0, 148.0, 11.0],
                 [Vz(5.), 152.0, 108.0, 5.0],
                 [Vz(75.), 152.0, 108.0, 75.0],
                 [Vz(11.), 192.0, 148.0, 11.0]])

Dh = gh.computeD(Vpts[:,1:])
Dr = gr.compute_D(Vpts[:,1:])
d = Dh - Dr
print('D', d.count_nonzero())


print('### 1 ###')
t1 = time.time()
tth, raysh, Mh = gh.raytrace(slowness, src, rcv_data)
# tth, raysh = gh.raytrace(slowness, src, rcv_data)

print('Compute time: ', time.time()-t1)
print('### 2 ###')
t1 = time.time()
ttr, raysr, Mr = gr.raytrace(slowness, src, rcv_data, compute_M=True,
                                return_rays=True)
# ttr, raysr = gr.raytrace(slowness, src, rcv_data, return_rays=True)
print('Compute time: ', time.time()-t1)

checkRes = True

if checkRes:

    d = np.abs(tth - ttr)
    e_rel = d / tth


    plt.figure()
    plt.plot(tth,'o')
    plt.plot(ttr,'*')
    plt.xlabel('Data number')
    plt.ylabel('Traveltime')
    plt.show(block=False)

    plt.figure()
    plt.plot(100*e_rel)
    plt.show()


    for n in range(len(raysh)):
        d = np.abs(raysh[n] - raysr[n]).flatten()
        if d.max() > 0.0:
            print(n, d.max())



    Mr = sp.vstack(Mr)
    Mh = sp.vstack(Mh)
    d = Mh - Mr
    print('M', d.count_nonzero())
    if d.count_nonzero() > 0:
        ind = d.nonzero()
        e_rel = np.array(d[ind]) / np.array(Mh[ind])

        plt.figure()
        plt.plot(np.array(d[ind]).flatten())
        plt.plot(np.array(Mh[ind]).flatten())
        plt.show(block=False)

        plt.figure()
        plt.plot(100*e_rel.flatten())
        plt.show()

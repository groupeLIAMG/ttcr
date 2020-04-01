#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from multiprocessing import Process, Queue

from ttcrpy.cgrid3d import Grid3Drn


def worker(idx, grid, istart, iend, vTx, Rx, iRx, t0, nout, tt_queue):
    # slowness must have been set before raytracing
    for n in range(istart, iend):
        t = grid.raytrace(None,
                          np.atleast_2d(vTx[n, :]),
                          np.atleast_2d(Rx[iRx[n], :]),
                          t0[n],
                          nout=nout,
                          thread_no=idx)
        tt_queue.put((t, iRx[n], n))


nx = 50
ny = 15
nz = 75

dx = 1.0
x0 = 0.0
y0 = 5.0
z0 = -5.0

num_threads = 2

g = Grid3Drn(nx, ny, nz, dx, x0, y0, z0, 1.e-15, 15, 1, num_threads)

slowness = np.ones(((nx+1)*(ny+1)*(nz+1)))

Tx = np.array([[x0+1, y0+1, z0+1],
               [x0+10*dx, y0+1, z0+1],
               [x0+20*dx, y0+1, z0+1],
               [x0+10*dx, y0+5*dx, z0+1],
               [x0+20*dx, y0+10*dx, z0+1]])

Rx = np.array([[x0+45*dx, y0+1, z0+1],
               [x0+45*dx, y0+dx, z0+dx],
               [x0+45*dx, y0+2*dx, z0+70*dx]])

t0 = np.arange(Tx.shape[0]).reshape(-1, 1)


src = np.kron(Tx, np.ones((Rx.shape[0], 1)))
t0 = np.kron(t0, np.ones((Rx.shape[0], 1)))
rcv = np.kron(np.ones((Tx.shape[0], 1)), Rx)





vTx, i0, ind = np.unique(src, axis=0, return_index=True, return_inverse=True)
nTx = vTx.shape[0]
iRx = []
uind = np.unique(ind)
t0 = t0[i0]
for i in uind:
    ii = ind == i
    iRx.append(ii)


for nr in range(2):

    print('n =', g.get_nthreads())

    g.set_slowness(slowness)
    tt = np.zeros((rcv.shape[0],))
    v0 = np.zeros((rcv.shape[0],))
    rays = [ [0.0] for n in range(rcv.shape[0])]

    if nTx < 1.5*num_threads or num_threads == 1:
        for n in range(nTx):
            t, r, v = g.raytrace(None,
                           np.atleast_2d(vTx[n, :]),
                           np.atleast_2d(rcv[iRx[n], :]),
                           t0[n],
                           nout=3)
            ind = np.where(iRx[n])[0]
            tt[iRx[n]] = t
            v0[iRx[n]] = v
            for nn in range(len(ind)):
                rays[ind[nn]] = r[nn]
    else:

        nout = 3
        blk_size = np.zeros((num_threads,), dtype=np.int64)
        nj = nTx
        while nj > 0:
            for n in range(num_threads):
                blk_size[n] += 1
                nj -= 1
                if nj == 0:
                    break

        processes = []
        blk_start = 0
        tt_queue = Queue(nTx)
        for n in range(num_threads):
            blk_end = blk_start + blk_size[n]
            p = Process(target=worker,
                        args=(n, g, blk_start, blk_end, vTx, rcv,
                              iRx, t0, nout, tt_queue))
            p.start()
            processes.append(p)
            blk_start += blk_size[n]

        for p in processes:
            p.join()

        while not tt_queue.empty():
            t, ind, n = tt_queue.get()
            if nout == 1:
                tt[ind] = t
            else:
                tt[ind] = t[0]
                v0[ind] = t[2]


    d = np.sqrt(np.sum((src - rcv)**2, axis=1))

    print(np.vstack((d, tt)).T)

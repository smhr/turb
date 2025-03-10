#!/usr/bin/env python3
import sys
import os

# Add the parent directory to the system path
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.insert(0, parent_dir)

import numpy as np
import matplotlib.pyplot as plt
from make_turb import TurbField

kmin = 2 # lambda_max, minimum populated turbulent wavenumber for Gaussian initial velocity field, in units of pi/R [default: 2]
turb_sol = 0.5 # fraction of turbulence in solenoidal modes [default: 0.5]
seed = 401
res = 64
N = 2*res # N is 2x cube resolution. E.g. N = 128 means a 64^3 velocity cube.

fname = "vturb%d_sol%g_seed%d.npy"%(kmin,turb_sol, seed)
fname_v = "vturb%d_sol%g_seed%d.bin"%(kmin,turb_sol, seed)

fname_vx = "vx.bin"
fname_vy = "vy.bin"
fname_vz = "vz.bin"

vt = TurbField(N=N, kmin=kmin, sol_weight=turb_sol, seed=seed)
print(vt.shape)
nmin, nmax = vt.shape[-1]// 4, 3*vt.shape[-1]//4

vt = vt[:, nmin:nmax, nmin:nmax, nmin:nmax]  # we take the central cube of size L/2 so that opposide sides of the cloud are not correlated
print(vt.shape)

vt.T.tofile(fname_v)
vtx = vt[0, 1, :, :]
vty = vt[1, 1, :, :]
vtz = vt[2, :, :, :]
vtx.T.tofile(fname_vx)
vty.T.tofile(fname_vy)
vtz.T.tofile(fname_vz)
np.save(fname, vt)

print(nmin, nmax) 
print(vt.shape)
print(vtx.shape)
print(vtx.size)
print("===================")
# print(vtx[:,:])
aa = np.linspace(-1, 1, N//2)
x,y = np.meshgrid(np.linspace(-1, 1, N//2), np.linspace(-1, 1, N//2))
#z = x*np.exp(-x**2 - y**2)
#v, u = np.gradient(z, .2, .2)
plt.quiver(x,y,vtx,vty)

plt.show()

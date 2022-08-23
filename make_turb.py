#!/usr/bin/env python3

import numpy as np
from scipy import fftpack

def TurbField(N=128, kmin=2, maxmode = 64, sol_weight=0.5, seed=42):
    freqs = fftpack.fftfreq(N)
    freq3d = np.array(np.meshgrid(freqs,freqs,freqs,indexing='ij'))
    intfreq = np.around(freq3d*N)
    kSqr = np.sum(np.abs(freq3d)**2,axis=0)
    intkSqr = np.sum(np.abs(intfreq)**2, axis=0)
    VK = []

    # apply ~k^-2 exp(-k^2/kmax^2) filter to white noise to get x, y, and z components of velocity field
    for i in range(3):
        np.random.seed(seed+i)
        rand_phase = fftpack.fftn(np.random.normal(size=kSqr.shape)) # fourier transform of white noise
        vk = rand_phase * (float(kmin)/N)**2 / (kSqr+1e-300)
        #vk[intkSqr < kmin**2] = 0.0     # freeze out modes lower than kmin
#        print(intkSqr[intkSqr < kmin**2])
        vk[intkSqr==0] = 0.0
#        vk[intkSqr>0] *= np.exp(-kmin**2/intkSqr)
        vk[intkSqr < kmin**2] *= intkSqr[intkSqr < kmin**2]**2/kmin**4 # smoother filter than mode-freezing; should give less "ringing" artifacts
        vk *= np.exp(-intkSqr/maxmode**2)

        VK.append(vk)
    VK = np.array(VK)
    
    vk_new = np.zeros_like(VK)
    
    # do projection operator to get the correct mix of compressive and solenoidal
    for i in range(3):
        for j in range(3):
            if i==j:
                vk_new[i] += sol_weight * VK[j]
            vk_new[i] += (1 - 2 * sol_weight) * freq3d[i]*freq3d[j]/(kSqr+1e-300) * VK[j]
    vk_new[:,kSqr==0] = 0.0
    VK = vk_new
    
    vel = np.array([fftpack.ifftn(vk).real for vk in VK]) # transform back to real space
    vel -= np.average(vel,axis=(1,2,3))[:,np.newaxis,np.newaxis,np.newaxis]
    vel = vel / np.sqrt(np.sum(vel**2,axis=0).mean()) # normalize so that RMS is 1
    return np.array(vel)


kmin = 2 # lambda_max, minimum populated turbulent wavenumber for Gaussian initial velocity field, in units of pi/R [default: 2]
turb_sol = 0.5 # fraction of turbulence in solenoidal modes [default: 0.5]
seed = 55
res = 128
N = 2*res # N is 2x cube resolution. E.g. N = 128 means a 64^3 velocity cube.

print(f" Resolution = {res}^3, kmin = {kmin}, turb_sol = {turb_sol}, seed = {seed}" )
#fname = "vturb%d_sol%g_seed%d.npy"%(kmin,turb_sol, seed) # currently not needed
#fname_v = "vturb%d_sol%g_seed%d.bin"%(kmin,turb_sol, seed) # currently not needed

fname_vx = "vx.bin"
fname_vy = "vy.bin"
fname_vz = "vz.bin"

vt = TurbField(N=N, kmin=kmin, sol_weight=turb_sol, seed=seed)

nmin, nmax = vt.shape[-1]// 4, 3*vt.shape[-1]//4

vt = vt[:, nmin:nmax, nmin:nmax, nmin:nmax]  # we take the central cube of size L/2 so that opposite sides of the cloud are not correlated
print(vt.shape)

#vt.T.tofile(fname_v) # currently not needed
vtx = vt[0, :, :, :]
vty = vt[1, :, :, :]
vtz = vt[2, :, :, :]
vtx.T.tofile(fname_vx)
vty.T.tofile(fname_vy)
vtz.T.tofile(fname_vz)
#np.save(fname, vt) # currently not needed

#print(nmin, nmax) 
#print(vt.shape)
#print(vtx.shape)
#print(vtx.size)
print("===================")

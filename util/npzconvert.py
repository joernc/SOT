import numpy as np
import h5py

fl = np.load("data/argo/dtaus_Nias_H08_argo.npz")
tau = fl["dtau"].T
freq = fl["f"]
n, l = np.shape(tau)
time = np.arange(n)

fl = np.load("data/argo/dTz_Nias_H08_argo.npz")
T = fl["dT"]
z = fl["z"]
nz = np.size(z)

with h5py.File("data/argo/nias_H08.h5", "w") as f:
    ftime = f.create_dataset("time", (n,), dtype='i8')
    ftau = f.create_dataset("tau", (n, l), dtype='f8')
    ffreq = f.create_dataset("freq", (l,), dtype='f8')
    fT = f.create_dataset("T", (n, nz,), dtype='f8')
    fz = f.create_dataset("z", (nz,), dtype='f8')
    ftime[:] = time
    ftau[:,:] = tau
    ffreq[:] = freq
    fT[:,:] = T
    fz[:] = z

fl = np.load("data/kernels/KTs_Nias_H08.npz")
x = fl["x"]
z = fl["z"]
K = fl["KT"]
freq = fl["f"]
l, nx, nz = np.shape(K)

with h5py.File("data/kernels/nias_H08.h5", "w") as f:
    fx = f.create_dataset("x", (nx,), dtype='f8')
    fz = f.create_dataset("z", (nz,), dtype='f8')
    fK = f.create_dataset("K", (l, nx, nz), dtype='f8')
    ffreq = f.create_dataset("freq", (l,), dtype='f8')
    fx[:] = x
    fz[:] = z
    fK[:,:,:] = K
    ffreq[:] = freq

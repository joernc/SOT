import numpy as np
import matplotlib.pyplot as plt
import dedalus.public as de
import pandas as pd
from scipy.interpolate import interp1d
import xarray as xr
import time

# 1d extrapolation
def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return np.array(list(map(pointwise, np.array(xs))))

    return ufunclike

def cdif(x):
    return np.insert(x[2:]-x[:-2],[0,-1],[x[1]-x[0],x[-1]-x[-2]])
def grad(x,y):
    return cdif(y)/cdif(x)

model = 'ecco'
lonmin,lonmax = 55,105
latmin,latmax = -10,10

if model == 'ecco':
    ds_ref = xr.open_dataset(f'data/{model}/sigma1_2005_mean_wof.nc')
else:
    ds_ref = xr.open_dataset(f'data/{model}/sigma1_2005_mean.nc')
lon_e = ds_ref.lon.values
lat_e = ds_ref.lat.values
dlon_e = np.diff(lon_e).max()
dlat_e = np.diff(lat_e).max()
iidx = np.where((lon_e<=lonmax+dlon_e) & (lon_e>=lonmin-dlon_e))[0]
jidx = np.where((lat_e<=latmax+dlat_e) & (lat_e>=latmin-dlat_e))[0]
nx,ny = len(iidx),len(jidx)

zp = ds_ref.z.values
nz = len(zp)
print(zp[-3:])

n2 = np.full((nx,ny,nz),np.nan)
t1 = time.time()
print('Calculate stratification...',flush=True)
for i in range(nx):
    print(f'lon = {lon_e[iidx[i]]}, t0 = {(time.time()-t1)/60:.2f} min',flush=True)
    for j in range(ny):
        n2[i,j,:] = -9.8*grad(zp,ds_ref.sigma1[iidx[i],jidx[j],:].values)/1e3

ds = xr.Dataset(
    data_vars=dict(
        n2 = (["x", "y", "z"], n2),),
    coords=dict(
        lon = lon_e[iidx],
        lat = lat_e[jidx],
        z = zp,),
        attrs=dict(description=f"mean stratification from {model.upper()} 2005 daily"),
)

# grid points
n = 256

# gravity acceleration
ga = 9.8

dtype = np.complex128

nmode = 15
hnll = np.full((nmode,nx,ny,nz),np.nan)
cnll = np.full((nmode,nx,ny),np.nan)
Hll = np.full((nx,ny),np.nan)

# set up basis
zcoord = de.Coordinate('z')
dist = de.Distributor(zcoord, dtype=dtype)

t1 = time.time()
mean = True
print('Calculate modes...',flush=True)
for i in range(nx):
    for j in range(ny):
        print(f'(lon,lat) = ({lon_e[iidx[i]]},{lat_e[jidx[j]]}), t1 = {(time.time()-t1)/60:.2f} min',flush=True)
        if mean:
            n2p = np.nanmean(ds.n2.values,axis=(0,1))
        else:
            n2p = ds.n2[i,j].values.flatten()
        if len(np.argwhere(np.isnan(n2p)))>0:
            k = np.argwhere(np.isnan(n2p)).max()+1
            if k+1<nz:
                H = -zp[k]
                Hll[i,j] = H
                zbasis = de.Chebyshev(zcoord, n, bounds=(-H, 0))
                
                # generate sound slowness field
                h = dist.Field(name='h', bases=zbasis)
                tau_1 = dist.Field(name='tau_1')
                tau_2 = dist.Field(name='tau_2')
                c = dist.Field(name='c')
                n2 = dist.Field(name='n2', bases=zbasis)
                z = dist.local_grid(zbasis)
                
                # interpolated sound speed profile from data
                knots = zp
                itn2 = interp1d(zp, n2p)
                etn2 = extrap1d(itn2)
                n2['g'] = etn2(z)
                
                # Substitutions
                dz = lambda A: de.Differentiate(A, zcoord)
                lift_basis = zbasis.derivative_basis(1)
                lift = lambda A: de.Lift(A, lift_basis, -1)
                hz = dz(h) + lift(tau_1) # First-order reduction
                
                # Problem
                problem = de.EVP([h, tau_1, tau_2], eigenvalue=c, namespace=locals())
                problem.add_equation("c*n2*h + dz(hz) + lift(tau_2) = 0")
                problem.add_equation("h(z=-H) = 0")
                problem.add_equation("h(z=0) = 0")
                #problem.add_equation("c*ga*h(z=0) - dz(h)(z=0)  = 0")
                
                # Solve
                solver = problem.build_solver()
                solver.solve_dense(solver.subproblems[0])
                evals = np.sort(solver.eigenvalues)
                ievls = np.argsort(solver.eigenvalues)[evals>0][:]
                evals = evals[evals>0][:]
                
                knots = z
                for m, idx in enumerate(ievls[:nmode]):
                    solver.set_state(idx, solver.subsystems[0])
                    hn = ((h / np.sqrt(de.Integrate(n2*h**2, zcoord)/H)).evaluate()['g']).real
                    #hn = hn*np.sign(hn[np.argmax(np.abs(hn))])
                    ithn = interp1d(z, hn)
                    ethn = extrap1d(ithn)
                    hnll[m,i,j] = ethn(zp)
                    hnll[m,i,j] = hnll[m,i,j]*np.sign(hnll[m,i,j,-3])
                    hnll[m,i,j,zp<-H] = np.nan
                    cnll[m,i,j] = 1/np.abs(evals[m])**.5
                    
        if mean:
            dst = xr.Dataset(
                data_vars=dict(
                    cn=(["n"], cnll[:,i,j]),
                    hn=(["n","z"], hnll[:,i,j]), ),
                coords=dict(
                    n = 1+np.arange(nmode),
                    z = zp,),
                attrs=dict(description=f"first {nmode} dynamical modes"),)

            dst.to_netcdf(f"data/{model}/modes_mean_{model}_3d.nc")
            break
    if mean:
        break
    else:
        dst = xr.Dataset(
            data_vars=dict(
                H = (["lon","lat"], Hll),
                cn=(["n","lon","lat"], cnll),
                hn=(["n","lon", "lat","z"], hnll), ),
            coords=dict(
                n = 1+np.arange(nmode),
                lon = lon_e[iidx],
                lat = lat_e[jidx],
                z = zp,),
            attrs=dict(description=f"first {nmode} dynamical modes"),)
        dst.to_netcdf(f"data/ecco/modes_3d_{model}.nc")

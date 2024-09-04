## Import external packages and modules
import numpy as np
from scipy.interpolate import interp1d
from scipy import arange, exp
import xarray as xr
import gsw as gsw
import great_circle_calculator.great_circle_calculator as gcc
from scipy.interpolate import RegularGridInterpolator as RGI
from datetime import datetime,timedelta,date,timezone
import itertools
import warnings
warnings.filterwarnings('ignore')
import time
from scipy.signal import butter, lfilter
from obspy.signal.filter import bandpass
import h5py

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y
    
# fill in coastal and bottom points (ignores points at the edge of the array)
def filledgesxy(a):
    # fill in all points that have a missing value with the average of the eight surrounding points (if there are any values there)
    _, nx, ny = a.shape
    b = np.copy(a)
    for i,j in itertools.product(range(1,nx-1), range(1,ny-1)):
        if (np.isnan(a[0,i,j]) and any(~np.isnan(a[0,i-1:i+2,j-1:j+2]).flatten())):
            b[:,i,j] = np.nanmean(a[:,i-1:i+2,j-1:j+2],axis=(1,2))
    return b

# fill in coastal and bottom points (ignores points at the edge of the array)
def filledges(a):
    # fill in all points that have a missing value with the average of the eight surrounding points (if there are any values there)
    _, nx, ny, nz = a.shape
    b = np.copy(a)
    for i,j,k in itertools.product(range(1,nx-1), range(1,ny-1), range(nz)):
        if (np.isnan(a[0,i,j,k]) and any(~np.isnan(a[0,i-1:i+2,j-1:j+2,k]).flatten())):
            b[:,i,j,k] = np.nanmean(a[:,i-1:i+2,j-1:j+2,k],axis=(1,2))
    # fill in bottom points
    for i,j,k in itertools.product(range(nx), range(ny), range(nz-1)):
        if (np.isnan(a[0,i,j,k]) and ~np.isnan(a[0,i,j,k+1])):
            b[:,i,j,k] = a[:,i,j,k+1]
    return b
    
# in situ temperature from potential temperature
def T(pT, S, lon, lat, z):
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    CT = gsw.CT_from_pt(SA, pT)
    T = gsw.t_from_CT(SA, CT, p)
    return T

# in situ density as a function of salinity and temperature
def sigma1(pT, S, lon, lat, z):
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    CT = gsw.CT_from_pt(SA, pT)
    prho1 = gsw.sigma1(SA, CT)+1e3
    return prho1

# absolute salinity from potential temperature
def SA(S, lon, lat, z):
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    return SA

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

model = 'ecco'
var,mvar = 'T','itemp'
nmode = 5
head = '/central/groups/oceanphysics/shirui/SOT/'
f = h5py.File(f'{head}/results/anomalies/nias_tm_H08_coarsen.h5', 'r')
dsK = xr.open_dataset('~/acoustics/data/knl/coarsen/KTs_Nias_H08_coarsen.nc')
dsn2 = xr.open_dataset(f'{head}/data/{model}/n2_Nias_H08_2005.nc')
dsm = xr.open_dataset(f'{head}/results/{model}/modes_mean_{model}.nc').isel(n=slice(nmode))
dsmT = xr.open_dataset(f'{head}/results/{model}/hnd{var}_Nias_H08_0506.nc').isel(n=slice(nmode))

cmdl = dsmT.hndsig1.isel(n=slice(3)).values#np.array(f['cecco'])
t = np.array(f['ti'])/24/3600/1000
tm = .5*(t[0]+t[-1])
X = np.ones((2,len(t)))
X[1,:]=t-tm
A = np.linalg.inv(X@X.T)
dx,dz = np.diff(dsK.x).mean(),np.diff(dsK.z).mean()
N0 = 1e-3
n2G = np.nan_to_num(np.nanmean(dsn2.n2.values,axis=0)[None,:]*dsm.hn.values)
Ak = np.nan_to_num(np.nansum(dsK.SEMkernels_T.values,axis=1))@n2G.T*dx*dz/N0
print(Ak.shape)
u,s,vh = np.linalg.svd(Ak)
u *= -1
vh *= -1
# Sample rate and desired cutoff frequencies (in Hz).
fs = 1.0
lowcuts = [1/17]#,1/18,1/12,1/9.5,1/7.5,1/5.5]
highcuts = [1/10]#,1/12,1/9.5,1/7.5,1/5.5,1/4.5]
nband = len(lowcuts)
    
ds_mode =  xr.open_dataset(f'{head}results/{model}/modes_mean_{model}_3d.nc')
ds_sig1 =  xr.open_dataset(f'{head}data/{model}/{mvar}_2005_mean_wof.nc')

lonmin,lonmax = 55,105
latmin,latmax = -10,10

lon_e = ds_sig1.lon.values
lat_e = ds_sig1.lat.values
ze = ds_sig1.z.values
dlon_e = np.diff(lon_e).max()
dlat_e = np.diff(lat_e).max()
iidx = np.where((lon_e<=lonmax+dlon_e) & (lon_e>=lonmin-dlon_e))[0]
jidx = np.where((lat_e<=latmax+dlat_e) & (lat_e>=latmin-dlat_e))[0]
nx,ny,nz = len(iidx),len(jidx),len(ze)

head = '/central/groups/oceanphysics/shirui/ecco/v4r4/'
ds_temp = xr.open_mfdataset(f'{head}nctiles_daily/theta/THETA_2005*.nc')
ds_salt = xr.open_mfdataset(f'{head}nctiles_daily/salt/SALT_2005*.nc')
ds_ssh = xr.open_mfdataset(f'{head}nctiles_daily/sshdyn/SSHDYN_2005*.nc')
ds_temp,ds_salt = np.transpose(ds_temp.THETA[:,::-1,4].values,(0,3,2,1)),np.transpose(ds_salt.SALT[:,::-1,4].values,(0,3,2,1))
dssh = np.transpose(ds_ssh.SSHDYN[:,4].values,(0,2,1))
ds_temp[ds_temp==0.] = np.nan 
ds_salt[ds_salt==0.] = np.nan
dssh[dssh==0.] = np.nan
if var=='T':
    dsig1 = T(ds_temp, ds_salt,lon_e.reshape((1,-1,1,1)),lat_e.reshape((1,1,-1,1)), ze.reshape((1,1,1,-1)))-ds_sig1.itemp.values.reshape((1,len(lon_e),-1,nz))
if var=='sig1':
    dsig1 = sigma1(ds_temp, ds_salt,lon_e.reshape((1,-1,1,1)),lat_e.reshape((1,1,-1,1)), ze.reshape((1,1,1,-1)))-ds_sig1.sigma1.values.reshape((1,len(lon_e),-1,nz))
dsig1 = dsig1[:,iidx][:,:,jidx]
dssh = dssh[:,iidx][:,:,jidx]
dz = np.diff(ds_sig1.z.values)[0]
nlag = 15
cll = np.full((2,2*nlag,nx,ny,nband),np.nan)
t0 = time.time()

hndsig1 = np.nansum(dsig1*ds_mode.hn[0].values.reshape((1,1,1,nz)),axis=-1)[None,:,:,:]
for i in range(1,nmode):
    hndsig1 = np.concatenate((hndsig1,np.nansum(dsig1*ds_mode.hn[i].values.reshape((1,1,1,nz)),axis=-1)[None,:,:,:]),axis=0)
kH = np.argwhere(np.isnan(ds_mode.hn[0].values)).max()+1
H = -ds_mode.z.values[kH] 
hndsig1 *= dz/H

def ip(frac):
    return gcc.intermediate_point(p1, p2, frac)
xmaxk,xmink = max(dsn2.x),min(dsn2.x)
nxk = len(dsn2.x)
p10,p20 = (72.49, -7.65),(96.92, 1.62)
d0 = gcc.distance_between_points(p10, p20, unit='meters')
crs1,crs2 = gcc.bearing_at_p1(p10, p20),gcc.bearing_at_p1(p20, p10)
# kernel end points
p2 = gcc.point_given_start_and_bearing(p10, crs1, xmaxk, unit='meters')
p1 = gcc.point_given_start_and_bearing(p20, crs2, d0-xmink, unit='meters')
lonT1,lonT2 = p1[0],p2[0]
latT1,latT2 = p1[1],p2[1]

# Get (lon,lat) for the grid points along the path
fr = np.linspace(0,1,nxk)
gridT = np.array(list(map(ip,fr)))
gridT = gridT[(gridT[:,0]>=p10[0]) & (gridT[:,0]<=p20[0]),:]

# interpolate horizontally onto great circle path
knots = (lon_e[iidx], lat_e[jidx])

for nt in range(365):
    for nh in range(3):
        ith = RGI(knots, hndsig1[nh,nt,:,:])
        cmdl[nh,nt] = np.nanmean(ith(gridT))

# Plot the frequency response for a few different orders.
cmdlf = np.zeros((cmdl.shape[1],cmdl.shape[0],nband))
order = 5
for i in range(3):
    ki = (A@(X@cmdl[i].reshape(-1,1))).flatten()
    for j in range(nband):
        cmdlf[:,i,j] = bandpass(cmdl[i]-ki@X, lowcuts[j], highcuts[j], fs, order, True)#butter_bandpass_filter(cmdl[:,i], lowcuts[j], highcuts[j], fs, order=6)
        cmdlf[:,i,j] = cmdlf[:,i,j]/np.std(cmdlf[:,i,j])
    
hnvn = False
if hnvn:
    ds0 = xr.Dataset(
        data_vars=dict(
            h1dsig1 =  (["t","lon", "lat"], h1dsig1),
            h2dsig1 =  (["t","lon", "lat"], h2dsig1),
            vndsig1 =  (["t","m"], cmdl[:,:2]),),
        coords=dict(
            t = np.arange(365)+1,
            m = np.arange(2)+1,
            lon = lon_e[iidx],
            lat = lat_e[jidx],),
            attrs=dict(description="dynamical and acoustic mode projections"),
    )
    ds0.to_netcdf(f'{head}hnvn_{var}_wof.nc')

for i in range(nx):
    print(f'lon = {lon_e[iidx[i]]}, time = {(time.time()-t0)/60:.2f} min')
    for j in range(ny):
        if ~np.isnan(dssh[0,i,j]):
            for k in range(nband):
                x1,x2 = cmdlf[:,0,k],cmdlf[:,2,k]
                y1,y2 = hndsig1[0,:,i,j],hndsig1[2,:,i,j]#vh[0,:nmode]@hndsig1[:,:,i,j],vh[1,:nmode]@hndsig1[:,:,i,j]
                k1,k2 = (A@(X@y1[:,None])).flatten(),(A@(X@y2[:,None])).flatten()
                y1 = bandpass(y1-k1@X, lowcuts[k], highcuts[k], fs, order, True).flatten()
                y2 = bandpass(y2-k2@X, lowcuts[k], highcuts[k], fs, order, True).flatten()

                for l in range(-nlag,nlag):
                    if l<=0:
                        mtx1 = np.concatenate((np.ones(x1[-l:,None].shape),x1[-l:,None]),axis=1)
                        mtx2 = np.concatenate((np.ones(x2[-l:,None].shape),x2[-l:,None]),axis=1)
                    else:
                        mtx1 = np.concatenate((np.ones(x1[:-l,None].shape),x1[:-l,None]),axis=1)
                        mtx2 = np.concatenate((np.ones(x2[:-l,None].shape),x2[:-l,None]),axis=1)

                    if l<0:
                        try:
                            cll[0,l+nlag,i,j,k]=np.linalg.solve(mtx1.T@mtx1,mtx1.T@y1[:l,None]).flatten()[1]
                            cll[1,l+nlag,i,j,k]=np.linalg.solve(mtx2.T@mtx2,mtx2.T@y2[:l,None]).flatten()[1]
                        except:
                            print(y1)
                    else:
                        try:
                            cll[0,l+nlag,i,j,k]=np.linalg.solve(mtx1.T@mtx1,mtx1.T@y1[l:,None]).flatten()[1]
                            cll[1,l+nlag,i,j,k]=np.linalg.solve(mtx2.T@mtx2,mtx2.T@y2[l:,None]).flatten()[1]
                        except:
                            print(y1)
                        
                
    ds = xr.Dataset(
        data_vars=dict(
            c = (["i","lag","lon", "lat","fband"], cll),
            lowcuts = (["fband"], lowcuts),
            highcuts = (["fband"], highcuts),),
            #hndsig1 =  (["n","lag2","lon", "lat","fband"], hdsll),
            #c2 = (["m","lag","lon", "lat","fband"], c2ll),
            #ssh =  (["lag2","lon", "lat","fband"], sshll),
            #ccr = (["n","m","lon", "lat","fband"], ccrll),),
        coords=dict(
            #n = np.arange(nmode+1)+1,
            i = [1,3],#np.arange(2)+1,
            lag = np.arange(-nlag,nlag),
            #lag2 = np.arange(2*nlag),
            lon = lon_e[iidx],
            lat = lat_e[jidx],
            fband = np.arange(nband)+1,),
            attrs=dict(description="regressed Ti against ci in the biweekly band"),
    )
    ds.to_netcdf(f'results/ecco/regression_{var}_Omodes.nc')

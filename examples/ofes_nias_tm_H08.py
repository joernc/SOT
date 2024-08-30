## Import external packages and modules
import numpy as np
from scipy.interpolate import interp1d
from scipy import arange, exp
import xarray as xr
import gsw as gsw
import great_circle_calculator.great_circle_calculator as gcc
from scipy.interpolate import RegularGridInterpolator as RGI
from scipy.interpolate import RectBivariateSpline as RBS
from datetime import datetime,timedelta,date,timezone
import itertools
import warnings
warnings.filterwarnings('ignore')
import time
from scipy.interpolate import interp1d

# output timestep, and output choices 
# for vertical anomalies, weighted anomalies
dtself,vert,wdtau,xvert = 1,0,0,0
vdT,KdT,Kuv = 0,0,0        # get vertical, kernel-weighted, or shifted anomalies
xzr1,mdT = 0,1
lat = np.nan#3
mvar,dvar = 'itemp','T'

########
# fill in coastal and bottom points (ignores points at the edge of the array)
def filledgeswt(a):
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
    
# fill in coastal and bottom points (ignores points at the edge of the array)
def filledges(a):
    # fill in all points that have a missing value with the average of the eight surrounding points (if there are any values there)
    nx, ny, nz = a.shape
    b = np.empty([nx, ny, nz])
    b[:,:,:] = a[:,:,:]
    for i,j,k in itertools.product(range(1,nx-1), range(1,ny-1), range(nz)):
        if (np.isnan(a[i,j,k]) and any(~np.isnan(a[i-1:i+2,j-1:j+2,k]).flatten())):
            b[i,j,k] = np.nanmean(a[i-1:i+2,j-1:j+2,k])
    # fill in bottom points
    for i,j,k in itertools.product(range(nx), range(ny), range(nz-1)):
        if (np.isnan(a[i,j,k]) and ~np.isnan(a[i,j,k+1])):
            b[i,j,k] = a[i,j,k+1]
    return b

# fill in bottom points along the path
def fillbtm(a,n):
    # fill in all points that have a missing value with the average of the eight surrounding points (if there are any values there)
    nx, nz = a.shape
    b = np.empty([nx, nz])
    b[:,:] = a[:,:]
    for i,k in itertools.product(range(nx), range(nz-1)):
        if (np.isnan(a[i,k]) and ~np.isnan(a[i,k+1])):
            b[i,k] = a[i,k+1]
            if n>0:
                b[i,max([0,k-n]):k] = a[i,k+1]
    return b

# fill in bottom points along the path
def fillbtm2(a):  
  # fill in all points that have a missing value with the average of the two surrounding points (if there are any values there)
    nx, nz = a.shape
    b = np.empty([nx, nz])
    b[:,:] = a[:,:]
    for i in range(nx):
        for k in range(nz-1):
            if (np.isnan(b[i,nz-k-2]) and ~np.isnan(b[i,nz-k-1])):
                b[i,nz-k-2] = b[i,nz-k-1]
    return b
    
def ip(frac):
    return gcc.intermediate_point(p1, p2, frac)

# sound speed as a function of in situ temperature and depth
def c(T, S, lon, lat, z):
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    c = gsw.sound_speed_t_exact(SA, T, p)
    return c

# in situ temperature from potential temperature
def T(pT, S, lon, lat, z):
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    CT = gsw.CT_from_pt(SA, pT)
    T = gsw.t_from_CT(SA, CT, p)
    return T

# in situ density as a function of salinity and temperature
def rho(T, S, lon, lat, z):
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    rho = gsw.rho_t_exact(SA, T, p)
    return rho
    
# in situ density as a function of salinity and temperature
def prho0(pT, S, lon, lat, z):
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    CT = gsw.CT_from_pt(SA, pT)
    prho0 = gsw.sigma0(SA, CT)+1e3
    return prho0

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

# temperature sensitivity of sound speed
def dcdT(T, S, lon, lat, z):
    eps = 1e-4
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    c1 = gsw.sound_speed_t_exact(SA, T-eps, p)
    c2 = gsw.sound_speed_t_exact(SA, T+eps, p)
    return (c2-c1)/(2*eps)

# salinity sensitivity of sound speed
def dcdS(T, S, lon, lat, z):
    eps = 1e-4
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    c1 = gsw.sound_speed_t_exact(SA-eps, T, p)
    c2 = gsw.sound_speed_t_exact(SA+eps, T, p)
    return (c2-c1)/(2*eps)

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

src,stn = 'Nias','H08'
name = src+'_'+stn
#resolution = 'Gaussian'

# Define path end points
if stn == 'H01':
    p10,p20 = (114.14, -34.88),(96.95, 1.12)
elif stn == 'H08':
    p10,p20 = (72.49, -7.65),(96.92, 1.62)
elif stn == 'DGAR':
    p10,p20 = (72.45, -7.41),(96.96, 1.63)
d0 = gcc.distance_between_points(p10, p20, unit='meters')
crs1,crs2 = gcc.bearing_at_p1(p10, p20),gcc.bearing_at_p1(p20, p10)

ds_ref = xr.open_dataset(f'data/ofes/{mvar}_2005_mean.nc')
lon_e = ds_ref.lon.values
lat_e = ds_ref.lat.values
dlon_e = np.diff(lon_e).max()
dlat_e = np.diff(lat_e).max()
ze = ds_ref.z.values
nze = len(ze)

grid = 'coarsen'#'Gaussian'

if mdT==1 or mdT>2:
    dsm = xr.open_dataset('results/ofes/modes_full_ofes.nc')
    idx_topo = np.where(np.isnan(dsm.hn[0].values))
    
    xk,zk = dsm.x,dsm.z
    xmaxk,xmink = max(xk),min(xk)
    nxk,nzk = len(xk),len(zk)
    dxk,dzk = np.ptp(xk.values)/(nxk-1),np.ptp(zk.values)/(nzk-1)

    mean,nmode = True,3#15
    ns=np.arange(nmode)+1
    if mean:
        dsmm = xr.open_dataset('results/ofes/modes_mean_ofes.nc')
        kH = np.argwhere(np.isnan(dsmm.hn[0].values)).max()+1
        H = -dsmm.z.values[kH] 

dsk = xr.open_dataset(f'~/acoustics/data/knl/coarsen/KTs_Nias_H08_{grid}.nc')
dxk,dzk = np.diff(dsk.x)[0],np.diff(dsk.z)[0]
if KdT or xzr1 or mdT>1:
    dsk = xr.open_dataset(f'data/knl/KTs_Nias_H08_{grid}.nc')
    fs = dsk.f
    idx_topo = np.where(dsk.SEMkernels_T[0].values==0)
    
    xk,zk = dsk.x,dsk.z
    xmaxk,xmink = max(xk),min(xk)
    nxk,nzk = len(xk),len(zk)
    dxk,dzk = np.ptp(xk.values)/(nxk-1),np.ptp(zk.values)/(nzk-1)

# kernel end points
p2 = gcc.point_given_start_and_bearing(p10, crs1, xmaxk, unit='meters')
p1 = gcc.point_given_start_and_bearing(p20, crs2, d0-xmink, unit='meters')
lonT1,lonT2 = p1[0],p2[0]
latT1,latT2 = p1[1],p2[1]

# Get (lon,lat) for the grid points along the path
fr = np.linspace(0,1,nxk)
gridT = np.array(list(map(ip,fr)))
lat_k, lon_k = gridT[:,1],gridT[:,0]    
if np.isnan(lat):
    i_k = np.where((lon_e<=max(p10[0],p20[0])+dlon_e) & (lon_e>=min(p10[0],p20[0])-dlon_e))[0]
    j_k = np.where((lat_e<=max(p10[1],p20[1])+dlat_e) & (lat_e>=min(p10[1],p20[1])-dlat_e))[0]
else:
    lat_k = -lat*np.ones(lat_k.shape)
    i_k = np.where((lon_e<=np.max(lon_k)+dlon_e) & (lon_e>=np.min(lon_k)-dlon_e))[0]
    j_k = np.where((lat_e<=np.max(lat_k)+dlat_e) & (lat_e>=np.min(lat_k)-dlat_e))[0]
if dvar == 'T':
    T_ref = ds_ref.itemp.isel(lon=i_k,lat=j_k).values
else:
    T_ref = ds_ref.sigma1.isel(lon=i_k,lat=j_k).values

idx_lonk = np.where((lon_k>=lon_e[i_k].min()) & (lon_k<=lon_e[i_k].max()))[0]
idx_latk = np.where((lat_k>=lat_e[j_k].min()) & (lat_k<=lat_e[j_k].max()))[0]
idx_k = np.intersect1d(idx_lonk,idx_latk)

#########
year = '2005'
ds_temp = xr.open_dataset('data/ofes/temp_ofes_'+year+'_dy.nc')
ds_salt = xr.open_dataset('data/ofes/salinity_ofes_'+year+'_dy.nc')
dt64_events = ds_temp.time.values
t0 = time.time()

# kernel coordinates
dTk = np.full([nxk, nzk],np.nan)

if mdT:
    if mean:
        dTGn = np.full([len(ns),len(dt64_events)],np.nan)
    else:
        dtauhn = np.full([len(ns),len(dt64_events),nxk],np.nan)
    e0 = 0
    
if vert:
    dTtv = np.zeros([len(dt64_events),nzk])
    e0 = 0
    
if KdT:
    e0 = 0
    dtauSEM = np.zeros([len(fs),len(dt64_events)])
    dtauMODE = np.zeros([len(fs),len(dt64_events)])

if xvert:
    prxz = np.zeros([nxk,nzk])
    e0 = 0
    
for e, dt64 in enumerate(dt64_events[e0:]):
    if (vert or KdT or xvert or mdT)==0:
        break
    print(dt64)
    print((time.time()-t0)/60,flush=True)
    
    pTe = np.flip(ds_temp.temp.isel(lon=i_k,lat=j_k)[e].values,axis=0).T
    Se = np.flip(ds_salt.salinity.isel(lon=i_k,lat=j_k)[e].values,axis=0).T
    pTe[pTe==0.] = np.nan
    Se[Se==0.] = np.nan
    nxe, nye, nze = pTe.shape
    # fill in the coastal and bottom points
    if dvar=='T':
        Te = T(pTe, Se, lon_e[i_k,None,None], lat_e[None,j_k,None], ze[None,None,:])
    if dvar=='sig1':
        Te = sigma1(pTe, Se, lon_e[i_k,None,None], lat_e[None,j_k,None], ze[None,None,:])
    dTe = Te-T_ref
    # interpolate horizontally onto great circle path
    knots = (lon_e[i_k], lat_e[j_k])
    locsk = np.empty([len(idx_k),2])
    locsk[:,0],locsk[:,1] = lon_k[idx_k],lat_k[idx_k]

    dTkze = np.full((nxk, nze),np.nan)
    for k in range(nze):
        itpT = RGI(knots, dTe[:,:,k])
        dTkze[idx_k,k] = itpT(locsk)
    # interpolate vertically onto kernel grid
    # fill in the coastal and bottom points
    dTkze = fillbtm2(dTkze)
    knots = ze
    for i in idx_k:
        itpT = interp1d(knots, dTkze[i,:])
        etpT = extrap1d(itpT)
        dTk[i,:] = etpT(zk)

    # calculate in situ temperature
    dTk[idx_topo] = np.nan
    
    if mdT:
        if mdT==1 or mdT>2:
            #dTGn[:,e+e0] = np.nansum(dsmm.hn*np.nanmean(dTk,axis=0)[None,:],axis=-1)*dzk/H
            dTGn[:,e+e0] = np.nansum(np.nansum(dsk.SEMkernels_T,axis=1)*np.nanmean(dTk,axis=0)[None,:],axis=-1)*dxk*dzk
            dst = xr.Dataset(
                data_vars=dict(
                    hndsig1=(["n", "t"], dTGn), ),
                coords=dict(
                    n = ns,#.values,
                    t = dt64_events,),
                attrs=dict(description="Daily anomalies projected onto mean modes"),)
        else:
            dtauhn[:,e+e0] = np.nansum(dsm.hn*dTk[None,:,:],axis=-1)*dzk/dsm.H.values[None,:]/1e3
            dst = xr.Dataset(
                data_vars=dict(
                    hndsig1=(["n", "t", "x"], dtauhn), ),
                coords=dict(
                    n = ns,#.values,
                    t = dt64_events,#),
                    x = xk,),
                attrs=dict(description="Daily anomalies projected onto local modes"),)
        #dst.to_netcdf(f'results/ofes/hnd{dvar}_{name}_{year}.nc')
        dst.to_netcdf(f'results/ofes/Kjd{dvar}_{name}_{year}.nc')
        if mdT>=2:
            dtauvn[:,e+e0] = np.nansum(v2[:,None,:]*dTk[None,:,:],axis=-1)*dzk/H0
            dst = xr.Dataset(
                data_vars=dict(
                    vndsig1=(["n", "t", "x"], dtauvn), ),
                coords=dict(
                    n = ns,#.values,
                    t = dt64_events,#),
                    x = xk,),
                attrs=dict(description="Daily anomalies projected onto singular vectors"),)
            dst.to_netcdf(f'results/ofes/vnd{dvar}_{name}_{year}.nc')
        
    if xvert:
            dav = xr.DataArray(prxz,[("x", xk),("z", zk)],)
            dav.to_netcdf('results/ofes/sig1xz_'+name+'_'+year+f'_{lat}S.nc')
            
    if vert:
            dTtv[e,:] = np.nanmean(dTk,axis=0)
            dav = xr.DataArray(dTtv,[("t", dt64_events),("z", zk)],)
            dav.to_netcdf('results/ofes/dTz_'+name+'_2005.nc')
            
    if KdT:
        dtauSEM[:,e+e0] = np.nansum(dsk.SEMkernels_T*dTk[None,:,:],axis=(1,2))*dxk*dzk
            #dtauMODE[j,e+e0] = np.nansum(ds.MODEkernels_T[j]*dTk)*dxk*dzk
        dst = xr.Dataset(
            data_vars=dict(
                SEMdtaus=(["f", "t"], dtauSEM), ),
                #MODEdtaus=(["f", "t"], dtauMODE), ),
            coords=dict(
                f = fs.values,
                t = dt64_events,),
            attrs=dict(description="Daily kernel-weighted anomalies estimated from OFES"),)
        dst.to_netcdf(f'results/ofes/kTd{dvar}_{name}_{year}_{grid}.nc')
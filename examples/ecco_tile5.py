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
import requests

# output choices 
eref = 0                   # get monthly temperature and salinity on tile 5
kcrho,k2K = 0,0            # get kernel reference or convert kernel  
vdT,KdT,Kuv = 0,0,0        # get vertical, kernel-weighted, or shifted anomalies
xzr0,xzr1,mdT = 0,0,1
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

########## 
if eref:
    print(f'Reference field calculation on ECCO grid...\n')
    ## Load temperature, salinity and pressure dataset
    tile = 4
    head = '/central/groups/oceanphysics/shirui/'
    if dvar=='v':
        mfdspT = xr.open_mfdataset(f'{head}ecco/v4r4/nctiles_daily/nvel/NVEL_2005*')
    else:
        mfdspT = xr.open_mfdataset('{head}ecco/v4r4/nctiles_daily/theta/THETA_2005*')
        mfdsS = xr.open_mfdataset('{head}ecco/v4r4/nctiles_daily/salt/SALT_2005*')
        
    lon_e = mfdspT.XC.isel(tile=tile)[0]
    lat_e = mfdspT.YC.isel(tile=tile)[:,0]
    # reverse depth coordinate
    ze = np.flip(mfdspT.Z)
    te = mfdspT.time
    i_tmean = np.where((te>=np.datetime64('2005-01-01'))&(te<np.datetime64('2006-01-01')))[0]
    
    mfdspT = mfdspT.isel(time = i_tmean)
    if dvar!='v':
        mfdsS = mfdsS.isel(time = i_tmean)
    te = mfdspT.time
    
    xm,ym = 95,0
    dlon,dlat = np.diff(lon_e).min(),np.diff(lat_e)[np.argmin(np.abs(lat_e.values-ym))]
    im,jm = np.where(np.abs(lon_e.values-85)<15)[0],np.where(np.abs(lat_e.values-ym)<dlat)[0]
    print(im,jm)
    
    lon_e = lon_e[im]
    lat_e = lat_e[jm]
    if dvar=='v':
        v_tile5 = np.transpose(mfdspT.NVEL.isel(i=im,j=jm)[:,::-1,tile].values,(0,3,2,1))  
        ds_var = xr.Dataset(
            data_vars=dict(
                v = (["t","lon", "z"], np.nanmean(v_tile5,axis=2)),),
            coords=dict(
                t = te.values,
                lon = lon_e.values,
                z = ze.values,),
                attrs=dict(description="meridional velocity from ECCOv4r4 daily 2005"),
        ) 
        ds_var.to_netcdf(f'{head}SOT/data/ecco/{mvar}_year05_lat{ym:.0f}.nc')
    else:
        pT_tile5 = np.transpose(mfdspT.THETA.isel(i=im,j=jm)[:,::-1,tile].values,(0,3,2,1))
        S_tile5 = np.transpose(mfdsS.SALT.isel(i=im,j=jm)[:,::-1,tile].values,(0,3,2,1))
    T_tile5 = np.full(pT_tile5.shape,np.nan)
    wot = False
    print(wot)

    t0 = time.time()
    if ~wot:
        pT_tile5[pT_tile5==0.] = np.nan
        S_tile5[S_tile5==0.] = np.nan
        #pT_tile5,S_tile5 = filledgeswt(pT_tile5),filledgeswt(S_tile5)
        if dvar=='T':
            T_tile5 = T(pT_tile5, S_tile5, lon_e.values.reshape((1,-1,1,1)), lat_e.values.reshape((1,1,-1,1)), ze.values.reshape((1,1,1,-1)))
        if dvar=='sig1':
            T_tile5 = sigma1(pT_tile5, S_tile5, lon_e.values.reshape((1,-1,1,1)), lat_e.values.reshape((1,1,-1,1)), ze.values.reshape((1,1,1,-1)))
    else:
        for e,dt in enumerate(te):
            print(dt.values)
            print((time.time()-t0)/60)
            # mask missing data and reverse depth coordinate
            pT_tile5[e][pT_tile5[e]==0.] = np.nan
            S_tile5[e][S_tile5[e]==0.] = np.nan
            pT_tile5[e] = filledges(pT_tile5[e])
            S_tile5[e] = filledges(S_tile5[e])
            if dvar=='T':
                T_tile5[e] = T(pT_tile5[e], S_tile5[e], lon_e.values.reshape((-1,1,1)), lat_e.values.reshape((1,-1,1)), ze.values.reshape((1,1,-1)))
            if dvar=='sig1':
                T_tile5[e] = sigma1(pT_tile5[e], S_tile5[e], lon_e.values.reshape((-1,1,1)), lat_e.values.reshape((1,-1,1)), ze.values.reshape((1,1,-1)))
    
    if dvar=='T':
        ds_var = xr.Dataset(
            data_vars=dict(
                itemp = (["lon", "lat", "z"], np.nanmean(T_tile5,axis=0)),),
            coords=dict(
                #t = te.values,
                lon = lon_e.values,
                lat = lat_e.values,
                z = ze.values,),
                attrs=dict(description="In-situ temperature (K) from ECCOv4r4 daily 2005, w/o filledges"),
        )
    if dvar=='sig1':
        ds_var = xr.Dataset(
            data_vars=dict(
                sigma1 = (["lon", "z"], np.nanmean(T_tile5,axis=(0,2))),),
            coords=dict(
                #t = te.values,
                lon = lon_e.values,
                #lat = lat_e.values,
                z = ze.values,),
                attrs=dict(description="sigma1 (kg/m^3) from ECCOv4r4 daily 2005, w/o filledges"),
        )       
    dsm = xr.Dataset(
        data_vars=dict(
            itemp = (["t", "lon", "lat", "z"], T_tile5),),
        coords=dict(
            t = te.values,
            x = lon_e.values,
            y = lat_e.values,
            z = ze.values,),
            attrs=dict(description="In-situ temperature (K) from ECCOv4r4 daily 2005, w/o filledges"),
    )
    #dsm.to_netcdf(f'{head}data/ecco/{mvar}_2005_{-ym:.0f}S_wof.nc')
    ds_var.to_netcdf(f'{head}SOT/data/ecco/{mvar}_year05_lat{ym:.0f}.nc')
    
##########
src,stn = 'Nias','H08'
name = src+'_'+stn

# Define path end points
if stn == 'H01':
    p10,p20 = (114.14, -34.88),(96.95, 1.12)
elif stn == 'H08':
    p10,p20 = (72.49, -7.65),(96.92, 1.62)
elif stn == 'DGAR':
    p10,p20 = (72.45, -7.41),(96.96, 1.63)
d0 = gcc.distance_between_points(p10, p20, unit='meters')
crs1,crs2 = gcc.bearing_at_p1(p10, p20),gcc.bearing_at_p1(p20, p10)

head = '/central/groups/oceanphysics/shirui/SOT/'
ds_ref = xr.open_dataset(f'{head}data/ecco/{mvar}_2005_mean.nc')
lon_e = ds_ref.lon.values
lat_e = ds_ref.lat.values
dlon_e = np.diff(lon_e).max()
dlat_e = np.diff(lat_e).max()
ze = ds_ref.z.values
nze = len(ze)

grid = 'coarsen'

if kcrho:
    # start and end of path
    if stn != 'H01':
        xmink,xmaxk,dxk = -100e3,3000e3, 10e3 
    else:
        xmink,xmaxk,dxk = -100e3,4500e3, 5e3 

    zmink,zmaxk,dzk = -6500.,0.,50.

    # size of kernel array
    nxk,nzk = int((xmaxk-xmink)/dxk)+1,int((zmaxk-zmink)/dzk)+1

    # kernel coordinates
    xk = np.linspace(xmink,xmaxk,nxk)
    zk = np.linspace(zmink,zmaxk,nzk)
elif mdT==1 or mdT>2:
    dsm = xr.open_dataset('/central/groups/oceanphysics/shirui/SOT/results/ecco/modes_full_ecco.nc')
    idx_topo = np.where(np.isnan(dsm.hn[0].values))
    
    xk,zk = dsm.x,dsm.z
    xmaxk,xmink = max(xk),min(xk)
    nxk,nzk = len(xk),len(zk)
    dxk,dzk = np.ptp(xk.values)/(nxk-1),np.ptp(zk.values)/(nzk-1)
    
    mean,nmode = True,15
    ns=np.arange(nmode)+1
    if mean:
        dsmm = xr.open_dataset('/central/groups/oceanphysics/shirui/SOT/results/ecco/modes_mean_ecco.nc')
        kH = np.argwhere(np.isnan(dsmm.hn[0].values)).max()+1
        H = -dsmm.z.values[kH] 
    
if KdT or Kuv or xzr1 or mdT>1:
    if Kuv:
        dsk = xr.open_dataset(f'data/knl/coarsen/KUs_KTs_'+name+'_'+grid+'.nc')
    else:
        dsk = xr.open_dataset(f'data/knl/coarsen/KTs_Nias_H08_{grid}.nc')
    fs = dsk.f
    idx_topo = np.where(dsk.SEMkernels_T[0].values==0)
    
    xk,zk = dsk.x,dsk.z
    xmaxk,xmink = max(xk),min(xk)
    nxk,nzk = len(xk),len(zk)
    dxk,dzk = np.ptp(xk.values)/(nxk-1),np.ptp(zk.values)/(nzk-1)
    
    if mdT>=2:
        kT = np.nansum(dsk.SEMkernels_T,axis=1)*dxk
        H0 = 5e3
        u, s, vh = np.linalg.svd(kT, full_matrices=True)
        v2 = (H0/dzk)**0.5*vh[:2,:]
        ns=[1,2]

# kernel end points
p2 = gcc.point_given_start_and_bearing(p10, crs1, xmaxk, unit='meters')
p1 = gcc.point_given_start_and_bearing(p20, crs2, d0-xmink, unit='meters')
lonT1,lonT2 = p1[0],p2[0]
latT1,latT2 = p1[1],p2[1]

# Get (lon,lat) for the grid points along the path
fr = np.linspace(0,1,nxk)
gridT = np.array(list(map(ip,fr)))
lat_k, lon_k = gridT[:,1],gridT[:,0]
i_k = np.where((lon_e<=max(p10[0],p20[0])+dlon_e) & (lon_e>=min(p10[0],p20[0])-dlon_e))[0]
j_k = np.where((lat_e<=max(p10[1],p20[1])+dlat_e) & (lat_e>=min(p10[1],p20[1])-dlat_e))[0]
if mvar=='itemp':
    T_ref = ds_ref.itemp.isel(lon=i_k,lat=j_k).values
if mvar=='sigma1':
    T_ref = ds_ref.sigma1.isel(lon=i_k,lat=j_k).values

idx_lonk = np.where((lon_k>=lon_e[i_k].min()) & (lon_k<=lon_e[i_k].max()))[0]
idx_latk = np.where((lat_k>=lat_e[j_k].min()) & (lat_k<=lat_e[j_k].max()))[0]
idx_k = np.intersect1d(idx_lonk,idx_latk)

if xzr1:
    # interpolate horizontally onto great circle path
    knots = (lon_e[i_k], lat_e[j_k])
    locsk = np.empty([len(idx_k),2])
    locsk[:,0],locsk[:,1] = lon_k[idx_k],lat_k[idx_k]
    
    sig1kze = np.full((nxk, nze),np.nan)
    for k in range(nze):
        itpT = RGI(knots, T_ref[:,:,k])
        sig1kze[idx_k,k] = itpT(locsk)
                
    sig1k = np.zeros([nxk,nzk])
    # interpolate vertically onto kernel grid
    # fill in the coastal and bottom points
    sig1kze = fillbtm2(sig1kze)
    knots = ze
    for i in idx_k:
        itpT = interp1d(knots, sig1kze[i,:])
        etpT = extrap1d(itpT)
        sig1k[i,:] = etpT(zk)
    sig1k[idx_topo] = np.nan
    dav = xr.DataArray(sig1k,[("x", xk),("z", zk)],)
    dav.to_netcdf('results/ecco/sig1xz_Nias_H08_2005.nc')
    
#########
if kcrho:
    print(f'Reference fields calculation for {name}...\n')
    ck = np.full((nxk, nzk),np.nan)
    rhok = np.full((nxk, nzk),np.nan)
    
    # load salinity data
    ds_psmean = xr.open_dataset(f'data/ref/ecco_psalt_tile5_0515_mymean.nc')
    S_ref = ds_psmean.psalt.isel(x=i_k,y=j_k).values

    # calculate sound speed
    ce = c(T_ref, S_ref, lon_e[i_k].values.reshape((-1,1,1)), lat_e[j_k].values.reshape((1,-1,1)), ze.values.reshape((1,1,-1)))

    # calculate sound speed
    rhoe = rho(T_ref, S_ref, lon_e[i_k].values.reshape((-1,1,1)), lat_e[j_k].values.reshape((1,-1,1)), ze.values.reshape((1,1,-1)))

    # interpolate horizontally onto great circle path
    knots = (lon_e[i_k], lat_e[j_k])
    locsk = np.empty([len(idx_k),2])
    locsk[:,0],locsk[:,1] = lon_k[idx_k],lat_k[idx_k]

    ckze = np.full((nxk, nze),np.nan)
    rkze = np.full((nxk, nze),np.nan)
    for k in range(nze):
        itc = RGI(knots, ce[:,:,k])
        itr = RGI(knots, rhoe[:,:,k])
        ckze[idx_k,k] = itc(locsk)
        rkze[idx_k,k] = itr(locsk)
    # interpolate vertically onto kernel grid
    # fill in the coastal and bottom points
    ckze = fillbtm2(ckze)
    rkze = fillbtm2(rkze)
    knots = ze
    for i in range(nxk):
        itc = interp1d(knots, ckze[i,:])
        etc = extrap1d(itc)
        ck[i,:] = etc(zk)   
        itr = interp1d(knots, rkze[i,:])
        etr = extrap1d(itr)
        rhok[i,:] = etr(zk)   

    head1 = f'{xmink:.2f} {zmink:.2f} {xmaxk:.2f} {zmaxk:.2f} #x_min(m) z_min(m) x_max(m) z_max(m)\n'
    head2 = f'{dxk:.4f} {dzk:.4f}                #dx(m) dz(m)\n'
    head3 = f'{nxk} {nzk}                           #nx nz\n'

    file = open(f'{name}.txt',"w") 

    file.write(head1) 
    file.write(head2) 
    file.write(head3) 
    for i,k in itertools.product(range(nxk),range(nzk)):
        file.write(f'{xk[i]:.2f} {abs(zk[-1-k]):.2f} {ck[i,-1-k]:.4f} 0 {rhok[i,-1-k]:.4f}\n') 

    file.close()
    
###########
if k2K:
    if stn != 'H01':
        frqs = [2,3,4]
    else:
        frqs = [2.5,3.25,4]
    
    # Define path end points
    if stn == 'H01':
        p10,p20 = (114.14, -34.88),(96.95, 1.12)
    elif stn == 'H08':
        p10,p20 = (72.49, -7.65),(96.92, 1.62)
    elif stn == 'DGAR':
        p10,p20 = (72.45, -7.41),(96.96, 1.63)
    
    d0 = gcc.distance_between_points(p10, p20, unit='meters')
    crs1,crs2 = gcc.bearing_at_p1(p10, p20),gcc.bearing_at_p1(p20, p10)
    
    for i,f in enumerate(frqs):
        frq = str(f)+'Hz'
        da = xr.open_dataset(f'data/knl/all_ncFiles/{stn}_{frq}.nc')
        if i<1:
            xk,zk = np.sort(1e3*da.x),np.sort(-1e3*da.depth)
            nxk,nzk = len(xk),len(zk)
            dxk,dzk = np.ptp(xk)/(nxk-1),np.ptp(zk)/(nzk-1)
            SEMks = np.empty([len(frqs),nxk, nzk])
            SEMKTs = np.empty([len(frqs),nxk, nzk])
            MODEKTs = np.empty([len(frqs),nxk, nzk])
        SEMks[i] = np.flip(da.SEMkernels_saved[0],axis=0).T.values
        if stn == 'H01':
            da2 = xr.open_dataset(f'data/ref/ref_{name}_path2_Oct21.nc')
        else:
            da2 = xr.open_dataset(f'data/ref/ref_{name}_Oct21.nc')
        da2 = da2.interp(x=xk, z=zk)
        SEMKTs[i,:,:] = np.flip(da.SEMkernels_saved[0],axis=0).T.values/da2.c*da2.dcdT
        MODEKTs[i,:,:] = np.flip(da.MODEkernels_saved[0],axis=0).T.values/da2.c*da2.dcdT
        
    ds = xr.Dataset(
        data_vars=dict(
            SEMkernels_T=(["f", "x", "z"], SEMKTs),
            MODEkernels_T=(["f", "x", "z"], MODEKTs),
        ),
        coords=dict(
            f = frqs,
            x = xk,
            z = zk,
        ),
        attrs=dict(description="Temperature kernels in s/m^-2/K"),
    )
    #ds.to_netcdf('data/knl/KTs_'+name+'.nc')
###########
timestep = 'daily'
year = '0017'
# all events
dt64_events = np.arange('2000-01', '2018-01', dtype='datetime64[D]')
#dt64_events = np.arange('2005-01', '2016-01', dtype='datetime64[M]')
  
if vdT or KdT or Kuv or mdT:
    print(f'Anomaly calculation for {name}: vdT{vdT}, KdT{KdT}, Kuv{Kuv}...\n')

    # kernel coordinates
    if vdT or KdT or mdT:
        dTk = np.full([nxk, nzk],np.nan)
    if Kuv:
        uk = np.empty([nxk, nzk])
        vk = np.empty([nxk, nzk])
        dlatk = np.deg2rad(np.insert(lat_k[2:]-lat_k[:-2],[1,-1],[lat_k[1]-lat_k[0],lat_k[-1]-lat_k[-2]]))*6371e3
        dlonk = np.deg2rad(np.insert(lon_k[2:]-lon_k[:-2],[1,-1],[lon_k[1]-lon_k[0],lon_k[-1]-lon_k[-2]]))
        dlonk = dlonk*np.cos(np.deg2rad(lat_k))*6371e3
        dlatk = dlatk/np.insert(2*dxk*np.ones(nxk-2),[1,-1],[dxk,dxk])
        dlonk = dlonk/np.insert(2*dxk*np.ones(nxk-2),[1,-1],[dxk,dxk])

    if mdT:
        if mean:
            dTGn = np.full([len(ns),len(dt64_events)],np.nan)
        else:
            dtauhn = np.full([len(ns),len(dt64_events),nxk],np.nan)
            dtauvn = np.full([len(ns),len(dt64_events),nxk],np.nan)
        e0 = 0

    if vdT:
        try:
            dTtv = xr.open_dataarray('result/ecco/dTz_'+name+'_ofes.nc').values
            e0 = np.nonzero(dTtv[:,0])[0].max()+1
        except FileNotFoundError:
            dTtv = np.zeros([len(dt64_events),nzk])
            e0 = 0
    if KdT:
        try:
            dst = xr.open_dataset('result/ecco/kTdT_'+name+'_'+year+'.nc')
            dtauSEM = dst.SEMdtaus.values
            #dtauMODE = dst.MODEdtaus.values
            e0 = np.nonzero(dtauSEM[0])[0].max()+1
        except FileNotFoundError:
            dtauSEM = np.zeros([len(fs),len(dt64_events)])
            #dtauMODE = np.zeros([len(fs),len(dt64_events)])
            e0=0
    if Kuv:
        dtauSEM = np.zeros([len(fs),len(dt64_events)])
        e0=0
    
    t0 = time.time()
    for e, dt64 in enumerate(dt64_events[e0:]):

        # find time stamps of ECCO data bracketing the event
        ts = (dt64 - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
        event = datetime.utcfromtimestamp(ts)
        sec_since_epoch = (event.date() - date(1970, 1, 1)).total_seconds()
        day0 = datetime.utcfromtimestamp(sec_since_epoch)
        if (timestep == 'daily') or (timestep == 'monthly'):
            day1 = day0 - timedelta(hours=12)
        else:
            if event.hour<12:
                day1 = day0 - timedelta(hours=12)
            else:
                day1 = day0 + timedelta(hours=12)

        day2 = day1 + timedelta(1)
        print(day2,flush=True)
        print((time.time()-t0)/60,flush=True)

        if day2.year>=2018:
            break

        head = '/central/groups/oceanphysics/shirui/ecco/v4r4/'
        webhead = 'https://data.nas.nasa.gov/ecco/download_data.php?file=/eccodata/llc_90/ECCOv4/Release4/'
        if vdT or KdT or xzr0 or mdT:
            # load ECCO data
            if timestep == 'monthly':
                dspT2 = xr.open_dataset(f'{head}nctiles_monthly/theta/THETA_{day2.year:4}_{day2.month:02}.nc')
                dsS2 = xr.open_dataset(f'{head}nctiles_monthly/salt/SALT_{day2.year:4}_{day2.month:02}.nc')
            else: 
                dspT2 = None
                while dspT2 is None:
                    try:
                        dspT2 = xr.open_dataset(f'{head}nctiles_daily/theta/THETA_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc')
                    except:
                        url = f'https://data.nas.nasa.gov/ecco/download_data.php?file=/eccodata/llc_90/ECCOv4/Release4/nctiles_daily/THETA/{day2.year:4}/{day2.timetuple().tm_yday:03}/THETA_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc'
                        print('Downloading data!')
                        print(url)
                        r = requests.get(url, allow_redirects=True)
                        open(f'{head}nctiles_daily/theta/THETA_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc', 'wb').write(r.content)
             
                dsS2 = None
                while dsS2 is None:
                    try:
                        dsS2 = xr.open_dataset(f'{head}nctiles_daily/salt/SALT_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc')
                    except:
                        url = f'https://data.nas.nasa.gov/ecco/download_data.php?file=/eccodata/llc_90/ECCOv4/Release4/nctiles_daily/SALT/{day2.year:4}/{day2.timetuple().tm_yday:03}/SALT_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc'
                        print('Downloading data!')
                        print(url)
                        r = requests.get(url, allow_redirects=True)
                        open(f'{head}nctiles_daily/salt/SALT_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc', 'wb').write(r.content)   
    
            # mask missing data and reverse depth coordinate
            pTe = np.flip(dspT2.THETA.isel(i=i_k,j=j_k)[0,:,4].values,axis=0).T
            pTe[pTe==0.] = np.nan
            Se = np.flip(dsS2.SALT.isel(i=i_k,j=j_k)[0,:,4].values,axis=0).T
            Se[Se==0.] = np.nan    
            pTe = filledges(pTe)
            Se = filledges(Se)
    
            # fill in the coastal and bottom points
            if dvar=='T':
                Te = T(pTe, Se, lon_e[i_k].reshape((-1,1,1)), lat_e[j_k].reshape((1,-1,1)), ze.reshape((1,1,-1)))
            if dvar=='sig1':
                Te = sigma1(pTe, Se, lon_e[i_k].reshape((-1,1,1)), lat_e[j_k].reshape((1,-1,1)), ze.reshape((1,1,-1)))
            dTe = Te-T_ref
    
            # interpolate horizontally onto great circle path
            knots = (lon_e[i_k], lat_e[j_k])
            locsk = np.empty([len(idx_k),2])
            locsk[:,0],locsk[:,1] = lon_k[idx_k],lat_k[idx_k]
    
            dTkze = np.zeros((nxk, nze))
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
                    if mean:
                        dTGn[:,e+e0] = np.nansum(dsmm.hn*np.nanmean(dTk,axis=0)[None,:],axis=-1)*dzk/H
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
                    dst.to_netcdf(f'results/ecco/hnd{dvar}_{name}_{year}.nc')
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
                    dst.to_netcdf(f'results/ecco/vnd{dvar}_{name}_{year}.nc')
            if vdT:
                dTtv[e+e0,:] = np.nanmean(dTk,axis=0)
                da = xr.DataArray(dTtv,[("time", dt64_events),("z",zk.values)],)
                da.to_netcdf('results/ecco/dTz_'+name+'_ecco_2005.nc')
            if KdT:
                for j in range(len(fs)):
                    dtauSEM[j,e+e0] = np.nansum(dsk.SEMkernels_T[j]*dTk)*dxk*dzk
                    #dtauMODE[j,e+e0] = np.nansum(ds.MODEkernels_T[j]*dTk)*dxk*dzk
                dst = xr.Dataset(
                    data_vars=dict(
                        SEMdtaus=(["f", "t"], dtauSEM), ),
                        #MODEdtaus=(["f", "t"], dtauMODE), ),
                    coords=dict(
                        f = fs.values,
                        t = dt64_events,),
                    attrs=dict(description="Daily kernel-weighted anomalies estimated from ECCOv4r4"),)
                dst.to_netcdf(f'results/ecco/kTd{dvar}_{name}_{year}.nc')
        else:
            dsU = None
            while dsU is None:
                try:
                    dsU = xr.open_dataset(f'{head}nctiles_daily/evel/EVEL_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc')
                except:
                    url = f'https://data.nas.nasa.gov/ecco/download_data.php?file=/eccodata/llc_90/ECCOv4/Release4/nctiles_daily/EVEL/{day2.year:4}/{day2.timetuple().tm_yday:03}/EVEL_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc'
                    print('Downloading data!')
                    print(url)
                    r = requests.get(url, allow_redirects=True)
                    open(f'{head}nctiles_daily/evel/EVEL_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc', 'wb').write(r.content)
         
            dsV = None
            while dsV is None:
                try:
                    dsV = xr.open_dataset(f'{head}nctiles_daily/nvel/NVEL_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc')
                except:
                    url = f'https://data.nas.nasa.gov/ecco/download_data.php?file=/eccodata/llc_90/ECCOv4/Release4/nctiles_daily/NVEL/{day2.year:4}/{day2.timetuple().tm_yday:03}/NVEL_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc'
                    print('Downloading data!')
                    print(url)
                    r = requests.get(url, allow_redirects=True)
                    open(f'{head}nctiles_daily/nvel/NVEL_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc', 'wb').write(r.content)  
                    
            # mask missing data and reverse depth coordinate
            ue = np.flip(dsU.EVEL.isel(i=i_k,j=j_k)[0,:,4].values,axis=0).T
            ue[ue==0.] = np.nan
            ve = np.flip(dsV.NVEL.isel(i=i_k,j=j_k)[0,:,4].values,axis=0).T
            ve[ve==0.] = np.nan    
            ue = filledges(ue)
            ve = filledges(ve)
            
            # interpolate horizontally onto great circle path
            knots = (lon_e[i_k], lat_e[j_k])
            locsk = np.empty([len(idx_k),2])
            locsk[:,0],locsk[:,1] = lon_k[idx_k],lat_k[idx_k]
    
            ukze = np.full((nxk, nze),np.nan)
            vkze = np.full((nxk, nze),np.nan)
            for k in range(nze):
                itu = RGI(knots, ue[:,:,k])
                ukze[idx_k,k] = itu(locsk)
                itv = RGI(knots, ve[:,:,k])
                vkze[idx_k,k] = itv(locsk)
            # interpolate vertically onto kernel grid
            # fill in the coastal and bottom points
            ukze = fillbtm2(ukze)
            vkze = fillbtm2(vkze)
            knots = ze
            for i in range(nxk):
                itu = interp1d(knots, ukze[i,:])
                etu = extrap1d(itu)
                uk[i,:] = etu(zk) 
                itv = interp1d(knots, vkze[i,:])
                etv = extrap1d(itv)
                vk[i,:] = etv(zk) 
    
            uk[idx_topo] = np.nan
            vk[idx_topo] = np.nan
            velk = -(uk*dlonk.reshape((-1,1))+vk*dlatk.reshape((-1,1)))
            for j in range(len(fs)):
                    dtauSEM[j,e+e0] = np.nansum(ds.SEMkernels_U[j]*velk)*dxk*dzk
            dst = xr.Dataset(
                data_vars=dict(
                    SEMdtaus=(["f", "t"], dtauSEM), ),
                coords=dict(
                    f = fs.values,
                    t = dt64_events,),
                attrs=dict(description="Daily Doppler shift (s) estimated from ECCOv4r4"),)
            dst.to_netcdf('results/ecco/doppler_'+name+'_ecco_'+grid+'KTs.nc')
## Import external packages and modules
import numpy as np
from scipy.interpolate import interp1d
from scipy import arange, exp
import xarray as xr
import h5py
import gsw as gsw 
import great_circle_calculator.great_circle_calculator as gcc
from scipy.interpolate import RegularGridInterpolator as RGI
# from scipy.interpolate import RectBivariateSpline as RBS
from datetime import datetime,timedelta,date,timezone
import itertools
import warnings
warnings.filterwarnings('ignore')
import time
import requests

# output choices 
eref = 1                   # get monthly temperature and salinity on tile 5
kcrho,k2K = 0,0            # get kernel reference or convert kernel  
vdT,KdT,Kuv = 0,0,0        # get vertical, kernel-weighted, or shifted anomalies
xzr,mdT,xtau = 0,0,0
mvar,dvar = 'itemp','T'
model = 'ecco'
lat = -3
########
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
def sigma0(pT, S, lon, lat, z):
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    CT = gsw.CT_from_pt(SA, pT)
    prho0 = gsw.sigma0(SA, CT)
    return prho0

# in situ density as a function of salinity and temperature
def sigma1(pT, S, lon, lat, z):
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    CT = gsw.CT_from_pt(SA, pT)
    prho1 = gsw.sigma1(SA, CT)
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
    print(f'Reference {mvar} field calculation on ECCO grid...\n')
    ## Load temperature, salinity and pressure dataset
    mfdspT = xr.open_mfdataset('/central/groups/oceanphysics/shirui/ecco/v4r4/nctiles_daily/theta/THETA_2005*')
    mfdsS = xr.open_mfdataset('/central/groups/oceanphysics/shirui/ecco/v4r4/nctiles_daily/salt/SALT_2005*')
    
    lon_e2 = mfdspT.XC[1][0]
    lon_e5 = mfdspT.XC[4][0]
    lat_e = mfdspT.YC[4][:,0].values
    dlat = np.diff(lat_e).max()
    # reverse depth coordinate
    ze = np.flip(mfdspT.Z).values
    te = mfdspT.time
    i_tmean = np.where((te>=np.datetime64('2005-01-01'))&(te<np.datetime64('2006-01-01')))[0]
    mfdspT = mfdspT.isel(time = i_tmean)
    mfdsS = mfdsS.isel(time = i_tmean)
    te = mfdspT.time.values

    zbe = np.zeros(1)
    for z in ze[::-1]:
        zbe = np.insert(zbe,0,2*z-zbe[0])
    dz = zbe[1:]-zbe[:-1]
    print(dz)

    # mask missing data and reverse depth coordinate
    laidx = np.argmin(np.abs(lat_e-lat))
    #np.where((lat_e>=lat-dlat)&(lat_e<=lat+dlat))[0]
    lidx2 = np.where((lon_e2>=35)&(lon_e2<=105))[0]
    lidx5 = np.where((lon_e5>=35)&(lon_e5<=105))[0]
    lon_e = np.append(lon_e2[lidx2],lon_e5[lidx5])
    print(lon_e)
    print(lat_e[laidx])
    
    pT_tile2 = np.transpose(mfdspT.isel(tile=1,i=lidx2,j=laidx).THETA[:,::-1].values,(0,2,1))
    pT_tile5 = np.transpose(mfdspT.isel(tile=4,i=lidx5,j=laidx).THETA[:,::-1].values,(0,2,1))
    pT_e = np.concatenate((pT_tile2,pT_tile5),axis=1)
    pT_e[pT_e==0.] = np.nan
    S_tile2 = np.transpose(mfdsS.isel(tile=1,i=lidx2,j=laidx).SALT[:,::-1].values,(0,2,1))
    S_tile5 = np.transpose(mfdsS.isel(tile=4,i=lidx5,j=laidx).SALT[:,::-1].values,(0,2,1))
    S_e = np.concatenate((S_tile2,S_tile5),axis=1)
    S_e[S_e==0.] = np.nan
    
    if dvar=='T':
        ds_var = T(pT_e, S_e, lon_e.reshape((1,-1,1)), lat_e[laidx], ze.reshape((1,1,-1)))
        
        ds_var = xr.Dataset(
            data_vars=dict(
                dT = (["time","lon", "z"], ds_var-np.nanmean(ds_var,axis=0)),),
            coords=dict(
                time = te,
                lon = lon_e,
                z = ze,),
                attrs=dict(description="In-situ temperature anomaly (K) from ECCOv4r4 daily 2005, w/o filledges"),
        )

        nmode = 15
        ns=np.arange(nmode)+1
        dsmm = xr.open_dataset('/central/groups/oceanphysics/shirui/SOT/results/ecco/modes_mean_ecco_-3lat.nc')
        kH = np.argwhere(np.isnan(dsmm.hn[0].values)).max()+1
        H = -dsmm.z.values[kH] 
        print(f'Project onto modes...\n') 
        for i in range(nmode):
            Gi = np.nan_to_num(dsmm.hn[i].values)
            if i==0:
                hndT = np.nansum(ds_var.dT.values*Gi[None,None,:]*dz[None,None,:],axis=-1)[None,:,:]
            else:
                hndT = np.concatenate((hndT,np.nansum(ds_var.dT.values*Gi[None,None,:]*dz[None,None,:],axis=-1)[None,:,:]),axis=0)
        hndT = hndT/H
        
        dst = xr.Dataset(
            data_vars=dict(
                hndT=(["n", "t", "x"], hndT), ),
            coords=dict(
                n = np.arange(nmode)+1,
                t = te,
                x = lon_e,),
            attrs=dict(description="Daily in-situ temperature anomalies (K) projected onto local modes"),)
        dst.to_netcdf(f"results/ecco/hndT_{model}_{-lat:.0f}S.nc")
    else:
        ds_var = sigma0(pT_e, S_e, lon_e.reshape((1,-1,1)), lat_e[laidx], ze.reshape((1,1,-1)))
        
        ds_var = xr.Dataset(
            data_vars=dict(
                sigma0 = (["lon", "z"], np.nanmean(ds_var,axis=0)),),
            coords=dict(
                lon = lon_e,
                z = ze,),
                attrs=dict(description="mean sigma0 (kg/m^3) from ECCOv4r4 daily 2005, w/o filledges"),
        )
        head = '/central/groups/oceanphysics/shirui/SOT/'
        ds_var.to_netcdf(f'{head}data/ecco/{mvar}_05_nias{-lat:.0f}S.nc')
        
#############
ds_itmean = xr.open_dataset(f'~/acoustics/data/ref/ecco_sigma1_tile2a5_0506_mymean.nc')
lon_e = ds_itmean.x
lat_e = ds_itmean.y
dlon_e = np.diff(lon_e).max()
dlat_e = np.diff(lat_e).max()
ze = ds_itmean.z
nze = len(ze)

grid = 'coarsen'
lat = 3
if mdT==1:
    ds = xr.open_dataset(f'/central/groups/oceanphysics/shirui/ofes/results/modes_full_ecco_{lat}S.nc')
    idx_topo = np.where(np.isnan(ds.hn[0].values))
    
    xk,zk = ds.x,ds.z
    xmaxk,xmink = max(xk),min(xk)
    nxk,nzk = len(xk),len(zk)
    dxk,dzk = np.ptp(xk.values)/(nxk-1),np.ptp(zk.values)/(nzk-1)
    
    mean,nmode = True,15
    ns=np.arange(nmode)+1
    if mean:
        dsmm = xr.open_dataset('/central/groups/oceanphysics/shirui/SOT/results/ecco/modes_mean_ecco_{lat}S.nc')
        kH = np.argwhere(np.isnan(dsmm.hn[0].values)).max()+1
        H = -dsmm.z.values[kH] 
else:
    ds = xr.open_dataset('data/knl/coarsen/KTs_Nias_H08_coarsen.nc')
    zk = ds.z
    nzk = len(zk)
    dzk = np.ptp(zk.values)/(nzk-1)
    
    xmink,xmaxk,dxk = 0,6100e3, 20e3 

    # size of kernel array
    nxk,nzk = int((xmaxk-xmink)/dxk)+1,len(zk)
    xk = np.linspace(xmink,xmaxk,nxk)
    if mdT==2:
        kT = np.nansum(ds.SEMkernels_T,axis=1)*dxk
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
lat_k = -lat*np.ones(lat_k.shape)
i_k = np.where((lon_e<=np.max(lon_k)+dlon_e) & (lon_e>=np.min(lon_k)-dlon_e))[0]
j_k = np.where((lat_e<=np.max(lat_k)+dlat_e) & (lat_e>=np.min(lat_k)-dlat_e))[0]
T_ref = ds_itmean.sigma1.isel(x=i_k,y=j_k).values

idx_lonk = np.where((lon_k>=lon_e[i_k].min().values) & (lon_k<=lon_e[i_k].max().values))[0]
idx_latk = np.where((lat_k>=lat_e[j_k].min().values) & (lat_k<=lat_e[j_k].max().values))[0]
idx_k = np.intersect1d(idx_lonk,idx_latk)

#########
timestep = 'daily'
year = '2005'
# all events
dt64_events = np.arange('2005-01', '2006-01', dtype='datetime64[D]')
#dt64_events = np.arange('2005-01', '2016-01', dtype='datetime64[M]')
  
if vdT or KdT or Kuv or xzr or mdT or xtau:
    print(f'Anomaly calculation for {name}: vdT{vdT}, KdT{KdT}, Kuv{Kuv}, mdT{mdT}...\n')
    print(f'Absolute calculation for {name}: xzr{xzr}, xtau{xtau}...\n')

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
        dtauhn = np.full([len(ns),len(dt64_events),nxk],np.nan)
        e0 = 0
    if xzr:
        prxz = np.zeros([nxk,nzk])
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
            dst = xr.open_dataset('result/ecco/dtaus_'+name+'_ecco_'+resolution+'KTs.nc')
            dtauSEM = dst.SEMdtaus.values
            dtauMODE = dst.MODEdtaus.values
            e0 = np.nonzero(dtauSEM[0])[0].max()+1
        except FileNotFoundError:
            dtauSEM = np.zeros([len(fs),len(dt64_events)])
            dtauMODE = np.zeros([len(fs),len(dt64_events)])
            e0=0
    if Kuv:
        dtauSEM = np.zeros([len(fs),len(dt64_events)])
        e0=0
    if xtau:
        tauxk = np.full([len(dt64_events),nxk],np.nan)
        tauyk = np.full([len(dt64_events),nxk],np.nan)
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
        if vdT or KdT or xzr or mdT:
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
    
            if e==0:
                lon_e2 = dspT2.XC[1][0].values
                lon_e5 = dspT2.XC[4][0].values
                dlon_e2,dlon_e5 = np.diff(lon_e2).max(),np.diff(lon_e5).max()
                i_k2 = np.where((lon_e2<=np.max(lon_k)+dlon_e2) & (lon_e2>=np.min(lon_k)-dlon_e2))[0]
                i_k5 = np.where((lon_e5<=np.max(lon_k)+dlon_e5) & (lon_e5>=np.min(lon_k)-dlon_e5))[0]
            # mask missing data and reverse depth coordinate
            pTe2 = np.flip(dspT2.THETA.isel(time=0,tile=1,i=i_k2,j=j_k).values,axis=0).T
            pTe5 = np.flip(dspT2.THETA.isel(time=0,tile=4,i=i_k5,j=j_k).values,axis=0).T
            pTe = np.concatenate((pTe2,pTe5),axis=0)
            pTe[pTe==0.] = np.nan
            Se2 = np.flip(dsS2.SALT.isel(time=0,tile=1,i=i_k2,j=j_k).values,axis=0).T
            Se5 = np.flip(dsS2.SALT.isel(time=0,tile=4,i=i_k5,j=j_k).values,axis=0).T
            Se = np.concatenate((Se2,Se5),axis=0)
            Se[Se==0.] = np.nan    
            pTe = filledges(pTe)
            Se = filledges(Se)
    
            # fill in the coastal and bottom points
            #pre = prho0(pTe, Se, lon_e[i_k].values.reshape((-1,1,1)), lat_e[j_k].values.reshape((1,-1,1)), ze.values.reshape((1,1,-1)))
            Te = sigma1(pTe, Se, lon_e[i_k].values.reshape((-1,1,1)), lat_e[j_k].values.reshape((1,-1,1)), ze.values.reshape((1,1,-1)))
            dTe = Te-(1-xzr)*T_ref
    
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
            #dTkze = fillbtm2(dTkze)
            knots = ze
            for i in idx_k:
                itpT = interp1d(knots, dTkze[i,:])
                etpT = extrap1d(itpT)
                if xzr:
                    prxz[i,:] = (e*prxz[i,:]+etpT(zk))/(e+1) 
                else:
                    dTk[i,:] = etpT(zk)        
    
            # calculate in situ temperature
            if xzr:
                #prxz[idx_topo] = np.nan
                dav = xr.DataArray(prxz,[("x", xk),("z", zk)],)
                dav.to_netcdf('results/ecco/sig1xz_'+name+'_'+year+f'_{lat}S.nc')
            #else:
            #    dTk[idx_topo] = np.nan
                
            if mdT:
                for j in range(len(ns)):
                    if mdT==1:
                        dtauhn[j,e+e0] = np.nansum(ds.hn[j]*dTk,axis=1)*dzk/ds.H
                    else:
                        dtauhn[j,e+e0] = np.nansum(v2[j].reshape((1,-1))*dTk,axis=1)*dzk/H0
                if mdT==1:
                    dst = xr.Dataset(
                        data_vars=dict(
                            hndsig1=(["n", "t", "x"], dtauhn), ),
                        coords=dict(
                            n = ns,#.values,
                            t = dt64_events,#),
                            x = xk,),
                        attrs=dict(description="Daily potential density anomalies (s) projected onto local modes"),)
                    dst.to_netcdf(f'results/ecco/hndsig1_{name}_{year}_{lat}S.nc')
                else:
                    dst = xr.Dataset(
                        data_vars=dict(
                            vndsig1=(["n", "t", "x"], dtauhn), ),
                        coords=dict(
                            n = ns,#.values,
                            t = dt64_events,#),
                            x = xk,),
                        attrs=dict(description="Daily potential density anomalies (s) projected onto singular vectors"),)
                    dst.to_netcdf(f'results/ecco/vndsig1_{name}_{year}_{lat}S.nc')

        if xtau:
            # load ECCO data
            dstaux = None
            while dstaux is None:
                try:
                    dstaux = xr.open_dataset(f'{head}nctiles_daily/taux/TAUX_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc')
                except:
                    url = f'https://data.nas.nasa.gov/ecco/download_data.php?file=/eccodata/llc_90/ECCOv4/Release4/nctiles_daily/EXFtaux/{day2.year:4}/{day2.timetuple().tm_yday:03}/EXFtaux_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc'
                    print('Downloading data!')
                    print(url)
                    r = requests.get(url, allow_redirects=True)
                    open(f'{head}nctiles_daily/taux/TAUX_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc', 'wb').write(r.content)
         
            dstauy = None
            while dstauy is None:
                try:
                    dstauy = xr.open_dataset(f'{head}nctiles_daily/tauy/TAUY_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc')
                except:
                    url = f'https://data.nas.nasa.gov/ecco/download_data.php?file=/eccodata/llc_90/ECCOv4/Release4/nctiles_daily/EXFtauy/{day2.year:4}/{day2.timetuple().tm_yday:03}/EXFtauy_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc'
                    print('Downloading data!')
                    print(url)
                    r = requests.get(url, allow_redirects=True)
                    open(f'{head}nctiles_daily/tauy/TAUY_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc', 'wb').write(r.content)  
    
            if e==0:
                lon_e2 = dstaux.XC[1][0].values
                lon_e5 = dstaux.XC[4][0].values
                dlon_e2,dlon_e5 = np.diff(lon_e2).max(),np.diff(lon_e5).max()
                i_k2 = np.where((lon_e2<=np.max(lon_k)+dlon_e2) & (lon_e2>=np.min(lon_k)-dlon_e2))[0]
                i_k5 = np.where((lon_e5<=np.max(lon_k)+dlon_e5) & (lon_e5>=np.min(lon_k)-dlon_e5))[0]
            # mask missing data and reverse depth coordinate
            tauxe2 = np.flip(dstaux.EXFtaux.isel(time=0,tile=1,i=i_k2,j=j_k).values,axis=0).T
            tauxe5 = np.flip(dstaux.EXFtaux.isel(time=0,tile=4,i=i_k5,j=j_k).values,axis=0).T
            tauxe = np.concatenate((tauxe2,tauxe5),axis=0)
            tauxe[tauxe==0.] = np.nan
            tauye2 = np.flip(dstauy.EXFtauy.isel(time=0,tile=1,i=i_k2,j=j_k).values,axis=0).T
            tauye5 = np.flip(dstauy.EXFtauy.isel(time=0,tile=4,i=i_k5,j=j_k).values,axis=0).T
            tauye = np.concatenate((tauye2,tauye5),axis=0)
            tauye[tauye==0.] = np.nan  
            tauxe = filledgesxy(tauxe)
            tauye = filledgesxy(tauye)
            
            # interpolate horizontally onto great circle path
            knots = (lon_e[i_k], lat_e[j_k])
            locsk = np.empty([len(idx_k),2])
            locsk[:,0],locsk[:,1] = lon_k[idx_k],lat_k[idx_k]
    
            itptx = RGI(knots, tauxe)
            tauxk[e+e0,idx_k] = itptx(locsk)
            itpty = RGI(knots, tauye)
            tauyk[e+e0,idx_k] = itpty(locsk)
            
            dst = xr.Dataset(
                data_vars=dict(
                    taux=(["t", "x"], tauxk), tauy=(["t", "x"], tauyk), ),
                coords=dict(
                    t = dt64_events,#),
                    x = xk,),
                attrs=dict(description="surface wind stress"),)
            dst.to_netcdf('results/ecco/tauxy_'+name+'_'+year+f'_{lat}S.nc')
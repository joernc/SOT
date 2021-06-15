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

name = 'H08S2'
frqs = ['2Hz','2Hz_attenuation','3Hz','4Hz']
#frqs = ['2.5Hz','3.25Hz','4Hz']

# Define path end points
#p10,p20 = (114.1361, -34.8832),(96.8, 1.)
p10,p20 = (72.49, -7.65),(97.0, 1.65)
#p10,p20 = (166.90986,19.71786),(142.11, 39) # (142.11, 39) (141,37.65)
d0 = gcc.distance_between_points(p10, p20, unit='meters')
crs1,crs2 = gcc.bearing_at_p1(p10, p20),gcc.bearing_at_p1(p20, p10)

# start and end of path
xmink,xmaxk = -100e3,3000e3
zmink,zmaxk = -6500.,0.

# kernel grid spacing
dxk,dzk = 10e3,100.

# size of kernel array
nxk,nzk = int((xmaxk-xmink)/dxk)+1,int((zmaxk-zmink)/dzk)+1

# kernel coordinates
xk = np.linspace(xmink,xmaxk,nxk)
zk = np.linspace(zmink,zmaxk,nzk)

# output timestep, and output choices 
# for vertical anomalies, weighted anomalies
dtself,vert,wdtau,refself = 1,0,1,0

# kernel end points
p2 = gcc.point_given_start_and_bearing(p10, crs1, xmaxk, unit='meters')
p1 = gcc.point_given_start_and_bearing(p20, crs2, d0-xmink, unit='meters')
lonT1,lonT2 = p1[0],p2[0]
latT1,latT2 = p1[1],p2[1]

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

def ip(frac):
    return gcc.intermediate_point(p1, p2, frac)

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

# sound speed as a function of in situ temperature and depth
def c(T, S, lon, lat, z):
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    c = gsw.sound_speed_t_exact(SA, T, p)
    return c

# temperature sensitivity of sound speed
def dcdT(T, S, lon, lat, z):
    eps = 1e-4
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    c1 = gsw.sound_speed_t_exact(SA, T-eps, p)
    c2 = gsw.sound_speed_t_exact(SA, T+eps, p)
    return (c2-c1)/(2*eps)

def dt_midmonth(year,month):
    dt_start = datetime(year,month,1)
    if month<12:
        dt_end = datetime(year,month+1,1)
    else:
        dt_end = datetime(year+1,1,1)
    dt_mid = dt_start+(dt_end-dt_start)/2
    return dt_mid

# Get (lon,lat) for the grid points along the path
fr = np.linspace(0,1,nxk)
gridxk = np.deg2rad(np.array(list(map(ip,fr))))
# Get (lon,lat) for the grid points along the path
lat_k, lon_k = gridxk[:,1],gridxk[:,0]
# mean pot. temperature and salinity
data_bar = xr.open_dataset(f'data/ref/ref_{name}_614.nc')
T_bar = data_bar.T
c_bar = data_bar.c
dcdT = data_bar.dcdT

KTk = np.empty([len(frqs),nxk, nzk])

if wdtau:
    if name[0]=="j":
        flnm = 'data/knl/kernel2D.txt'
        a = np.loadtxt(flnm, skiprows=7)
        # read kernel and mask missing values
        Kc = 1e-7*a[:,2].reshape((nzk, nxk)).T
        Kc[Kc==0.] = np.nan
    else:    
        for i,f in enumerate(frqs):
            flnm = 'data/knl/kernel_'+name+'_'+f+'.txt'
            a = np.loadtxt(flnm, skiprows=4)
            # read kernel and mask missing values
            Kc = 1e-7*a[:,2].reshape((nzk, nxk)).T
            Kc[Kc==0.] = np.nan
            # temperature kernel
            KTk[i] = Kc/c_bar*dcdT

head = 'data/argo/'
t_argo = xr.open_mfdataset(f'{head}RG_ArgoClim_Temperature_2019.nc',decode_times=False,combine='by_coords')
s_argo = xr.open_mfdataset(f'{head}RG_ArgoClim_Salinity_2019.nc',decode_times=False,combine='by_coords')
lon_sorti = np.append(np.array(range(160,360)),np.array(range(0,160)))
lon_a = np.heaviside(t_argo.LONGITUDE-180,0)*(t_argo.LONGITUDE-360)+np.heaviside(180-t_argo.LONGITUDE,1)*t_argo.LONGITUDE

t_a = t_argo.ARGO_TEMPERATURE_MEAN+t_argo.ARGO_TEMPERATURE_ANOMALY
s_a = s_argo.ARGO_SALINITY_MEAN+s_argo.ARGO_SALINITY_ANOMALY

lonMin = min(lonT1,lonT2)
lonMax = max(lonT1,lonT2)
latMin = min(latT1,latT2)
latMax = max(latT2,latT2)
lonT_r = (lon_a>=lonMin-3.)&(lon_a<=lonMax+3.)
latT_r = (t_argo.LATITUDE>=latMin-3.)&(t_argo.LATITUDE<=latMax+3.)
# convert coordinates to radiuas
lon_a = np.deg2rad(lon_a[lonT_r])
lat_a = np.deg2rad(t_argo.LATITUDE.values[latT_r])
# reverse depth coordinate
za = -np.flip(t_argo.PRESSURE)

t_inda = t_a[:,latT_r,lonT_r]
s_inda = s_a[:,latT_r,lonT_r]
T_a = np.flip(t_inda.values,axis=0).T
T_a[T_a==0.] = np.nan
S_a = np.flip(s_inda.values,axis=0).T
S_a[S_a==0.] = np.nan
# all events
if dtself:
    dt64_events = t_argo.TIME
else:
    dt64_events = np.load('dt64_H01W3_101520.npy')
    dt64_events = dt64_events[dt64_events>np.datetime64('2004-01-01T00:00:00Z')]

Tk = np.empty([len(dt64_events), nxk, nzk])
Sk = np.empty([len(dt64_events), nxk, nzk])
if vert:
    dTtv = np.zeros([len(dt64_events),nzk])
    
if wdtau:
    dtaut = np.zeros([len(frqs),len(dt64_events)])
    
for e, dt64 in enumerate(dt64_events):
    if dtself:
        print(f'{e} of {len(dt64_events)}')
        Ta = T_a[e]
        Sa = S_a[e]
    else:
        # find time stamps of ECCO data bracketing the event
        ts = (dt64 - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
        event = datetime.utcfromtimestamp(ts)
        dt_mid = dt_midmonth(event.year,event.month)
        nt_a = 12*(event.year-2004)+event.month-1
        if event<dt_mid:
            nt_day1 = nt_a-1
            if event.month>1:
                day1 = dt_midmonth(event.year,event.month-1)
            else:
                day1 = dt_midmonth(event.year-1,12)
            day2 = dt_mid
        else:
            day1 = dt_mid
            nt_day1 = nt_a
            if event.month<12:
                day2 = dt_midmonth(event.year,event.month+1)
            else:
                day2 = dt_midmonth(event.year+1,1)

        nt_day2 = nt_day1+1

        print(day2)

        # mask missing data and reverse depth coordinate
        T1,T2 = T_a[nt_day1],T_a[nt_day2]
       # S1,S2 = S_a[nt_day1],S_a[nt_day2]

        # interpolate to the time of the event
        ddt = day2-day1
        w1 = (day2 - event)/ddt
        w2 = (event - day1)/ddt
        Ta = w1*T1 + w2*T2
        print(f'{w1} {w2}')
        
    # Sa = w1*S1 + w2*S2
    # ECCO data array size
    nxa, nya, nza = Ta.shape
    # fill in the coastal and bottom points
    Ta = filledges(Ta)
    Sa = filledges(Sa)
    # Get (lon,lat) for the grid points along the path
    # interpolate horizontally onto great circle path
    knots = (lon_a, lat_a)
    locsk = np.empty([nxk,2])
    locsk[:,0],locsk[:,1] = lon_k,lat_k
    Tkza = np.empty([nxk, nza])
    Skza = np.empty([nxk, nza])
    for k in range(nza):
        itpT = RGI(knots, Ta[:,:,k])
        itpS = RGI(knots, Sa[:,:,k])
        Tkza[:,k] = itpT(locsk)
        Skza[:,k] = itpS(locsk)
    # interpolate vertically onto kernel grid
    # fill in the coastal and bottom points
    Tkza = fillbtm(Tkza,2)
    Skza = fillbtm(Skza,2)
    knots = za
    for i in range(nxk):
        itpT = interp1d(knots, Tkza[i,:])
        itpS = interp1d(knots, Skza[i,:])
        etpT = extrap1d(itpT)
        etpS = extrap1d(itpS)
        Tk[e,i,:] = etpT(zk)
        Tk[e,i,zk<-2e3] = np.nan
        Sk[e,i,:] = etpS(zk)       
        Sk[e,i,zk<-2e3] = np.nan
    if refself != 1:
        # calculate in situ temperature
        dT = Tk[e]-T_bar
        if vert:
            dTtv[e,:] = np.nanmean(dT,axis=0)
            np.savez('result/ecco_argo/dTz_'+name+'_argo.npz',dT = dTtv,time=dt64_events,z=zk)
        if wdtau:
            for i in range(len(frqs)):
                dtaut[i,e] = np.nansum(KTk[i]*dT)*dxk*dzk
            da = xr.DataArray(dtaut,[("freq", frqs),("time", dt64_events)],)
            da.to_netcdf('result/ecco_argo/dtaus_'+name+'_argo.nc')
       
if refself == 1:
    Tk_bar = np.nanmean(Tk,axis=0)
    Sk_bar = np.nanmean(Sk,axis=0)
    dT = Tk-Tk_bar.reshape((1,nxk,nzk))
    # calculate sound speed
    c_bar = c(Tk_bar, Sk_bar, np.rad2deg(lon_k).reshape((-1,1)), np.rad2deg(lat_k).reshape((-1,1)), zk.reshape((1,-1)))
    dcdTk = dcdT(Tk_bar, Sk_bar, np.rad2deg(lon_k).reshape((-1,1)), np.rad2deg(lat_k).reshape((-1,1)), zk.reshape((1,-1)))
    KTk[0,:,:] = Kc/c_bar*dcdTk
    dtaut[0,:] = np.nansum(KTk[0,:,:]*dT,axis=(1,2))*dxk*dzk
    np.savez('result/ecco_argo/dtaus_'+name+'_argo.npz',dtau = dtaut,time=dt64_events,f=frqs)
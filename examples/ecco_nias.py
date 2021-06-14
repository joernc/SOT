## Import external packages and modules
import numpy as np
from scipy.interpolate import interp1d
from scipy import arange, array, exp
import xarray as xr
import gsw as gsw
import great_circle_calculator.great_circle_calculator as gcc
from scipy.interpolate import RegularGridInterpolator as RGI
from datetime import datetime,timedelta,date,timezone
import itertools
import requests
import warnings
warnings.filterwarnings('ignore')

name = 'Nias_H01'
frqs = ['2.5Hz','3.25Hz','4Hz']#['2Hz','2Hz_attenuation','3Hz','4Hz']

# Define path end points
p10,p20 = (114.1361, -34.8832),(96.8, 1.)
#p10,p20 = (72.49, -7.65),(97.0, 1.65)
#p10,p20 = (166.90986,19.71786),(142.11, 39) # (142.11, 39) (141,37.65)
d0 = gcc.distance_between_points(p10, p20, unit='meters')
crs1,crs2 = gcc.bearing_at_p1(p10, p20),gcc.bearing_at_p1(p20, p10)

# start and end of path
xmink,xmaxk = -100e3,4500e3
zmink,zmaxk = -6500.,0.

# kernel grid spacing
dxk,dzk = 5e3,50.

# kernel end points
p2 = gcc.point_given_start_and_bearing(p10, crs1, xmaxk, unit='meters')
p1 = gcc.point_given_start_and_bearing(p20, crs2, d0-xmink, unit='meters')
#p1,p2 = (140,40),(168,18)

# size of kernel array
nxk,nzk = int((xmaxk-xmink)/dxk)+1,int((zmaxk-zmink)/dzk)+1

# kernel coordinates
xk = np.linspace(xmink,xmaxk,nxk)
zk = np.linspace(zmink,zmaxk,nzk)

ref,crho = 0,0
# output timestep, and output choices 
# for 1d vertical anomalies, scalar weighted anomalies and temperature kernels
timestep = 'daily'
vert,wdtau,ktemp,wvel = 0,0,1,0
    
########
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

# sound speed as a function of in situ temperature and depth
def c(pT, S, lon, lat, z):
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    CT = gsw.CT_from_pt(SA, pT)
    c = gsw.sound_speed(SA, CT, p)
    return c

# in situ temperature from potential temperature
def T(pT, S, lon, lat, z):
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    CT = gsw.CT_from_pt(SA, pT)
    T = gsw.t_from_CT(SA, CT, p)
    return T

# in situ density as a function of salinity and temperature
def rho(pT, S, lon, lat, z):
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    CT = gsw.CT_from_pt(SA, pT)
    rho = gsw.rho(SA, CT, p)
    return rho

# absolute salinity from potential temperature
def SA(S, lon, lat, z):
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    return SA

# temperature sensitivity of sound speed
def dcdT(pT, S, lon, lat, z):
    eps = 1e-4
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    CT = gsw.CT_from_pt(SA, pT)
    T = gsw.t_from_CT(SA, CT, p)
    CT1 = gsw.CT_from_t(SA,T-eps,p)
    CT2 = gsw.CT_from_t(SA,T+eps,p)
    c1 = gsw.sound_speed(SA, CT1, p)
    c2 = gsw.sound_speed(SA, CT2, p)
    return (c2-c1)/(2*eps)

# salinity sensitivity of sound speed
def dcdS(pT, S, lon, lat, z):
    eps = 1e-4
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(S, p, lon, lat)
    CT1 = gsw.CT_from_pt(SA-eps, pT)
    CT2 = gsw.CT_from_pt(SA+eps, pT)
    c1 = gsw.sound_speed(SA-eps, CT1, p)
    c2 = gsw.sound_speed(SA+eps, CT2, p)
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

#########
if ref or crho:
    head = 'D:\\Caltech\\Acoustics\\ecco\\v4r4\\'
    # load ECCO data
    dspT = xr.open_mfdataset(f'{head}inputFig2_Shirui\\THETA_2006_10_20.nc')

    # convert coordinates to radiuas
    lon_e = np.deg2rad(dspT.XC[4][0])
    lat_e = np.deg2rad(dspT.YC[4][:,0])
    # reverse depth coordinate
    ze = np.flip(dspT.Z)

    data = np.load(f'{head}subset\\fle.npz')
    pTe = data['pTe']
    Se = data['Se']
    # ECCO data array size
    nte, nxe, nye, nze = pTe.shape

    # Get (lon,lat) for the grid points along the path
    fr = np.linspace(0,1,nxk)
    gridT = np.deg2rad(np.array(list(map(ip,fr))))
    lat_k, lon_k = gridT[:,1],gridT[:,0]

    # interpolate horizontally onto great circle path
    locsk = np.empty([nxk,2])
    locsk[:,0],locsk[:,1] = lon_k,lat_k
    pTkze = np.empty([nte, nxk, nze])
    Skze = np.empty([nte, nxk, nze])
    pTk = np.empty([nte, nxk, nzk])
    Sk = np.empty([nte, nxk, nzk])

    rngl = range(12,144)

    for l in rngl:
        print(f'tl = {l}')
        knots = (lon_e, lat_e)
        pTel,Sel = filledges(pTe[l]),filledges(Se[l])
        for k in range(nze):
            itppT = RGI(knots, pTel[:,:,k])
            itpS = RGI(knots, Sel[:,:,k])
            pTkze[l,:,k] = itppT(locsk)
            Skze[l,:,k] = itpS(locsk)
        pTkze[l,:,:] = fillbtm(pTkze[l,:,:],2)
        Skze[l,:,:] = fillbtm(Skze[l,:,:],2)

        # interpolate vertically onto kernel grid
        knots = ze
        for i in range(nxk):
            itppT = interp1d(knots, pTkze[l,i,:])
            itpS = interp1d(knots, Skze[l,i,:])
            etppT = extrap1d(itppT)
            etpS = extrap1d(itpS)
            pTk[l,i,:] = etppT(zk)
            Sk[l,i,:] = etpS(zk)

    pT_bar,S_bar = pTk[12:144].mean(axis=0),Sk[12:144].mean(axis=0)

    # calculate sound speed
    ck = c(pT_bar, S_bar, np.rad2deg(lon_k).reshape((-1,1)), np.rad2deg(lat_k).reshape((-1,1)), zk.reshape((1,-1)))

    if ref:
        # calculate in-situ temperature
        Tk = T(pT_bar, S_bar, np.rad2deg(lon_k).reshape((-1,1)), np.rad2deg(lat_k).reshape((-1,1)), zk.reshape((1,-1)))
        # calculate dcdT
        dcdTk = dcdT(pT_bar, S_bar, np.rad2deg(lon_k).reshape((-1,1)), np.rad2deg(lat_k).reshape((-1,1)), zk.reshape((1,-1)))
        ds = xr.Dataset(data_vars=dict(c=(["x", "z"], ck),T=(["x", "z"], Tk),dcdT=(["x", "z"], dcdTk),),
                        coords=dict(x=xk,z=zk,),)   
        ds.to_netcdf(f'result/ref_{name}_Jun21.nc')         

    if crho:
        # calculate sound speed
        rhok = rho(pT_bar, S_bar, np.rad2deg(lon_k).reshape((-1,1)), np.rad2deg(lat_k).reshape((-1,1)), zk.reshape((1,-1)))

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
else:
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
    # all events
    if timestep == 'daily':
        dt64_events = np.arange('2000-01', '2019-12', dtype='datetime64[D]')
    elif timestep == 'monthly':
        dt64_events = np.arange('2000-01', '2019-12', dtype='datetime64[M]')
    else:
        dt64_flnm = 'dt64_H01W3_101520.npy'
        # all events
        dt64_events = np.load(dt64_flnm)

    # kernel coordinates
    xk = np.linspace(xmink,xmaxk,nxk)
    zk = np.linspace(zmink,zmaxk,nzk)
    pTk = np.empty([nxk, nzk])
    Sk = np.empty([nxk, nzk])
    KTk = np.empty([len(frqs),nxk, nzk])
    if wvel:
        uk = np.empty([nxk, nzk])
        wut = np.zeros([len(frqs),len(dt64_events)])
    if vert:
        dTtv = np.zeros([len(dt64_events),nzk])
    if wdtau:
        dtaut = np.zeros([len(frqs),len(dt64_events)])

    for i,f in enumerate(frqs):
        flnm = 'data/knl/kernel_'+name+'_'+f+'.txt'
        a = np.loadtxt(flnm, skiprows=4)
        # read kernel and mask missing values
        Kc = 1e-7*a[:,2].reshape((nzk, nxk)).T
        Kc[Kc==0.] = np.nan
        # temperature kernel
        KTk[i,:,:] = Kc/c_bar*dcdT
    if ktemp:
        da = xr.DataArray(KTk,[("freq", frqs),("x",xk),("z",zk)],)
        da.to_netcdf('data/knl/KTs_'+name+'.nc')

    for e, dt64 in enumerate(dt64_events):
        if (vert and wdtau and wvel)==0:
            break

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
        print(day2)

        if day2.year>=2018:
            break

        head = '/central/groups/oceanphysics/shirui/ecco/v4r4/'
        webhead = 'https://data.nas.nasa.gov/ecco/download_data.php?file=/eccodata/llc_90/ECCOv4/Release4/'
        # load ECCO data
        if timestep == 'monthly':
            u = None
            fileloc = f'{head}nctiles_monthly/EVEL/EVEL_{day2.year:4}_{day2.month:02}.nc'
            while u is None:
                try: u = xr.open_dataset(fileloc)
                except:
                    url = f'{webhead}nctiles_monthly/EVEL/{day2.year:4}/EVEL_{day2.year:4}_{day2.month:02}.nc'
                    print('Downloading data!')
                    print(url)
                    r = requests.get(url, allow_redirects=True)
                    open(fileloc,'wb').write(r.content)
        else: 
            if timestep != 'daily':
                dspT1 = None
                while dspT1 is None:
                    try:
                        dspT1 = xr.open_mfdataset(f'{head}theta/daily/*THETA_{day1.year:4}_{day1.month:02}_{day1.day:02}.nc',combine='by_coords')
                    except:
                        url = f'https://data.nas.nasa.gov/ecco/download_data.php?file=/eccodata/llc_90/ECCOv4/Release4/nctiles_daily/THETA/{day1.year:4}/{day1.timetuple().tm_yday:03}/THETA_{day1.year:4}_{day1.month:02}_{day1.day:02}.nc'
                        print('Downloading data!')
                        print(url)
                        r = requests.get(url, allow_redirects=True)
                        open(f'{head}theta/daily/download_THETA_{day1.year:4}_{day1.month:02}_{day1.day:02}.nc', 'wb').write(r.content)

            dspT2 = None
            while dspT2 is None:
                try:
                    dspT2 = xr.open_mfdataset(f'{head}theta/daily/*THETA_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc',combine='by_coords')
                except:
                    url = f'https://data.nas.nasa.gov/ecco/download_data.php?file=/eccodata/llc_90/ECCOv4/Release4/nctiles_daily/THETA/{day2.year:4}/{day2.timetuple().tm_yday:03}/THETA_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc'
                    print('Downloading data!')
                    print(url)
                    r = requests.get(url, allow_redirects=True)
                    open(f'{head}theta/daily/download_THETA_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc', 'wb').write(r.content)

            if timestep != 'daily':
                dsS1 = None
                while dsS1 is None:
                    try:
                        dsS1 = xr.open_mfdataset(f'{head}salt/daily/*SALT_{day1.year:4}_{day1.month:02}_{day1.day:02}.nc',combine='by_coords')
                    except:
                        url = f'https://data.nas.nasa.gov/ecco/download_data.php?file=/eccodata/llc_90/ECCOv4/Release4/nctiles_daily/SALT/{day1.year:4}/{day1.timetuple().tm_yday:03}/SALT_{day1.year:4}_{day1.month:02}_{day1.day:02}.nc'
                        print('Downloading data!')
                        print(url)
                        r = requests.get(url, allow_redirects=True)
                        open(f'{head}salt/daily/download_SALT_{day1.year:4}_{day1.month:02}_{day1.day:02}.nc', 'wb').write(r.content)

            dsS2 = None
            while dsS2 is None:
                try:
                    dsS2 = xr.open_mfdataset(f'{head}salt/daily/*SALT_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc',combine='by_coords')
                except:
                    url = f'https://data.nas.nasa.gov/ecco/download_data.php?file=/eccodata/llc_90/ECCOv4/Release4/nctiles_daily/SALT/{day2.year:4}/{day2.timetuple().tm_yday:03}/SALT_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc'
                    print('Downloading data!')
                    print(url)
                    r = requests.get(url, allow_redirects=True)
                    open(f'{head}salt/daily/download_SALT_{day2.year:4}_{day2.month:02}_{day2.day:02}.nc', 'wb').write(r.content)   


        # convert coordinates to radiuas
        lon_e = np.deg2rad(dspT2.XC[4][0])
        lat_e = np.deg2rad(dspT2.YC[4][:,0])
        # reverse depth coordinate
        ze = np.flip(dsS2.Z)
        nze = len(ze)
        # mask missing data and reverse depth coordinate
        if timestep == 'monthly':
            ue = np.swapaxes(np.flip(u.EVEL[0,:,4].values,axis=0),0,-1)
            ue[ue==0.]=np.nan
            ue = filledges(ue)
        else:
            if timestep != 'daily':
                pT1 = np.swapaxes(np.flip(dspT1.THETA[0,:,4].values,axis=0),0,-1)
                pT1[pT1==0.] = np.nan
                S1 = np.swapaxes(np.flip(dsS1.SALT[0,:,4].values,axis=0),0,-1) 
                S1[S1==0.] = np.nan

            pT2 = np.swapaxes(np.flip(dspT2.THETA[0,:,4].values,axis=0),0,-1)
            pT2[pT2==0.] = np.nan
            S2 = np.swapaxes(np.flip(dsS2.SALT[0,:,4].values,axis=0),0,-1) 
            S2[S2==0.] = np.nan    

            # interpolate to the time of the event
            if timestep == 'daily':
                pTe,Se = 1*pT2,1*S2
            else: 
                w1 = (day2 - event)/timedelta(1)
                w2 = (event - day1)/timedelta(1)
                pTe = w1*pT1 + w2*pT2
                Se = w1*S1 + w2*S2

            # ECCO data array size
            nxe, nye, nze = pTe.shape
            # fill in the coastal and bottom points
            pTe = filledges(pTe)
            Se = filledges(Se)       

        # interpolate horizontally onto great circle path
        knots = (lon_e, lat_e)
        locsk = np.empty([nxk,2])
        locsk[:,0],locsk[:,1] = lon_k,lat_k
        if timestep == 'monthly':
            ukze = np.empty([nxk, nze])
        else:
            pTkze = np.empty([nxk, nze])
            Skze = np.empty([nxk, nze])
        for k in range(nze):
            if monthly:
                itpu = RGI(knots, ue[:,:,k])
                ukze[:,k] = itpu(locsk)
            else:
                itppT = RGI(knots, pTe[:,:,k])
                itpS = RGI(knots, Se[:,:,k])
                pTkze[:,k] = itppT(locsk)
                Skze[:,k] = itpS(locsk)
        # interpolate vertically onto kernel grid
        # fill in the coastal and bottom points
        if timestep == 'monthly':
            ukze = fillbtm(ukze,2)
        else:
            pTkze = fillbtm(pTkze,2)
            Skze = fillbtm(Skze,2)
        knots = ze
        for i in range(nxk):
            if monthly:
                itpu = interp1d(knots, ukze[i,:])
                etpu = extrap1d(itpu)
                uk[i,:] = etpu(zk)
            else:
                itppT = interp1d(knots, pTkze[i,:])
                itpS = interp1d(knots, Skze[i,:])
                etppT = extrap1d(itppT)
                etpS = extrap1d(itpS)
                pTk[i,:] = etppT(zk)
                Sk[i,:] = etpS(zk)             

        if timestep == 'monthly':
            for i in range(len(frqs)):
                intK = np.nansum(KTk[i])*dxk*dzk
                wut[i,e] = np.nansum(KTk[i]*uk)*dxk*dzk/intK
            da = xr.DataArray(wut,[("freq", frqs),("time", dt64_events)],)
            da.to_netcdf('result/wgtU_'+name+'.nc')
        else:
            # calculate in situ temperature
            Tk = T(pTk, Sk, lon_k.reshape((-1,1)), lat_k.reshape((-1,1)), zk.reshape((1,-1)))
            dT = Tk-Tk_bar
            if vert:
                dTtv[e,:] = np.nanmean(dT,axis=0)
                da = xr.DataArray(dTtv,[("time", dt64_events),("z",zk)],)
                da.to_netcdf('result/dTz_'+name+'.nc')
            if wdtau:
                for i in range(len(frqs)):
                    dtaut[i,e] = np.nansum(KTk[i]*dT)*dxk*dzk
                da = xr.DataArray(dtaut,[("frq",frqs,"time", dt64_events)],)
                da.to_netcdf('result/dtaus_'+name+'.nc')

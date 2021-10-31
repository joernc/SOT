# TODO:
# - Check out-of-memory errors
# - Add caculation for repeating pairs 

### Gauss-Markov mapping for Argo float data
using PyCall, DelimitedFiles, Distances, CSV, DataFrames, Dates,LinearAlgebra, NCDatasets, Printf, Statistics, Interpolations, PyPlot, ProgressBars,HDF5 #,BenchmarkTools
np = pyimport("numpy")
gsw = pyimport("gsw")
gcc = pyimport("great_circle_calculator.great_circle_calculator")
cmo = pyimport("cmocean.cm")
apy = pyimport("argopy")
argo_loader = apy.DataFetcher()
#smp = pyimport("sklearn.metrics.pairwise") 

df = DateFormat("y-m-dTH:M:S.s") 
const earthradius = 6371e3
const p10 = ((166.90986,19.71786) .+ (166.89109,19.71356) .+ (166.89677,19.73115))./3

# start-end of path
name = "jpnth_H11"

# frequency/ies 
frq = ["2.5Hz"] #,"3.25Hz","4Hz"]
#frq = ["2Hz","3Hz","4Hz"]

# start and end of path
xmink = 0
xmaxk = 3200e3
zmink = -4000.
zminkK = -11e3
zmaxk = 0.

# kernel grid spacing
Δxk = 20e3
Δzk = 100.

# kernel end points
#p20 = (96.8, 1.) 
p20 = (142.11, 39) #(97.0, 1.65) (142.11,39) (140.57,37)
d0 = gcc.distance_between_points(p10, p20, unit="meters")
crs1 = gcc.bearing_at_p1(p10, p20)
crs2 = gcc.bearing_at_p1(p20, p10)
p2 = gcc.point_given_start_and_bearing(p10, crs1, xmaxk, unit="meters")
p1 = gcc.point_given_start_and_bearing(p20, crs2, d0-xmink, unit="meters")
λ1 = p1[1]
θ1 = p1[2]
λ2 = p2[1]
θ2 = p2[2]

function getTauErr(filename,fmttime,KT,K1d,Δxk,Δzk,region_params,dttm,Δday,ll_path,dist_max,mp_lvl,xyz_grds,lxy,lz,lt,omg,sign,ecco_ref)
    if !isfile(filename)
        grd_inv,P1 = getInvxP(region_params,dttm,Δday,ll_path,dist_max,mp_lvl,xyz_grds,lxy,lz,lt,omg,sign,ecco_ref)
        # write to file
        h5open(filename, "w") do fid
          # start time in microseconds since 1970-01-01T00:00:00
          fid["time"] = fmttime
          
          # station latitude
          fid["xinv"] = grd_inv
        
          # station longitude
          fid["P"] = P1
        end
    else
        # read data
        fid = h5open(filename, "r")
        grd_inv = read(fid, "xinv")
        P1 = read(fid, "P")
        close(fid)
    end
    offset = sum(skipmissing(KT[end:-1:1,end+1-nzk:end].*grd_inv))*Δxk*Δzk
    sig_dtau = (Δxk*Δzk)*(K1d' * P1 * K1d)^0.5
    return offset,sig_dtau
end

function getInvxP(region_params,dttm,Δday,ll_path,dist_max,mp_lvl,xyz_grds,lxy,lz,lt,omg,sign,ecco_ref)
    rxx1,rnn1,E1,y1,nflts1 = inversion(region_params,dttm,Δday,ll_path,dist_max,mp_lvl,xyz_grds,lxy,lz,lt,omg,sign,ecco_ref)
    @printf("Calculating 1st inversion...\n")
    grd_inv = reshape(getInversion(E1, rxx1, rnn1, y1)[1:length(z_grds)],:,nzk)
    @printf("Calculating 1st matrix P...\n")
    P1 = getP(E1, rxx1, rnn1)[1:length(z_grds),1:length(z_grds)]
    return grd_inv,P1
end

function getRyy(E::AbstractMatrix, rxx::AbstractMatrix, rnn)
    return rnn + E * rxx * E'
end

function getP(E::AbstractMatrix, rxx::AbstractMatrix, rnn)
    P = rxx - rxx * E' * inv(E * rxx * E' + rnn) * E * rxx
    return P
end

function getInversion(
        E::AbstractMatrix, 
        rxx::AbstractMatrix, 
        rnn, 
        y::AbstractArray,
    )
	return rxx * E' * inv(E * rxx * E' + rnn) * y
end

function getSigx(z)
    s_argo = DataFrame(CSV.File("result/sigx_"*name*".csv"))
    sigx = zeros(size(z))
    for (k,z_argo) in enumerate(s_argo.z)
        idx = findall(x->z_argo-1.5 < x ≤ z_argo+1.5,z)
        sigx[idx] .= s_argo.sigx[k]
    end
    kmin,kmax = argmin(s_argo.z),argmax(s_argo.z)
    idx = findall(x->x ≤ s_argo.z[kmin],z)
    sigx[idx] .= s_argo.sigx[kmin]
    idx = findall(x->s_argo.z[kmax] ≤ x,z)
    sigx[idx] .= s_argo.sigx[kmax]
    return sigx
end

# Calculate great-circle distance
function dist(lon1, lat1, lon2, lat2)
  Δσ = acos(min(1,sind(lat1)*sind(lat2) + cosd(lat1)*cosd(lat2)*cosd(lon1-lon2)))
  return earthradius*Δσ
end

function dist_flt2pth(ll_path,ll_flt)
    global d = dist(140, 18, 168, 40)
    for i = 1:length(ll_path[:,1])
        di = dist(ll_path[i,1], ll_path[i,2], ll_flt[1], ll_flt[2])
        if di < d
            global d = di
        else
            break
        end
    end
    return d
end

function getSTdist(xyz_all,t_all)
    n_all = length(t_all)
    dist_xy = pairwise(Haversine(1e-3earthradius), xyz_all[:,1:2]', dims=2)
    dist_z = pairwise(Euclidean(), xyz_all[:,3]', dims=2)
    lag_t = (repeat(t_all,1,n_all).-repeat(permutedims(t_all),n_all,1))
    #@printf("size and max lag_t: %s, %s\n",size(lag_t),maximum(lag_t))
    lag_t = 1e-3Dates.value.(lag_t)/(24*3600)
    return dist_xy,dist_z,lag_t
end

function inversion(region_params,dttm,Δday,ll_path,dist_max,mp_lvl,xyz_grds,lxy,lz,lt,omg,sign,ecco_ref;season=false,plot=false)
    dt1,dt2 = dttm - Day(Δday), dttm + Day(Δday)
    tstart = @sprintf("%4d-%02d-%02d",year(dt1),month(dt1),day(dt1))
    tend = @sprintf("%4d-%02d-%02d",year(dt2),month(dt2),day(dt2))
    append!(region_params,[tstart,tend])
    
    @printf("From %s to %s\n",tstart,tend)
    @printf("Loading data...\n")
    ds = nothing
    while isnothing(ds)
        try
            # connect
            ds = argo_loader.region(region_params).to_xarray()
        catch y
            continue
        end
    end
    @printf("Download finished!\n")
    
    ll_flts = hcat(Array{Float64,1}(ds["LONGITUDE"].values),Array{Float64,1}(ds["LATITUDE"].values))
    ll_eflts = unique(ll_flts,dims=1)
    dist_f2p = Array{Float64,1}(undef, length(ll_eflts[:,1]))
    for i = 1:length(ll_eflts[:,1])
        dist_f2p[i] = dist_flt2pth(ll_path,ll_eflts[i,:])/1e3
    end
    idx_near = findall(x -> x < dist_max,dist_f2p)
    @printf("Number of near-by floats: %d\n",length(idx_near))
    
    for idx in idx_near
        idx_lon,idx_lat = findall(x -> x==ll_eflts[idx,1],ll_flts[:,1]),findall(x -> x==ll_eflts[idx,2],ll_flts[:,2])
        if idx == idx_near[1]
            global iall_nrflts = intersect(idx_lat,idx_lon)
        else
            global iall_nrflts = vcat(iall_nrflts,intersect(idx_lat,idx_lon))
        end
    end

    xyzT_flts = hcat(ll_flts[iall_nrflts,:],-Array{Float64,1}(ds["PRES"].values)[iall_nrflts],Array{Float64,1}(ds["TEMP"].values)[iall_nrflts])
    t_flts = (Nanosecond.(ds["TIME"].values.astype(np.int64)) + DateTime(1970))[iall_nrflts]
    
    @printf("Mean field interpolation...\n")
    λe = Array{Union{Missing, Float64}}(ecco_ref["lon"]); replace!(λe,NaN=>missing)
    θe = Array{Union{Missing, Float64}}(ecco_ref["lat"]); replace!(θe,NaN=>missing)
    ze = Array{Union{Missing, Float64}}(ecco_ref["z"]); replace!(ze,NaN=>missing)
    Te = Array{Union{Missing, Float64}}(ecco_ref["__xarray_dataarray_variable__"]); replace!(Te,NaN=>missing)
    Te = permutedims(Te, [3, 2, 1])
    knots = (λe, θe, ze)
    itpT = interpolate(knots, Te, Gridded(Constant()))
    etpT = extrapolate(itpT, Flat())
    yref = etpT.(xyzT_flts[:,1],xyzT_flts[:,2],xyzT_flts[:,3])
    @printf("Mean field interpolation finished!\n")
    
    idxy = findall(x->!ismissing(x),yref)
    n_grds,n_flts,xyzT_flts,t_flts = length(xyz_grds[:,1]),length(idxy),xyzT_flts[idxy,:],t_flts[idxy]
    n_all = n_grds+n_flts
    xyz_all = vcat(xyz_grds,xyzT_flts[:,1:3])
    t_all =  vcat(repeat([dttm],n_grds),t_flts) 
    
    dist_xy,dist_z,lag_t = getSTdist(xyz_all,t_all)
    
    @printf("Spatial and temporal distances finished!\n")
    
    if season
        A_t = exp.(-abs.(lag_t)/lt) + .5*cos.(omg*lag_t) + .5*cos.(2*omg*lag_t)
        A_t = .5*A_t
    else
        A_t = exp.(-abs.(lag_t)/lt)
    end
    sigx = getSigx(xyz_all[:,3])
    rxx = sigx*sigx'.*exp.(-dist_xy/lxy-dist_z/lz).*A_t
    rnn = sign^2 * I
    E = hcat(zeros(Float64,n_flts,n_grds),I)
    y = xyzT_flts[:,end].-yref[idxy]
    
    return rxx,rnn,E,y,length(idx_near)
end

# lon. and lat. of n equally spaced points along the great circle path that originates at (λ1, θ1) and passes through (λ2, θ2)
function findpath(λ1d, θ1d, λ2d, θ2d, n)
  λ1 = deg2rad(λ1d)
  θ1 = deg2rad(θ1d)
  λ2 = deg2rad(λ2d)
  θ2 = deg2rad(θ2d)
  λ2p = atan(cos(θ2)*sin(λ2-λ1), sin(θ1)*cos(θ2)*cos(λ2-λ1) - cos(θ1)*sin(θ2))
  θ2p = asin(cos(θ1)*cos(θ2)*cos(λ2-λ1) + sin(θ1)*sin(θ2))
  R2 = [cos(λ2p) -sin(λ2p) 0; sin(λ2p) cos(λ2p) 0; 0 0 1]
  Ry = [sin(θ1) 0 cos(θ1); 0 1 0; -cos(θ1) 0 sin(θ1)]
  Rz = [cos(λ1) -sin(λ1) 0; sin(λ1) cos(λ1) 0; 0 0 1]
  θp = π/2 .- (0:n-1)/(n-1)*(π/2 - θ2p)
  p = [cos.(θp)'; zeros(n)'; sin.(θp)']
  q = Rz*Ry*R2*p
  λq = atan.(q[2,:], q[1,:])
  θq = asin.(q[3,:])
  return rad2deg.(λq), rad2deg.(θq)
end

# sound speed as a function of in situ temperature and depth
function c(t, S, λ, θ, z)
  if ismissing(t)
    return missing
  else
    p = gsw.p_from_z(z, θ)
    SA = gsw.SA_from_SP(S, p, λ, θ)
    c = gsw.sound_speed_t_exact(SA, t, p)
    return c
  end
end
    
# temperature sensitivity of sound speed
function dcdT(t, S, λ, θ, z)
  if ismissing(t)
    return missing
  else
    eps = 1e-4
    p = gsw.p_from_z(z, θ)
    SA = gsw.SA_from_SP(S, p, λ, θ)
    c1 = gsw.sound_speed_t_exact(SA, t-eps, p)
    c2 = gsw.sound_speed_t_exact(SA, t+eps, p)
    return (c2-c1)/2eps
  end
end


# kernel coordinates
xk = xmink : Δxk : xmaxk
zk = zmink : Δzk : zmaxk
zkK = zminkK : Δzk : zmaxk

# size of kernel array
nxk = length(xk)
nzk = length(zk)
nzkK = length(zkK)

# get T-wave path path coordinates
λk, θk = findpath(λ1, θ1, λ2, θ2, nxk)

# read data
ds = Dataset(@sprintf("data/ref/ref_%s_061921.nc", name))
T̄ = Array{Union{Missing, Float64}}(ds["T"]'); replace!(T̄, NaN=>missing)
c̄ = Array{Union{Missing, Float64}}(ds["c"]'); replace!(c̄, NaN=>missing)
dcdT0 =  Array{Union{Missing, Float64}}(ds["dcdT"]'); replace!(dcdT0, NaN=>missing)

# read kernel file
if size(frq,1) == 1
    kernelfile = "data/knl/kernel2D.txt"
    a = readdlm(kernelfile, skipstart=7)
    # read kernel and mask missing values
    Kc = Array{Union{Missing, Float64}}(1e-7reshape(a[:,3], (nxk, nzkK))); Kc[Kc.==0.] .= missing
    # temperature kernel
    KT = Kc[:,zkK.>=zmink]./c̄.*dcdT0
    #KT = Kc./c̄.*dcdT.(Θ̄, S̄, λk, θk, zk')
else
    KTs = Array{Union{Missing, Float64}}(undef, size(frq,1), nxk, nzk)
    for (i, f) in enumerate(frq)
        kernelfile = "data/knl/kernel_"*name*"_"*f*".txt"
        a = readdlm(kernelfile, skipstart=4)
        Kc = Array{Union{Missing, Float64}}(1e-7reshape(a[:,3], (nxk, nzk))); Kc[Kc.==0.] .= missing; replace!(Kc, NaN=>missing)
        # temperature kernel
        KTs[i,:,:] = Kc./c̄.*dcdT0
        #KTs[i] = Kc./c̄.*dcdT.(Θ̄, S̄, λk, θk, zk')
    end 
end

K1d = vec(KT[end:-1:1,end+1-nzk:end]); replace!(K1d, missing=>0)
ll_path = hcat(λk, θk)
xy_grds = repeat(ll_path,nzk)
z_grds = vec(repeat(zk,1,nxk)')
xyz_grds = hcat(xy_grds,z_grds)
lt,Δday = 30,5 # correlation days
region_params = Array{Union{String, Float64}}([140,168,18.,40.,0.,2000.])
dist_max = 500
mp_lvl = 50:50:2000
lxy,lz,omg,sign = 200,1e3,2*pi/365.2425,1e-3
ecco_ref = Dataset("data/ref/ref_jpwk_map_May21.nc")
head = "/central/groups/oceanphysics/shirui/SOT"
events = DataFrame(CSV.File(@sprintf("%s/results/dtau_jpn_H11N.csv", head)))
#pairs = DataFrame(CSV.File(@sprintf("%s/pairs_jpn_H11N.csv", head)))
# create directory if needed
mkpath(@sprintf("%s/data/argo/%s", head, name))
argos = DataFrame(time=String[], dtau=Float64[], err=Float64[], dtau0=Float64[], err0=Float64[])
# loop over events
for e in eachrow(events)
    dttm = e.time 
    #dttm = DateTime("2012-08-20T14:09:31.83", df)
    # formatted date and time
    fmttime = Dates.format(dttm, "yyyy-mm-ddTHH:MM:SS.ss")
    # filename
    filename = @sprintf("%s/data/argo/%s/%s.h5", head, name, fmttime)
    
    offset,sig_dtau = getTauErr(filename,fmttime,KT,K1d,Δxk,Δzk,region_params,dttm,Δday,ll_path,dist_max,mp_lvl,xyz_grds,lxy,lz,lt,omg,sign,ecco_ref)
    @printf("Offset: %.4f, Error: %.4f s\n", offset,sig_dtau)
    
    # record in catalog
    push!(argos, [fmttime,offset,sig_dtau^0.5,e.dtauT,e.err])
    
    @printf("\n")
end

# save catalog to file
CSV.write(@sprintf("%s/data/catalogs/argo_%s.csv", head, name), argos)
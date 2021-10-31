### ECCO/RG Clim estimate calculation

using PyCall, DelimitedFiles, Dates, NCDatasets, Printf, Statistics, Interpolations, PyPlot, ProgressBars,HDF5
np = pyimport("numpy")
gsw = pyimport("gsw")
gcc = pyimport("great_circle_calculator.great_circle_calculator")
cmo = pyimport("cmocean.cm")

const earthradius = 6371e3
const p10 = ((166.90986,19.71786) .+ (166.89109,19.71356) .+ (166.89677,19.73115))./3

data = "ecco"      # data source from Argo, RGClim or ECCO
date = "Sep14"     # Date of computation
timestep = "month" # time series in every month, day, or t-wave event
getcrho = true
getref = true      # compute reference state (S,T) from ECCO
getKTs = false     # save 2D temperature kernel profiles
getanl = false      # compute travel time anomaly estimate from ARGO or ECCO
missdeep = false    # ignore deep extrapolation below 2km for ARGO
plot = false        # plot intersection profiles

# all events
if data == "rgclim"
    events = Date(2004,1,16):Month(1):Date(2018,12,16)
elseif data == "ecco"
    if timestep == "month"
        events = Date(2008,1):Month(1):Date(2017,12)
    elseif timestep == "day"
        events = DateTime(2008, 03, 01, 12):Day(1):DateTime(2008, 03, 10, 12)
    end
end

# start-end of path
name = "H11N3"
nend = "jpn39"

# frequency/ies 
#frq = ["2.5Hz","3.25Hz","4Hz"]
frq = ["2Hz","3Hz","4Hz"]

# start and end of path
xmink = -300e3
xmaxk = 3500e3
zmink = -8000.
zmaxk = 0.

# kernel grid spacing
Δxk = 10e3
Δzk = 100.

# kernel end points
#p10 = (114.1361, -34.8832) 
#p10 = (166.90986,19.71786)#(72.49, -7.65)
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

# fill in coastal and bottom points (ignores points at the edge of the array)
function filledges!(a)
  # fill in all points that have a missing value with the average of the eight surrounding points (if there are any values there)
  nx, ny, nz = size(a)
  a[2:nx-1,2:ny-1,1:nz] = [ismissing(a[i,j,k]) & !all(ismissing.(a[i-1:i+1,j-1:j+1,k])) ? mean(skipmissing(a[i-1:i+1,j-1:j+1,k])) :
                           a[i,j,k] for i = 2:nx-1, j = 2:ny-1, k = 1:nz]
  # fill in bottom points
  a[:,:,1:nz-1] = [ismissing(a[i,j,k]) & !ismissing(a[i,j,k+1]) ? a[i,j,k+1] : a[i,j,k] for i = 1:nx, j = 1:ny, k = 1:nz-1]
end

# fill in coastal and bottom points (ignores points at the edge of the array)
function filledgesxz!(a)
  # fill in all points that have a missing value with the average of the eight surrounding points (if there are any values there)
  nx, nz = size(a)
  a[2:nx-1,1:nz-1] = [ismissing(a[i,k]) & !all(ismissing.(a[i-1:i+1,k:k+1])) ? mean(skipmissing(a[i-1:i+1,k:k+1])) :
                      a[i,k] for i = 2:nx-1, k = 1:nz-1]
                      
  # fill in bottom points
  a[:,1:nz-1] = [ismissing(a[i,k]) & !ismissing(a[i,k+1]) ? a[i,k+1] : a[i,k] for i = 1:nx, k = 1:nz-1]
end

# fill in bottom points along the path
function fillbottom!(a)
  # fill in all points that have a missing value with the average of the two surrounding points (if there are any values there)
  nx, nz = size(a)
  while sum([ismissing(a[i,k]) for i = 1:nx, k = 1:nz-1])>0
    # fill in bottom points
    a[:,1:nz-1] = [ismissing(a[i,k]) & !ismissing(a[i,k+1]) ? a[i,k+1] : a[i,k] for i = 1:nx, k = 1:nz-1]
  end
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

# in situ temperature from potential temperature
function T(Θ, S, λ, θ, z)
  if ismissing(Θ)
    return missing
  else
    p = gsw.p_from_z(z, θ)
    SA = gsw.SA_from_SP(S, p, λ, θ)
    CT = gsw.CT_from_pt(SA, Θ)
    T = gsw.t_from_CT(SA, CT, p)
    return T
  end
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

# in situ density as a function of salinity and temperature
function rho(t, S, λ, θ, z)
  if ismissing(t)
    return missing
  else
    p = gsw.p_from_z(z, θ)
    SA = gsw.SA_from_SP(S, p, λ, θ)
    rho = gsw.rho_t_exact(SA, t, p)
    return rho
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

function getJP(ds,var;ll=false)
  λe = ds["XC"][:,:,6]
  θe = ds["YC"][:,:,6]
  idx = intersect(findall(x->x>120, λe),findall(x->50>x>10, θe))
  λjp1 = unique(λe[idx])
  θjp = unique(θe[idx])
  dt1 = Array{Union{Missing, Float64}}(ds[var][unique([x[1] for x in idx]),unique([x[2] for x in idx]),6,:,1])[:,:,end:-1:1]; dt1[dt1.==0] .= missing
  λe = ds["XC"][:,:,8]
  θe = ds["YC"][:,:,8]
  idx = intersect(findall(x->x>120, λe),findall(x->50>x>10, θe))
  λjp2 = unique(λe[idx])
  λjp = vcat(λjp1,λjp2) 
  dt2 = Array{Union{Missing, Float64}}(ds[var][unique([x[1] for x in idx]),unique([x[2] for x in idx]),8,:,1])[end:-1:1,:,end:-1:1]; dt2[dt2.==0] .= missing
  dt = vcat(dt1,permutedims(dt2, [2, 1, 3]))
  if ll == true
    return λjp,θjp,dt
  else
    return dt
  end
end

# kernel coordinates
xk = xmink : Δxk : xmaxk
zk = zmink : Δzk : zmaxk

# size of kernel array
nxk = length(xk)
nzk = length(zk)

# get T-wave path path coordinates
λk, θk = findpath(λ1, θ1, λ2, θ2, nxk)

if getanl
    # read data
    ds = Dataset(@sprintf("data/ref/ref_%s_%s.nc", name, date))
    T̄ = Array{Union{Missing, Float64}}(ds["temp"]); replace!(T̄, NaN=>missing)
    S̄ = Array{Union{Missing, Float64}}(ds["salt"]); replace!(S̄, NaN=>missing)
    c̄ = c.(T̄, S̄, λk, θk, zk')
    dcdT0 = dcdT.(T̄, S̄, λk, θk, zk')
    
    # read kernel file
    if size(frq,1) == 1
        kernelfile = "data/knl/kernel2D.txt"
        a = readdlm(kernelfile, skipstart=7)
        # read kernel and mask missing values
        Kc = Array{Union{Missing, Float64}}(1e-7reshape(a[:,3], (nxk, nzk))); Kc[Kc.==0.] .= missing
        # temperature kernel
        KT = Kc./c̄.*dcdT0
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
end

if getKTs
    filename = @sprintf("data/knl/KTs_%s_%s.nc", name, date)
    ds = Dataset(filename,"c")
    
    # Define the dimension 
    defDim(ds,"f",size(frq,1))
    defDim(ds,"x",nxk)
    defDim(ds,"z",nzk)
    
    # Define the variables 
    vKT = defVar(ds,"KT",Float64,("f","x","z"))
    vKT[:,:,:] = KTs
    close(ds)
end

if data == "rgclim"
    head = "/central/groups/oceanphysics/shirui/argo/"
    dsTa = Dataset(@sprintf("%stemp/RG_ArgoClim_Temperature_2019.nc", head))
    #dsSa = Dataset(@sprintf("%sRG_ArgoClim_Salinity_2019.nc", head))
    λa = [dsTa["LONGITUDE"][i]>180.0 ? dsTa["LONGITUDE"][i]-360.0 : dsTa["LONGITUDE"][i] for i = 1:length(dsTa["LONGITUDE"])]
    Δλ = 1.
    Δθ = 1.
    iλ = findfirst(min(λ1,λ2)-3Δλ .< λa)
    iθ = findfirst(min(θ1,θ2)-3Δθ .< dsTa["LATITUDE"])
    nλ = Int(round((abs(λ1-λ2)+6Δλ)/Δλ))
    nθ = Int(round((abs(θ1-θ2)+6Δθ)/Δθ))
    λa = λa[iλ:iλ+nλ-1]
    θa = dsTa["LATITUDE"][iθ:iθ+nθ-1]
    # reverse depth coordinate
    za = -Array(dsTa["PRESSURE"])[end:-1:1]
    T̄a = Array{Union{Missing, Float64}}(dsTa["ARGO_TEMPERATURE_MEAN"][iλ:iλ+nλ-1,iθ:iθ+nθ-1,:])[:,:,end:-1:1]; replace!(T̄a, NaN=>missing)
    ΔTa = Array{Union{Missing, Float64}}(dsTa["ARGO_TEMPERATURE_ANOMALY"][iλ:iλ+nλ-1,iθ:iθ+nθ-1,:,:])[:,:,end:-1:1,:]; replace!(ΔTa, NaN=>missing)
    #S̄a = Array{Union{Missing, Float64}}(dsTa["ARGO_SALINITY_MEAN"][:,iθ:iθ+nθ-1,iλ:iλ+nλ-1])[end:-1:1,:,:]; S̄a[S̄a.==0] .= missing
    #ΔSa = Array{Union{Missing, Float64}}(dsTa["ARGO_SALINITY_ANOMALY"][:,:,iθ:iθ+nθ-1,iλ:iλ+nλ-1])[:,end:-1:1,:,:]; ΔSa[ΔSa.==0] .= missing
end

let
    if getanl
        global T̄
        Δτ = Array{Union{Missing, Float64}}(undef, length(events),size(frq,1))
    end
    if getref || getcrho
        Θsum = zeros(Union{Missing, Float64},nxk, nzk)
        Ssum = zeros(Union{Missing, Float64},nxk, nzk)
    end
    @printf("Interpolation begin...\n")
    for (e, event) in enumerate(tqdm(events))
        if data == "rgclim"
            Te = T̄a .+ ΔTa[:,:,:,e]
            #Sa = permutedims(S̄a .+ ΔSa[e],[3,2,1])
            # ECCO data array size
            nxe, nye, nze = size(Te)
            # fill in the coastal and bottom points
            filledges!(Te)
            #filledges!(Sa)
            # interpolate horizontally onto great circle path
            knots = (λa, θa,)
            Tkze = Array{Union{Missing, Float64}}(undef, nxk, nze)
            #Skze = Array{Union{Missing, Float64}}(undef, nxk, nze)
            for k = 1:nze
                itpT = interpolate(knots, Te[:,:,k], Gridded(Constant()))
                #itpS = interpolate(knots, Se[:,:,k], Gridded(Constant()))
                etpT = extrapolate(itpT, Flat())
                #etpS = extrapolate(itpS, Flat())
                Tkze[:,k] = etpT.(λk, θk)
                #Skze[:,k] = etpS.(λk, θk)
            end
            # interpolate vertically onto kernel grid
            knots = (za,)
            Tk = Array{Union{Missing, Float64}}(undef, nxk, nzk)
            for i = 1:nxk
                itpT = interpolate(knots, Tkze[i,:], Gridded(Constant()))
                #itpS = interpolate(knots, Skze[i,:], Gridded(Constant()))
                etpT = extrapolate(itpT, Flat())
                #etpS = extrapolate(itpS, Flat())
                Tk[i,:] = etpT.(zk)
                #Sk[i,:] = etpS.(zk)
            end
        elseif data == "ecco"
            head = "/central/groups/oceanphysics/shirui/ecco/v4r4/"
            if timestep == "month"
                dsΘ = Dataset(@sprintf("%snctiles_monthly/theta/THETA_%04d_%02d.nc",head,year(event), month(event)))
                dsS = Dataset(@sprintf("%snctiles_monthly/salt/SALT_%04d_%02d.nc", head,year(event), month(event)))
            elseif timestep == "day"
                dsΘ = Dataset(@sprintf("%snctiles_daily/theta/THETA_%04d_%02d_%02d.nc",head,year(event), month(event), day(event)))
                dsS = Dataset(@sprintf("%snctiles_daily/salt/SALT_%04d_%02d_%02d.nc", head,year(event), month(event), day(event)))
            end
            if timestep == "twave"
                # find time stamps of ECCO data bracketing the event
                day1 = round(event - Day(1), Day) + Hour(12)
                day2 = round(event, Day) + Hour(12)
                # load ECCO data
                dsΘ1 = Dataset(@sprintf("%stheta/THETA_%04d_%02d_%02d.nc",head, year(day1), month(day1), day(day1)))
                dsΘ2 = Dataset(@sprintf("%stheta/THETA_%04d_%02d_%02d.nc",head, year(day2), month(day2), day(day2)))
                dsS1 = Dataset(@sprintf("%ssalt/SALT_%04d_%02d_%02d.nc", head, year(day1), month(day1), day(day1)))
                dsS2 = Dataset(@sprintf("%ssalt/SALT_%04d_%02d_%02d.nc", head, year(day2), month(day2), day(day2)))
                # convert coordinates to radiuas
                λe = dsΘ1["XC"][:,:,5]
                θe = dsΘ1["YC"][:,:,5]
                # reverse depth coordinate
                ze = Array(dsΘ1["Z"])[end:-1:1]
                # mask missing data and reverse depth coordinate
                Θ1 = Array{Union{Missing, Float64}}(dsΘ1["THETA"][:,:,5,:,1])[:,:,end:-1:1]; Θ1[Θ1.==0] .= missing
                Θ2 = Array{Union{Missing, Float64}}(dsΘ2["THETA"][:,:,5,:,1])[:,:,end:-1:1]; Θ2[Θ2.==0] .= missing
                S1 = Array{Union{Missing, Float64}}(dsS1["SALT"][:,:,5,:,1])[:,:,end:-1:1]; S1[S1.==0] .= missing
                S2 = Array{Union{Missing, Float64}}(dsS2["SALT"][:,:,5,:,1])[:,:,end:-1:1]; S2[S2.==0] .= missing
                # interpolate to the time of the event
                w1 = (day2 - event)/Millisecond(Day(1))
                w2 = (event - day1)/Millisecond(Day(1))
                Θe = w1*Θ1 + w2*Θ2
                Se = w1*S1 + w2*S2
            else
                # reverse depth coordinate
                ze = Array(dsΘ["Z"])[end:-1:1]
                # mask missing data and reverse depth coordinate
                if name[1:3] == "H11"
                    @printf("Japan interpolation for %s\n",event)
                    λe,θe,Θe = getJP(dsΘ,"THETA";ll=true)
                    Se = getJP(dsS,"SALT")
                else
                    # coordinates
                    λe = dsΘ["XC"][:,1,5]
                    θe = dsΘ["YC"][1,:,5]
                    Θe = Array{Union{Missing, Float64}}(dsΘ["THETA"][:,:,5,:,1])[:,:,end:-1:1]; Θe[Θe.==0] .= missing
                    Se = Array{Union{Missing, Float64}}(dsS["SALT"][:,:,5,:,1])[:,:,end:-1:1]; Se[Se.==0] .= missing
                end
            end
            # ECCO data array size
            nxe, nye, nze = size(Θe)
            # fill in the coastal and bottom points
            filledges!(Θe)
            filledges!(Se)
            # interpolate horizontally onto great circle path
            knots = (λe, θe,)
            Θkze = Array{Union{Missing, Float64}}(undef, nxk, nze)
            Skze = Array{Union{Missing, Float64}}(undef, nxk, nze)
            for k = 1:nze
                itpΘ = interpolate(knots, Θe[:,:,k], Gridded(Constant()))
                itpS = interpolate(knots, Se[:,:,k], Gridded(Constant()))
                etpΘ = extrapolate(itpΘ, Flat())
                etpS = extrapolate(itpS, Flat())
                Θkze[:,k] = etpΘ.(λk, θk)
                Skze[:,k] = etpS.(λk, θk)
            end
            # fill in the bottom points
            filledgesxz!(Θkze)
            filledgesxz!(Skze)
            fillbottom!(Θkze)
            fillbottom!(Skze)
            # interpolate vertically onto kernel grid
            knots = (ze,)
            Θk = Array{Union{Missing, Float64}}(undef, nxk, nzk)
            Sk = Array{Union{Missing, Float64}}(undef, nxk, nzk)
            for i = 1:nxk
                itpΘ = interpolate(knots, Θkze[i,:], Gridded(Constant()))
                itpS = interpolate(knots, Skze[i,:], Gridded(Constant()))
                etpΘ = extrapolate(itpΘ, Flat())
                etpS = extrapolate(itpS, Flat())
                Θk[i,:] = etpΘ.(zk)
                Sk[i,:] = etpS.(zk)
            end
            if getanl
                Tk = [T(Θk[i,k], Sk[i,k], λk[i], θk[i], zk[k]) for i = 1:nxk, k = 1:nzk]
            end
            if getref || getcrho
                Θsum = Θsum .+ Θk
                Ssum = Ssum .+ Sk
            end
        end
    
        if getanl
            ΔT = Array{Union{Missing, Float64}}(Tk - T̄)
            if missdeep
                ΔT[:,zk.<-2e3] .= missing
            end
            for i = 1:size(frq,1)
                Δτ[e,i] = sum(skipmissing(KTs[i,:,:].*ΔT))*Δxk*Δzk
            end
            if (e == length(events) ÷ 3) && plot && (data == "rgargo")
                fig = figure()
                pcolormesh(1e-3xk, zk, replace(T̄, missing=>NaN)', vmin=0, vmax=25, cmap="Reds", shading="auto")
                plt.colorbar()
                fig.tight_layout()
                fig.savefig(string("result/rgargo_test_Tbar.pdf"))
                fig = figure()
                pcolormesh(1e-3xk, zk, replace(asinh.(ΔT/.1), missing=>NaN)', vmin=asinh(-100), vmax=asinh(100), cmap="RdBu_r", shading="auto")
                fig.tight_layout()
                fig.savefig(string("result/rgargo_test_dT.pdf"))
                fig = figure()
                pcolormesh(1e-3xk, zk, replace(asinh.(KTs[1,:,:].*ΔT/1e-11), missing=>NaN)', vmin=asinh(-100), vmax=asinh(100), cmap="RdBu", shading="auto")
                fig.savefig(string("result/rgargo_test_dtau.pdf"))
            end
        end
    end
    @printf("Interpolation finished!\n")
    if getref
        replace!(Θsum, missing=>NaN) 
        replace!(Ssum, missing=>NaN)
        filename = @sprintf("data/ref/ref_%s_%s.nc", name, nend)
        ds = Dataset(filename,"c")
        
        # Define the dimension "lon" and "lat" 
        defDim(ds,"x",nxk)
        defDim(ds,"z",nzk)
        
        # Define the variables temperature and salinity with the attribute units
        vT = defVar(ds,"temp",Float64,("x","z"))#, attrib = OrderedDict("units" => "degree Celsius"))
        vS = defVar(ds,"salt",Float64,("x","z"))#, attrib = OrderedDict("units" => "psu")))
        vT[:,:] = T.(Θsum/Float64(length(events)), Ssum/Float64(length(events)), λk, θk, zk')
        vS[:,:] = Ssum/Float64(length(events))
        close(ds)
    end
    if getcrho
        T̄ = T.(Θsum/Float64(length(events)), Ssum/Float64(length(events)), λk, θk, zk')
        S̄ = Ssum/Float64(length(events))
        c̄ = c.(T̄, S̄, λk, θk, zk')
        arho = rho.(T̄, S̄, λk, θk, zk')
        filename = @sprintf("data/ref/%s_%s.txt",name,nend)
        open(filename, "w") do io
           write(io, @sprintf("%.2f %.2f %.2f %.2f #x_min(m) z_min(m) x_max(m) z_max(m)\n",xmink,zmink,xmaxk,zmaxk))
           write(io, @sprintf("%.4f %.4f                #dx(m) dz(m)\n",Δxk,Δzk))
           write(io, @sprintf("%d %d                           #nx nz\n",nxk,nzk))
           for i = 1:nxk, k = 1:nzk
               write(io, @sprintf("%.2f %.2f %.4f 0 %.4f\n",xk[i],abs(zk[end+1-k]),c̄[i,end+1-k],arho[i,end+1-k]))
           end
        end
    end
    if getanl
        filename = @sprintf("result/dtaus_%s_%s.nc", name, date)
        ds = Dataset(filename,"c")
        
        # Define the dimension 
        defDim(ds,"t",length(events))
        defDim(ds,"f",size(frq,1))
        
        # Define the variables 
        vτ = defVar(ds,"dtau",Float64,("t","f"))#, attrib = OrderedDict("units" => "second"))
        vτ[:,:] = Δτ
        close(ds)
        
        fig, ax = subplots(1, 1, figsize=(16, 6.4))
        for i = 1:size(frq,1)
            ax.plot(events, Δτ[:,i], linewidth=1.0,label=frq[i])
        end
        ax.invert_yaxis()
        ax.legend(frameon=false, loc=1)
        ax.set_ylabel("travel time anomaly (s)")
        fig.tight_layout()
        fig.savefig(string("result/ts_",name,"_",data,".pdf"))
    end
end

using SOT, PyPlot, Printf, Dates, LinearAlgebra, Statistics, SparseArrays, HDF5, Interpolations, DataFrames, CSV
using Random, Distributions,StatsBase, NCDatasets, SeisNoise,StatsBase

source = "nias_tm"
station = "H08"
tstations = [station]
evtpos = [1.62, 96.92]
ppfile = @sprintf("results/%s_%s_ppairs.csv",source,station)
tpvfile = @sprintf("results/%s_%s_tpairs.csv",source,station)

tstalocs = [-7.65, 72.49]
pstations = ["PS.PSI..BHZ","MY.KUM..BHZ","II.WRAB.00.BHZ", "GE.GSI..BHZ"]
pstnlats,pstnlons = [2.69,5.2902,-19.934,1.3039],[98.92,100.6492,134.36,97.5755]
pstalocs = DataFrame(station=pstations,slat=pstnlats,slon=pstnlons)

# frequencies used in inversion
tinvfreq = 2.0:1.0:4.0

# number of frequencies
l = length(tinvfreq)

# manually exclude pairs
tpairs = CSV.read(tpvfile, DataFrame)
ppairs = CSV.read(ppfile, DataFrame)
#t0,t1 = DateTime(2005,1,1),DateTime(2006,1,1)
#tpairs = tpairs[(t1.> tpairs.event1.≥ t0),:]
#tpairs = tpairs[(tpairs.event2.< t1),:]
ppairs = innerjoin(ppairs, select(tpairs, [:event1, :event2]), on=[:event1, :event2])

# number of good T- and P-wave pairs
nt = size(tpairs, 1)
np = size(ppairs, 1)

# number of unique events
m = length(unique([tpairs.event1; tpairs.event2]))

tpairs.Δτl = [[parse(Float64, s) for s in split(tpairs.Δτl[i][2:end-1],", ")] for i = 1:nt]
tpairs.Δτc = [[parse(Float64, s) for s in split(tpairs.Δτc[i][2:end-1],", ")] for i = 1:nt]
tpairs.Δτr = [[parse(Float64, s) for s in split(tpairs.Δτr[i][2:end-1],", ")] for i = 1:nt]
tpairs.ccc = [[parse(Float64, s) for s in split(tpairs.ccc[i][2:end-1],", ")] for i = 1:nt]
tpairs.ccl = [[parse(Float64, s) for s in split(tpairs.ccl[i][2:end-1],", ")] for i = 1:nt]
tpairs.ccr = [[parse(Float64, s) for s in split(tpairs.ccr[i][2:end-1],", ")] for i = 1:nt]

###
model,Kgrid = "ecco","Gaussian"
fl = @sprintf("data/%s/covtau0_%s.h5",model,Kgrid)
lags = h5read(fl, "lags")
ctau = h5read(fl, "covs")

# solution standard deviation for covariance (s)
σc = ones(6)#[1.686,0.508,0,1.672,1.498,1.642]

# noise (s)
σn,σnp = 0.018,3e-3

# origin time correction standard deviation (s)
σp = 1.3

# trend priors for coefficients of travel time anomalies (s/day)
σtrend = 1e-2*ones(3)/SOT.meanyear

σs,σh = 0.03,0.03

# get inversion matrices
t, E, R, N, P, D, invR = SOT.invert(tpairs, ppairs, 1, σc, σn, σp, I, I, 1; cov_method=:complex, tstation=tstalocs, pstations=pstalocs, σs, σnp,σh,lags,ctau, σtrend)

# make cycle-skipping correction
returnz = true
tpairs.Δτ,tpairs.Δτp,zint,zinp = SOT.correctcycleskipping(tpairs, ppairs, E, R, N, P, m; source, tstations, returnz)

#CSV.write(@sprintf("results/pairs/%s_%s_tpairs.csv",source,tstations[1]), tpairs)

# collect delays into data vector
y = [reshape(vcat([(tpairs.Δτ[i])' for i = 1:nt]...), l*nt); ppairs.Δτ]

# invert
a = P*E'*inv(Array(N))*y

# extract trends
trends = a[(l+1)*m+1:(l+1)*m+l]
Ptrends = diag(P[(l+1)*m+1:(l+1)*m+l,(l+1)*m+1:(l+1)*m+l])
@printf("\ntrend = %s s/yr,\n ptrend = %s s/yr\n",trends*SOT.meanyear,2*sqrt.(Ptrends)*SOT.meanyear)

# plot residuals with and without smoothing
n1 = y - E*a
n2 = y - E*(E\y)
tpairs.n1,tpairs.n2 = collect(eachrow(reshape(n1[1:l*nt],(nt,l)))),collect(eachrow(reshape(n2[1:l*nt],(nt,l))))
ppairs.n1,ppairs.n2 = n1[l*nt+1:l*nt+np],n2[l*nt+1:l*nt+np]

tn1maxi = argmax(abs.(reshape(n1[1:l*nt],(nt,l))),dims=1)
tn2maxi = argmax(abs.(reshape(n2[1:l*nt],(nt,l))),dims=1)
for i = 1:3
  @printf("\ntpair outlier: %s %s,\n Dt=%s, Dtp=%.3fs, n1=%s s\n", tpairs.event1[tn1maxi[i][1]],tpairs.event2[tn1maxi[i][1]],tpairs.Δτ[tn1maxi[i][1]],tpairs.Δτp[tn1maxi[i][1]],tpairs.n1[tn1maxi[i][1]])
  @printf("tpair outlier: %s %s,\n Dt=%s, Dtp=%.3fs, n2=%s s\n", tpairs.event1[tn2maxi[i][1]],tpairs.event2[tn2maxi[i][1]],tpairs.Δτ[tn2maxi[i][1]],tpairs.Δτp[tn2maxi[i][1]],tpairs.n2[tn2maxi[i][1]])
end
pn1maxi = argmax(abs.(ppairs.n1))
@printf("\nppair n1 std: %.3f; outlier: %s %s, n1=%.2f s\n", std(ppairs.n1), ppairs.event1[pn1maxi],ppairs.event2[pn1maxi],ppairs.n1[pn1maxi])
pn2maxi = argmax(abs.(ppairs.n2))
@printf("ppair n2 std: %.3f; outlier: %s %s, n2=%.2f s\n", std(ppairs.n2), ppairs.event1[pn2maxi],ppairs.event2[pn2maxi],ppairs.n2[pn2maxi])

b = -0.25:0.005:0.25
#@assert all((b[1] .< n1 .< b[end]) .& (b[1] .< n2 .< b[end]))
fig, ax = subplots(2, l+1; sharex=true, sharey=true)
for i = 1:l
  ax[1,i].hist(n1[(i-1)*nt+1:i*nt]; bins=b)
  ax[2,i].hist(n2[(i-1)*nt+1:i*nt]; bins=b)
  ax[1,i].set_title("T wave, freq. $i")
  ax[2,i].set_title("T wave, freq. $i")
end
ax[1,l+1].hist(n1[l*nt+1:l*nt+np]; bins=b)
ax[2,l+1].hist(n2[l*nt+1:l*nt+np]; bins=b)
ax[1,l+1].set_title("P wave")
ax[2,l+1].set_title("P wave")
ax[1,1].set_yscale("log")
ax[1,1].set_ylabel("count")
ax[2,1].set_ylabel("count")
ax[2,1].set_xlabel("residual (s)")
ax[2,2].set_xlabel("residual (s)")
ax[2,3].set_xlabel("residual (s)")
ax[2,4].set_xlabel("residual (s)")
fig.tight_layout()
fig.savefig(@sprintf("results/plots/%s_%s_residual_hist_%s.pdf",source,tstations[1],Kgrid))

td = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24

# reconstruct full travel time anomalies
τ = reshape(D*a, (m, l))
eτ = reshape(sqrt.(diag(D*P*D')), (m, l))

# interpolate onto regular grid
ti = DateTime(2005, 1, 1, 12, 0, 0) : Day(1) : DateTime(2005, 12, 31, 12, 0, 0)
τi,Pi,ΔCi,zix,zix2 = SOT.regulargrid_gpcov(td, ti, a, R, invR, P, σc, ctau, lags)

n = length(ti)
  
###
nmode = 5
Lx = 29e5

nanmean(x) = mean(filter(!isnan,x))
nanmean(x,y) = mapslices(nanmean,Array{Float64}(x);dims=y)
fl = Dataset("data/ecco/n2_Nias_H08_2005.nc", "r")
range = Array{Float64}(fl["x"][:])
depth = -Array{Float64}(fl["z"][:])
idx = 0 .<= range .<= 29e5
nx = length(range[idx])
Δx = range[2] - range[1]
Δz = depth[1] - depth[2]
giofes = NCDataset("results/ofes/modes_mean_ofes.nc")["hn"][:,1:nmode]
n2ofes = nanmean(replace(NCDataset("data/ofes/n2_Nias_H08_2005.nc")["n2"][:,idx], missing=>NaN), 2)
giecco = NCDataset("results/ecco/modes_mean_ecco.nc")["hn"][:,1:nmode]
n2ecco = nanmean(replace(NCDataset("data/ecco/n2_Nias_H08_2005.nc")["n2"][:,idx], missing=>NaN), 2)

if Kgrid == "Gaussian"
  # kernels
  fl = Dataset("data/temperature/KTs_Nias_H08_Gaussian.nc", "r")
  K0 = Array{Float64}(replace(fl["SEMkernels_T"][:], missing=>0))
  xk = Array{Float64}(fl["x"][:])
  zk = Array{Float64}(fl["z"][:])
  K = Array{Float64, 3}(undef, length(depth), length(range), 3)
  for i = 1:3
    interp_linear = linear_interpolation((zk, xk), K0[:,:,i])
    for j = 1:length(range)
      K[:,j,i] = interp_linear.(-depth, range[j])
    end
  end
else
  K = h5read("data/temperature/nias_H08.h5", "K")
end
K = sum(K, dims=2)[:,1,:]*Δx

Aofes = Array{Float64, 2}(undef, 3, nmode)
Aecco = Array{Float64, 2}(undef, 3, nmode)
N0 = 1e-3
for i = 1:3
  for j = 1:nmode
    Aofes[i,j] = sum(replace(K[:,i].*n2ofes.*giofes[:,j],missing => 0))*Δz/N0
    Aecco[i,j] = sum(replace(K[:,i].*n2ecco.*giecco[:,j],missing => 0))*Δz/N0
  end
end

# SVD
Uo, Λo, Vo = svd(Aofes)
Ue, Λe, Ve = svd(Aecco)

# project travel time anomalies and trend reconstructions onto singular vectors
T = kron(sparse(inv(Diagonal(Λe))*Ue'), I(n))

fl = Dataset(@sprintf("data/ofes/kTdT_Nias_H08_2005_%s.nc",Kgrid), "r")
τofes = Array{Float64}(replace(fl["SEMdtaus"][:,:], missing=>0))
fl = Dataset(@sprintf("data/ecco/kTdT_Nias_H08_2005_%s.nc",Kgrid), "r")
τecco = Array{Float64}(replace(fl["SEMdtaus"][:,:], missing=>0))

tm = td[1] + (td[m] - td[1])/2
tid = Dates.value.(ti .- DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
for i = 1:l
  τofes[:,i] = detrend(τofes[:,i] .- mean(τofes[:,i]))
  τecco[:,i] = detrend(τecco[:,i] .- mean(τecco[:,i]))
end

# project Argo and ECCO travel time anomalies onto singular vectors
cofes = (τofes)*(Uo*inv(Diagonal(Λo)))
ci = (τi)*(Ue*inv(Diagonal(Λe)))
cecco = (τecco)*(Ue*inv(Diagonal(Λe)))
ec = reshape(sqrt.(max.(diag(T*Pi*T'),0).+max.(diag(T*ΔCi*T'),0)), (n, l))
  
qplot = true
if qplot
    p = 0.5:0.5:99.5
    α = 0.05
    qs4,qs5 = percentile(zix,p),percentile(zint,p)
    qs2,qs3 = percentile(zinp,p),percentile(zix2,p)
    qt = quantile.(Normal(), p/100)
    lb,ub = zeros(length(p),3),zeros(length(p),3)
    for i = 1:length(p)
        for (j,ns) in enumerate([n, l*nt, np])
            d = NoncentralT(ns-1,-sqrt(ns)*qt[i])
            lb[i,j],ub[i,j] = -quantile(d, 1-α/2)/sqrt(ns),-quantile(d, α/2)/sqrt(ns)
        end
    end
    # small font
    rc("font", size=8)
    rc("axes", titlesize="medium")
    fig,ax=plt.subplots(figsize=(4.8,3.6))
    ax.plot(qt,qs3-qt,label="2005 daily anomaly",c="tab:blue")
    #ax.plot(qt,qs3-qt,c="tab:blue",ls="--")
    ax.fill_between(qt, (lb[:,1]-qt), (ub[:,1]-qt), alpha=.2, zorder=3, color="tab:blue", linewidth=0)
    ax.plot(qt,qs5-qt,label="\$T\$-wave residual",c="tab:orange")
    ax.fill_between(qt, (lb[:,2]-qt), (ub[:,2]-qt), alpha=.2, zorder=3, color="tab:orange", linewidth=0)
    ax.plot(qt,qs2-qt,label="\$P\$-wave residual",c="tab:green")
    ax.fill_between(qt, (lb[:,3]-qt), (ub[:,3]-qt), alpha=.2, zorder=3, color="tab:green", linewidth=0)
    ax.axhline(0,color="black",lw=.8)
    ax.set_xlabel("theoretical quantile")
    ax.set_ylabel("sample quantile \$-\$ theoretical quantile")
    ax.legend(frameon=false)
    fig.tight_layout()
    fig.savefig(@sprintf("results/ofes/quantile_%s.pdf",Kgrid))
end

# save catalogs of used pairs
#CSV.write("results/pairs/nias_tm_H08_tpairs.csv", tpairs)
#CSV.write("results/pairs/nias_tm_H08_ppairs.csv", ppairs)

# save to file
h5open(@sprintf("results/anomalies/nias_tm_H08_%s.h5",Kgrid), "w") do file
  write(file, "t", Dates.value.(t .- DateTime(2000, 1, 1, 0, 0, 0)))
  write(file, "ti", Dates.value.(ti .- DateTime(2000, 1, 1, 0, 0, 0)))
  write(file, "τ", τ)
  write(file, "eτ", eτ)
  write(file, "τi", τi)
  write(file, "eτi", reshape(sqrt.(max.(diag(Pi),0).+max.(diag(ΔCi),0)), (n, l)))
  write(file, "ci", ci)
  write(file, "cofes", cofes)
  write(file, "cecco", cecco)
  write(file, "ec", ec)
  write(file, "zix", zix)
  write(file, "zix2", zix2)
  write(file, "zint", zint)
  write(file, "zinp", zinp)
  write(file,"Pi",Pi)
  write(file,"DCi",ΔCi)
end
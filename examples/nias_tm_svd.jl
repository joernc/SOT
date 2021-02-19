using .SOT
using PyPlot, Printf, Dates
using Interpolations, Statistics
using HDF5, DelimitedFiles, NPZ
using LinearAlgebra, SparseArrays

# small font
rc("font", size=8)

# excluded time periods for H08: before 2004-12-01 and periods with uncorrected clock error
excludetimes08 = [[Date(2001, 1, 1) Date(2004, 12, 1)],
                  [DateTime("2010-03-16T00:00:00") DateTime("2010-05-17T02:06:17.760")],
                  [Date(2017, 6, 1) Date(2018, 1, 1)]]

# number of frequencies
l = 3

# read H08 data
h5open("results/nias_tm_H08.h5", "r") do file
  global tr08 = read(file, "t")
  global τ08 = read(file, "tau")
  global τerr08 = read(file, "tauerr")
  global δτ08 = read(file, "dtau")
  global δτerr08 = read(file, "dtauerr")
  global τargo08 = read(file, "tauargo")
  global τecco08 = read(file, "tauecco")
  global Targo08 = read(file, "Targo")
  global Tecco08 = read(file, "Tecco")
end
t08 = DateTime(2000, 1, 1, 0, 0, 0) .+ Millisecond.(Int.(round.(tr08*1000*3600*24)))

# read H01 data
h5open("results/nias_tm_H01.h5", "r") do file
  global tr01 = read(file, "t")
  global τ01 = read(file, "tau")
  global τerr01 = read(file, "tauerr")
  global δτ01 = read(file, "dtau")
  global δτerr01 = read(file, "dtauerr")
  global τargo01 = read(file, "tauargo")
  global τecco01 = read(file, "tauecco")
  global Targo01 = read(file, "Targo")
  global Tecco01 = read(file, "Tecco")
end
t01 = DateTime(2000, 1, 1, 0, 0, 0) .+ Millisecond.(Int.(round.(tr01*1000*3600*24)))

# offset H08 based on trends
tr08 = Dates.value.(t08 .- t08[1])/1000/3600/24
for i = 1:l
  c, _ = SOT.lineartrend(t08, τ08[:,i]; fitannual=true, fitsemiannual=true)
  τ08[:,i] .-= c[2] + c[1]*tr08[end]/2
  if i == 1
    δτ08[:,:] .-= c[2] + c[1]*tr08[end]/2
  else
    δτ08[:,i-1] .+= c[2] + c[1]*tr08[end]/2
  end
  c, _ = SOT.lineartrend(t08, τecco08[:,i]; fitannual=true, fitsemiannual=true)
  τecco08[:,i] .-= c[2] + c[1]*tr08[end]/2
  c, _ = SOT.lineartrend(t08, τargo08[:,i]; fitannual=true, fitsemiannual=true)
  τargo08[:,i] .-= c[2] + c[1]*tr08[end]/2
end

# offset H01 based on trends
tr01 = Dates.value.(t01 .- t01[1])/1000/3600/24
for i = 1:l
  c, _ = SOT.lineartrend(t01, τ01[:,i]; fitannual=true, fitsemiannual=true)
  τ01[:,i] .-= c[2] + c[1]*tr01[end]/2
  if i == 1
    δτ01[:,:] .-= c[2] + c[1]*tr01[end]/2
  else
    δτ01[:,i-1] .+= c[2] + c[1]*tr01[end]/2
  end
  c, _ = SOT.lineartrend(t01, τecco01[:,i]; fitannual=true, fitsemiannual=true)
  τecco01[:,i] .-= c[2] + c[1]*tr01[end]/2
  c, _ = SOT.lineartrend(t01, τargo01[:,i]; fitannual=true, fitsemiannual=true)
  τargo01[:,i] .-= c[2] + c[1]*tr01[end]/2
end

# estimate H08 trends
cτ08, _ = SOT.lineartrend(t08, τ08[:,1]; fitannual=true, fitsemiannual=true)
cδτ08, _ = SOT.lineartrend(t08, δτ08[:,l-1]; fitannual=true, fitsemiannual=true)
cτecco08, _ = SOT.lineartrend(t08, τecco08[:,1]; fitannual=true, fitsemiannual=true)
cδτecco08, _ = SOT.lineartrend(t08, τecco08[:,1] - τecco08[:,l]; fitannual=true, fitsemiannual=true)
cτargo08, _ = SOT.lineartrend(t08, τargo08[:,1]; fitannual=true, fitsemiannual=true)
cδτargo08, _ = SOT.lineartrend(t08, τargo08[:,1] - τargo08[:,2]; fitannual=true, fitsemiannual=true)

# estimate H01 trends
cτ01, _ = SOT.lineartrend(t01, τ01[:,1]; fitannual=true, fitsemiannual=true)
cδτ01, _ = SOT.lineartrend(t01, δτ01[:,l-1]; fitannual=true, fitsemiannual=true)
cτecco01, _ = SOT.lineartrend(t01, τecco01[:,1]; fitannual=true, fitsemiannual=true)
cδτecco01, _ = SOT.lineartrend(t01, τecco01[:,1] - τecco01[:,2]; fitannual=true, fitsemiannual=true)
cτargo01, _ = SOT.lineartrend(t01, τargo01[:,1]; fitannual=true, fitsemiannual=true)
cδτargo01, _ = SOT.lineartrend(t01, τargo01[:,1] - τargo01[:,2]; fitannual=true, fitsemiannual=true)

###

Δx = 5e3
Δz = 50.
c0 = 1.5e3
α = 4.6
depth = 6.5:-1e-3Δz:0
depthi = [depth; -1e-3Δz] .+ 1e-3Δz/2

K08 = h5read("data/kernels/nias_H08.h5", "K")
replace!(K08, NaN=>0)
K08 = sum(K08; dims=2)[:,1,:]'*Δx*Δz

K01 = h5read("data/kernels/nias_H01.h5", "K")
replace!(K01, NaN=>0)
K01 = sum(K01; dims=2)[:,1,:]'*Δx*Δz

nz = size(K08, 2)

###

# SVD
U08, Σ08, V08 = svd(K08)
U01, Σ01, V01 = svd(K01)

# infer temperature profiles
T08 = V08[:,1:2]*diagm(Σ08[1:2].^-1)*U08[:,1:2]'*τ08'
Targoi08 = V08[:,1:2]*diagm(Σ08[1:2].^-1)*U08[:,1:2]'*τargo08'
Teccoi08 = V08[:,1:2]*diagm(Σ08[1:2].^-1)*U08[:,1:2]'*τecco08'
T01 = V01[:,1:2]*diagm(Σ01[1:2].^-1)*U01[:,1:2]'*τ01'
Targoi01 = V01[:,1:2]*diagm(Σ01[1:2].^-1)*U01[:,1:2]'*τargo01'
Teccoi01 = V01[:,1:2]*diagm(Σ01[1:2].^-1)*U01[:,1:2]'*τecco01'

# get τ trends at all frequencies
τtrendtwave08 = [SOT.lineartrend(t08, τ08[:,i]; fitannual=true, fitsemiannual=true)[1][1] for i = 1:l]
τtrendargo08 = [SOT.lineartrend(t08, τargo08[:,i]; fitannual=true, fitsemiannual=true)[1][1] for i = 1:l]
τtrendecco08 = [SOT.lineartrend(t08, τecco08[:,i]; fitannual=true, fitsemiannual=true)[1][1] for i = 1:l]
τtrendtwave01 = [SOT.lineartrend(t01, τ01[:,i]; fitannual=true, fitsemiannual=true)[1][1] for i = 1:l]
τtrendargo01 = [SOT.lineartrend(t01, τargo01[:,i]; fitannual=true, fitsemiannual=true)[1][1] for i = 1:l]
τtrendecco01 = [SOT.lineartrend(t01, τecco01[:,i]; fitannual=true, fitsemiannual=true)[1][1] for i = 1:l]

# infer temperature trend profiles
Ttrendtwave08 = V08[:,1:2]*diagm(Σ08[1:2].^-1)*U08[:,1:2]'*τtrendtwave08
Ttrendargo08 = V08[:,1:2]*diagm(Σ08[1:2].^-1)*U08[:,1:2]'*τtrendargo08
Ttrendecco08 = V08[:,1:2]*diagm(Σ08[1:2].^-1)*U08[:,1:2]'*τtrendecco08
Ttrendtwave01 = V01[:,1:2]*diagm(Σ01[1:2].^-1)*U01[:,1:2]'*τtrendtwave01
Ttrendargo01 = V01[:,1:2]*diagm(Σ01[1:2].^-1)*U01[:,1:2]'*τtrendargo01
Ttrendecco01 = V01[:,1:2]*diagm(Σ01[1:2].^-1)*U01[:,1:2]'*τtrendecco01

Ttrendargof08 = [SOT.lineartrend(t08, Targo08[:,i]; fitannual=true, fitsemiannual=true)[1][1] for i = 1:nz]
Ttrendeccof08 = [SOT.lineartrend(t08, Tecco08[:,i]; fitannual=true, fitsemiannual=true)[1][1] for i = 1:nz]
Ttrendargof01 = [SOT.lineartrend(t01, Targo01[:,i]; fitannual=true, fitsemiannual=true)[1][1] for i = 1:nz]
Ttrendeccof01 = [SOT.lineartrend(t01, Tecco01[:,i]; fitannual=true, fitsemiannual=true)[1][1] for i = 1:nz]

Tmeanargof08 = [SOT.lineartrend(t08, Targo08[:,i]; fitannual=true, fitsemiannual=true)[1][2] for i = 1:nz]
Tmeaneccof08 = [SOT.lineartrend(t08, Tecco08[:,i]; fitannual=true, fitsemiannual=true)[1][2] for i = 1:nz]
Tmeanargof01 = [SOT.lineartrend(t01, Targo01[:,i]; fitannual=true, fitsemiannual=true)[1][2] for i = 1:nz]
Tmeaneccof01 = [SOT.lineartrend(t01, Tecco01[:,i]; fitannual=true, fitsemiannual=true)[1][2] for i = 1:nz]

Targo08 .-= Tmeanargof08' + Ttrendargof08'.*tr08[end]/2
Tecco08 .-= Tmeaneccof08' + Ttrendeccof08'.*tr08[end]/2
Targo01 .-= Tmeanargof01' + Ttrendargof01'.*tr01[end]/2
Tecco01 .-= Tmeaneccof01' + Ttrendeccof01'.*tr01[end]/2

# insert gaps
function insertgap(t, ys, tgap)
  idx = findlast(t .< tgap)
  m = length(t)
  if !isnothing(idx) && idx < m
    tg = [t[1:idx]; tgap; t[idx+1:m]]
    ysg = []
    for y in ys
      l = size(y, 2)
      push!(ysg, [y[1:idx,:]; [NaN for i=1:l]'; y[idx+1:m,:]])
    end
  else
    tg = t
    ysg = ys
  end
  return tg, ysg...
end

# hydrophone gaps
tg08, τg08, τerrg08, δτg08, δτerrg08, τargog08, τeccog08 = insertgap(t08, [τ08, τerr08, δτ08, δτerr08, τargo08, τecco08], Date(2007, 1, 1))
tg01, τg01, τerrg01, δτg01, δτerrg01, τargog01, τeccog01 = insertgap(t01, [τ01, τerr01, δτ01, δτerr01, τargo01, τecco01], Date(2011, 10, 27))

# excluded periods
for i = 1:length(excludetimes08)
  global tg08, τg08, τerrg08, δτg08, δτerrg08, τargog08, τeccog08
  tg08, τg08, τerrg08, δτg08, δτerrg08, τargog08, τeccog08 = insertgap(tg08, [τg08, τerrg08, δτg08, δτerrg08, τargog08, τeccog08], excludetimes08[i][1])
end

# seismic station gap
tg08, τg08, τerrg08, δτg08, δτerrg08, τargog08, τeccog08 = insertgap(tg08, [τg08, τerrg08, δτg08, δτerrg08, τargog08, τeccog08], Date(2009, 7, 1))
tg01, τg01, τerrg01, δτg01, δτerrg01, τargog01, τeccog01 = insertgap(tg01, [τg01, τerrg01, δτg01, δτerrg01, τargog01, τeccog01], Date(2009, 7, 1))

###

tr = Dates.value.(t01 .- DateTime(2000, 1, 1))/1000/3600/24
itp = interpolate((tr,), τ01[:,1], Gridded(Linear()))
thyd = [DateTime(2007, 4, 1); DateTime(2016, 4, 1)]
trhyd = Dates.value.(thyd .- DateTime(2000, 1, 1))/1000/3600/24
τhyd = itp.(trhyd)

depthctd = 1e-3h5read("data/ctd/i08.h5", "depth")
Ttrendctdf = h5read("data/ctd/i08.h5", "dT") / h5read("data/ctd/i08.h5", "years")
Ttrendctdf = [mean(Ttrendctdf[depth[i]-.025 .< depthctd .< depth[i]+.025]) for i = 1:nz]

###

# plot H08 timeseries
colors = matplotlib.rcParams["axes.prop_cycle"].by_key()["color"]
fig, ax = subplots(3, 1, figsize=(190/25.4, 190/25.4), sharex=true)
ax[1].plot(tg08, τg08[:,1], color="tab:blue", linewidth=1, zorder=1, label=L"$T$ waves")
ax[1].scatter(t08, τ08[:,1], s=1, c="tab:blue", zorder=1)
ax[1].fill_between(tg08, τg08[:,1] - 2τerrg08[:,1], τg08[:,1] + 2τerrg08[:,1], alpha=.25,
                   color="tab:blue", linewidths=0, zorder=1)
ax[2].plot(tg08, δτg08[:,l-1], color="tab:blue", linewidth=1, zorder=1, label=L"$T$ waves")
ax[2].scatter(t08, δτ08[:,l-1], s=1, color="tab:blue", zorder=1)
ax[2].fill_between(tg08, δτg08[:,l-1] - 2δτerrg08[:,l-1], δτg08[:,l-1] + 2δτerrg08[:,l-1], alpha=.25,
                   color="tab:blue", linewidths=0, zorder=1)
ax[1].plot(tg08, τargog08[:,1], color="tab:orange", linewidth=1, zorder=0, label="Argo")
ax[1].plot(tg08, τeccog08[:,1], color="tab:green", linewidth=1, zorder=0, label="ECCO")
ax[2].plot(tg08, τargog08[:,1] - τargog08[:,l], color="tab:orange", linewidth=1, zorder=0,
           label="Argo")
ax[2].plot(tg08, τeccog08[:,1] - τeccog08[:,l], color="tab:green", linewidth=1, zorder=0,
           label="ECCO")
ax[1].plot(t08, cτargo08[1]*tr08 .+ cτargo08[2], color="tab:orange", linewidth=1, zorder=0)
ax[1].plot(t08, cτecco08[1]*tr08 .+ cτecco08[2], color="tab:green", linewidth=1, zorder=0)
ax[1].plot(t08, cτ08[1]*tr08 .+ cτ08[2], color="tab:blue", linewidth=1, zorder=1)
ax[2].plot(t08, cδτargo08[1]*tr08 .+ cδτargo08[2], color="tab:orange", linewidth=1, zorder=0)
ax[2].plot(t08, cδτecco08[1]*tr08 .+ cδτecco08[2], color="tab:green", linewidth=1, zorder=0)
ax[2].plot(t08, cδτ08[1]*tr08 .+ cδτ08[2], color="tab:blue", linewidth=1, zorder=1)
idx = findall(isnan.(τg08[:,1]))
for i = 1 : length(idx)+1
  i0 = i == 1 ? 1 : idx[i-1]+1
  i1 = i == length(idx)+1 ? length(t08) : idx[i]-1
  ti = tg08[i0:i1-1] + .5(tg08[i0+1:i1] - tg08[i0:i1-1])
  ti = [tg08[i0] - (ti[1]-tg08[i0]); ti; tg08[i1] + (tg08[i1]-ti[end])]
  ax[3].pcolormesh(ti, depthi, T08[:,i0:i1], cmap="RdBu_r", vmin=-.11, vmax=.11)
end
ax[1].invert_yaxis()
ax[1].legend(frameon=false, loc=4, ncol=3)
ax[1].set_xlim(t08[1] - (t08[end] - t08[1])÷100, t08[end] + (t08[end] - t08[1])÷100)
ax[3].set_ylim(5, 0)
ax[1].set_ylabel("travel time anomaly (s)")
ax[2].set_ylabel("travel time difference (s)")
ax[3].set_ylabel("depth (km)")
axr = ax[1].twinx()
axr.set_yticks(-.06:.030:.06)
axr.set_ylim(ax[1].get_ylim()./sum(K08[1,:]))
axr.set_ylabel("temperature anomaly (K)")
fig.align_ylabels()
fig.tight_layout()

# plot H01 timeseries
colors = matplotlib.rcParams["axes.prop_cycle"].by_key()["color"]
fig, ax = subplots(3, 1, figsize=(190/25.4, 190/25.4), sharex=true)
#ax[1].scatter(thyd, τhyd, s=20, color="tab:red", zorder=10)
ax[1].plot(tg01, τg01[:,1], color="tab:blue", linewidth=1, zorder=1, label=L"$T$ waves")
ax[1].scatter(t01, τ01[:,1], s=1, c="tab:blue", zorder=1)
ax[1].fill_between(tg01, τg01[:,1] - 2τerrg01[:,1], τg01[:,1] + 2τerrg01[:,1], alpha=.25,
                   color="tab:blue", linewidths=0, zorder=1)
ax[2].plot(tg01, δτg01[:,l-1], color="tab:blue", linewidth=1, zorder=1, label=L"$T$ waves")
ax[2].scatter(t01, δτ01[:,l-1], s=1, color="tab:blue", zorder=1)
ax[2].fill_between(tg01, δτg01[:,l-1] - 2δτerrg01[:,l-1], δτg01[:,l-1] + 2δτerrg01[:,l-1], alpha=.25,
                   color="tab:blue", linewidths=0, zorder=1)
ax[1].plot(tg01, τargog01[:,1], color="tab:orange", linewidth=1, zorder=0, label="ECCO")
ax[2].plot(tg01, τargog01[:,1] - τargog01[:,l], color="tab:orange", linewidth=1, zorder=0,
           label="Argo")
ax[1].plot(tg01, τeccog01[:,1], color="tab:green", linewidth=1, zorder=0, label="ECCO")
ax[2].plot(tg01, τeccog01[:,1] - τeccog01[:,l], color="tab:green", linewidth=1, zorder=0,
           label="ECCO")
ax[1].plot(t01, cτargo01[1]*tr01 .+ cτargo01[2], color="tab:orange", linewidth=1, zorder=0)
ax[1].plot(t01, cτecco01[1]*tr01 .+ cτecco01[2], color="tab:green", linewidth=1, zorder=0)
ax[1].plot(t01, cτ01[1]*tr01 .+ cτ01[2], color="tab:blue", linewidth=1, zorder=1)
ax[2].plot(t01, cδτargo01[1]*tr01 .+ cδτargo01[2], color="tab:orange", linewidth=1, zorder=0)
ax[2].plot(t01, cδτecco01[1]*tr01 .+ cδτecco01[2], color="tab:green", linewidth=1, zorder=0)
ax[2].plot(t01, cδτ01[1]*tr01 .+ cδτ01[2], color="tab:blue", linewidth=1, zorder=1)
idx = findall(isnan.(τg01[:,1]))
for i = 1 : length(idx)+1
  i0 = i == 1 ? 1 : idx[i-1]+1
  i1 = i == length(idx)+1 ? length(t01) : idx[i]-1
  ti = tg01[i0:i1-1] + .5(tg01[i0+1:i1] - tg01[i0:i1-1])
  ti = [tg01[i0] - (ti[1]-tg01[i0]); ti; tg01[i1] + (tg01[i1]-ti[end])]
  ax[3].pcolormesh(ti, depthi, T01[:,i0:i1], cmap="RdBu_r", vmin=-.17, vmax=.17)
end
ax[1].invert_yaxis()
ax[1].legend(frameon=false, loc=4, ncol=3)
ax[1].set_xlim(t01[1] - (t01[end] - t01[1])÷100, t01[end] + (t01[end] - t01[1])÷100)
ax[3].set_ylim(5, 0)
ax[1].set_ylabel("travel time anomaly (s)")
ax[2].set_ylabel("travel time difference (s)")
ax[3].set_ylabel("depth (km)")
ax[1].set_yticks(-.8:.4:.8)
axr = ax[1].twinx()
axr.set_yticks(-.1:.05:.1)
axr.set_ylim(ax[1].get_ylim()./sum(K01[1,:]))
axr.set_ylabel("temperature anomaly (K)")
fig.align_ylabels()
fig.tight_layout()

# H08 depth–time plots
fig, ax = subplots(3, 1, figsize=(190/25.4, 190/25.4), sharex=true, sharey=true)
idx = findall(isnan.(τg08[:,1]))
for i = 1 : length(idx)+1
  i0 = i == 1 ? 1 : idx[i-1]+1
  i1 = i == length(idx)+1 ? length(t08) : idx[i]-1
  ti = tg08[i0:i1-1] + .5(tg08[i0+1:i1] - tg08[i0:i1-1])
  ti = [tg08[i0] - (ti[1]-tg08[i0]); ti; tg08[i1] + (tg08[i1]-ti[end])]
  ax[1].pcolormesh(ti, depthi, T08[:,i0:i1], cmap="RdBu_r", vmin=-.11, vmax=.11)
  ax[2].pcolormesh(ti, depthi, Targoi08[:,i0:i1], cmap="RdBu_r", vmin=-.11, vmax=.11)
  ax[3].pcolormesh(ti, depthi, Teccoi08[:,i0:i1], cmap="RdBu_r", vmin=-.11, vmax=.11)
end
ax[1].set_ylim(5, 0)
ax[1].set_ylabel("depth (km)")
ax[2].set_ylabel("depth (km)")
ax[3].set_ylabel("depth (km)")
fig.tight_layout()

# H01 depth–time plots
fig, ax = subplots(3, 1, figsize=(190/25.4, 190/25.4), sharex=true, sharey=true)
idx = findall(isnan.(τg01[:,1]))
for i = 1 : length(idx)+1
  i0 = i == 1 ? 1 : idx[i-1]+1
  i1 = i == length(idx)+1 ? length(t01) : idx[i]-1
  ti = tg01[i0:i1-1] + .5(tg01[i0+1:i1] - tg01[i0:i1-1])
  ti = [tg01[i0] - (ti[1]-tg01[i0]); ti; tg01[i1] + (tg01[i1]-ti[end])]
  ax[1].pcolormesh(ti, depthi, T01[:,i0:i1], cmap="RdBu_r", vmin=-.17, vmax=.17)
  ax[2].pcolormesh(ti, depthi, Targoi01[:,i0:i1], cmap="RdBu_r", vmin=-.17, vmax=.17)
  ax[3].pcolormesh(ti, depthi, Teccoi01[:,i0:i1], cmap="RdBu_r", vmin=-.17, vmax=.17)
end
ax[1].set_ylim(5, 0)
ax[1].set_ylabel("depth (km)")
ax[2].set_ylabel("depth (km)")
ax[3].set_ylabel("depth (km)")
fig.tight_layout()

# plot kernels, singular vectors, resolution matrix
fig, ax = subplots(2, 3, sharex="col", sharey=true, figsize=(190/25.4, 150/25.4),
                   gridspec_kw=Dict("width_ratios"=>(1, 1, 1.9)))
ax[1,1].plot(1e3K08[1,:]/Δz, depth, linewidth=1, label="2.00 Hz")
ax[1,1].plot(1e3K08[2,:]/Δz, depth, linewidth=1, label="3.00 Hz")
ax[1,1].plot(1e3K08[3,:]/Δz, depth, linewidth=1, label="4.00 Hz")
ax[1,2].plot(V08[:,1], depth, linewidth=1, color="tab:olive")
ax[1,2].plot(V08[:,2], depth, linewidth=1, color="tab:purple")
ax[1,2].axvline(0, color="black", linewidth=.8)
img = ax[1,3].imshow(V08[:,1:2]*V08[:,1:2]', cmap="RdBu_r", origin="lower", vmin=-.07,
                   vmax=.07, extent=[depthi[1], depthi[end], depthi[1], depthi[end]])
cb = colorbar(img, ax=ax[1,3])
#cb.outline.set_visible(false)
ax[2,1].plot(1e3K01[1,:]/Δz, depth, linewidth=1, label="2.50 Hz")
ax[2,1].plot(1e3K01[2,:]/Δz, depth, linewidth=1, label="3.25 Hz")
ax[2,1].plot(1e3K01[3,:]/Δz, depth, linewidth=1, label="4.00 Hz")
ax[2,2].plot(V01[:,1], depth, linewidth=1, color="tab:olive")
ax[2,2].plot(V01[:,2], depth, linewidth=1, color="tab:purple")
ax[2,2].axvline(0, color="black", linewidth=.8)
img = ax[2,3].imshow(V01[:,1:2]*V01[:,1:2]', cmap="RdBu_r", origin="lower", vmin=-.07,
                   vmax=.07, extent=[depthi[1], depthi[end], depthi[1], depthi[end]])
cb = colorbar(img, ax=ax[2,3])
#cb.outline.set_visible(false)
ax[1,1].set_ylim(5, 0)
ax[1,1].set_ylabel("depth (km)")
ax[2,1].set_ylabel("depth (km)")
ax[1,1].set_xlim(-8, 0)
ax[1,2].set_xlim(-.25, .25)
ax[1,3].set_xlim(0, 5)
ax[2,1].set_xlabel(L"sensitivity (s$\,$K$^{-1}\,$km$^{-1}$)")
ax[2,3].set_xlabel("depth (km)")
ax[1,1].legend(frameon=false, loc=3)
ax[2,1].legend(frameon=false, loc=3)
ax[1,2].legend([L"$v_1$", L"$v_2$"], frameon=false)
ax[2,2].legend([L"$v_1$", L"$v_2$"], frameon=false)
ax[1,1].set_title("H08 sensitivities")
ax[1,2].set_title("H08 singular vectors")
ax[1,3].set_title("H08 resolution matrix")
ax[2,1].set_title("H01 sensitivities")
ax[2,2].set_title("H01 singular vectors")
ax[2,3].set_title("H01 resolution matrix")
fig.tight_layout()

# plot trend profiles
fig, ax = subplots(2, 2, sharex=true, sharey=true, figsize=(95/25.4, 150/25.4))
ax[1,1].plot(NaN*depth, depth, color="tab:blue", linewidth=1, label=L"$T$ wave")
ax[1,1].plot(1e3Ttrendargof08*SOT.meanyear, depth, linewidth=1, color="tab:orange",
             label="Argo")
ax[1,1].plot(1e3Ttrendeccof08*SOT.meanyear, depth, linewidth=1, color="tab:green",
             label="ECCO")
ax[1,1].plot(1e3Ttrendctdf[depth.≥2], depth[depth.≥2], linewidth=1, color="tab:red",
             zorder=1, label="Hydr.")
ax[1,2].plot(1e3Ttrendtwave08*SOT.meanyear, depth, zorder=2, linewidth=1, color="tab:blue")
ax[1,2].plot(1e3Ttrendargo08*SOT.meanyear, depth, zorder=1, linewidth=1, color="tab:orange")
ax[1,2].plot(1e3Ttrendecco08*SOT.meanyear, depth, zorder=1, linewidth=1, color="tab:green")
ax[2,1].plot(1e3Ttrendargof01*SOT.meanyear, depth, linewidth=1, color="tab:orange")
ax[2,1].plot(1e3Ttrendeccof01*SOT.meanyear, depth, linewidth=1, color="tab:green")
ax[2,1].plot(1e3Ttrendctdf[depth.≥2], depth[depth.≥2], linewidth=1, color="tab:red",
             zorder=1)
ax[2,2].plot(1e3Ttrendtwave01*SOT.meanyear, depth, zorder=2, linewidth=1, color="tab:blue")
ax[2,2].plot(1e3Ttrendargo01*SOT.meanyear, depth, zorder=1, linewidth=1, color="tab:orange")
ax[2,2].plot(1e3Ttrendecco01*SOT.meanyear, depth, zorder=1, linewidth=1, color="tab:green")
ax[1].axvline(0, color="black", linewidth=.8, zorder=0)
ax[2].axvline(0, color="black", linewidth=.8, zorder=0)
ax[3].axvline(0, color="black", linewidth=.8, zorder=0)
ax[4].axvline(0, color="black", linewidth=.8, zorder=0)
ax[1].axhline(2, color="lightgray", linewidth=.8, zorder=0)
ax[2].axhline(2, color="lightgray", linewidth=.8, zorder=0)
ax[3].axhline(2, color="lightgray", linewidth=.8, zorder=0)
ax[4].axhline(2, color="lightgray", linewidth=.8, zorder=0)
ax[1].set_xlim(-2, 8)
ax[1].set_ylim(5, 0)
ax[1,1].set_ylabel("depth (km)")
ax[2,1].set_ylabel("depth (km)")
ax[2,1].set_xlabel(L"trend (mK$\,$yr$^{-1}$)")
ax[2,2].set_xlabel(L"trend (mK$\,$yr$^{-1}$)")
ax[1,1].set_title("H08 trends")
ax[1,2].set_title("H08 projected trends")
ax[2,1].set_title("H01 trends")
ax[2,2].set_title("H01 projected trends")
ax[1,1].legend(frameon=false, loc=4)
fig.tight_layout(w_pad=1.5)

using .SOT
using PyPlot, Printf, Dates
using Interpolations, Statistics
using HDF5, DelimitedFiles, NPZ
using LinearAlgebra, SparseArrays
using PyCall

# small font
rc("font", size=10)
rcParams = PyDict(matplotlib["rcParams"])
rcParams["font.sans-serif"] = "Fira Sans"
rcParams["mathtext.fontset"] = "custom"
rcParams["mathtext.sf"] = "Fira Sans"

# excluded time periods for H08: before 2004-12-01 and periods with uncorrected clock error
excludetimes08 = [[Date(2001, 1, 1) Date(2004, 12, 1)],
                  [DateTime("2010-03-16T00:00:00") DateTime("2010-05-17T02:06:17.760")],
                  [Date(2017, 6, 1) Date(2018, 1, 1)]]

# excluded time period for DGAR: before 2004-12-01
excludetimesdgar = [[Date(2001, 1, 1) Date(2004, 12, 1)]]

# number of frequencies
l01 = 3
l08 = 3
ldgar = 2

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

# read DGAR data
h5open("results/nias_DGAR.h5", "r") do file
  global trdgar = read(file, "t")
  global τdgar = read(file, "tau")
  global τerrdgar = read(file, "tauerr")
  global δτdgar = read(file, "dtau")
  global δτerrdgar = read(file, "dtauerr")
  global τargodgar = read(file, "tauargo")
  global τeccodgar = read(file, "tauecco")
end
tdgar = DateTime(2000, 1, 1, 0, 0, 0) .+ Millisecond.(Int.(round.(trdgar*1000*3600*24)))

# offset H01 based on trends
tr01 = Dates.value.(t01 .- t01[1])/1000/3600/24
for i = 1:l01
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

# offset H08 based on trends
tr08 = Dates.value.(t08 .- t08[1])/1000/3600/24
for i = 1:l08
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

# offset DGAR based on trends
trdgar = Dates.value.(tdgar .- tdgar[1])/1000/3600/24
for i = 1:ldgar
  c, _ = SOT.lineartrend(tdgar, τdgar[:,i]; fitannual=true, fitsemiannual=true)
  τdgar[:,i] .-= c[2] + c[1]*trdgar[end]/2
  if i == 1
    δτdgar[:,:] .-= c[2] + c[1]*trdgar[end]/2
  else
    δτdgar[:,i-1] .+= c[2] + c[1]*trdgar[end]/2
  end
  c, _ = SOT.lineartrend(tdgar, τeccodgar[:,i]; fitannual=true, fitsemiannual=true)
  τeccodgar[:,i] .-= c[2] + c[1]*trdgar[end]/2
  c, _ = SOT.lineartrend(tdgar, τargodgar[:,i]; fitannual=true, fitsemiannual=true)
  τargodgar[:,i] .-= c[2] + c[1]*trdgar[end]/2
end

# estimate H01 trends
cτ01, _ = SOT.lineartrend(t01, τ01[:,1]; fitannual=true, fitsemiannual=true)
cδτ01, _ = SOT.lineartrend(t01, δτ01[:,l01-1]; fitannual=true, fitsemiannual=true)
cτecco01, _ = SOT.lineartrend(t01, τecco01[:,1]; fitannual=true, fitsemiannual=true)
cδτecco01, _ = SOT.lineartrend(t01, τecco01[:,1] - τecco01[:,l01]; fitannual=true, fitsemiannual=true)
cτargo01, _ = SOT.lineartrend(t01, τargo01[:,1]; fitannual=true, fitsemiannual=true)
cδτargo01, _ = SOT.lineartrend(t01, τargo01[:,1] - τargo01[:,2]; fitannual=true, fitsemiannual=true)

# estimate H08 trends
cτ08, _ = SOT.lineartrend(t08, τ08[:,1]; fitannual=true, fitsemiannual=true)
cδτ08, _ = SOT.lineartrend(t08, δτ08[:,l08-1]; fitannual=true, fitsemiannual=true)
cτecco08, _ = SOT.lineartrend(t08, τecco08[:,1]; fitannual=true, fitsemiannual=true)
cδτecco08, _ = SOT.lineartrend(t08, τecco08[:,1] - τecco08[:,l08]; fitannual=true, fitsemiannual=true)
cτargo08, _ = SOT.lineartrend(t08, τargo08[:,1]; fitannual=true, fitsemiannual=true)
cδτargo08, _ = SOT.lineartrend(t08, τargo08[:,1] - τargo08[:,2]; fitannual=true, fitsemiannual=true)

# estimate H01 trends
cτdgar, _ = SOT.lineartrend(tdgar, τdgar[:,1]; fitannual=true, fitsemiannual=true)
cδτdgar, _ = SOT.lineartrend(tdgar, δτdgar[:,ldgar-1]; fitannual=true, fitsemiannual=true)
cτeccodgar, _ = SOT.lineartrend(tdgar, τeccodgar[:,1]; fitannual=true, fitsemiannual=true)
cδτeccodgar, _ = SOT.lineartrend(tdgar, τeccodgar[:,1] - τeccodgar[:,2]; fitannual=true, fitsemiannual=true)
cτargodgar, _ = SOT.lineartrend(tdgar, τargodgar[:,1]; fitannual=true, fitsemiannual=true)
cδτargodgar, _ = SOT.lineartrend(tdgar, τargodgar[:,1] - τargodgar[:,2]; fitannual=true, fitsemiannual=true)

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
tg01, τg01, τerrg01, δτg01, δτerrg01, τargog01, τeccog01 = insertgap(t01, [τ01, τerr01, δτ01, δτerr01, τargo01, τecco01], Date(2011, 10, 27))
tg08, τg08, τerrg08, δτg08, δτerrg08, τargog08, τeccog08 = insertgap(t08, [τ08, τerr08, δτ08, δτerr08, τargo08, τecco08], Date(2007, 1, 1))

tgdgar = copy(tdgar)
τgdgar = copy(τdgar)
τerrgdgar = copy(τerrdgar)
δτgdgar = copy(δτdgar)
δτerrgdgar = copy(δτerrdgar)
τargogdgar = copy(τargodgar)
τeccogdgar = copy(τeccodgar)

# excluded periods
for i = 1:length(excludetimes08)
  global tg08, τg08, τerrg08, δτg08, δτerrg08, τargog08, τeccog08
  tg08, τg08, τerrg08, δτg08, δτerrg08, τargog08, τeccog08 = insertgap(tg08, [τg08, τerrg08, δτg08, δτerrg08, τargog08, τeccog08], excludetimes08[i][1])
end
for i = 1:length(excludetimesdgar)
  global tgdgar, τgdgar, τerrgdgar, δτgdgar, δτerrgdgar, τargogdgar, τeccogdgar
  tgdgar, τgdgar, τerrgdgar, δτgdgar, δτerrgdgar, τargogdgar, τeccogdgar = insertgap(tgdgar, [τgdgar, τerrgdgar, δτgdgar, δτerrgdgar, τargogdgar, τeccogdgar], excludetimesdgar[i][1])
end

# seismic station gap
tg01, τg01, τerrg01, δτg01, δτerrg01, τargog01, τeccog01 = insertgap(tg01, [τg01, τerrg01, δτg01, δτerrg01, τargog01, τeccog01], Date(2009, 7, 1))
tg08, τg08, τerrg08, δτg08, δτerrg08, τargog08, τeccog08 = insertgap(tg08, [τg08, τerrg08, δτg08, δτerrg08, τargog08, τeccog08], Date(2009, 7, 1))

tgdgar, τgdgar, τerrgdgar, δτgdgar, δτerrgdgar, τargogdgar, τeccogdgar = insertgap(tgdgar, [τgdgar, τerrgdgar, δτgdgar, δτerrgdgar, τargogdgar, τeccogdgar], Date(2007, 4, 1))
tgdgar, τgdgar, τerrgdgar, δτgdgar, δτerrgdgar, τargogdgar, τeccogdgar = insertgap(tgdgar, [τgdgar, τerrgdgar, δτgdgar, δτerrgdgar, τargogdgar, τeccogdgar], Date(2007, 8, 1))
tgdgar, τgdgar, τerrgdgar, δτgdgar, δτerrgdgar, τargogdgar, τeccogdgar = insertgap(tgdgar, [τgdgar, τerrgdgar, δτgdgar, δτerrgdgar, τargogdgar, τeccogdgar], Date(2007, 11, 1))

Δx = 5e3
Δz = 50.
c0 = 1.5e3
α = 4.6
depth = 6.5:-1e-3Δz:0
depthi = [depth; -1e-3Δz] .+ 1e-3Δz/2

K01 = h5read("data/kernels/nias_H01.h5", "K")
replace!(K01, NaN=>0)
K01 = sum(K01; dims=2)[:,1,:]'*Δx*Δz

K08 = h5read("data/kernels/nias_H08.h5", "K")
replace!(K08, NaN=>0)
K08 = sum(K08; dims=2)[:,1,:]'*Δx*Δz

# SVD
U01, Σ01, V01 = svd(K01)
U08, Σ08, V08 = svd(K08)

# infer temperature profiles
T01 = V01[:,1:2]*diagm(Σ01[1:2].^-1)*U01[:,1:2]'*τ01'
T08 = V08[:,1:2]*diagm(Σ08[1:2].^-1)*U08[:,1:2]'*τ08'

fig, ax = subplots(1, 2, sharey=true)
colors = PyPlot.cm.inferno([.15, .5, .85])
for i = 1:l08
  ax[1].plot(K08[i,:], depth, color=colors[i,:])
end
ax[2].plot(V08[:,1:2], depth)
ax[1].legend(["2 Hz", "3 Hz", "4 Hz"], frameon=false, loc=3)
ax[2].legend([L"$v_1$", L"$v_2$"], frameon=false, loc=4)
ax[1].set_xlim(-.2, 0)
ax[2].set_xlim(-.25, .25)
ax[1].set_ylim(5, 0)
ax[1].set_ylabel("depth (km)")
ax[1].set_xlabel(L"sensitivity (s$\,$K$^{-1}\,$km$^{-1}$)")
ax[1].set_title("kernels")
ax[2].set_title("singular vectors")
fig.tight_layout(w_pad=2)

fig, ax = subplots(1, 1, figsize=(10, 4))
ax.plot(tgdgar, τgdgar[:,1], zorder=1, label="DGAR")
ax.plot(tg08, τg08[:,1], color="tab:red", zorder=2, label="CTBTO H08")
ax.scatter(tgdgar, τgdgar[:,1], s=3, c="tab:blue", zorder=1)
ax.scatter(tg08, τg08[:,1], s=3, c="tab:red", zorder=2)
ax.fill_between(tgdgar, τgdgar[:,1] - 2τerrgdgar[:,1], τgdgar[:,1] + 2τerrgdgar[:,1], alpha=.25, color="tab:blue", linewidths=0, zorder=1)
ax.fill_between(tg08, τg08[:,1] - 2τerrg08[:,1], τg08[:,1] + 2τerrg08[:,1], alpha=.25, color="tab:red", linewidths=0, zorder=2)
ax.set_xlim(t08[1] - (t08[end] - t08[1])÷100, t08[end] + (t08[end] - t08[1])÷100)
ax.legend(frameon=false, ncol=2, loc=4)
ax.set_ylim(-.09*sum(K08[1,:]), .09*sum(K08[1,:]))
ax.set_ylabel("travel time anomaly (s)")
ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y"))
ax.set_xticks(Date(2005, 1, 1) : Year(2) : Date(2017, 1, 1))
ax.set_xticks(Date(2006, 1, 1) : Year(2) : Date(2016, 1, 1), "minor")
axr = ax.twinx()
axr.set_yticks(-90:30:90)
axr.set_ylim(1e3.*ax.get_ylim()./sum(K08[1,:]))
axr.set_ylabel("temperature anomaly (mK)")
fig.tight_layout()

fig, ax = subplots(2, 1, figsize=(10, 5.75), sharex=true)
ax[1].plot(tg08, τg08[:,1], color="tab:blue", zorder=2, label="CTBTO H08")
ax[1].scatter(tg08, τg08[:,1], s=3, c="tab:blue", zorder=2)
ax[1].fill_between(tg08, τg08[:,1] - 2τerrg08[:,1], τg08[:,1] + 2τerrg08[:,1], alpha=.25, color="tab:blue", linewidths=0, zorder=2)
ax[2].plot(tg08, δτg08[:,l08-1], color="tab:blue", zorder=2, label="CTBTO H08")
ax[2].scatter(tg08, δτg08[:,l08-1], s=3, c="tab:blue", zorder=2)
ax[2].fill_between(tg08, δτg08[:,l08-1] - 2δτerrg08[:,l08-1], δτg08[:,l08-1] + 2δτerrg08[:,l08-1], alpha=.25, color="tab:blue", linewidths=0, zorder=2)
ax[1].set_xlim(t08[1] - (t08[end] - t08[1])÷100, t08[end] + (t08[end] - t08[1])÷100)
ax[1].set_ylim(-.09*sum(K08[1,:]), .09*sum(K08[1,:]))
ax[2].set_ylim(-.132, .132)
ax[1].set_ylabel("travel time anomaly (s)")
ax[2].set_ylabel("travel time difference (s)")
ax[1].text(0.015, 0.95, "2 Hz", va="top", transform=ax[1].transAxes)
ax[2].text(0.015, 0.95, "4 Hz – 2 Hz", va="top", transform=ax[2].transAxes)
ax[1].xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y"))
ax[1].set_xticks(Date(2005, 1, 1) : Year(2) : Date(2017, 1, 1))
ax[1].set_xticks(Date(2006, 1, 1) : Year(2) : Date(2016, 1, 1), "minor")
axr = ax[1].twinx()
axr.set_yticks(-90:30:90)
axr.set_ylim(1e3.*ax[1].get_ylim()./sum(K08[1,:]))
axr.set_ylabel("temperature anomaly (mK)")
fig.align_ylabels()
fig.tight_layout()

fig, ax = subplots(2, 1, figsize=(10, 5.75), sharex=true)
ax[1].plot(tg08, τg08[:,1], color="tab:blue", zorder=2, label="CTBTO H08")
ax[1].plot(tg08, τargog08[:,1], color="tab:orange", zorder=1, label="Argo (SIO)")
ax[1].plot(tg08, τeccog08[:,1], color="tab:green", zorder=1, label="ECCO")
ax[1].scatter(tg08, τg08[:,1], s=3, c="tab:blue", zorder=2)
ax[1].fill_between(tg08, τg08[:,1] - 2τerrg08[:,1], τg08[:,1] + 2τerrg08[:,1], alpha=.25, color="tab:blue", linewidths=0, zorder=2)
ax[2].plot(tg08, δτg08[:,l08-1], color="tab:blue", zorder=2, label="CTBTO H08")
ax[2].plot(tg08, τargog08[:,1] - τargog08[:,l08], color="tab:orange", zorder=1, label="Argo (SIO)")
ax[2].plot(tg08, τeccog08[:,1] - τeccog08[:,l08], color="tab:green", zorder=1, label="ECCO")
ax[2].scatter(tg08, δτg08[:,l08-1], s=3, c="tab:blue", zorder=2)
ax[2].fill_between(tg08, δτg08[:,l08-1] - 2δτerrg08[:,l08-1], δτg08[:,l08-1] + 2δτerrg08[:,l08-1], alpha=.25, color="tab:blue", linewidths=0, zorder=2)
ax[1].set_xlim(t08[1] - (t08[end] - t08[1])÷100, t08[end] + (t08[end] - t08[1])÷100)
ax[1].legend(frameon=false, ncol=3, loc=4)
ax[1].set_ylim(-.09*sum(K08[1,:]), .09*sum(K08[1,:]))
ax[2].set_ylim(-.132, .132)
ax[1].set_ylabel("travel time anomaly (s)")
ax[2].set_ylabel("travel time difference (s)")
ax[1].text(0.015, 0.95, "2 Hz", va="top", transform=ax[1].transAxes)
ax[2].text(0.015, 0.95, "4 Hz – 2 Hz", va="top", transform=ax[2].transAxes)
ax[1].xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y"))
ax[1].set_xticks(Date(2005, 1, 1) : Year(2) : Date(2017, 1, 1))
ax[1].set_xticks(Date(2006, 1, 1) : Year(2) : Date(2016, 1, 1), "minor")
axr = ax[1].twinx()
axr.set_yticks(-90:30:90)
axr.set_ylim(1e3.*ax[1].get_ylim()./sum(K08[1,:]))
axr.set_ylabel("temperature anomaly (mK)")
fig.align_ylabels()
fig.tight_layout()

fig, ax = subplots(1, 1, figsize=(10, 4))
idx = findall(isnan.(τg08[:,1]))
for i = 1 : length(idx)+1
  i0 = i == 1 ? 1 : idx[i-1]+1
  i1 = i == length(idx)+1 ? length(t08) : idx[i]-1
  ti = tg08[i0:i1-1] + .5(tg08[i0+1:i1] - tg08[i0:i1-1])
  ti = [tg08[i0] - (ti[1]-tg08[i0]); ti; tg08[i1] + (tg08[i1]-ti[end])]
  img = ax.pcolormesh(ti, depthi, 1e3*T08[:,i0:i1], cmap="RdBu_r", vmin=-120, vmax=120, rasterized=true)
  if i == 1
    colorbar(img, label="temperature anomaly (mK)", fraction=.05, pad=.025)
  end
end
ax.set_xlim(t08[1] - (t08[end] - t08[1])÷100, t08[end] + (t08[end] - t08[1])÷100)
ax.legend(frameon=false, ncol=2, loc=4)
ax.set_ylim(5, 0)
ax.set_ylabel("depth (km)")
ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y"))
ax.set_xticks(Date(2005, 1, 1) : Year(2) : Date(2017, 1, 1))
ax.set_xticks(Date(2006, 1, 1) : Year(2) : Date(2016, 1, 1), "minor")
fig.tight_layout()

fig, ax = subplots(1, 1, figsize=(10, 4))
idx = findall(isnan.(τg08[:,1]))
for i = 1 : length(idx)+1
  i0 = i == 1 ? 1 : idx[i-1]+1
  i1 = i == length(idx)+1 ? length(t08) : idx[i]-1
  ti = tg08[i0:i1-1] + .5(tg08[i0+1:i1] - tg08[i0:i1-1])
  ti = [tg08[i0] - (ti[1]-tg08[i0]); ti; tg08[i1] + (tg08[i1]-ti[end])]
  img = ax.pcolormesh(ti, depthi, 1e3*T08[:,i0:i1], cmap="RdBu_r", vmin=-120, vmax=120, rasterized=true)
  if i == 1
    colorbar(img, label="temperature anomaly (mK)", fraction=.05, pad=.025)
  end
end
ax.set_xlim(Date(2005, 1, 1), Date(2006, 10, 1))
ax.legend(frameon=false, ncol=2, loc=4)
ax.set_ylim(5, 0)
ax.set_ylabel("depth (km)")
ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y-%m"))
ax.set_xticks(Date(2005, 1, 1) : Month(3) : Date(2006, 10, 1))
ax.set_xticks(Date(2005, 1, 1) : Month(1) : Date(2006, 10, 1), "minor")
fig.tight_layout()

fig, ax = subplots(1, 1, figsize=(10, 4))
idx = findall(isnan.(τg01[:,1]))
for i = 1 : length(idx)+1
  i0 = i == 1 ? 1 : idx[i-1]+1
  i1 = i == length(idx)+1 ? length(t01) : idx[i]-1
  ti = tg01[i0:i1-1] + .5(tg01[i0+1:i1] - tg01[i0:i1-1])
  ti = [tg01[i0] - (ti[1]-tg01[i0]); ti; tg01[i1] + (tg01[i1]-ti[end])]
  img = ax.pcolormesh(ti, depthi, 1e3*T01[:,i0:i1], cmap="RdBu_r", vmin=-170, vmax=170, rasterized=true)
  if i == 1
    colorbar(img, label="temperature anomaly (mK)", fraction=.05, pad=.025)
  end
end
ax.set_xlim(t01[1] - (t01[end] - t01[1])÷100, t01[end] + (t01[end] - t01[1])÷100)
ax.legend(frameon=false, ncol=2, loc=4)
ax.set_ylim(5, 0)
ax.set_ylabel("depth (km)")
ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y"))
ax.set_xticks(Date(2005, 1, 1) : Year(2) : Date(2017, 1, 1))
ax.set_xticks(Date(2006, 1, 1) : Year(2) : Date(2016, 1, 1), "minor")
fig.tight_layout()

fig, ax = subplots(1, 1, figsize=(10, 4))
ax.plot(tg01, τg01[:,1], color="tab:blue", zorder=2, label="CTBTO H01")
ax.scatter(tg01, τg01[:,1], s=3, c="tab:blue", zorder=2)
ax.fill_between(tg01, τg01[:,1] - 2τerrg01[:,1], τg01[:,1] + 2τerrg01[:,1], alpha=.25, color="tab:blue", linewidths=0, zorder=2)
ax.set_xlim(t01[1] - (t01[end] - t01[1])÷100, t01[end] + (t01[end] - t01[1])÷100)
ax.set_yticks(-1:.5:1)
ax.set_ylim(1, -1)
ax.set_ylabel("travel time anomaly (s)")
ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y"))
ax.set_xticks(Date(2005, 1, 1) : Year(2) : Date(2017, 1, 1))
ax.set_xticks(Date(2006, 1, 1) : Year(2) : Date(2016, 1, 1), "minor")
axr = ax.twinx()
axr.set_yticks(-120:40:120)
axr.set_ylim(1e3.*ax.get_ylim()./sum(K01[1,:]))
axr.set_ylabel("temperature anomaly (mK)")
fig.tight_layout()

fig, ax = subplots(1, 1, figsize=(10, 4))
ax.plot(tg01, τg01[:,1], color="tab:blue", zorder=2, label="CTBTO H01")
ax.plot(tg01, τargog01[:,1], color="tab:orange", zorder=1, label="Argo (SIO)")
ax.plot(tg01, τeccog01[:,1], color="tab:green", zorder=1, label="ECCO")
ax.scatter(tg01, τg01[:,1], s=3, c="tab:blue", zorder=2)
ax.fill_between(tg01, τg01[:,1] - 2τerrg01[:,1], τg01[:,1] + 2τerrg01[:,1], alpha=.25, color="tab:blue", linewidths=0, zorder=2)
ax.set_xlim(t01[1] - (t01[end] - t01[1])÷100, t01[end] + (t01[end] - t01[1])÷100)
ax.legend(frameon=false, ncol=3, loc=4)
ax.set_yticks(-1:.5:1)
ax.set_ylim(1, -1)
ax.set_ylabel("travel time anomaly (s)")
ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y"))
ax.set_xticks(Date(2005, 1, 1) : Year(2) : Date(2017, 1, 1))
ax.set_xticks(Date(2006, 1, 1) : Year(2) : Date(2016, 1, 1), "minor")
axr = ax.twinx()
axr.set_yticks(-120:40:120)
axr.set_ylim(1e3.*ax.get_ylim()./sum(K01[1,:]))
axr.set_ylabel("temperature anomaly (K)")
fig.tight_layout()

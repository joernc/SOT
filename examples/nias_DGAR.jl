using .SOT, PyPlot, Printf, Dates, LinearAlgebra, Statistics, SparseArrays
using HDF5, Interpolations, DataFrames, CSV

# identifier for experiment
eqname = "nias"

# P-wave (reference) stations
pstations = ["PS.PSI..BHZ", "MY.KUM..BHZ", "II.WRAB.00.BHZ", "GE.GSI..BHZ"]

# intervals to which to cut P waveforms
pintervals = [[-3, 47], [-3, 47], [-3, 47], [-3, 47]]

# frequency bands to which to filter P waveforms
pfreqbands = [[1, 3], [1, 3], [1.5, 2.5], [1, 3]]

# T-wave station
tstations = ["II.DGAR.00.BHZ"]

# T-wave time window around predicted arrival time
tintervals = [[-10, 70]]

# T-wave filtering window width
tavgwidth = 0.5

# T-wave reference frequency at which to find first max CC
treffreq = 2.0

# frequencies used in inversion
tinvfreq = [2.0, 3.0]

# minimum CCs for T-wave pairs (at inversion frequencies)
tmincc = [0.6, 0.5]

# excluded time period: before 2004-12-01
excludetimes = [[Date(2001, 1, 1) Date(2004, 12, 1)]]

# manually exclude pairs
excludepairs = CSV.read("data/catalogs/nias_DGAR_exclude.csv", DataFrame)

# download P-wave data
SOT.downloadseisdata(eqname, pstations)

# cut and filter P waveforms
SOT.cutpwaves(eqname, pstations, pintervals, pfreqbands)

# find P-wave pairs
SOT.findpairs(eqname, pstations, pintervals, pfreqbands, saveplot=true)

# download T-wave data
SOT.downloadseisdata(eqname, tstations)

# measure T-wave lags Δτ
SOT.twavepick(eqname, tstations, tintervals, tavgwidth, treffreq, pstations, pintervals,
              pfreqbands, saveplot=true)

# collect usable pairs
tpairs, ppairs = SOT.collectpairs(eqname, tstations, tintervals, tavgwidth, treffreq,
                                  tinvfreq, tmincc, pstations, pintervals, pfreqbands;
                                  excludetimes, excludepairs)

# DGAR clock error correction
breaktime = DateTime(2012, 3, 17)
idx = (tpairs.event1 .< breaktime) .& (tpairs.event2 .> breaktime)
for i = findall(idx)
  tpairs.Δτl[i] .-= 1.0
  tpairs.Δτc[i] .-= 1.0
  tpairs.Δτr[i] .-= 1.0
end

# perform inversion
t, E, S, P, D = SOT.invert(tpairs, ppairs)

# number of good T- and P-wave pairs
nt = size(tpairs, 1)
np = size(ppairs, 1)

# number of unique events
m = length(t)

# number of frequencies
l = length(tinvfreq)

# invert for P-wave delays without smoothing
Δτp = E[1:nt,1:m]*(E[l*nt+1:l*nt+np,l*m+1:(l+1)*m]\ppairs.Δτ)

# make cycle-skipping correction
tpairs.Δτ = SOT.correctcycleskipping(tpairs, ppairs, E, S, P)

# data vector
y = [reshape(vcat(tpairs.Δτ'...), l*nt); ppairs.Δτ]

# get least-squares solution
x = P*(E'*y)

# variance estimate
n = y - E*x
σt = sqrt(mean(n[1:l*nt].^2)*(l*nt+np)/(l*nt+np-l*m-1))
σp = sqrt(mean(n[l*nt+1:l*nt+np].^2)*(l*nt+np)/(l*nt+np-l*m-1))
R = spdiagm(0 => [σt^2*ones(l*nt); σp^2*ones(np)])

# matrix for error propagation
A = P*E'

# get travel time anomalies and their errors
τ = reshape(D*x, (m, l))
τerr = reshape(sqrt.(diag(D*A*R*A'*D')), (m, l))

# get frequency differences and their errors
δD = [vcat([I(m) for i = 1:l-1]...) -I((l-1)*m) spzeros((l-1)*m, m)]
δτ = reshape(δD*x, (m, l-1))
δτerr = reshape(sqrt.(diag(δD*A*R*A'*δD')), (m, l-1))

# read and interpolate Argo data
targo = Dates.value.(Date(2004, 1, 15) .+ Month.(h5read("data/argo/nias_H08.h5", "time"))
                     .- Date(2000, 1, 1))
tr = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
τargo = hcat([extrapolate(interpolate((targo,), h5read("data/argo/nias_H08.h5", "tau")[i,:],
                                      Gridded(Linear())), Flat())(tr) for i = 1:l]...)

# read and interpolate ECCO data
tecco = h5read("data/ecco/nias_H08.h5", "time")
tr = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
τecco = hcat([extrapolate(interpolate((tecco,), h5read("data/ecco/nias_H08.h5", "tau")[i,:],
                                      Gridded(Linear())), Flat())(tr) for i = 1:2]...)

# save timeseries to file
tr = Dates.value.(t - DateTime(2000, 1, 1, 0, 0, 0))/1000/3600/24
h5open("results/nias_DGAR.h5", "w") do file
  write(file, "t", tr)
  write(file, "tau", τ)
  write(file, "tauerr", τerr)
  write(file, "dtau", δτ)
  write(file, "dtauerr", δτerr)
  write(file, "tauargo", τargo)
  write(file, "tauecco", τecco)
end

# offset based on trends
tr = Dates.value.(t .- t[1])/1000/3600/24
for i = 1:l
  c, _ = SOT.lineartrend(t, τ[:,i]; fitannual=true, fitsemiannual=true)
  τ[:,i] .-= c[2] + c[1]*tr[end]/2
  if i == 1
    δτ[:,:] .-= c[2] + c[1]*tr[end]/2
  else
    δτ[:,i-1] .+= c[2] + c[1]*tr[end]/2
  end
  c, _ = SOT.lineartrend(t, τargo[:,i]; fitannual=true, fitsemiannual=true)
  τargo[:,i] .-= c[2] + c[1]*tr[end]/2
  c, _ = SOT.lineartrend(t, τecco[:,i]; fitannual=true, fitsemiannual=true)
  τecco[:,i] .-= c[2] + c[1]*tr[end]/2
end

# estimate trends
cτ, _ = SOT.lineartrend(t, τ[:,1]; fitannual=true, fitsemiannual=true)
cτargo, _ = SOT.lineartrend(t, τargo[:,1]; fitannual=true, fitsemiannual=true)
cτecco, _ = SOT.lineartrend(t, τecco[:,1]; fitannual=true, fitsemiannual=true)
cδτ, _ = SOT.lineartrend(t, δτ[:,l-1]; fitannual=true, fitsemiannual=true)
cδτargo, _ = SOT.lineartrend(t, τargo[:,1] - τargo[:,l]; fitannual=true, fitsemiannual=true)
cδτecco, _ = SOT.lineartrend(t, τecco[:,1] - τecco[:,l]; fitannual=true, fitsemiannual=true)

# plot measured vs. inverted T-wave delays (lowest freq.)
fig, ax = subplots(1, 1)
ax.scatter(y[1:nt,1], E[1:nt,:]*x[:,1], s=5)
ax.set_aspect(1)
xl = ax.get_xlim()
yl = ax.get_ylim()
a = [-20, 20]
ax.plot(a, a, color="black", linewidth=.8)
ax.plot(a, a .+ 1/tinvfreq[1], color="black", linewidth=.8, zorder=0)
ax.plot(a, a .- 1/tinvfreq[1], color="black", linewidth=.8, zorder=0)
ax.plot(a, a .+ 1/2tinvfreq[1], color="black", linewidth=.8, zorder=0)
ax.plot(a, a .- 1/2tinvfreq[1], color="black", linewidth=.8, zorder=0)
ax.set_xlim(xl)
ax.set_ylim(yl)
ax.set_title("T waves")
ax.set_xlabel("measured delay (s)")
ax.set_ylabel("inverted delay (s)")

# plot measured vs. inverted P-wave delays (lowest freq.)
fig, ax = subplots(1, 1)
ax.scatter(y[nt+1:nt+np], E[nt+1:nt+np,:]*x[:,1], s=5)
ax.set_aspect(1)
xl = ax.get_xlim()
yl = ax.get_ylim()
a = [-20, 20]
ax.plot(a, a, color="black", linewidth=.8)
ax.plot(a, a .+ 1/tinvfreq[1], color="black", linewidth=.8, zorder=0)
ax.plot(a, a .- 1/tinvfreq[1], color="black", linewidth=.8, zorder=0)
ax.plot(a, a .+ 1/2tinvfreq[1], color="black", linewidth=.8, zorder=0)
ax.plot(a, a .- 1/2tinvfreq[1], color="black", linewidth=.8, zorder=0)
ax.set_xlim(xl)
ax.set_ylim(yl)
ax.set_title("P waves")
ax.set_xlabel("measured delay (s)")
ax.set_ylabel("inverted delay (s)")

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

tg = copy(t)
τg = copy(τ)
τerrg = copy(τerr)
δτg = copy(δτ)
δτerrg = copy(δτerr)
τargog = copy(τargo)
τeccog = copy(τecco)
# excluded periods
for i = 1:length(excludetimes)
  global tg, τg, τerrg, δτg, δτerrg, τargog, τeccog
  tg, τg, τerrg, δτg, δτerrg, τargog, τeccog = insertgap(tg, [τg, τerrg, δτg, δτerrg,
                                                              τargog, τeccog],
                                                         excludetimes[i][1])
end

# plot timeseries
fig, ax = subplots(2, 1, figsize=(16, 6.4), sharex=true)
ax[1].plot(tg, τg[:,1], color="tab:blue", zorder=3, label=L"$T$ waves")
ax[1].scatter(t, τ[:,1], s=2, c="tab:blue", zorder=3)
ax[1].fill_between(tg, τg[:,1] - 2τerrg[:,1], τg[:,1] + 2τerrg[:,1], alpha=.25,
                   color="tab:blue", linewidths=0, zorder=3)
ax[2].plot(tg, δτg[:,l-1], color="tab:blue", zorder=3, label=L"$T$ waves")
ax[2].scatter(t, δτ[:,l-1], s=2, color="tab:blue", zorder=3)
ax[2].fill_between(tg, δτg[:,l-1] - 2δτerrg[:,l-1], δτg[:,l-1] + 2δτerrg[:,l-1], alpha=.25,
                   color="tab:blue", linewidths=0, zorder=3)
ax[1].plot(tg, τargog[:,1], color="tab:orange", zorder=1, label="Argo")
ax[2].plot(tg, τargog[:,1] - τargog[:,l], color="tab:orange", zorder=1, label="Argo")
ax[1].plot(tg, τeccog[:,1], color="tab:green", zorder=2, label="ECCO")
ax[2].plot(tg, τeccog[:,1] - τeccog[:,l], color="tab:green", zorder=2, label="ECCO")
ax[1].plot(t, cτargo[1]*tr .+ cτargo[2], color="tab:orange", zorder=1)
ax[1].plot(t, cτecco[1]*tr .+ cτecco[2], color="tab:green", zorder=2)
ax[1].plot(t, cτ[1]*tr .+ cτ[2], color="tab:blue", zorder=3)
ax[2].plot(t, cδτargo[1]*tr .+ cδτargo[2], color="tab:orange", zorder=1)
ax[2].plot(t, cδτecco[1]*tr .+ cδτecco[2], color="tab:green", zorder=2)
ax[2].plot(t, cδτ[1]*tr .+ cδτ[2], color="tab:blue", zorder=3)
ax[1].invert_yaxis()
ax[1].legend(frameon=false, loc=4, ncol=3)
ax[2].legend(frameon=false, loc=4, ncol=3)
ax[1].set_xlim(t[1] - (t[end] - t[1])÷100, t[end] + (t[end] - t[1])÷100)
ax[1].set_ylabel("travel time anomaly (s)")
ax[2].set_ylabel("travel time difference (s)")
fig.align_ylabels()
fig.tight_layout()

# plot measured vs. inverted T-wave delays with P-wave delays estimated from inversion
# without smoothing
fig, ax = subplots()
ax.set_aspect(1)
ax.scatter(y[1:nt,1] - Δτp[:,1], E[1:nt,:]*x[:,1] - Δτp[:,1], s=5)
xl = ax.get_xlim()
yl = ax.get_ylim()
a = [-20, 20]
ax.plot(a, a, color="black", linewidth=.8)
ax.plot(a, a .+ 1/tinvfreq[1], color="black", linewidth=.8, zorder=0)
ax.plot(a, a .- 1/tinvfreq[1], color="black", linewidth=.8, zorder=0)
ax.set_xlim(xl)
ax.set_ylim(yl)
ax.set_xlabel("measured delay (s)")
ax.set_ylabel("inverted delay (s)")

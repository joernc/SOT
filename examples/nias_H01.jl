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
tstations = ["H01W3"]

# T-wave time window around predicted arrival time
tintervals = [[-10, 70]]

# T-wave filtering window width
tavgwidth = 0.5

# T-wave reference frequency at which to find first max CC
treffreq = 2.5

# frequencies used in inversion
tinvfreq = 2.5:0.75:4.0

# minimum CCs for T-wave pairs (at inversion frequencies)
tmincc = 0.6:-0.1:0.4

# manually exclude pairs
excludepairs = CSV.read("data/catalogs/nias_H01_exclude.csv", DataFrame)

## download P-wave data
#SOT.downloadseisdata(eqname, pstations)
#
## cut and filter P waveforms
#SOT.cutpwaves(eqname, pstations, pintervals, pfreqbands)
#
## find P-wave pairs
#SOT.findpairs(eqname, pstations, pintervals, pfreqbands, saveplot=true)
#
## measure T-wave lags Δτ
#SOT.twavepick(eqname, tstations, tintervals, tavgwidth, treffreq, pstations, pintervals,
#              pfreqbands, saveplot=true)

# collect usable pairs
tpairs, ppairs = SOT.collectpairs(eqname, tstations, tintervals, tavgwidth, treffreq,
                                  tinvfreq, tmincc, pstations, pintervals, pfreqbands;
                                  excludepairs)

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
targo = Dates.value.(Date(2004, 1, 15) .+ Month.(h5read("data/argo/nias_H01.h5", "time"))
                     .- Date(2000, 1, 1))
tr = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
τargo = hcat([interpolate((targo,), h5read("data/argo/nias_H01.h5", "tau")[i,:],
                          Gridded(Linear()))(tr) for i = 1:l]...)
nz = size(h5read("data/argo/nias_H01.h5", "T"), 1)
Targo = hcat([interpolate((targo,), h5read("data/argo/nias_H01.h5", "T")[i,:],
                          Gridded(Linear()))(tr) for i = 1:nz]...)

# read and interpolate ECCO data
tecco = h5read("data/ecco/nias_H01.h5", "time")
tr = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
τecco = hcat([interpolate((tecco,), h5read("data/ecco/nias_H01.h5", "tau")[i,:],
                          Gridded(Linear()))(tr) for i = 1:l]...)
nz = size(h5read("data/ecco/nias_H01.h5", "T"), 1)
Tecco = hcat([interpolate((tecco,), h5read("data/ecco/nias_H01.h5", "T")[i,:],
                          Gridded(Linear()))(tr) for i = 1:nz]...)

# save timeseries to file
tr = Dates.value.(t - DateTime(2000, 1, 1, 0, 0, 0))/1000/3600/24
h5open("results/nias_H01.h5", "w") do file
  write(file, "t", tr)
  write(file, "tau", τ)
  write(file, "tauerr", τerr)
  write(file, "dtau", δτ)
  write(file, "dtauerr", δτerr)
  write(file, "tauargo", τargo)
  write(file, "tauecco", τecco)
  write(file, "Targo", Targo)
  write(file, "Tecco", Tecco)
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

# plot timeseries
fig, ax = subplots(2, 1, figsize=(16, 6.4), sharex=true)
ax[1].plot(t, τ[:,1], color="tab:blue", zorder=3, label=L"$T$ waves")
ax[1].scatter(t, τ[:,1], s=2, c="tab:blue", zorder=3)
ax[1].fill_between(t, τ[:,1] - 2τerr[:,1], τ[:,1] + 2τerr[:,1], alpha=.25,
                   color="tab:blue", linewidths=0, zorder=3)
ax[2].plot(t, δτ[:,l-1], color="tab:blue", zorder=3, label=L"$T$ waves")
ax[2].scatter(t, δτ[:,l-1], s=2, color="tab:blue", zorder=3)
ax[2].fill_between(t, δτ[:,l-1] - 2δτerr[:,l-1], δτ[:,l-1] + 2δτerr[:,l-1], alpha=.25,
                   color="tab:blue", linewidths=0, zorder=3)
ax[1].plot(t, τargo[:,1], color="tab:orange", zorder=1, label="Argo")
ax[2].plot(t, τargo[:,1] - τargo[:,l], color="tab:orange", zorder=1, label="Argo")
ax[1].plot(t, τecco[:,1], color="tab:green", zorder=2, label="ECCO")
ax[2].plot(t, τecco[:,1] - τecco[:,l], color="tab:green", zorder=2, label="ECCO")
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

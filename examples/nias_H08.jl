using .SOT, PyPlot, Printf, Dates, LinearAlgebra, HDF5, Interpolations

# identifier for experiment
eqname = "nias"

# P-wave (reference) stations
pstations = ["PS.PSI..BHZ", "MY.KUM..BHZ", "II.WRAB.00.BHZ", "GE.GSI..BHZ"]

# intervals to which to cut P waveforms
pintervals = [[-3, 47], [-3, 47], [-3, 47], [-3, 47]]

# frequency bands to which to filter P waveforms
pfreqbands = [[1, 3], [1, 3], [1.5, 2.5], [1, 3]]

# T-wave station
tstations = ["H08S2"]

# T-wave time window around predicted arrival time
tintervals = [[-10, 70]]

# T-wave filtering window width
tavgwidth = 0.5

# T-wave reference frequency at which to find first max CC
treffreq = 2.0

# frequencies used in inversion
tinvfreq = [2.0, 3.0, 4.0]

# minimum CCs for T-wave pairs (at inversion frequencies)
tmincc = [0.6, 0.5, 0.4]

# excluded time periods: before 2004-08-01 and periods with uncorrected clock error
excludetimes = [[Date(2001, 1, 1) Date(2004, 8, 1)],
                [DateTime("2010-03-16T00:00:00") DateTime("2010-05-17T02:06:17.760")],
                [Date(2017, 6, 1) Date(2018, 1, 1)]]

# download P-wave data
SOT.downloadpwaves(eqname, pstations)

# cut and filter P waveforms
SOT.cutpwaves(eqname, pstations, pintervals, pfreqbands)

# find P-wave pairs
SOT.findpairs(eqname, pstations, pintervals, pfreqbands, saveplot=true)

# measure T-wave lags Δτ
SOT.twavepick(eqname, tstations, tintervals, tavgwidth, treffreq, pstations, pintervals,
              pfreqbands, saveplot=true)

# collect usable pairs
tpairs, ppairs = SOT.collectpairs(eqname, tstations, tintervals, tavgwidth, treffreq,
                                  tinvfreq, tmincc, pstations, pintervals, pfreqbands;
                                  excludetimes)

# H08 clock error correction
function timingcorrection!(tpairs, starttime, endtime, c)
  idx1 = starttime .< tpairs.event1 .< endtime
  idx2 = starttime .< tpairs.event2 .< endtime
  for i = findall(idx1)
    tpairs.Δτl[i] .-= c*((tpairs.event1[i] .- starttime)./(endtime - starttime))
    tpairs.Δτc[i] .-= c*((tpairs.event1[i] .- starttime)./(endtime - starttime))
    tpairs.Δτr[i] .-= c*((tpairs.event1[i] .- starttime)./(endtime - starttime))
  end
  for i = findall(idx2)
    tpairs.Δτl[i] .+= c*((tpairs.event2[i] .- starttime)./(endtime - starttime))
    tpairs.Δτc[i] .+= c*((tpairs.event2[i] .- starttime)./(endtime - starttime))
    tpairs.Δτr[i] .+= c*((tpairs.event2[i] .- starttime)./(endtime - starttime))
  end
end

# correct based on comparison with DGAR
timingcorrection!(tpairs, DateTime("2010-01-23T08:13:52.666"),
                  DateTime("2010-03-16T00:00:00"), 2.311)
timingcorrection!(tpairs, DateTime("2010-05-17T02:06:17.760"),
                  DateTime("2010-10-09T23:59:52.747"), 7.060)
timingcorrection!(tpairs, DateTime("2010-10-09T23:59:52.747"),
                  DateTime("2010-12-16T02:03:23.040"), 3.427)
timingcorrection!(tpairs, DateTime("2010-12-16T02:03:23.040"),
                  DateTime("2011-02-22T09:43:45.858"), 3.495)
timingcorrection!(tpairs, DateTime("2011-02-22T09:43:45.858"),
                  DateTime("2011-03-24T07:55:37.542"), 1.593)
timingcorrection!(tpairs, DateTime("2011-03-24T07:55:37.542"),
                  DateTime("2011-04-18T01:59:00.042"), 1.276)
timingcorrection!(tpairs, DateTime("2011-04-18T01:59:00.042"),
                  DateTime("2011-08-27T23:59:52.317"), 6.678)
timingcorrection!(tpairs, DateTime("2011-08-27T23:59:52.317"),
                  DateTime("2011-12-24T23:59:52.453"), 6.786)
timingcorrection!(tpairs, DateTime("2011-12-24T23:59:52.453"),
                  DateTime("2012-01-20T00:00:00.000"), 1.495)

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
Δτp = E[1:nt,1:m]*(E[nt+1:nt+np,m+1:2m]\ppairs.Δτ)

# make cycle-skipping correction
tpairs.Δτ = SOT.correctcycleskipping(tpairs, ppairs, E, S, P)

# get travel time anomalies
y = [vcat(tpairs.Δτ'...); repeat(ppairs.Δτ, 1, l)]
x = P*(E'*y)
τ = D*x

# calculate error
σ2 = 1/(nt+np-m-1)*sum((y - E*x).^2, dims=1)
A = D*P*E'
τerr = sqrt.(σ2.*diag(A*A'))

# get travel time differences between frequencies
δy = y[:,1] .- y[:,2:l]
δx = x[:,1] .- x[:,2:l]
δτ = D*δx

# calculate error
σ2 = 1/(nt+np-m-1)*sum((δy - E*δx).^2, dims=1)
δτerr = sqrt.(σ2.*diag(A*A'))

# read and interpolate ECCO data
tecco = h5read("data/ecco/nias_H08.h5", "time")
tr = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
τecco = hcat([interpolate((tecco,), h5read("data/ecco/nias_H08.h5", "tau")[:,i],
                          Gridded(Linear()))(tr) for i = 1:2]...)

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
colors = matplotlib.rcParams["axes.prop_cycle"].by_key()["color"]
fig, ax = subplots(2, 1, figsize=(16, 6.4), sharex=true)
ax[1].plot(t, τecco[:,2], color="tab:gray", zorder=0)
ax[2].plot(t, τecco[:,1] - τecco[:,2], color="tab:gray", zorder=0)
for i = 1:l
  ax[1].plot(t, τ[:,i], color=colors[i], label=@sprintf("%3.1f Hz", tinvfreq[i]))
  ax[1].scatter(t, τ[:,i], s=5, c=colors[i])
  ax[1].fill_between(t, τ[:,i] - 2τerr[:,i], τ[:,i] + 2τerr[:,i], alpha=.25,
                     color=colors[i], linewidths=0)
  if i > 1
    ax[2].plot(t, δτ[:,i-1], color=colors[i], label=@sprintf("%3.1f Hz – %3.1f Hz",
                                                             tinvfreq[1], tinvfreq[i]))
    ax[2].scatter(t, δτ[:,i-1], s=5, color=colors[i])
    ax[2].fill_between(t, δτ[:,i-1] - 2δτerr[:,i-1], δτ[:,i-1] + 2δτerr[:,i-1], alpha=.25,
                       color=colors[i], linewidths=0)
  end
end
ax[1].invert_yaxis()
ax[1].legend(frameon=false)
ax[2].legend(frameon=false)
ax[1].axhline(0, color="black", linewidth=.8, zorder=0)
ax[2].axhline(0, color="black", linewidth=.8, zorder=0)
ax[1].set_xlim(t[1] - (t[end] - t[1])÷100, t[end] + (t[end] - t[1])÷100)
ax[1].set_ylabel("travel time anomaly (s)")
ax[2].set_ylabel("travel time difference (s)")
fig.tight_layout()

# plot measured vs. inverted delays for all possible combination of T and P pairs
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

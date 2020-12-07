# TODO:
# - Understand large anomalies in 2017 — clock error?

using .SOT, PyPlot, Printf, Dates

# identifier for experiment
eqname = "nias"

# P-wave (reference) stations
pstations = ["PSI", "KUM", "WRAB"]

# T-wave station and channel
tstations = ["H08S2..EDH"]

# time window (seconds before and after predicted arrival time)
starttime = 10
endtime = 70

# frequencies at which to measure travel time change
frequencies = 0.1:0.1:10

# frequency averaging width
avgwidth = 0.5

# reference frequency at which to pick max CC
reffreq = 2.0

# frequencies used in inversion
invfreq = [2.0, 4.0]

# minimum requirement for CC at the frequencies used for inversion
mincc = [0.6, 0.4]

# maximum allowable T-wave |Δτ| (discard outliers)
maxΔτt = 20.

# excluded time periods: before 2004-12-01 and period with clock error
excludetimes = [[Date(2001, 1, 1) Date(2004, 12, 1)], [Date(2010, 1, 1) Date(2012, 1, 20)]]

# measure T-wave lags Δτ
for s in tstations
  SOT.twavepick(eqname, s, pstations, starttime, endtime, frequencies, avgwidth, reffreq;
                saveplot=true)
end

# invert for travel time anomalies τ
t, τ, τerr, tpairs, ppairs = SOT.invert(eqname, tstations, pstations, invfreq, mincc;
                                        maxΔτt, excludetimes, csc=true)

# number of used T- and P-wave pairs
nt = size(tpairs, 1)
np = size(ppairs, 1)

# plot measured vs. inverted T-wave delays (lowest freq.)
fig, ax = subplots(1, 1)
ax.scatter([tpairs.Δτ[i,1] for i = 1:nt], [tpairs.Δτi[i,1] for i = 1:nt], s=5)
ax.set_aspect(1)
xl = ax.get_xlim()
yl = ax.get_ylim()
x = [-2maxΔτt, 2maxΔτt]
ax.plot(x, x, color="black", linewidth=.8)
ax.plot(x, x .+ 1/invfreq[1], color="black", linewidth=.8, zorder=0)
ax.plot(x, x .- 1/invfreq[1], color="black", linewidth=.8, zorder=0)
ax.plot(x, x .+ 1/2invfreq[1], color="black", linewidth=.8, zorder=0)
ax.plot(x, x .- 1/2invfreq[1], color="black", linewidth=.8, zorder=0)
ax.set_xlim(xl)
ax.set_ylim(yl)
ax.set_title("T waves")
ax.set_xlabel("measured delay (s)")
ax.set_ylabel("inverted delay (s)")

# plot measured vs. inverted P-wave delays (lowest freq.)
fig, ax = subplots(1, 1)
ax.scatter(ppairs.Δτ, [ppairs.Δτi[i][1] for i = 1:np], s=5)
ax.set_aspect(1)
xl = ax.get_xlim()
yl = ax.get_ylim()
x = [-2maxΔτt, 2maxΔτt]
ax.plot(x, x, color="black", linewidth=.8)
ax.plot(x, x .+ 1/invfreq[1], color="black", linewidth=.8, zorder=0)
ax.plot(x, x .- 1/invfreq[1], color="black", linewidth=.8, zorder=0)
ax.plot(x, x .+ 1/2invfreq[1], color="black", linewidth=.8, zorder=0)
ax.plot(x, x .- 1/2invfreq[1], color="black", linewidth=.8, zorder=0)
ax.set_xlim(xl)
ax.set_ylim(yl)
ax.set_title("P waves")
ax.set_xlabel("measured delay (s)")
ax.set_ylabel("inverted delay (s)")

# plot timeseries
colors = matplotlib.rcParams["axes.prop_cycle"].by_key()["color"]
fig, ax = subplots(2, 1, figsize=(16, 6.4), sharex=true)
for i = 1:length(invfreq)
  ax[1].plot(t, τ[:,i], color=colors[i], label=@sprintf("%3.1f Hz", invfreq[i]))
  ax[1].scatter(t, τ[:,i], s=5, c=colors[i])
  ax[1].fill_between(t, τ[:,i] - 2τerr[:,i], τ[:,i] + 2τerr[:,i], alpha=.25,
                     color=colors[i], linewidths=0)
  if i > 1
    δτ = τ[:,1] - τ[:,i]
    δτerr = sqrt.(τerr[:,1].^2 + τerr[:,i].^2)
    ax[2].plot(t, δτ, color=colors[i], label=@sprintf("%3.1f Hz – %3.1f Hz", invfreq[1],
                                                      invfreq[i]))
    ax[2].scatter(t, δτ, s=5, color=colors[i])
    ax[2].fill_between(t, δτ - 2δτerr, δτ + 2δτerr, alpha=.25, color=colors[i],
                       linewidths=0)
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

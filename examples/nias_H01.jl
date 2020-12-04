# TODO:
# - remove pairs that cause |Δτt| > 20

using .SOT, PyPlot, Printf, Dates

# identifier for experiment
eqname = "nias"

# P-wave (reference) stations
pstations = ["PSI", "KUM", "WRAB"]

# T-wave station and channel
tstation = "H01W3..EDH"

# time window (seconds before and after predicted arrival time)
starttime = 10
endtime = 70

# frequencies at which to measure travel time change
frequencies = 0.1:0.1:10

# frequency averaging width
avgwidth = 0.5

# reference frequency at which to pick max CC
reffreq = 2.5

# frequencies used in inversion
invfreq = [2.5, 4.0]

# minimum requirement for CC at the frequencies used for inversion
mincc = [0.6, 0.4]

# maximum allowable T-wave |Δτ| (discard outliers)
maxΔτt = 20.

# excluded time periods: before 2004-12-01 and period with clock error
excludetimes = [[Date(2001, 1, 1) Date(2004, 12, 1)]]

# measure T-wave lags Δτ
SOT.twavepick(eqname, tstation, pstations, starttime, endtime, frequencies, avgwidth,
              reffreq; saveplot=true)

# invert for travel time anomalies τ
t, τ, τerr = SOT.invert(eqname, tstation, pstations, invfreq, mincc; maxΔτt, excludetimes)

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

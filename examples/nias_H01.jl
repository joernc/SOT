# TODO:
# - remove pairs that cause |Δτt| > 20
# - update to process P waves

using .SOT, PyPlot, Printf, Dates, LinearAlgebra

# identifier for experiment
eqname = "nias"

# P-wave (reference) stations
pstations = ["PSI", "KUM", "WRAB"]

# T-wave station and channel
tstations = ["H01W1..EDH", "H01W2..EDH", "H01W3..EDH"]

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
maxΔτ = 20.

# excluded time periods: before 2004-12-01 and period with clock error
excludetimes = [[Date(2001, 1, 1) Date(2004, 12, 1)]]

# number of frequencies
l = length(invfreq)

# measure T-wave lags Δτ
for s in tstations
  SOT.twavepick(eqname, s, pstations, starttime, endtime, frequencies, avgwidth, reffreq;
                saveplot=true)
end

# collect usable pairs
tpairs, ppairs = SOT.collectpairs(eqname, tstations, pstations, invfreq, mincc;
                                  maxΔτ, excludetimes)

# perform inversion
t, E, S, P, D = SOT.invert(tpairs, ppairs)

# number of good T- and P-wave pairs
nt = size(tpairs, 1)
np = size(ppairs, 1)

# number of unique events
m = length(t)

# get travel time anomalies
y = [vcat(tpairs.Δτc'...); repeat(ppairs.Δτ, 1, l)]
x = P*(E'*y)
τ = D*x

# calculate error
σ2 = 1/(nt+np-m-1)*sum((y - E*x).^2, dims=1)
A = D*P*E'
τerr = sqrt.(σ2.*diag(A*A'))

# record inverted delays
tpairs.Δτi = collect(eachrow(E[1:nt,:]*x))
ppairs.Δτi = collect(eachrow(E[nt+1:nt+np,:]*x))

# get travel time differences between frequencies
δy = y[:,1] - y[:,2:l]
δx = x[:,1] - x[:,2:l]
δτ = D*δx

# calculate error
σ2 = 1/(nt+np-m-1)*sum((δy - E*δx).^2, dims=1)
δτerr = sqrt.(σ2.*diag(A*A'))

# plot measured vs. inverted T-wave delays (lowest freq.)
fig, ax = subplots(1, 1)
ax.scatter(y[1:nt,1], E[1:nt,:]*x[:,1], s=5)
ax.set_aspect(1)
xl = ax.get_xlim()
yl = ax.get_ylim()
a = [-2maxΔτ, 2maxΔτ]
ax.plot(a, a, color="black", linewidth=.8)
ax.plot(a, a .+ 1/invfreq[1], color="black", linewidth=.8, zorder=0)
ax.plot(a, a .- 1/invfreq[1], color="black", linewidth=.8, zorder=0)
ax.plot(a, a .+ 1/2invfreq[1], color="black", linewidth=.8, zorder=0)
ax.plot(a, a .- 1/2invfreq[1], color="black", linewidth=.8, zorder=0)
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
a = [-2maxΔτ, 2maxΔτ]
ax.plot(a, a, color="black", linewidth=.8)
ax.plot(a, a .+ 1/invfreq[1], color="black", linewidth=.8, zorder=0)
ax.plot(a, a .- 1/invfreq[1], color="black", linewidth=.8, zorder=0)
ax.plot(a, a .+ 1/2invfreq[1], color="black", linewidth=.8, zorder=0)
ax.plot(a, a .- 1/2invfreq[1], color="black", linewidth=.8, zorder=0)
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
    ax[2].plot(t, δτ[:,i-1], color=colors[i], label=@sprintf("%3.1f Hz – %3.1f Hz",
                                                             invfreq[1], invfreq[i]))
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

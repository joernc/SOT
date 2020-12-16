# TODO:
# - cut T waveforms
# - allow different processing parameters for different T-wave stations
# - save T-wave data in folders specific to processing parameters

using .SOT, PyPlot, Printf, Dates, LinearAlgebra

# identifier for experiment
eqname = "nias"

# P-wave (reference) stations
pstations = ["PS.PSI..BHZ", "MY.KUM..BHZ", "II.WRAB.00.BHZ"]

# intervals to which to cut P waveforms
pintervals = [[-3, 47], [-3, 47], [-3, 47]]

# frequency bands to which to filter P waveforms
pfreqbands = [[1, 3], [1, 3], [1.5, 2.5]]

# T-wave station and channel
tstations = ["H08S2..EDH"]

# frequencies used in inversion
tinvfreq = [2.0, 4.0]

# excluded time periods: before 2004-12-01 and period with clock error
excludetimes = [[Date(2001, 1, 1) Date(2004, 12, 1)],
                [Date(2010, 1, 1) Date(2012, 1, 20)],
                [Date(2017, 6, 1) Date(2018, 1, 1)]]

# download P-wave data
SOT.downloadpwaves(eqname, pstations)

# cut and filter P waveforms
SOT.cutpwaves(eqname, pstations, pintervals, pfreqbands)

# find P-wave pairs
SOT.findpairs(eqname, pstations, pintervals, pfreqbands, saveplot=true)

# measure T-wave lags Δτ
for s in tstations
  SOT.twavepick(eqname, s, pstations, [-10, 70], 0.1:0.1:10, 0.5, 2.0; saveplot=true)
end

# collect usable pairs
tpairs, ppairs = SOT.collectpairs(eqname, tstations, pstations, tinvfreq, [0.6, 0.4];
                                  excludetimes)

# perform inversion
t, E, S, P, D = SOT.invert(tpairs, ppairs)

# number of good T- and P-wave pairs
nt = size(tpairs, 1)
np = size(ppairs, 1)

# number of unique events
m = length(t)

# Apply cycle-skipping correction
tpairs.Δτ = SOT.correctcycleskipping(tpairs, ppairs, E, S, P)

# number of frequencies
l = length(tinvfreq)

# get travel time anomalies
y = [vcat(tpairs.Δτ'...); repeat(ppairs.Δτ, 1, l)]
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

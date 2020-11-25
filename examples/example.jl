using SOT, CSV, DataFrames, PyPlot

# identifier for experiment
eqname = "nias"

# P-wave (reference) stations
pstations = ["PSI", "GSI"]

# T-wave station and channel
tstation = "H08S2..EDH"

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

# maximum allowable |Δτ| (discard outliers)
maxΔτ = 0.65

# excluded time periods: before 2004-12-01 and period with clock error
excludetimes = [[Date(2001, 1, 1) Date(2004, 12, 1)], [Date(2010, 1, 1) Date(2012, 1, 20)]]

# measure T-wave lags Δτ
for s in pstations
  SOT.twavepick(eqname, s, tstation, starttime, endtime, frequencies, avgwidth, reffreq;
                saveplot=true)
end

# invert for travel time anomalies τ
t, τ, τerr = SOT.invert(eqname, pstations, tstation, invfreq, mincc; maxΔτ=maxΔτ,
                        excludetimes=excludetimes)

# simple plot
colors = matplotlib.rcParams["axes.prop_cycle"].by_key()["color"]
figure()
for i = 1:length(invfreq)
  plot(t, τ[:,i], color=colors[i])
  scatter(t, τ[:,i], s=5, c=colors[i])
  fill_between(t, τ[:,i] - 2τerr[:,i], τ[:,i] + 2τerr[:,i], alpha=.25, color=colors[i], linewidths=0)
end

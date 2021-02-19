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
tstations = ["II.DGAR.00.BHZ", "H08S2"]

# T-wave time window around predicted arrival time
tintervals = [[-10, 70], [-10, 70]]

# T-wave filtering window width
tavgwidth = 0.5

# T-wave reference frequency at which to find first max CC
treffreq = 2.0

# frequencies used in inversion
tinvfreq = [2.0, 3.0]

# minimum CCs for T-wave pairs (at inversion frequencies)
tmincc = [0.6, 0.5]

# excluded time period: before 2004-12-01
excludetimes = [[Date(2001, 1, 1) Date(2004, 12, 1)],
                [DateTime("2010-03-16T00:00:00") DateTime("2010-05-17T02:06:17.760")],
                [Date(2017, 6, 1) Date(2021, 1, 1)]]

# manually exclude pairs
excludepairs = [CSV.read("data/catalogs/nias_DGAR_exclude.csv", DataFrame),
                CSV.read("data/catalogs/nias_H08_exclude.csv", DataFrame)]
append!(excludepairs[2], DataFrame(pstation="PS.PSI..BHZ",
                                   event1=DateTime("2008-03-03T19:53:32.85"),
                                   event2=DateTime("2013-12-27T04:18:41.55"),
                                   reason="wrong CSC"))
append!(excludepairs[2], DataFrame(pstation="MY.KUM..BHZ",
                                   event1=DateTime("2008-03-03T19:53:32.85"),
                                   event2=DateTime("2013-12-27T04:18:41.55"),
                                   reason="wrong CSC"))
append!(excludepairs[2], DataFrame(pstation="GE.GSI..BHZ",
                                   event1=DateTime("2008-03-03T19:53:32.85"),
                                   event2=DateTime("2013-12-27T04:18:41.55"),
                                   reason="wrong CSC"))

# collect usable pairs
tpairs1, ppairs1 = SOT.collectpairs(eqname, [tstations[1]], tintervals, tavgwidth, treffreq,
                                    tinvfreq, tmincc, pstations, pintervals, pfreqbands;
                                    excludetimes, excludepairs=excludepairs[1])
tpairs2, ppairs2 = SOT.collectpairs(eqname, [tstations[2]], tintervals, tavgwidth, treffreq,
                                    tinvfreq, tmincc, pstations, pintervals, pfreqbands;
                                    excludetimes, excludepairs=excludepairs[2])

# DGAR clock error correction
breaktime = DateTime(2012, 3, 17)
idx = (tpairs1.event1 .< breaktime) .& (tpairs1.event2 .> breaktime)
for i = findall(idx)
  tpairs1.Δτl[i] .-= 1.02
  tpairs1.Δτc[i] .-= 1.02
  tpairs1.Δτr[i] .-= 1.02
end

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
timingcorrection!(tpairs2, DateTime("2010-01-23T08:13:52.666"),
                  DateTime("2010-03-16T00:00:00.000"), 2.276)
timingcorrection!(tpairs2, DateTime("2010-05-17T02:06:17.760"),
                  DateTime("2010-10-09T23:59:52.747"), 7.045)
timingcorrection!(tpairs2, DateTime("2010-10-09T23:59:52.747"),
                  DateTime("2010-12-16T02:03:23.040"), 3.454)
timingcorrection!(tpairs2, DateTime("2010-12-16T02:03:23.040"),
                  DateTime("2011-02-22T09:43:45.858"), 3.422)
timingcorrection!(tpairs2, DateTime("2011-02-22T09:43:45.858"),
                  DateTime("2011-03-24T07:55:37.542"), 1.552)
timingcorrection!(tpairs2, DateTime("2011-03-24T07:55:37.542"),
                  DateTime("2011-04-18T01:59:00.042"), 1.278)
timingcorrection!(tpairs2, DateTime("2011-04-18T01:59:00.042"),
                  DateTime("2011-08-27T23:59:52.317"), 6.714)
timingcorrection!(tpairs2, DateTime("2011-08-27T23:59:52.317"),
                  DateTime("2011-12-24T23:59:52.453"), 6.683)
timingcorrection!(tpairs2, DateTime("2011-12-24T23:59:52.453"),
                  DateTime("2012-01-20T00:00:00.000"), 1.374)

# find common pairs
tidx1 = [any([(tpairs1.event1[i] .== tpairs2.event1[j]) .& (tpairs1.event2[i] .== tpairs2.event2[j]) for j = 1 : size(tpairs2, 1)]) for i = 1 : size(tpairs1, 1)]
pidx1 = [any([(ppairs1.event1[i] .== ppairs2.event1[j]) .& (ppairs1.event2[i] .== ppairs2.event2[j]) for j = 1 : size(ppairs2, 1)]) for i = 1 : size(ppairs1, 1)]
tidx2 = [any([(tpairs2.event1[i] .== tpairs1.event1[j]) .& (tpairs2.event2[i] .== tpairs1.event2[j]) for j = 1 : size(tpairs1, 1)]) for i = 1 : size(tpairs2, 1)]
pidx2 = [any([(ppairs2.event1[i] .== ppairs1.event1[j]) .& (ppairs2.event2[i] .== ppairs1.event2[j]) for j = 1 : size(ppairs1, 1)]) for i = 1 : size(ppairs2, 1)]

# select common pairs
tpairs1 = tpairs1[tidx1,:]
ppairs1 = ppairs1[pidx1,:]
tpairs2 = tpairs2[tidx2,:]
ppairs2 = ppairs2[pidx2,:]

# perform inversion
t, E, S, P, D = SOT.invert(tpairs1, ppairs1)

# number of frequencies
l = length(tinvfreq)

# number of good T- and P-wave pairs
nt = size(tpairs1, 1)
np = size(ppairs1, 1)

# number of unique events
m = length(t)

# make cycle-skipping correction
tpairs1.Δτ = SOT.correctcycleskipping(tpairs1, ppairs1, E, S, P)
tpairs2.Δτ = SOT.correctcycleskipping(tpairs2, ppairs2, E, S, P)

# data vector
y1 = [reshape(vcat(tpairs1.Δτ'...), l*nt); ppairs1.Δτ]
y2 = [reshape(vcat(tpairs2.Δτ'...), l*nt); ppairs2.Δτ]

# get least-squares solution
x1 = P*(E'*y1)
x2 = P*(E'*y2)

# get travel time anomalies
τ1 = reshape(D*x1, (m, l))
τ2 = reshape(D*x2, (m, l))

# plot differences
figure()
plot(t, τ1[:,1] - τ2[:,1])
scatter(t, τ1[:,1] - τ2[:,1], s=5)
xlim(t[1]-Day(30), t[m]+Day(30))
ylabel("DGAR minus H08 (s)")

# calculate trends
cτ1, _ = SOT.lineartrend(t, τ1[:,1]; fitannual=true, fitsemiannual=true)
cτ2, _ = SOT.lineartrend(t, τ2[:,1]; fitannual=true, fitsemiannual=true)
tr = Dates.value.(t .- t[1])/1000/3600/24

@printf("DGAR trend: %3.1f mK/yr", -cτ1[1]/5.4*SOT.meanyear*1e3)
@printf("H08 trend:  %3.1f mK/yr", -cτ2[1]/5.4*SOT.meanyear*1e3)

# plot both time series
figure()
plot(t, τ1[:,1], zorder=1)
plot(t, τ2[:,1], zorder=2)
scatter(t, τ1[:,1], s=5, zorder=1)
scatter(t, τ2[:,1], s=5, zorder=2)
plot(t, cτ1[1]*tr .+ cτ1[2], zorder=3, color="tab:blue")
plot(t, cτ2[1]*tr .+ cτ2[2], zorder=3, color="tab:orange")

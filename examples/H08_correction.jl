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
                [Date(2017, 6, 1) Date(2021, 1, 1)]]
                #[Date(2016, 7, 1) Date(2021, 1, 1)]]

# manually exclude pairs
excludepairs = CSV.read("data/catalogs/nias_DGAR_exclude.csv", DataFrame)
append!(excludepairs, CSV.read("data/catalogs/nias_H08_exclude.csv", DataFrame))
unique!(excludepairs)

# collect usable pairs
tpairs, ppairs = SOT.collectpairs(eqname, tstations, tintervals, tavgwidth, treffreq,
                                  tinvfreq, tmincc, pstations, pintervals, pfreqbands;
                                  excludetimes, excludepairs)

# DGAR clock error correction
breaktime = DateTime(2012, 3, 17)
idx = (tpairs.station .== "II.DGAR.00.BHZ") .& (tpairs.event1 .< breaktime) .& (tpairs.event2 .> breaktime)
for i = findall(idx)
  tpairs.Δτl[i] .-= 1.02
  tpairs.Δτc[i] .-= 1.02
  tpairs.Δτr[i] .-= 1.02
end

# find common pairs
gdf = groupby(tpairs, :station)
common = innerjoin(DataFrame(gdf[1]), DataFrame(gdf[2]); on=[:event1, :event2], makeunique=true)

# lowest-freq delay difference
d = [common.Δτc[i][1] - common.Δτc_1[i][1] for i = 1 : size(common, 1)]

# find pairs that have at least one member outside of corrupt period
idx1 = common.event2 .> Date(2012, 1, 20)
idx2 = common.event1 .< Date(2010, 1, 23)

# collect these pairs and the delay differences
e = [common.event1[idx1]; common.event2[idx2]]
d = [-d[idx1]; d[idx2]]

# sort
idx = sortperm(e)
e = e[idx]
d = d[idx]

# time periods over which to correct (based on gap data)
periods = [[DateTime("2010-01-23T08:13:52.666"), DateTime("2010-03-16T00:00:00.000")],
           [DateTime("2010-05-17T02:06:17.760"), DateTime("2010-10-09T23:59:52.747")],
           [DateTime("2010-10-09T23:59:52.747"), DateTime("2010-12-16T02:03:23.040")],
           [DateTime("2010-12-16T02:03:23.040"), DateTime("2011-02-22T09:43:45.858")],
           [DateTime("2011-02-22T09:43:45.858"), DateTime("2011-03-24T07:55:37.542")],
           [DateTime("2011-03-24T07:55:37.542"), DateTime("2011-04-18T01:59:00.042")],
           [DateTime("2011-04-18T01:59:00.042"), DateTime("2011-08-27T23:59:52.317")],
           [DateTime("2011-08-27T23:59:52.317"), DateTime("2011-12-24T23:59:52.453")],
           [DateTime("2011-12-24T23:59:52.453"), DateTime("2012-01-20T00:00:00.000")]]

# exclude pairs affected by cycle skipping
cc = [1148, 1182, 1214, 1292]
deleteat!(e, cc)
deleteat!(d, cc)

# fit one single drift

A = Float64[]
b = Float64[]
for p = periods
  idx = p[1] .≤ e .≤ p[2]
  append!(A, Float64.(Dates.value.(e[idx] - p[1])))
  append!(b, d[idx])
end

s = A\b

figure()
for p = periods
  idx = p[1] .≤ e .≤ p[2]
  t = Float64.(Dates.value.(e[idx] - p[1]))
  scatter(e[idx], d[idx], color="tab:blue")
  scatter([p[1]], [0], color="tab:orange")
  plot([p[1]; e[idx]], [0; s*t], color="tab:blue", zorder=0)
end

# fit independent drift for each period

figure()
for p = periods
  idx = p[1] .≤ e .≤ p[2]
  t = Float64.(Dates.value.(e[idx] - p[1]))
  s = t\d[idx]
  @printf("%6.3e %5.3f\n", s, s*Dates.value(p[2] - p[1]))
  scatter(e[idx], d[idx], color="tab:blue")
  scatter([p[1]], [0], color="tab:orange")
  plot([p[1]; e[idx]], [0; s*t], color="tab:blue", zorder=0)
end

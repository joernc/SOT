using HDF5, LinearAlgebra, PyPlot, NCDatasets, Statistics, Dates, Printf

# load time series
t8 = DateTime(2000, 1, 1) .+ Millisecond.(h5read("results/nias_tm_H08.h5", "t"))
td = DateTime(2000, 1, 1) .+ Millisecond.(h5read("results/nias_DGAR.h5", "t"))
τ8 = h5read("results/nias_tm_H08.h5", "τ")
τd = h5read("results/nias_DGAR.h5", "τ")
eτ8 = h5read("results/nias_tm_H08.h5", "eτ")
eτd = h5read("results/nias_DGAR.h5", "eτ")

# gaps
gaps8 = [Date(2004, 1, 1), Date(2007, 1, 1), Date(2009, 7, 1), Date(2011, 1, 1), Date(2020, 1, 1)]
gapsd = [Date(2004, 1, 1), Date(2007, 3, 1), Date(2020, 1, 1)]
ir8 = [findfirst(t8 .> gaps8[i]) : findlast(t8 .< gaps8[i+1]) for i = 1:length(gaps8)-1]
ird = [findfirst(td .> gapsd[i]) : findlast(td .< gapsd[i+1]) for i = 1:length(gapsd)-1]

# plot time series
fig, ax = subplots(1, 1, sharex=true, figsize=(190/25.4, 90/25.4))
ax.scatter(t8, τ8[:,1], s=2, zorder=2, color="tab:blue")
ax.scatter(td, τd[:,1], s=2, zorder=2, color="tab:orange")
for j = ir8
  global p1
  p1, = ax.plot(t8[j], τ8[j,1], zorder=3, color="tab:blue", linewidth=1)
  ax.fill_between(t8[j], τ8[j,1]-2eτ8[j,1], τ8[j,1]+2eτ8[j,1], alpha=.2, zorder=3, color="tab:blue", linewidth=0)
end
for j = ird
  global p2
  p2, = ax.plot(td[j], τd[j,1], zorder=3, color="tab:orange", linewidth=1)
  ax.fill_between(td[j], τd[j,1]-2eτd[j,1], τd[j,1]+2eτd[j,1], alpha=.2, zorder=3, color="tab:orange", linewidth=0)
end
ax.set_ylabel("travel time anomaly (s)")
ylim = maximum(abs.(ax.get_ylim()))
ax.set_ylim(-ylim, ylim)
ax.set_xlim(td[1]-Day(30), td[end]+Day(30))
ax.legend([p1, p2], ["H08", "DGAR"], ncol=3, loc="upper center", frameon=false)
ax.set_xticks(ceil(td[1], Year) : Year(1) : floor(td[end], Year), minor=true)
fig.tight_layout()

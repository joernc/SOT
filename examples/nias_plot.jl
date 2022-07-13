using HDF5, LinearAlgebra, PyPlot, NCDatasets, Statistics, Dates, Printf

# load time series
t8 = DateTime(2000, 1, 1) .+ Millisecond.(h5read("results/nias_tm_H08.h5", "t"))
td = DateTime(2000, 1, 1) .+ Millisecond.(h5read("results/nias_DGAR.h5", "t"))
τ8 = h5read("results/nias_tm_H08.h5", "τ")
τd = h5read("results/nias_DGAR.h5", "τ")
eτ8 = h5read("results/nias_tm_H08.h5", "eτ")
eτd = h5read("results/nias_DGAR.h5", "eτ")

# load ECCO time series
τecco8 = h5read("results/nias_tm_H08.h5", "τecco")
τeccod = h5read("results/nias_DGAR.h5", "τecco")

# gaps
gaps8 = [Date(2004, 1, 1), Date(2007, 1, 1), Date(2009, 7, 1), Date(2011, 1, 1), Date(2020, 1, 1)]
gapsd = [Date(2004, 1, 1), Date(2007, 3, 1), Date(2020, 1, 1)]
ir8 = [findfirst(t8 .> gaps8[i]) : findlast(t8 .< gaps8[i+1]) for i = 1:length(gaps8)-1]
ird = [findfirst(td .> gapsd[i]) : findlast(td .< gapsd[i+1]) for i = 1:length(gapsd)-1]

# plot time series
fig, ax = subplots(2, 1, sharex=true, figsize=(190/25.4, 150/25.4))
ax[1].scatter(t8, τ8[:,1], s=2, zorder=2, color="tab:blue")
ax[1].scatter(td, τd[:,1], s=2, zorder=2, color="tab:orange")
ax[2].scatter(t8, τecco8[:,1], s=2, zorder=2, color="tab:blue")
ax[2].scatter(td, τeccod[:,1], s=2, zorder=2, color="tab:orange")
for j = ir8
  global p1
  p1, = ax[1].plot(t8[j], τ8[j,1], zorder=3, color="tab:blue", linewidth=1)
  ax[1].fill_between(t8[j], τ8[j,1]-2eτ8[j,1], τ8[j,1]+2eτ8[j,1], alpha=.2, zorder=3, color="tab:blue", linewidth=0)
  p1, = ax[2].plot(t8[j], τecco8[j,1], zorder=3, color="tab:blue", linewidth=1)
end
for j = ird
  global p2
  p2, = ax[1].plot(td[j], τd[j,1], zorder=3, color="tab:orange", linewidth=1)
  ax[1].fill_between(td[j], τd[j,1]-2eτd[j,1], τd[j,1]+2eτd[j,1], alpha=.2, zorder=3, color="tab:orange", linewidth=0)
  p2, = ax[2].plot(td[j], τeccod[j,1], zorder=3, color="tab:orange", linewidth=1)
end
ax[1].set_ylabel("travel time anomaly (s)")
ax[2].set_ylabel("travel time anomaly (s)")
ax[1].set_title("\$T\$ waves")
ax[2].set_title("ECCO")
ylim = maximum(abs.(ax[1].get_ylim()))
ax[1].set_ylim(-ylim, ylim)
ax[2].set_ylim(-ylim, ylim)
ax[1].set_xlim(td[1]-Day(30), td[end]+Day(30))
ax[1].legend([p1, p2], ["H08", "DGAR"], ncol=3, loc="upper center", frameon=false)
ax[2].legend([p1, p2], ["H08", "DGAR"], ncol=3, loc="upper center", frameon=false)
ax[1].set_xticks(ceil(td[1], Year) : Year(1) : floor(td[end], Year), minor=true)
fig.tight_layout()

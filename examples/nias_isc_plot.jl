using HDF5, LinearAlgebra, PyPlot, NCDatasets, Statistics, Dates, Printf

h = 5e3

# read DGAR kernel
x = h5read("data/temperature/nias_DGAR.h5", "x")
z = h5read("data/temperature/nias_DGAR.h5", "z")
Δx = x[2] - x[1]
Δz = z[2] - z[1]
K = sum(h5read("data/temperature/nias_DGAR.h5", "K"), dims=2)[:,1,:]'*Δx
Km = sum(h5read("data/temperature/nias_DGAR.h5", "Km"), dims=2)[:,1,:]'*Δx

# load DGAR time series
t = DateTime(2000, 1, 1) .+ Millisecond.(h5read("results/nias_isc_DGAR.h5", "t"))
τ = h5read("results/nias_isc_DGAR.h5", "τ")
τargo = h5read("results/nias_isc_DGAR.h5", "τargo")
τecco = h5read("results/nias_isc_DGAR.h5", "τecco")
eτ = h5read("results/nias_isc_DGAR.h5", "eτ")

# load DGAR pairs
tpairs = CSV.read("results/nias_isc_DGAR_tpairs.csv", DataFrame)

# gaps
gaps = [Date(2004, 12, 1), Date(2007, 3, 1), Date(2007, 11, 1), Date(2018, 1, 1)]
ir = [findfirst(t .> gaps[i]) : findlast(t .< gaps[i+1]) for i = 1:length(gaps)-1]

# plot DGAR time series
fig, ax = subplots(1, 1, sharex=true, figsize=(8.4, 4.2))
ax.scatter(t, τ[:,1], s=2, zorder=2, color="tab:blue")
for i = 1:size(tpairs, 1)
  i1 = findfirst(t .== tpairs.event1[i])
  i2 = findfirst(t .== tpairs.event2[i])
  plot([t[i1], t[i2]], [τ[i1], τ[i2]]; color="black", alpha=0.2, zorder=0, linewidth=.75)
end
for j = ir
  #global p1, p2, p3
  #p1, = ax.plot(t[j], τ[j,1], zorder=3, color="tab:blue", linewidth=1)
  ax.fill_between(t[j], τ[j,1]-2eτ[j,1], τ[j,1]+2eτ[j,1], alpha=0, zorder=3, color="tab:blue", linewidth=0)
end
#  ax.plot(t, 1e3ctrend8[:,i], zorder=1, linewidth=1, color="tab:blue")
#  ax.fill_between(t8, 1e3*(ctrend8[:,i]-2ectrend8[:,i]), 1e3*(ctrend8[:,i]+2ectrend8[:,i]), alpha=.2, zorder=0, color="tab:blue", linewidth=0)
#  ax.plot(t8, 1e3ctrendargo8[:,i], zorder=0, linewidth=1, color="tab:orange")
#  ax.fill_between(t8, 1e3*(ctrendargo8[:,i]-2ectrendargo8[:,i]), 1e3*(ctrendargo8[:,i]+2ectrendargo8[:,i]), alpha=.2, zorder=0, color="tab:orange", linewidth=0)
#  ax.plot(t8, 1e3ctrendecco8[:,i], zorder=0, linewidth=1, color="tab:green")
#  ax.fill_between(t8, 1e3*(ctrendecco8[:,i]-2ectrendecco8[:,i]), 1e3*(ctrendecco8[:,i]+2ectrendecco8[:,i]), alpha=.2, zorder=0, color="tab:green", linewidth=0)
ax.set_ylabel("travel time anomaly (s)")
#  ax.set_title("projection onto the $i$(ending(i)) singular vector")
ylim = maximum(abs.(ax.get_ylim()))
ax.set_ylim(ylim, -ylim)
ax.set_xlim(t[1]-Day(30), t[end]+Day(30))
axs = ax.twinx()
axs.set_ylim(1e3ylim/(sum(Km[1,:])*Δz), -1e3ylim/(sum(Km[1,:])*Δz))
axs.set_ylabel("temperature anomaly (mK)")
ax.set_xticks(ceil(t[1], Year) : Year(1) : floor(t[end], Year), minor=true)
fig.align_ylabels()
fig.tight_layout()
fig.savefig("../figures/pres/nias_isc_DGAR_timeseries-1.pdf", dpi=300)

# plot DGAR time series
fig, ax = subplots(1, 1, sharex=true, figsize=(8.4, 4.2))
ax.scatter(t, τ[:,1], s=2, zorder=2, color="tab:blue")
for i = 1:size(tpairs, 1)
  i1 = findfirst(t .== tpairs.event1[i])
  i2 = findfirst(t .== tpairs.event2[i])
  plot([t[i1], t[i2]], [τ[i1], τ[i2]]; color="black", alpha=0.2, zorder=0, linewidth=.75)
end
for j = ir
  global p1, p2, p3
  p1, = ax.plot(t[j], τ[j,1], zorder=3, color="tab:blue", linewidth=1)
  ax.fill_between(t[j], τ[j,1]-2eτ[j,1], τ[j,1]+2eτ[j,1], alpha=0, zorder=3, color="tab:blue", linewidth=0)
end
#  ax.plot(t, 1e3ctrend8[:,i], zorder=1, linewidth=1, color="tab:blue")
#  ax.fill_between(t8, 1e3*(ctrend8[:,i]-2ectrend8[:,i]), 1e3*(ctrend8[:,i]+2ectrend8[:,i]), alpha=.2, zorder=0, color="tab:blue", linewidth=0)
#  ax.plot(t8, 1e3ctrendargo8[:,i], zorder=0, linewidth=1, color="tab:orange")
#  ax.fill_between(t8, 1e3*(ctrendargo8[:,i]-2ectrendargo8[:,i]), 1e3*(ctrendargo8[:,i]+2ectrendargo8[:,i]), alpha=.2, zorder=0, color="tab:orange", linewidth=0)
#  ax.plot(t8, 1e3ctrendecco8[:,i], zorder=0, linewidth=1, color="tab:green")
#  ax.fill_between(t8, 1e3*(ctrendecco8[:,i]-2ectrendecco8[:,i]), 1e3*(ctrendecco8[:,i]+2ectrendecco8[:,i]), alpha=.2, zorder=0, color="tab:green", linewidth=0)
ax.set_ylabel("travel time anomaly (s)")
#  ax.set_title("projection onto the $i$(ending(i)) singular vector")
ylim = maximum(abs.(ax.get_ylim()))
ax.set_ylim(ylim, -ylim)
ax.set_xlim(t[1]-Day(30), t[end]+Day(30))
axs = ax.twinx()
axs.set_ylim(1e3ylim/(sum(Km[1,:])*Δz), -1e3ylim/(sum(Km[1,:])*Δz))
axs.set_ylabel("temperature anomaly (mK)")
ax.set_xticks(ceil(t[1], Year) : Year(1) : floor(t[end], Year), minor=true)
fig.align_ylabels()
fig.tight_layout()
fig.savefig("../figures/pres/nias_isc_DGAR_timeseries-2.pdf", dpi=300)

# plot DGAR time series
fig, ax = subplots(1, 1, sharex=true, figsize=(8.4, 4.2))
ax.scatter(t, τ[:,1], s=2, zorder=2, color="tab:blue")
for j = ir
  global p1, p2, p3
  p1, = ax.plot(t[j], τ[j,1], zorder=3, color="tab:blue", linewidth=1)
  ax.fill_between(t[j], τ[j,1]-2eτ[j,1], τ[j,1]+2eτ[j,1], alpha=.2, zorder=3, color="tab:blue", linewidth=0)
end
#  ax.plot(t, 1e3ctrend8[:,i], zorder=1, linewidth=1, color="tab:blue")
#  ax.fill_between(t8, 1e3*(ctrend8[:,i]-2ectrend8[:,i]), 1e3*(ctrend8[:,i]+2ectrend8[:,i]), alpha=.2, zorder=0, color="tab:blue", linewidth=0)
#  ax.plot(t8, 1e3ctrendargo8[:,i], zorder=0, linewidth=1, color="tab:orange")
#  ax.fill_between(t8, 1e3*(ctrendargo8[:,i]-2ectrendargo8[:,i]), 1e3*(ctrendargo8[:,i]+2ectrendargo8[:,i]), alpha=.2, zorder=0, color="tab:orange", linewidth=0)
#  ax.plot(t8, 1e3ctrendecco8[:,i], zorder=0, linewidth=1, color="tab:green")
#  ax.fill_between(t8, 1e3*(ctrendecco8[:,i]-2ectrendecco8[:,i]), 1e3*(ctrendecco8[:,i]+2ectrendecco8[:,i]), alpha=.2, zorder=0, color="tab:green", linewidth=0)
ax.set_ylabel("travel time anomaly (s)")
#  ax.set_title("projection onto the $i$(ending(i)) singular vector")
ylim = maximum(abs.(ax.get_ylim()))
ax.set_ylim(ylim, -ylim)
ax.set_xlim(t[1]-Day(30), t[end]+Day(30))
axs = ax.twinx()
axs.set_ylim(1e3ylim/(sum(Km[1,:])*Δz), -1e3ylim/(sum(Km[1,:])*Δz))
axs.set_ylabel("temperature anomaly (mK)")
ax.set_xticks(ceil(t[1], Year) : Year(1) : floor(t[end], Year), minor=true)
fig.align_ylabels()
fig.tight_layout()
fig.savefig("../figures/pres/nias_isc_DGAR_timeseries-3.pdf", dpi=300)

# plot DGAR time series
fig, ax = subplots(1, 1, sharex=true, figsize=(8.4, 4.2))
ax.scatter(t, τ[:,1], s=2, zorder=1, color="tab:blue")
for j = ir
  global p1, p2, p3
  p1, = ax.plot(t[j], τ[j,1], zorder=1, color="tab:blue", linewidth=1)
  ax.fill_between(t[j], τ[j,1]-2eτ[j,1], τ[j,1]+2eτ[j,1], alpha=.2, zorder=1, color="tab:blue", linewidth=0)
  p2, = ax.plot(t[j], τargo[j,1], zorder=2, linewidth=1, color="tab:orange")
  p3, = ax.plot(t[j], τecco[j,1], zorder=3, linewidth=1, color="tab:green")
end
#  ax.plot(t, 1e3ctrend8[:,i], zorder=1, linewidth=1, color="tab:blue")
#  ax.fill_between(t8, 1e3*(ctrend8[:,i]-2ectrend8[:,i]), 1e3*(ctrend8[:,i]+2ectrend8[:,i]), alpha=.2, zorder=0, color="tab:blue", linewidth=0)
#  ax.plot(t8, 1e3ctrendargo8[:,i], zorder=0, linewidth=1, color="tab:orange")
#  ax.fill_between(t8, 1e3*(ctrendargo8[:,i]-2ectrendargo8[:,i]), 1e3*(ctrendargo8[:,i]+2ectrendargo8[:,i]), alpha=.2, zorder=0, color="tab:orange", linewidth=0)
#  ax.plot(t8, 1e3ctrendecco8[:,i], zorder=0, linewidth=1, color="tab:green")
#  ax.fill_between(t8, 1e3*(ctrendecco8[:,i]-2ectrendecco8[:,i]), 1e3*(ctrendecco8[:,i]+2ectrendecco8[:,i]), alpha=.2, zorder=0, color="tab:green", linewidth=0)
ax.set_ylabel("travel time anomaly (s)")
#  ax.set_title("projection onto the $i$(ending(i)) singular vector")
ylim = maximum(abs.(ax.get_ylim()))
ax.set_ylim(ylim, -ylim)
ax.set_xlim(t[1]-Day(30), t[end]+Day(30))
axs = ax.twinx()
axs.set_ylim(1e3ylim/(sum(Km[1,:])*Δz), -1e3ylim/(sum(Km[1,:])*Δz))
axs.set_ylabel("temperature anomaly (mK)")
ax.legend([p1, p2, p3], ["\$T\$ waves", "Argo (SIO)", "ECCO"], ncol=3, loc="lower center", frameon=false)
ax.set_xticks(ceil(t[1], Year) : Year(1) : floor(t[end], Year), minor=true)
fig.align_ylabels()
fig.tight_layout()
fig.savefig("../figures/pres/nias_isc_DGAR_timeseries-4.pdf", dpi=300)

xlim1 = [Date(2005, 3, 15), Date(2006, 9, 15)]
ylim1 = [0.5, -0.35]

# plot DGAR time series
fig, ax = subplots(1, 1, sharex=true, figsize=(8.4, 4.2))
ax.scatter(t, τ[:,1], s=2, zorder=1, color="tab:blue")
for j = ir
  global p1, p2, p3
  p1, = ax.plot(t[j], τ[j,1], zorder=1, color="tab:blue", linewidth=1)
  ax.fill_between(t[j], τ[j,1]-2eτ[j,1], τ[j,1]+2eτ[j,1], alpha=.2, zorder=1, color="tab:blue", linewidth=0)
  p2, = ax.plot(t[j], τargo[j,1], zorder=2, linewidth=1, color="tab:orange")
  p3, = ax.plot(t[j], τecco[j,1], zorder=3, linewidth=1, color="tab:green")
end
#  ax.plot(t, 1e3ctrend8[:,i], zorder=1, linewidth=1, color="tab:blue")
#  ax.fill_between(t8, 1e3*(ctrend8[:,i]-2ectrend8[:,i]), 1e3*(ctrend8[:,i]+2ectrend8[:,i]), alpha=.2, zorder=0, color="tab:blue", linewidth=0)
#  ax.plot(t8, 1e3ctrendargo8[:,i], zorder=0, linewidth=1, color="tab:orange")
#  ax.fill_between(t8, 1e3*(ctrendargo8[:,i]-2ectrendargo8[:,i]), 1e3*(ctrendargo8[:,i]+2ectrendargo8[:,i]), alpha=.2, zorder=0, color="tab:orange", linewidth=0)
#  ax.plot(t8, 1e3ctrendecco8[:,i], zorder=0, linewidth=1, color="tab:green")
#  ax.fill_between(t8, 1e3*(ctrendecco8[:,i]-2ectrendecco8[:,i]), 1e3*(ctrendecco8[:,i]+2ectrendecco8[:,i]), alpha=.2, zorder=0, color="tab:green", linewidth=0)
ax.set_ylabel("travel time anomaly (s)")
#  ax.set_title("projection onto the $i$(ending(i)) singular vector")
ylim = maximum(abs.(ax.get_ylim()))
ax.plot([xlim1[1], xlim1[2], xlim1[2], xlim1[1], xlim1[1]], [ylim1[1], ylim1[1], ylim1[2], ylim1[2], ylim1[1]]; color="black", linewidth=0.8)
ax.set_ylim(ylim, -ylim)
ax.set_xlim(t[1]-Day(30), t[end]+Day(30))
axs = ax.twinx()
axs.set_ylim(1e3ylim/(sum(Km[1,:])*Δz), -1e3ylim/(sum(Km[1,:])*Δz))
axs.set_ylabel("temperature anomaly (mK)")
ax.legend([p1, p2, p3], ["\$T\$ waves", "Argo (SIO)", "ECCO"], ncol=3, loc="lower center", frameon=false)
ax.set_xticks(ceil(t[1], Year) : Year(1) : floor(t[end], Year), minor=true)
fig.align_ylabels()
fig.tight_layout()
fig.savefig("../figures/pres/nias_isc_DGAR_timeseries-5.pdf", dpi=300)

# plot DGAR time series
fig, ax = subplots(1, 1, sharex=true, figsize=(8.4, 4.2))
ax.scatter(t, τ[:,1], s=2, zorder=1, color="tab:blue")
for j = ir
  global p1, p2, p3
  p1, = ax.plot(t[j], τ[j,1], zorder=1, color="tab:blue", linewidth=1)
  ax.fill_between(t[j], τ[j,1]-2eτ[j,1], τ[j,1]+2eτ[j,1], alpha=.2, zorder=1, color="tab:blue", linewidth=0)
  p2, = ax.plot(t[j], τargo[j,1], zorder=2, linewidth=1, color="tab:orange")
  p3, = ax.plot(t[j], τecco[j,1], zorder=3, linewidth=1, color="tab:green")
end
#  ax.plot(t, 1e3ctrend8[:,i], zorder=1, linewidth=1, color="tab:blue")
#  ax.fill_between(t8, 1e3*(ctrend8[:,i]-2ectrend8[:,i]), 1e3*(ctrend8[:,i]+2ectrend8[:,i]), alpha=.2, zorder=0, color="tab:blue", linewidth=0)
#  ax.plot(t8, 1e3ctrendargo8[:,i], zorder=0, linewidth=1, color="tab:orange")
#  ax.fill_between(t8, 1e3*(ctrendargo8[:,i]-2ectrendargo8[:,i]), 1e3*(ctrendargo8[:,i]+2ectrendargo8[:,i]), alpha=.2, zorder=0, color="tab:orange", linewidth=0)
#  ax.plot(t8, 1e3ctrendecco8[:,i], zorder=0, linewidth=1, color="tab:green")
#  ax.fill_between(t8, 1e3*(ctrendecco8[:,i]-2ectrendecco8[:,i]), 1e3*(ctrendecco8[:,i]+2ectrendecco8[:,i]), alpha=.2, zorder=0, color="tab:green", linewidth=0)
ax.set_ylabel("travel time anomaly (s)")
#  ax.set_title("projection onto the $i$(ending(i)) singular vector")
ax.set_ylim(ylim1...)
ax.set_xlim(xlim1...)
ax.set_yticks(0.4 : -0.2 : ylim1[2])
axs = ax.twinx()
axs.set_ylim((1e3ylim1/(sum(Km[1,:])*Δz))...)
axs.set_ylabel("temperature anomaly (mK)")
ax.legend([p1, p2, p3], ["\$T\$ waves", "Argo (SIO)", "ECCO"], ncol=3, loc="lower center", frameon=false)
ax.set_xticks(Date(2005, 4, 1) : Month(2) : xlim1[2])
ax.set_xticks(ceil(xlim1[1], Month) : Month(1) : xlim1[2]; minor=true)
ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y-%m"))
fig.align_ylabels()
fig.tight_layout()
fig.savefig("../figures/pres/nias_isc_DGAR_timeseries-6.pdf", dpi=300)

xlim2 = [Date(2005, 3, 28), Date(2005, 4, 25)]
ylim2 = [0.12, -0.275]

# plot DGAR time series
fig, ax = subplots(1, 1, sharex=true, figsize=(8.4, 4.2))
ax.scatter(t, τ[:,1], s=2, zorder=1, color="tab:blue")
for j = ir
  global p1, p2, p3
  p1, = ax.plot(t[j], τ[j,1], zorder=1, color="tab:blue", linewidth=1)
  ax.fill_between(t[j], τ[j,1]-2eτ[j,1], τ[j,1]+2eτ[j,1], alpha=.2, zorder=1, color="tab:blue", linewidth=0)
  p2, = ax.plot(t[j], τargo[j,1], zorder=2, linewidth=1, color="tab:orange")
  p3, = ax.plot(t[j], τecco[j,1], zorder=3, linewidth=1, color="tab:green")
end
#  ax.plot(t, 1e3ctrend8[:,i], zorder=1, linewidth=1, color="tab:blue")
#  ax.fill_between(t8, 1e3*(ctrend8[:,i]-2ectrend8[:,i]), 1e3*(ctrend8[:,i]+2ectrend8[:,i]), alpha=.2, zorder=0, color="tab:blue", linewidth=0)
#  ax.plot(t8, 1e3ctrendargo8[:,i], zorder=0, linewidth=1, color="tab:orange")
#  ax.fill_between(t8, 1e3*(ctrendargo8[:,i]-2ectrendargo8[:,i]), 1e3*(ctrendargo8[:,i]+2ectrendargo8[:,i]), alpha=.2, zorder=0, color="tab:orange", linewidth=0)
#  ax.plot(t8, 1e3ctrendecco8[:,i], zorder=0, linewidth=1, color="tab:green")
#  ax.fill_between(t8, 1e3*(ctrendecco8[:,i]-2ectrendecco8[:,i]), 1e3*(ctrendecco8[:,i]+2ectrendecco8[:,i]), alpha=.2, zorder=0, color="tab:green", linewidth=0)
ax.set_ylabel("travel time anomaly (s)")
#  ax.set_title("projection onto the $i$(ending(i)) singular vector")
ax.plot([xlim2[1], xlim2[2], xlim2[2], xlim2[1], xlim2[1]], [ylim2[1], ylim2[1], ylim2[2], ylim2[2], ylim2[1]]; color="black", linewidth=0.8)
ax.set_ylim(ylim1...)
ax.set_xlim(xlim1...)
ax.set_yticks(0.4 : -0.2 : ylim1[2])
axs = ax.twinx()
axs.set_ylim((1e3ylim1/(sum(Km[1,:])*Δz))...)
axs.set_ylabel("temperature anomaly (mK)")
ax.legend([p1, p2, p3], ["\$T\$ waves", "Argo (SIO)", "ECCO"], ncol=3, loc="lower center", frameon=false)
ax.set_xticks(Date(2005, 4, 1) : Month(2) : xlim1[2])
ax.set_xticks(ceil(xlim1[1], Month) : Month(1) : xlim1[2]; minor=true)
ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y-%m"))
fig.align_ylabels()
fig.tight_layout()
fig.savefig("../figures/pres/nias_isc_DGAR_timeseries-7.pdf", dpi=300)

# plot DGAR time series
fig, ax = subplots(1, 1, sharex=true, figsize=(8.4, 4.2))
ax.scatter(t, τ[:,1], s=2, zorder=1, color="tab:blue")
for j = ir
  global p1, p2, p3
  p1, = ax.plot(t[j], τ[j,1], zorder=1, color="tab:blue", linewidth=1)
  ax.fill_between(t[j], τ[j,1]-2eτ[j,1], τ[j,1]+2eτ[j,1], alpha=.2, zorder=1, color="tab:blue", linewidth=0)
  p2, = ax.plot(t[j], τargo[j,1], zorder=2, linewidth=1, color="tab:orange")
  p3, = ax.plot(t[j], τecco[j,1], zorder=3, linewidth=1, color="tab:green")
end
#  ax.plot(t, 1e3ctrend8[:,i], zorder=1, linewidth=1, color="tab:blue")
#  ax.fill_between(t8, 1e3*(ctrend8[:,i]-2ectrend8[:,i]), 1e3*(ctrend8[:,i]+2ectrend8[:,i]), alpha=.2, zorder=0, color="tab:blue", linewidth=0)
#  ax.plot(t8, 1e3ctrendargo8[:,i], zorder=0, linewidth=1, color="tab:orange")
#  ax.fill_between(t8, 1e3*(ctrendargo8[:,i]-2ectrendargo8[:,i]), 1e3*(ctrendargo8[:,i]+2ectrendargo8[:,i]), alpha=.2, zorder=0, color="tab:orange", linewidth=0)
#  ax.plot(t8, 1e3ctrendecco8[:,i], zorder=0, linewidth=1, color="tab:green")
#  ax.fill_between(t8, 1e3*(ctrendecco8[:,i]-2ectrendecco8[:,i]), 1e3*(ctrendecco8[:,i]+2ectrendecco8[:,i]), alpha=.2, zorder=0, color="tab:green", linewidth=0)
ax.set_ylabel("travel time anomaly (s)")
#  ax.set_title("projection onto the $i$(ending(i)) singular vector")
ax.set_ylim(ylim2...)
ax.set_xlim(xlim2...)
ax.set_yticks(round(ylim2[1]; sigdigits=1) : -0.1 : ylim2[2])
axs = ax.twinx()
axs.set_ylim((1e3ylim2/(sum(Km[1,:])*Δz))...)
axs.set_ylabel("temperature anomaly (mK)")
ax.legend([p1, p2, p3], ["\$T\$ waves", "Argo (SIO)", "ECCO"], ncol=3, loc="lower center", frameon=false)
ax.set_xticks(Date(2005, 3, 28) : Week(1) : xlim2[2])
ax.set_xticks(xlim2[1] : Day(1) : xlim2[2]; minor=true)
fig.align_ylabels()
fig.tight_layout()
fig.savefig("../figures/pres/nias_isc_DGAR_timeseries-8.pdf", dpi=300)

# read H08 kernel
x = h5read("data/temperature/nias_H08.h5", "x")
z = h5read("data/temperature/nias_H08.h5", "z")
Δx = x[2] - x[1]
Δz = z[2] - z[1]
K8 = sum(h5read("data/temperature/nias_H08.h5", "K"), dims=2)[:,1,:]'*Δx
Km8 = sum(h5read("data/temperature/nias_H08.h5", "Km"), dims=2)[:,1,:]'*Δx

# read H01 kernel
x = h5read("data/temperature/nias_H01.h5", "x")
z = h5read("data/temperature/nias_H01.h5", "z")
Δx = x[2] - x[1]
Δz = z[2] - z[1]
K1 = sum(h5read("data/temperature/nias_H01.h5", "K"), dims=2)[:,1,:]'*Δx
Km1 = sum(h5read("data/temperature/nias_H01.h5", "Km"), dims=2)[:,1,:]'*Δx

# SVD
U8, Λ8, V8 = svd(K8)
U1, Λ1, V1 = svd(K1)

# rescale to be independent of resolution (∫v'vdz = h)
V8 /= sqrt(Δz/h)
V1 /= sqrt(Δz/h)
Λ8 *= sqrt(Δz/h)
Λ1 *= sqrt(Δz/h)

# load H08 time series
t8 = DateTime(2000, 1, 1) .+ Millisecond.(h5read("results/nias_tm_H08.h5", "t"))
τ8 = h5read("results/nias_tm_H08.h5", "τ")
τ8argo = h5read("results/nias_tm_H08.h5", "τargo")
τ8ecco = h5read("results/nias_tm_H08.h5", "τecco")
eτ8 = h5read("results/nias_tm_H08.h5", "eτ")
c8 = h5read("results/nias_tm_H08.h5", "c")
c8argo = h5read("results/nias_tm_H08.h5", "cargo")
c8ecco = h5read("results/nias_tm_H08.h5", "cecco")
ec8 = h5read("results/nias_tm_H08.h5", "ec")

# gaps
gaps8 = [Date(2004, 1, 1), Date(2007, 1, 1), Date(2011, 1, 1), Date(2020, 1, 1)]
ir8 = [findfirst(t8 .> gaps8[i]) : findlast(t8 .< gaps8[i+1]) for i = 1:length(gaps8)-1]

function ending(i)
  if i == 1
    return "st"
  elseif i == 2
    return "nd"
  elseif i == 3
    return "rd"
  else
    return "th"
  end
end

# plot DGAG–H08 comparison
fig, ax = subplots(1, 1, sharex=true, figsize=(8.4, 4.2))
ax.scatter(t, τ[:,1], s=2, zorder=1, color="tab:blue")
ax.scatter(t8, τ8[:,1], s=2, zorder=2, color="tab:red")
for j = ir
  global p1
  p1, = ax.plot(t[j], τ[j,1], zorder=1, color="tab:blue", linewidth=1)
  ax.fill_between(t[j], τ[j,1]-2eτ[j,1], τ[j,1]+2eτ[j,1], alpha=.2, zorder=1, color="tab:blue", linewidth=0)
end
for j = ir8
  global p2
  p2, = ax.plot(t8[j], τ8[j,1], zorder=1, color="tab:red", linewidth=1)
  ax.fill_between(t8[j], τ8[j,1]-2eτ8[j,1], τ8[j,1]+2eτ8[j,1], alpha=.2, zorder=1, color="tab:red", linewidth=0)
end
#  ax.plot(t, 1e3ctrend8[:,i], zorder=1, linewidth=1, color="tab:blue")
#  ax.fill_between(t8, 1e3*(ctrend8[:,i]-2ectrend8[:,i]), 1e3*(ctrend8[:,i]+2ectrend8[:,i]), alpha=.2, zorder=0, color="tab:blue", linewidth=0)
#  ax.plot(t8, 1e3ctrendargo8[:,i], zorder=0, linewidth=1, color="tab:orange")
#  ax.fill_between(t8, 1e3*(ctrendargo8[:,i]-2ectrendargo8[:,i]), 1e3*(ctrendargo8[:,i]+2ectrendargo8[:,i]), alpha=.2, zorder=0, color="tab:orange", linewidth=0)
#  ax.plot(t8, 1e3ctrendecco8[:,i], zorder=0, linewidth=1, color="tab:green")
#  ax.fill_between(t8, 1e3*(ctrendecco8[:,i]-2ectrendecco8[:,i]), 1e3*(ctrendecco8[:,i]+2ectrendecco8[:,i]), alpha=.2, zorder=0, color="tab:green", linewidth=0)
ax.set_ylabel("travel time anomaly (s)")
ax.legend([p1, p2], ["DGAR", "CTBTO H08"], ncol=2, loc="lower center", frameon=false)
#  ax.set_title("projection onto the $i$(ending(i)) singular vector")
ylim = maximum(abs.(ax.get_ylim()))
ax.set_ylim(ylim, -ylim)
ax.set_xlim(t[1]-Day(30), t[end]+Day(30))
ax.set_xticks(ceil(t[1], Year) : Year(1) : floor(t8[end], Year), minor=true)
fig.align_ylabels()
fig.tight_layout()
fig.savefig("../figures/pres/nias_DGAR_H08_comparison.pdf", dpi=300)

# plot H08 time series
fig, ax = subplots(2, 1, sharex=true, figsize=(8.4, 4.6))
for i = 1:2
  ax[i].scatter(t8, 1e3c8[:,i], s=2, zorder=2, color="tab:blue")
  for j = ir8
    global p1, p2, p3
    p1, = ax[i].plot(t8[j], 1e3c8[j,i], zorder=3, color="tab:blue", linewidth=1)
    ax[i].fill_between(t8[j], 1e3*(c8[j,i]-2ec8[j,i]), 1e3*(c8[j,i]+2ec8[j,i]), alpha=.2, zorder=3, color="tab:blue", linewidth=0)
#    p2, = ax[i].plot(t8[j], 1e3c8argo[j,i], zorder=1, linewidth=1, color="tab:orange")
#    p3, = ax[i].plot(t8[j], 1e3c8ecco[j,i], zorder=2, linewidth=1, color="tab:green")
  end
#  ax[i].plot(t8, 1e3ctrend8[:,i], zorder=1, linewidth=1, color="tab:blue")
#  ax[i].fill_between(t8, 1e3*(ctrend8[:,i]-2ectrend8[:,i]), 1e3*(ctrend8[:,i]+2ectrend8[:,i]), alpha=.2, zorder=0, color="tab:blue", linewidth=0)
#  ax[i].plot(t8, 1e3ctrendargo8[:,i], zorder=0, linewidth=1, color="tab:orange")
#  ax[i].fill_between(t8, 1e3*(ctrendargo8[:,i]-2ectrendargo8[:,i]), 1e3*(ctrendargo8[:,i]+2ectrendargo8[:,i]), alpha=.2, zorder=0, color="tab:orange", linewidth=0)
#  ax[i].plot(t8, 1e3ctrendecco8[:,i], zorder=0, linewidth=1, color="tab:green")
#  ax[i].fill_between(t8, 1e3*(ctrendecco8[:,i]-2ectrendecco8[:,i]), 1e3*(ctrendecco8[:,i]+2ectrendecco8[:,i]), alpha=.2, zorder=0, color="tab:green", linewidth=0)
  ax[i].set_ylabel("temperature anomaly (mK)")
  ax[i].set_title("$i$(ending(i)) singular vector")
  ylim = maximum(abs.(ax[i].get_ylim()))
  ax[i].set_ylim(-ylim, ylim)
  ax[i].set_xlim(t8[1]-Day(30), t8[end]+Day(30))
end
#ax[1].legend([p1, p2, p3], ["\$T\$ waves", "Argo (SIO)", "ECCO"], ncol=3, loc="lower center", frameon=false)
ax[1].set_xticks(ceil(t8[1], Year) : Year(1) : floor(t8[end], Year), minor=true)
fig.align_ylabels()
fig.tight_layout()
savefig("../figures/pres/nias_tm_H08_timeseries-1.pdf", dpi=300)

# plot H08 time series
fig, ax = subplots(2, 1, sharex=true, figsize=(8.4, 4.6))
for i = 1:2
  ax[i].scatter(t8, 1e3c8[:,i], s=2, zorder=1, color="tab:blue")
  for j = ir8
    global p1, p2, p3
    p1, = ax[i].plot(t8[j], 1e3c8[j,i], zorder=1, color="tab:blue", linewidth=1)
    ax[i].fill_between(t8[j], 1e3*(c8[j,i]-2ec8[j,i]), 1e3*(c8[j,i]+2ec8[j,i]), alpha=.2, zorder=1, color="tab:blue", linewidth=0)
    p2, = ax[i].plot(t8[j], 1e3c8argo[j,i], zorder=2, linewidth=1, color="tab:orange")
    p3, = ax[i].plot(t8[j], 1e3c8ecco[j,i], zorder=3, linewidth=1, color="tab:green")
  end
#  ax[i].plot(t8, 1e3ctrend8[:,i], zorder=1, linewidth=1, color="tab:blue")
#  ax[i].fill_between(t8, 1e3*(ctrend8[:,i]-2ectrend8[:,i]), 1e3*(ctrend8[:,i]+2ectrend8[:,i]), alpha=.2, zorder=0, color="tab:blue", linewidth=0)
#  ax[i].plot(t8, 1e3ctrendargo8[:,i], zorder=0, linewidth=1, color="tab:orange")
#  ax[i].fill_between(t8, 1e3*(ctrendargo8[:,i]-2ectrendargo8[:,i]), 1e3*(ctrendargo8[:,i]+2ectrendargo8[:,i]), alpha=.2, zorder=0, color="tab:orange", linewidth=0)
#  ax[i].plot(t8, 1e3ctrendecco8[:,i], zorder=0, linewidth=1, color="tab:green")
#  ax[i].fill_between(t8, 1e3*(ctrendecco8[:,i]-2ectrendecco8[:,i]), 1e3*(ctrendecco8[:,i]+2ectrendecco8[:,i]), alpha=.2, zorder=0, color="tab:green", linewidth=0)
  ax[i].set_ylabel("temperature anomaly (mK)")
  ax[i].set_title("$i$(ending(i)) singular vector")
  ylim = maximum(abs.(ax[i].get_ylim()))
  ax[i].set_ylim(-ylim, ylim)
  ax[i].set_xlim(t8[1]-Day(30), t8[end]+Day(30))
end
ax[1].legend([p1, p2, p3], ["\$T\$ waves", "Argo (SIO)", "ECCO"], ncol=3, loc="lower center", frameon=false)
ax[1].set_xticks(ceil(t8[1], Year) : Year(1) : floor(t8[end], Year), minor=true)
fig.align_ylabels()
fig.tight_layout()
fig.savefig("../figures/pres/nias_tm_H08_timeseries-2.pdf", dpi=300)

# plot H08 time series
fig, ax = subplots(1, 1, sharex=true, figsize=(8.4, 3.8))
for j = ir8
  global img
  img = ax.pcolormesh(t8[j], -1e-3z, 1e3V8[:,1:2]*c8[j,1:2]', cmap="RdBu_r", vmin=-170, vmax=170, shading="nearest", rasterized="true")
end
colorbar(img, label="temperature anomaly (mK)")
ax.set_ylim(5, 0)
ax.set_ylabel("depth (km)")
ax.set_title("temperature reconstruction")
ax.set_xticks(ceil(t8[1], Year) : Year(1) : floor(t8[end], Year), minor=true)
fig.align_ylabels()
fig.tight_layout()
fig.savefig("../figures/pres/nias_tm_H08_timeseries-3.pdf", dpi=300)

# plot H08 time series
fig, ax = subplots(1, 1, sharex=true, figsize=(8.4, 3.8))
for j = ir8
  global img
  img = ax.pcolormesh(t8[j], -1e-3z, 1e3V8[:,1:2]*c8[j,1:2]', cmap="RdBu_r", vmin=-170, vmax=170, shading="nearest", rasterized="true")
end
colorbar(img, label="temperature anomaly (mK)")
ax.set_ylim(5, 0)
ax.set_ylabel("depth (km)")
ax.set_title("temperature reconstruction")
#ax.set_xticks(ceil(t8[1], Year) : Year(1) : floor(t8[end], Year), minor=true)
ax.set_xlim(xlim1...)
ax.set_xticks(Date(2005, 4, 1) : Month(2) : xlim1[2])
ax.set_xticks(ceil(xlim1[1], Month) : Month(1) : xlim1[2]; minor=true)
ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y-%m"))
fig.align_ylabels()
fig.tight_layout()
fig.savefig("../figures/pres/nias_tm_H08_timeseries-4.pdf", dpi=300)

# load H01 time series
t1 = DateTime(2000, 1, 1) .+ Millisecond.(h5read("results/nias_tm_H01.h5", "t"))
τ1 = h5read("results/nias_tm_H01.h5", "τ")
τ1argo = h5read("results/nias_tm_H01.h5", "τargo")
τ1ecco = h5read("results/nias_tm_H01.h5", "τecco")
eτ1 = h5read("results/nias_tm_H01.h5", "eτ")
c1 = h5read("results/nias_tm_H01.h5", "c")
c1argo = h5read("results/nias_tm_H01.h5", "cargo")
c1ecco = h5read("results/nias_tm_H01.h5", "cecco")
ec1 = h5read("results/nias_tm_H01.h5", "ec")

# gaps
gaps1 = [Date(2004, 1, 1), Date(2020, 1, 1)]
ir1 = [findfirst(t1 .> gaps1[i]) : findlast(t1 .< gaps1[i+1]) for i = 1:length(gaps1)-1]

# plot H01 time series
fig, ax = subplots(1, 1, sharex=true, figsize=(8.4, 4.2))
ax.scatter(t1, τ1[:,1], s=2, zorder=2, color="tab:blue")
for j = ir1
  global p1, p2, p3
  p1, = ax.plot(t1[j], τ1[j,1], zorder=3, color="tab:blue", linewidth=1)
  ax.fill_between(t1[j], τ1[j,1]-2eτ1[j,1], τ1[j,1]+2eτ1[j,1], alpha=.2, zorder=3, color="tab:blue", linewidth=0)
end
#  ax.plot(t, 1e3ctrend8[:,i], zorder=1, linewidth=1, color="tab:blue")
#  ax.fill_between(t8, 1e3*(ctrend8[:,i]-2ectrend8[:,i]), 1e3*(ctrend8[:,i]+2ectrend8[:,i]), alpha=.2, zorder=0, color="tab:blue", linewidth=0)
#  ax.plot(t8, 1e3ctrendargo8[:,i], zorder=0, linewidth=1, color="tab:orange")
#  ax.fill_between(t8, 1e3*(ctrendargo8[:,i]-2ectrendargo8[:,i]), 1e3*(ctrendargo8[:,i]+2ectrendargo8[:,i]), alpha=.2, zorder=0, color="tab:orange", linewidth=0)
#  ax.plot(t8, 1e3ctrendecco8[:,i], zorder=0, linewidth=1, color="tab:green")
#  ax.fill_between(t8, 1e3*(ctrendecco8[:,i]-2ectrendecco8[:,i]), 1e3*(ctrendecco8[:,i]+2ectrendecco8[:,i]), alpha=.2, zorder=0, color="tab:green", linewidth=0)
ax.set_ylabel("travel time anomaly (s)")
#  ax.set_title("projection onto the $i$(ending(i)) singular vector")
ylim = maximum(abs.(ax.get_ylim()))
ax.set_ylim(ylim, -ylim)
ax.set_yticks(1 : -0.5 : -1)
ax.set_xlim(t1[1]-Day(30), t1[end]+Day(30))
axs = ax.twinx()
axs.set_ylim(1e3ylim/(sum(Km1[1,:])*Δz), -1e3ylim/(sum(Km1[1,:])*Δz))
axs.set_ylabel("temperature anomaly (mK)")
ax.set_xticks(ceil(t[1], Year) : Year(1) : floor(t1[end], Year), minor=true)
fig.align_ylabels()
fig.tight_layout()
fig.savefig("../figures/pres/nias_tm_H01_timeseries-1.pdf", dpi=300)

# plot H01 time series
fig, ax = subplots(1, 1, sharex=true, figsize=(8.4, 4.2))
ax.scatter(t1, τ1[:,1], s=2, zorder=1, color="tab:blue")
for j = ir1
  global p1, p2, p3
  p1, = ax.plot(t1[j], τ1[j,1], zorder=1, color="tab:blue", linewidth=1)
  ax.fill_between(t1[j], τ1[j,1]-2eτ1[j,1], τ1[j,1]+2eτ1[j,1], alpha=.2, zorder=1, color="tab:blue", linewidth=0)
  p2, = ax.plot(t1[j], τ1argo[j,1], zorder=2, linewidth=1, color="tab:orange")
  p3, = ax.plot(t1[j], τ1ecco[j,1], zorder=3, linewidth=1, color="tab:green")
end
#  ax.plot(t, 1e3ctrend8[:,i], zorder=1, linewidth=1, color="tab:blue")
#  ax.fill_between(t8, 1e3*(ctrend8[:,i]-2ectrend8[:,i]), 1e3*(ctrend8[:,i]+2ectrend8[:,i]), alpha=.2, zorder=0, color="tab:blue", linewidth=0)
#  ax.plot(t8, 1e3ctrendargo8[:,i], zorder=0, linewidth=1, color="tab:orange")
#  ax.fill_between(t8, 1e3*(ctrendargo8[:,i]-2ectrendargo8[:,i]), 1e3*(ctrendargo8[:,i]+2ectrendargo8[:,i]), alpha=.2, zorder=0, color="tab:orange", linewidth=0)
#  ax.plot(t8, 1e3ctrendecco8[:,i], zorder=0, linewidth=1, color="tab:green")
#  ax.fill_between(t8, 1e3*(ctrendecco8[:,i]-2ectrendecco8[:,i]), 1e3*(ctrendecco8[:,i]+2ectrendecco8[:,i]), alpha=.2, zorder=0, color="tab:green", linewidth=0)
ax.set_ylabel("travel time anomaly (s)")
#  ax.set_title("projection onto the $i$(ending(i)) singular vector")
ylim = maximum(abs.(ax.get_ylim()))
ax.set_ylim(ylim, -ylim)
ax.set_yticks(1 : -0.5 : -1)
ax.set_xlim(t[1]-Day(30), t[end]+Day(30))
axs = ax.twinx()
axs.set_ylim(1e3ylim/(sum(Km1[1,:])*Δz), -1e3ylim/(sum(Km1[1,:])*Δz))
axs.set_ylabel("temperature anomaly (mK)")
ax.legend([p1, p2, p3], ["\$T\$ waves", "Argo (SIO)", "ECCO"], ncol=3, loc="lower center", frameon=false)
ax.set_xticks(ceil(t1[1], Year) : Year(1) : floor(t1[end], Year), minor=true)
fig.align_ylabels()
fig.tight_layout()
fig.savefig("../figures/pres/nias_tm_H01_timeseries-2.pdf", dpi=300)

# plot H08 time series
fig, ax = subplots(1, 1, sharex=true, figsize=(8.4, 3.8))
for j = ir1
  global img
  img = ax.pcolormesh(t1[j], -1e-3z, 1e3V1[:,1:2]*c1[j,1:2]', cmap="RdBu_r", vmin=-200, vmax=200, shading="nearest", rasterized="true")
end
colorbar(img, label="temperature anomaly (mK)")
ax.set_ylim(5, 0)
ax.set_ylabel("depth (km)")
ax.set_title("temperature reconstruction")
ax.set_xticks(ceil(t1[1], Year) : Year(1) : floor(t1[end], Year), minor=true)
ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y-%m"))
fig.align_ylabels()
fig.tight_layout()
fig.savefig("../figures/pres/nias_tm_H01_timeseries-3.pdf", dpi=300)

stop

ctrend8 = h5read("results/nias_tm_H08.h5", "ctrend")
ctrend1 = h5read("results/nias_tm_H01.h5", "ctrend")
ctrendargo8 = h5read("results/nias_tm_H08.h5", "ctrendargo")
ctrendargo1 = h5read("results/nias_tm_H01.h5", "ctrendargo")
ctrendecco8 = h5read("results/nias_tm_H08.h5", "ctrendecco")
ctrendecco1 = h5read("results/nias_tm_H01.h5", "ctrendecco")
ectrend8 = h5read("results/nias_tm_H08.h5", "ectrend")
ectrend1 = h5read("results/nias_tm_H01.h5", "ectrend")
ectrendargo8 = h5read("results/nias_tm_H08.h5", "ectrendargo")
ectrendargo1 = h5read("results/nias_tm_H01.h5", "ectrendargo")
ectrendecco8 = h5read("results/nias_tm_H08.h5", "ectrendecco")
ectrendecco1 = h5read("results/nias_tm_H01.h5", "ectrendecco")

# calculate standard deviations
σT8 = std(V8[:,1:2]*c8[:,1:2]'; dims=2)
σT1 = std(V1[:,1:2]*c1[:,1:2]'; dims=2)
σTargo8 = std(V8[:,1:2]*cargo8[:,1:2]'; dims=2)
σTargo1 = std(V1[:,1:2]*cargo1[:,1:2]'; dims=2)
σTecco8 = std(V8[:,1:2]*cecco8[:,1:2]'; dims=2)
σTecco1 = std(V1[:,1:2]*cecco1[:,1:2]'; dims=2)

fl = Dataset("../shirui/argo/dTz_Nias_H08_argo.nc", "r")
targo8 = DateTime(2004, 1, 15, 0, 0, 0) .+ Month.(0:length(fl["t"][:])-1)
Targo8 = Array{Float64}(replace(fl["__xarray_dataarray_variable__"][:], missing=>NaN))
Targo8 .-= mean(Targo8[:,t8[1].≤targo8.≤t8[end]]; dims=2)

fl = Dataset("../shirui/argo/dTz_Nias_H01_argo.nc", "r")
targo1 = DateTime(2004, 1, 15, 0, 0, 0) .+ Month.(0:length(fl["t"][:])-1)
Targo1 = Array{Float64}(replace(fl["__xarray_dataarray_variable__"][:], missing=>NaN))
Targo1 .-= mean(Targo1[:,t1[1].≤targo1.≤t1[end]]; dims=2)

fl = Dataset("../shirui/ecco/dTz_Nias_H08_ecco.nc", "r")
tecco8 = fl["time"][:] .+ Hour(12)
Tecco8 = Array{Float64}(replace(fl["__xarray_dataarray_variable__"][:], missing=>NaN))
Tecco8 .-= mean(Tecco8[:,t8[1].≤tecco8.≤t8[end]]; dims=2)

fl = Dataset("../shirui/ecco/dTz_Nias_H01_ecco.nc", "r")
tecco1 = fl["time"][:] .+ Hour(12)
Tecco1 = Array{Float64}(replace(fl["__xarray_dataarray_variable__"][:], missing=>NaN))[1:120,:]
Tecco1 .-= mean(Tecco1[:,t1[1].≤tecco1.≤t1[end]]; dims=2)

# small font
rc("font", size=8)
rc("axes", titlesize="medium")

# plot kernels, singular vectors, and resolution matrices
fig, ax = subplots(2, 3, sharex="col", sharey=true, figsize=(6.4, 5.2), gridspec_kw=Dict("width_ratios"=>(1, 1, 1.3)))
ax[1,1].plot(1e3K8', -1e-3z)
ax[2,1].plot(1e3K1', -1e-3z)
ax[1,2].set_prop_cycle(matplotlib.cycler("color", ["tab:olive", "tab:purple"]))
ax[2,2].set_prop_cycle(matplotlib.cycler("color", ["tab:olive", "tab:purple"]))
ax[1,2].plot(V8[:,1:2], -1e-3z)
ax[2,2].plot(V1[:,1:2], -1e-3z)
img8 = ax[1,3].imshow(V8[end:-1:1,1:2]*V8[end:-1:1,1:2]', extent=1e-3*[-z[end]+Δz/2, -z[1]-Δz/2, -z[1]-Δz/2, -z[end]+Δz/2], cmap="RdBu_r", vmin=-8, vmax=8)
img1 = ax[2,3].imshow(V1[end:-1:1,1:2]*V1[end:-1:1,1:2]', extent=1e-3*[-z[end]+Δz/2, -z[1]-Δz/2, -z[1]-Δz/2, -z[end]+Δz/2], cmap="RdBu_r", vmin=-8, vmax=8)
#cb8 = colorbar(img8, ax=ax[1,3], label="resolution matrix", ticks=-0.04:0.02:0.04)
#cb1 = colorbar(img1, ax=ax[2,3], label="resolution matrix", ticks=-0.04:0.02:0.04)
ax[1,1].legend(["2.00 Hz", "3.00 Hz", "4.00 Hz"], frameon=false)
ax[2,1].legend(["2.50 Hz", "3.25 Hz", "4.00 Hz"], frameon=false)
ax[1,2].legend([L"$v_1$", L"$v_2$"], frameon=false, loc=3)
ax[2,2].legend([L"$v_1$", L"$v_2$"], frameon=false, loc=3)
ax[1,1].set_ylim(5, 0)
ax[1,3].set_xlim(0, 5)
#ax[1,2].set_xlim(-0.25, 0.25)
ax[1,1].set_ylabel("depth (km)")
ax[2,1].set_ylabel("depth (km)")
ax[2,1].set_xlabel(L"sensitivity (s$\,$K$^{-1}\,$km$^{-1}$)")
ax[2,2].set_xlabel("singular vector")
ax[2,3].set_xlabel("depth (km)")
ax[1,1].set_title("Diego Garcia")
ax[1,2].set_title("Diego Garcia")
ax[1,3].set_title("Diego Garcia")
ax[2,1].set_title("Cape Leeuwin")
ax[2,2].set_title("Cape Leeuwin")
ax[2,3].set_title("Cape Leeuwin")
for i = 1:2, j=1:3
  ax[i,j].set_title("($(('a':'z')[(i-1)*3+j]))", loc="left")
end
fig.tight_layout(w_pad=2)
fig.savefig("../ms/2021_grl2/fig/kernels.pdf")

fig, ax = subplots(2, 2, sharex=true, sharey=true, figsize=(4.8, 4.8))
ax[1,1].plot(1e3K8', -1e-3z)
ax[1,2].plot(1e3Km8', -1e-3z)
ax[1,1].legend(["2.00 Hz", "3.00 Hz", "4.00 Hz"], frameon=false, loc=3)
ax[2,1].plot(1e3K1', -1e-3z)
ax[2,2].plot(1e3Km1', -1e-3z)
ax[2,1].legend(["2.50 Hz", "3.25 Hz", "4.00 Hz"], frameon=false, loc=3)
ax[1,1].set_ylim(5, 0)
ax[1,1].set_title("Diego Garcia")
ax[1,2].set_title("Diego Garcia")
ax[2,1].set_title("Cape Leeuwin")
ax[2,2].set_title("Cape Leeuwin")
ax[1,1].set_title("SEM", loc="right")
ax[1,2].set_title("modal", loc="right")
ax[2,1].set_title("SEM", loc="right")
ax[2,2].set_title("modal", loc="right")
ax[1,1].set_ylabel("depth (km)")
ax[2,1].set_ylabel("depth (km)")
ax[2,1].set_xlabel(L"sensitivity (s$\,$K$^{-1}\,$km$^{-1}$)")
ax[2,2].set_xlabel(L"sensitivity (s$\,$K$^{-1}\,$km$^{-1}$)")
for i = 1:2, j=1:2
  ax[i,j].set_title("($(('a':'z')[(i-1)*2+j]))", loc="left")
end
fig.tight_layout(w_pad=2)
fig.savefig("../ms/2021_grl2/fig/modekernels.pdf")

fig, ax = subplots(1, 2, sharex=true, sharey=true)
ax[1].plot(σT8, -1e-3z)
ax[1].plot(σTargo8, -1e-3z)
ax[1].plot(σTecco8, -1e-3z)
ax[2].plot(σT1, -1e-3z)
ax[2].plot(σTargo1, -1e-3z)
ax[2].plot(σTecco1, -1e-3z)
ax[1].set_ylim(5, 0)

fig, ax = subplots(1, 2, sharex=true, sharey=true, figsize=(4.8, 4.))
ax[1].plot(1e3*SOT.meanyear*Ttrend8, -1e-3z; zorder=6)
ax[1].fill_betweenx(-1e-3z, 1e3*SOT.meanyear*(Ttrend8 - 2eTtrend8), 1e3*SOT.meanyear*(Ttrend8 + 2eTtrend8); alpha=0.2, zorder=3)
ax[1].plot(1e3*SOT.meanyear*Ttrendargo8, -1e-3z; zorder=5)
ax[1].fill_betweenx(-1e-3z, 1e3*SOT.meanyear*(Ttrendargo8 - 2eTtrendargo8), 1e3*SOT.meanyear*(Ttrendargo8 + 2eTtrendargo8), alpha=0.2, zorder=2)
ax[1].plot(1e3*SOT.meanyear*Ttrendecco8, -1e-3z; zorder=4)
ax[1].fill_betweenx(-1e-3z, 1e3*SOT.meanyear*(Ttrendecco8 - 2eTtrendecco8), 1e3*SOT.meanyear*(Ttrendecco8 + 2eTtrendecco8), alpha=0.2, zorder=1)
ax[2].plot(1e3*SOT.meanyear*Ttrend1, -1e-3z; zorder=6)
ax[2].fill_betweenx(-1e-3z, 1e3*SOT.meanyear*(Ttrend1 - 2eTtrend1), 1e3*SOT.meanyear*(Ttrend1 + 2eTtrend1), alpha=0.2, zorder=3)
ax[2].plot(1e3*SOT.meanyear*Ttrendargo1, -1e-3z; zorder=5)
ax[2].fill_betweenx(-1e-3z, 1e3*SOT.meanyear*(Ttrendargo1 - 2eTtrendargo1), 1e3*SOT.meanyear*(Ttrendargo1 + 2eTtrendargo1), alpha=0.2, zorder=2)
ax[2].plot(1e3*SOT.meanyear*Ttrendecco1, -1e-3z; zorder=4)
ax[2].fill_betweenx(-1e-3z, 1e3*SOT.meanyear*(Ttrendecco1 - 2eTtrendecco1), 1e3*SOT.meanyear*(Ttrendecco1 + 2eTtrendecco1), alpha=0.2, zorder=1)
ax[1].set_ylim(5, 0)
ax[1].set_ylabel("depth (km)")
ax[1].set_xlabel(L"trend (mK$\,$yr$^{-1}$)")
ax[2].set_xlabel(L"trend (mK$\,$yr$^{-1}$)")
ax[1].set_title("Diego Garcia")
ax[2].set_title("Cape Leeuwin")
fig.tight_layout()

# gaps
gaps8 = [Date(2004, 1, 1), Date(2007, 1, 1), Date(2009, 7, 1), Date(2011, 1, 1), Date(2020, 1, 1)]
gaps1 = [Date(2004, 1, 1), Date(2006, 4, 1), Date(2009, 7, 1), Date(2020, 1, 1)]
ir8 = [findfirst(t8 .> gaps8[i]) : findlast(t8 .< gaps8[i+1]) for i = 1:length(gaps8)-1]
ir1 = [findfirst(t1 .> gaps1[i]) : findlast(t1 .< gaps1[i+1]) for i = 1:length(gaps1)-1]

# lengths of time series
l8 = size(c8, 2)
l1 = size(c1, 2)

function ending(i)
  if i == 1
    return "st"
  elseif i == 2
    return "nd"
  elseif i == 3
    return "rd"
  else
    return "th"
  end
end

# plot H08 time series
fig, ax = subplots(l8, 1, sharex=true, figsize=(190/25.4, 190/25.4))
for i = 1:l8-1
  ax[i].scatter(t8, 1e3c8[:,i], s=2, zorder=2, color="tab:blue")
  for j = ir8
    global p1, p2, p3
    p1, = ax[i].plot(t8[j], 1e3c8[j,i], zorder=3, color="tab:blue", linewidth=1)
    ax[i].fill_between(t8[j], 1e3*(c8[j,i]-2ec8[j,i]), 1e3*(c8[j,i]+2ec8[j,i]), alpha=.2, zorder=3, color="tab:blue", linewidth=0)
    p2, = ax[i].plot(t8[j], 1e3cargo8[j,i], zorder=1, linewidth=1, color="tab:orange")
    p3, = ax[i].plot(t8[j], 1e3cecco8[j,i], zorder=2, linewidth=1, color="tab:green")
  end
  ax[i].plot(t8, 1e3ctrend8[:,i], zorder=1, linewidth=1, color="tab:blue")
  ax[i].fill_between(t8, 1e3*(ctrend8[:,i]-2ectrend8[:,i]), 1e3*(ctrend8[:,i]+2ectrend8[:,i]), alpha=.2, zorder=0, color="tab:blue", linewidth=0)
  ax[i].plot(t8, 1e3ctrendargo8[:,i], zorder=0, linewidth=1, color="tab:orange")
  ax[i].fill_between(t8, 1e3*(ctrendargo8[:,i]-2ectrendargo8[:,i]), 1e3*(ctrendargo8[:,i]+2ectrendargo8[:,i]), alpha=.2, zorder=0, color="tab:orange", linewidth=0)
  ax[i].plot(t8, 1e3ctrendecco8[:,i], zorder=0, linewidth=1, color="tab:green")
  ax[i].fill_between(t8, 1e3*(ctrendecco8[:,i]-2ectrendecco8[:,i]), 1e3*(ctrendecco8[:,i]+2ectrendecco8[:,i]), alpha=.2, zorder=0, color="tab:green", linewidth=0)
  ax[i].set_ylabel("temperature anomaly (mK)")
  ax[i].set_title("projection onto the $i$(ending(i)) singular vector")
  ylim = maximum(abs.(ax[i].get_ylim()))
  ax[i].set_ylim(-ylim, ylim)
  ax[i].set_xlim(t8[1]-Day(30), t8[end]+Day(30))
end
ax[1].legend([p1, p2, p3], ["\$T\$ waves", "Argo", "ECCO"], ncol=3, loc="lower center", frameon=false)
for j = ir8
  ax[l8].pcolormesh(t8[j], -1e-3z, V8[:,1:2]*c8[j,1:2]', cmap="RdBu_r", vmin=-.2, vmax=.2, shading="nearest", rasterized="true")
end
ax[l8].set_ylim(5, 0)
ax[l8].set_ylabel("depth (km)")
ax[l8].set_title("temperature reconstruction")
ax[1].set_xticks(ceil(t8[1], Year) : Year(1) : floor(t8[end], Year), minor=true)
for i = 1:l1
  ax[i].set_title("($(('a':'z')[i]))", loc="left")
end
fig.align_ylabels()
fig.tight_layout()
savefig("../ms/2021_grl2/fig/nias_H08_timeseries.pdf", dpi=300)

# plot H01 time series
fig, ax = subplots(l1, 1, sharex=true, figsize=(190/25.4, 190/25.4))
for i = 1:l1-1
  ax[i].scatter(t1, 1e3c1[:,i], s=2, zorder=2, color="tab:blue")
  for j = ir1
    global p1, p2, p3
    p1, = ax[i].plot(t1[j], 1e3c1[j,i], zorder=3, color="tab:blue", linewidth=1)
    ax[i].fill_between(t1[j], 1e3*(c1[j,i]-2ec1[j,i]), 1e3*(c1[j,i]+2ec1[j,i]), alpha=.2, zorder=3, color="tab:blue", linewidth=0)
    p2, = ax[i].plot(t1[j], 1e3cargo1[j,i], zorder=1, linewidth=1, color="tab:orange")
    p3, = ax[i].plot(t1[j], 1e3cecco1[j,i], zorder=2, linewidth=1, color="tab:green")
  end
  ax[i].plot(t1, 1e3ctrend1[:,i], zorder=1, linewidth=1, color="tab:blue")
  ax[i].fill_between(t1, 1e3*(ctrend1[:,i]-2ectrend1[:,i]), 1e3*(ctrend1[:,i]+2ectrend1[:,i]), alpha=.2, zorder=0, color="tab:blue", linewidth=0)
  ax[i].plot(t1, 1e3ctrendargo1[:,i], zorder=0, linewidth=1, color="tab:orange")
  ax[i].fill_between(t1, 1e3*(ctrendargo1[:,i]-2ectrendargo1[:,i]), 1e3*(ctrendargo1[:,i]+2ectrendargo1[:,i]), alpha=.2, zorder=0, color="tab:orange", linewidth=0)
  ax[i].plot(t1, 1e3ctrendecco1[:,i], zorder=0, linewidth=1, color="tab:green")
  ax[i].fill_between(t1, 1e3*(ctrendecco1[:,i]-2ectrendecco1[:,i]), 1e3*(ctrendecco1[:,i]+2ectrendecco1[:,i]), alpha=.2, zorder=0, color="tab:green", linewidth=0)
  ax[i].set_ylabel("temperature anomaly (mK)")
  ax[i].set_title("projection onto the $i$(ending(i)) singular vector")
  ylim = maximum(abs.(ax[i].get_ylim()))
  ax[i].set_ylim(-ylim, ylim)
  ax[i].set_xlim(t1[1]-Day(30), t1[end]+Day(30))
end
ax[1].legend([p1, p2, p3], ["\$T\$ waves", "Argo", "ECCO"], ncol=3, loc="lower center", frameon=false)
for j = ir1
  ax[l1].pcolormesh(t1[j], -1e-3z, V1[:,1:2]*c1[j,1:2]', cmap="RdBu_r", vmin=-.2, vmax=.2, shading="nearest", rasterized="true")
end
ax[l1].set_ylim(5, 0)
ax[l1].set_ylabel("depth (km)")
ax[l1].set_title("temperature reconstruction")
ax[1].set_xticks(ceil(t1[1], Year) : Year(1) : floor(t1[end], Year), minor=true)
for i = 1:l1
  ax[i].set_title("($(('a':'z')[i]))", loc="left")
end
fig.align_ylabels()
fig.tight_layout()
savefig("../ms/2021_grl2/fig/nias_H01_timeseries.pdf", dpi=300)

cnorm(x) = asinh(asinh(asinh(x/0.1)))

m = 5
major = [-1; -0.1; 0; 0.1; 1]
minor = [-5:1; -0.9:0.1:-0.1; -0.09:0.01:-0.01; 0; 0.01:0.01:0.09; 0.1:0.1:0.9; 1:5]

# H08 temperature plots
fig, ax = subplots(3, 2, sharex=true, sharey=true, figsize=(190/25.4, 190/25.4))
for j = ir8
  ax[1,2].pcolormesh(t8[j], -1e-3z, cnorm.(V8[:,1:2]*c8[j,1:2]'), cmap="RdBu_r", vmin=cnorm(-m), vmax=cnorm(m), shading="nearest", rasterized="true")
end
img = ax[2,1].pcolormesh(targo8, -1e-3z, cnorm.(Targo8), cmap="RdBu_r", vmin=cnorm(-m), vmax=cnorm(m), shading="nearest", rasterized="true")
img = ax[3,1].pcolormesh(tecco8, -1e-3z, cnorm.(Tecco8), cmap="RdBu_r", vmin=cnorm(-m), vmax=cnorm(m), shading="nearest", rasterized="true")
img = ax[2,2].pcolormesh(targo8, -1e-3z, cnorm.(V8[:,1:2]*V8[:,1:2]'*replace(Targo8, NaN=>0)*Δz/h), cmap="RdBu_r", vmin=cnorm(-m), vmax=cnorm(m), shading="nearest", rasterized="true")
img = ax[3,2].pcolormesh(tecco8, -1e-3z, cnorm.(V8[:,1:2]*V8[:,1:2]'*replace(Tecco8, NaN=>0)*Δz/h), cmap="RdBu_r", vmin=cnorm(-m), vmax=cnorm(m), shading="nearest", rasterized="true")
ax[1,1].set_visible(false)
cb = colorbar(img, ax=ax[1,1], location="bottom", label="temperature anomaly (K)")
cb.ax.xaxis.set_ticks_position("top")
cb.ax.xaxis.set_label_position("top")
cb.set_ticks(cnorm.(major))
cb.set_ticks(cnorm.(minor); minor=true)
cb.set_ticklabels([@sprintf("\$%g\$", l) for l = major])
for i = 1:5
  ax[i÷2+1,i%2+1].set_title("($(('a':'z')[i]))", loc="left")
end
ax[1,1].set_xlim(Date(2005, 4, 1), Date(2006, 10, 1))
ax[1,1].set_xticks(Date(2005, 7, 1) : Month(6) : Date(2006, 7, 1))
ax[1,1].set_xticks(Date(2005, 4, 1) : Month(1) : Date(2006, 10, 1), minor=true)
ax[1,1].xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y-%m"))
ax[1,1].set_ylim(5, 0)
ax[2,1].set_title("Argo")
ax[3,1].set_title("ECCO")
ax[1,2].set_title(L"$T$ waves")
ax[2,2].set_title("Argo")
ax[3,2].set_title("ECCO")
ax[2,1].set_title("original", loc="right")
ax[3,1].set_title("original", loc="right")
ax[1,2].set_title("projected", loc="right")
ax[2,2].set_title("projected", loc="right")
ax[3,2].set_title("projected", loc="right")
ax[1,1].set_ylabel("depth (km)")
ax[2,1].set_ylabel("depth (km)")
ax[3,1].set_ylabel("depth (km)")
fig.tight_layout(w_pad=2)
savefig("../ms/2021_grl2/fig/nias_H08_depth.pdf", dpi=300)

# H01 temperature plots
fig, ax = subplots(3, 2, sharex=true, sharey=true, figsize=(190/25.4, 190/25.4))
for j = ir1
  ax[1,2].pcolormesh(t1[j], -1e-3z, cnorm.(V1[:,1:2]*c1[j,1:2]'), cmap="RdBu_r", vmin=cnorm(-m), vmax=cnorm(m), shading="nearest", rasterized="true")
end
img = ax[2,1].pcolormesh(targo1, -1e-3z, cnorm.(Targo1), cmap="RdBu_r", shading="nearest", rasterized="true", vmin=cnorm(-m), vmax=cnorm(m))
img = ax[3,1].pcolormesh(tecco1, -1e-3z, cnorm.(Tecco1), cmap="RdBu_r", shading="nearest", rasterized="true", vmin=cnorm(-m), vmax=cnorm(m))
img = ax[2,2].pcolormesh(targo1, -1e-3z, cnorm.(V1[:,1:2]*V1[:,1:2]'*replace(Targo1, NaN=>0)*Δz/h), cmap="RdBu_r", vmin=cnorm(-m), vmax=cnorm(m), shading="nearest", rasterized="true")
img = ax[3,2].pcolormesh(tecco1, -1e-3z, cnorm.(V1[:,1:2]*V1[:,1:2]'*replace(Tecco1, NaN=>0)*Δz/h), cmap="RdBu_r", vmin=cnorm(-m), vmax=cnorm(m), shading="nearest", rasterized="true")
ax[1,1].set_visible(false)
cb = colorbar(img, ax=ax[1,1], location="bottom", label="temperature anomaly (K)")
cb.ax.xaxis.set_ticks_position("top")
cb.ax.xaxis.set_label_position("top")
cb.set_ticks(cnorm.(major))
cb.set_ticks(cnorm.(minor); minor=true)
cb.set_ticklabels([@sprintf("\$%g\$", l) for l = major])
for i = 1:5
  ax[i÷2+1,i%2+1].set_title("($(('a':'z')[i]))", loc="left")
end
ax[1,1].set_xlim(Date(2005, 4, 1), Date(2006, 10, 1))
ax[1,1].set_xticks(Date(2005, 7, 1) : Month(6) : Date(2006, 7, 1))
ax[1,1].set_xticks(Date(2005, 4, 1) : Month(1) : Date(2006, 10, 1), minor=true)
ax[1,1].xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y-%m"))
ax[1,1].set_ylim(5, 0)
ax[2,1].set_title("Argo")
ax[3,1].set_title("ECCO")
ax[1,2].set_title(L"$T$ waves")
ax[2,2].set_title("Argo")
ax[3,2].set_title("ECCO")
ax[2,1].set_title("original", loc="right")
ax[3,1].set_title("original", loc="right")
ax[1,2].set_title("projected", loc="right")
ax[2,2].set_title("projected", loc="right")
ax[3,2].set_title("projected", loc="right")
ax[1,1].set_ylabel("depth (km)")
ax[2,1].set_ylabel("depth (km)")
ax[3,1].set_ylabel("depth (km)")
fig.tight_layout(w_pad=2)
savefig("../ms/2021_grl2/fig/nias_H01_depth.pdf", dpi=300)

# H08 inverted temperature
fig, ax = subplots(3, 1, sharex=true, sharey=true, figsize=(190/25.4, 190/25.4))
for j = ir8
  ax[1].pcolormesh(t8[j], -1e-3z, V8[:,1:2]*c8[j,1:2]', cmap="RdBu_r", vmin=-0.2, vmax=0.2, shading="nearest", rasterized="true")
end
ax[2].pcolormesh(targo8, -1e-3z, V8[:,1:2]*V8[:,1:2]'*replace(Targo8, NaN=>0)*Δz/h, cmap="RdBu_r", vmin=-0.2, vmax=0.2, shading="nearest", rasterized="true")
ax[3].pcolormesh(tecco8, -1e-3z, V8[:,1:2]*V8[:,1:2]'*replace(Tecco8, NaN=>0)*Δz/h, cmap="RdBu_r", vmin=-0.2, vmax=0.2, shading="nearest", rasterized="true")
ax[1].set_xticks(ceil(t8[1], Year) : Year(1) : floor(t8[end], Year), minor=true)
for i = 1:3
  ax[i].set_title("($(('a':'z')[i]))", loc="left")
end
ax[1].set_xlim(t8[1], t8[end])
ax[1].set_ylim(5, 0)
ax[1].set_title(L"$T$ waves")
ax[2].set_title("Argo")
ax[3].set_title("ECCO")
ax[1].set_ylabel("depth (km)")
ax[2].set_ylabel("depth (km)")
ax[3].set_ylabel("depth (km)")
fig.tight_layout()
savefig("../ms/2021_grl2/fig/nias_H08_inv_depth.pdf", dpi=300)

# H01 inverted temperature
fig, ax = subplots(3, 1, sharex=true, sharey=true, figsize=(190/25.4, 190/25.4))
for j = ir1
  ax[1].pcolormesh(t1[j], -1e-3z, V1[:,1:2]*c1[j,1:2]', cmap="RdBu_r", vmin=-0.2, vmax=0.2, shading="nearest", rasterized="true")
end
ax[2].pcolormesh(targo1, -1e-3z, V1[:,1:2]*V1[:,1:2]'*replace(Targo1, NaN=>0)*Δz/h, cmap="RdBu_r", vmin=-0.2, vmax=0.2, shading="nearest", rasterized="true")
ax[3].pcolormesh(tecco1, -1e-3z, V1[:,1:2]*V1[:,1:2]'*replace(Tecco1, NaN=>0)*Δz/h, cmap="RdBu_r", vmin=-0.2, vmax=0.2, shading="nearest", rasterized="true")
ax[1].set_xticks(ceil(t1[1], Year) : Year(1) : floor(t1[end], Year), minor=true)
for i = 1:3
  ax[i].set_title("($(('a':'z')[i]))", loc="left")
end
ax[1].set_xlim(t1[1], t1[end])
ax[1].set_ylim(5, 0)
ax[1].set_title(L"$T$ waves")
ax[2].set_title("Argo")
ax[3].set_title("ECCO")
ax[1].set_ylabel("depth (km)")
ax[2].set_ylabel("depth (km)")
ax[3].set_ylabel("depth (km)")
fig.tight_layout()
savefig("../ms/2021_grl2/fig/nias_H01_inv_depth.pdf", dpi=300)

using HDF5, LinearAlgebra, PyPlot, NCDatasets, Statistics, Dates, Printf

h = 5e3

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

# read trends and their uncertainties
ctrends8 = h5read("results/nias_tm_H08.h5", "ctrends")
ctrends1 = h5read("results/nias_tm_H01.h5", "ctrends")
ctrendsargo8 = h5read("results/nias_tm_H08.h5", "ctrendsargo")
ctrendsargo1 = h5read("results/nias_tm_H01.h5", "ctrendsargo")
ctrendsecco8 = h5read("results/nias_tm_H08.h5", "ctrendsecco")
ctrendsecco1 = h5read("results/nias_tm_H01.h5", "ctrendsecco")
ectrends8 = h5read("results/nias_tm_H08.h5", "ectrends")
ectrends1 = h5read("results/nias_tm_H01.h5", "ectrends")
ectrendsargo8 = h5read("results/nias_tm_H08.h5", "ectrendsargo")
ectrendsargo1 = h5read("results/nias_tm_H01.h5", "ectrendsargo")
ectrendsecco8 = h5read("results/nias_tm_H08.h5", "ectrendsecco")
ectrendsecco1 = h5read("results/nias_tm_H01.h5", "ectrendsecco")

# estimate trend profiles
Ttrend8 = V8[:,1:2]*ctrends8[1:2]
Ttrend1 = V1[:,1:2]*ctrends1[1:2]
Ttrendargo8 = V8[:,1:2]*ctrendsargo8[1:2]
Ttrendargo1 = V1[:,1:2]*ctrendsargo1[1:2]
Ttrendecco8 = V8[:,1:2]*ctrendsecco8[1:2]
Ttrendecco1 = V1[:,1:2]*ctrendsecco1[1:2]

# estimate trend profile uncertainties
eTtrend8 = sqrt.(diag(V8[:,1:2]*Diagonal(ectrends8[1:2].^2)*V8[:,1:2]'))
eTtrend1 = sqrt.(diag(V1[:,1:2]*Diagonal(ectrends1[1:2].^2)*V1[:,1:2]'))
eTtrendargo8 = sqrt.(diag(V8[:,1:2]*Diagonal(ectrendsargo8[1:2].^2)*V8[:,1:2]'))
eTtrendargo1 = sqrt.(diag(V1[:,1:2]*Diagonal(ectrendsargo1[1:2].^2)*V1[:,1:2]'))
eTtrendecco8 = sqrt.(diag(V8[:,1:2]*Diagonal(ectrendsecco8[1:2].^2)*V8[:,1:2]'))
eTtrendecco1 = sqrt.(diag(V1[:,1:2]*Diagonal(ectrendsecco1[1:2].^2)*V1[:,1:2]'))

# load time series
t8 = DateTime(2000, 1, 1) .+ Millisecond.(h5read("results/nias_tm_H08.h5", "t"))
t1 = DateTime(2000, 1, 1) .+ Millisecond.(h5read("results/nias_tm_H01.h5", "t"))
c8 = h5read("results/nias_tm_H08.h5", "c")
c1 = h5read("results/nias_tm_H01.h5", "c")
cargo8 = h5read("results/nias_tm_H08.h5", "cargo")
cargo1 = h5read("results/nias_tm_H01.h5", "cargo")
cecco8 = h5read("results/nias_tm_H08.h5", "cecco")
cecco1 = h5read("results/nias_tm_H01.h5", "cecco")
ec8 = h5read("results/nias_tm_H08.h5", "ec")
ec1 = h5read("results/nias_tm_H01.h5", "ec")
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

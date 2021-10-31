using .SOT, PyPlot, Printf, Dates, LinearAlgebra, Statistics, SparseArrays
using HDF5, Interpolations, DataFrames, CSV

# TODO: resolve clock errors in H08
# TODO: update kernels
# TODO: resolve difference with DGAR

# identifier for experiment
eqname = "nias"

# P-wave (reference) stations
pstations = ["PS.PSI..BHZ", "MY.KUM..BHZ", "II.WRAB.00.BHZ", "GE.GSI..BHZ"]

# intervals to which to cut P waveforms
pintervals = [[-3, 47], [-3, 47], [-3, 47], [-3, 47]]

# frequency bands to which to filter P waveforms
pfreqbands = [[1, 3], [1, 3], [1.5, 2.5], [1, 3]]

# T-wave station
tstations = ["H08S2"]

# T-wave time window around predicted arrival time
tintervals = [[-10, 70]]

# T-wave filtering window width
tavgwidth = 0.5

# T-wave reference frequency at which to find first max CC
treffreq = 2.0

# frequencies used in inversion
tinvfreq = 2.0:1.0:4.0

# minimum CCs for T-wave pairs (at inversion frequencies)
tmincc = 0.6:-0.1:0.4

# excluded time periods: before 2004-12-01 and periods with uncorrected clock error
excludetimes = [[Date(2001, 1, 1) Date(2004, 12, 1)],
                [DateTime("2010-01-23T00:00:00") DateTime("2012-01-20T00:00:00")],
                [Date(2017, 6, 1) Date(2018, 1, 1)]]

# manually exclude pairs
excludepairs = CSV.read("data/catalogs/nias_H08_exclude.csv", DataFrame)

# download P-wave data
SOT.downloadseisdata(eqname, pstations)

# cut and filter P waveforms
SOT.cutpwaves(eqname, pstations, pintervals, pfreqbands)

# find P-wave pairs
SOT.findpairs(eqname, pstations, pintervals, pfreqbands, saveplot=true)

# measure T-wave lags Δτ
SOT.twavepick(eqname, tstations, tintervals, tavgwidth, treffreq, pstations, pintervals,
              pfreqbands, saveplot=true)

# collect usable pairs
tpairs, ppairs = SOT.collectpairs(eqname, tstations, tintervals, tavgwidth, treffreq,
                                  tinvfreq, tmincc, pstations, pintervals, pfreqbands;
                                  excludetimes, excludepairs)

# number of good T- and P-wave pairs
nt = size(tpairs, 1)
np = size(ppairs, 1)

# number of unique events
m = length(unique([tpairs.event1; tpairs.event2]))

# number of frequencies
l = length(tinvfreq)

###

range = h5read("data/kernels/nias_H08.h5", "range")
depth = h5read("data/kernels/nias_H08.h5", "depth")
Δx = range[2] - range[1]
Δz = depth[1] - depth[2]
K = sum(h5read("data/kernels/nias_H08.h5", "K"), dims=2)[:,1,:]'*Δx*Δz

# SVD
U, Σ, V = svd(K)

# correlation time (days)
λ = 30

# solution standard deviation (K)
σx = [0.1, 0.05, 0.05]

# noise (s)
σn = 0.02

# origin time correction standard deviation (s)
σp = 5.0

# trend priors (K/day)
σtrend = 5e-4*σx

# annual cycle prior (K)
σannual = σx

# semi-annual cycle prior (K)
σsemiannual = σx

# get inversion matrices
t, E, Rxx, Rnn, P, D = SOT.invert(tpairs, ppairs, λ, σx, σn, σp, U, Σ; σtrend, σannual, σsemiannual)

tpairs.Δτ = SOT.correctcycleskipping(tpairs, ppairs, E, Rxx, Rnn, P, m)

# collect delays into data vector
y = [reshape(vcat([(tpairs.Δτ[i])' for i = 1:nt]...), l*nt); ppairs.Δτ]

# invert
x = P*E'*inv(Rnn)*y

# extract trends
trends = x[(l+1)*m+1:(l+1)*m+l]
etrends = sqrt.(diag(P[(l+1)*m+1:(l+1)*m+l,(l+1)*m+1:(l+1)*m+l]))

# trend and uncertainty of lowest frequency
@printf("T-wave trend: %+3.1f ± %3.1f mK/yr\n", SOT.meanyear/sum(K[1,:])*1e3*[trends[1] -etrends[1]]...)

ω = 2π/SOT.meanyear
td = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
tm = td[1]+(td[m]-td[1])/2

# reconstruct full travel time anomalies
τ = reshape(D*x, (m, l))
e = reshape(sqrt.(diag(D*P*D')), (m, l))

# reconstruct trends
Dtrend = [zeros(l*m, (l+1)*m) D[:,(l+1)*m+1:(l+1)*m+l] zeros(l*m, 4l)]
τtrend = reshape(Dtrend*x, (m, l))
etrend = reshape(sqrt.(diag(Dtrend*P*Dtrend')), (m, l))

# reconstruct seasonal signal
Dseasonal = [zeros(l*m, (l+1)*m+l) D[:,(l+1)*m+l+1:(l+1)*m+5l]]
τseasonal = reshape(Dseasonal*x, (m, l))
eseasonal = reshape(sqrt.(diag(Dseasonal*P*Dseasonal')), (m, l))

# read and interpolate Argo data
targo = h5read("data/argo/nias_H08.h5", "t")
τargo = hcat([interpolate((targo,), h5read("data/argo/nias_H08.h5", "tau")[:,i], Gridded(Linear()))(td) for i = 1:l]...)

# read and interpolate ECCO data
tecco = h5read("data/ecco/nias_H08.h5", "t")
τecco = hcat([interpolate((tecco,), h5read("data/ecco/nias_H08.h5", "tau")[:,i], Gridded(Linear()))(td) for i = 1:l]...)

# project travel time anomalies onto singular vectors
M = kron(sparse((U/Diagonal(Σ))'), I(m))
τs = reshape(M*D*x, (m, l))
es = reshape(sqrt.(diag(M*D*P*D'*M')), (m, l))
τstrend = reshape(M*Dtrend*x, (m, l))
estrend = reshape(sqrt.(diag(M*Dtrend*P*Dtrend'*M')), (m, l))

# project Argo travel time anomalies onto singular vectors
τsargo = τargo*U/Diagonal(Σ)

# project ECCO travel time anomalies onto singular vectors
τsecco = τecco*U/Diagonal(Σ)

# Argo/ECCO trend setup
E = [kron(I(l), td.-tm) kron(I(l), cos.(ω*td)) kron(I(l), sin.(ω*td)) kron(I(l), cos.(2ω*td)) kron(I(l), sin.(2ω*td))]
Rxx = blockdiag(spdiagm(σtrend.^2), kron(I(2), spdiagm(σannual.^2/2)), kron(I(2), spdiagm(σsemiannual.^2/2)))
A = exp.(-abs.(td.-td')/λ)
Rnn = kron(spdiagm(σx.^2), A)
Ps = inv(inv(Array(Rxx)) + E'*inv(Array(Rnn))*E)
T = kron(I(5), U*Diagonal(Σ))
P = T*Ps*T'
M = zeros(size(E))
M[:,1:3] = E[:,1:3]

# Argo trends
ys = reshape(τsargo, l*m)
xs = Ps*E'*inv(Array(Rnn))*ys
x = T*xs
τargotrend = reshape(M*x, (m, l))
eargotrend = reshape(sqrt.(diag(M*P*M')), (m, l))
τsargotrend = reshape(M*xs, (m, l))
esargotrend = reshape(sqrt.(diag(M*Ps*M')), (m, l))
trendsargo = x[1:3]
etrendsargo = sqrt.(diag(P[1:3,1:3]))

# trend and uncertainty of lowest frequency
@printf("Argo trend: %+3.1f ± %3.1f mK/yr\n", SOT.meanyear/sum(K[1,:])*1e3*[trendsargo[1] -etrendsargo[1]]...)

# ECCO trends
ys = reshape(τsecco, l*m)
xs = Ps*E'*inv(Array(Rnn))*ys
x = T*xs
τeccotrend = reshape(M*x, (m, l))
eeccotrend = reshape(sqrt.(diag(M*P*M')), (m, l))
τseccotrend = reshape(M*xs, (m, l))
eseccotrend = reshape(sqrt.(diag(M*Ps*M')), (m, l))
trendsecco = x[1:3]
etrendsecco = sqrt.(diag(P[1:3,1:3]))

# trend and uncertainty of lowest frequency
@printf("ECCO trend: %+3.1f ± %3.1f mK/yr\n", SOT.meanyear/sum(K[1,:])*1e3*[trendsecco[1] -etrendsecco[1]]...)

# small font
rc("font", size=8)

# plot
fig, ax = subplots(l, 1, sharex=true, figsize=(8, 6.4))
for i = 1:l
  ax[i].plot(t, τ[:,i], color="tab:blue")
  ax[i].fill_between(t, τ[:,i]-2e[:,i], τ[:,i]+2e[:,i], alpha=0.2, linewidth=0, color="tab:blue")
  ax[i].plot(t, τtrend[:,i], color="tab:blue")
  ax[i].fill_between(t, τtrend[:,i]-2etrend[:,i], τtrend[:,i]+2etrend[:,i], alpha=0.2, linewidth=0, color="tab:blue")
  ax[i].plot(t, τargo[:,i], zorder=0, color="tab:orange")
  ax[i].plot(t, τargotrend[:,i], zorder=0, color="tab:orange")
  ax[i].fill_between(t, τargotrend[:,i]-2eargotrend[:,i], τargotrend[:,i]+2eargotrend[:,i], alpha=0.2, zorder=0, linewidth=0, color="tab:orange")
  ax[i].plot(t, τecco[:,i], zorder=0, color="tab:green")
  ax[i].plot(t, τeccotrend[:,i], zorder=0, color="tab:green")
  ax[i].fill_between(t, τeccotrend[:,i]-2eeccotrend[:,i], τeccotrend[:,i]+2eeccotrend[:,i], alpha=0.2, zorder=0, linewidth=0, color="tab:green")
end
fig.tight_layout()

ir = [1:796, 797:997, 998:m]
# plot
fig, ax = subplots(l, 1, sharex=true, figsize=(190/25.4, 190/25.4))
for i = 1:l-1
  ax[i].scatter(t, τs[:,i], s=2, zorder=2, color="tab:blue")
  for j = ir
    global p1, p2, p3
    p1, = ax[i].plot(t[j], τs[j,i], zorder=3, color="tab:blue", linewidth=1)
    ax[i].fill_between(t[j], τs[j,i]-2es[j,i], τs[j,i]+2es[j,i], alpha=.2, zorder=3, color="tab:blue", linewidth=0)
    p2, = ax[i].plot(t[j], τsargo[j,i], zorder=1, linewidth=1, color="tab:orange")
    p3, = ax[i].plot(t[j], τsecco[j,i], zorder=2, linewidth=1, color="tab:green")
  end
  ax[i].plot(t, τstrend[:,i], zorder=1, color="tab:blue")
  ax[i].fill_between(t, τstrend[:,i]-2estrend[:,i], τstrend[:,i]+2estrend[:,i], alpha=.2, zorder=0, color="tab:blue", linewidth=0)
  ax[i].plot(t, τsargotrend[:,i], zorder=0, color="tab:orange")
  ax[i].fill_between(t, τsargotrend[:,i]-2esargotrend[:,i], τsargotrend[:,i]+2esargotrend[:,i], alpha=.2, zorder=0, color="tab:orange", linewidth=0)
  ax[i].plot(t, τseccotrend[:,i], zorder=0, color="tab:green")
  ax[i].fill_between(t, τseccotrend[:,i]-2eseccotrend[:,i], τseccotrend[:,i]+2eseccotrend[:,i], alpha=.2, zorder=0, color="tab:green", linewidth=0)
  ax[i].set_ylabel("\$T_$i\$ (K)")
  ylim = maximum(abs.(ax[i].get_ylim()))
  ax[i].set_ylim(-ylim, ylim)
  ax[i].set_xlim(t[1]-Day(30), t[m]+Day(30))
end
ax[1].set_yticks(-.6:.3:.6)
ax[1].legend([p1, p2, p3], ["\$T\$ waves", "Argo", "ECCO"], ncol=3, loc="lower center", frameon=false)
for j = ir
  ax[l].pcolormesh(t[j], 1e-3depth, V[:,1:2]*τs[j,1:2]', cmap="RdBu_r", vmin=-.15, vmax=.15, shading="nearest", rasterized="true")
end
ax[l].set_ylim(5, 0)
ax[l].set_ylabel("depth (km)")
fig.align_ylabels()
fig.tight_layout()

fig, ax = subplots(1, 3, sharey=true)
ax[1].plot(1e3K'/Δz, 1e-3depth)
ax[1].set_ylim(5, 0)
ax[1].set_xlabel(L"kernel (s$\,$K$^{-1}\,$km$^{-1}$)")
ax[1].set_ylabel("depth (km)")
ax[1].legend(["2 Hz", "3 Hz", "4 Hz"], frameon=false)
ax[2].plot(V, 1e-3depth)
ax[2].legend([L"$v_1$", L"$v_2$", L"$v_3$"], frameon=false)
ax[3].plot(1e3*SOT.meanyear*V[:,1:2]*inv(Diagonal(Σ[1:2]))*U[:,1:2]'*trends, 1e-3depth)
ax[3].plot(1e3*SOT.meanyear*V[:,1:2]*inv(Diagonal(Σ[1:2]))*U[:,1:2]'*trendsargo, 1e-3depth)
ax[3].plot(1e3*SOT.meanyear*V[:,1:2]*inv(Diagonal(Σ[1:2]))*U[:,1:2]'*trendsecco, 1e-3depth)
ax[3].set_xlabel(L"trend (mK$\,$yr$^{-1}$)")
ax[3].legend(["\$T\$ waves", "Argo", "ECCO"], frameon=false)
fig.tight_layout()

using .SOT, PyPlot, Printf, Dates, LinearAlgebra, Statistics, SparseArrays
using HDF5, Interpolations, DataFrames, CSV

# identifier for experiment
eqname = "nias_tm"

# P-wave (reference) stations
pstations = ["PS.PSI..BHZ", "GE.GSI..BHZ"]

# intervals to which to cut P waveforms
pintervals = [[-3, 47], [-3, 47]]

# frequency bands to which to filter P waveforms
pfreqbands = [[1, 3], [1, 3]]

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
#excludetimes = [[Date(2001, 1, 1) Date(2004, 12, 1)],
#                [DateTime("2010-03-16T00:00:00") DateTime("2010-05-17T02:06:17.760")],
#                [Date(2017, 6, 1) Date(2018, 1, 1)]]

# manually exclude pairs
excludepairs = CSV.read("data/catalogs/nias_tm_H08_exclude.csv", DataFrame)

# measure T-wave lags Δτ
SOT.twavepick(eqname, tstations, tintervals, tavgwidth, treffreq, pstations, pintervals,
              pfreqbands, saveplot=true)

# collect usable pairs
tpairs, ppairs = SOT.collectpairs(eqname, tstations, tintervals, tavgwidth, treffreq,
                                  tinvfreq, tmincc, pstations, pintervals, pfreqbands;
                                  excludetimes, excludepairs)

## H08 clock error correction
#function timingcorrection!(tpairs, starttime, endtime, c)
#  idx1 = starttime .< tpairs.event1 .< endtime
#  idx2 = starttime .< tpairs.event2 .< endtime
#  for i = findall(idx1)
#    tpairs.Δτl[i] .-= c*((tpairs.event1[i] .- starttime)./(endtime - starttime))
#    tpairs.Δτc[i] .-= c*((tpairs.event1[i] .- starttime)./(endtime - starttime))
#    tpairs.Δτr[i] .-= c*((tpairs.event1[i] .- starttime)./(endtime - starttime))
#  end
#  for i = findall(idx2)
#    tpairs.Δτl[i] .+= c*((tpairs.event2[i] .- starttime)./(endtime - starttime))
#    tpairs.Δτc[i] .+= c*((tpairs.event2[i] .- starttime)./(endtime - starttime))
#    tpairs.Δτr[i] .+= c*((tpairs.event2[i] .- starttime)./(endtime - starttime))
#  end
#end
#
## correct based on comparison with DGAR
#timingcorrection!(tpairs, DateTime("2010-01-23T08:13:52.666"),
#                  DateTime("2010-03-16T00:00:00.000"), 2.276)
#timingcorrection!(tpairs, DateTime("2010-05-17T02:06:17.760"),
#                  DateTime("2010-10-09T23:59:52.747"), 7.045)
#timingcorrection!(tpairs, DateTime("2010-10-09T23:59:52.747"),
#                  DateTime("2010-12-16T02:03:23.040"), 3.454)
#timingcorrection!(tpairs, DateTime("2010-12-16T02:03:23.040"),
#                  DateTime("2011-02-22T09:43:45.858"), 3.422)
#timingcorrection!(tpairs, DateTime("2011-02-22T09:43:45.858"),
#                  DateTime("2011-03-24T07:55:37.542"), 1.552)
#timingcorrection!(tpairs, DateTime("2011-03-24T07:55:37.542"),
#                  DateTime("2011-04-18T01:59:00.042"), 1.278)
#timingcorrection!(tpairs, DateTime("2011-04-18T01:59:00.042"),
#                  DateTime("2011-08-27T23:59:52.317"), 6.714)
#timingcorrection!(tpairs, DateTime("2011-08-27T23:59:52.317"),
#                  DateTime("2011-12-24T23:59:52.453"), 6.683)
#timingcorrection!(tpairs, DateTime("2011-12-24T23:59:52.453"),
#                  DateTime("2012-01-20T00:00:00.000"), 1.374)

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

###

CSV.write("results/nias_tm_H08_tpairs.csv", tpairs)
CSV.write("results/nias_tm_H08_ppairs.csv", ppairs)

## perform inversion
#t, E, S, P, D = SOT.invert(tpairs, ppairs)
#
## number of good T- and P-wave pairs
#nt = size(tpairs, 1)
#np = size(ppairs, 1)
#
## number of unique events
#m = length(t)
#
## number of frequencies
#l = length(tinvfreq)
#
## invert for P-wave delays without smoothing
#Δτp = E[1:nt,1:m]*(E[l*nt+1:l*nt+np,l*m+1:(l+1)*m]\ppairs.Δτ)
#
## make cycle-skipping correction
#tpairs.Δτ = SOT.correctcycleskipping(tpairs, ppairs, E, S, P)
#
## data vector
#y = [reshape(vcat(tpairs.Δτ'...), l*nt); ppairs.Δτ]
#
## get least-squares solution
#x = P*(E'*y)
#
## variance estimate
#n = y - E*x
#σt = sqrt(mean(n[1:l*nt].^2)*(l*nt+np)/(l*nt+np-l*m-1))
#σp = sqrt(mean(n[l*nt+1:l*nt+np].^2)*(l*nt+np)/(l*nt+np-l*m-1))
#R = spdiagm(0 => [σt^2*ones(l*nt); σp^2*ones(np)])
#
## matrix for error propagation
#A = P*E'
#
## get travel time anomalies and their errors
#τ = reshape(D*x, (m, l))
#τerr = reshape(sqrt.(diag(D*A*R*A'*D')), (m, l))
#
## get frequency differences and their errors
#δD = [vcat([I(m) for i = 1:l-1]...) -I((l-1)*m) spzeros((l-1)*m, m)]
#δτ = reshape(δD*x, (m, l-1))
#δτerr = reshape(sqrt.(diag(δD*A*R*A'*δD')), (m, l-1))
#
## read and interpolate Argo data
#targo = Dates.value.(Date(2004, 1, 15) .+ Month.(h5read("data/argo/nias_H08.h5", "time"))
#                     .- Date(2000, 1, 1))
#tr = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
#τargo = hcat([interpolate((targo,), h5read("data/argo/nias_H08.h5", "tau")[i,:],
#                          Gridded(Linear()))(tr) for i = 1:l]...)
#nz = size(h5read("data/argo/nias_H08.h5", "T"), 1)
#Targo = hcat([interpolate((targo,), h5read("data/argo/nias_H08.h5", "T")[i,:],
#                          Gridded(Linear()))(tr) for i = 1:nz]...)
#
## read and interpolate ECCO data
#tecco = h5read("data/ecco/nias_H08.h5", "time")
#tr = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
#τecco = hcat([interpolate((tecco,), h5read("data/ecco/nias_H08.h5", "tau")[i,:],
#                          Gridded(Linear()))(tr) for i = 1:l]...)
#nz = size(h5read("data/ecco/nias_H08.h5", "T"), 1)
#Tecco = hcat([interpolate((tecco,), h5read("data/ecco/nias_H08.h5", "T")[i,:],
#                          Gridded(Linear()))(tr) for i = 1:nz]...)
#
## save timeseries to file
#tr = Dates.value.(t - DateTime(2000, 1, 1, 0, 0, 0))/1000/3600/24
#h5open("results/nias_tm_H08.h5", "w") do file
#  write(file, "t", tr)
#  write(file, "tau", τ)
#  write(file, "tauerr", τerr)
#  write(file, "dtau", δτ)
#  write(file, "dtauerr", δτerr)
#  write(file, "tauargo", τargo)
#  write(file, "tauecco", τecco)
#  write(file, "Targo", Targo)
#  write(file, "Tecco", Tecco)
#end
#
## offset based on trends
#tr = Dates.value.(t .- t[1])/1000/3600/24
#for i = 1:l
#  c, _ = SOT.lineartrend(t, τ[:,i]; fitannual=true, fitsemiannual=true)
#  τ[:,i] .-= c[2] + c[1]*tr[end]/2
#  if i == 1
#    δτ[:,:] .-= c[2] + c[1]*tr[end]/2
#  else
#    δτ[:,i-1] .+= c[2] + c[1]*tr[end]/2
#  end
#  c, _ = SOT.lineartrend(t, τargo[:,i]; fitannual=true, fitsemiannual=true)
#  τargo[:,i] .-= c[2] + c[1]*tr[end]/2
#  c, _ = SOT.lineartrend(t, τecco[:,i]; fitannual=true, fitsemiannual=true)
#  τecco[:,i] .-= c[2] + c[1]*tr[end]/2
#end
#
## estimate trends
#cτ, _ = SOT.lineartrend(t, τ[:,1]; fitannual=true, fitsemiannual=true)
#cτargo, _ = SOT.lineartrend(t, τargo[:,1]; fitannual=true, fitsemiannual=true)
#cτecco, _ = SOT.lineartrend(t, τecco[:,1]; fitannual=true, fitsemiannual=true)
#cδτ, _ = SOT.lineartrend(t, δτ[:,l-1]; fitannual=true, fitsemiannual=true)
#cδτargo, _ = SOT.lineartrend(t, τargo[:,1] - τargo[:,l]; fitannual=true, fitsemiannual=true)
#cδτecco, _ = SOT.lineartrend(t, τecco[:,1] - τecco[:,l]; fitannual=true, fitsemiannual=true)
#
## plot measured vs. inverted T-wave delays (lowest freq.)
#fig, ax = subplots(1, 1)
#ax.scatter(y[1:nt,1], E[1:nt,:]*x[:,1], s=5)
#ax.set_aspect(1)
#xl = ax.get_xlim()
#yl = ax.get_ylim()
#a = [-20, 20]
#ax.plot(a, a, color="black", linewidth=.8)
#ax.plot(a, a .+ 1/tinvfreq[1], color="black", linewidth=.8, zorder=0)
#ax.plot(a, a .- 1/tinvfreq[1], color="black", linewidth=.8, zorder=0)
#ax.plot(a, a .+ 1/2tinvfreq[1], color="black", linewidth=.8, zorder=0)
#ax.plot(a, a .- 1/2tinvfreq[1], color="black", linewidth=.8, zorder=0)
#ax.set_xlim(xl)
#ax.set_ylim(yl)
#ax.set_title("T waves")
#ax.set_xlabel("measured delay (s)")
#ax.set_ylabel("inverted delay (s)")
#
## plot measured vs. inverted P-wave delays (lowest freq.)
#fig, ax = subplots(1, 1)
#ax.scatter(y[nt+1:nt+np], E[nt+1:nt+np,:]*x[:,1], s=5)
#ax.set_aspect(1)
#xl = ax.get_xlim()
#yl = ax.get_ylim()
#a = [-20, 20]
#ax.plot(a, a, color="black", linewidth=.8)
#ax.plot(a, a .+ 1/tinvfreq[1], color="black", linewidth=.8, zorder=0)
#ax.plot(a, a .- 1/tinvfreq[1], color="black", linewidth=.8, zorder=0)
#ax.plot(a, a .+ 1/2tinvfreq[1], color="black", linewidth=.8, zorder=0)
#ax.plot(a, a .- 1/2tinvfreq[1], color="black", linewidth=.8, zorder=0)
#ax.set_xlim(xl)
#ax.set_ylim(yl)
#ax.set_title("P waves")
#ax.set_xlabel("measured delay (s)")
#ax.set_ylabel("inverted delay (s)")
#
## insert gaps
#function insertgap(t, ys, tgap)
#  idx = findlast(t .< tgap)
#  m = length(t)
#  if !isnothing(idx) && idx < m
#    tg = [t[1:idx]; tgap; t[idx+1:m]]
#    ysg = []
#    for y in ys
#      l = size(y, 2)
#      push!(ysg, [y[1:idx,:]; [NaN for i=1:l]'; y[idx+1:m,:]])
#    end
#  else
#    tg = t
#    ysg = ys
#  end
#  return tg, ysg...
#end
#
## hydrophone gap
#tg, τg, τerrg, δτg, δτerrg, τargog, τeccog = insertgap(t, [τ, τerr, δτ, δτerr, τargo,
#                                                           τecco], Date(2007, 1, 1))
#
## excluded periods
#for i = 1:length(excludetimes)
#  global tg, τg, τerrg, δτg, δτerrg, τargog, τeccog
#  tg, τg, τerrg, δτg, δτerrg, τargog, τeccog = insertgap(tg, [τg, τerrg, δτg, δτerrg,
#                                                              τargog, τeccog],
#                                                         excludetimes[i][1])
#end
#
## seismic station gap
#tg, τg, τerrg, δτg, δτerrg, τargog, τeccog = insertgap(tg, [τg, τerrg, δτg, δτerrg, τargog,
#                                                            τeccog], Date(2009, 7, 1))
#
## plot timeseries
#fig, ax = subplots(2, 1, figsize=(16, 6.4), sharex=true)
#ax[1].plot(tg, τg[:,1], color="tab:blue", zorder=3, label=L"$T$ waves")
#ax[1].scatter(t, τ[:,1], s=2, c="tab:blue", zorder=3)
#ax[1].fill_between(tg, τg[:,1] - 2τerrg[:,1], τg[:,1] + 2τerrg[:,1], alpha=.25,
#                   color="tab:blue", linewidths=0, zorder=3)
#ax[2].plot(tg, δτg[:,l-1], color="tab:blue", zorder=3, label=L"$T$ waves")
#ax[2].scatter(t, δτ[:,l-1], s=2, color="tab:blue", zorder=3)
#ax[2].fill_between(tg, δτg[:,l-1] - 2δτerrg[:,l-1], δτg[:,l-1] + 2δτerrg[:,l-1], alpha=.25,
#                   color="tab:blue", linewidths=0, zorder=3)
#ax[1].plot(tg, τargog[:,1], color="tab:orange", zorder=1, label="Argo")
#ax[2].plot(tg, τargog[:,1] - τargog[:,l], color="tab:orange", zorder=1, label="Argo")
#ax[1].plot(tg, τeccog[:,1], color="tab:green", zorder=2, label="ECCO")
#ax[2].plot(tg, τeccog[:,1] - τeccog[:,l], color="tab:green", zorder=2, label="ECCO")
#ax[1].plot(t, cτargo[1]*tr .+ cτargo[2], color="tab:orange", zorder=1)
#ax[1].plot(t, cτecco[1]*tr .+ cτecco[2], color="tab:green", zorder=2)
#ax[1].plot(t, cτ[1]*tr .+ cτ[2], color="tab:blue", zorder=3)
#ax[2].plot(t, cδτargo[1]*tr .+ cδτargo[2], color="tab:orange", zorder=1)
#ax[2].plot(t, cδτecco[1]*tr .+ cδτecco[2], color="tab:green", zorder=2)
#ax[2].plot(t, cδτ[1]*tr .+ cδτ[2], color="tab:blue", zorder=3)
#ax[1].invert_yaxis()
#ax[1].legend(frameon=false, loc=4, ncol=3)
#ax[2].legend(frameon=false, loc=4, ncol=3)
#ax[1].set_xlim(t[1] - (t[end] - t[1])÷100, t[end] + (t[end] - t[1])÷100)
#ax[1].set_ylabel("travel time anomaly (s)")
#ax[2].set_ylabel("travel time difference (s)")
#fig.align_ylabels()
#fig.tight_layout()
#
## plot measured vs. inverted T-wave delays with P-wave delays estimated from inversion
## without smoothing
#fig, ax = subplots()
#ax.set_aspect(1)
#ax.scatter(y[1:nt,1] - Δτp[:,1], E[1:nt,:]*x[:,1] - Δτp[:,1], s=5)
#xl = ax.get_xlim()
#yl = ax.get_ylim()
#a = [-20, 20]
#ax.plot(a, a, color="black", linewidth=.8)
#ax.plot(a, a .+ 1/tinvfreq[1], color="black", linewidth=.8, zorder=0)
#ax.plot(a, a .- 1/tinvfreq[1], color="black", linewidth=.8, zorder=0)
#ax.set_xlim(xl)
#ax.set_ylim(yl)
#ax.set_xlabel("measured delay (s)")
#ax.set_ylabel("inverted delay (s)")

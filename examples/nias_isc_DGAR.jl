using SOT, PyPlot, Printf, Dates, LinearAlgebra, Statistics, SparseArrays, HDF5, Interpolations, DataFrames, CSV

# identifier for experiment
eqname = "nias_isc"

# P-wave (reference) stations
pstations = ["PS.PSI..BHZ", "MY.KUM..BHZ", "II.WRAB.00.BHZ", "GE.GSI..BHZ"]

# P-wave download source
psrc = ["IRIS", "IRIS", "IRIS", "GFZ"]

# intervals to which to cut P waveforms
pintervals = [[-3, 47], [-3, 47], [-3, 47], [-3, 47]]

# frequency bands to which to filter P waveforms
pfreqbands = [[1, 3], [1, 3], [1.5, 2.5], [1, 3]]

# T-wave station
tstations = ["II.DGAR.00.BHZ"]

# T-wave download source
tsrc = ["IRIS"]

# T-wave time window around predicted arrival time
tintervals = [[-10, 70]]

# T-wave filtering window width
tavgwidth = 0.5

# T-wave reference frequency at which to find first max CC
treffreq = 2.0

# frequencies used in inversion
tinvfreq = [2.0, 3.0]

# minimum CCs for T-wave pairs (at inversion frequencies)
tmincc = [0.6, 0.5]

# excluded time period: before 2004-12-01
# TODO: include more recent data?
excludetimes = [[Date(2001, 1, 1) Date(2004, 12, 1)], [Date(2018, 1, 1) Date(2022, 1, 1)]]

# manually exclude pairs
excludepairs = CSV.read("data/catalogs/nias_isc_DGAR_exclude.csv", DataFrame)

# reference depth (for normalization of singular vectors)
h = 5e3

# download P-wave data
downloadseisdata(eqname, pstations; src=psrc)

# cut and filter P waveforms
cutpwaves(eqname, pstations, pintervals, pfreqbands)

# find P-wave pairs
findpairs(eqname, pstations, pintervals, pfreqbands)

# download T-wave data
downloadseisdata(eqname, tstations; src=tsrc)

# measure T-wave lags Δτ
twavepick(eqname, tstations, tintervals, tavgwidth, treffreq, pstations, pintervals, pfreqbands)

# collect usable pairs
tpairs, ppairs = collectpairs(eqname, tstations, tintervals, tavgwidth, treffreq, tinvfreq, tmincc, pstations, pintervals, pfreqbands; excludetimes, excludepairs)

# DGAR clock error correction
breaktime = DateTime(2012, 3, 17)
idx = (tpairs.event1 .< breaktime) .& (tpairs.event2 .> breaktime)
for i = findall(idx)
  tpairs.Δτl[i] .-= 1.0
  tpairs.Δτc[i] .-= 1.0
  tpairs.Δτr[i] .-= 1.0
end

# number of good T- and P-wave pairs
nt = size(tpairs, 1)
np = size(ppairs, 1)

# number of unique events
m = length(unique([tpairs.event1; tpairs.event2]))

# number of frequencies
l = length(tinvfreq)

###

range = h5read("data/temperature/nias_DGAR.h5", "x")
depth = -h5read("data/temperature/nias_DGAR.h5", "z")
Δx = range[2] - range[1]
Δz = depth[1] - depth[2]
K = sum(h5read("data/temperature/nias_DGAR.h5", "K"), dims=2)[:,1,1:l]'*Δx

# SVD
U, Λ, V = svd(K)

# rescale to be independent of resolution (∫v'vdz = h)
V /= sqrt(Δz/h)
Λ *= sqrt(Δz/h)

# correlation time (days)
λ = 15

# solution standard deviation for coefficients of singular vectors (K)
σc = 1e-3*[15, 10]

# noise (s)
σn = 0.01

# origin time correction standard deviation (s)
σp = 5.0

# trend priors for coefficients of singular vectors (K/day)
σtrend = 1e-3*[4, 2]/SOT.meanyear

# annual cycle prior (K)
σannual = 1e-3*[15, 10]

# semi-annual cycle prior (K)
σsemiannual = 1e-3*[15, 10]

# get inversion matrices
t, E, R, N, P, D = invert(tpairs, ppairs, λ, σc, σn, σp, U, Λ, Δz, h; σtrend, σannual, σsemiannual)

tpairs.Δτ = correctcycleskipping(tpairs, ppairs, E, R, N, P, m)

# collect delays into data vector
y = [reshape(vcat([(tpairs.Δτ[i])' for i = 1:nt]...), l*nt); ppairs.Δτ]

# invert
a = P*E'*inv(N)*y

# plot residuals with and without smoothing
n1 = y - E*a
n2 = y - E*(E\y)
b = -0.16:0.005:0.16
@assert all((b[1] .< n1 .< b[end]) .& (b[1] .< n2 .< b[end]))
fig, ax = subplots(2, l+1; sharex=true, sharey=true)
for i = 1:l
  ax[1,i].hist(n1[(i-1)*nt+1:i*nt]; bins=b)
  ax[2,i].hist(n2[(i-1)*nt+1:i*nt]; bins=b)
  ax[1,i].set_title("T wave, freq. $i")
  ax[2,i].set_title("T wave, freq. $i")
end
ax[1,l+1].hist(n1[l*nt+1:l*nt+np]; bins=b)
ax[2,l+1].hist(n2[l*nt+1:l*nt+np]; bins=b)
ax[1,l+1].set_title("P wave")
ax[2,l+1].set_title("P wave")
ax[1,1].set_yscale("log")
ax[1,1].set_ylabel("count")
ax[2,1].set_ylabel("count")
for i = 1:l+1
  ax[2,l].set_xlabel("residual (s)")
end
fig.tight_layout()

# extract trends
trends = a[(l+1)*m+1:(l+1)*m+l]
Ptrends = P[(l+1)*m+1:(l+1)*m+l,(l+1)*m+1:(l+1)*m+l]

# extract annual amplitude
annuals = hypot.(a[(l+1)*m+l+1:(l+1)*m+2l], a[(l+1)*m+2l+1:(l+1)*m+3l])

# extract semi-annual amplitude
semiannuals = hypot.(a[(l+1)*m+3l+1:(l+1)*m+4l], a[(l+1)*m+4l+1:(l+1)*m+5l])

ω = 2π/SOT.meanyear
td = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24

# reconstruct full travel time anomalies
τ = reshape(D*a, (m, l))
eτ = reshape(sqrt.(diag(D*P*D')), (m, l))

# reconstruct random signal
Drand = [D[:,1:(l+1)*m] zeros(l*m, 5l)]
τrand = reshape(Drand*a, (m, l))

# reconstruct trends
Dtrend = [zeros(l*m, (l+1)*m) D[:,(l+1)*m+1:(l+1)*m+l] zeros(l*m, 4l)]
τtrend = reshape(Dtrend*a, (m, l))
eτtrend = reshape(sqrt.(diag(Dtrend*P*Dtrend')), (m, l))

# reconstruct seasonal signal
Dseasonal = [zeros(l*m, (l+1)*m+l) D[:,(l+1)*m+l+1:(l+1)*m+5l]]
τseasonal = reshape(Dseasonal*a, (m, l))
eτseasonal = reshape(sqrt.(diag(Dseasonal*P*Dseasonal')), (m, l))

# project travel time anomalies and trend reconstructions onto singular vectors
T = h^-1*kron(sparse(inv(Diagonal(Λ))*U'), I(m))
c = reshape(T*D*a, (m, l))
ec = reshape(sqrt.(diag(T*D*P*D'*T')), (m, l))
ctrend = reshape(T*Dtrend*a, (m, l))
ectrend = reshape(sqrt.(diag(T*Dtrend*P*Dtrend'*T')), (m, l))
crand = reshape(T*Drand*a, (m, l))

# project trends onto singular vectors
T = h^-1*inv(Diagonal(Λ))*U'
ctrends = T*trends
ectrends = sqrt.(diag(T*Ptrends*T'))
tannuals = kron(I(2), T)*a[(l+1)*m+l+1:(l+1)*m+3l]
cannuals = hypot.(tannuals[1:l], tannuals[l+1:2l])
J = spzeros(l, 2l)
for i = 1:l
  J[i,2(i-1)+1:2i] = abs.(tannuals[i:l:i+l])/cannuals[i]
end
ecannuals = sqrt.(diag(J*kron(I(2), T)*P[(l+1)*m+l+1:(l+1)*m+3l,(l+1)*m+l+1:(l+1)*m+3l]*kron(I(2), T)'*J'))
tsemiannuals = kron(I(2), T)*a[(l+1)*m+3l+1:(l+1)*m+5l]
csemiannuals = hypot.(tsemiannuals[1:l], tsemiannuals[l+1:2l])
J = spzeros(l, 2l)
for i = 1:l
  J[i,2(i-1)+1:2i] = abs.(tsemiannuals[i:l:i+l])/csemiannuals[i]
end
ecsemiannuals = sqrt.(diag(J*kron(I(2), T)*P[(l+1)*m+l+1:(l+1)*m+3l,(l+1)*m+l+1:(l+1)*m+3l]*kron(I(2), T)'*J'))

# interpolate onto regular grid
ti = DateTime(2005, 1, 1, 12, 0, 0) : Day(1) : DateTime(2005, 12, 31, 12, 0, 0)
τi, ci = regulargrid(td, ti, a, R, λ, h, U, Λ, σc)

# read and interpolate Argo data
targo = h5read("data/temperature/nias_DGAR.h5", "targo")
τargo = h5read("data/temperature/nias_DGAR.h5", "τargo")
Targo = h5read("data/temperature/nias_DGAR.h5", "Targo")
τargo = hcat([interpolate((targo,), τargo[:,i], Gridded(Linear()))(td) for i = 1:l]...)
Targo = hcat([interpolate((targo,), Targo[i,:], Gridded(Linear()))(td) for i = 1:length(depth)]...)

# read and interpolate ECCO data
tecco = h5read("data/temperature/nias_DGAR.h5", "tecco")
τecco = h5read("data/temperature/nias_DGAR.h5", "τecco")
#Tecco = h5read("data/temperature/nias_DGAR.h5", "Tecco")
τecco = hcat([interpolate((tecco,), τecco[:,i], Gridded(Linear()))(td) for i = 1:l]...)
#Tecco = hcat([interpolate((tecco,), Tecco[i,:], Gridded(Linear()))(td) for i = 1:length(depth)]...)

# project Argo and ECCO travel time anomalies onto singular vectors
cargo = τargo*T'
cecco = τecco*T'

# estimate ECCO and Argo trends
Et, Rxxt, Rnnt, Pt, Mt = estimatetrend(td, λ, σc, σc, σtrend, σannual, σsemiannual)
Tt = h*kron(I(6), U*Diagonal(Λ))

# Argo trends
yt = reshape(cargo, l*m)
xt = Pt*Et'*(Rnnt\yt)
# travel time anomalies associated with the trend, plus standard errors
τargotrend = reshape(Mt*Tt*xt, (m, l))
eτargotrend = reshape(sqrt.(diag(Mt*Tt*Pt*Tt'*Mt')), (m, l))
# travel time anomalies associated with the trend, projected onto the sinular vectors, plus standard errors
ctrendargo = reshape(Mt*xt, (m, l))
ectrendargo = reshape(sqrt.(diag(Mt*Pt*Mt')), (m, l))
# trends for travel time anomalies
trendsargo = Tt[l+1:2l,:]*xt
# trends for travel time anomalies projected onto singular vectors, plus uncertainties
ctrendsargo = xt[l+1:2l]
ectrendsargo = sqrt.(diag(Pt[l+1:2l,l+1:2l]))
# annuals for travel time anomalies projected onto singular vectors, plus uncertainties
tannualsargo = xt[2l+1:4l]
cannualsargo = hypot.(tannualsargo[1:l], tannualsargo[l+1:2l])
J = spzeros(l, 2l)
for i = 1:l
  J[i,2(i-1)+1:2i] = abs.(tannualsargo[i:l:i+l])/cannualsargo[i]
end
ecannualsargo = sqrt.(diag(J*Pt[2l+1:4l,2l+1:4l]*J'))
# semiannuals for travel time anomalies projected onto singular vectors, plus uncertainties
tsemiannualsargo = xt[4l+1:6l]
csemiannualsargo = hypot.(tsemiannualsargo[1:l], tsemiannualsargo[l+1:2l])
J = spzeros(l, 2l)
for i = 1:l
  J[i,2(i-1)+1:2i] = abs.(tsemiannualsargo[i:l:i+l])/csemiannualsargo[i]
end
ecsemiannualsargo = sqrt.(diag(J*Pt[4l+1:6l,4l+1:6l]*J'))

# remove fitted means
cargo .-= xt[1:l]'
τargo = cargo*inv(T)'

# ECCO trends
yt = reshape(cecco, l*m)
xt = Pt*Et'*(Rnnt\yt)
# travel time anomalies associated with the trend, plus standard errors
τeccotrend = reshape(Mt*Tt*xt, (m, l))
eτeccotrend = reshape(sqrt.(diag(Mt*Tt*Pt*Tt'*Mt')), (m, l))
# travel time anomalies associated with the trend, projected onto the sinular vectors, plus standard errors
ctrendecco = reshape(Mt*xt, (m, l))
ectrendecco = reshape(sqrt.(diag(Mt*Pt*Mt')), (m, l))
# trends for travel time anomalies
trendsecco = Tt[l+1:2l,:]*xt
# trends for travel time anomalies projected onto singular vectors, plus uncertainties
ctrendsecco = xt[l+1:2l]
ectrendsecco = sqrt.(diag(Pt[l+1:2l,l+1:2l]))
# annuals for travel time anomalies projected onto singular vectors, plus uncertainties
tannualsecco = xt[2l+1:4l]
cannualsecco = hypot.(tannualsecco[1:l], tannualsecco[l+1:2l])
J = spzeros(l, 2l)
for i = 1:l
  J[i,2(i-1)+1:2i] = abs.(tannualsecco[i:l:i+l])/cannualsecco[i]
end
ecannualsecco = sqrt.(diag(J*Pt[2l+1:4l,2l+1:4l]*J'))
# semiannuals for travel time anomalies projected onto singular vectors, plus uncertainties
tsemiannualsecco = xt[4l+1:6l]
csemiannualsecco = hypot.(tsemiannualsecco[1:l], tsemiannualsecco[l+1:2l])
J = spzeros(l, 2l)
for i = 1:l
  J[i,2(i-1)+1:2i] = abs.(tsemiannualsecco[i:l:i+l])/csemiannualsecco[i]
end
ecsemiannualsecco = sqrt.(diag(J*Pt[4l+1:6l,4l+1:6l]*J'))

# remove fitted means
cecco .-= xt[1:l]'
τecco = cecco*inv(T)'

# print trends and seasonal amplitudes
@printf("         %-17s %-9s %-9s\n", "      trend", "  12 mo.", "  6 mo.")
@printf("         %-17s %-9s %-9s\n", "     (mK/yr)", "   (mK)", "  (mK)")
@printf("         %-8s %-8s %-4s %-4s %-4s %-4s\n", "    1", "    2", "  1", "  2", "  1", "  2")
@printf("T-waves: %+3.1f±%3.1f %+3.1f±%3.1f %2.0f±%1.0f %2.0f±%1.0f %2.0f±%1.0f %2.0f±%1.0f\n", SOT.meanyear*1e3*reshape([ctrends[1:2]';2ectrends[1:2]'], 4)..., 1e3*reshape([cannuals[1:2]'; 2ecannuals[1:2]'], 4)..., 1e3*reshape([csemiannuals[1:2]'; 2ecsemiannuals[1:2]'], 4)...,)
@printf("Argo:    %+3.1f±%3.1f %+3.1f±%3.1f %2.0f±%1.0f %2.0f±%1.0f %2.0f±%1.0f %2.0f±%1.0f\n", SOT.meanyear*1e3*reshape([ctrendsargo[1:2]';2ectrendsargo[1:2]'], 4)..., 1e3*reshape([cannualsargo[1:2]'; 2ecannualsargo[1:2]'], 4)..., 1e3*reshape([csemiannualsargo[1:2]'; 2ecsemiannualsargo[1:2]'], 4)...,)
@printf("ECCO:    %+3.1f±%3.1f %+3.1f±%3.1f %2.0f±%1.0f %2.0f±%1.0f %2.0f±%1.0f %2.0f±%1.0f\n", SOT.meanyear*1e3*reshape([ctrendsecco[1:2]';2ectrendsecco[1:2]'], 4)..., 1e3*reshape([cannualsecco[1:2]'; 2ecannualsecco[1:2]'], 4)..., 1e3*reshape([csemiannualsecco[1:2]'; 2ecsemiannualsecco[1:2]'], 4)...,)

# temperature trend profile
ETt, RxxTt, RnnTt, PTt, MTt = estimatetrend(td, λ, [σc[1]], [σc[1]], [σtrend[1]], [σannual[1]], [σsemiannual[1]])
Ttrendargo = [(PTt*ETt'*(RnnTt\Targo[:,i]))[2] for i = 1:size(Targo, 2)]

# save catalogs of used pairs
CSV.write("results/nias_isc_DGAR_tpairs.csv", tpairs)
CSV.write("results/nias_isc_DGAR_ppairs.csv", ppairs)

# save to file
h5open("results/nias_isc_DGAR.h5", "w") do file
  write(file, "t", Dates.value.(t .- DateTime(2000, 1, 1, 0, 0, 0)))
  write(file, "ti", Dates.value.(ti .- DateTime(2000, 1, 1, 0, 0, 0)))
  write(file, "z", -depth)
  write(file, "τ", τ)
  write(file, "τi", τi)
  write(file, "τargo", τargo)
  write(file, "τecco", τecco)
  write(file, "eτ", eτ)
  write(file, "c", c)
  write(file, "ci", ci)
  write(file, "cargo", cargo)
  write(file, "cecco", cecco)
  write(file, "ec", ec)
  write(file, "ctrend", ctrend)
  write(file, "ctrendargo", ctrendargo)
  write(file, "ctrendecco", ctrendecco)
  write(file, "ectrend", ectrend)
  write(file, "ectrendargo", ectrendargo)
  write(file, "ectrendecco", ectrendecco)
  write(file, "ctrends", ctrends)
  write(file, "ctrendsargo", ctrendsargo)
  write(file, "ctrendsecco", ctrendsecco)
  write(file, "ectrends", ectrends)
  write(file, "ectrendsargo", ectrendsargo)
  write(file, "ectrendsecco", ectrendsecco)
  write(file, "Targo", Targo)
  #write(file, "Tecco", Tecco)
end

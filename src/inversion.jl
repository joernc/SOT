# TODO:
# - Use sparse QR decomposition if inversion matrix becomes too big or ill-conditioned
#   (qr([E;S])\[y;zeros(l*m)]).

# mean year in days
global const meanyear = 365.2425

"""
    collectpairs(eqname, tstations, pstations, invfreq, mincc; maxΔτ=Inf, excludetimes=[], excludepairs=[])

Collect all *T*- and *P*-wave pairs from the specified `tstations` and `pstations` that meet
certain criteria and record the measurements at the frequencies specified by `invfreq`. The
*T*-wave measurements are required to have a minimum cross-correlation coefficient of
`mincc`. *T*-wave outliers can be excluded with `maxΔτ`, time periods can be excluded with
`excludetimes`, and specific pairs can be excluded with `excludepairs`.

# Examples
```
julia> tpairs, ppairs = SOT.collectpairs("nias", ["H01W1..EDH", "H01W3..EDH"], ["PSI", "KUM", "WRAB"], [2.5, 4], [0.6, 0.4], maxΔτ=20)
[...]

julia> tpairs, ppairs = SOT.collectpairs("nias", ["H08S2..EDH"], ["PSI"], [2, 4], [0.6, 0.4]; excludetimes=[[Date(2001, 1, 1) Date(2004, 12, 1)], [Date(2010, 1, 1) Date(2012, 1, 20)]])
[...]
```
"""
function collectpairs(eqname, tstations, tintervals, tavgwidth, treffreq, tinvfreq, tmincc, pstations, pintervals, pfreqbands; tmaxΔτ=Inf, excludetimes=[], excludepairs=DataFrame(pstation=String[], event1=DateTime[], event2=DateTime[]))

  # number of frequencies at which to perform inversion
  l = length(tinvfreq)

  # load and combine catalogs of P-wave pairs
  ppairs = DataFrame[]
  for i = 1 : size(pstations, 1)
    filename = paircatfile(eqname, pstations[i], pintervals[i], pfreqbands[i])
    push!(ppairs, DataFrame(CSV.File(filename, select=1:6)))
    # exclude specified pairs
    for e = eachrow(excludepairs[excludepairs.pstation.==pstations[i],:])
      exclude = (ppairs[i].event1 .== e.event1) .& (ppairs[i].event2 .== e.event2)
      ppairs[i] = ppairs[i][.!exclude,:]
    end
  end
  ppairs = sort(unique(vcat(ppairs...)))

  # exclude events in specified time periods
  for i = 1 : length(excludetimes)
    exclude1 = excludetimes[i][1] .< ppairs.event1 .< excludetimes[i][2]
    exclude2 = excludetimes[i][1] .< ppairs.event2 .< excludetimes[i][2]
    ppairs = ppairs[.!(exclude1 .| exclude2),:]
  end

  # list of all selected T-wave pairs
  tpairs = DataFrame(station=String[], event1=DateTime[], event2=DateTime[], Δτl=Array{Float64,1}[], Δτc=Array{Float64,1}[], Δτr=Array{Float64,1}[], ccl=Array{Float64,1}[], ccc=Array{Float64,1}[], ccr=Array{Float64,1}[])

  # collect all T-wave pairs that meet the criteria
  @showprogress for i = 1 : length(tstations), j = 1 : size(ppairs, 1)

    # file with T-wave data
    tdelayfile = @sprintf("%s/%s_%s.h5", tdelaydir(eqname, tstations[i], tintervals[i], tavgwidth, treffreq), fmttime(ppairs[j,:event1]), fmttime(ppairs[j,:event2]))

    # read data if present
    if isfile(tdelayfile) && filesize(tdelayfile) > 0

      # open file
      fid = h5open(tdelayfile, "r")

      # frequencies
      tfreq = read(fid, "freq")

      # find frequency indices
      idx = [argmin(abs.(tfreq .- tinvfreq[i])) for i = 1:l]

      # read central delay and CC
      Δτc = read(fid, "Δτc")[idx]
      ccc = read(fid, "ccc")[idx]

      # check if criteria are met
      if all(ccc .> tmincc) && all(abs.(Δτc) .< tmaxΔτ)

        # record measurements
        push!(tpairs, [tstations[i], ppairs[j,:event1], ppairs[j,:event2], read(fid, "Δτl")[idx], Δτc, read(fid, "Δτr")[idx], read(fid, "ccl")[idx], ccc, read(fid, "ccr")[idx]])

      end

      # close file
      close(fid)

    end

  end

  # list of all selected P-wave pairs
  ppairs = DataFrame(station=String[], event1=DateTime[], latitude1=Float64[], longitude1=Float64[], event2=DateTime[], latitude2=Float64[], longitude2=Float64[], Δτ=Float64[], cc=Float64[])

  # find all pairs that were selected
  for i = 1 : length(pstations)

    # P-wave catalog
    filename = paircatfile(eqname, pstations[i], pintervals[i], pfreqbands[i])
    pairs = DataFrame(CSV.File(filename))

    # find pairs that were selected based on T-wave measurements
    selpairs = innerjoin(pairs, select(tpairs, [:event1, :event2]), on=[:event1, :event2])

    # exclude specified pairs
    for e = eachrow(excludepairs[excludepairs.pstation.==pstations[i],:])
      exclude = (selpairs.event1 .== e.event1) .& (selpairs.event2 .== e.event2)
      selpairs = selpairs[.!exclude,:]
    end

    # add station information
    selpairs.station = repeat([pstations[i]], size(selpairs, 1))

    # record
    append!(ppairs, selpairs)

  end

  # return pairs
  return tpairs, ppairs

end

"""
    invert(tpairs, ppairs; timescale=NaN)

Set up inversion of the travel time changes between events ``Δτ`` for travel time anomalies
``τ`` relative to an arbitrary but common reference. Returns the unique events `t`, the
design matrix `E`, the smoothing matrix `S`, the pseudoinverse `P = pinv(E'*E + S'*S)`, and
the difference matrix `D` that take the difference between *T*- and *P*-wave anomalies. If
`timescale` is set (in days), smoothing is applied at that scale. By default
(`timescale=NaN`), smoothing is applied at the sampling scale. (See Wunsch: Discrete Inverse
and State Estimation Problems, p. 56.)

From these matrices, the travel time anomalies can be computed as
```
# number of T- and P-wave pairs
nt = size(tpairs, 1)
np = size(ppairs, 1)
# number of unique events
m = length(t)
# number of frequencies
l = length(tpairs.Δτc[1])
# data vector
y = [reshape(vcat(tpairs.Δτ'...), l*nt); ppairs.Δτ]
# least-square solution
x = P*(E'*y)
# travel time anomalies as differences between T- and P-wave anomalies
τ = reshape(D*x, (m, l))
```
The standard error of the solution can then be obtained using
```
# residuals
n = y - E*x
# estimate variances for T- and P-wave delays
σt = sqrt(mean(n[1:l*nt].^2)*(l*nt+np)/(l*nt+np-l*m-1))
σp = sqrt(mean(n[l*nt+1:l*nt+np].^2)*(l*nt+np)/(l*nt+np-l*m-1))
R = spdiagm(0 => [σt^2*ones(l*nt); σp^2*ones(np)])
# propagate error
A = P*E'
τerr = reshape(sqrt.(diag(D*A*R*A'*D')), (m, l))
```
  λ: correlation time (days)
  σx: solution standard deviation
  σn: measurement noise (s)
  σp: origin time error (s)
"""
function invert(tpairs, ppairs, λ, σc, σn, σp, U, Λ, Δz, h; σtrend=0, σannual=0, σsemiannual=0)

  # number of frequencies
  l = length(tpairs.Δτc[1])

  # find unique events
  t = sort(unique([tpairs.event1; tpairs.event2]))

  # real time (days)
  tr = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24

  # number of T- and P-wave pairs
  nt = size(tpairs, 1)
  np = size(ppairs, 1)

  # number of unique events
  m = length(t)

  # T-wave pair matrix
  tidx1 = indexin(tpairs.event1, t)
  tidx2 = indexin(tpairs.event2, t)
  Xt = sparse([1:nt; 1:nt], [tidx1; tidx2], [-ones(nt); ones(nt)])

  # P-wave pair matrix
  pidx1 = indexin(ppairs.event1, t)
  pidx2 = indexin(ppairs.event2, t)
  Xp = sparse([1:np; 1:np], [pidx1; pidx2], [-ones(np); ones(np)])

  # full design matrix
  E = blockdiag([Xt for i = 1:l]..., Xp)

  # difference matrix (T- minus P-anomalies)
  D = [I(l*m) vcat([-I(m) for i = 1:l]...)]

  # solution covariance in time
  C = exp.(-abs.(tr.-tr')/λ)

  # transformation from singular vectors to observed frequencies
  T = h*kron(sparse(U*Diagonal(Λ)), I(m))

  # covariance matrix assuming no correlation between singular vector expansion coefficients
  R = [T*kron(spdiagm(σc.^2), C)*T' zeros(l*m, m); zeros(m, (l+1)*m)] + σp^2*kron(sparse(ones(l+1, l+1)), I(m))

  # add trend if desired
  if !(σtrend == 0)
    tm = tr[1]+(tr[m]-tr[1])/2
    E = [E [kron(I(l), Xt*(tr.-tm)); zeros(np, l)]]
    T = h*U*Diagonal(Λ)
    R = [R zeros(size(R, 1), l); zeros(l, size(R, 2)) T*spdiagm(σtrend.^2)*T']
    D = [D kron(I(l), tr.-tm)]
  end

  # add annual cycle if desired
  if !(σannual == 0)
    ω = 2π/SOT.meanyear
    E = [E [kron(I(l), Xt*cos.(ω*tr)) kron(I(l), Xt*sin.(ω*tr)); zeros(np, 2l)]]
    T = h*U*Diagonal(Λ)
    R = [R zeros(size(R, 1), 2l); zeros(2l, size(R, 2)) kron(I(2), T*spdiagm(σannual.^2/2)*T')]
    D = [D kron(I(l), cos.(ω*tr)) kron(I(l), sin.(ω*tr))]
  end

  # add semi-annual cycle if desired
  if !(σsemiannual == 0)
    ω = 4π/SOT.meanyear
    E = [E [kron(I(l), Xt*cos.(ω*tr)) kron(I(l), Xt*sin.(ω*tr)); zeros(np, 2l)]]
    T = h*U*Diagonal(Λ)
    R = [R zeros(size(R, 1), 2l); zeros(2l, size(R, 2)) kron(I(2), T*spdiagm(σsemiannual.^2/2)*T')]
    D = [D kron(I(l), cos.(ω*tr)) kron(I(l), sin.(ω*tr))]
  end

  # noise covariance
  N = σn^2*I

  # Gauss–Markov solution covariance matrix
  P = inv(inv(Array(R)) + E'*inv(N)*E)

  # Best estimate, subtracting T- and P-delays, adding trend and seasonal cycles:
  #   τ = D*P*E'*inv(N)*y
  # Uncertainty:
  #   e = sqrt.(diag(D*P*D'))

  return t, E, R, N, P, D

end

### with model covariance
function invert_gpcov(tpairs, ppairs, tstation, pstations, evtpos, σc,σp,σs,σn,σnp,σh,lags,ctau; σtrend=0, σannual=0, σsemiannual=0)

  # number of frequencies
  l = length(tpairs.Δτc[1])

  # find unique events
  t = sort(unique([tpairs.event1; tpairs.event2]))

  # real time (days)
  tr = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24

  # number of T- and P-wave pairs
  nt = size(tpairs, 1)
  np = size(ppairs, 1)

  # number of unique events
  m = length(t)

  # T-wave pair matrix
  tidx1 = indexin(tpairs.event1, t)
  tidx2 = indexin(tpairs.event2, t)
  Xt = sparse([1:nt; 1:nt], [tidx1; tidx2], [-ones(nt); ones(nt)])

  # P-wave pair matrix
  pidx1 = indexin(ppairs.event1, t)
  pidx2 = indexin(ppairs.event2, t)
  Xp = sparse([1:np; 1:np], [pidx1; pidx2], [-ones(np); ones(np)])
  
  ppairs.avglat = 0.5(ppairs.latitude1+ppairs.latitude2)
  ppairs.avglon = 0.5(ppairs.longitude1+ppairs.longitude2)
  leftjoin!(ppairs, pstations, on=:station)
  ppazms = azimuth.(ppairs.slat,ppairs.slon,ppairs.avglat,ppairs.avglon)
  #@printf("selpairs: %s\n%s\n%s\n%s\n",selpairs.slat[1:5],selpairs.slon[1:5],selpairs.avglat[1:5],selpairs.avglon[1:5])
  Xpcs = sparse([1:np; 1:np], [pidx1; pidx2], [-cosd.(ppazms); cosd.(ppazms)])
  Xpsn = sparse([1:np; 1:np], [pidx1; pidx2], [-sind.(ppazms); sind.(ppazms)])
  ppairs = select(ppairs, [:event1, :event2, :avglat, :avglon])
  unique!(ppairs)
  tpairs = select(tpairs, [:event1, :event2])
  leftjoin!(tpairs, ppairs, on=[:event1, :event2])
  tpazms = azimuth.(tstation[1],tstation[2],tpairs.avglat,tpairs.avglon)
  Xtcs = sparse([1:nt; 1:nt], [tidx1; tidx2], [-cosd.(tpazms); cosd.(tpazms)])
  Xtsn = sparse([1:nt; 1:nt], [tidx1; tidx2], [-sind.(tpazms); sind.(tpazms)])

  Rcsn = [[Xtcs for i = 1:l]...;Xpcs]*[[Xtcs for i = 1:l]...;Xpcs]'
  Rcsn = Rcsn + [[Xtsn for i = 1:l]...;Xpsn]*[[Xtsn for i = 1:l]...;Xpsn]'
  Rcsnh = [[Xtcs for i = 1:l]...;0*Xpcs]*[[Xtcs for i = 1:l]...;0*Xpcs]'
  Rcsnh = Rcsnh + [[Xtsn for i = 1:l]...;0*Xpsn]*[[Xtsn for i = 1:l]...;0*Xpsn]'

  # full design matrix
  E = blockdiag([Xt for i = 1:l]..., Xp)

  # difference matrix (T- minus P-anomalies)
  D = [I(l*m) vcat([-I(m) for i = 1:l]...)]

  # solution covariance in time
  Δtrd = tr.-tr'
  C = []
  for i = 1:6
    cij = ctau[:,i]
    interp_linear = linear_interpolation(lags, cij, extrapolation_bc=0)
    push!(C,interp_linear.(Δtrd))
  end
  C = Symmetric([σc[1]^2*C[1] σc[2]^2*C[2] σc[3]^2*C[3]; σc[2]^2*C[2]' σc[4]^2*C[4] σc[5]^2*C[5]; σc[3]^2*C[3]' σc[5]^2*C[5]' σc[6]^2*C[6]])
  #C = Matrix([C[1] C[2] C[3]; C[2]' C[4] C[5]; C[3]' C[5]' C[6]])
  #C = (C+C')/2

  # covariance matrix assuming no correlation between singular vector expansion coefficients
  R = [C zeros(l*m, m); zeros(m, (l+1)*m)] + σp^2*kron(sparse(ones(l+1, l+1)), I(m))

  # add trend if desired
  if !(σtrend == 0)
    tm = tr[1]+(tr[m]-tr[1])/2
    E = [E [kron(I(l), Xt*(tr.-tm)); zeros(np, l)]]
    R = [R zeros(size(R, 1), l); zeros(l, size(R, 2)) spdiagm(σtrend.^2)]
    D = [D kron(I(l), tr.-tm)]
  end

  # add annual cycle if desired
  if !(σannual == 0)
    ω = 2π/SOT.meanyear
    E = [E [kron(I(l), Xt*cos.(ω*tr)) kron(I(l), Xt*sin.(ω*tr)); zeros(np, 2l)]]
    R = [R zeros(size(R, 1), 2l); zeros(2l, size(R, 2)) kron(I(2), spdiagm(σannual.^2/2))]
    D = [D kron(I(l), cos.(ω*tr)) kron(I(l), sin.(ω*tr))]
  end

  # add semi-annual cycle if desired
  if !(σsemiannual == 0)
    ω = 4π/SOT.meanyear
    E = [E [kron(I(l), Xt*cos.(ω*tr)) kron(I(l), Xt*sin.(ω*tr)); zeros(np, 2l)]]
    R = [R zeros(size(R, 1), 2l); zeros(2l, size(R, 2)) kron(I(2), spdiagm(σsemiannual.^2/2))]
    D = [D kron(I(l), cos.(ω*tr)) kron(I(l), sin.(ω*tr))]
  end

  # noise covariance
  N = σn^2*spdiagm([ones(l*nt);zeros(np)])+σs^2*Rcsn+σh^2*Rcsnh
  N += σnp^2*spdiagm([zeros(l*nt);ones(np)])

  invR = inv(Array(R))
  # Gauss–Markov solution covariance matrix
  P = inv(invR + E'*inv(Array(N))*E)
  #P1 = E'*inv((Array(N)) + E*inv(Array(R))*E')*E

  # Best estimate, subtracting T- and P-delays, adding trend and seasonal cycles:
  #   τ = D*P*E'*inv(N)*y
  # Uncertainty:
  #   e = sqrt.(diag(D*P*D'))

  return t, E, R, N, P, D, invR

end

### with model covariance and 1 frequency
function invert_gpcov1f(tpairs, ppairs, tstation, pstations, evtpos, λt, λθ, στ, σx, σh, σn, σnp, σp; σtrend=0, σannual=0, σsemiannual=0, σθtrend=0)

  # number of frequencies
  l = 1

  stalon,stalat,evtlon,evtlat = tstation[2],tstation[1],evtpos[2],evtpos[1]

  # find unique events
  t = sort(unique([tpairs.event1; tpairs.event2]))

  # real time (days)
  tr = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24

  # number of T- and P-wave pairs
  nt = size(tpairs, 1)
  np = size(ppairs, 1)

  # number of unique events
  m = length(t)
  
  # T-wave pair matrix
  tidx1 = indexin(tpairs.event1, t)
  tidx2 = indexin(tpairs.event2, t)
  Xt = sparse([1:nt; 1:nt], [tidx1; tidx2], [-ones(nt); ones(nt)])

  # P-wave pair matrix
  pidx1 = indexin(ppairs.event1, t)
  pidx2 = indexin(ppairs.event2, t)
  Xp = sparse([1:np; 1:np], [pidx1; pidx2], [-ones(np); ones(np)])

  selpairs = select(ppairs, [:station, :event1, :event2])
  selpairs.avglat = 0.5(ppairs.latitude1+ppairs.latitude2)
  selpairs.avglon = 0.5(ppairs.longitude1+ppairs.longitude2)
  leftjoin!(selpairs, pstations, on=:station)
  ppazms = azimuth.(selpairs.slat,selpairs.slon,selpairs.avglat,selpairs.avglon)
  Xpcs = sparse([1:np; 1:np], [pidx1; pidx2], [-cosd.(ppazms); cosd.(ppazms)])
  Xpsn = sparse([1:np; 1:np], [pidx1; pidx2], [-sind.(ppazms); sind.(ppazms)])
  selpairs = select(selpairs, [:event1, :event2, :avglat, :avglon])
  unique!(selpairs)
  tpairs = select(tpairs, [:event1, :event2])
  leftjoin!(tpairs, selpairs, on=[:event1, :event2])
  tpazms = azimuth.(tstation[1],tstation[2],tpairs.avglat,tpairs.avglon)
  Xtcs = sparse([1:nt; 1:nt], [tidx1; tidx2], [-cosd.(tpazms); cosd.(tpazms)])
  Xtsn = sparse([1:nt; 1:nt], [tidx1; tidx2], [-sind.(tpazms); sind.(tpazms)])
  Ecs,Esn = [Xtcs;Xpcs],[Xtsn;Xpsn]
  Ecsh,Esnh = [Xtcs;zeros(size(Xp))],[Xtsn;zeros(size(Xp))]
  
  # full design matrix
  E = blockdiag([Xt for i = 1:l]..., Xp)
  
  # difference matrix (T- minus P-anomalies)
  D = [I(l*m) vcat([-I(m) for i = 1:l]...)]

  # solution covariance in time
  if λθ>0
    xd,yd,θd = getazimuth(tpairs,ppairs,stalon,stalat,evtlon,evtlat)
    C = exp.(-abs.(tr.-tr')/λt-(θd.-θd').^2/λθ.^2/2)
  else
    C = exp.(-abs.(tr.-tr')/λt)
  end
  
  # transformation from singular vectors to observed frequencies
  T = I#h*kron(sparse(U*Diagonal(Λ)), I(m))

  # covariance matrix assuming no correlation between singular vector expansion coefficients
  R = [T*kron((στ.^2), C)*T' zeros(l*m, m); zeros(m, (l+1)*m)] + σp^2*kron(sparse(ones(l+1, l+1)), I(m))

  # add trend if desired
  if !(σtrend == 0)
    tm = tr[1]+(tr[m]-tr[1])/2
    E = [E [kron(I(l), Xt*(tr.-tm)); zeros(np, l)]]
    T = I#h*U*Diagonal(Λ)
    R = [R zeros(size(R, 1), l); zeros(l, size(R, 2)) T*(σtrend.^2)*T']
    D = [D kron(I(l), tr.-tm)]
  end

  # add annual cycle if desired
  if !(σannual == 0)
    ω = 2π/SOT.meanyear
    E = [E [kron(I(l), Xt*cos.(ω*tr)) kron(I(l), Xt*sin.(ω*tr)); zeros(np, 2l)]]
    T = I#h*U*Diagonal(Λ)
    R = [R zeros(size(R, 1), 2l); zeros(2l, size(R, 2)) I(2)*(σannual.^2/2)]
    D = [D kron(I(l), cos.(ω*tr)) kron(I(l), sin.(ω*tr))]
  end

  # add semi-annual cycle if desired
  if !(σsemiannual == 0)
    ω = 4π/SOT.meanyear
    E = [E [kron(I(l), Xt*cos.(ω*tr)) kron(I(l), Xt*sin.(ω*tr)); zeros(np, 2l)]]
    T = I#h*U*Diagonal(Λ)
    R = [R zeros(size(R, 1), 2l); zeros(2l, size(R, 2)) I(2)*(σsemiannual.^2/2)]
    D = [D kron(I(l), cos.(ω*tr)) kron(I(l), sin.(ω*tr))]
  end
  
  # add azm trend if desired
  if !(σθtrend == 0)
    E = [E [kron(I(l), Xt*θd); zeros(np, l)]]
    T = I
    R = [R zeros(size(R, 1), l); zeros(l, size(R, 2)) T*(σθtrend.^2)*T']
    D = [D kron(I(l), θd)]
  end

  # noise covariance
  RcsnEa = Ecs*Ecs'+Esn*Esn'
  Rcsnh = Ecsh*Ecsh'+Esnh*Esnh'
  N = σn^2*spdiagm([ones(nt);zeros(np)])+σx^2*RcsnEa+σh^2*Rcsnh
  N += σnp^2*spdiagm([zeros(nt);ones(np)])

  #invR = inv(Array(R))
  invR = inv(E*Array(R)*E'+Array(N))
  # Gauss–Markov solution covariance matrix
  P = inv(inv(Array(R)) + E'*inv(Array(N))*E)

  # Best estimate, subtracting T- and P-delays, adding trend and seasonal cycles:
  #   τ = D*P*E'*inv(N)*y
  # Uncertainty:
  #   e = sqrt.(diag(D*P*D'))
  if λθ>0
    return t, xd, yd, θd, E, R, N, P, D, invR
  else
    return t, E, R, N, P, D, invR
  end

end

"""
    correctcycleskipping(tpairs, ppairs, E, S, P)

Find cycle-skipping corrections. Applies corrections to adjacent maxima of the
cross-correlation function whenever they reduce the cost function. Returns the corrected
*T*-wave delays `Δτ`.
"""
function correctcycleskipping(tpairs, ppairs, E, Rxx, Rnn, P, m)
  
  # number of frequencies
  l = length(tpairs.Δτc[1])

  # number of T- and P-wave pairs
  nt = size(tpairs, 1)
  np = size(ppairs, 1)

  # collect the three delays and CCs
  Δτa = cat(vcat(tpairs.Δτl'...), vcat(tpairs.Δτc'...), vcat(tpairs.Δτr'...), dims=3)
  cca = cat(vcat(tpairs.ccl'...), vcat(tpairs.ccc'...), vcat(tpairs.ccr'...), dims=3)

  # invert for P-wave delays without smoothing
  Δτp = E[1:nt,1:m]*(E[l*nt+1:l*nt+np,l*m+1:(l+1)*m]\ppairs.Δτ)

  # use Gaussian mixture model with three members and shared covariance to find three
  # clusters of pairs, estimate parameters using an EM algorithm, perform initial cycle-
  # skipping correction by shifting left cluster to right and right cluster to left
  x1 = [tpairs.Δτc[i][1] - Δτp[i] for i = 1:nt]
  x2 = [tpairs.Δτc[i][1] - tpairs.Δτc[i][l] for i = 1:nt]
  x = [x1 x2]
  μ = [1 0.1; 0 0; -1 -0.1]
  Λ = [0.2^2 0; 0 0.1^2]
  τ = [0.05, 0.9, 0.05]
  f(x, μ, Λ) = 1/sqrt((2π)^length(μ)*norm(Λ))*exp(-1/2*(x-μ)'*(Λ\(x-μ)))
  T = nothing
  for n = 1:100
    T = [τ[j].*f(x[i,:], μ[j,:], Λ) for j=1:3, i=1:nt]
    T ./= sum(T; dims=1)
    τ = sum(T; dims=2)/nt
    μ = sum([T[j,i]*x[i,k] for j=1:3, k=1:2, i=1:nt]; dims=3)[:,:,1]./sum(T; dims=2)
    Λ = sum([T[j,i]*(x[i,:]-μ[j,:])*(x[i,:]-μ[j,:])' for i = 1:nt, j = 1:3])/nt
  end
  cs = [argmax(T[:,i]) for i = 1:nt]

  # plot clusters used for initial correction
  rc("font", size=8)
  rc("axes", titlesize="medium")
  fig, ax = subplots(1, 2, sharex=true, sharey=true, figsize=(190/25.4, 95/25.4))
  for i = 1:3
    ax[1].scatter(x[cs.==i,1], x[cs.==i,2], s=5)
  end
  xlim = round(maximum(abs.(ax[1].get_xlim())); digits=2)
  ylim = round(maximum(abs.(ax[1].get_ylim())); digits=2)
  ax[1].set_xlim(-xlim, xlim)
  ax[1].set_ylim(-ylim, ylim)
  ax[1].set_xlabel("low-frequency delay (s)")
  ax[1].set_ylabel("differential delay (s)")
  ax[1].set_title("after clustering")
  ax[1].set_title("(a)", loc="left")

  for i = 1:nt
    if cs[i] != 2
      s = @sprintf("%s %s\n", tpairs.event1[i], tpairs.event2[i])
      printstyled(s, color=(cs[i]==1) ? :blue : :green)
    end
  end

  # data vector using initial correction
  y = [reshape([Δτa[i,j,cs[i]] for i = 1:nt, j = 1:l], l*nt); ppairs.Δτ]

  # inverses
  iRnn = inv(Rnn)
  iRxx = inv(Array(Rxx))

  # solution operator (x̃ = Ay)
  A = P*E'*iRnn

#  # residual operator
#  R = I - E*A

  # initial residuals and cost
  x = A*y
  n = y - E*x
  r = n.^2
  J = .5n'*iRnn*n + .5x'*iRxx*x

  # cycle through T-wave pairs until no further corrections are made
  i = 1
  idx = nothing
  while i ≤ nt

    @printf("\r%4d", i)

    # sort unique pairs by size of residuals
    if i == 1
      idx = sortperm(r[1:nt], rev=true)
    end

    # corrected this pair?
    corrected = false

    # trial correction(s)
    if cs[idx[i]] == 1 || cs[idx[i]] == 3
      dir = [2]
    else
      dir = [1, 3]
    end

    # go through trials
    for j = dir

      # check if probability of belonging to different cluster is at least 0.1%
      if T[j,idx[i]] ≥ 1e-3

        # swap out Δτ
        y[idx[i].+(0:l-1)*nt] = Δτa[idx[i],:,j]

        # new residuals and cost
        x = A*y
        n = y - E*x
        rp = n.^2
        Jp = .5n'*iRnn*n + .5x'*iRxx*x

        # record if cost is reduced, revert otherwise
        if Jp < J
          cs[idx[i]] = j
          J = Jp
          r = copy(rp)
          corrected = true
          @printf("\r%4d %4d %d %7.5f %s %s\n", i, idx[i], j, 1e3J/(nt+np+m-2)/l,
                  tpairs.event1[idx[i]], tpairs.event2[idx[i]])
        else
          y[idx[i].+(0:l-1)*nt] = Δτa[idx[i],:,cs[idx[i]]]
        end

      end
    end

    # move back to first pair if correction was made, advance otherwise
    i = corrected ? 1 : i+1

  end

  @printf("\n")
  @printf("Total number of T-wave pairs:            %4d\n", nt)
  @printf("Number of pairs corrected to left max.:  %4d\n", sum(cs.==1))
  @printf("Number of uncorrected pairs:             %4d\n", sum(cs.==2))
  @printf("Number of pairs corrected to right max.: %4d\n", sum(cs.==3))

  # plot clusters used for initial correction
  x1 = [tpairs.Δτc[i][1] - Δτp[i] for i = 1:nt]
  x2 = [tpairs.Δτc[i][1] - tpairs.Δτc[i][l] for i = 1:nt]
  x = [x1 x2]
  for i = 1:3
    ax[2].scatter(x[cs.==i,1], x[cs.==i,2], s=5)
  end
  ax[2].set_xlabel("low-frequency delay (s)")
  ax[2].set_title("after additional corrections")
  ax[2].set_title("(b)", loc="left")
  fig.tight_layout(w_pad=2)

  # return optimal delays
  return collect(eachrow(reshape(y[1:l*nt], (nt, l))))

end

### return normalized residuals z
function correctcycleskipping_returnz(eqname, tstations, tpairs, ppairs, E, Rxx, Rnn, invR, P, m)
  
  # number of frequencies
  l = length(tpairs.Δτc[1])

  # number of T- and P-wave pairs
  nt = size(tpairs, 1)
  np = size(ppairs, 1)

  # collect the three delays and CCs
  Δτa = cat(vcat(tpairs.Δτl'...), vcat(tpairs.Δτc'...), vcat(tpairs.Δτr'...), dims=3)
  cca = cat(vcat(tpairs.ccl'...), vcat(tpairs.ccc'...), vcat(tpairs.ccr'...), dims=3)
  
  # inverses
  iRnn = inv(Array(Rnn))
  #iRxx = inv(Symmetric(Rxx))
  #iRyy = inv(Array(E*Rxx*E' + Rnn))

  # solution operator (x̃ = Ay)
  A = P*E'*iRnn

  # invert for P-wave delays without smoothing
  @printf("origin correction w/ covariance...\n") 
  invRp = inv(Array(Rxx[l*m+1:(l+1)*m,l*m+1:(l+1)*m]))
  invNp = inv(Array(Rnn[l*nt+1:end,l*nt+1:end]))
  # Gauss–Markov solution covariance matrix
  Pp = inv(invRp + E[l*nt+1:end,l*m+1:(l+1)*m]'*invNp*E[l*nt+1:end,l*m+1:(l+1)*m])
  τp = Pp*E[l*nt+1:end,l*m+1:(l+1)*m]'*invNp*ppairs.Δτ
  Δτp = E[1:nt,1:m]*τp
#  Δτp = E[1:nt,1:m]*(E[l*nt+1:l*nt+np,l*m+1:(l+1)*m]\ppairs.Δτ)

  # use Gaussian mixture model with three members and shared covariance to find three
  # clusters of pairs, estimate parameters using an EM algorithm, perform initial cycle-
  # skipping correction by shifting left cluster to right and right cluster to left
  x1 = [tpairs.Δτc[i][1] - Δτp[i] for i = 1:nt]
  x2 = [tpairs.Δτc[i][1] - tpairs.Δτc[i][l] for i = 1:nt]
  x = [x1 x2]
  μ = [1 0.1; 0 0; -1 -0.1]
  Λ = [0.2^2 0; 0 0.1^2]
  τ = [0.05, 0.9, 0.05]
  f(x, μ, Λ) = 1/sqrt((2π)^length(μ)*norm(Λ))*exp(-1/2*(x-μ)'*(Λ\(x-μ)))
  T = nothing
  for n = 1:100
    T = [τ[j].*f(x[i,:], μ[j,:], Λ) for j=1:3, i=1:nt]
    T ./= sum(T; dims=1)
    τ = sum(T; dims=2)/nt
    μ = sum([T[j,i]*x[i,k] for j=1:3, k=1:2, i=1:nt]; dims=3)[:,:,1]./sum(T; dims=2)
    Λ = sum([T[j,i]*(x[i,:]-μ[j,:])*(x[i,:]-μ[j,:])' for i = 1:nt, j = 1:3])/nt
  end
  cs = [argmax(T[:,i]) for i = 1:nt]

  # plot clusters used for initial correction
  rc("font", size=8)
  rc("axes", titlesize="medium")
  fig, ax = subplots(1, 2, sharex=true, sharey=true, figsize=(190/25.4, 95/25.4))
  for i = 1:3
    ax[1].scatter(x[cs.==i,1], x[cs.==i,2], s=5)
  end
  xlim = round(maximum(abs.(ax[1].get_xlim())); digits=2)
  ylim = round(maximum(abs.(ax[1].get_ylim())); digits=2)
  ax[1].set_xlim(-xlim, xlim)
  ax[1].set_ylim(-ylim, ylim)
  ax[1].set_xlabel("low-frequency delay (s)")
  ax[1].set_ylabel("differential delay (s)")
  ax[1].set_title("after clustering")
  ax[1].set_title("(a)", loc="left")

  for i = 1:nt
    if cs[i] != 2
      s = @sprintf("%s %s\n", tpairs.event1[i], tpairs.event2[i])
      printstyled(s, color=(cs[i]==1) ? :blue : :green)
    end
  end

  # data vector using initial correction
  y = [reshape([Δτa[i,j,cs[i]] for i = 1:nt, j = 1:l], l*nt); ppairs.Δτ]

#  # residual operator
#  R = I - E*A

  # initial residuals and cost
  x = A*y
  n = y - E*x
  r = n.^2
  J = .5n'*iRnn*n + .5x'*invR*x

  # cycle through T-wave pairs until no further corrections are made
  i = 1
  idx = nothing
  while i ≤ nt

    @printf("\r%4d", i)

    # sort unique pairs by size of residuals
    if i == 1
      idx = sortperm(r[1:nt], rev=true)
    end

    # corrected this pair?
    corrected = false

    # trial correction(s)
    if cs[idx[i]] == 1 || cs[idx[i]] == 3
      dir = [2]
    else
      dir = [1, 3]
    end

    # go through trials
    for j = dir

      # check if probability of belonging to different cluster is at least 0.1%
      if T[j,idx[i]] ≥ 1e-3

        # swap out Δτ
        y[idx[i].+(0:l-1)*nt] = Δτa[idx[i],:,j]

        # new residuals and cost
        x = A*y
        n = y - E*x
        rp = n.^2
        Jp = .5n'*iRnn*n + .5x'*invR*x

        # record if cost is reduced, revert otherwise
        if Jp < J
          cs[idx[i]] = j
          J = Jp
          r = copy(rp)
          corrected = true
          @printf("\r%4d %4d %d %7.5f %s %s\n", i, idx[i], j, 1e3J/(nt+np+m-2)/l,
                  tpairs.event1[idx[i]], tpairs.event2[idx[i]])
        else
          y[idx[i].+(0:l-1)*nt] = Δτa[idx[i],:,cs[idx[i]]]
        end

      end
    end

    # move back to first pair if correction was made, advance otherwise
    i = corrected ? 1 : i+1

  end

  @printf("\n")
  @printf("Total number of T-wave pairs:            %4d\n", nt)
  @printf("Number of pairs corrected to left max.:  %4d\n", sum(cs.==1))
  @printf("Number of uncorrected pairs:             %4d\n", sum(cs.==2))
  @printf("Number of pairs corrected to right max.: %4d\n", sum(cs.==3))

  # plot clusters used for initial correction
  Δτp = E[1:nt,1:m]*((A*y)[l*m+1:(l+1)*m])
  x1 = [tpairs.Δτc[i][1] - Δτp[i] for i = 1:nt]
  x2 = [tpairs.Δτc[i][1] - tpairs.Δτc[i][l] for i = 1:nt]
  x = [x1 x2]
  for i = 1:3
    ax[2].scatter(x[cs.==i,1], x[cs.==i,2], s=5)
  end
  ax[2].set_xlabel("low-frequency delay (s)")
  ax[2].set_title("after additional corrections")
  ax[2].set_title("(b)", loc="left")
  fig.tight_layout(w_pad=2)
  fig.savefig(@sprintf("results/plots/%s_%s_cluster.pdf",eqname,tstations[1]))
  
  #Δy = y[1:l*nt] .- repeat(Δτp, outer = l)
  #yh = E*A*y
  #Δτph = E[1:nt,1:m]*(E[l*nt+1:l*nt+np,l*m+1:(l+1)*m]\yh[(l*nt+1):end])
  #Δr = (Δy .- (yh[1:l*nt] .- repeat(Δτph, outer = l))).^2
  x = A*y
  n = y - E*x
  #T = I(l*nt+np)-E*A
  
  #Cn = T*(Rnn+E*Rxx*E')*T'
  iCnt = inv(Symmetric(Matrix(Rnn[1:l*nt,1:l*nt])))
  iCnp = inv(Symmetric(Matrix(Rnn[l*nt+1:end,l*nt+1:end])))
  #U,S,Vt = svd(Symmetric(Rxx-P)[1+l*m:(l+1)*m,1+l*m:(l+1)*m])
  #iCo = inv(Matrix(Diagonal(diag(Rxx-P)[1+l*m:(l+1)*m])))
  @printf("iCnt posdef: %s, iCnp posdef: %s\n",isposdef(iCnt),isposdef(iCnp))
  zint = isposdef(iCnt) ? cholesky(iCnt).U*n[1:l*nt] : zeros(l*nt)
  zinp = isposdef(iCnp) ? cholesky(iCnp).U*n[1+l*nt:end] : zeros(np)
  #zio = Diagonal(sqrt.(S.^-1))*U'*x[1+l*m:(l+1)*m]

  # return optimal delays
  return collect(eachrow(reshape(y[1:l*nt], (nt, l)))), Δτp, zint, zinp

end

### return normalized residuals z with 1 frequency and 2 possible right cycle skips
function correctcycleskipping_1f2r(eqname, tstations, tpairs, ppairs, E, Rxx, Rnn, P, m; l2=false)
  
  # number of frequencies
  l = 1#length(tpairs.Δτc[1])

  # number of T- and P-wave pairs
  nt = size(tpairs, 1)
  np = size(ppairs, 1)

  # collect the three delays and CCs
  if l2
    Δτa = cat(vcat(tpairs.Δτl'...), vcat(tpairs.Δτc'...), vcat(tpairs.Δτr'...), vcat(tpairs.Δτr2'...), vcat(tpairs.Δτl2'...), dims=3)
  else
    Δτa = cat(vcat(tpairs.Δτl'...), vcat(tpairs.Δτc'...), vcat(tpairs.Δτr'...), vcat(tpairs.Δτr2'...), dims=3)
  end
  cca = cat(vcat(tpairs.ccl'...), vcat(tpairs.ccc'...), vcat(tpairs.ccr'...), dims=3)

  # invert for P-wave delays without smoothing
  @printf("origin correction w/ covariance...\n") 
  invRp = inv(Array(Rxx[l*m+1:(l+1)*m,l*m+1:(l+1)*m]))
  invNp = inv(Array(Rnn[nt+1:end,nt+1:end]))
  # Gauss–Markov solution covariance matrix
  Pp = inv(invRp + E[nt+1:end,l*m+1:(l+1)*m]'*invNp*E[nt+1:end,l*m+1:(l+1)*m])
  τp = Pp*E[nt+1:end,l*m+1:(l+1)*m]'*invNp*ppairs.Δτ
  Δτp = E[1:nt,1:m]*τp
  Δτp0 = E[nt+1:end,m+1:2m]*τp
  fig, ax = subplots(1, 2)
  ax[1].hist(ppairs.Δτ-Δτp0;bins=50)
  ax[2].hist(Δτp;bins=50)
  fig.tight_layout()
  fig.savefig(@sprintf("results/plots/japan_%s_origincorrection_hist.pdf",tstations[1]))
  #Δτp = E[1:nt,1:m]*(E[l*nt+1:l*nt+np,l*m+1:(l+1)*m]\ppairs.Δτ)

  # use Gaussian mixture model with three members and shared covariance to find three
  # clusters of pairs, estimate parameters using an EM algorithm, perform initial cycle-
  # skipping correction by shifting left cluster to right and right cluster to left
  x1 = [tpairs.Δτc[i][1] - Δτp[i] for i = 1:nt]
  x2 = [tpairs.Δτc[i][1] - tpairs.Δτc[i][2] for i = 1:nt]
  x = [x1 x2]
  if l2
    jmax = 5
    μ = [0.9 0.01; 0 0; -0.9 -0.01; -1.8 -0.02; 1.8 0.02]
    τ = [0.05, 0.68, 0.23, 0.02, 0.02]
  else
    jmax = 4
    μ = [0.9 0.01; 0 0; -0.9 -0.01; -1.8 -0.02]
    τ = [0.05, 0.70, 0.23, 0.02]
  end
  Λ = [0.2^2 0; 0 0.01^2]
  f(x, μ, Λ) = 1/sqrt((2π)^length(μ)*norm(Λ))*exp(-1/2*(x-μ)'*(Λ\(x-μ)))
  T = nothing
  for n = 1:100
    T = [τ[j].*f(x[i,:], μ[j,:], Λ) for j=1:jmax, i=1:nt]
    T ./= sum(T; dims=1)
    τ = sum(T; dims=2)/nt
    μ = sum([T[j,i]*x[i,k] for j=1:jmax, k=1:2, i=1:nt]; dims=3)[:,:,1]./sum(T; dims=2)
    Λ = sum([T[j,i]*(x[i,:]-μ[j,:])*(x[i,:]-μ[j,:])' for i = 1:nt, j = 1:jmax])/nt
  end
  cs = [argmax(T[:,i]) for i = 1:nt]
    
  # plot clusters used for initial correction
  rc("font", size=8)
  rc("axes", titlesize="medium")
  fig, ax = subplots(1, 2, sharex=true, sharey=true, figsize=(190/25.4, 95/25.4))
  for i = 1:jmax
    ax[1].scatter(x[cs.==i,1], x[cs.==i,2], s=5)
  end
  xlim = round(maximum(abs.(ax[1].get_xlim())); digits=2)
  ylim = round(maximum(abs.(ax[1].get_ylim())); digits=2)
  ax[1].set_xlim(-xlim, xlim)
  ax[1].set_ylim(-ylim, ylim)
  ax[1].set_xlabel("low-frequency delay (s)")
  ax[1].set_ylabel("differential delay (s)")
  ax[1].set_title("after clustering")
  ax[1].set_title("(a)", loc="left")

  for i = 1:nt
    if cs[i] != 2
      s = @sprintf("%s %s\n", tpairs.event1[i], tpairs.event2[i])
      printstyled(s, color=(cs[i]==1) ? :blue : :green)
    end
  end

  # data vector using initial correction
  y = [reshape([Δτa[i,j,cs[i]] for i = 1:nt, j = 1:l], l*nt); ppairs.Δτ]
    
  # inverses
  iRnn = inv(Array(Rnn))
  #iRxx = inv(Array(Rxx))
  iRyy = inv(Array(E*Rxx*E' + Rnn))

  # solution operator (x̃ = Ay)
  A = P*E'*iRnn

#  # residual operator
#  R = I - E*A

  # initial residuals and cost
  x = A*y
  n = y - E*x
  r = n.^2
  J = .5*(y'*iRyy*y)#.5n'*iRnn*n + .5x'*iRxx*x

  # cycle through T-wave pairs until no further corrections are made
  i = 1
  idx = nothing
  while i ≤ nt

    @printf("\r%4d", i)

    # sort unique pairs by size of residuals
    if i == 1
      idx = sortperm(r[1:nt], rev=true)
    end
 
    # corrected this pair?
    corrected = false

    # trial correction(s)
    if cs[idx[i]] == 1
      if l2
        dir = [2, 5]
      else
        dir = [2]
      end
    elseif cs[idx[i]] == 3
      dir = [2, 4]
    elseif cs[idx[i]] == 4
      dir = [3]
    elseif cs[idx[i]] == 2
      dir = [1, 3]
    else
      dir = [1]
    end

    # go through trials
    for j = dir

      # check if probability of belonging to different cluster is at least 0.1%
      if T[j,idx[i]] ≥ 1e-3

        # swap out Δτ
        y[idx[i]] = Δτa[idx[i],1,j]

        # new residuals and cost
        x = A*y
        n = y - E*x
        rp = n.^2
        Jp = .5*(y'*iRyy*y)#.5n'*iRnn*n + .5x'*iRxx*x

        # record if cost is reduced, revert otherwise
        if Jp < J
          cs[idx[i]] = j
          J = Jp
          r = copy(rp)
          corrected = true
          @printf("\r%4d %4d %d %7.5f %s %s\n", i, idx[i], j, 1e3J/(nt+np+m-2)/l,
                  tpairs.event1[idx[i]], tpairs.event2[idx[i]])
        else
          y[idx[i]] = Δτa[idx[i],1,cs[idx[i]]]
        end

      end
    end

    # move back to first pair if correction was made, advance otherwise
    i = corrected ? 1 : i+1

  end

  @printf("\n")
  @printf("Total number of T-wave pairs:            %4d\n", nt)
  @printf("Number of pairs corrected to left max.:  %4d\n", sum(cs.==1))
  @printf("Number of uncorrected pairs:             %4d\n", sum(cs.==2))
  @printf("Number of pairs corrected to right max.: %4d\n", sum(cs.==3))
  @printf("Number of pairs corrected to second right max.: %4d\n", sum(cs.==4))
  if l2
    @printf("Number of pairs corrected to second left max.: %4d\n", sum(cs.==5))
  end
  
  # plot clusters used for initial correction
  Δτp = E[1:nt,1:m]*((A*y)[l*m+1:(l+1)*m])
  x1 = [tpairs.Δτc[i][1] - Δτp[i] for i = 1:nt]
  x2 = [tpairs.Δτc[i][1] - tpairs.Δτc[i][2] for i = 1:nt]
  x = [x1 x2]
  for i = 1:jmax
    ax[2].scatter(x[cs.==i,1], x[cs.==i,2], s=5)
  end
  ax[2].set_xlabel("low-frequency delay (s)")
  ax[2].set_title("after additional corrections")
  ax[2].set_title("(b)", loc="left")
  fig.tight_layout(w_pad=2)
  fig.savefig(@sprintf("results/plots/japan_%s_cluster.pdf",tstations[1]))
  
  x1outi = argmax(abs.(x1))
  ppairsexclude = ppairs[(ppairs.event1 .== tpairs.event1[x1outi]) .&& (ppairs.event2 .== tpairs.event2[x1outi]), :]
  ppairsexclude = ppairsexclude[:, [:station, :event1, :event2]]
  CSV.write(@sprintf("results/pairs/japan_%s_ppairsexclude_cs.csv",tstations[1]), ppairsexclude)
  
  @printf("\nx1 outlier: %s %s, x1 = %.3f s\n", tpairs.event1[x1outi],tpairs.event2[x1outi],x1[x1outi])
  
  x = A*y
  n = y - E*x
  Cn = (I(nt+np)-E*A)*(Rnn+E*Rxx*E')*(I(nt+np)-A'*E')
  iCnt = inv(Symmetric(Matrix(Cn[1:l*nt,1:l*nt])))
  iCnp = inv(Symmetric(Matrix(Cn[l*nt+1:end,l*nt+1:end])))
  U,S,Vt = svd(Symmetric(Rxx-P)[1+l*m:(l+1)*m,1+l*m:(l+1)*m])
  #iCo = inv(Matrix(Diagonal(diag(Rxx-P)[1+l*m:(l+1)*m])))
  @printf("iCnt posdef: %s, iCnp posdef: %s\n",isposdef(iCnt),isposdef(iCnp))
  zint = isposdef(iCnt) ? cholesky(iCnt).U*n[1:l*nt] : zeros(l*nt)
  zinp = isposdef(iCnp) ? cholesky(iCnp).U*n[1+l*nt:end] : zeros(np)
  zio = Diagonal(sqrt.(S.^-1))*U'*x[1+l*m:(l+1)*m]#isposdef(iCo) ? cholesky(iCo).U*x[1+l*m:(l+1)*m] : zeros(m)
  
  Δτ2 = [Δτa[i,2,cs[i]] for i = 1:nt]
  # return optimal delays
  return collect(eachrow(reshape(y[1:l*nt], (nt, l)))), Δτp, cs, x1,x2,Δτ2,zint,zinp,zio

end

"""
    c, ỹ = lineartrend(t, y; fitannual=false, fitsemiannual=false)

Performs least-squares regression ``y = c_1 t + c_2 + ε``. If `fitannual`,
``c_4 sin ωt + c_5 sin ωt`` is also included in the regression, where ``ω`` is the annual
frequency. If `fitsemiannual`, ``c_6 sin 2ωt + c_7 sin 2ωt`` is additionally included. The
actual numbering of the output coefficients ``c_i`` (`c[i]`) depends on which cycles are
included; the coefficients are in the above order. Returns the coefficient vector `c` and
the adjusted data `b`.
"""
function lineartrend(t, y; fitannual=false, fitsemiannual=false)

  # length of time vector
  m = length(t)

  # convert time to real numbers (in days)
  tr = float.(Dates.value.(t .- t[1]))/1000/3600/24

  # annual frequency
  ω = 2π/meanyear

  # design matrix for linear regression
  E = [tr ones(m)]

  # add annual cycle
  if fitannual
    E = hcat(E, sin.(ω*tr), cos.(ω*tr))
  end

  # add semiannual cycle
  if fitsemiannual
    E = hcat(E, sin.(2ω*tr), cos.(2ω*tr))
  end

  # invert
  x = E\y

  # return coefficients and adjusted data
  return x, E*x

end

# Argo/ECCO trend estimation
function estimatetrend(td, λ, σc, σm, σtrend, σannual, σsemiannual)
  # number of time points
  m = length(td)
  # number of frequencies
  l = length(σc)
  # annual frequency
  ω = 2π/SOT.meanyear
  # time mid-point
  tm = td[1]+(td[m]-td[1])/2
  # design matrix
  E = [kron(I(l), ones(m)) kron(I(l), td.-tm) kron(I(l), cos.(ω*td)) kron(I(l), sin.(ω*td)) kron(I(l), cos.(2ω*td)) kron(I(l), sin.(2ω*td))]
  # prior solution uncertainty
  Rxx = blockdiag(spdiagm(σm.^2), spdiagm(σtrend.^2), kron(I(2), spdiagm(σannual.^2/2)), kron(I(2), spdiagm(σsemiannual.^2/2)))
  # noise = anomalies correlated in time with correlation time λ
  C = exp.(-abs.(td.-td')/λ)
  Rnn = kron(spdiagm(σc.^2), C)
  # uncertainty matrix
  P = inv(inv(Array(Rxx)) + E'*inv(Array(Rnn))*E)
  # matrix to extract projected travel time anomalies associated with the trends
  M = zeros(size(E))
  M[:,l+1:2l] = E[:,l+1:2l]
  return E, Rxx, Rnn, P, M
end

""" interpolation onto regular grid of averages """
function regulargrid(td, ti, a, R, λ, h, U, Λ, σc)
  m = length(td)
  l = length(Λ)
  tavg = Dates.value(ti[2] - ti[1])/86400000
  T = h*U*Diagonal(Λ)
  tid = Dates.value.(ti .- DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24 # TODO why 12:00:00?
  tm = td[1] + (td[m] - td[1])/2
  ω = 2π/SOT.meanyear
  p = R\a
  corr(t1, t2) = abs(t1-t2)>tavg/2 ? 2λ/tavg*exp(-abs(t1-t2)/λ)*sinh(tavg/2λ) : λ/tavg*(2-exp((t1-t2-tavg/2)/λ)-exp(-(t1-t2+tavg/2)/λ))
  C = corr.(td', tid)
  F = [kron(T, I(length(ti)))*kron(spdiagm(σc.^2), C)*kron(T', I(m)) zeros(l*length(ti), m) zeros(l*length(ti), 5l)]
  τi = reshape(F*p, (length(ti), l)) + a[(l+1)*m+1:(l+1)*m+l]'.*(tid.-tm) + a[(l+1)*m+l+1:(l+1)*m+2l]'.*cos.(ω*tid) + a[(l+1)*m+2l+1:(l+1)*m+3l]'.*sin.(ω*tid) + a[(l+1)*m+3l+1:(l+1)*m+4l]'.*cos.(2ω*tid) + a[(l+1)*m+4l+1:(l+1)*m+5l]'.*sin.(2ω*tid)
  ci = Array((inv(T)*τi')')
  return τi, ci
end

""" interpolation onto regular grid of averages and return normalized z-scores """
function regulargrid_returnz(td, ti, a, R, λ, h, U, Λ, σc, P)
  m = length(td)
  l = length(Λ)
  tavg = Dates.value(ti[2] - ti[1])/86400000
  T = h*U*Diagonal(Λ)
  tid = Dates.value.(ti .- DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24 # TODO why 12:00:00?
  iR = inv(Symmetric(Matrix(R)))
  corr(t1, t2) = abs(t1-t2)>tavg/2 ? 2λ/tavg*exp(-abs(t1-t2)/λ)*sinh(tavg/2λ) : λ/tavg*(2-exp((t1-t2-tavg/2)/λ)-exp(-(t1-t2+tavg/2)/λ))
  C = corr.(td', tid)
  F = [kron(T, I(length(ti)))*kron(spdiagm(σc.^2), C)*kron(T', I(m)) zeros(l*length(ti), m) zeros(l*length(ti), 5l)]
  Ti = h^-1*kron(sparse(inv(Diagonal(Λ))*U'), I(length(ti)))
  ci = Ti*F*iR*a 
  Pi = Ti*F*iR*P*iR'*F'*Ti'
  Ci = kron(spdiagm(σc.^2), corr.(tid', tid))
  iCi = inv(Symmetric(Matrix(Ci-Pi)))
  @printf("is iCi posdef: %s\n",isposdef(iCi))
  zix = isposdef(iCi) ? cholesky(iCi).U*ci : zeros(length(ti)*l)
  return ci, zix
end

""" interpolation onto regular grid for 1 frequency """
function regulargrid_1f(td, θd, ti, θi, a, R, iRyy, E, λt, λθ, στ)
  @printf("\n")
  @printf("interpolation onto regular grid\n")
  m = length(td)
  l = 1

  T = I(1)
  tid = Dates.value.(ti .- DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24 # TODO why 12:00:00?
  tidv = vec(tid' .* ones(length(θi)))
  θiv = vec(ones(length(ti))' .* θi)
  tm = td[1] + (td[m] - td[1])/2
  ω = 2π/SOT.meanyear
  p = R\a
  corr(t1, θ1, t2, θ2) = exp(-abs(t1-t2)/λt-(θ1-θ2)^2/λθ^2/2)
  C = corr.(td', θd', tidv, θiv)
  @printf("covariance matrix size: %s\n",size(C))
  F = [kron(T, I(length(tidv)))*kron(spdiagm(στ.^2), C)*kron(T', I(m)) zeros(l*length(tidv), m) zeros(l*length(tidv), 5l)]
  F = [F; zeros(5l, (l+1)*m) R[end-4:end,end-4:end]]
  Di = [I(length(tidv)) (tidv.-tm) cos.(ω*tidv) sin.(ω*tidv) cos.(2ω*tidv) sin.(2ω*tidv)]
  @printf("calculate regular anomalies...\n")
  τi = reshape(Di*F*p, (length(θi),length(ti)))
  @printf("regular anomalies shape: %s\n",size(τi))
  Riy = Di*F*E'
  Rii = corr.(tidv', θiv', tidv, θiv)
  Rii = Di*[kron(spdiagm(στ.^2), Rii) zeros(length(tidv),5); zeros(5,length(tidv)) R[end-4:end,end-4:end]]*Di'
  ei = reshape(sqrt.(max.(0,diag(Rii-Riy*iRyy*Riy'))), (length(θi),length(ti))) 
  ei0 = reshape(sqrt.(diag(Rii)), (length(θi),length(ti))) 
  return τi, ei, ei0
end

""" interpolation onto regular grid with model covariance """
function regulargrid_gpcov(td, ti, a, y, R, E, N, invR, P, σc, ctau, lags)
  @printf("\n")
  @printf("interpolation onto regular grid\n")
  m = length(td)
  mi = length(ti)
  l = 3
  tid = Dates.value.(ti .- DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24 # TODO why 12:00:00?

  # solution covariance in time
  Δtrd = tid.-td'
  C = []
  for i = 1:6
    etpf = linear_interpolation(lags, ctau[:,i], extrapolation_bc=0)
    push!(C,etpf.(Δtrd))
    if i in [2,3,5]
      push!(C,etpf.(-Δtrd))
    end
  end
  C = Matrix([σc[1]^2*C[1] σc[2]^2*C[2] σc[3]^2*C[4]; σc[2]^2*C[3] σc[4]^2*C[6] σc[5]^2*C[7]; σc[3]^2*C[5] σc[5]^2*C[8] σc[6]^2*C[9]])
  A = I(l*mi)-kron(I(l), ones(mi,mi)./mi)

  F = [C zeros(l*mi, m) zeros(l*mi, l)]
  ΔC = F*invR*F'

  τi = F*E'*inv(Matrix(N+E*R*E'))*y 
  F = F*invR
  τi2 = F*a 
  Pi = F*P*F'
  iCi = inv(Symmetric(Matrix(ΔC-Pi)))
  
  # solution covariance in time
  Δtrd = tid.-tid'
  C = []
  for i = 1:6
    etpf = linear_interpolation(lags, ctau[:,i], extrapolation_bc=0)
    push!(C,etpf.(Δtrd))
  end
  C = Symmetric([σc[1]^2*C[1] σc[2]^2*C[2] σc[3]^2*C[3]; σc[2]^2*C[2]' σc[4]^2*C[4] σc[5]^2*C[5]; σc[3]^2*C[3]' σc[5]^2*C[5]' σc[6]^2*C[6]])

  @printf("is iCi posdef: %s\n",isposdef(iCi))
  zix = cholesky(inv(Symmetric(Matrix(C)))).U*τi2
  zix2 = cholesky(inv(Symmetric(Matrix(C)))).U*τi

  return reshape(A*τi, (mi, l)),A*Pi*A',C,zix,zix2
end
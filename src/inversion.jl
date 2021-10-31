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
function collectpairs(eqname, tstations, tintervals, tavgwidth, treffreq, tinvfreq, tmincc,
                      pstations, pintervals, pfreqbands; tmaxΔτ=Inf, excludetimes=[],
                      excludepairs=DataFrame(pstation=String[], event1=DateTime[],
                                             event2=DateTime[]))

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
  tpairs = DataFrame(station=String[], event1=DateTime[], event2=DateTime[],
                     Δτl=Array{Float64,1}[], Δτc=Array{Float64,1}[],
                     Δτr=Array{Float64,1}[], ccl=Array{Float64,1}[],
                     ccc=Array{Float64,1}[], ccr=Array{Float64,1}[])

  # collect all T-wave pairs that meet the criteria
  @showprogress for i = 1 : length(tstations), j = 1 : size(ppairs, 1)

    # file with T-wave data
    tdelayfile = @sprintf("%s/%s_%s.h5",
                          tdelaydir(eqname, tstations[i], tintervals[i], tavgwidth,
                                    treffreq),
                          fmttime(ppairs[j,:event1]), fmttime(ppairs[j,:event2]))

    # read data if present
    if isfile(tdelayfile)

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
        push!(tpairs, [tstations[i], ppairs[j,:event1], ppairs[j,:event2],
                       read(fid, "Δτl")[idx], Δτc, read(fid, "Δτr")[idx],
                       read(fid, "ccl")[idx], ccc, read(fid, "ccr")[idx]])

      end

      # close file
      close(fid)

    end

  end

  # list of all selected P-wave pairs
  ppairs = DataFrame(station=String[],
                     event1=DateTime[], latitude1=Float64[], longitude1=Float64[],
                     event2=DateTime[], latitude2=Float64[], longitude2=Float64[],
                     Δτ=Float64[], cc=Float64[])

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
function invert(tpairs, ppairs, λ, σx, σn, σp, U, Σ; σtrend=0, σannual=0, σsemiannual=0)

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
  A = exp.(-abs.(tr.-tr')/λ)

  # transformation from singular vectors to observed frequencies
  T = kron(sparse(U*Diagonal(Σ)), I(m))

  # covariance matrix assuming no correlation between singular vectors expansion coefficients
  Rxx = [T*kron(spdiagm(σx.^2), A)*T' zeros(l*m, m); zeros(m, (l+1)*m)] + σp^2*kron(sparse(ones(l+1, l+1)), I(m))

  # add trend if desired
  if !(σtrend == 0)
    tm = tr[1]+(tr[m]-tr[1])/2
    E = [E E*[kron(I(l), tr.-tm); zeros(m, 3)]]
    T = U*Diagonal(Σ)
    Rxx = [Rxx zeros(size(Rxx, 1), l); zeros(l, size(Rxx, 2)) T*spdiagm(σtrend.^2)*T']
    D = [D kron(I(l), tr.-tm)]
  end

  # add annual cycle if desired
  if !(σannual == 0)
    ω = 2π/365.2425
    E = [E E[:,1:l*m]*[kron(I(l), cos.(ω*tr)) kron(I(l), sin.(ω*tr))]]
    Rxx = [Rxx zeros(size(Rxx, 1), 2l); zeros(2l, size(Rxx, 2)) kron(I(2), T*spdiagm(σannual.^2/2)*T')]
    D = [D kron(I(l), cos.(ω*tr)) kron(I(l), sin.(ω*tr))]
  end

  # add semi-annual cycle if desired
  if !(σsemiannual == 0)
    ω = 2π/365.2425
    E = [E E[:,1:l*m]*[kron(I(l), cos.(2ω*tr)) kron(I(l), sin.(2ω*tr))]]
    Rxx = [Rxx zeros(size(Rxx, 1), 2l); zeros(2l, size(Rxx, 2)) kron(I(2), T*spdiagm(σannual.^2/2)*T')]
    D = [D kron(I(l), cos.(2ω*tr)) kron(I(l), sin.(2ω*tr))]
  end

  # noise covariance
  Rnn = σn^2*I

  # Gauss–Markov solution covariance matrix
  P = inv(inv(Array(Rxx)) + E'*inv(Rnn)*E)

  # Best estimate, subtracting T- and P-delays, adding trend and seasonal cycles:
  #   τ = D*P*E'*inv(Rnn)*y
  # Uncertainty:
  #   e = sqrt.(diag(D*P*D'))

  return t, E, Rxx, Rnn, P, D

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
  Σ = [0.2^2 0; 0 0.1^2]
  τ = [0.05, 0.9, 0.05]
  f(x, μ, Σ) = 1/sqrt((2π)^length(μ)*norm(Σ))*exp(-1/2*(x-μ)'*(Σ\(x-μ)))
  T = nothing
  for n = 1:100
    T = [τ[j].*f(x[i,:], μ[j,:], Σ) for j=1:3, i=1:nt]
    T ./= sum(T; dims=1)
    τ = sum(T; dims=2)/nt
    μ = sum([T[j,i]*x[i,k] for j=1:3, k=1:2, i=1:nt]; dims=3)[:,:,1]./sum(T; dims=2)
    Σ = sum([T[j,i]*(x[i,:]-μ[j,:])*(x[i,:]-μ[j,:])' for i = 1:nt, j = 1:3])/nt
  end
  cs = [argmax(T[:,i]) for i = 1:nt]

  # plot clusters used for initial correction
  fig, ax = subplots(1, 2, sharex=true, sharey=true, figsize=(190/25.4, 95/25.4))
  for i = 1:3
    ax[1].scatter(x[cs.==i,1], x[cs.==i,2], s=5)
  end
  ax[1].set_xlabel("low-frequency delay (s)")
  ax[1].set_ylabel("differential delay (s)")
  ax[1].set_title("after clustering")

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
  fig.tight_layout()

  # return optimal delays
  return collect(eachrow(reshape(y[1:l*nt], (nt, l))))

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

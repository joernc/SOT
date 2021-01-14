# TODO:
# - Use sparse QR decomposition if inversion matrix becomes too big (qr([E; S])\y)?

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
    selpairs.station = pstations[i]

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
"""
function invert(tpairs, ppairs; timescale=NaN)

  # number of frequencies
  l = length(tpairs.Δτc[1])

  # find unique events
  t = sort(unique([tpairs.event1; tpairs.event2]))

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

  # smoothing matrix
  tr = Dates.value.(t .- t[1])/1000/3600/24
  Δ = isnan(timescale) ? (tr[3:m] - tr[1:m-2])/2 : timescale*ones(m-2)
  F = spdiagm(m-2, m,
              0 => Δ.*(tr[2:m-1] - tr[1:m-2]).^-1,
              1 => -Δ.*((tr[2:m-1] - tr[1:m-2]).^-1 + (tr[3:m] - tr[2:m-1]).^-1),
              2 => Δ.*(tr[3:m] - tr[2:m-1]).^-1)
  F = [F[1,:]'; F; F[m-2,:]']

  # difference matrix to get τ from T- and P-wave anomalies
  D = [I(l*m) vcat([-I(m) for i = 1:l]...)]

  # full smoothing matrix
  S = blockdiag([F for i = 1:l]...)*D

  # calculate pseudoinverse
  P = pinv(Array(E'*E + S'*S))

  return t, E, S, P, D

end

"""
    correctcycleskipping(tpairs, ppairs, E, S, P)

Find cycle-skipping corrections. Applies corrections to adjacent maxima of the
cross-correlation function whenever they reduce the cost function. Returns the corrected
*T*-wave delays `Δτ`.
"""
function correctcycleskipping(tpairs, ppairs, E, S, P)
  
  # residual operator
  H = E*P*E'
  K = S*P*E'

  # number of frequencies
  l = length(tpairs.Δτc[1])

  # number of unique events
  m = size(E, 2)÷(l+1)

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
  x1 = [tpairs.Δτc[i][1] - Δτp[i,1] for i = 1:nt]
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

  # data vector using initial correction
  y = [reshape([Δτa[i,j,cs[i]] for i = 1:nt, j = 1:l], l*nt); ppairs.Δτ]

  # initial residuals and cost
  r = (y - H*y).^2
  J = .5sum(r) + .5sum((K*y).^2)

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
        rp = (y - H*y).^2
        Jp = .5sum(rp) + .5sum((K*y).^2)

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
    lineartrend(t, b; fitannual=false, fitsemiannual=false)

Performs least-squares regression ``b = c_1 t + c_2 + ε``. If `fitannual`,
``c_4 sin ωt + c_5 sin ωt`` is also included in the regression, where ``ω`` is the annual
frequency. If `fitsemiannual`, ``c_6 sin 2ωt + c_7 sin 2ωt`` is additionally included. The
actual numbering of the output coefficients ``c_i`` (`c[i]`) depends on which cycles are
included; the coefficients are in the above order. Returns the coefficient vector `c` and
the adjusted data `b`.
"""
function lineartrend(t, b; fitannual=false, fitsemiannual=false)

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
  c = E\b

  # return coefficients and adjusted data
  return c, E*c

end

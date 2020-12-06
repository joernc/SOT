# TODO:
# - Generate catalog of used P-wave pairs and stations

# reference time
global const tref = DateTime(2000, 1, 1, 0, 0, 0)

# mean year in days
global const meanyear = 365.2425

"""
    invert(eqname, tstations, pstations, invfreq, mincc; maxΔτ=Inf, excludetimes=[], timscale=NaN)

Invert measurements of the travel time changes between events ``Δτ`` for travel time
anomalies ``τ`` relative to an arbitrary but common reference.

# Arguments
- `eqname::String`: earthquake name to identify experiment.
- `tstations::String`: list of *T*-wave station designations.
- `pstations::String`: list of *P*-wave station designations.
- `invfreq::Array{Float64,1}`: frequencies at which to perform inversion.
- `mincc::Array{Float64,1}`: minimum CC requirement at the difference `invfreq`.
- `maxΔτ::Float64=Inf`: maximum allowable T-wave ``Δτ``; larger measurements are discarded.
- `excludetimes::Array{Array{Date,2},1}`: time periods to exclude from inversion.
- `timescale::Float64=NaN`: time scale at which to apply smoothing; default is to use time
  scale set by the sampling rather than a fixed scale.

# Examples
```
julia> t, τ, τerr, tpairs, ppairs = SOT.invert("nias", "H08S2..EDH", ["PSI"], [2, 4], [0.6, 0.4])
[...]

julia> t, τ, τerr, tpairs, ppairs = SOT.invert("nias", "H08S2..EDH", ["PSI"], [2, 4], [0.6, 0.4]; excludetimes=[[Date(2001, 1, 1) Date(2004, 12, 1)], [Date(2010, 1, 1) Date(2012, 1, 20)]])
[...]
```
"""
function invert(eqname, tstations, pstations, invfreq, mincc; maxΔτt=Inf, excludetimes=[],
               timescale=NaN)

  # number of frequencies at which to perform inversion
  l = length(invfreq)

  # load and combine catalogs of P-wave pairs
  catalog = Array{DataFrame,1}(undef, 0)
  for s in pstations
    filename = @sprintf("catalogs/%s_%s.csv", eqname, s)
    push!(catalog, DataFrame(CSV.File(filename, select=[1, 2], comment="#")))
  end
  catalog = sort(unique(vcat(catalog...)))

  # exclude events in specified time periods
  for i = 1:length(excludetimes)
    exclude1 = excludetimes[i][1] .< catalog[:event1] .< excludetimes[i][2]
    exclude2 = excludetimes[i][1] .< catalog[:event2] .< excludetimes[i][2]
    catalog = catalog[.!(exclude1 .| exclude2),:]
  end

  # list of all selected T-wave pairs
  tpairs = DataFrame(station=String[], event1=DateTime[], event2=DateTime[],
                     Δτl=Array{Float64,1}[], Δτc=Array{Float64,1}[],
                     Δτr=Array{Float64,1}[], ccl=Array{Float64,1}[],
                     ccc=Array{Float64,1}[], ccr=Array{Float64,1}[])

  # collect all T-wave pairs that meet the criteria
  @showprogress for pair in eachrow(catalog)

    # loop over T-wave stations
    for s in tstations

      # file with T-wave data
      filename = @sprintf("twavedelays/%s_%s/%s_%s.h5", eqname, s,
                          pair.event1, pair.event2)

      # read data if present
      if isfile(filename)

        # open file
        fid = h5open(filename, "r")

        # frequencies
        frequencies = read(fid, "freq")

        # find frequency indices
        idx = [argmin(abs.(frequencies .- invfreq[i])) for i = 1:l]

        # read central delay and CC
        Δτc = read(fid, "lagc")[idx]
        ccc = read(fid, "ccc")[idx]

        # check if criteria are met
        if all(ccc .> mincc) && all(abs.(Δτc) .< maxΔτt)

          # record measurements
          push!(tpairs, [s, pair.event1, pair.event2, read(fid, "lagl")[idx], Δτc,
                         read(fid, "lagr")[idx], read(fid, "ccl")[idx], ccc,
                         read(fid, "ccr")[idx]])

        end

        # close file
        close(fid)

      end

    end

  end

  # list of all selected P-wave pairs
  ppairs = DataFrame(station=String[], event1=DateTime[], event2=DateTime[], Δτ=Float64[],
                     cc=Float64[])

  # find all pairs that were selected
  for s in pstations

    # P-wave catalog
    filename = @sprintf("catalogs/%s_%s.csv", eqname, s)
    pairs = DataFrame(CSV.File(filename, comment="#"))

    # find pairs that were selected based on T-wave measurements
    selpairs = innerjoin(pairs, select(tpairs, [:event1, :event2]), on=[:event1, :event2])

    # add station information
    selpairs.station = s

    # record
    append!(ppairs, selpairs)

  end

  # perform inversion
  t, E, S, P = invmatrix(tpairs, ppairs; timescale)

  # number of good T- and P-wave pairs
  nt = size(tpairs, 1)
  np = size(ppairs, 1)

  # number of unique events
  m = length(t)

  @printf("Number of T-wave pairs:  %4d\n", nt)
  @printf("Number of P-wave pairs:  %4d\n", np)
  @printf("Number of unique events: %4d\n", m)

  # cycle skipping correction
  tpairs.Δτ = correctcycleskipping(tpairs, ppairs, E, S, P)

  # get travel time anomalies
  y = [vcat(tpairs.Δτ'...); repeat(ppairs.Δτ, 1, l)]
  x = P*(E'*y)
  τ = x[1:m,:] - x[m+1:2m,:]

  # calculate error
  if nt + np > m + 1
    H = diag(P)
    se = 1/(nt+np-m-1)*sum((y - E*x).^2, dims=1)
    xerr = sqrt.(se.*H)
    τerr = sqrt.(xerr[1:m,:].^2 + xerr[m+1:2m,:].^2)
  else
    τerr = Inf
  end

  # record inverted delays
  tpairs.Δτi = collect(eachrow(E[1:nt,:]*x))
  ppairs.Δτi = collect(eachrow(E[nt+1:nt+np,:]*x))

  # return inversion results, used pairs
  return t, τ, τerr, tpairs, ppairs

end

"""
    invmatrix(t1t, t2t, t1p, t2p; timescale=NaN)

Construct inversion matrix for T-wave pairs (`t1t`, `t2t`) and P-wave pairs (`t1p`, `t2p`).
Returns common times `t`, the design matrix `E`, the smoothing matrix `S`, and the
pseudoinverse `P = pinv(E'*E + S'*S)`. The inversion can then be performed with
`P*E'*[Δτt; Δτp]`. (See Wunsch: Discrete Inverse and State Estimation Problems, p. 56.) If
`timescale` is set (in days), smoothing is applied at that scale. By default
(`timescale=NaN`), smoothing is applied at the sampling scale.
"""
function invmatrix(tpairs, ppairs; timescale=NaN)

  # find unique events
  t = sort(unique([tpairs.event1; tpairs.event2]))
  tidx1 = indexin(tpairs.event1, t)
  tidx2 = indexin(tpairs.event2, t)
  pidx1 = indexin(ppairs.event1, t)
  pidx2 = indexin(ppairs.event2, t)

  # number of T- and P-wave pairs
  nt = size(tpairs, 1)
  np = size(ppairs, 1)

  # number of unique events
  m = length(t)

  # time in days
  tr = float.(Dates.value.(t .- tref))/1000/3600/24

  # T-wave pair matrix
  Xt = zeros(nt, m)
  for i = 1:nt
    Xt[i,tidx1[i]] = -1
    Xt[i,tidx2[i]] = 1
  end

  # P-wave pair matrix
  Xp = zeros(np, m)
  for i = 1:np
    Xp[i,pidx1[i]] = -1
    Xp[i,pidx2[i]] = 1
  end

  # smoothing matrix
  S = zeros(m-2, m)
  for i = 1:m-2
    Δ = isnan(timescale) ? (tr[i+2] - tr[i])/2 : timescale
    S[i,i] = Δ/(tr[i+1] - tr[i])
    S[i,i+1] = -Δ*(1/(tr[i+1] - tr[i]) + 1/(tr[i+2] - tr[i+1]))
    S[i,i+2] = Δ/(tr[i+2] - tr[i+1])
  end

  # full design matrix
  E = [Xt zeros(nt, m); zeros(np, m) Xp]

  # full smoothing matrix
  S = [S/2 -S/2]

  # calculate pseudoinverse
  P = pinv(E'*E + S'*S)

  return t, E, S, P

end

"""
    correctcycleskipping(Δτtl, Δτtc, Δτtr, ccl, ccc, ccr, Δτp, A, Ap)

Find cycle-skipping corrections. Applies corrections to adjacent maxima of the
cross-correlation function whenever they reduce the cost function. Returns the corrected
T-wave delays ``Δτt``.
"""
function correctcycleskipping(tpairs, ppairs, E, S, P)
  
  # residual operator
  H = E*P*E'
  K = S*P*E'

  # number of unique events
  m = size(E, 2)÷2

  # number of T- and P-wave pairs
  nt = size(tpairs, 1)
  np = size(ppairs, 1)

  # number of frequencies
  l = length(tpairs.Δτc[1])

  # collect the three delays and CCs
  Δτa = cat(vcat(tpairs.Δτl'...), vcat(tpairs.Δτc'...), vcat(tpairs.Δτr'...), dims=3)
  cca = cat(vcat(tpairs.ccl'...), vcat(tpairs.ccc'...), vcat(tpairs.ccr'...), dims=3)

  # data vector using central measurements
  y = [vcat(tpairs.Δτc'...); repeat(ppairs.Δτ, 1, l)]

  # initial residuals and cost
  r = sum((y - H*y).^2, dims=2)[:,1]
  J = .5sum(r) + .5sum((K*y).^2)

  # index used (start with central)
  cs = 2ones(Int, nt)

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

      # check if neighboring CCs are close to max CCs
      if all(cca[idx[i],:,j] .≥ cca[idx[i],:,2] .- 0.15)

        # swap out Δτ
        y[idx[i],:] = Δτa[idx[i],:,j]

        # new residuals and cost
        rp = sum((y - H*y).^2, dims=2)[:,1]
        Jp = .5sum(rp) + .5sum((K*y).^2)

        # record if cost is reduced, revert otherwise
        if Jp < J
          cs[idx[i]] = j
          J = Jp
          r = copy(rp)
          corrected = true
          @printf("\r%4d %4d %d %7.5f\n", i, idx[i], j, 1e3J/(nt+np+m-2)/l)
        else
          y[idx[i],:] = Δτa[idx[i],:,cs[idx[i]]]
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

  # return optimal delays
  return collect(eachrow(y[1:nt,:]))

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
  tr = float.(Dates.value.(t .- tref))/1000/3600/24

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

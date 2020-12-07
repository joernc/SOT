# TODO:
# - Use sparse QR decomposition if inversion matrix becomes too big (qr([E; S])\y)

# reference time
global const tref = DateTime(2000, 1, 1, 0, 0, 0)

# mean year in days
global const meanyear = 365.2425

"""
    collectpairs(eqname, tstations, pstations, invfreq, mincc; maxΔτ=Inf, excludetimes=[])

Collect all *T*- and *P*-wave pairs from the specified `tstations` and `pstations` that meet
certain criteria and record the measurements at the frequencies specified by `invfreq`. The
*T*-wave measurements are required to have a minimum cross-correlation coefficient of
`mincc`. *T*-wave outliers can be excluded with `maxΔτ`, and time periods can be excluded
with `excludetimes`.

# Examples
```
julia> tpairs, ppairs = SOT.collectpairs("nias", ["H01W1..EDH", "H01W3..EDH"], ["PSI", "KUM", "WRAB"], [2.5, 4], [0.6, 0.4], maxΔτ=20)
[...]

julia> tpairs, ppairs = SOT.collectpairs("nias", ["H08S2..EDH"], ["PSI"], [2, 4], [0.6, 0.4]; excludetimes=[[Date(2001, 1, 1) Date(2004, 12, 1)], [Date(2010, 1, 1) Date(2012, 1, 20)]])
[...]
```
"""
function collectpairs(eqname, tstations, pstations, invfreq, mincc;
                      maxΔτ=Inf, excludetimes=[])

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
        if all(ccc .> mincc) && all(abs.(Δτc) .< maxΔτ)

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
# number of unique events
m = length(t)
# data vector
y = [vcat(tpairs.Δτ'...); repeat(ppairs.Δτ, 1, l)]
# least-square solution
x = P*(E'*y)
# difference between T- and P-wave anomalies
τ = D*x
```
The standard error of the solution can then be obtained using
```
# number of T- and P-wave pairs
nt = size(tpairs, 1)
np = size(ppairs, 1)
# estimate variance
σ2 = 1/(nt+np-m-1)*sum((y - E*x).^2, dims=1)
# propagate error
A = D*P*E'
τerr = sqrt.(σ2.*diag(A*A'))
```
"""
function invert(tpairs, ppairs; timescale=NaN)

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
  Xt = sparse([1:nt; 1:nt], [tidx1; tidx2], [-ones(nt); ones(nt)])

  # P-wave pair matrix
  Xp = sparse([1:np; 1:np], [pidx1; pidx2], [-ones(np); ones(np)])

  # smoothing matrix
  Δ = isnan(timescale) ? (tr[3:m] - tr[1:m-2])/2 : timescale*ones(m-2)
  S = spdiagm(m-2, m, 0 => Δ ./ (tr[2:m-1] - tr[1:m-2]),
              1 => -Δ.*(1 ./ (tr[2:m-1] - tr[1:m-2]) + 1 ./ (tr[3:m] - tr[2:m-1])),
              2 => Δ ./ (tr[3:m] - tr[2:m-1]))

  # full design matrix
  E = [Xt spzeros(nt, m); spzeros(np, m) Xp]

  # full smoothing matrix
  S = [S/2 -S/2]

  # calculate pseudoinverse
  P = pinv(Array(E'*E + S'*S))

  # difference matrix to get τ from T- and P-wave anomalies
  D = [spdiagm(0 => ones(m)) spdiagm(0 => -ones(m))]

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

# reference time
global const tref = DateTime(2000, 1, 1, 0, 0, 0)

# mean year in days
global const meanyear = 365.2425

"""
    invert(eqname, pstations, tstation, invfreq, mincc; maxΔτ=Inf, excludetimes=[], timscale=NaN)

Invert measurements of the travel time changes between events ``Δτ`` for travel time
anomalies ``τ`` relative to an arbitrary but common reference.

# Arguments
- `eqname::String`: earthquake name to identify experiment.
- `pstations::String`: list of *P*-wave station designations.
- `tstation::String`: *T*-wave station designation.
- `invfreq::Array{Float64,1}`: frequencies at which to perform inversion.
- `mincc::Array{Float64,1}`: minimum CC requirement at the difference `invfreq`.
- `maxΔτ::Float64=Inf`: maximum allowable ``Δτ``; larger measurements are discarded.
- `excludetimes::Array{Array{Date,2},1}`: time periods to exclude from inversion.
- `timescale::Float64=NaN`: time scale at which to apply smoothing; default is to use time
  scale set by the sampling rather than a fixed scale.

# Examples
```
julia> t, τ, τerr = SOT.invert("nias", catalog, ["PSI"], "H08S2..EDH", [2, 4], [0.6, 0.4]; Δτmax=0.65)
[...]

julia> t, τ, τerr = SOT.invert("nias", catalog, ["PSI"], "H08S2..EDH", [2, 4], [0.6, 0.4]; Δτmax=0.65, excludetimes=[[Date(2001, 1, 1) Date(2004, 12, 1)], [Date(2010, 1, 1) Date(2012, 1, 20)]])
[...]
```
"""
function invert(eqname, pstations, tstation, invfreq, mincc; maxΔτ=Inf, excludetimes=[],
               timescale=NaN)

  # number of frequencies at which to perform inversion
  l = length(invfreq)

  # initialize event times
  t1 = Array{DateTime}(undef, 0)
  t2 = Array{DateTime}(undef, 0)

  # initialize delay measurements
  Δτc = Array{Array{Float64,1}}(undef, 0)
  Δτr = Array{Array{Float64,1}}(undef, 0)
  Δτl = Array{Array{Float64,1}}(undef, 0)

  # initialize CC measurements
  ccc = Array{Array{Float64,1}}(undef, 0)
  ccr = Array{Array{Float64,1}}(undef, 0)
  ccl = Array{Array{Float64,1}}(undef, 0)

  # iterate over stations
  for s in pstations

    # load catalog of P-wave pairs
    catalog = DataFrame(CSV.File("catalogs/$(eqname)_$s.csv"))

    # exclude events in specified time periods
    for i = 1:length(excludetimes)
      exclude1 = excludetimes[i][1] .< catalog[:event1] .< excludetimes[i][2]
      exclude2 = excludetimes[i][1] .< catalog[:event2] .< excludetimes[i][2]
      catalog = catalog[.!(exclude1 .| exclude2),:]
    end

    # loop over all pairs
    @showprogress for (i, pair) in enumerate(eachrow(catalog))

      # file with T-wave data
      filename = @sprintf("twavedelays/%s_%s_%s/%s_%s.h5", eqname, s, tstation,
                          pair.event1, pair.event2)

      # read data if present
      if isfile(filename)

        # frequencies
        frequencies = h5read(filename, "freq")
        # find frequency indices
        idx = [argmin(abs.(frequencies.-invfreq[i])) for i = 1:l]
        # read data
        push!(t1, pair.event1)
        push!(t2, pair.event2)
        push!(Δτc, h5read(filename, "lagc")[idx])
        push!(Δτr, h5read(filename, "lagr")[idx])
        push!(Δτl, h5read(filename, "lagl")[idx])
        push!(ccc, h5read(filename, "ccc")[idx])
        push!(ccr, h5read(filename, "ccr")[idx])
        push!(ccl, h5read(filename, "ccl")[idx])
      end
    end

  end

  # concatenate into arrays
  Δτc = hcat(Δτc...)'
  Δτr = hcat(Δτr...)'
  Δτl = hcat(Δτl...)'
  ccc = hcat(ccc...)'
  ccr = hcat(ccr...)'
  ccl = hcat(ccl...)'

  # select good pairs
  idx = all(ccc .≥ mincc'; dims=2)[:,1] .& (abs.(Δτc[:,1]) .≤ maxΔτ)

  # number of good pairs
  n = sum(idx)

  # perform inversion
  t, X, A, Ap = invmatrix(t1[idx], t2[idx], timescale=timescale)

  # number of unique events
  m = length(t)

  @printf("Number of used pairs:    %4d\n", n)
  @printf("Number of unique events: %4d\n", m)

  # cycle skipping correction
  Δτ = correctcycleskipping(t1[idx], t2[idx], Δτc[idx,:], Δτr[idx,:], Δτl[idx,:],
                                ccc[idx,:], ccr[idx,:], ccl[idx,:], X, A, Ap)

  # get travel times
  τ = Ap*[Δτ; zeros(m-2, l)]

  # calculate error
  if n > m + 1
    H = diag(pinv(A'*A))
    τerr = sqrt.(1/(n-m-1)*sum((Δτ - X*τ).^2, dims=1).*H)
  else
    τerr = Inf
  end

  return t, τ, τerr

end

"""
Construct inversion matrix for pairs (`t1`, `t2`). Returns common times `t`, the pair
matrix `X`, the full design matrix `A`, and the pseudoinverse `Ap`. The inversion can then
be performed with `Ap*Δτ`. If `timescale` is set (in days), smoothing is applied at that
scale. By default (`timescale=NaN`), smoothing is applied at the sampling scale.
"""
function invmatrix(t1, t2; timescale=NaN)

  # find unique events
  t = sort(unique([t1; t2]))
  idx1 = indexin(t1, t)
  idx2 = indexin(t2, t)

  # number of pairs
  n = length(t1)

  # number of unique events
  m = length(t)

  # time in days
  tr = float.(Dates.value.(t .- tref))/1000/3600/24

  # pair matrix
  X = zeros(n, m)
  for i = 1:n
    X[i,idx1[i]] = -1
    X[i,idx2[i]] = 1
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
  A = [X; S]

  # calculate pseudoinverse
  Ap = pinv(A)

  return t, X, A, Ap

end

"""
Find cycle skipping corrections. Applies corrections to adjacent maxima of the
cross-correlation function whenever they reduce the loss function. Returns the corrected
``Δτ``.
"""
function correctcycleskipping(t1, t2, Δτc, Δτr, Δτl, ccc, ccr, ccl, X, A, Ap)
  
  # residual operator
  R = I - A*Ap

  # number of pairs and unique events
  n, m = size(X)

  # number of frequencies
  l = size(Δτc)[2]

  # collect the three delays and CCs
  Δτa = cat(Δτl, Δτc, Δτr, dims=3)
  cca = cat(ccl, ccc, ccr, dims=3)

  # data vector using central measurements
  b = [Δτc; zeros(m-2, l)]

  # initial residuals and cost
  r = R*b
  J = .5*sum(r.^2)

  # find unique pairs
  pairs = unique([findall((t1 .== t1[i]) .& (t2 .== t2[i])) for i = 1:n])

  # number of unique pairs
  nu = length(pairs)

  # index used (start with central)
  cs = 2ones(Int, nu)

  # cycle through unique pairs until no further corrections are made
  i = 1
  while i ≤ nu

    @printf("\r%4d", i)

    # sort unique pairs by size of residuals
    if i == 1
      res = [mean(r[pairs[i],:].^2) for i = 1:nu]
      idx = sortperm(res, rev=true)
      pairs = pairs[idx]
      cs = cs[idx]
    end

    # corrected this pair?
    corrected = false

    # trial correction(s)
    if cs[i] == 1 || cs[i] == 3
      dir = [2]
    else
      dir = [1, 3]
    end

    # go through trials
    for j = dir

      # check if neighboring CCs are close to max CCs
      if all(cca[pairs[i],:,j] .≥ ccc[pairs[i],:] .- .15)

        # swap out Δτ
        b[pairs[i],:] = Δτa[pairs[i],:,j]

        # new residuals and cost
        rp = R*b
        Jp = .5*sum(rp.^2)

        # record if cost is reduced, revert otherwise
        if Jp < J
          cs[i] = j
          J = Jp
          r = copy(rp)
          corrected = true
          @printf("\r%4d %d %7.5f %s\n", i, j, 1e3J/(n+m-2)/l, string(pairs[i]))
        else
          b[pairs[i],:] = Δτa[pairs[i],:,cs[i]]
        end

      end
    end

    # move back to first pair if correction was made, advance otherwise
    i = corrected ? 1 : i+1

  end

  @printf("\n")
  @printf("Total number of unique pairs:            %4d\n", nu)
  @printf("Number of pairs corrected to left max.:  %4d\n", sum(cs.==1))
  @printf("Number of uncorrected pairs:             %4d\n", sum(cs.==2))
  @printf("Number of pairs corrected to right max.: %4d\n", sum(cs.==3))

  # return optimal delays
  return b[1:n,:]

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
  A = [tr ones(m)]

  # add annual cycle
  if fitannual
    A = hcat(A, sin.(ω*tr), cos.(ω*tr))
  end

  # add semiannual cycle
  if fitsemiannual
    A = hcat(A, sin.(2ω*tr), cos.(2ω*tr))
  end

  # invert
  c = A\b

  # return coefficients and adjusted data
  return c, A*c

end

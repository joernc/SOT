# TODO:
# - Generate catalog of used P-wave pairs and stations

# reference time
global const tref = DateTime(2000, 1, 1, 0, 0, 0)

# mean year in days
global const meanyear = 365.2425

"""
    invert(eqname, tstation, pstations, invfreq, mincc; maxΔτ=Inf, excludetimes=[], timscale=NaN)

Invert measurements of the travel time changes between events ``Δτ`` for travel time
anomalies ``τ`` relative to an arbitrary but common reference.

# Arguments
- `eqname::String`: earthquake name to identify experiment.
- `tstation::String`: *T*-wave station designation.
- `pstations::String`: list of *P*-wave station designations.
- `invfreq::Array{Float64,1}`: frequencies at which to perform inversion.
- `mincc::Array{Float64,1}`: minimum CC requirement at the difference `invfreq`.
- `maxΔτ::Float64=Inf`: maximum allowable T-wave ``Δτ``; larger measurements are discarded.
- `excludetimes::Array{Array{Date,2},1}`: time periods to exclude from inversion.
- `timescale::Float64=NaN`: time scale at which to apply smoothing; default is to use time
  scale set by the sampling rather than a fixed scale.

# Examples
```
julia> t, τ, τerr = SOT.invert("nias", "H08S2..EDH", ["PSI"], [2, 4], [0.6, 0.4])
[...]

julia> t, τ, τerr = SOT.invert("nias", "H08S2..EDH", ["PSI"], [2, 4], [0.6, 0.4]; excludetimes=[[Date(2001, 1, 1) Date(2004, 12, 1)], [Date(2010, 1, 1) Date(2012, 1, 20)]])
[...]
```
"""
function invert(eqname, tstation, pstations, invfreq, mincc; maxΔτt=Inf, excludetimes=[],
               timescale=NaN)

  # number of frequencies at which to perform inversion
  l = length(invfreq)

  # initialize event times
  t1t = Array{DateTime}(undef, 0)
  t2t = Array{DateTime}(undef, 0)

  # initialize T-wave delay measurements
  Δτtl = Array{Array{Float64,1}}(undef, 0)
  Δτtr = Array{Array{Float64,1}}(undef, 0)
  Δτtc = Array{Array{Float64,1}}(undef, 0)

  # initialize T-wave CC measurements
  ccc = Array{Array{Float64,1}}(undef, 0)
  ccr = Array{Array{Float64,1}}(undef, 0)
  ccl = Array{Array{Float64,1}}(undef, 0)

  # load and combine catalogs of P-wave pairs
  catalog = Array{DataFrame,1}(undef, 0)
  for s in pstations
    push!(catalog, DataFrame(CSV.File("catalogs/$(eqname)_$s.csv", select=[1, 2],
                                      comment="#")))
  end
  catalog = sort(unique(vcat(catalog...)))

  # exclude events in specified time periods
  for i = 1:length(excludetimes)
    exclude1 = excludetimes[i][1] .< catalog[:event1] .< excludetimes[i][2]
    exclude2 = excludetimes[i][1] .< catalog[:event2] .< excludetimes[i][2]
    catalog = catalog[.!(exclude1 .| exclude2),:]
  end

  exclude = Array{Int}(undef, 0)

  # loop over all pairs
  @showprogress for (i, pair) in enumerate(eachrow(catalog))

    # file with T-wave data
    filename = @sprintf("twavedelays/%s_%s/%s_%s.h5", eqname, tstation,
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
      Δτt = read(fid, "lagc")[idx]
      cc = read(fid, "ccc")[idx]

      # check if criteria are met
      if all(cc .> mincc) && all(abs.(Δτt) .< maxΔτt)

        # record event times
        push!(t1t, pair.event1)
        push!(t2t, pair.event2)

        # read data
        push!(Δτtl, read(fid, "lagl")[idx])
        push!(Δτtc, Δτt)
        push!(Δτtr, read(fid, "lagr")[idx])
        push!(ccl, read(fid, "ccl")[idx])
        push!(ccc, cc)
        push!(ccr, read(fid, "ccr")[idx])

      else

        # record excluded event
        push!(exclude, i)

      end

      # close file
      close(fid)

    else

      # record missing file
      push!(exclude, i)

    end
  end

  # concatenate into arrays
  Δτtc = hcat(Δτtc...)'
  Δτtr = hcat(Δτtr...)'
  Δτtl = hcat(Δτtl...)'
  ccc = hcat(ccc...)'
  ccr = hcat(ccr...)'
  ccl = hcat(ccl...)'

  # remove pairs from catalog that were excluded
  delete!(catalog, sort(unique(exclude)))

  # initialize event times
  t1p = Array{DateTime}(undef, 0)
  t2p = Array{DateTime}(undef, 0)

  # initialize P-wave delay measurements
  Δτp = Array{Float64}(undef, 0)

  # iterate over P-wave stations
  for s in pstations

    # load catalog of P-wave pairs
    pcatalog = DataFrame(CSV.File("catalogs/$(eqname)_$s.csv", comment="#"))

    # iterate over all pairs in P-wave catalog
    for (i, pair) in enumerate(eachrow(pcatalog))

      # record if T-wave measurement is present
      if any((pair.event1 .== catalog.event1) .& (pair.event2 .== catalog.event2))

        # record event times
        push!(t1p, pair.event1)
        push!(t2p, pair.event2)

        # record measurement
        push!(Δτp, pair.origincorrection)

      end

    end

  end

#  # select good pairs
#  idx = all(ccc .≥ mincc'; dims=2)[:,1] .& (abs.(Δτc[:,1]) .≤ maxΔτ)

  # perform inversion
  t, E, S, P = invmatrix(t1t, t2t, t1p, t2p; timescale)

  # number of good T-wave pairs
  nt = length(t1t)

  # number of good P-wave pairs
  np = length(t1p)

  # number of unique events
  m = length(t)

  @printf("Number of T-wave pairs:  %4d\n", nt)
  @printf("Number of P-wave pairs:  %4d\n", np)
  @printf("Number of unique events: %4d\n", m)

  # cycle skipping correction
  Δτt = correctcycleskipping(Δτtl, Δτtc, Δτtr, ccl, ccc, ccr, Δτp, E, P)

  # get travel times
  τa = P*E'*[Δτt; repeat(Δτp, 1, l)]
  τ = τa[1:m,:] - τa[m+1:2m,:]

  # calculate error
  if nt + np > m + 1
    H = diag(P)
    se = 1/(nt+np-m-1)*sum(([Δτt; repeat(Δτp, 1, l)] - E*τa).^2, dims=1)
    τaerr = sqrt.(se.*H)
    τerr = sqrt.(τaerr[1:m,:].^2 + τaerr[m+1:2m,:].^2)
  else
    τerr = Inf
  end

  # inverted delays
  Δτti = E[1:nt,:]*τa
  Δτpi = E[nt+1:nt+np,:]*τa

  idx = sortperm(abs.(Δτt[:,1] - Δτti[:,1]), rev=true)
  println("Largest five T-wave residuals:")
  println(catalog[idx,:][1:5,:])

  # plot measured vs. inverted T-wave delays
  fig, ax = subplots(1, 1)
  ax.scatter(Δτt[:,1], Δτti[:,1], s=5)
  ax.set_aspect(1)
  xlim = ax.get_xlim()
  ylim = ax.get_ylim()
  x = [-2maxΔτt, 2maxΔτt]
  ax.plot(x, x, color="black", linewidth=.8)
  ax.plot(x, x .+ 1/invfreq[1], color="black", linewidth=.8, zorder=0)
  ax.plot(x, x .- 1/invfreq[1], color="black", linewidth=.8, zorder=0)
  ax.plot(x, x .+ 1/2invfreq[1], color="black", linewidth=.8, zorder=0)
  ax.plot(x, x .- 1/2invfreq[1], color="black", linewidth=.8, zorder=0)
  ax.set_xlim(xlim)
  ax.set_ylim(ylim)
  ax.set_title("T waves")
  ax.set_xlabel("measured delay (s)")
  ax.set_ylabel("inverted delay (s)")

  # plot measured vs. inverted P-wave delays
  fig, ax = subplots(1, 1)
  scatter(Δτp[:,1], Δτpi[:,1], s=5)
  ax.set_aspect(1)
  xlim = ax.get_xlim()
  ylim = ax.get_ylim()
  x = [-2maxΔτt, 2maxΔτt]
  ax.plot(x, x, color="black", linewidth=.8)
  ax.plot(x, x .+ 1/invfreq[1], color="black", linewidth=.8, zorder=0)
  ax.plot(x, x .- 1/invfreq[1], color="black", linewidth=.8, zorder=0)
  ax.plot(x, x .+ 1/2invfreq[1], color="black", linewidth=.8, zorder=0)
  ax.plot(x, x .- 1/2invfreq[1], color="black", linewidth=.8, zorder=0)
  ax.set_xlim(xlim)
  ax.set_ylim(ylim)
  ax.set_title("P waves")
  ax.set_xlabel("measured delay (s)")
  ax.set_ylabel("inverted delay (s)")

  return t, τ, τerr

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
function invmatrix(t1t, t2t, t1p, t2p; timescale=NaN)

  # find unique events
  t = sort(unique([t1t; t2t]))
  idx1t = indexin(t1t, t)
  idx2t = indexin(t2t, t)
  idx1p = indexin(t1p, t)
  idx2p = indexin(t2p, t)

  # number of T- and P-wave pairs
  nt = length(t1t)
  np = length(t1p)

  # number of unique events
  m = length(t)

  # time in days
  tr = float.(Dates.value.(t .- tref))/1000/3600/24

  # T-wave pair matrix
  Xt = zeros(nt, m)
  for i = 1:nt
    Xt[i,idx1t[i]] = -1
    Xt[i,idx2t[i]] = 1
  end

  # P-wave pair matrix
  Xp = zeros(np, m)
  for i = 1:np
    Xp[i,idx1p[i]] = -1
    Xp[i,idx2p[i]] = 1
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
function correctcycleskipping(Δτtl, Δτtc, Δτtr, ccl, ccc, ccr, Δτp, E, P)
  
  # residual operator
  R = I - E*P*E'

  # number of unique events
  m = size(E, 2)÷2

  # number of T-wave pairs and frequencies
  nt, l = size(Δτtc)

  # number of P-wave pairs
  np = length(Δτp)

  # collect the three delays and CCs
  Δτa = cat(Δτtl, Δτtc, Δτtr, dims=3)
  cca = cat(ccl, ccc, ccr, dims=3)

  # data vector using central measurements
  b = [Δτtc; repeat(Δτp, 1, l)]

  # initial residuals and cost
  r = R*b
  J = .5*sum(r.^2)

  # index used (start with central)
  cs = 2ones(Int, nt)

  # cycle through T-wave pairs until no further corrections are made
  i = 1
  idx = nothing
  while i ≤ nt

    @printf("\r%4d", i)

    # sort unique pairs by size of residuals
    if i == 1
      idx = sortperm(r[1:nt].^2, rev=true)
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
      if all(cca[idx[i],:,j] .≥ ccc[idx[i],:] .- 0.15)

        # swap out Δτ
        b[idx[i],:] = Δτa[idx[i],:,j]

        # new residuals and cost
        rp = R*b
        Jp = .5*sum(rp.^2)

        # record if cost is reduced, revert otherwise
        if Jp < J
          cs[idx[i]] = j
          J = Jp
          r = copy(rp)
          corrected = true
          @printf("\r%4d %4d %d %7.5f\n", i, idx[i], j, 1e3J/(nt+np+m-2)/l)
        else
          b[idx[i],:] = Δτa[idx[i],:,cs[idx[i]]]
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
  return b[1:nt,:]

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

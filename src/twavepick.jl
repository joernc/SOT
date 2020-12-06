# TODO:
# - Allow for clock error corrections (before cutting).

# Earth's radius
global const earthradius = 6371e3

# sound speed in water (for arrival time estimation)
global const soundspeed = 1.51e3

"""
    twavepick(eqname, tstation, pstations, starttime, endtime, frequencies, avgwidth, reffreq; saveplot=false)

Picks the lags at which the cross-correlation function between pairs is maximized, saving
the results to file.

# Arguments
- `eqname::String`: earthquake name to identify experiment.
- `tstation::String`: *T*-wave station designation.
- `pstations::String`: list of *P*-wave station designations.
- `starttime::Float64`: beginning of time window in number of seconds before predicted
  arrival time.
- `endtime::Float64`: end of time window in number of seconds after predicted arrival time.
- `frequencies::Array{Float64,1}`: frequencies at which to measure.
- `avgwidth::Float64`: width of Gaussian frequency window.
- `reffreq::Float64`: reference frequency at which to find max cross correlation.
- `saveplot::Boolean`: whether to save a plot.

# Example
```
julia> twavepick("nias", "H08S2..EDH", ["PSI", "KUM", "WRAB"], 10, 70, 0.1:0.1:10, 0.5, 2.0)
[...]
```
"""
function twavepick(eqname, tstation, pstations, starttime, endtime, frequencies, avgwidth,
                   reffreq; saveplot=false)

  # load and combine catalogs of P-wave pairs
  catalog = Array{DataFrame,1}(undef, 0)
  for s in pstations
    push!(catalog, DataFrame(CSV.File("catalogs/$(eqname)_$s.csv", select=[1, 2],
                                      comment="#")))
  end
  catalog = sort(unique(vcat(catalog...)))

  # iterate over pairs
  for (i, pair) in enumerate(eachrow(catalog))

    @printf("%d/%d: %s %s\n", i, size(catalog)[1], pair.event1, pair.event2)

    # paths to T waveforms
    path1 = @sprintf("twaveforms/%s/%s", Dates.format(pair.event1, "yyyymmddHHMMSS"),
                     tstation)
    path2 = @sprintf("twaveforms/%s/%s", Dates.format(pair.event2, "yyyymmddHHMMSS"),
                     tstation)

    # file to which measurements are saved
    savepath = @sprintf("twavedelays/%s_%s", eqname, tstation)
    !isdir(savepath) && mkdir(savepath)
    filename = @sprintf("%s/%s_%s.h5", savepath, pair.event1, pair.event2)
    figpath = @sprintf("twaveplots/%s_%s", eqname, tstation)
    !isdir(figpath) && mkdir(figpath)

    # check whether all data is present
    if !isfile(filename) && isfile(path1) && isfile(path2)

      # load both waveforms
      trace1 = read_sac(path1)
      trace2 = read_sac(path2)

      # estimate arrival time from geodetic distance (relative to event time)
      evtlon = .5(trace1.evt.lon + trace2.evt.lon)
      evtlat = .5(trace1.evt.lat + trace2.evt.lat)
      stalon = .5(trace1.sta.lon + trace2.sta.lon)
      stalat = .5(trace1.sta.lat + trace2.sta.lat)
      range = dist(evtlon, evtlat, stalon, stalat)
      arrivaltime = range/soundspeed

      # sampling frequency
      fs = 1/Float64(trace1.delta)

      # sampling times
      time1 = trace1.b .+ (0:length(trace1.t)-1)/fs
      time2 = trace2.b .+ (0:length(trace2.t)-1)/fs

      # define window
      idx1 = arrivaltime - starttime .< time1 .< arrivaltime + endtime
      idx2 = arrivaltime - starttime .< time2 .< arrivaltime + endtime

      # check whether data is present and has same length in both traces
      if any(idx1) && sum(idx1) == sum(idx2)

        # cut signals
        s1 = Float64.(trace1.t[idx1])
        s2 = Float64.(trace2.t[idx2])

        # sample length
        n = length(s1)

        # cut last data point if length is odd
        if mod(n, 2) == 1
          n -= 1
          s1 = s1[1:n]
          s2 = s2[1:n]
        end

        # remove means
        s1 .-= mean(s1)
        s2 .-= mean(s2)

        # Hann window
        hann = .5*(1 .- cos.(2π*(1:n)/n))

        # Fourier transform
        st1 = rfft(hann.*s1)
        st2 = rfft(hann.*s2)

        # sample frequencies
        ω = (0:n÷2)*fs/n

        # Gaussian for freq. averages
        G = exp.(-(ω.-frequencies').^2/2avgwidth^2)

        # filter
        st1f = st1.*G
        st2f = st2.*G

        # normalization
        n1 = (abs.(st1f[1,:]).^2 + 2sum(abs.(st1f[2:n÷2,:]).^2, dims=1)[1,:]
              + abs.(st1f[n÷2+1,:]).^2)/n
        n2 = (abs.(st2f[1,:]).^2 + 2sum(abs.(st2f[2:n÷2,:]).^2, dims=1)[1,:]
              + abs.(st2f[n÷2+1,:]).^2)/n

        # cross-correlation function
        cc = circshift(irfft(conj(st1f).*st2f, n, 1)./sqrt.(n1.*n2)', (n÷2, 0))

        # offsets
        lags = collect(-n÷2:n÷2-1)/fs

        # correct for difference in initial times
        lags .+= time2[1] - time1[1]

        # pick max CC and adjacent max
        lagc, lagr, lagl, ccc, ccr, ccl = findmaxcc(cc, frequencies, lags, reffreq, fs)
        i0 = argmin(abs.(frequencies .- reffreq))

        # save figure if CC ≥ 0.6
        if saveplot && ccc[i0] ≥ 0.6
          fig = figure()
          imshow(cc', cmap="RdBu_r", vmin=-1, vmax=1, origin="lower", aspect="auto",
                 extent=[lags[1] - 1/2fs, lags[end] + 1/2fs,
                         frequencies[1] - .5(frequencies[2]-frequencies[1]),
                         frequencies[end] + .5(frequencies[end]-frequencies[end-1])])
          plot(lagc, frequencies, color="black")
          plot(lagr, frequencies, color="black")
          plot(lagl, frequencies, color="black")
          xlim(lagc[i0]-1, lagc[i0]+1)
          xlabel("lag (s)")
          ylabel("frequency (Hz)")
          savefig(@sprintf("%s/%s_%s.pdf", figpath, pair.event1, pair.event2), dpi=200)
          close(fig)
        end

        # save to file
        h5open(filename, "w") do fid
          write(fid, "freq", collect(frequencies))
          write(fid, "lagc", lagc)
          write(fid, "lagr", lagr)
          write(fid, "lagl", lagl)
          write(fid, "ccc", ccc)
          write(fid, "ccr", ccr)
          write(fid, "ccl", ccl)
        end

      end

    end

  end

end

"Calculate great-circle distance"
function dist(lon1, lat1, lon2, lat2)
  Δσ = acos(sind(lat1)*sind(lat2) + cosd(lat1)*cosd(lat2)*cosd(lon1-lon2))
  return earthradius*Δσ
end

"Find maxima in CC starting at a reference frequency"
function findmaxcc(cc, frequencies, lags, reffreq, fs)

  # sample length
  n = length(lags)

  # number of central frequencies
  m = length(frequencies)

  # initialize lags, CCs
  lagc = Array{Float64,1}(undef, m)
  lagr = Array{Float64,1}(undef, m)
  lagl = Array{Float64,1}(undef, m)
  ccc = Array{Float64,1}(undef, m)
  ccr = Array{Float64,1}(undef, m)
  ccl = Array{Float64,1}(undef, m)

  # frequency closest to the reference frequency
  i0 = argmin(abs.(frequencies.-reffreq))

  # index of max CC at reference frequency
  imaxc = argmax(cc[:,i0])

  # find max by interpolation
  lagc[i0] = lags[imaxc] + argmaxquad(cc[mod1.(imaxc-1:imaxc+1, n),i0], 1/fs)
  ccc[i0] = maxquad(cc[mod1.(imaxc-1:imaxc+1, n),i0])

  # sweep frequencies
  sweep!(lagc, ccc, cc, frequencies, lags, i0, imaxc, fs)

  # index of right secondary max CC at reference frequency
  idx = mod1.(imaxc .+ (Int(round(fs/2reffreq)) : Int(round(3fs/2reffreq))), n)
  imaxr = idx[argmax(cc[idx,i0])]

  # find max by interpolation
  lagr[i0] = lags[imaxr] + argmaxquad(cc[mod1.(imaxr-1:imaxr+1, n),i0], 1/fs)
  ccr[i0] = maxquad(cc[mod1.(imaxr-1:imaxr+1, n),i0])

  # sweep frequencies
  sweep!(lagr, ccr, cc, frequencies, lags, i0, imaxr, fs)

  # index of left secondary max CC at reference frequency
  idx = mod1.(imaxc .+ (-Int(round(3fs/2reffreq)) : -Int(round(fs/2reffreq))), n)
  imaxl = idx[argmax(cc[idx,i0])]
  lagl[i0] = lags[imaxl] + argmaxquad(cc[mod1.(imaxl-1:imaxl+1, n),i0], 1/fs)
  ccl[i0] = maxquad(cc[mod1.(imaxl-1:imaxl+1, n),i0])

  # sweep frequencies
  sweep!(lagl, ccl, cc, frequencies, lags, i0, imaxl, fs)

  return lagc, lagr, lagl, ccc, ccr, ccl

end

"Get maximum of CC function by fitting a quadratic to three points"
maxquad(c) = c[2] + (c[1]-c[3])^2/(16c[2]-8(c[1]+c[3]))

"Get index of maximum of CC function by fitting a quadratic to three points"
argmaxquad(c, h) = h*(c[1]-c[3])/2(c[1]-2c[2]+c[3])

"Sweep across frequencies to trace max CC"
function sweep!(maxlags, ccs, cc, frequencies, lags, i0, imax, fs)

  # sample length
  n = length(lags)

  # number of central frequencies
  m = length(frequencies)

  # save reference index of maximum
  imaxref = imax

  # sweep up
  for i = i0+1:m

    # search region: ± half a period
    idx = mod1.(imax .+ (-Int(round(fs/2frequencies[i])) : Int(round(fs/2frequencies[i]))), n)

    # maximum on grid
    imax = idx[argmax(cc[idx,i])]

    # interpolated max
    maxlags[i] = lags[imax] + argmaxquad(cc[mod1.(imax-1:imax+1, n),i], 1/fs)
    ccs[i] = maxquad(cc[mod1.(imax-1:imax+1, n),i])

  end

  # restore reference index
  imax = imaxref

  # sweep down
  for i = i0-1:-1:1

    # search region: ± half a period
    idx = mod1.(imax .+ (-Int(round(fs/2frequencies[i])) : Int(round(fs/2frequencies[i]))), n)

    # maximum on grid
    imax = idx[argmax(cc[idx,i])]

    # interpolated max
    maxlags[i] = lags[imax] + argmaxquad(cc[mod1.(imax-1:imax+1, n),i], 1/fs)
    ccs[i] = maxquad(cc[mod1.(imax-1:imax+1, n),i])

  end

end

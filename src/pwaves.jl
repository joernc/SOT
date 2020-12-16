# TODO:
# - Write empty file when data is missing?
# - Prevent pairing with neighboring event?

# initialize model for ray prediction
global taup = pyimport("obspy.taup")
global model = taup.TauPyModel(model="prem")

"""
    downloadpwaves(eqname, stations)

Download the P waveforms for the event catalog specified by `eqname` and the stations (and
channels) specified in `stations`. Currently, one full hour after the event is downloaded.
"""
function downloadpwaves(eqname, stations; src="IRIS")

  # load event catalog
  events = DataFrame(CSV.read(@sprintf("data/catalogs/%s.csv", eqname)))

  # loop over stations
  for s in stations

    # path to P-waveform directory
    path = @sprintf("data/pdata/%s_%s", eqname, s)

    # create directory if needed
    mkpath(path)

    # loop over events
    for e in eachrow(events)

      # formatted date and time
      fmttime = Dates.format(e.time, "yyyy-mm-ddTHH:MM:SS.ss")

      # filename
      filename = @sprintf("%s/%s.h5", path, fmttime)

      # check if file exists already
      if !isfile(filename)

        @printf("Downloading %s ... ", fmttime)

        # download data
        try

          S = get_data("FDSN", s; s=string(e.time), t=string(e.time + Hour(1)), src)

          # check if there are no gaps
          if length(S) == 1

            # write to file
            h5open(filename, "w") do fid

              # station latitude
              fid["latitude"] = S[1].loc.lat

              # station longitude
              fid["longitude"] = S[1].loc.lon

              # start time in microseconds since 1970-01-01T00:00:00
              fid["starttime"] = S[1].t[1,2]

              # sampling frequency (in Hz)
              fid["fs"] = S[1].fs
              
              # trace
              fid["trace"] = S[1].x

            end

          end

          @printf("done\n")

        catch y

          @printf("failed\n")

        end
    
      end
  
    end
  
  end

end

"""
    cutpwaves(eqname, stations, interval, freqband)

Cut out and filter P-waveform sections used for cross-correlation. The data are read from
file based on `eqname` and `stations`, and are cut around the predicted P-wave arrival
time. The `interval` specifies what interval around that predicted arrival should be cut,
e.g. `interval = [-3, 47]` cuts a 50 s interval that starts 3 s before the predicted
arrival. The traces are filtered (prior to cutting) with a bandpass specified in
`freqband`, e.g. `freqband = [1, 3]` will filter to 1 to 3 Hz.
"""
function cutpwaves(eqname, stations, interval, freqband)

  # load event catalog
  events = DataFrame(CSV.read(@sprintf("data/catalogs/%s.csv", eqname)))

  for s in stations

    # path to P-waveform directory
    inpath = @sprintf("data/pdata/%s_%s", eqname, s)
    outpath = @sprintf("data/pwaves/%s_%s", eqname, s)

    # create directory if needed
    mkpath(outpath)

    for e in eachrow(events)

      # formatted date and time
      fmttime = Dates.format(e.time, "yyyy-mm-ddTHH:MM:SS.ss")

      # filename
      infile = @sprintf("%s/%s.h5", inpath, fmttime)
      outfile = @sprintf("%s/%s.h5", outpath, fmttime)

      println(fmttime)

      # check if input file is present and output file is not
      if isfile(infile) && !isfile(outfile)

        # read data
        fid = h5open(infile, "r")
        latitude = read(fid, "latitude")
        longitude = read(fid, "longitude")
        starttime = read(fid, "starttime")
        fs = read(fid, "fs")
        trace = read(fid, "trace")
        close(fid)

        # time stamps
        time = (1:length(trace))/fs

        # predict P-wave travel time
        d = dist(e.longitude, e.latitude, longitude, latitude)
        traveltime = model.get_travel_times(source_depth_in_km=0,
                                            distance_in_degree=rad2deg(d/earthradius),
                                            phase_list=["P"])[1].time

        # band pass filter
        responsetype = Bandpass(freqband[1], freqband[2]; fs)
        designmethod = Butterworth(4)
        filttrace = filtfilt(digitalfilter(responsetype, designmethod), trace)

        # time shift due to discrepancy between start time of time series and event time
        timeshift = (starttime - 10^3*(e.time - DateTime(1970, 1, 1)).value)

        # find times
        idx = interval[1] .< time .+ timeshift/1e6 .- traveltime .< interval[2]

        # check if data is present
        if sum(idx) > 1

          # cut
          cutstarttime = starttime + Int(round(10^6*time[idx][1]))
          cuttrace = filttrace[idx]

          # save to file
          h5open(outfile, "w") do fid
            fid["starttime"] = cutstarttime
            fid["trace"] = cuttrace
            fid["fs"] = fs
          end

        end

      end

    end

  end

end

"""
    findpairs(eqname, stations; saveplot=false)

Cross-correlate P waveforms to find repeating pairs. The events are read from catalogs
according to `eqname` and `stations`. Set `saveplot=true` if plots of the cross-correlation
functions are desired.
"""
function findpairs(eqname, stations; saveplot=false)

  # load event catalog
  events = DataFrame(CSV.read(@sprintf("data/catalogs/%s.csv", eqname)))

  # loop over stations
  for s in stations

    # path to P-wave directory
    datapath = @sprintf("data/pwaves/%s_%s", eqname, s)

    # path to plot directory
    plotpath = @sprintf("data/pplots/%s_%s", eqname, s)
    if saveplot
      mkpath(plotpath)
    end

    # read traces, sampling rates, and start times from file
    traces = Array{Float64,1}[]
    fs = Float64[]
    starttimes = Int[]
    present = falses(size(events, 1))
    for i = 1 : size(events, 1)
      fmttime = Dates.format(events[i,:time], "yyyy-mm-ddTHH:MM:SS.ss")
      filename = @sprintf("%s/%s.h5", datapath, fmttime)
      if isfile(filename)
        present[i] = true
        push!(traces, h5read(filename, "trace"))
        push!(fs, h5read(filename, "fs"))
        push!(starttimes, h5read(filename, "starttime"))
      end
    end

    # delete events without file
    events = events[present,:]

    # initialize pair catalog
    pairs = DataFrame(event1=DateTime[], event2=DateTime[], Δτ=Float64[], cc=Float64[])

    # loop over event pairs
    for i = 1 : size(events, 1) - 1, j = i + 1 : size(events, 1)

      # lengths of waveforms
      n1 = length(traces[i])
      n2 = length(traces[j])

      # pad with zeros
      padtrace1 = [traces[i]; zeros(n2)]
      padtrace2 = [traces[j]; zeros(n1)]

      # calculate cross-correlation function
      norm = sqrt(sum(traces[i].^2)*sum(traces[j].^2))
      cc = irfft(conj.(rfft(padtrace1)).*rfft(padtrace2), n1+n2)/norm

      # find index of max cross-correlation
      maxidx = argmax(cc)

      # calculate max CC using quadratic interpolation
      maxcc = maxquad(cc[mod1.(maxidx-1:maxidx+1, n1+n2)])

      # record if CC is above 0.9, sampling rate is identical, and pairs is not a repeater
      if maxcc ≥ 0.9 && fs[i] == fs[j] && starttimes[j] - starttimes[i] > n1/fs[i]*1e6

        # formatted event times
        fmttime1 = Dates.format(events[i,:time], "yyyy-mm-ddTHH:MM:SS.ss")
        fmttime2 = Dates.format(events[j,:time], "yyyy-mm-ddTHH:MM:SS.ss")

        # integer lags (depending on whether n1 + n2 is even or odd)
        lags = iseven(n1+n2) ? (-(n1+n2)÷2 : (n1+n2)÷2-1) : (-(n1+n2)÷2 : (n1+n2)÷2)

        # calculate lag from grid max CC
        Δτ = mod(maxidx-1, lags)/fs[i]

        # adjust using quadratic interpolation
        Δτ += argmaxquad(cc[mod1.(maxidx-1:maxidx+1, n1+n2)], 1/fs[i])

        # adjust for starttime recorded in waveforms
        originadjustment = starttimes[j] - starttimes[i] - 10^3*(events[j,:time]
                                                                 - events[i,:time]).value
        Δτ += originadjustment/1e6

        @printf("%s %s %4.2f %+6.3f\n", fmttime1, fmttime2, maxcc, Δτ)

        # record in catalog
        push!(pairs, [events[i,:time], events[j,:time], Δτ, maxcc])

        # plot cross-correlation function if desired
        if saveplot
          fig = figure()
          plot(lags/fs[i] .+ originadjustment/1e6, circshift(cc, (n1+n2)÷2))
          xlim(Δτ-5, Δτ+5)
          xlabel("lag (s)")
          ylabel("cross-correlation")
          savefig(@sprintf("%s/%s_%s.pdf", plotpath, fmttime1, fmttime2))
          close(fig)
        end

      end

    end
    
    # save catalog to file
    CSV.write(@sprintf("data/catalogs/%s_%s.csv", eqname, s), pairs)

  end

end

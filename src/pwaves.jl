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
function downloadseisdata(eqname, stations; src="IRIS")

  # load event catalog
  events = DataFrame(CSV.File(@sprintf("data/catalogs/%s.csv", eqname)))

  # loop over stations
  for s in stations

    # path to P-waveform directory
    path = seisdatadir(eqname, s)

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
    cutpwaves(eqname, stations, intervals, freqbands)

Cut out and filter P-waveform sections used for cross-correlation. The data are read from
file based on `eqname` and `stations`, and are cut around the predicted P-wave arrival
time. The `intervals` specify what intervals around that predicted arrivals should be cut,
e.g. `intervals = [[-3, 47], [-3, 47]]` cuts a 50 s interval that starts 3 s before the
predicted arrival for both P-wave stations to be processed. The traces are filtered (prior
to cutting) with a bandpass specified in `freqbands`, e.g. `freqbands = [[1, 3], [1.5, 2.5]`
will filter data from the first station to 1 to 3 Hz and that from the second station to
1.5 to 2.5 Hz.
"""
function cutpwaves(eqname, stations, intervals, freqbands)

  # load event catalog
  events = DataFrame(CSV.File(@sprintf("data/catalogs/%s.csv", eqname)))

  for i in 1 : length(stations)

    # path to P-waveform directory
    datapath = seisdatadir(eqname, stations[i])
    wavepath = seiswavedir(eqname, stations[i], intervals[i], freqbands[i])

    # create directory if needed
    mkpath(wavepath)

    for e in eachrow(events)

      # filename
      datafile = @sprintf("%s/%s.h5", datapath, fmttime(e.time))
      wavefile = @sprintf("%s/%s.h5", wavepath, fmttime(e.time))

      @printf("%s %s\n", stations[i], fmttime(e.time))

      # check if input file is present and output file is not
      if isfile(datafile) && !isfile(wavefile)

        # read data
        fid = h5open(datafile, "r")
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
        responsetype = Bandpass(freqbands[i][1], freqbands[i][2]; fs)
        designmethod = Butterworth(4)
        filttrace = filtfilt(digitalfilter(responsetype, designmethod), trace)

        # time shift due to discrepancy between start time of time series and event time
        timeshift = (starttime - 10^3*(e.time - DateTime(1970, 1, 1)).value)

        # find times
        idx = intervals[i][1] .< time .+ timeshift/1e6 .- traveltime .< intervals[i][2]

        # check if data is present
        if sum(idx) > 1

          # cut
          cutstarttime = starttime + Int(round(10^6*time[idx][1]))
          cuttrace = filttrace[idx]

          # save to file
          h5open(wavefile, "w") do fid
            fid["latitude"] = latitude
            fid["longitude"] = longitude
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
    findpairs(eqname, stations, intervals, freqbands; saveplot=false)

Cross-correlate P waveforms to find repeating pairs. The events are read from catalogs
according to `eqname` and `stations`, and the waveforms are read from folders based on
`eqname`, `stations`, `intervals`, `freqbands`. Set `saveplot=true` if plots of the
cross-correlation functions should be saved.
"""
function findpairs(eqname, stations, intervals, freqbands; saveplot=false)

  # load event catalog
  allevents = DataFrame(CSV.File(@sprintf("data/catalogs/%s.csv", eqname)))

  # loop over stations
  for i = 1 : length(stations)

    # path to P-wave directory
    wavepath = pwavedir(eqname, stations[i], intervals[i], freqbands[i])

    # path to plot directory
    plotpath = pplotdir(eqname, stations[i], intervals[i], freqbands[i])
    if saveplot
      mkpath(plotpath)
    end

    # read traces, sampling rates, and start times from file
    latitudes = Float64[]
    longitudes = Float64[]
    starttimes = Int[]
    fs = Float64[]
    traces = Array{Float64,1}[]
    present = falses(size(allevents, 1))
    for i = 1 : size(allevents, 1)
      filename = @sprintf("%s/%s.h5", wavepath, fmttime(allevents[i,:time]))
      if isfile(filename)
        present[i] = true
        push!(latitudes, h5read(filename, "latitude"))
        push!(longitudes, h5read(filename, "longitude"))
        push!(starttimes, h5read(filename, "starttime"))
        push!(fs, h5read(filename, "fs"))
        push!(traces, h5read(filename, "trace"))
      end
    end

    # delete events without file
    events = allevents[present,:]

    # initialize pair catalog
    pairs = DataFrame(event1=DateTime[], latitude1=Float64[], longitude1=Float64[],
                      event2=DateTime[], latitude2=Float64[], longitude2=Float64[],
                      Δτ=Float64[], cc=Float64[])

    # loop over event pairs
    for j = 1 : size(events, 1) - 1, k = j + 1 : size(events, 1)

      # lengths of waveforms
      n1 = length(traces[j])
      n2 = length(traces[k])

      # pad with zeros
      padtrace1 = [traces[j]; zeros(n2)]
      padtrace2 = [traces[k]; zeros(n1)]

      # calculate cross-correlation function
      norm = sqrt(sum(traces[j].^2)*sum(traces[k].^2))
      cc = irfft(conj.(rfft(padtrace1)).*rfft(padtrace2), n1+n2)/norm

      # find index of max cross-correlation
      maxidx = argmax(cc)

      # calculate max CC using quadratic interpolation
      maxcc = maxquad(cc[mod1.(maxidx-1:maxidx+1, n1+n2)])

      # record if CC is above 0.9, sampling rate is identical, and pairs is not a repeater
      if maxcc ≥ 0.9 && fs[j] == fs[k] && starttimes[k] - starttimes[j] > n1/fs[j]*1e6

        # integer lags (depending on whether n1 + n2 is even or odd)
        lags = iseven(n1+n2) ? (-(n1+n2)÷2 : (n1+n2)÷2-1) : (-(n1+n2)÷2 : (n1+n2)÷2)

        # calculate lag from grid max CC
        Δτ = mod(maxidx-1, lags)/fs[j]

        # adjust using quadratic interpolation
        Δτ += argmaxquad(cc[mod1.(maxidx-1:maxidx+1, n1+n2)], 1/fs[j])

        # adjust for starttime recorded in waveforms
        originadjustment = starttimes[k] - starttimes[j] - 10^3*(events[k,:time]
                                                                 - events[j,:time]).value
        Δτ += originadjustment/1e6

        @printf("%s %s %s %4.2f %+6.3f\n", stations[i], fmttime(events[j,:time]),
                fmttime(events[k,:time]), maxcc, Δτ)

        # record in catalog
        push!(pairs, [events[j,:time], events[j,:latitude], events[j,:longitude],
                      events[k,:time], events[k,:latitude], events[k,:longitude],
                      Δτ, maxcc])

        # plot cross-correlation function if desired
        if saveplot
          fig = figure()
          plot(lags/fs[j] .+ originadjustment/1e6, circshift(cc, (n1+n2)÷2))
          xlim(Δτ-5, Δτ+5)
          xlabel("lag (s)")
          ylabel("cross-correlation")
          savefig(@sprintf("%s/%s_%s.pdf", plotpath, fmttime(events[j,:time]),
                           fmttime(events[k,:time])))
          close(fig)
        end

      end

    end
    
    # save catalog to file
    CSV.write(paircatfile(eqname, stations[i], intervals[i], freqbands[i]), pairs)

  end

end

"Time format"
fmttime(time) = Dates.format(time, "yyyy-mm-ddTHH:MM:SS.ss")

"Directory for raw P-wave data"
seisdatadir(eqname, station) = @sprintf("data/seisdata/%s_%s", eqname, station)

"Directory for processed P waveforms"
pwavedir(eqname, station, interval, freqband) = @sprintf("data/pwaves/%s_%s_%+03d_%+03d_%3.1f_%3.1f", eqname, station, interval[1], interval[2], freqband[1], freqband[2])

"Directory for P-wave cross-correlation plots"
pplotdir(eqname, station, interval, freqband) = @sprintf("data/pplots/%s_%s_%+03d_%+03d_%3.1f_%3.1f", eqname, station, interval[1], interval[2], freqband[1], freqband[2])

"File name for P-wave pair catalogs"
paircatfile(eqname, station, interval, freqband) = @sprintf("data/catalogs/%s_%s_%+03d_%+03d_%3.1f_%3.1f.csv", eqname, station, interval[1], interval[2], freqband[1], freqband[2])

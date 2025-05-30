# TODO:
# - Prevent pairing with neighboring event?
# - Allow for gaps!

"""
    downloadpwaves(eqname, stations)

Download the P waveforms for the event catalog specified by `eqname` and the stations (and
channels) specified in `stations`. Currently, one full hour after the event is downloaded.
"""
function downloadseisdata(eqname, stations; src="IRIS", paircat=false)

  obspy = pyimport("obspy")
  fdsn = pyimport("obspy.clients.fdsn")
  client = fdsn.Client("IRIS")

  if paircat
    # load pair catalog
    pairs = DataFrame(CSV.File(@sprintf("data/catalogs/%s.csv", eqname)))
    # build catalog of all unique events in pair catalog
    events1 = pairs[:,[1;3:5]]
    events2 = pairs[:,[2;6:8]]
    rename!(events1, :event1=>:time, :latitude1=>:latitude, :longitude1=>:longitude, :depth1=>:depth)
    rename!(events2, :event2=>:time, :latitude2=>:latitude, :longitude2=>:longitude, :depth2=>:depth)
    events = sort(unique(vcat(events1, events2)))
  else
    # load event catalog
    events = DataFrame(CSV.File(@sprintf("data/catalogs/%s.csv", eqname)))
  end

  # loop over stations
  for (i, s) in enumerate(stations)
    # select source (station-dependent if given as vector)
    if typeof(src) == Vector{String}
      source = src[i]
    else
      source = src
    end
    # get station info
    sta = client.get_stations(network=split(s, ".")[1], station=split(s, ".")[2])[1][1]
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
        @printf("Downloading %s %s ... ", s, fmttime)
        try
          # download data
          st = client.get_waveforms(split(s, ".")..., obspy.UTCDateTime(string(e.time)), obspy.UTCDateTime(string(e.time + Hour(1))); attach_response=true)
          # check if there are no gaps
          if length(st.traces) == 1 && size(st[1].data, 1) > 0
            # remove instrument response
            st.remove_response()
            # write to file
            h5open(filename, "w") do fid
              # station latitude
              fid["latitude"] = sta.latitude
              # station longitude
              fid["longitude"] = sta.longitude
              # start time in microseconds since 1970-01-01T00:00:00
              fid["starttime"] = Int(round(1e6*(st[1].stats.starttime - obspy.UTCDateTime(1970, 1, 1, 0, 0, 0)).real))
              # sampling frequency (in Hz)
              fid["fs"] = 1/st[1].stats.delta
              # trace
              fid["trace"] = st[1].data
            end
            @printf("done\n")
          else
            # create empty file, so download is not attempted again
            touch(filename)
            @printf("%d gap(s)\n", length(st.traces)-1)
          end
        catch y
          # create empty file, so download is not attempted again
          if isa(y, PyCall.PyError) && y.T == obspy.clients.fdsn.header.FDSNNoDataException
            touch(filename)
          end
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
function cutpwaves(eqname, stations, intervals, freqbands; paircat=false)

  # initialize model for ray prediction
  taup = pyimport("obspy.taup")
  model = taup.TauPyModel(model="prem")

  if paircat

    # load pair catalog
    pairs = DataFrame(CSV.File(@sprintf("data/catalogs/%s.csv", eqname)))

    # build catalog of all unique events in pair catalog
    events1 = pairs[:,[1;3:5]]
    events2 = pairs[:,[2;6:8]]
    rename!(events1, :event1=>:time, :latitude1=>:latitude, :longitude1=>:longitude, :depth1=>:depth)
    rename!(events2, :event2=>:time, :latitude2=>:latitude, :longitude2=>:longitude, :depth2=>:depth)
    events = sort(unique(vcat(events1, events2)))

  else

    # load event catalog
    events = DataFrame(CSV.File(@sprintf("data/catalogs/%s.csv", eqname)))

  end

  for i in 1 : length(stations)

    # path to P-waveform directory
    datapath = seisdatadir(eqname, stations[i])
    wavepath = pwavedir(eqname, stations[i], intervals[i], freqbands[i])

    # create directory if needed
    mkpath(wavepath)

    for e in eachrow(events)

      # filename
      datafile = @sprintf("%s/%s.h5", datapath, fmttime(e.time))
      wavefile = @sprintf("%s/%s.h5", wavepath, fmttime(e.time))

      @printf("%s %s\n", stations[i], fmttime(e.time))

      # check if input file is present and output file is not
      if isfile(datafile) && !isfile(wavefile) && filesize(datafile) > 0

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
        traveltime = model.get_travel_times(source_depth_in_km=e.depth, distance_in_degree=rad2deg(d/earthradius), phase_list=["P", "p"])[1].time

        # time shift due to discrepancy between start time of time series and event time
        timeshift = (starttime - 10^3*(e.time - DateTime(1970, 1, 1)).value)

        # find times
        idx = intervals[i][1] .< time .+ timeshift/1e6 .- traveltime .< intervals[i][2]

        # check if data is present (25 data points needed for bandpass filter)
        if sum(idx) > 25

          # band pass filter
          responsetype = Bandpass(freqbands[i][1], freqbands[i][2])
          designmethod = Butterworth(4)
          filttrace = filtfilt(digitalfilter(responsetype, designmethod; fs), trace)

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
    pairs = DataFrame(event1=DateTime[], latitude1=Float64[], longitude1=Float64[], event2=DateTime[], latitude2=Float64[], longitude2=Float64[], Δτ=Float64[], cc=Float64[])

    # loop over event pairs
    for j = 1 : size(events, 1) - 1, k = j + 1 : size(events, 1)

      # ensure sampling rate is identical
      if fs[j] == fs[k]

        # cross-correlation measurement
        maxcc, Δτ = xcorr(traces[j], traces[k], fs[j])

        # record if CC is above 0.9 and pair is not a self-repeater
        if maxcc ≥ 0.9 && starttimes[k] - starttimes[j] > length(traces[j])/fs[j]*1e6

          # adjust for starttime recorded in waveforms
          originadjustment = starttimes[k] - starttimes[j] - 10^3*(events[k,:time] - events[j,:time]).value
          Δτ += originadjustment/1e6

          @printf("%s %s %s %5.3f %+6.3f\n", stations[i], fmttime(events[j,:time]), fmttime(events[k,:time]), maxcc, Δτ)

          # record in catalog
          push!(pairs, [events[j,:time], events[j,:latitude], events[j,:longitude], events[k,:time], events[k,:latitude], events[k,:longitude], Δτ, maxcc])

        end

      end

    end

    # save catalog to file
    CSV.write(paircatfile(eqname, stations[i], intervals[i], freqbands[i]), pairs)

  end

end

"""
    measurepairs(eqname, stations, intervals, freqbands; saveplot=false)

Cross-correlate P waveforms for cataloged pairs. The event pairs are read from catalogs
according to `eqname` and `stations`, and the waveforms are read from folders based on
`eqname`, `stations`, `intervals`, `freqbands`. Set `saveplot=true` if plots of the
cross-correlation functions should be saved.
"""
function measurepairs(eqname, stations, intervals, freqbands)

  # load event pair catalog
  allpairs = DataFrame(CSV.File(@sprintf("data/catalogs/%s.csv", eqname)))

  # loop over stations
  for i = 1 : length(stations)

    # path to P-wave directory
    wavepath = pwavedir(eqname, stations[i], intervals[i], freqbands[i])

    # initialize pair catalog
    pairs = DataFrame(event1=DateTime[], latitude1=Float64[], longitude1=Float64[], event2=DateTime[], latitude2=Float64[], longitude2=Float64[], Δτ=Float64[], cc=Float64[])

    # loop over event pairs
    for pair in eachrow(allpairs)

      # data file names
      filename1 = @sprintf("%s/%s.h5", wavepath, fmttime(pair.event1))
      filename2 = @sprintf("%s/%s.h5", wavepath, fmttime(pair.event2))

      @printf("%s %s %s ", stations[i], fmttime(pair.event1), fmttime(pair.event2))

      # ensure files are present
      if !isfile(filename1)
        println("file1 missing")
      elseif !isfile(filename2)
        println("file2 missing")
      else

        # read data from files
        latitude1 = h5read(filename1, "latitude")
        latitude2 = h5read(filename2, "latitude")
        longitude1 = h5read(filename1, "longitude")
        longitude2 = h5read(filename2, "longitude")
        starttime1 = h5read(filename1, "starttime")
        starttime2 = h5read(filename2, "starttime")
        fs1 = h5read(filename1, "fs")
        fs2 = h5read(filename2, "fs")
        trace1 = h5read(filename1, "trace")
        trace2 = h5read(filename2, "trace")

        # ensure sampling rate is identical
        if fs1 != fs2
          println("fs1 ≠ fs2")
        else

          # cross-correlation measurement
          maxcc, Δτ = xcorr(trace1, trace2, fs1)

          # record if CC is above 0.9 and pair is not a self-repeater
          if maxcc < 0.9
            @printf("%5.3f\n", maxcc)
          elseif starttime2 - starttime1 ≤ length(trace1)/fs1*1e6
            println("suspected self-repeater")
          else

            # adjust for starttime recorded in waveforms
            originadjustment = starttime2 - starttime1 - 10^3*(pair.event2 - pair.event1).value
            Δτ += originadjustment/1e6

            @printf("%5.3f %+6.3f\n", maxcc, Δτ)

            # record in catalog
            push!(pairs, [pair.event1, pair.latitude1, pair.longitude1, pair.event2, pair.latitude2, pair.longitude2, Δτ, maxcc])

          end

        end

      end

    end
    
    # save catalog to file
    CSV.write(paircatfile(eqname, stations[i], intervals[i], freqbands[i]), pairs)

  end

end

"Cross-correlation measurement"
function xcorr(trace1, trace2, fs; maxlag=12)

  # lengths of waveforms
  n1 = length(trace1)
  n2 = length(trace2)

  # pad with zeros
  padtrace1 = [trace1; zeros(n2)]
  padtrace2 = [trace2; zeros(n1)]

  # calculate covariance function
  cov = irfft(conj.(rfft(padtrace1)).*rfft(padtrace2), n1+n2)

  # indicator functions
  ind1 = [ones(n1); zeros(n2)]
  ind2 = [ones(n2); zeros(n1)]

  # calculate normalization factors
  norm1 = irfft(conj.(rfft(padtrace1.^2)).*rfft(ind2), n1+n2)
  norm2 = irfft(conj.(rfft(ind1)).*rfft(padtrace2.^2), n1+n2)

  # exclude negative and zero normalization factors (due to roundoff error and non-overlapping segments)
  norm1[norm1.≤0] .= Inf
  norm2[norm2.≤0] .= Inf

  # cross-correlation
  cc = cov./sqrt.(norm1.*norm2)

  # find index of max covariance in allowable range (±maxlag, in seconds)
  k = Int(round(maxlag*fs))
  cca = copy(cc); cca[k+2:n1+n2-k] .= 0
  maxidx = argmax(cca)

  # check whether CC = 0
  if cc[maxidx] == 0

    return 0.0, NaN

  else

    # integer lags (depending on whether n1 + n2 is even or odd)
    lags = iseven(n1+n2) ? (-(n1+n2)÷2 : (n1+n2)÷2-1) : (-(n1+n2)÷2 : (n1+n2)÷2)

    # index shift from grid max covariance
    Δi = mod(maxidx-1, lags)

    # calculate max CC using quadratic interpolation
    maxcc = maxquad(cc[mod1.(maxidx-1:maxidx+1, n1+n2)])

    # time shift adjusted using quadratic interpolation
    Δτ = Δi/fs + argmaxquad(cc[mod1.(maxidx-1:maxidx+1, n1+n2)], 1/fs)

    return maxcc, Δτ

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

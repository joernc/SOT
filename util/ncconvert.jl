# H08 ECCO

using NCDatasets, HDF5

fl = Dataset("../shirui/dtaus_H08S2.nc", "r")

t = fl["time"][:]
freq = [2 3 4]
τ = float.(fl["__xarray_dataarray_variable__"][:,[1,3,4]])

tr = Dates.value.(t .- DateTime(2000, 1, 1))/1000/3600/24

h5open("data/ecco/nias_H08.h5", "w") do file
  write(file, "t", tr)
  write(file, "freq", freq)
  write(file, "tau", τ)
end

# H08 Argo

fl = Dataset("../shirui/dtaus_H08S2_argo.nc", "r")

t = DateTime(2004, 1, 15) .+ Month.(floor.(fl["time"][:]))
freq = [2 3 4]
τ = float.(fl["__xarray_dataarray_variable__"][:,[1,3,4]])

tr = Dates.value.(t .- DateTime(2000, 1, 1))/1000/3600/24

h5open("data/argo/nias_H08.h5", "w") do file
  write(file, "t", tr)
  write(file, "freq", freq)
  write(file, "tau", τ)
end

# H08 kernels

fl = Dataset("../shirui/KTs_H08S2.nc", "r")

range = float.(fl["x"][:])
depth = -fl["z"][:]
freq = [2 3 4]
K = float.(replace(fl["__xarray_dataarray_variable__"][:,:,[1,3,4]], missing=>0))

h5open("data/kernels/nias_H08.h5", "w") do file
  write(file, "range", range)
  write(file, "depth", depth)
  write(file, "freq", freq)
  write(file, "K", K)
end

# H01 ECCO

fl = Dataset("../shirui/dtaus_Nias_H01.nc", "r")

t = fl["time"][:]
freq = [2.5 3.25 4.0]
τ = float.(fl["__xarray_dataarray_variable__"][:,:])

tr = Dates.value.(t .- DateTime(2000, 1, 1))/1000/3600/24

h5open("data/ecco/nias_H01.h5", "w") do file
  write(file, "t", tr)
  write(file, "freq", freq)
  write(file, "tau", τ)
end

# H01 Argo

fl = Dataset("../shirui/dtaus_Nias_H01_argo.nc", "r")

t = DateTime(2004, 1, 15) .+ Month.(floor.(fl["time"][:]))
freq = [2.5 3.25 4.0]
τ = float.(fl["__xarray_dataarray_variable__"][:,:])

tr = Dates.value.(t .- DateTime(2000, 1, 1))/1000/3600/24

h5open("data/argo/nias_H01.h5", "w") do file
  write(file, "t", tr)
  write(file, "freq", freq)
  write(file, "tau", τ)
end

# H01 kernels

fl = Dataset("../shirui/KTs_Nias_H01.nc", "r")

range = float.(fl["x"][:])
depth = -fl["z"][:]
freq = [2.5 3.75 4.0]
K = float.(replace(fl["__xarray_dataarray_variable__"][:,:,:], missing=>0))

h5open("data/kernels/nias_H01.h5", "w") do file
  write(file, "range", range)
  write(file, "depth", depth)
  write(file, "freq", freq)
  write(file, "K", K)
end

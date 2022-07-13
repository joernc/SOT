using NCDatasets, HDF5, Dates

# H08

# kernels
fl = Dataset("../shirui/kernels/KTs_Nias_H08_coarsen.nc", "r")
K = Array{Float64}(replace(fl["SEMkernels_T"][:], missing=>0))
Km = Array{Float64}(replace(fl["MODEkernels_T"][:], missing=>0))
ω = Float64.(fl["f"][:])
x = Array{Float64}(fl["x"][:])
z = Array{Float64}(fl["z"][:])

# Argo travel time anomalies
fl = Dataset("../shirui/argo/dtaus_Nias_H08_argo_coarsenKTs.nc", "r")
targo = DateTime(2004, 1, 15) .+ Month.(floor.(fl["t"][:]))
τargo = float.(fl["SEMdtaus"][:,:,2])

# ECCO travel time anomalies
fl = Dataset("../shirui/ecco/dtaus_Nias_H08_ecco_coarseKTs.nc", "r")
tecco = fl["t"][:]
τecco = float.(fl["SEMdtaus"][:,:,2])

# Argo range-averaged temperature anomalies
fl = Dataset("../shirui/argo/dTz_Nias_H08_argo.nc", "r")
Targo = Array{Float64}(replace(fl["__xarray_dataarray_variable__"][:], missing=>NaN))
@assert DateTime(2004, 1, 15) .+ Month.(floor.(fl["t"][:])) == targo
@assert fl["z"][:] == z

## ECCO range-averaged temperature anomalies
#fl = Dataset("../shirui/ecco/dTz_Nias_H08_ecco.nc", "r")
#Tecco = Array{Float64}(replace(fl["__xarray_dataarray_variable__"][:], missing=>NaN))
#@assert fl["time"][:] == tecco
#@assert fl["z"][:] == z

# time in days since 2000-01-01T00:00:00
targo = Dates.value.(targo .- DateTime(2000, 1, 1))/1000/3600/24
tecco = Dates.value.(tecco .- DateTime(2000, 1, 1))/1000/3600/24

# save to file
h5open("data/temperature/nias_H08.h5", "w") do file
  write(file, "K", K)
  write(file, "Km", Km)
  write(file, "x", x)
  write(file, "z", z)
  write(file, "ω", ω)
  write(file, "targo", targo)
  write(file, "τargo", τargo)
  write(file, "Targo", Targo)
  write(file, "tecco", tecco)
  write(file, "τecco", τecco)
#  write(file, "Tecco", Tecco)
end

# H01

# kernels
fl = Dataset("../shirui/kernels/KTs_Nias_H01_coarsen.nc", "r")
K = Array{Float64}(replace(fl["SEMkernels_T"][:], missing=>0))
Km = Array{Float64}(replace(fl["MODEkernels_T"][:], missing=>0))
ω = Float64.(fl["f"][:])
x = Array{Float64}(fl["x"][:])
z = Array{Float64}(fl["z"][:])

# Argo travel time anomalies
fl = Dataset("../shirui/argo/dtaus_Nias_H01_argo_coarsenKTs.nc", "r")
targo = DateTime(2004, 1, 15) .+ Month.(floor.(fl["t"][:]))
τargo = float.(fl["SEMdtaus"][:,:,2])

# ECCO travel time anomalies
fl = Dataset("../shirui/ecco/dtaus_Nias_H01_ecco_coarseKTs.nc", "r")
tecco = fl["t"][:]
τecco = float.(fl["SEMdtaus"][:,:,2])

# Argo range-averaged temperature anomalies
fl = Dataset("../shirui/argo/dTz_Nias_H01_argo.nc", "r")
Targo = Array{Float64}(replace(fl["__xarray_dataarray_variable__"][:], missing=>NaN))
@assert DateTime(2004, 1, 15) .+ Month.(floor.(fl["t"][:])) == targo
@assert fl["z"][:] == z

## ECCO range-averaged temperature anomalies
#fl = Dataset("../shirui/ecco/dTz_Nias_H01_ecco.nc", "r")
#Tecco = Array{Float64}(replace(fl["__xarray_dataarray_variable__"][:], missing=>NaN))
#@assert fl["time"][:] == tecco
#@assert fl["z"][:] == z

# time in days since 2000-01-01T00:00:00
targo = Dates.value.(targo .- DateTime(2000, 1, 1))/1000/3600/24
tecco = Dates.value.(tecco .- DateTime(2000, 1, 1))/1000/3600/24

# save to file
h5open("data/temperature/nias_H01.h5", "w") do file
  write(file, "K", K)
  write(file, "Km", Km)
  write(file, "x", x)
  write(file, "z", z)
  write(file, "ω", ω)
  write(file, "targo", targo)
  write(file, "τargo", τargo)
  write(file, "Targo", Targo)
  write(file, "tecco", tecco)
  write(file, "τecco", τecco)
#  write(file, "Tecco", Tecco)
end

# DGAR

# kernels
fl = Dataset("../shirui/kernels/KTs_Nias_DGAR_coarsen.nc", "r")
K = Array{Float64}(replace(fl["SEMkernels_T"][:], missing=>0))
Km = Array{Float64}(replace(fl["MODEkernels_T"][:], missing=>0))
ω = Float64.(fl["f"][:])
x = Array{Float64}(fl["x"][:])
z = Array{Float64}(fl["z"][:])

# Argo travel time anomalies
fl = Dataset("../shirui/argo/dtaus_Nias_DGAR_argo_coarsenKTs.nc", "r")
targo = DateTime(2004, 1, 15) .+ Month.(floor.(fl["t"][:]))
τargo = float.(fl["SEMdtaus"][:,:,2])

# ECCO travel time anomalies
fl = Dataset("../shirui/ecco/dtaus_Nias_DGAR_ecco_coarseKTs.nc", "r")
tecco = fl["t"][:]
τecco = float.(fl["SEMdtaus"][:,:,2])

# Argo range-averaged temperature anomalies
fl = Dataset("../shirui/argo/dTz_Nias_DGAR_argo.nc", "r")
Targo = Array{Float64}(replace(fl["__xarray_dataarray_variable__"][:], missing=>NaN))
@assert DateTime(2004, 1, 15) .+ Month.(floor.(fl["t"][:])) == targo
@assert fl["z"][:] == z

## ECCO range-averaged temperature anomalies
#fl = Dataset("../shirui/ecco/dTz_Nias_DGAR_ecco.nc", "r")
#Tecco = Array{Float64}(replace(fl["__xarray_dataarray_variable__"][:], missing=>NaN))
#@assert fl["time"][:] == tecco
#@assert fl["z"][:] == z

# time in days since 2000-01-01T00:00:00
targo = Dates.value.(targo .- DateTime(2000, 1, 1))/1000/3600/24
tecco = Dates.value.(tecco .- DateTime(2000, 1, 1))/1000/3600/24

# save to file
h5open("data/temperature/nias_DGAR.h5", "w") do file
  write(file, "K", K)
  write(file, "Km", Km)
  write(file, "x", x)
  write(file, "z", z)
  write(file, "ω", ω)
  write(file, "targo", targo)
  write(file, "τargo", τargo)
  write(file, "Targo", Targo)
  write(file, "tecco", tecco)
  write(file, "τecco", τecco)
#  write(file, "Tecco", Tecco)
end

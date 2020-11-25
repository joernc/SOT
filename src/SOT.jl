module SOT

using LinearAlgebra, Dates, DataFrames, Seis, Printf, HDF5, FFTW, PyPlot, Statistics, ProgressMeter, CSV

export twavepick, invert, lineartrend

include("twavepick.jl")
include("inversion.jl")

end # module

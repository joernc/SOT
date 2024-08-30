using SOT, Printf, LinearAlgebra, HDF5, DataFrames, CSV, Dates, PyPlot,SparseArrays,Random, Distributions,NCDatasets
using Interpolations, Statistics, StatsBase
source = "nias_tm"
station = "H08"

# load data
#h5file = @sprintf("data/japan_%s_x2y_%dex.h5",station,nexclude)
ppfile = @sprintf("results/%s_%s_ppairs.csv",source,station)
tpvfile = @sprintf("results/%s_%s_tpairs.csv",source,station)

# frequencies used in inversion
tinvfreq = 2.0:1.0:4.0

# number of frequencies
l = length(tinvfreq)
i = 0
# manually exclude pairs
tpairs = CSV.read(tpvfile, DataFrame)
ppairs = CSV.read(ppfile, DataFrame)
#t0,t1 = DateTime(2005,1,1),DateTime(2006,1,1)
#tpairs = tpairs[(t1.> tpairs.event1.≥ t0),:]
#tpairs = tpairs[(tpairs.event2.< t1),:]
ppairs = innerjoin(ppairs, select(tpairs, [:event1, :event2]), on=[:event1, :event2])

tstation = [-7.65, 72.49]
pstations = ["PS.PSI..BHZ","MY.KUM..BHZ","II.WRAB.00.BHZ", "GE.GSI..BHZ"]
pstnlats,pstnlons = [2.69,5.2902,-19.934,1.3039],[98.92,100.6492,134.36,97.5755]
pstations = DataFrame(station=pstations,slat=pstnlats,slon=pstnlons)

# number of good T- and P-wave pairs
nt = size(tpairs, 1)
np = size(ppairs, 1)

# number of unique events
m = length(unique([tpairs.event1; tpairs.event2]))

###
model,grid = "ecco","Gaussian"
fl = @sprintf("results/covtau0_%s_%s.h5",model,grid)
lags = h5read(fl, "lags")
ctau = h5read(fl, "covs")

σc = ones(6)
# find unique events
t = sort(unique([tpairs.event1; tpairs.event2]))

# real time (days)
tr = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
#tr = 1.0:1.0:365.0
C = []
for i = 1:6
  cij = ctau[:,i]
  etpf = linear_interpolation(lags, cij, extrapolation_bc=0)
  push!(C,etpf.( tr.-tr'))
end
C = Matrix([σc[1]^2*C[1] σc[2]^2*C[2] σc[3]^2*C[3]; σc[2]^2*C[2]' σc[4]^2*C[4] σc[5]^2*C[5]; σc[3]^2*C[3]' σc[5]^2*C[5]' σc[6]^2*C[6]])
isposdef((inv(Symmetric(C))))

# noise (s)
σn,σnp = 0.014175,1.7e-3

# origin time correction standard deviation (s)
σp = 1.27

# trend priors for coefficients of travel time anomalies (s/day)
σtrend = 1e-2*ones(l)/SOT.meanyear

σs,σh = 0.0263628899,0.021

# get inversion matrices
t, E, C, Rcsn, Rcsnh = SOT.getEIO(tpairs, ppairs, tstation, pstations; l,σtrend,ctau,lags,getC=true,dcov=false)

C = Symmetric([σc[1]^2*C[1] σc[2]^2*C[2] σc[3]^2*C[3]; σc[2]^2*C[2]' σc[4]^2*C[4] σc[5]^2*C[5]; σc[3]^2*C[3]' σc[5]^2*C[5]' σc[6]^2*C[6]])
# transformation from singular vectors to observed frequencies
#T = kron(sparse(V'*Diagonal(sqrt.(Γ))),I(m))
#C = T*C*T'

#y = [reshape(vcat([ [parse(Float64, ss) for ss in split(tpairs.Δτ[i][2:end-1],", ")][1] for i = 1:nt]...), nt); ppairs.Δτ]
y = [reshape(vcat([ [parse(Float64, ss) for ss in split(tpairs.Δτ[i][2:end-1],", ")]' for i = 1:nt]...), l*nt); ppairs.Δτ]

C=C[i*m+1:(i+l)*m,i*m+1:(i+l)*m]

@printf("%d pairs, %d events\n",nt,m)
f(x) = -1 .* SOT.loglikelihoodIO(x, y, t, E, Rcsn, Rcsnh, l, nt, np;σtrend,C,hydro=true,finds=true)

## BFGS
x0 = log.([σp,σn,σnp,σh, 1.0, σs])
xk,Hk = SOT.BFGS(f,x0,5e-5,100; xrtol=1e-6)

@printf("MLE: %s\n",exp.(xk))
@printf("lb: %s\n ub: %s\n",exp.(xk.-2*sqrt.(diag(Hk))),exp.(xk.+2*sqrt.(diag(Hk))))


#= ### sediment kernel
MLE: [1.2719791639525246, 0.014130256015823685, 0.0017231361425545742, 0.021207243614479142, 1.0, 0.026040194418779264]
lb: [1.2157419667817706, 0.013842195425683029, 0.0016270663773979804, 0.017010835562479565, 0.1353352832366127, 0.02442966103781115]
ub: [1.330817753879337, 0.014424311240560976, 0.001824878325201598, 0.026438864808961175, 7.38905609893065, 0.027756902738777355]
### sedimentless kernel
MLE: [1.2721241577018088, 0.01351474939046164, 0.0016532357804016843, 0.026643055221222522, 1.0, 0.028080211154085197]
lb: [1.2180612477552393, 0.013230526748541948, 0.001565515069819472, 0.022539491393714892, 0.1353352832366127, 0.026408543763459687]
ub: [1.3285866171269265, 0.013805077799121776, 0.0017458717570285316, 0.031493718253068216, 7.38905609893065, 0.02985769550644516]
 =#
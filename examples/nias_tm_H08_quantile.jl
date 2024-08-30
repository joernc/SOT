using SOT, Printf, LinearAlgebra, HDF5, DataFrames, CSV, Dates, PyPlot,SparseArrays,Random, Distributions,NCDatasets
using Interpolations, Statistics, StatsBase
source = "nias_tm"
station = "H08"

grid = "coarsen"
fl = @sprintf("results/%s_%s_%s.h5",source,station,grid)
zix = h5read(fl, "zix2")
println(std(zix))
#zix ./= std(zix)
zint = h5read(fl, "zint")
zinp = h5read(fl, "zinp")
println(std(zint))
println(std(zinp))
grid = "Gaussian"
fl = @sprintf("results/%s_%s_%s.h5",source,station,grid)
zixg = h5read(fl, "zix2")
#zixg ./= std(zixg)
println(std(zixg))
zintg = h5read(fl, "zint")
zinpg = h5read(fl, "zinp")
println(std(zintg))
println(std(zinpg))
nx,nt,np = length(zix),length(zint),length(zinp)

p = 0.5:0.5:99.5
α = 0.05
qs4,qs5,qs6 = percentile(zix,p),percentile(zint,p),percentile(zinp,p)
qs1,qs2,qs3 = percentile(zixg,p),percentile(zintg,p),percentile(zinpg,p)
qt = quantile.(Normal(), p/100)
lb,ub = zeros(length(p),3),zeros(length(p),3)
for i = 1:length(p)
    for (j,ns) in enumerate([nx, nt, np])
        d = NoncentralT(ns-1,-sqrt(ns)*qt[i])
        lb[i,j],ub[i,j] = -quantile(d, 1-α/2)/sqrt(ns),-quantile(d, α/2)/sqrt(ns)
    end
end
# small font
rc("font", size=10)
rc("axes", titlesize="medium")
PyPlot.matplotlib[:rc]("mathtext",fontset="cm")        #computer modern font 
PyPlot.matplotlib[:rc]("font",family="STIXGeneral")
fig,ax=plt.subplots(figsize=(4.8,3.6))
ax.plot(qt,qs4-qt,label="2005 daily anomaly",c="tab:blue")
ax.plot(qt,qs1-qt,ls="--",c="tab:blue")
ax.fill_between(qt, (lb[:,1]-qt), (ub[:,1]-qt), alpha=.2, zorder=3, color="tab:blue", linewidth=0)
ax.plot(qt,qs5-qt,label="\$T\$-wave residual",c="tab:orange")
ax.plot(qt,qs2-qt,ls="--",c="tab:orange")
ax.fill_between(qt, (lb[:,2]-qt), (ub[:,2]-qt), alpha=.2, zorder=3, color="tab:orange", linewidth=0)
ax.plot(qt,qs6-qt,label="\$P\$-wave residual",c="tab:green")
ax.plot(qt,qs3-qt,ls="--",c="tab:green")
ax.fill_between(qt, (lb[:,3]-qt), (ub[:,3]-qt), alpha=.2, zorder=3, color="tab:green", linewidth=0)
ax.axhline(0,color="black",ls=":",lw=1)
ax.set_xlabel("theoretical quantile")
ax.set_ylabel("sample quantile \$-\$ theoretical quantile")
ax.legend(frameon=false)
fig.tight_layout()
fig.savefig(@sprintf("results/%s_%s_quantile.pdf",source,station))
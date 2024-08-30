using PyPlot, Printf, Dates, LinearAlgebra, Statistics, SparseArrays, DataFrames, CSV,FFTW,NCDatasets
using HDF5, Interpolations, Distributions, SeisNoise, StatsBase

nmode = 5
Lx = 29e5

nanmean(x) = mean(filter(!isnan,x))
nanmean(x,y) = mapslices(nanmean,Array{Float64}(x);dims=y)
fl = Dataset("data/ecco/n2_Nias_H08_2005.nc", "r")
range = Array{Float64}(fl["x"][:])
depth = -Array{Float64}(fl["z"][:])
xidx = 0 .<= range .<= Lx
nx = length(range[xidx])
Δx = range[2] - range[1]
Δz = depth[1] - depth[2]
giofes = NCDataset("results/ofes/modes_mean_ofes.nc")["hn"][:,1:nmode]
n2ofes = nanmean(replace(NCDataset("data/ofes/n2_Nias_H08_2005.nc")["n2"][:,xidx], missing=>NaN), 2)
giecco = NCDataset("results/ecco/modes_mean_ecco.nc")["hn"][:,1:nmode]
n2ecco = nanmean(replace(NCDataset("data/ecco/n2_Nias_H08_2005.nc")["n2"][:,xidx], missing=>NaN), 2)

K = sum(h5read("data/temperature/nias_H08.h5", "K"), dims=2)[:,1,:]*Δx

Aofes = Array{Float64, 2}(undef, 3, nmode)
Aecco = Array{Float64, 2}(undef, 3, nmode)
N0 = 1e-3
for i = 1:3
  for j = 1:nmode
    Aofes[i,j] = sum(replace(K[:,i].*n2ofes.*giofes[:,j],missing => 0))*Δz/N0
    Aecco[i,j] = sum(replace(K[:,i].*n2ecco.*giecco[:,j],missing => 0))*Δz/N0
  end
end

# SVD
Uo, Λo, Vo = svd(Aofes)
Ue, Λe, Ve = svd(Aecco)
println(Uo)
println(Λo)
println(Vo)


fl = Dataset("results/ofes/KjdT_Nias_H08_2005.nc", "r")
τofes = Array{Float64}(replace(fl["hndsig1"][:,:], missing=>0))

for i = 1:3
  τofes[:,i] = detrend(τofes[:,i] .- mean(τofes[:,i]))
end

# project Argo and ECCO travel time anomalies onto singular vectors
#cofes = (τofes)*(Uo*inv(Diagonal(Λo)))

PyPlot.matplotlib[:rc]("mathtext",fontset="cm")        #computer modern font 
PyPlot.matplotlib[:rc]("font",family="STIXGeneral")#,size=14)  #font similar to LaTeX
# small font
rc("font", size=10)
rc("axes", titlesize="medium")

function solvequadratic(a, b, c)
    d = sqrt(b^2 - 4a*c)
    return [(-b - d) / 2a, (-b + d) / 2a]
end

function make_FFT1D(field,n,Δ)
    signal=field
    FFT=rfft(signal)
    freq1=rfftfreq(n,Δ)
    
    om=freq1
    PSD=abs.(2*FFT.*conj(FFT)/n)
    
    return PSD,FFT,om
end

#N0 = 1e3
Kgrid = "coarsen"
flnm = @sprintf("results/anomalies/nias_tm_H08_%s.h5",Kgrid)
# load time series
t = DateTime(2000, 1, 1) .+ Millisecond.(h5read(flnm, "ti"))
ci = h5read(flnm, "ci")
cofes = h5read(flnm, "cofes")
cecco = h5read(flnm, "cecco")
ec = h5read(flnm, "ec")

#ci = ci .- mean(ci,dims=1)
#cofes = cofes .- mean(cofes,dims=1)
#cecco = cecco .- mean(cecco,dims=1)

# estimate ECCO and Argo trends
#Et, Rxxt, Rnnt, Pt, Mt = SOT.estimatetrend(t, λ, σc, σc, σtrend, σannual, σsemiannual)

n = length(t)
# sampling interval
Δ = 1

clrs = ["#a6cee3","tab:blue", "tab:orange", "tab:green"]

fig0, ax = subplots(2, 1, sharex=true, figsize=(5.6, 3.6))
for i = 1:2
    @printf("c%d ECCO overlap T-wave fraction: %.2f\n",i,sum((ci[:,i] - cecco[:,i]) .< 2*ec[:,i])/n)
    global p2,p3,p4
    #p1, = ax[i].plot(t, cix[:,i], zorder=1, color=clrs[1])
    p2, = ax[i].plot(t, 1e3sign(1.5-i)*ci[:,i], zorder=3, linewidth=1, color=clrs[2])
    ax[i].fill_between(t, 1e3(sign(1.5-i)*ci[:,i]-2*ec[:,i]), 1e3(sign(1.5-i)*ci[:,i]+2*ec[:,i]), alpha=.2, zorder=3, color=clrs[2], linewidth=0)
    p3, = ax[i].plot(t, 1e3sign(1.5-i)*cofes[:,i], zorder=1, linewidth=1, color=clrs[3])
    p4, = ax[i].plot(t, 1e3sign(1.5-i)*cecco[:,i], zorder=2, linewidth=1, color=clrs[4])
    ax[i].set_ylabel("\$a_$i~(\\mathrm{mK})\$")
    ylim = maximum(abs.(ax[i].get_ylim()))
    ax[i].set_ylim(-ylim, ylim)
    ax[i].set_title("($(('a':'z')[i]))", loc="left")
    ax[i].axhline(0, linestyle="--", linewidth=1, color="k", alpha=0.2)
end
ax[1].legend([p2, p3, p4], ["\$T\$ waves", "OFES", "ECCO"], ncol=3, loc="lower center", frameon=false)
ax[2].set_xlim(t[1]-Day(5), t[end]+Day(5))
fig0.align_ylabels()
fig0.tight_layout()
fig0.savefig(@sprintf("results/ofes/timeseries_ci_%s.pdf",Kgrid))

# convert time of samples to seconds since event time
time = 1:Δ:n

# Hann window
hann = .5*(1 .- cos.(2π*(0:n-1)/n))

n1 = n ÷ 2 + 1
idx1 = time .< n1+1
idx2 = time .> n-n1
idx3 = n1 ÷ 2 .< time .≤ n1 ÷ 2 + n1
idxsub = [idx1,idx2,idx3]

# Hann window
hann1 = .5*(1 .- cos.(2π*(1:n1)/n1))
edf,α = 3*36/19,1-0.95
d = Chisq(edf)
cl95 = 1-α^(1/(edf-1))
ci95 = [edf/quantile(d,1-α/2),edf/quantile(d,α/2)]

fig1, ax1 = subplots(2, 2, figsize=(5.6, 4.8), sharex="col")    
fig2, ax2 = subplots(2, 2, figsize=(190/25.4, 190/25.4*2/3), sharex=true)
# Fourier transform
for i = 1:2
    si = detrend(ci[:,i])
    #six = cix[:,i]
    sofes = detrend(cofes[:,i])
    secco = detrend(cecco[:,i])

    giisub = Array{Float64, 2}(undef, 3, 1+n1÷2)
    stisub = Array{Complex{Float64}, 2}(undef, 3, 1+n1÷2)
    #giixsub = Array{Float64, 2}(undef, 3, 1+n1÷2)
    #stixsub = Array{Complex{Float64}, 2}(undef, 3, 1+n1÷2)
    goosub = Array{Float64, 2}(undef, 3, 1+n1÷2)
    stosub = Array{Complex{Float64}, 2}(undef, 3, 1+n1÷2)
    geesub = Array{Float64, 2}(undef, 3, 1+n1÷2)
    stesub = Array{Complex{Float64}, 2}(undef, 3, 1+n1÷2)
    freqsub = Array{Float64, 2}(undef, 3, 1+n1÷2)
    for j = 1:3
        giisub[j,:],stisub[j,:],freqsub[j,:] = make_FFT1D(sqrt(8/3)*hann1.*si[idxsub[j]],n1,Δ)
        #giixsub[j,:],stixsub[j,:],freqsub[j,:] = make_FFT1D(sqrt(8/3)*hann1.*six[idxsub[j]],n1,Δ)
        goosub[j,:],stosub[j,:],freqsub[j,:] = make_FFT1D(sqrt(8/3)*hann1.*sofes[idxsub[j]],n1,Δ)
        geesub[j,:],stesub[j,:],freqsub[j,:] = make_FFT1D(sqrt(8/3)*hann1.*secco[idxsub[j]],n1,Δ)
    end

    idx = freqsub[1,:] .> 0.02
    #ax1[1,i].loglog(freqsub[1,idx],mean(giixsub,dims=1)[idx],zorder=1, color=clrs[1],label="\$T\$ waves, \$C_m\$")
    ax1[i,1].loglog(freqsub[1,idx],1e6*mean(giisub,dims=1)[idx],zorder=1, color=clrs[2],label="\$T\$ waves")
    ax1[i,1].loglog(freqsub[1,idx],1e6*mean(goosub,dims=1)[idx],zorder=0, color=clrs[3],label="OFES")
    ax1[i,1].loglog(freqsub[1,idx],1e6*mean(geesub,dims=1)[idx],zorder=0, color=clrs[4],label="ECCO")
    #ax1[1,i].axvline(1/4.5, linestyle=":", color="k", alpha=0.2)
    #ax1[1,i].axvline(1/5.5, linestyle=":", color="k", alpha=0.2)
    #ax1[1,i].axvline(1/7.5, linestyle=":", color="k", alpha=0.2)
    #ax1[1,i].axvline(1/9.5, linestyle="--", linewidth=1, color="k", alpha=0.2)
    #ax1[1,i].axvline(1/12, linestyle="-.", linewidth=1, color="k", alpha=0.2)
    for j = 1:2
      ax1[i,j].axvline(1/10, linestyle="-", linewidth=0.8, color="k",zorder=0)
      ax1[i,j].axvline(1/17, linestyle="-", linewidth=0.8, color="k",zorder=0)
    end
    #ax1[1,i].axvline(1/18, linestyle="--", linewidth=1, color="k", alpha=0.2)
    #ax1[1,i].axvline(1/26, linestyle=":", color="k", alpha=0.2)
    ax1[i,1].set_ylabel("PSD of \$a_$i~(\\mathrm{mK}^2~\\mathrm{cpd}^{-1})\$")
    #ax1[i,1].set_title("spectra for \$ c^{($i)} \$")
    #ax1[2,i].plot(freqsub[1,idx],freqsub[1,idx].*mean(giixsub,dims=1)[idx], color=clrs[1])
    fgii = 1*mean(giisub,dims=1)[:]
    fgoo = 1*mean(goosub,dims=1)[:]
    fgee = 1*mean(geesub,dims=1)[:]
    ax1[i,2].plot(freqsub[1,idx],1e6fgii[idx], color=clrs[2],label="\$T\$ waves")
    ax1[i,2].plot(freqsub[1,idx],1e6fgoo[idx], color=clrs[3],label="OFES")
    ax1[i,2].plot(freqsub[1,idx],1e6fgee[idx], color=clrs[4],label="ECCO")

    i0 = argmax(freqsub[1,idx].*fgii[idx]) 
    println(1/freqsub[1,idx][i0])
    println(fgii[idx][i0]/fgoo[idx][i0])
    println(fgii[idx][i0]/fgee[idx][i0])
    ie = argmax(freqsub[1,idx].*fgee[idx]) 
    println(1/freqsub[1,idx][ie])
    #ax1[i,2].axvline(1/10, linestyle="-", linewidth=0.8, color="k")
    #ax1[i,2].axvline(1/17, linestyle="-", linewidth=0.8, color="k")
    #ax1[2,i].axvline(1/18, linestyle="--", linewidth=1, color="k", alpha=0.2)
    #ax1[2,i].axvline(1/26, linestyle=":", color="k", alpha=0.2)
    #ax1[i,2].set_xscale("log")
    #ax1[i,2].set_ylabel("\$c_$i\$ variance \$(\\mathrm{mK})^2\$")
    ax1[i,2].set_ylim([0,100+50i])
    ax1[i,2].set_xlim([1/20,0.11])
    ax1[2,i].set_xlabel("frequency (cpd)")
    ax1[1,i].set_title("($(('a':'z')[i]))", loc="left",y=1.0)
    ax1[2,i].set_title("($(('a':'z')[2+i]))", loc="left")
    
    g2io = (abs.(sum(stisub.*conj(stosub),dims=1))).^2 ./ 
           (sum((abs.(stisub)).^2,dims=1).*sum((abs.(stosub)).^2,dims=1))
    dphio = angle.(sum(stisub.*conj(stosub),dims=1))*180/pi 
    dphio[dphio.>180] = dphio[dphio.>180] .- 360
    dphio[dphio.<-180] = dphio[dphio.<-180] .+ 360

    ax2[1,i].plot(freqsub[1,idx],abs.(g2io[idx]),"k")
    ax2[1,i].axhline(cl95, c="k", ls="--", linewidth=.8)
    ax2[1,i].axvline(1/10, linewidth=.8, color="k",zorder=0)
    ax2[1,i].axvline(1/17, linewidth=.8, color="k",zorder=0)
    #ax2[1,i].axvline(1/18, linestyle="--", linewidth=1, color="k", alpha=0.2)
    axt = ax2[1,i].twinx()
    axt.plot(freqsub[1,idx],dphio[idx],"tab:red")
    axt.axhline(0, c="tab:red", ls="--", linewidth=.8)
    ax2[1,i].set_xscale("log")
    ax2[1,i].set_ylabel("coherence\$^2\$")
    ax2[1,i].set_ylim([-1,1])
    ax2[1,i].yaxis.set_ticks(0:0.25:1)
    ax2[1,i].set_title(@sprintf("\$ a_%d \$: \$T\$ waves\$-\$OFES",i))
    axt.set_ylabel("phase",color="tab:red")
    axt.set_ylim([-180,270])
    axt.yaxis.set_ticks(-180:90:180)
    ax2[1,i].yaxis.set_label_coords(-0.15,0.75)
    axt.yaxis.set_label_coords(1.15,0.4)
    ax2[1,i].set_title("($(('a':'z')[i]))", loc="left")
    
    g2ie = (abs.(sum(stisub.*conj(stesub),dims=1))).^2 ./ 
           (sum((abs.(stisub)).^2,dims=1).*sum((abs.(stesub)).^2,dims=1))
    dphie = angle.(sum(stisub.*conj(stesub),dims=1))*180/pi 
    dphie[dphie.>180] = dphie[dphie.>180] .- 360
    dphie[dphie.<-180] = dphie[dphie.<-180] .+ 360

    ax2[2,i].plot(freqsub[1,idx],abs.(g2ie[idx]),"k")
    axt = ax2[2,i].twinx()
    axt.plot(freqsub[1,idx],dphie[idx],"tab:red")
    ax2[2,i].axhline(cl95, c="k", ls="--", linewidth=.8)
    ax2[2,i].axvline(1/10, linewidth=.8, color="k",zorder=0)
    ax2[2,i].axvline(1/17, linewidth=.8, color="k",zorder=0)
    #ax2[2,i].axvline(1/9.5, linestyle="--", linewidth=1, color="k", alpha=0.2)
    #ax2[2,i].axvline(1/12, linestyle="-.", linewidth=1, color="k", alpha=0.2)
    #ax2[2,i].axvline(1/18, linestyle="--", linewidth=1, color="k", alpha=0.2)
    axt.axhline(0, c="tab:red", ls="--", linewidth=0.8)
    ax2[2,i].set_xscale("log")
    ax2[2,i].set_ylabel("coherence\$^2\$")
    #ax2[2,i].set_xlabel("Frequency (day\$^{-1}\$)")
    ax2[2,i].set_ylim([-1,1])
    ax2[2,i].yaxis.set_ticks(0:0.25:1)
    ax2[2,i].set_title(@sprintf("\$ a_%d \$: \$T\$ waves\$-\$ECCO",i))
    axt.set_ylabel("phase",color="tab:red")
    axt.set_ylim([-180,270])
    axt.yaxis.set_ticks(-180:90:180)
    ax2[2,i].yaxis.set_label_coords(-0.15,0.75)
    axt.yaxis.set_label_coords(1.15,0.4)
    ax2[2,i].set_title("($(('a':'z')[2+i]))", loc="left")
    ax2[2,i].set_xlabel("frequency (cpd)")
    
    g2oe = (abs.(sum(stosub.*conj(stesub),dims=1))).^2 ./ 
           (sum((abs.(stosub)).^2,dims=1).*sum((abs.(stesub)).^2,dims=1))
    dphoe = angle.(sum(stosub.*conj(stesub),dims=1))*180/pi 
    dphoe[dphoe.>180] = dphoe[dphoe.>180] .- 360
    dphoe[dphoe.<-180] = dphoe[dphoe.<-180] .+ 360

end
for p in [11,13]
  ax1[1,2].axvline(1/p, linestyle="-", linewidth=0.8, color="k",zorder=3, ymin=0.1,ymax=0.2)
end
for p in [15,13]
  ax1[2,2].axvline(1/p, linestyle="-", linewidth=0.8, color="k",zorder=3, ymin=0.1,ymax=0.2)
end
ax1[2,2].ticklabel_format(axis="y",style="plain")
ax1[1,1].legend(frameon=false,loc=3)
#ax1[1,2].legend(frameon=false)
ax1[2,2].set_xlim([0.05,0.11])
yerr = reshape([1-ci95[1], ci95[2]-1],(2,1))
ax1[1,1].errorbar(2.8e-2, 1, yerr,fmt=".",capsize=2,ecolor="k",mfc="k", mec="k")
ax1[1,1].text(3.2e-2, 0.6, "\$95 \\%\$",c="k")
ax1t = ax1[1,1].twiny()
ax1t.set_xlim(ax1[1,1].get_xlim())
ax1t.set_xscale("log")
new_tick_locations = [1/20, 1/10, 1/5]
ax1t.set_xticks(new_tick_locations)
ax1t.set_xticklabels(new_tick_locations.^-1)
ax1t.set_xlabel("period (d)")
ax1t = ax1[1,2].twiny()
ax1t.set_xlim(ax1[1,2].get_xlim())
new_tick_locations = [1/17, 1/15, 1/13, 1/11,1/10]
ax1t.set_xticks(new_tick_locations)
ax1t.set_xticklabels(new_tick_locations.^-1)
ax1t.set_xlabel("period (d)")
fig1.tight_layout(w_pad=2)
fig1.savefig(@sprintf("results/ofes/spectra_ci_%s.pdf",Kgrid))
fig2.tight_layout(w_pad=2)
fig2.savefig(@sprintf("results/ofes/coherence_ci_%s.pdf",Kgrid))
    
            
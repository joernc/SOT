using PyPlot, Printf, Dates, LinearAlgebra, Statistics, SparseArrays, DataFrames, CSV,FFTW,NCDatasets
using HDF5, Interpolations, Distributions, SeisNoise

PyPlot.matplotlib[:rc]("mathtext",fontset="cm")        #computer modern font 
PyPlot.matplotlib[:rc]("font",family="STIXGeneral")#,size=14)  #font similar to LaTeX
# small font
rc("font", size=8)
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

Kgrid = "Gaussian"
flnm = @sprintf("results/anomalies/nias_tm_H08_%s.h5",Kgrid)
# load time series
t = DateTime(2000, 1, 1) .+ Millisecond.(h5read(flnm, "ti"))
ci = h5read(flnm, "ci")
cofes = h5read(flnm, "cofes")
cecco = h5read(flnm, "cecco")
ec = h5read(flnm, "ec")

ci = ci .- mean(ci,dims=1)
cofes = cofes .- mean(cofes,dims=1)
cecco = cecco .- mean(cecco,dims=1)

# lengths of traces
n = length(t)

# sampling interval
Δ = 1
clrs = ["#a6cee3","tab:blue", "tab:orange", "tab:green"]

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

fig1, ax1 = subplots(3, 2, figsize=(190/25.4, 190/25.4), sharex=true)    
#fig2, ax2 = subplots(3, 2, figsize=(190/25.4, 190/25.4), sharex=true)
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
        goosub[j,:],stosub[j,:],freqsub[j,:] = make_FFT1D(sqrt(8/3)*hann1.*sofes[idxsub[j]],n1,Δ)
        geesub[j,:],stesub[j,:],freqsub[j,:] = make_FFT1D(sqrt(8/3)*hann1.*secco[idxsub[j]],n1,Δ)
    end

    idx = freqsub[1,:] .> 0.02
    fgii = 1*mean(giisub,dims=1)[:]
    fgoo = 1*mean(goosub,dims=1)[:]
    fgee = 1*mean(geesub,dims=1)[:]
    ax1[1,i].plot(freqsub[1,idx],1e6fgii[idx], color=clrs[2],label="\$T\$ waves")
    ax1[1,i].plot(freqsub[1,idx],1e6fgoo[idx], color=clrs[3],label="OFES")
    ax1[1,i].plot(freqsub[1,idx],1e6fgee[idx], color=clrs[4],label="ECCO")
    i0 = argmax(freqsub[1,idx].*fgii[idx]) 
    println(1/freqsub[1,idx][i0])
    println(fgii[idx][i0]/fgoo[idx][i0])
    println(fgii[idx][i0]/fgee[idx][i0])
    ie = argmax(freqsub[1,idx].*fgee[idx]) 
    println(1/freqsub[1,idx][ie])
    for j = 1:3
      ax1[j,i].axvline(1/10, linestyle="-", linewidth=0.8, color="k",zorder=0)
      ax1[j,i].axvline(1/17, linestyle="-", linewidth=0.8, color="k",zorder=0)
    end
    ax1[1,i].set_ylabel("PSD of \$c_$i~(\\mathrm{mK}^2~\\mathrm{cpd}^{-1})\$")
    ax1[1,i].set_title("($(('a':'z')[i]))", loc="left")
    ax1[1,i].set_ylim([0,100+50*i])
    
    g2io = (abs.(sum(stisub.*conj(stosub),dims=1))).^2 ./ 
           (sum((abs.(stisub)).^2,dims=1).*sum((abs.(stosub)).^2,dims=1))
    dphio = angle.(sum(stisub.*conj(stosub),dims=1))*180/pi 
    dphio[dphio.>180] = dphio[dphio.>180] .- 360
    dphio[dphio.<-180] = dphio[dphio.<-180] .+ 360

    ax1[2,i].plot(freqsub[1,idx],abs.(g2io[idx]),"k")
    ax1[2,i].axhline(cl95, c="k", ls="--", linewidth=.8)
    ax1[2,i].axvline(1/10, linewidth=.8, color="k",zorder=0)
    ax1[2,i].axvline(1/17, linewidth=.8, color="k",zorder=0)
    axt = ax1[2,i].twinx()
    axt.plot(freqsub[1,idx],dphio[idx],"tab:red")
    axt.axhline(0, c="tab:red", ls="--", linewidth=.8)
    #ax1[2,i].set_xscale("log")
    ax1[2,i].set_ylabel("coherence\$^2\$")
    ax1[2,i].set_ylim([-1,1])
    ax1[2,i].yaxis.set_ticks(0:0.25:1)
    ax1[2,i].set_title(@sprintf("\$ c_%d \$: \$T\$ waves--OFES",i))
    axt.set_ylabel("phase",color="tab:red")
    axt.set_ylim([-180,270])
    axt.yaxis.set_ticks(-180:90:180)
    ax1[2,i].yaxis.set_label_coords(-0.15,0.75)
    axt.yaxis.set_label_coords(1.15,0.4)
    ax1[2,i].set_title("($(('a':'z')[i+2]))", loc="left")
    
    g2ie = (abs.(sum(stisub.*conj(stesub),dims=1))).^2 ./ 
           (sum((abs.(stisub)).^2,dims=1).*sum((abs.(stesub)).^2,dims=1))
    dphie = angle.(sum(stisub.*conj(stesub),dims=1))*180/pi 
    dphie[dphie.>180] = dphie[dphie.>180] .- 360
    dphie[dphie.<-180] = dphie[dphie.<-180] .+ 360

    ax1[3,i].plot(freqsub[1,idx],abs.(g2ie[idx]),"k")
    axt = ax1[3,i].twinx()
    axt.plot(freqsub[1,idx],dphie[idx],"tab:red")
    ax1[3,i].axhline(cl95, c="k", ls="--", linewidth=.8)
    ax1[3,i].axvline(1/10, linewidth=.8, color="k",zorder=0)
    ax1[3,i].axvline(1/17, linewidth=.8, color="k",zorder=0)
    axt.axhline(0, c="tab:red", ls="--", linewidth=0.8)
    #ax1[3,i].set_xscale("log")
    ax1[3,i].set_ylabel("coherence\$^2\$")
    ax1[3,i].set_ylim([-1,1])
    ax1[3,i].yaxis.set_ticks(0:0.25:1)
    ax1[3,i].set_title(@sprintf("\$ c_%d \$: \$T\$ waves--ECCO",i))
    axt.set_ylabel("phase",color="tab:red")
    axt.set_ylim([-180,270])
    axt.yaxis.set_ticks(-180:90:180)
    ax1[3,i].yaxis.set_label_coords(-0.15,0.75)
    axt.yaxis.set_label_coords(1.15,0.4)
    ax1[3,i].set_title("($(('a':'z')[4+i]))", loc="left")
    ax1[3,i].set_xlabel("frequency (cpd)")
    ax1[3,i].set_xlim([0.05,0.11])

end
for p in [11,13]
  ax1[1,1].axvline(1/p, linestyle="-", linewidth=0.8, color="k",zorder=3, ymin=0.1,ymax=0.2)
end
for p in [15,13]
  ax1[1,2].axvline(1/p, linestyle="-", linewidth=0.8, color="k",zorder=3, ymin=0.1,ymax=0.2)
end
ax1[1,1].legend(frameon=false)
fig1.tight_layout()
fig1.savefig(@sprintf("results/ofes/spectra_ci_%s.pdf",Kgrid))
    
            
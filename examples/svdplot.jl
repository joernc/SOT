using PyPlot
PyPlot.matplotlib[:rc]("mathtext",fontset="cm")        #computer modern font 
PyPlot.matplotlib[:rc]("font",family="STIXGeneral")
#rc("font", size=10)
#rc("axes", titlesize="large")
rc("font", size=8)
rc("axes", titlesize="medium")

λ = [29; 28; 30; 29;; 1.6; 1.5; 1.9; 1.9]
U = [-0.53; -0.59; -0.61;; -0.53; -0.59; -0.61;; -0.48; -0.60; -0.64;; -0.48; -0.60; -0.64;;; -0.81; 0.15; 0.55;; -0.82; 0.15; 0.56;; -0.84; 0.10; 0.54;; -0.84; 0.10; 0.54]
V = [0.85; -0.49; -0.02; 0.16; -0.14;; 0.84; -0.50; -0.03; 0.17; -0.12;; 0.83; -0.49; 0.01; 0.19; -0.18;; 0.83; -0.50; 0.00; 0.20; -0.16;;; -0.09; -0.32; 0.88; -0.32; 0.05;; -0.10; -0.34; 0.87; -0.34; -0.02;; 0.11; -0.20; 0.60; -0.61; 0.47;; 0.11; -0.21; 0.59; -0.67; 0.38]

handles = [plt.Rectangle((0, 0), 1, 1; color, linewidth=0) for color in ["tab:orange", "tab:green"]]

fig, ax = subplots(2, 3; gridspec_kw=Dict("width_ratios"=>[1, 3, 5]), sharex="col", figsize=(8.0, 4.0))
p = ax[1,1].bar(1:4, λ[:,1]; color=repeat(["tab:orange", "tab:green"], 2))#, alpha=[1.0, 1.0, 0.8, 0.8])
q = ax[2,1].bar(1:4, λ[:,2]; color=repeat(["tab:orange", "tab:green"], 2))
for i = 3:4
  p[i].set_alpha(0.5)
  q[i].set_alpha(0.5)
end
p = ax[1,2].bar([1:4; 6:9; 11:14], reshape(U[:,:,1]', 12); color=repeat(["tab:orange", "tab:green"], 6))
q = ax[2,2].bar([1:4; 6:9; 11:14], reshape(U[:,:,2]', 12); color=repeat(["tab:orange", "tab:green"], 6))
for i = [3, 4, 7, 8, 11, 12]
  p[i].set_alpha(0.5)
  q[i].set_alpha(0.5)
end
p = ax[1,3].bar([1:4; 6:9; 11:14; 16:19; 21:24], reshape(V[:,:,1]', 20); color=repeat(["tab:orange", "tab:green"], 10))
q = ax[2,3].bar([1:4; 6:9; 11:14; 16:19; 21:24], reshape(V[:,:,2]', 20); color=repeat(["tab:orange", "tab:green"], 10))
for i = [3, 4, 7, 8, 11, 12, 15, 16, 19, 20]
  p[i].set_alpha(0.5)
  q[i].set_alpha(0.5)
end
for i = 1:6
  ax[i].spines["right"].set_visible(false)
  ax[i].spines["top"].set_visible(false)
end
for i = 1:2, j = 2:3
  ax[i,j].spines["bottom"].set_visible(false)
  ax[i,j].axhline(0; zorder=2, color="black", linewidth=0.8)
  ax[i,j].set_yticks(-1:1)
end
for i = 1:2, j = 1:3
  l = ('a':'z')[(i-1)*3+j]
  ax[i,j].set_title("($l)"; loc="left")
end
ax[1,1].set_xlim(0, 5)
ax[1,2].set_xlim(0, 15)
ax[1,3].set_xlim(0, 25)
ax[1,2].set_ylim(-1, 1)
ax[2,2].set_ylim(-1, 1)
ax[1,3].set_ylim(-1, 1)
ax[2,3].set_ylim(-1, 1)
ax[2,1].set_xticks([])
ax[2,2].set_xticks(2.5:5:12.5, 1:3)
ax[2,3].set_xticks(2.5:5:22.5, 1:5)
ax[1,1].set_ylabel(L"$λ_1$ (s K$^{-1}$)")
ax[2,1].set_ylabel(L"$λ_2$ (s K$^{-1}$)")
ax[1,2].set_title("left singular vector")
ax[2,2].set_title("left singular vector")
ax[1,3].set_title("right singular vector")
ax[2,3].set_title("right singular vector")
ax[2,2].set_xlabel(L"frequency $j$")
ax[2,3].set_xlabel(L"baroclinic mode $n$")
ax[1,3].legend(handles, ["OFES", "ECCO"]; frameon=false, ncol=2)
fig.align_ylabels()
fig.tight_layout()
fig.savefig("results/ofes/svd.pdf")

"""Same as cutjust_plot_pi0match_dist.py but with the raw m_gg invariant mass of the selected
pair (not the |m_gg-m_pi0| difference), mirroring the fig1 -> fig3 relation. One panel per gen
type (1/2/3), pi0-window cut shaded. Uses cutjust_pi0match_tables.npz (with best_m saved)."""
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.facecolor":"white","axes.facecolor":"white","axes.grid":True,
                      "grid.alpha":0.3,"font.size":10})
MPI0=0.13498; WIN=0.05
colors={"same_pi0":"#1baf7a","cross_pi0":"#eda100","to_spurious":"#eb6834"}
labels={"same_pi0":"correct: same $\\pi^0$","cross_pi0":"wrong: other $\\pi^0$ of same tau",
        "to_spurious":"wrong: non-$\\pi^0$ photon"}

d=np.load("cutjust_pi0match_tables.npz")
gtype,truepi0,has_partner,same_pi0,cross_pi0,to_spurious,best_m=[d[k] for k in
    ["gtype","truepi0","has_partner","same_pi0","cross_pi0","to_spurious","best_m"]]

cats=[(1,"gen 1 ($\\pi\\pi^0$, 1 $\\pi^0$)"),(2,"gen 2 ($\\pi 2\\pi^0$, 2 $\\pi^0$'s)"),(3,"gen 3 ($3\\pi^0$+, $\\geq$3 $\\pi^0$'s)")]
bins=np.logspace(-3,0.7,50)

fig,axs=plt.subplots(1,3,figsize=(14,4.6),sharey=True)
for ax,(g,lab) in zip(axs,cats):
    base=(gtype==g)&truepi0&has_partner
    ax.axvspan(MPI0-WIN,MPI0+WIN,color="gray",alpha=0.15,zorder=0,
               label=f"pi0 window (±{WIN*1000:.0f} MeV)" if ax is axs[0] else None)
    for key in ["same_pi0","cross_pi0","to_spurious"]:
        sel=base&d[key]
        if sel.sum()==0: continue
        x=best_m[sel]; x=x[np.isfinite(x)&(x>0)]
        ax.hist(x,bins=bins,histtype="step",density=True,lw=1.8,color=colors[key],
                 label=f"{labels[key]} (N={sel.sum()})")
    ax.axvline(MPI0,color="k",ls="--",lw=1.2)
    ax.set_xscale("log")
    ax.set_xlabel(r"$m_{\gamma\gamma}$ of the selected pair [GeV]")
    ax.set_title(lab,fontsize=10)
    ax.legend(fontsize=7.3,loc="upper left",framealpha=0.6,facecolor="lightgray")
axs[0].set_ylabel("true $\\pi^0$ photons (normalized)")
fig.suptitle(r"$m_{\gamma\gamma}$ of the selected pair relative to the $\pi^0$-window cut, by pairing outcome",fontsize=10.5)
fig.tight_layout(rect=[0,0,1,0.94]); fig.savefig("cutjust_fig7_pi0match_mass.png",dpi=140); plt.close(fig)
print("done")

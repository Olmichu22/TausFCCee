"""Repeats the bottom-right panel of figs/discrim_fig1_distribuciones.png (min |m_gg - m_pi0|
with the best partner in the cone) but without subtracting m_pi0: the raw m_gg invariant mass,
for the same 4 populations. Uses cutjust_mgg_tables.npz (cutjust_extract_mgg.py)."""
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.facecolor":"white","axes.facecolor":"white","axes.grid":True,
                      "grid.alpha":0.3,"font.size":10})
C=["#2a78d6","#eb6834","#1baf7a","#eda100"]
MPI0=0.13498

d=np.load("cutjust_mgg_tables.npz")
gtype,rtype,truepi0,mgg=[d[k] for k in ["gtype","rtype","truepi0","mgg_best"]]

g1=(gtype==1)&np.isin(rtype,[1,2])
g2=(gtype==2)
spur=np.isin(gtype,[1,2,3])&~truepi0
# "bad: gen0->reco1" omitted: 100% of these photons have no partner (see cutjust_fig1),
# so there is nothing to histogram for that category here.

pops=[("good: gen1->reco1/2",g1,C[0]),("good: gen2",g2,C[2]),("spurious in gen1,2,3",spur,C[3])]

fig,ax=plt.subplots(figsize=(6.4,4.8))
bins=np.logspace(-3,0.7,50)
WIN=0.05
ax.axvspan(MPI0-WIN,MPI0+WIN,color="gray",alpha=0.15,zorder=0,label=f"pi0 window (±{WIN*1000:.0f} MeV)")
for lab,sel,c in pops:
    x=mgg[sel]; x=x[np.isfinite(x)&(x>0)]
    nopart=np.sum(sel&~np.isfinite(mgg))
    ax.hist(x,bins=bins,histtype="step",density=True,color=c,lw=1.8,
            label=f"{lab} (N={sel.sum()}, {nopart/sel.sum():.0%} no partner)")
ax.axvline(MPI0,color="k",ls="--",lw=1.2)
ax.text(MPI0*1.05,ax.get_ylim()[1]*0.99,f"$m_{{\\pi^0}}$ = {MPI0:.3f} GeV",rotation=90,va="top",ha="left",fontsize=8.5,
         bbox=dict(facecolor="lightgray",alpha=0.6,edgecolor="none",pad=1.5))
ax.set_xscale("log")
ax.set_xlabel(r"$m_{\gamma\gamma}$ with the cone partner that minimizes $|m_{\gamma\gamma}-m_{\pi^0}|$ [GeV]")
ax.set_ylabel("photons (normalized)")
ax.set_title(r"$m_{\gamma\gamma}$ of the best partner",fontsize=10)
ax.legend(fontsize=7.2,loc="upper left",framealpha=0.6,facecolor="lightgray",handlelength=1.5,
          labelspacing=0.3,borderpad=0.5)
fig.tight_layout(); fig.savefig("cutjust_fig3_mgg_raw.png",dpi=140); plt.close(fig)
print("done")

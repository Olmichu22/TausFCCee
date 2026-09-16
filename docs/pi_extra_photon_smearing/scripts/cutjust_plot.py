"""Two figures justifying the final recommended photon-rejection cut (REPORT.md, Section 4):
reject an orphan photon (no pi0 partner in the cone) if P_gamma/P_pion < 0.05, or if
P_gamma > 2 GeV and m(pion+gamma) > 1.2 GeV. Uses cutjust_tables.npz (cutjust_extract.py)."""
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.facecolor":"white","axes.facecolor":"white","axes.grid":True,
                      "grid.alpha":0.3,"font.size":10})
C_BAD="#eb6834"; C_TRUE="#2a78d6"; C_CUT="#1a1a1a"

d=np.load("cutjust_tables.npz")
gtype,rtype,phP,Ppi,mpig,dmpi0=[d[k] for k in ["gtype","rtype","phP","Ppi","mpig","dmpi0"]]
orphan=~np.isfinite(dmpi0)|(dmpi0>0.05)

bad=(gtype==0)&(rtype==1)                 # extra photon: tau->pi nu migrating to reco 1 (the pathology the cut targets)
true=np.isin(gtype,[1,2,3])&(dmpi0<=0.05)&np.isfinite(dmpi0)  # genuine pi0 photon (has a pi0 partner)
# background the recommended cut must spare: orphan photons from real pi0-bearing decays (partner lost)
truth_orphan=np.isin(gtype,[1,2,3])&orphan

fpi=phP/Ppi

# ------------------------------------------------------------------ Figure 1 (1D): P_gamma/P_pion
fig,ax=plt.subplots(figsize=(6.2,4.6))
bins=np.logspace(-3,1.5,60)
for lab,sel,c in [("orphan photon, gen 1/2/3 (true $\\pi^0$, partner lost — must be kept)",truth_orphan,C_TRUE),
                   ("extra photon, gen 0 $\\to$ reco 1 (FSR / shower fragment — should be rejected)",bad,C_BAD)]:
    x=fpi[sel]; x=x[np.isfinite(x)&(x>0)]
    ax.hist(x,bins=bins,density=True,histtype="step",lw=2.0,color=c,label=f"{lab}  (N={sel.sum()})")
ax.axvline(0.05,color=C_CUT,ls="--",lw=1.5)
ax.text(0.055,ax.get_ylim()[1]*0.02,"cut: $P_\\gamma/P_\\pi$ < 0.05",rotation=90,va="bottom",ha="left",fontsize=8.5)
ax.set_xscale("log")
ax.set_xlabel(r"$P_\gamma / P_\pi$")
ax.set_ylabel("photons (normalized)")
ax.set_title("Soft-photon cut: orphan photons only",fontsize=10.5)
ax.legend(fontsize=7.8,loc="upper left")
fig.tight_layout(); fig.savefig("cutjust_fig1_soft_ratio_1d.png",dpi=140); plt.close(fig)

# ------------------------------------------------------------------ Figure 2 (2D): P_gamma vs m(pi+gamma)
fig,axs=plt.subplots(1,2,figsize=(11.5,4.8),sharex=True,sharey=True)
xbins=np.linspace(0,3,61); ybins=np.logspace(-1,1.3,61)
pops=[("orphan photon, gen 1/2/3\n(true $\\pi^0$, partner lost)",truth_orphan,axs[0]),
      ("extra photon, gen 0 $\\to$ reco 1\n(FSR of the tau line / shower fragment)",bad,axs[1])]
for lab,sel,ax in pops:
    h=ax.hist2d(mpig[sel&orphan],phP[sel&orphan],bins=[xbins,ybins],cmap="viridis",cmin=1)
    ax.set_yscale("log")
    ax.axvline(1.2,color="white",ls="--",lw=1.3,zorder=5)
    ax.axhline(2.0,color="white",ls="--",lw=1.3,zorder=5)
    ax.fill_betweenx([2.0,ybins[-1]],1.2,xbins[-1],facecolor="none",edgecolor="white",
                      hatch="////",lw=0,zorder=4)
    ax.text(1.25,ybins[-1]*0.75,"rejected",color="white",fontsize=8.5,ha="left",va="top",zorder=6)
    ax.set_xlabel(r"$m(\pi+\gamma)$ [GeV]")
    ax.set_title(f"{lab}\n(orphan photons, N={sel.sum()})",fontsize=9.5)
    fig.colorbar(h[3],ax=ax,shrink=0.85,pad=0.02)
axs[0].set_ylabel(r"$P_\gamma$ [GeV]")
fig.suptitle(r"Hard-photon cut: reject if $P_\gamma$ > 2 GeV and $m(\pi+\gamma)$ > 1.2 GeV (dashed lines; shaded = rejected)",fontsize=10.5)
fig.tight_layout(rect=[0,0,1,0.93]); fig.savefig("cutjust_fig2_hard_mass_2d.png",dpi=140); plt.close(fig)
print("done")

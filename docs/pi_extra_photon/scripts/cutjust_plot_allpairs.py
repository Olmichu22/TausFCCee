"""Figure 2: compares, per tau category (gen1->reco1/2, gen2), the m_gg of ALL gamma-gamma
combinations in the cone (combinatorial background) against the m_gg of ONLY the selected
(minimum |m_gg-m_pi0|) pair per tau. Shows whether picking the minimum actually pulls a clean
pi0 peak out of the combinatorics. Uses cutjust_allpairs_tables.npz (cutjust_extract_allpairs.py)."""
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.facecolor":"white","axes.facecolor":"white","axes.grid":True,
                      "grid.alpha":0.3,"font.size":10})
MPI0=0.13498; WIN=0.05
C_ALL="#b0b0b0"; C_BEST="#2a78d6"

d=np.load("cutjust_allpairs_tables.npz")
pg,pr,pm=[d[k] for k in ["pair_gtype","pair_rtype","pair_mgg"]]
tg,tr,tb=[d[k] for k in ["tau_gtype","tau_rtype","tau_best_mgg"]]

cats=[("gen1 -> reco1/2",(pg==1)&np.isin(pr,[1,2]),(tg==1)&np.isin(tr,[1,2])),
      ("gen2",(pg==2),(tg==2))]

fig,axs=plt.subplots(1,2,figsize=(11,4.6),sharey=False)
bins=np.logspace(-3,0.7,50)
for ax,(lab,psel,tsel) in zip(axs,cats):
    xall=pm[psel]; xall=xall[np.isfinite(xall)&(xall>0)]
    xbest=tb[tsel]; xbest=xbest[np.isfinite(xbest)&(xbest>0)]
    ax.axvspan(MPI0-WIN,MPI0+WIN,color="gray",alpha=0.15,zorder=0)
    ax.hist(xall,bins=bins,density=True,color=C_ALL,alpha=0.6,label=f"all $\\gamma\\gamma$ pairs (N={len(xall)})")
    ax.hist(xbest,bins=bins,histtype="step",density=True,color=C_BEST,lw=2.0,
            label=f"selected pair only (N={len(xbest)})")
    ax.axvline(MPI0,color="k",ls="--",lw=1.2)
    ax.set_xscale("log")
    ax.set_xlabel(r"$m_{\gamma\gamma}$ [GeV]")
    ax.set_title(lab,fontsize=10.5)
    ax.legend(fontsize=8,loc="upper left",framealpha=0.6,facecolor="lightgray")
axs[0].set_ylabel("pairs / taus (normalized)")
fig.suptitle(r"Does the minimum $|m_{\gamma\gamma}-m_{\pi^0}|$ selection pull a clean $\pi^0$ peak out of the combinatorics?",fontsize=10.5)
fig.tight_layout(rect=[0,0,1,0.94]); fig.savefig("cutjust_fig4_allpairs_vs_best.png",dpi=140); plt.close(fig)
print("done")

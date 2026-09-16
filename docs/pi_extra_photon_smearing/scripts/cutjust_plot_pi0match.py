"""For real pi0 photons (truepi0), classifies what the 'minimum |m_gg-m_pi0|' selection
actually picked as partner: the true sibling photon of the same pi0 (same_pi0, correct),
a photon from a DIFFERENT pi0 of the same tau (cross_pi0, wrong pairing but both photons are
still genuine pi0 photons), a non-pi0 photon (to_spurious, wrong), or nothing (no_partner).
Uses cutjust_pi0match_tables.npz (cutjust_extract_pi0match.py)."""
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.facecolor":"white","axes.facecolor":"white","axes.grid":True,
                      "grid.alpha":0.3,"axes.axisbelow":True,"font.size":10})

d=np.load("cutjust_pi0match_tables.npz")
gtype,rtype,truepi0,has_partner,same_pi0,cross_pi0,to_spurious=[d[k] for k in
    ["gtype","rtype","truepi0","has_partner","same_pi0","cross_pi0","to_spurious"]]

cats=[(1,"gen 1\n($\\pi\\pi^0$, 1 $\\pi^0$)"),(2,"gen 2\n($\\pi 2\\pi^0$, 2 $\\pi^0$'s)"),(3,"gen 3\n($3\\pi^0$+, $\\geq$3 $\\pi^0$'s)")]
labels=["correct: same $\\pi^0$","wrong: other $\\pi^0$ of same tau","wrong: non-$\\pi^0$ photon","no partner found"]
colors=["#1baf7a","#eda100","#eb6834","#6b6b6b"]

fracs=[]; Ns=[]
for g,_ in cats:
    sel=(gtype==g)&truepi0
    N=sel.sum(); Ns.append(N)
    fracs.append([np.mean(same_pi0[sel]),np.mean(cross_pi0[sel]),np.mean(to_spurious[sel]),np.mean(~has_partner[sel])])
fracs=np.array(fracs)

fig,ax=plt.subplots(figsize=(6.6,4.8))
x=np.arange(len(cats)); bottom=np.zeros(len(cats))
for k in range(4):
    ax.bar(x,fracs[:,k],bottom=bottom,color=colors[k],label=labels[k],width=0.55,edgecolor="white",lw=0.8)
    for xi,(f,b) in enumerate(zip(fracs[:,k],bottom)):
        if f>0.03: ax.text(xi,b+f/2,f"{f:.0%}",ha="center",va="center",fontsize=8.5,
                             color="white" if k in (0,2) else "black")
    bottom+=fracs[:,k]
ax.set_xticks(x); ax.set_xticklabels([f"{lab}\n(N={n})" for (g,lab),n in zip(cats,Ns)],fontsize=9)
ax.set_ylabel("fraction of true $\\pi^0$ photons")
ax.set_ylim(0,1.02)
ax.set_title("Does the 'minimum $|m_{\\gamma\\gamma}-m_{\\pi^0}|$' selection pick the true sibling photon?",fontsize=10)
ax.legend(fontsize=8,loc="upper left",bbox_to_anchor=(1.01,1.02),framealpha=0.6,facecolor="lightgray")
fig.tight_layout(); fig.savefig("cutjust_fig5_pi0_pairing_purity.png",dpi=140); plt.close(fig)
print("done")

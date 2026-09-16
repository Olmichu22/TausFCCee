"""Discriminacion del foton extra en tau->pi nu: distribuciones, ROC, puntos de trabajo,
matrices de migracion y eficiencias vs P. Produce figuras png y RESULTS.md."""
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, collections
from discrim_common import *
C=["#2a78d6","#eb6834","#1baf7a","#eda100","#e87ba4","#008300","#7a5cc7","#6b6b6b"]
plt.rcParams.update({"figure.facecolor":"white","axes.facecolor":"white","axes.grid":True,"grid.alpha":0.3,"font.size":9})
t,p,A=load()
NT=len(t["gtype"]); NP=len(p["P"])
md=[]  # lineas del RESULTS.md
def W(s=""): md.append(s)

# ------------------------------------------------------------------ categorias de fotones
bad=(p["gtype"]==0)&(p["rtype"]==1)
g1=(p["gtype"]==1)&np.isin(p["rtype"],[1,2])
g2=(p["gtype"]==2)
spur=np.isin(p["gtype"],[1,2,3])&~p["truepi0"]          # fotones espurios dentro de taus con pi0
true=np.isin(p["gtype"],[1,2,3])&p["truepi0"]            # fotones verdaderos de pi0 del propio tau
W("# Discriminacion del foton extra en tau -> pi nu (gen 0 -> reco 1)")
W(); W(f"Muestra: fichero completo (2M eventos, {len(A['gtype'])} gen taus, {NT} con reco emparejado, {NP} fotones constituyentes).")
W(); W("## 1. Poblaciones de fotones")
W(); W("| Poblacion | N fotones | frac. foton de pi0 del propio tau |"); W("|---|---|---|")
for n,s in [("malo: gen 0 -> reco 1",bad),("bueno: gen 1 -> reco 1/2",g1),("bueno: gen 2 (todos)",g2),
            ("verdadero pi0 (gen 1,2,3)",true),("espurio en gen 1,2,3",spur)]:
    W(f"| {n} | {s.sum()} | {p['truepi0'][s].mean():.3f} |")
W(); W("Origen gen del foton malo (RecoPhotonGenMatchIdx -> GenPhotonOrigin; -9 = sin enlace a foton gen):")
W(); W("| origen | N | frac | P mediana [GeV] | dR mediana |"); W("|---|---|---|---|---|")
names={1:"FSR del tau",-9:"sin match gen (fragmento/cluster partido)",2:"radiacion de cargada",0:"pi0 (otro tau)",3:"otro",-99:"no en bloque RecoPhoton"}
names_en={1:"tau FSR",-9:"no gen match (split fragment/cluster)",2:"charged radiation",0:"pi0 (other tau)",3:"other",-99:"not in RecoPhoton block"}
for o,n in collections.Counter(p["origin"][bad].tolist()).most_common():
    s=bad&(p["origin"]==o); W(f"| {names.get(o,o)} | {n} | {n/bad.sum():.3f} | {np.median(p['P'][s]):.2f} | {np.median(p['dR'][s]):.3f} |")

# ------------------------------------------------------------------ (a) distribuciones
vars_=[("P","photon P [GeV]",np.logspace(-1.3,1.7,45),True),("dR","dR(photon, leading pion) [rad]",np.linspace(0,0.4,41),False),
       ("fpi","P_gamma / P_pion",np.logspace(-2.5,1.3,45),True),("ftau","P_gamma / P_tau(reco)",np.logspace(-2.5,0,45),True),
       ("mpig","m(pion + photon) [GeV]",np.linspace(0,2,51),False),("dmpi0","min |m_gg - m_pi0| with another photon in the cone [GeV]",np.logspace(-3,0.5,45),True)]
pops=[("bad: gen0->reco1",bad,C[1]),("good: gen1->reco1/2",g1,C[0]),("good: gen2",g2,C[2]),("spurious in gen1,2,3",spur,C[3])]
fig,axs=plt.subplots(2,3,figsize=(13,7.5))
for ax,(v,lab,bins,logx) in zip(axs.flat,vars_):
    for n,s,c in pops:
        x=p[v][s]; x=x[np.isfinite(x)]
        ax.hist(x,bins=bins,histtype="step",density=True,color=c,lw=1.6,label=n if ax is axs.flat[0] else None)
    if logx: ax.set_xscale("log")
    ax.set_xlabel(lab); ax.set_ylabel("density (norm. to 1)")
    if v=="dmpi0": ax.axvline(0.05,color="k",ls="--",lw=1); ax.text(0.055,ax.get_ylim()[1]*0.8,"50 MeV window",fontsize=8)
axs.flat[0].legend(fontsize=8); fig.suptitle("Per-photon variables: extra photon (gen 0 -> reco 1) vs pi0 photons")
fig.tight_layout(); fig.savefig("fig1_distribuciones.png",dpi=130); plt.close(fig)

# distribucion 2D-ish del foton malo por origen
fig,axs=plt.subplots(1,3,figsize=(13,3.8))
for (o,c) in [(1,C[0]),(-9,C[1]),(2,C[2])]:
    s=bad&(p["origin"]==o)
    axs[0].hist(p["P"][s],bins=np.logspace(-1.3,1.7,40),histtype="step",density=True,color=c,lw=1.6,label=names_en[o])
    axs[1].hist(p["dR"][s],bins=np.linspace(0,0.4,41),histtype="step",density=True,color=c,lw=1.6)
    axs[2].hist(p["mpig"][s],bins=np.linspace(0,2,51),histtype="step",density=True,color=c,lw=1.6)
axs[0].set_xscale("log"); axs[0].set_xlabel("photon P [GeV]"); axs[1].set_xlabel("dR(photon, pion) [rad]"); axs[2].set_xlabel("m(pion+photon) [GeV]")
for ax in axs: ax.set_ylabel("density")
axs[0].legend(fontsize=8); fig.suptitle("Extra photon in gen 0 -> reco 1, by photon origin")
fig.tight_layout(); fig.savefig("fig2_malo_por_origen.png",dpi=130); plt.close(fig)

# ------------------------------------------------------------------ (b) ROC a nivel de foton
W(); W("## 2. ROC a nivel de foton")
W(); W("Senal = foton malo (gen 0 -> reco 1); fondo a preservar = foton verdadero de pi0 del propio tau (gen 1,2,3).")
W("Se corta 'rechazar si var < x' (P, P/Ppi, P/Ptau, m(pi+gamma)) o 'var > x' (dR, |m_gg-m_pi0|).")
fig,axs=plt.subplots(1,2,figsize=(11,4.5))
def roc(v,sig,bkg,upper):
    xs=np.quantile(p[v][sig][np.isfinite(p[v][sig])],np.linspace(0.005,0.995,150))
    e_s=np.array([(p[v][sig]<x).mean() if upper else (p[v][sig]>x).mean() for x in xs])
    e_b=np.array([(p[v][bkg]<x).mean() if upper else (p[v][bkg]>x).mean() for x in xs])
    return xs,e_s,e_b
W(); W("| variable | rechazo malo al perder 1% verdaderos | al perder 2% | al perder 5% | corte (1%) |"); W("|---|---|---|---|---|")
for i,(v,lab,up) in enumerate([("P","P",True),("fpi","P/Ppi",True),("ftau","P/Ptau",True),("mpig","m(pi+g)",True),("dR","dR",False)]):
    xs,es,eb=roc(v,bad,true,up); axs[0].plot(eb,es,color=C[i],lw=1.8,label=lab)
    if not up: xs,es,eb=xs[::-1],es[::-1],eb[::-1]
    f=lambda l: es[np.searchsorted(eb,l)] if eb[-1]>=l else np.nan
    W(f"| {lab} (todos los fotones) | {f(0.01):.3f} | {f(0.02):.3f} | {f(0.05):.3f} | {xs[min(np.searchsorted(eb,0.01),len(xs)-1)]:.3g} |")
nop=~p["haspartner"]|(p["dmpi0"]>0.05)
for i,(v,lab,up) in enumerate([("P","P",True),("fpi","P/Ppi",True),("ftau","P/Ptau",True),("mpig","m(pi+g)",True),("dR","dR",False)]):
    # rechazo solo si ademas no tiene pareja pi0: la perdida se mide sobre todos los verdaderos
    xs=np.quantile(p[v][bad],np.linspace(0.005,0.995,150))
    es=np.array([((p[v][bad]<x)&nop[bad]).mean() if up else ((p[v][bad]>x)&nop[bad]).mean() for x in xs])
    eb=np.array([((p[v][true]<x)&nop[true]).mean() if up else ((p[v][true]>x)&nop[true]).mean() for x in xs])
    axs[1].plot(eb,es,color=C[i],lw=1.8,label=lab+" & no pi0 partner")
    if not up: xs,es,eb=xs[::-1],es[::-1],eb[::-1]
    f=lambda l: es[np.searchsorted(eb,l)] if eb[-1]>=l else np.nan
    W(f"| {lab} & sin pareja pi0 | {f(0.01):.3f} | {f(0.02):.3f} | {f(0.05):.3f} | {xs[min(np.searchsorted(eb,0.01),len(xs)-1)]:.3g} |")
for ax,tt in zip(axs,["simple cut on the variable","cut only if the photon does not form m_gg ~ m_pi0 (|dm|>50 MeV)"]):
    ax.set_xlabel("fraction of true pi0 photons rejected"); ax.set_ylabel("fraction of bad photons rejected"); ax.set_xlim(0,0.3); ax.set_ylim(0,1); ax.legend(fontsize=8); ax.set_title(tt,fontsize=9)
fig.suptitle("Per-photon ROC: extra photon (gen 0 -> reco 1) vs pi0 photons (gen 1,2,3)")
fig.tight_layout(); fig.savefig("fig3_roc_foton.png",dpi=130); plt.close(fig)

# ------------------------------------------------------------------ eficiencias a nivel de tau
GENS=[0,1,2,3,10,11,12]
Ntot={g:np.sum(A["gtype"]==g) for g in GENS}
def eff_all(dm,lo=None,hi=None):
    out={}
    for g in GENS:
        if lo is None: N=Ntot[g]; s=(t["gtype"]==g)&(dm==g)
        else: N=np.sum((A["gtype"]==g)&(A["gvisP"]>=lo)&(A["gvisP"]<hi)); s=(t["gtype"]==g)&(dm==g)&(t["gvisP"]>=lo)&(t["gvisP"]<hi)
        out[g]=s.sum()/N
    return out
base=eff_all(t["rdm"]); base2040=eff_all(t["rdm"],20,40)

# ROC a nivel de tau por familias
fams={"P < x":lambda x:p["P"]<x,"P/Ppi < x":lambda x:p["fpi"]<x,"P/Ptau < x":lambda x:p["ftau"]<x,
      "1 foton en cono & P < x":lambda x:(p["nph"]==1)&(p["P"]<x),
      "sin pareja pi0 & P < x":lambda x:(p["dmpi0"]>0.05)&(p["P"]<x),
      "sin pareja pi0 & P/Ppi < x":lambda x:(p["dmpi0"]>0.05)&(p["fpi"]<x),
      "sin pareja pi0 & P/Ptau < x":lambda x:(p["dmpi0"]>0.05)&(p["ftau"]<x),
      "sin pareja pi0 & P < x*sqrt(Ppi)":lambda x:(p["dmpi0"]>0.05)&(p["P"]<x*np.sqrt(p["Ppi"]))}
grids={"P < x":[0.2,0.3,0.5,0.7,1,1.5,2,3,5],"P/Ppi < x":[0.01,0.02,0.03,0.05,0.07,0.1,0.15,0.2,0.3],"P/Ptau < x":[0.01,0.02,0.03,0.05,0.07,0.1,0.15,0.2],
       "1 foton en cono & P < x":[0.3,0.5,1,1.5,2,3,5,10],"sin pareja pi0 & P < x":[0.3,0.5,0.7,1,1.5,2,3,5],
       "sin pareja pi0 & P/Ppi < x":[0.01,0.02,0.03,0.05,0.07,0.1,0.15,0.2,0.3],"sin pareja pi0 & P/Ptau < x":[0.01,0.02,0.03,0.05,0.07,0.1,0.15,0.2],
       "sin pareja pi0 & P < x*sqrt(Ppi)":[0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.7]}
fams_en={"P < x":"P < x","P/Ppi < x":"P/Ppi < x","P/Ptau < x":"P/Ptau < x",
      "1 foton en cono & P < x":"1 photon in cone & P < x",
      "sin pareja pi0 & P < x":"no pi0 partner & P < x",
      "sin pareja pi0 & P/Ppi < x":"no pi0 partner & P/Ppi < x",
      "sin pareja pi0 & P/Ptau < x":"no pi0 partner & P/Ptau < x",
      "sin pareja pi0 & P < x*sqrt(Ppi)":"no pi0 partner & P < x*sqrt(Ppi)"}
n01=np.sum((t["gtype"]==0)&(t["rtype"]==1))
fig,axs=plt.subplots(1,3,figsize=(14,4.6))
W(); W("## 3. Barrido a nivel de tau (fichero completo)")
W(); W("Delta de eficiencia gen X -> DM X (denominador: todos los gen taus de ese tipo, incluidos los no emparejados) respecto al baseline.")
W(f"Baseline: " + ", ".join(f"e{g}={base[g]:.4f}" for g in GENS) + f". gen 0 -> reco 1: {n01} taus ({n01/Ntot[0]:.4f} de gen 0).")
W(); W("| criterio | x | de0 | de1 | de2 | de3 | de10 | de11 | de12 | gen0->reco1 restantes |"); W("|---|---|---|---|---|---|---|---|---|---|")
for i,(name,f) in enumerate(fams.items()):
    pts=[]
    for x in grids[name]:
        rt,dm,M,P,n=emulate(t,p,~f(x)); e=eff_all(dm); rem=np.sum((t["gtype"]==0)&(rt==1))
        pts.append([e[g]-base[g] for g in GENS]); W(f"| {name} | {x} | "+" | ".join(f"{e[g]-base[g]:+.4f}" for g in GENS)+f" | {rem} ({rem/n01:.2f}) |")
    pts=np.array(pts)
    for ax,j,gl in zip(axs,[1,2,3],["gen 1 -> DM 1","gen 2 -> DM 2","gen 3 -> DM 3"]):
        ax.plot(-pts[:,j],pts[:,0],"-o",ms=4,color=C[i],lw=1.6,label=fams_en[name])
        ax.set_xlabel(f"efficiency loss {gl}"); ax.set_ylabel("efficiency gain gen 0 -> DM 0"); ax.set_xlim(-0.012,0.08); ax.set_ylim(0,0.035); ax.axvline(0,color="k",lw=0.8)
axs[0].legend(fontsize=7); fig.suptitle("Tau-level cost-benefit of each criteria family (full file)")
fig.tight_layout(); fig.savefig("fig4_roc_tau.png",dpi=130); plt.close(fig)

# ------------------------------------------------------------------ (c) puntos de trabajo
WPS=collections.OrderedDict([
 ("baseline","(sin corte)"),
 ("P<0.5 GeV",lambda:p["P"]<0.5),("P<1 GeV",lambda:p["P"]<1.0),("P<2 GeV",lambda:p["P"]<2.0),
 ("WP-A: sin pareja pi0 & P<1 GeV",lambda:(p["dmpi0"]>0.05)&(p["P"]<1.0)),
 ("WP-B: sin pareja pi0 & P<1.5 GeV",lambda:(p["dmpi0"]>0.05)&(p["P"]<1.5)),
 ("WP-C: sin pareja pi0 & P/Ppi<0.05",lambda:(p["dmpi0"]>0.05)&(p["fpi"]<0.05)),
 ("WP-D: sin pareja pi0 & P<0.3*sqrt(Ppi)",lambda:(p["dmpi0"]>0.05)&(p["P"]<0.3*np.sqrt(p["Ppi"]))),
 ("WP-E: P<0.3 GeV | (sin pareja pi0 & P<1 GeV)",lambda:(p["P"]<0.3)|((p["dmpi0"]>0.05)&(p["P"]<1.0))),
])
WPS_en={"baseline":"baseline","P<0.5 GeV":"P<0.5 GeV","P<1 GeV":"P<1 GeV","P<2 GeV":"P<2 GeV",
 "WP-A: sin pareja pi0 & P<1 GeV":"WP-A: no pi0 partner & P<1 GeV",
 "WP-B: sin pareja pi0 & P<1.5 GeV":"WP-B: no pi0 partner & P<1.5 GeV",
 "WP-C: sin pareja pi0 & P/Ppi<0.05":"WP-C: no pi0 partner & P/Ppi<0.05",
 "WP-D: sin pareja pi0 & P<0.3*sqrt(Ppi)":"WP-D: no pi0 partner & P<0.3*sqrt(Ppi)",
 "WP-E: P<0.3 GeV | (sin pareja pi0 & P<1 GeV)":"WP-E: P<0.3 GeV | (no pi0 partner & P<1 GeV)"}
DMS=[0,1,2,3,10,11,12,-11,-13,-20]; GROWS=[0,1,2,3,10,11,12,-11,-13]
def matrix(dm):
    Mx=np.zeros((len(GROWS),len(DMS)+2))
    for i,g in enumerate(GROWS):
        N=np.sum(A["gtype"]==g); sg=t["gtype"]==g
        Mx[i,:len(DMS)]=[np.sum(sg&(dm==x))/N for x in DMS]
        Mx[i,-2]=np.sum(sg&~np.isin(dm,DMS))/N; Mx[i,-1]=1-sg.sum()/N
    return Mx
bins=np.arange(0,55,5)
def eff_bins(dm,g):
    out=[]
    for lo,hi in zip(bins[:-1],bins[1:]):
        N=np.sum((A["gtype"]==g)&(A["gvisP"]>=lo)&(A["gvisP"]<hi)); s=np.sum((t["gtype"]==g)&(dm==g)&(t["gvisP"]>=lo)&(t["gvisP"]<hi)); out.append(s/N)
    return np.array(out)
res={}
for name,f in WPS.items():
    keep=np.ones(NP,bool) if name=="baseline" else ~f()
    rt,dm,M,P,n=emulate(t,p,keep)
    res[name]=dict(rt=rt,dm=dm,M=M,P=P,mat=matrix(dm),eff=eff_all(dm),eff2040=eff_all(dm,20,40),
                   effb={g:eff_bins(dm,g) for g in (0,1,2,10,11)},nrej=(~keep).sum(),
                   rej_bad=(~keep&bad).sum()/bad.sum(),rej_true=(~keep&true).sum()/true.sum())
W(); W("## 4. Puntos de trabajo")
W(); W("Eficiencias gen X -> DM X (todo el rango / 20-40 GeV) y fotones rechazados.")
W(); W("| punto de trabajo | fotones rechazados | frac. malos rechazados | frac. pi0 verdaderos rechazados | "+" | ".join(f"e{g}" for g in GENS)+" | e0 (20-40) | e1 (20-40) | e2 (20-40) |")
W("|---|---|---|---|"+"---|"*(len(GENS)+3))
for name,r in res.items():
    W(f"| {name} | {r['nrej']} | {r['rej_bad']:.3f} | {r['rej_true']:.4f} | "+" | ".join(f"{r['eff'][g]:.4f}" for g in GENS)+" | "+" | ".join(f"{r['eff2040'][g]:.4f}" for g in (0,1,2))+" |")
W(); W("Diferencias respecto al baseline:")
W(); W("| punto de trabajo | "+" | ".join(f"de{g}" for g in GENS)+" | de0 (20-40) | de1 (20-40) | de2 (20-40) |"); W("|---|"+"---|"*(len(GENS)+3))
for name,r in res.items():
    if name=="baseline": continue
    W(f"| {name} | "+" | ".join(f"{r['eff'][g]-base[g]:+.4f}" for g in GENS)+" | "+" | ".join(f"{r['eff2040'][g]-base2040[g]:+.4f}" for g in (0,1,2))+" |")
W(); W("### 4.1 Matrices de migracion GenTauType x RecoTauDM (fracciones por fila; 'otro' = DM no listada, 'sin reco' = gen tau sin reco emparejado)")
hdr="| gen \\ DM | "+" | ".join(str(x) for x in DMS)+" | otro | sin reco |"
for name,r in res.items():
    W(); W(f"**{name}**"); W(); W(hdr); W("|---|"+"---|"*(len(DMS)+2))
    for i,g in enumerate(GROWS): W(f"| {g} | "+" | ".join(f"{x:.3f}" for x in r["mat"][i])+" |")
W(); W("### 4.2 Eficiencias en bins de GenVisTauP (GeV)")
for g,lab in [(0,"gen 0 -> DM 0"),(1,"gen 1 -> DM 1"),(2,"gen 2 -> DM 2"),(10,"gen 10 -> DM 10"),(11,"gen 11 -> DM 11")]:
    W(); W(f"**{lab}**"); W(); W("| punto de trabajo | "+" | ".join(f"{lo}-{hi}" for lo,hi in zip(bins[:-1],bins[1:]))+" |"); W("|---|"+"---|"*(len(bins)-1))
    for name,r in res.items(): W(f"| {name} | "+" | ".join(f"{x:.3f}" for x in r["effb"][g])+" |")

# figura eficiencia vs P
sel_wp=["baseline","P<0.5 GeV","P<1 GeV","P<2 GeV","WP-A: sin pareja pi0 & P<1 GeV","WP-B: sin pareja pi0 & P<1.5 GeV","WP-C: sin pareja pi0 & P/Ppi<0.05","WP-D: sin pareja pi0 & P<0.3*sqrt(Ppi)","WP-E: P<0.3 GeV | (sin pareja pi0 & P<1 GeV)"]
cols={n:c for n,c in zip(sel_wp,["k",C[7],C[3],C[1],C[0],C[2],C[4],C[6],C[5]])}
xc=(bins[:-1]+bins[1:])/2
fig,axs=plt.subplots(2,3,figsize=(14,8))
for ax,(g,lab) in zip(axs.flat,[(0,"gen 0 -> DM 0"),(1,"gen 1 -> DM 1"),(2,"gen 2 -> DM 2"),(10,"gen 10 -> DM 10"),(11,"gen 11 -> DM 11")]):
    for n in sel_wp: ax.plot(xc,res[n]["effb"][g],"-o",ms=3.5,lw=1.5,color=cols[n],ls="--" if n.startswith("P<") else "-",label=WPS_en[n])
    ax.set_xlabel("GenVisTauP [GeV]"); ax.set_ylabel("efficiency"); ax.set_title(lab,fontsize=10)
axd=axs.flat[5]
for n in sel_wp[1:]: axd.plot(xc,res[n]["effb"][0]-res["baseline"]["effb"][0],"-o",ms=3.5,lw=1.5,color=cols[n],ls="--" if n.startswith("P<") else "-",label=WPS_en[n])
for n in sel_wp[1:]: axd.plot(xc,res[n]["effb"][2]-res["baseline"]["effb"][2],":",lw=1.2,color=cols[n])
axd.axhline(0,color="k",lw=0.8); axd.set_xlabel("GenVisTauP [GeV]"); axd.set_ylabel("delta efficiency"); axd.set_title("delta: gen0->DM0 (solid) and gen2->DM2 (dotted)",fontsize=10)
axs.flat[0].legend(fontsize=6.5); fig.suptitle("Mode identification efficiency per visible-P bin, before and after photon rejection")
fig.tight_layout(); fig.savefig("fig5_eff_vs_P.png",dpi=130); plt.close(fig)

# figura matrices antes/despues para WP-A y P<1
def show(ax,Mx,tt):
    im=ax.imshow(Mx,cmap="Blues",vmin=0,vmax=1); ax.set_xticks(range(Mx.shape[1])); ax.set_xticklabels([str(x) for x in DMS]+["other","no reco"],rotation=60,fontsize=8)
    ax.set_yticks(range(len(GROWS))); ax.set_yticklabels([str(g) for g in GROWS],fontsize=8); ax.set_xlabel("RecoTauDM"); ax.set_ylabel("GenTauType"); ax.set_title(tt,fontsize=9); ax.grid(False)
    for i in range(Mx.shape[0]):
        for j in range(Mx.shape[1]):
            if Mx[i,j]>=0.005: ax.text(j,i,f"{Mx[i,j]:.2f}",ha="center",va="center",fontsize=6.5,color="white" if Mx[i,j]>0.5 else "black")
fig,axs=plt.subplots(1,3,figsize=(15,4.8))
show(axs[0],res["baseline"]["mat"],"baseline (no cut)"); show(axs[1],res["P<1 GeV"]["mat"],"P<1 GeV on all photons"); show(axs[2],res["WP-A: sin pareja pi0 & P<1 GeV"]["mat"],"WP-A: no pi0 partner & P<1 GeV")
fig.suptitle("Migration matrix GenTauType x RecoTauDM (fraction per row)"); fig.tight_layout(); fig.savefig("fig6_migracion.png",dpi=130); plt.close(fig)

# figura masa del tau reco (1 prong) antes/despues
fig,axs=plt.subplots(1,2,figsize=(11,4.2))
for ax,g,tt in zip(axs,[0,1],["gen 0 (tau->pi nu): reco 1-prong tau mass","gen 1 (tau->pi pi0 nu): reco 1-prong tau mass"]):
    for n in ["baseline","P<1 GeV","WP-A: sin pareja pi0 & P<1 GeV"]:
        s=(t["gtype"]==g)&(res[n]["rt"]>=0)&(res[n]["rt"]<10)
        ax.hist(res[n]["M"][s],bins=np.linspace(0,2,81),histtype="step",lw=1.6,color=cols[n],label=WPS_en[n],log=True)
    ax.set_xlabel("reco tau mass [GeV]"); ax.set_ylabel("taus"); ax.set_title(tt,fontsize=9); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig("fig7_masa_tau.png",dpi=130); plt.close(fig)

# restantes tras WP-A
rA=res["WP-A: sin pareja pi0 & P<1 GeV"]; keepA=~WPS["WP-A: sin pareja pi0 & P<1 GeV"]()
rem=(t["gtype"]==0)&(rA["rt"]==1); s=rem[p["tau"]]&keepA
W(); W("## 5. Lo que queda tras WP-A")
W(); W(f"gen 0 -> reco 1: {n01} -> {rem.sum()} taus ({rem.sum()/Ntot[0]:.4f} de gen 0). Origen de los fotones restantes: "+", ".join(f"{names.get(o,o)}: {n}" for o,n in collections.Counter(p['origin'][s].tolist()).most_common()))
W(f"P de los fotones restantes: percentiles 10/25/50/75/90 = {np.percentile(p['P'][s],[10,25,50,75,90]).round(2).tolist()} GeV; m(pi+gamma) mediana {np.median(p['mpig'][s]):.2f} GeV.")
# donde van los gen 0 -> reco 1 corregidos
fixed=(t["gtype"]==0)&(t["rtype"]==1)&(rA["rt"]==0)
W(f"gen 0 -> reco 1 corregidos a reco 0 con WP-A: {fixed.sum()} ({fixed.sum()/n01:.3f}).")
# efecto en gen 1: de donde sale la ganancia/perdida
for g in (1,2):
    up=(t["gtype"]==g)&(t["rdm"]==g+1)&(rA["dm"]==g); dn=(t["gtype"]==g)&(t["rdm"]==g)&(rA["dm"]==g-1)
    W(f"gen {g}: recuperados DM{g+1}->DM{g}: {up.sum()} ({up.sum()/Ntot[g]:+.4f}); perdidos DM{g}->DM{g-1}: {dn.sum()} ({-dn.sum()/Ntot[g]:+.4f}).")

W(); W("## 6. Conclusiones y criterio recomendado")
W(); W("""1. **El foton extra no es un fragmento colineal del shower.** El 54 % es FSR real del tau (tau -> pi nu gamma; GenTauType lo ignora),
   con P mediana 1.5 GeV y dR mediana 0.18; el 43 % no enlaza a ningun foton gen (fragmento del shower / cluster partido), con P mediana 0.57 GeV
   y dR mediana 0.07. Por eso dR NO discrimina (los fotones malos estan en promedio *mas lejos* del pion que los de pi0) y un corte "dR<X & P/Ppi<Y"
   solo recupera un tercio de lo que recupera el corte en P (barrido en la seccion 3).
2. **Las variables utiles son P y, sobre todo, P/P_pion**, pero un corte simple en cualquiera de ellas se paga carisimo en gen 2 y gen 3
   (los fotones verdaderos de pi0 son blandos: 13 % tienen P<1 GeV): P<1 GeV cuesta -0.061 en e2 y -0.157 en e3; P<2 GeV, -0.196 y -0.345.
   Incluso P<0.5 GeV cuesta -0.049 en e3.
3. **La pareja pi0 es lo que separa.** Un foton de pi0 casi siempre tiene otro foton en el cono con |m_gg - m_pi0| < 50 MeV (mediana 9 MeV);
   el foton malo, por construccion, esta solo. Exigir "sin pareja" antes de aplicar el corte reduce la perdida de fotones verdaderos de pi0 en un
   factor ~10 a igual rechazo del foton malo (ROC de la figura 3: al 1 % de perdida se pasa de rechazar el 11-15 % de los malos a rechazar el 44-47 %).
4. **Puntos de trabajo (fichero completo, seccion 4).** Todos recuperan ~la mitad de los gen 0 -> reco 1 (que son el 3.9 % de gen 0); la otra mitad
   son FSR duros (P mediana 3.5 GeV, hasta 18 GeV) indistinguibles de un pi0 con un foton perdido y no conviene tocarlos.
   - WP-A (sin pareja & P<1 GeV): de0 +0.022, de1 +0.000, de2 -0.006, de3 -0.005, de10 +0.031, de11 +0.003, de12 -0.001. En 20-40 GeV: +0.023 / +0.008 / -0.006.
     Su unico coste visible esta en gen 1 por debajo de 10 GeV (-0.10 en 0-5 GeV, -0.06 en 5-10), porque alli el pi0 con un foton perdido deja un foton solo y blando.
   - WP-C (sin pareja & P/Ppi<0.05): de0 +0.020, de1 -0.003, de2 -0.004, de3 -0.003, de10 +0.017, de11 +0.002, de12 +0.001. En 20-40 GeV: +0.026 / -0.000 / -0.004,
     y no toca gen 1 a baja P (0.409 frente a 0.415 en 0-5 GeV). Es el que mas gana en 20-40 GeV y el que menos cuesta en todas partes; gana menos que WP-A
     por debajo de 10 GeV y en 3 prongs.
   - WP-B (P<1.5) y WP-D (P<0.3 sqrt(Ppi)) ganan +0.003/+0.002 mas en e0 a cambio de -0.011/-0.007 en e2; WP-E (anadir P<0.3 a todos) mejora e1/e11 (+0.008/+0.012) pero cuesta -0.022 en e3.
   - Las filas -11/-13 apenas cambian (electron: DM0 0.017 -> 0.018-0.021 segun WP; muon: 0.073 -> 0.074-0.075); el reco -20 (pi->neutron) no se toca.
5. **Recomendacion:** rechazar un foton del cono si **no tiene pareja pi0** (ningun otro foton del cono con |m_gg - 0.135| < 0.05 GeV)
   **y P_gamma / P_pion < 0.05** (WP-C). Si se prefiere maximizar la recuperacion de gen 0 y 3 prongs a baja P y se acepta el coste en gen 1 por debajo
   de 10 GeV, usar **sin pareja & P_gamma < 1 GeV** (WP-A). En ambos casos la masa del tau reco de gen 0 vuelve a la masa del pion en la mitad de los
   casos migrados (figura 7) y el pico de la rho en gen 1 no se altera.
6. Nota tecnica: para 3 prongs el P4 del tau reco del tree no coincide con la suma de sus constituyentes en un 21 % de los casos (RecoTauP > suma),
   asi que la masa recalculada solo se usa para 1 prong; el recuento de fotones (y por tanto RecoTauType/DM) si es consistente al 100 %.""")
W(); W("## Figuras")
for f,d in [("fig1_distribuciones.png","distribuciones normalizadas de P, dR, P/Ppi, P/Ptau, m(pi+gamma) y min|m_gg-m_pi0| para foton malo / bueno / espurio"),
            ("fig2_malo_por_origen.png","foton malo separado por origen gen (FSR del tau, sin match, radiacion de cargada)"),
            ("fig3_roc_foton.png","ROC por foton: corte simple y corte condicionado a 'sin pareja pi0'"),
            ("fig4_roc_tau.png","ganancia en e0 frente a perdida en e1, e2, e3 para cada familia de criterios"),
            ("fig5_eff_vs_P.png","eficiencias gen0->DM0, gen1->DM1, gen2->DM2, gen10->DM10, gen11->DM11 en bins de GenVisTauP para cada punto de trabajo"),
            ("fig6_migracion.png","matrices de migracion baseline / P<1 GeV / WP-A"),
            ("fig7_masa_tau.png","masa del tau reco 1-prong en gen 0 y gen 1 antes y despues")]:
    W(f"- `{f}`: {d}")
open("RESULTS.md","w").write("\n".join(md)+"\n")
print("\n".join(md[:60]))

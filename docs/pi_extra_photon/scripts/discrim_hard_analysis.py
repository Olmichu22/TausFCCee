"""FSR duro: gen tau->pi nu con GenTauP 20-40 GeV migra a reco 1 por un foton FSR de la linea del tau (P ~ E_haz - P_tau).
Criterios sobre el foton sin pareja pi0 y duro (P>2 GeV): dR(g,pi), m(pi+g), dR(g, tau reco). Salidas con prefijo hard_."""
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, collections
from discrim_common import *
C=["#2a78d6","#eb6834","#1baf7a","#eda100","#e87ba4","#008300","#7a5cc7","#6b6b6b","#00a3b4","#b5651d"]
plt.rcParams.update({"figure.facecolor":"white","axes.facecolor":"white","axes.grid":True,"grid.alpha":0.3,"font.size":9})
t,p,A=load(); h=np.load("hard_tables.npz"); t["gP"]=h["tau_gP"]; t["gtheta"]=h["tau_gtheta"]; t["gphi"]=h["tau_gphi"]; A["gP"]=h["all_gP"]
NT=len(t["gtype"]); NP=len(p["P"]); ti=p["tau"]
def dR_tp(th1,ph1,th2,ph2):
    d=ph1-ph2; d=np.where(d>np.pi,2*np.pi-d,d); d=np.where(d<-np.pi,2*np.pi+d,d); return np.sqrt((th1-th2)**2+d**2)
# direccion del tau reco (cargado + todos los fotones del cono) y dR del foton a ella
px=t["chpx"]+np.bincount(ti,weights=p["px"],minlength=NT); py=t["chpy"]+np.bincount(ti,weights=p["py"],minlength=NT); pz=t["chpz"]+np.bincount(ti,weights=p["pz"],minlength=NT)
tth=np.arctan2(np.hypot(px,py),pz); tph=np.arctan2(py,px)
p["dRtau"]=dR_tp(p["theta"],p["phi"],tth[ti],tph[ti])
p["dRgen"]=dR_tp(p["theta"],p["phi"],t["gtheta"][ti],t["gphi"][ti])
p["gP"]=t["gP"][ti]
md=[]; W=lambda s="": md.append(s)
W("## 7. FSR duro (GenTauP 20-40 GeV)"); W()
W("Complemento: la caida de eficiencia en 20-40 GeV de **P total** del tau (GenTauP) la produce un foton FSR de la linea del tau con P ~ E_haz - P_tau. Ficheros `hard_*`; fichero completo (2M eventos).")
# --- 7.1 migracion gen 0 en bins de GenTauP
PB=[0,10,15,20,25,30,35,40,42,44,45,46,47]
W(); W("### 7.1 gen 0 por bin de GenTauP (baseline)"); W()
W("| GenTauP [GeV] | N gen 0 | reco 0 | reco 1 | reco -20 | sin reco | foton malo: FSR tau | sin match gen | P_gamma mediana [GeV] | dR(g,pi) mediana | m(pi+g) mediana [GeV] |"); W("|---|---|---|---|---|---|---|---|---|---|---|")
bad=(p["gtype"]==0)&(p["rtype"]==1)
for lo,hi in zip(PB[:-1],PB[1:]):
    N=np.sum((A["gtype"]==0)&(A["gP"]>=lo)&(A["gP"]<hi)); s=(t["gtype"]==0)&(t["gP"]>=lo)&(t["gP"]<hi)
    w=bad&(p["gP"]>=lo)&(p["gP"]<hi); c=collections.Counter(p["origin"][w].tolist())
    W(f"| {lo}-{hi} | {N} | {np.sum(s&(t['rtype']==0))/N:.3f} | {np.sum(s&(t['rtype']==1))/N:.3f} | {np.sum(s&(t['rtype']==-20))/N:.3f} | {1-s.sum()/N:.3f} | {c[1]/max(w.sum(),1):.2f} | {c[-9]/max(w.sum(),1):.2f} | {np.median(p['P'][w]):.2f} | {np.median(p['dR'][w]):.3f} | {np.median(p['mpig'][w]):.2f} |")
# --- 7.2 poblaciones duras sin pareja
hard=(p["dmpi0"]>0.05)&(p["P"]>2)
sig=bad&(p["origin"]==1)&hard                    # FSR duro en gen0->reco1
sigw=sig&(p["gP"]>=20)&(p["gP"]<40)
comp=(p["gtype"]==1)&(p["rtype"]==1)&p["truepi0"]&hard   # rho con un foton perdido (competidor directo)
comp_all=np.isin(p["gtype"],[1,2,3])&p["truepi0"]&hard    # cualquier foton de pi0 sin pareja y duro
W(); W("### 7.2 Foton duro (P>2 GeV) sin pareja pi0: FSR frente a foton de pi0 huerfano"); W()
W("| poblacion | N | dR(g,pi) p10/50/90 | m(pi+g) p10/50/90 [GeV] | dR(g,tau reco) p10/50/90 | dR(g,tau gen) p10/50/90 | P p10/50/90 [GeV] |"); W("|---|---|---|---|---|---|---|")
q=lambda x: "/".join(f"{v:.3f}" for v in np.percentile(x,[10,50,90]))
for n,s in [("FSR duro en gen0->reco1 (todo P_tau)",sig),("FSR duro en gen0->reco1, GenTauP 20-40",sigw),("foton pi0 huerfano gen1->reco1 (rho, un foton perdido)",comp),("foton pi0 huerfano duro, gen 1,2,3",comp_all)]:
    W(f"| {n} | {s.sum()} | {q(p['dR'][s])} | {q(p['mpig'][s])} | {q(p['dRtau'][s])} | {q(p['dRgen'][s])} | {q(p['P'][s])} |")
# tabla por foton de los criterios pedidos
crit=collections.OrderedDict([("dR>0.2",p["dR"]>0.2),("dR>0.25",p["dR"]>0.25),("dR>0.3",p["dR"]>0.3),
      ("m>1.0",p["mpig"]>1.0),("m>1.2",p["mpig"]>1.2),("m>1.5",p["mpig"]>1.5),("dR>0.2 & m>1.0",(p["dR"]>0.2)&(p["mpig"]>1.0)),
      ("dR>0.2 | m>1.2",(p["dR"]>0.2)|(p["mpig"]>1.2)),("dRtau>0.1",p["dRtau"]>0.1),("dRtau>0.15",p["dRtau"]>0.15),("dRtau>0.2",p["dRtau"]>0.2)])
W(); W("Rechazo del foton (solo se aplica a fotones sin pareja pi0 y P>2 GeV):"); W()
W("| criterio | rechazo FSR duro (todo P_tau) | rechazo FSR duro (GenTauP 20-40) | perdida foton pi0 huerfano gen1->reco1 | perdida foton pi0 huerfano gen 1,2,3 |"); W("|---|---|---|---|---|")
for n,c in crit.items():
    W(f"| {n} | {np.mean(c[sig]):.3f} | {np.mean(c[sigw]):.3f} | {np.mean(c[comp]):.3f} | {np.mean(c[comp_all]):.3f} |")
# ROC
fig,axs=plt.subplots(1,2,figsize=(11,4.6))
for ax,bk,tt in zip(axs,[comp,comp_all],["background: orphan pi0 photon in gen1->reco1","background: orphan pi0 photon in gen 1,2,3"]):
    for i,(v,lab) in enumerate([("dR","dR(g, pion)"),("mpig","m(pion+g)"),("dRtau","dR(g, reco tau)"),("dRgen","dR(g, gen tau) [diagnostic]"),("P","photon P"),("fpi","P_g/P_pi")]):
        xs=np.quantile(p[v][sigw],np.linspace(0.002,0.998,200))
        es=np.array([(p[v][sigw]>x).mean() for x in xs]); eb=np.array([(p[v][bk]>x).mean() for x in xs])
        ax.plot(eb,es,color=C[i],lw=1.8,label=lab)
    ax.set_xlabel("fraction of orphan pi0 photons rejected"); ax.set_ylabel("fraction of hard FSR rejected (GenTauP 20-40)"); ax.set_xlim(0,0.5); ax.set_ylim(0,1); ax.legend(fontsize=8); ax.set_title(tt,fontsize=9)
fig.suptitle("ROC restricted to photons with no pi0 partner and P>2 GeV: reject if variable > x")
fig.tight_layout(); fig.savefig("hard_fig_roc_duro.png",dpi=130); plt.close(fig)
# distribuciones
fig,axs=plt.subplots(1,3,figsize=(13,3.8))
for ax,(v,lab,bins) in zip(axs,[("dR","dR(photon, pion) [rad]",np.linspace(0,0.4,41)),("mpig","m(pion + photon) [GeV]",np.linspace(0,5,51)),("dRtau","dR(photon, reco tau) [rad]",np.linspace(0,0.4,41))]):
    for n,s,c in [("hard FSR gen0->reco1, GenTauP 20-40",sigw,C[1]),("hard FSR gen0->reco1, all P_tau",sig,C[3]),("orphan pi0 photon gen1->reco1",comp,C[0]),("orphan pi0 photon gen 1,2,3",comp_all,C[2])]:
        ax.hist(p[v][s],bins=bins,histtype="step",density=True,lw=1.6,color=c,label=n)
    ax.set_xlabel(lab); ax.set_ylabel("density (norm. to 1)")
axs[0].legend(fontsize=7); fig.suptitle("Photons with no pi0 partner and P>2 GeV: tau-line FSR vs orphan pi0 photon")
fig.tight_layout(); fig.savefig("hard_fig_dist_duro.png",dpi=130); plt.close(fig)
# --- 7.3 nivel de tau
GENS=[0,1,2]
def eff(dm,g,var,lo,hi):
    N=np.sum((A["gtype"]==g)&(A[var]>=lo)&(A[var]<hi)); s=np.sum((t["gtype"]==g)&(dm==g)&(t[var]>=lo)&(t[var]<hi)); return s/N
def effg(dm,g): return np.sum((t["gtype"]==g)&(dm==g))/np.sum(A["gtype"]==g)
wpc=(p["dmpi0"]>0.05)&(p["fpi"]<0.05)
WPS=collections.OrderedDict([("baseline",np.zeros(NP,bool)),("WP-C",wpc)])
for n,c in crit.items():
    if n.startswith("dRtau>0.1") : continue
    WPS["H: "+n]=hard&c
for n in ["m>1.0","m>1.2","dR>0.2 & m>1.0","dR>0.2 | m>1.2"]: WPS["WP-C + H: "+n]=wpc|(hard&crit[n])
VB=np.arange(0,55,5)
res={}
for n,rej in WPS.items():
    rt,dm,M,P,nn=emulate(t,p,~rej)
    res[n]=dict(dm=dm,g=[effg(dm,g) for g in GENS],gP={g:[eff(dm,g,"gP",lo,hi) for lo,hi in zip(PB[:-1],PB[1:])] for g in GENS},
                vis={g:[eff(dm,g,"gvisP",lo,hi) for lo,hi in zip(VB[:-1],VB[1:])] for g in GENS},
                w2040=[eff(dm,g,"gP",20,40) for g in GENS],v2040=[eff(dm,g,"gvisP",20,40) for g in GENS])
b=res["baseline"]
W(); W("### 7.3 Nivel de tau: eficiencias gen X -> DM X"); W()
W("Los criterios H se aplican solo a fotones sin pareja pi0 con P>2 GeV; 'WP-C + H' anade el corte blando recomendado antes (sin pareja & P/Ppi<0.05)."); W()
W("| criterio | e0 | e1 | e2 | e0 GenTauP 20-40 | e1 GenTauP 20-40 | e2 GenTauP 20-40 | e0 VisP 20-40 | e1 VisP 20-40 | e2 VisP 20-40 |"); W("|---|"+"---|"*9)
for n,r in res.items():
    d=lambda a,c: f"{a:.4f}" if n=="baseline" else f"{a:.4f} ({a-c:+.4f})"
    W(f"| {n} | "+" | ".join(d(r['g'][i],b['g'][i]) for i in range(3))+" | "+" | ".join(d(r['w2040'][i],b['w2040'][i]) for i in range(3))+" | "+" | ".join(d(r['v2040'][i],b['v2040'][i]) for i in range(3))+" |")
for g in GENS:
    W(); W(f"**gen {g} -> DM {g} en bins de GenTauP (GeV)**"); W(); W("| criterio | "+" | ".join(f"{lo}-{hi}" for lo,hi in zip(PB[:-1],PB[1:]))+" |"); W("|---|"+"---|"*(len(PB)-1))
    for n,r in res.items(): W(f"| {n} | "+" | ".join(f"{x:.3f}" for x in r["gP"][g])+" |")
for g in (1,2):
    W(); W(f"**gen {g} -> DM {g} en bins de GenVisTauP (GeV)**"); W(); W("| criterio | "+" | ".join(f"{lo}-{hi}" for lo,hi in zip(VB[:-1],VB[1:]))+" |"); W("|---|"+"---|"*(len(VB)-1))
    for n,r in res.items(): W(f"| {n} | "+" | ".join(f"{x:.3f}" for x in r["vis"][g])+" |")
# figura eff vs GenTauP
show=["baseline","WP-C","H: dR>0.2","H: dR>0.3","H: m>1.0","H: m>1.2","H: m>1.5","H: dR>0.2 & m>1.0","WP-C + H: m>1.2","WP-C + H: dR>0.2 | m>1.2"]
cols=dict(zip(show,["k",C[7],C[0],C[8],C[1],C[3],C[9],C[2],C[4],C[6]]))
xc=[(lo+hi)/2 for lo,hi in zip(PB[:-1],PB[1:])]
fig,axs=plt.subplots(2,3,figsize=(14,8))
for ax,g in zip(axs[0],GENS):
    for n in show: ax.plot(xc,res[n]["gP"][g],"-o",ms=3.5,lw=1.5,color=cols[n],ls="--" if n.startswith("H:") else "-",label=n)
    ax.set_xlabel("GenTauP [GeV]"); ax.set_ylabel("efficiency"); ax.set_title(f"gen {g} -> DM {g}",fontsize=10)
for ax,g in zip(axs[1],GENS):
    for n in show[1:]: ax.plot(xc,np.array(res[n]["gP"][g])-np.array(b["gP"][g]),"-o",ms=3.5,lw=1.5,color=cols[n],ls="--" if n.startswith("H:") else "-",label=n)
    ax.axhline(0,color="k",lw=0.8); ax.set_xlabel("GenTauP [GeV]"); ax.set_ylabel("delta efficiency vs baseline"); ax.set_title(f"delta gen {g} -> DM {g}",fontsize=10)
axs[0][0].legend(fontsize=6.5); fig.suptitle("Mode efficiency in bins of total gen tau P: hard-FSR criteria (full file)")
fig.tight_layout(); fig.savefig("hard_fig_eff_vs_GenTauP.png",dpi=130); plt.close(fig)
# figura eff vs VisP para gen1/gen2
fig,axs=plt.subplots(1,2,figsize=(11,4.3))
for ax,g in zip(axs,(1,2)):
    for n in show: ax.plot((VB[:-1]+VB[1:])/2,res[n]["vis"][g],"-o",ms=3.5,lw=1.5,color=cols[n],ls="--" if n.startswith("H:") else "-",label=n)
    ax.set_xlabel("GenVisTauP [GeV]"); ax.set_ylabel("efficiency"); ax.set_title(f"gen {g} -> DM {g}",fontsize=10)
axs[0].legend(fontsize=6.5); fig.suptitle("Cost of hard-FSR criteria in gen 1 and gen 2, vs visible P")
fig.tight_layout(); fig.savefig("hard_fig_eff_vs_VisP.png",dpi=130); plt.close(fig)
open("hard_RESULTS_section.md","w").write("\n".join(md)+"\n")
print("\n".join(md))

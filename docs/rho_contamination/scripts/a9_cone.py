"""Cono dependiente de P: fotones libres (sin tau) a 0.4<dR<X del pion lider de taus 1-prong."""
import numpy as np, glob, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from common import *
T,R,G,C=load(); out=open("a9_cone.md","w")
def P(*a): print(*a); print(*a,file=out)
c=cls(T["rtype"]); nt=len(T["ev"]); gt=T["gtype"]
U={}
for f in sorted(glob.glob("uncone_*.npz")):
    d=np.load(f)
    for k in d: U.setdefault(k,[]).append(d[k])
U={k:np.concatenate(v) for k,v in U.items()}
# fila del gen tau
tk=T["ev"].astype(np.int64)*1000+T["key"]; o=np.argsort(tk); tk_s=tk[o]; tr_s=np.arange(nt)[o]
q=U["ev"].astype(np.int64)*1000+U["gkey"]; i=np.clip(np.searchsorted(tk_s,q),0,nt-1); ok=tk_s[i]==q
U["tau"]=np.where(ok,tr_s[i],-1); assert ok.all()
ti=U["tau"]
U["dRpi"]=ang(U["th"] if "th" in U else np.nan,0,0,0) if False else None
# no guardamos theta/phi del foton: usamos dR a la direccion del reco tau. Para 0g el eje es el pion; para 1g/2g es el tau. Aceptable como aproximacion.
one=np.isin(T["rtype"],np.arange(0,10))
s=one[ti]
U["own"]=(U["gorigin"]==0)&(U["gtaukey"]==T["key"][ti])
P("## Fotones libres (sin tau) a dR<1 del eje del reco tau 1-prong: que son, por tipo gen del tau y anillo de dR")
P("| gen | anillo dR | N fotones | π0 propio | FSR | ISR/rad | otro tau π0 | sin match | P p50 | P p50 π0 propio |")
P("|---|---|---|---|---|---|---|---|---|---|")
for g,nm in [(1,"1 ρ"),(0,"0 π"),(2,"2"),(-11,"e")]:
    for lo,hi in [(0,0.4),(0.4,0.6),(0.6,0.8),(0.8,1.0)]:
        ss=s&(gt[ti]==g)&(U["dR"]>=lo)&(U["dR"]<hi); N=ss.sum()
        if N<50: continue
        oth=(U["gorigin"][ss]==0)&~U["own"][ss]
        P("| %s | %g-%g | %d | %.3f | %.3f | %.3f | %.3f | %.3f | %.2f | %.2f |"%(nm,lo,hi,N,U["own"][ss].mean(),(U["gorigin"][ss]==1).mean(),(U["gorigin"][ss]==2).mean(),oth.mean(),(U["gorigin"][ss]==-2).mean(),np.median(U["P"][ss]),np.median(U["P"][ss&U["own"]]) if (ss&U["own"]).any() else np.nan))
# emulacion: cono X para reco taus con P < Pcut (P visible reco)
P("\n## Emulacion: anadir al tau los fotones libres con dR<X si P_reco(tau) < Pcut. Fracciones de gen 1 / gen 0 en cada clase.")
P("| X | Pcut | gen1 0γ | gen1 1γ | gen1 2γ | gen1 3γ+ | gen1 2γ (P_vis<10) | gen0 0γ | gen0 1γ | gen0 0γ (P_vis<10) | gen0 1γ (P_vis<10) | e→2γ |")
P("|---|---|---|---|---|---|---|---|---|---|---|---|")
def row(nm,rt):
    cc=cls(rt); s1=gt==1; s0=gt==0; lo1=s1&(T["gvisP"]<10); lo0=s0&(T["gvisP"]<10); se=gt==-11
    P("| %s | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f |"%(nm,(cc[s1]=="0g").mean(),(cc[s1]=="1g").mean(),(cc[s1]=="2g").mean(),(cc[s1]=="3g+").mean(),(cc[lo1]=="2g").mean(),(cc[s0]=="0g").mean(),(cc[s0]=="1g").mean(),(cc[lo0]=="0g").mean(),(cc[lo0]=="1g").mean(),(cc[se]=="2g").mean()))
row("0.4 | - (baseline)",T["rtype"])
for X in (0.6,0.8,1.0):
    for Pc in (5,10,20,100):
        add=s&(U["dR"]<X)&(T["rP"][ti]<Pc)
        nadd=np.bincount(ti[add],minlength=nt)
        rt=T["rtype"].copy(); rt[one]=np.minimum(T["nph"][one]+nadd[one],9)
        row("%g | %g"%(X,Pc),rt)
# variante: solo fotones con P > 0.3 GeV
P("\nVariante: solo fotones libres con P>0.3 GeV")
P("| X | Pcut | gen1 0γ | gen1 1γ | gen1 2γ | gen1 3γ+ | gen1 2γ (P_vis<10) | gen0 0γ | gen0 1γ | gen0 0γ (P_vis<10) | gen0 1γ (P_vis<10) | e→2γ |")
P("|---|---|---|---|---|---|---|---|---|---|---|---|")
for X in (0.6,0.8):
    for Pc in (5,10,20):
        add=s&(U["dR"]<X)&(T["rP"][ti]<Pc)&(U["P"]>0.3)
        nadd=np.bincount(ti[add],minlength=nt)
        rt=T["rtype"].copy(); rt[one]=np.minimum(T["nph"][one]+nadd[one],9)
        row("%g | %g"%(X,Pc),rt)
# figura: dR de fotones propios vs ajenos para taus de P<10
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
bins=np.linspace(0,1,41)
for g,nm in [(1,"gen 1: π0 propio"),]:
    ss=s&(gt[ti]==1)&U["own"]&(T["rP"][ti]<10); ax[0].hist(U["dR"][ss],bins,histtype="step",label=nm+" (P_reco<10)")
ss=s&(gt[ti]==1)&~U["own"]&(T["rP"][ti]<10); ax[0].hist(U["dR"][ss],bins,histtype="step",label="gen 1: otros fotones (P_reco<10)")
ss=s&(gt[ti]==0)&(T["rP"][ti]<10); ax[0].hist(U["dR"][ss],bins,histtype="step",label="gen 0: todos (P_reco<10)")
ax[0].axvline(0.4,color="k",ls="--"); ax[0].set_xlabel("dR(foton libre, eje del reco tau)"); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3); ax[0].set_yscale("log")
ss=s&(gt[ti]==1)&U["own"]&(U["dR"]>0.4)&(U["dR"]<0.8); ax[1].hist(np.clip(U["P"][ss],0,5),np.linspace(0,5,51),histtype="step",label="gen 1: π0 propio, 0.4<dR<0.8")
ss=s&(gt[ti]==0)&(U["dR"]>0.4)&(U["dR"]<0.8); ax[1].hist(np.clip(U["P"][ss],0,5),np.linspace(0,5,51),histtype="step",label="gen 0: todos, 0.4<dR<0.8")
ax[1].set_xlabel("P del foton libre [GeV]"); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3); ax[1].set_yscale("log")
plt.tight_layout(); plt.savefig("../figs/cone_fig1_fotones_libres.png",dpi=120)
out.close()

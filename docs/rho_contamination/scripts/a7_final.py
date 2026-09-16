"""Candidatos finales: corte vigente (partes) + recombinacion en 3g+. Eficiencias por P, transicion, pureza, ROC.
Ademas: 1g como rho (fusion) y su pureza tras el corte vigente con una ventana en m(pi+gamma)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from common import *
T,R,G,C=load(); R=sort_by_tau(R); T,R=photon_pairing(T,R); out=open("a7_final.md","w")
def P(*a): print(*a); print(*a,file=out)
c=cls(T["rtype"]); ti=R["tau"]; nt=len(T["ev"]); gt=T["gtype"]; n=len(ti)
d5=np.load("a5_tau.npz"); best_d=d5["best_d"]; pairP=d5["pairP"]; restP=d5["restP"]; gpi0P=d5["gpi0P"]
paired=np.isfinite(R["dmpi0"])&(R["dmpi0"]<0.05); one=np.isin(T["rtype"][ti],np.arange(0,10))&(T["nneu"][ti]==0)
drop_soft=one&~paired&(R["fpi"]<0.05); drop_hard=one&~paired&(R["P"]>2)&(R["mpig"]>1.2)
three=(c=="3g+"); ratio=restP/np.maximum(pairP,1e-6)
def apply(keep_mask,recomb):
    """keep_mask por foton (corte vigente); recomb por tau: los 3g+ que cumplen se reclasifican como 2g (fusionando fotones)."""
    rt=emulate(T,R,keep_mask)
    cc=cls(rt); s=recomb&(cc=="3g+"); rt=rt.copy(); rt[s]=2
    return rt
def recomb_after(keep_mask,X,f):
    """recomputa m(all) y mejor pareja con los fotones que sobreviven al corte: aproximacion = usar valores baseline si no se quito nada del tau, si se quito algo usar solo criterio B sobre... (simplificacion: recalcular masa total)"""
    w=keep_mask.astype(float)
    E=np.bincount(ti,weights=R["P"]*w,minlength=nt); px=np.bincount(ti,weights=R["px"]*w,minlength=nt); py=np.bincount(ti,weights=R["py"]*w,minlength=nt); pz=np.bincount(ti,weights=R["pz"]*w,minlength=nt)
    mall=np.sqrt(np.maximum(E**2-px**2-py**2-pz**2,0))
    # pareja: si no se quito nada, la de baseline; si se quito, exigimos que la pareja baseline sobreviva (ambos fotones conservados)
    removed=np.bincount(ti,weights=~keep_mask,minlength=nt)>0
    restP2=np.bincount(ti,weights=R["P"]*w,minlength=nt)-pairP
    B=(best_d<0.05)&(restP2/np.maximum(pairP,1e-6)<f)
    A=mall<X
    return A|B, A, B, mall
cands={}
cands["baseline"]=T["rtype"].copy()
cands["vigente (blanda+dura)"]=apply(~(drop_soft|drop_hard),np.zeros(nt,bool))
cands["solo dura"]=apply(~drop_hard,np.zeros(nt,bool))
for X,f in [(0.25,0.1),(0.25,0.2),(0.3,0.2),(0.3,0.3),(0.4,0.3)]:
    rc,_,_,_=recomb_after(~drop_hard,X,f); cands["dura + recomb(m<%.2f | Prest/Ppar<%.1f)"%(X,f)]=apply(~drop_hard,rc)
    rc,_,_,_=recomb_after(~(drop_soft|drop_hard),X,f); cands["vigente + recomb(m<%.2f | Prest/Ppar<%.1f)"%(X,f)]=apply(~(drop_soft|drop_hard),rc)
rc,A,B,_=recomb_after(~drop_hard,0.25,0.0); cands["dura + recomb(solo m<0.25)"]=apply(~drop_hard,A)
rc,A,B,_=recomb_after(~drop_hard,0.0,0.2); cands["dura + recomb(solo Prest/Ppar<0.2)"]=apply(~drop_hard,B)
n1=(gt==1).sum()
P("## Candidatos: eficiencias globales y pureza de la clase 2g")
P("| version | gen1->2g | gen1->1g | gen1->3g+ | gen1->0g | gen2->2g | gen2->3g+ | gen3->3g+ | gen0->0g | gen0->1g | e->2g | pureza 2g | N 2g |")
P("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for nm,rt in cands.items():
    cc=cls(rt); s1=gt==1; s2=gt==2; s3=gt==3; s0=gt==0; se=gt==-11
    P("| %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.4f | %.3f | %d |"%(nm,(cc[s1]=="2g").mean(),(cc[s1]=="1g").mean(),(cc[s1]=="3g+").mean(),(cc[s1]=="0g").mean(),(cc[s2]=="2g").mean(),(cc[s2]=="3g+").mean(),(cc[s3]=="3g+").mean(),(cc[s0]=="0g").mean(),(cc[s0]=="1g").mean(),(cc[se]=="2g").mean(),((cc=="2g")&s1).sum()/(cc=="2g").sum(),(cc=="2g").sum()))
# eficiencia gen1->2g por bin de P visible y de P total
for var,bins in [("gvisP",[0,5,10,15,20,25,30,35,40,46]),("gP",[0,20,30,40,44,46])]:
    P("\n### gen1 -> 2g por bin de %s"%var)
    names=["baseline","solo dura","dura + recomb(m<0.25 | Prest/Ppar<0.2)","dura + recomb(m<0.30 | Prest/Ppar<0.3)","vigente + recomb(m<0.25 | Prest/Ppar<0.2)"]
    P(eff_table(T,gt==1,[cands[k] for k in names],names,bins,var,2))
P("\n### pureza de la clase 2g por bin de P visible reco")
names=["baseline","solo dura","dura + recomb(m<0.25 | Prest/Ppar<0.2)","dura + recomb(m<0.30 | Prest/Ppar<0.3)"]
P("| rP | "+" | ".join(names)+" |"); P("|"+"---|"*(len(names)+1))
for lo,hi in [(0,10),(10,20),(20,30),(30,40),(40,50),(0,50)]:
    s=(T["rP"]>=lo)&(T["rP"]<hi)
    P("| %g-%g | "%(lo,hi)+" | ".join("%.3f"%(((cls(cands[k])=="2g")&s&(gt==1)).sum()/max(((cls(cands[k])=="2g")&s).sum(),1)) for k in names)+" |")
# ROC: barrido fino de X y f sobre "dura + recomb"
P("\n## ROC de la recombinacion (sobre 'solo dura'): ganancia gen1->2g frente a contaminacion anadida a la clase 2g")
base=cands["solo dura"]; cb=cls(base); n2g=(cb=="2g").sum(); pts=[]
for X in [0,0.2,0.25,0.3,0.35,0.4]:
    for f in [0,0.1,0.2,0.3,0.5]:
        rc,_,_,_=recomb_after(~drop_hard,X,f); rt=apply(~drop_hard,rc); cc=cls(rt)
        gain=((cc=="2g")&(gt==1)).sum()-((cb=="2g")&(gt==1)).sum(); cont=((cc=="2g")&(gt!=1)).sum()-((cb=="2g")&(gt!=1)).sum()
        pts.append((X,f,gain/n1,cont/n2g,((cc=="2g")&(gt==1)).sum()/(cc=="2g").sum()))
P("| m(all)< | Prest/Ppar< | Δeff gen1->2g | contaminacion anadida/clase 2g | pureza 2g |"); P("|---|---|---|---|---|")
for p in pts: P("| %g | %g | +%.4f | %.4f | %.3f |"%p)
# ---- 1g como rho: m(pi+g) por poblacion tras el corte vigente
P("\n## Clase 1g tras el corte vigente: ¿es un rho con el pi0 fusionado?")
rt=cands["vigente (blanda+dura)"]; cc=cls(rt); s=(cc=="1g")
# m(pi+g) por tau con el foton superviviente (solo uno): usar mpig del foton conservado
keep=~(drop_soft|drop_hard); kept=np.full(nt,-1); w=np.where(keep&(cc[ti]=="1g"))[0]; kept[ti[w]]=w
ok=s&(kept>=0); m=R["mpig"][np.maximum(kept,0)]
P("| gen | N en 1g | m(pi+g) p10/50/90 | frac 0.5<m<1.1 | frac m<0.3 |"); P("|---|---|---|---|---|")
for g,nm in [(1,"1 rho"),(0,"0 pi"),(2,"2"),(-11,"e"),(-13,"mu"),(10,"10")]:
    ss=ok&(gt==g)
    if ss.sum()<100: continue
    P("| %s | %d | %.2f/%.2f/%.2f | %.3f | %.3f |"%(nm,ss.sum(),*np.percentile(m[ss],[10,50,90]),((m[ss]>0.5)&(m[ss]<1.1)).mean(),(m[ss]<0.3).mean()))
win=ok&(m>0.5)&(m<1.1)
P("\n1g con 0.5<m(pi+g)<1.1: N=%d, pureza gen1 = %.3f (sin ventana: %.3f). Fraccion de todos los gen1 que caen ahi: %.3f"%(win.sum(),(win&(gt==1)).sum()/win.sum(),(ok&(gt==1)).sum()/ok.sum(),(win&(gt==1)).sum()/n1))
# fusionados vs perdido en 1g gen1: usar P_reco/P_gen del foton (>1.15 ~ fusionado)
gi=G["tau"]
kk=kept[ok&(gt==1)]; r=R["P"][kk]/R["gP"][kk]; fused=r>1.15
P("gen1 en 1g: frac P_reco/P_gen>1.15 (fusion): %.3f; m(pi+g) p50 fusionados %.2f, no fusionados %.2f; frac en ventana: fusionados %.3f, no fusionados %.3f"%(fused.mean(),np.median(m[ok&(gt==1)][fused]),np.median(m[ok&(gt==1)][~fused]),((m[ok&(gt==1)][fused]>0.5)&(m[ok&(gt==1)][fused]<1.1)).mean(),((m[ok&(gt==1)][~fused]>0.5)&(m[ok&(gt==1)][~fused]<1.1)).mean()))
np.savez_compressed("a7_cands.npz",**{k.replace(" ","_").replace("|","or").replace("(","").replace(")","").replace("<","lt").replace("/","_").replace("+","p").replace(".",""):v for k,v in cands.items()})
# ---- figuras
fig,ax=plt.subplots(1,3,figsize=(16,4.5))
bins=[0,5,10,15,20,25,30,35,40,46]; xc=0.5*(np.array(bins[1:])+np.array(bins[:-1]))
for k,st in [("baseline","k-"),("solo dura","C0--"),("dura + recomb(m<0.25 | Prest/Ppar<0.2)","C1-"),("dura + recomb(m<0.30 | Prest/Ppar<0.3)","C2-"),("vigente (blanda+dura)","C3:")]:
    cc=cls(cands[k]); e=[((cc=="2g")&(gt==1)&(T["gvisP"]>=lo)&(T["gvisP"]<hi)).sum()/max(((gt==1)&(T["gvisP"]>=lo)&(T["gvisP"]<hi)).sum(),1) for lo,hi in zip(bins[:-1],bins[1:])]
    ax[0].plot(xc,e,st,marker="o",ms=3,label=k)
ax[0].set_xlabel("P visible gen del tau [GeV]"); ax[0].set_ylabel("eff gen 1 -> reco 2 (rho puro)"); ax[0].legend(fontsize=7); ax[0].grid(alpha=.3)
for k,st in [("baseline","k-"),("solo dura","C0--"),("dura + recomb(m<0.25 | Prest/Ppar<0.2)","C1-"),("dura + recomb(m<0.30 | Prest/Ppar<0.3)","C2-")]:
    cc=cls(cands[k]); p=[(((cc=="2g")&(gt==1)&(T["rP"]>=lo)&(T["rP"]<hi)).sum()/max(((cc=="2g")&(T["rP"]>=lo)&(T["rP"]<hi)).sum(),1)) for lo,hi in zip(bins[:-1],bins[1:])]
    ax[1].plot(xc,p,st,marker="o",ms=3,label=k)
ax[1].set_xlabel("P visible reco del tau [GeV]"); ax[1].set_ylabel("pureza gen 1 de la clase reco 2"); ax[1].legend(fontsize=7); ax[1].grid(alpha=.3)
pts=np.array(pts); 
for f in [0,0.1,0.2,0.3,0.5]:
    m_=pts[:,1]==f; ax[2].plot(pts[m_,3],pts[m_,2],marker="o",ms=3,label="Prest/Ppar<%g"%f)
for p in pts: ax[2].annotate("%g"%p[0],(p[3],p[2]),fontsize=6)
ax[2].set_xlabel("contaminacion anadida / clase 2g"); ax[2].set_ylabel("Δeff gen1->2g"); ax[2].legend(fontsize=7); ax[2].grid(alpha=.3); ax[2].set_title("etiqueta = corte en m(all fotones)")
plt.tight_layout(); plt.savefig("../figs/final_fig1_eff_pureza_roc.png",dpi=120)
# figura migracion por P baseline vs recomendado (gen1 -> clases)
fig,ax=plt.subplots(1,2,figsize=(12,4.5))
for a,k in zip(ax,["baseline","dura + recomb(m<0.25 | Prest/Ppar<0.2)"]):
    cc=cls(cands[k])
    for cl,st in [("2g","C0"),("1g","C1"),("3g+","C2"),("0g","C3"),("pi->n","C4"),("nomatch","C7")]:
        e=[((cc==cl)&(gt==1)&(T["gvisP"]>=lo)&(T["gvisP"]<hi)).sum()/max(((gt==1)&(T["gvisP"]>=lo)&(T["gvisP"]<hi)).sum(),1) for lo,hi in zip(bins[:-1],bins[1:])]
        a.plot(xc,e,color=st,marker="o",ms=3,label=cl)
    a.set_title("gen 1 (rho): "+k,fontsize=9); a.set_xlabel("P visible gen [GeV]"); a.set_ylabel("fraccion"); a.grid(alpha=.3); a.legend(fontsize=8); a.set_ylim(0,0.7)
plt.tight_layout(); plt.savefig("../figs/final_fig2_migracion_vs_P.png",dpi=120)
out.close()

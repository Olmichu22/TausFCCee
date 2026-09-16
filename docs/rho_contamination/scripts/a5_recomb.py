"""Criterios de recombinacion de fotones en 3g+ para devolver el rho a la clase 2g, y calidad del pi0 resultante."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from common import *
T,R,G,C=load(); R=sort_by_tau(R); T,R=photon_pairing(T,R); out=open("a5_recomb.md","w")
def P(*a): print(*a); print(*a,file=out)
c=cls(T["rtype"]); ti=R["tau"]; nt=len(T["ev"]); gt=T["gtype"]; gi=G["tau"]
# pi0 gen por tau
px,py,pz=cart(G["P"],G["th"],G["phi"]); sx=np.bincount(gi,weights=px,minlength=nt); sy=np.bincount(gi,weights=py,minlength=nt); sz=np.bincount(gi,weights=pz,minlength=nt); gpi0P=np.sqrt(sx**2+sy**2+sz**2)
# --- por tau: mejor pareja (min |m-mpi0| entre todas las parejas), P de la pareja, P y numero del resto
n=len(ti); best_d=np.full(nt,np.inf); bi=np.full(nt,-1); bj=np.full(nt,-1)
kmax=int(T["nph"].max())
for k in range(1,kmax):
    i=np.arange(0,n-k); j=i+k; same=ti[i]==ti[j]; i=i[same]; j=j[same]
    if len(i)==0: break
    cos=(R["px"][i]*R["px"][j]+R["py"][i]*R["py"][j]+R["pz"][i]*R["pz"][j])/(R["P"][i]*R["P"][j])
    m=np.sqrt(np.maximum(2*R["P"][i]*R["P"][j]*(1-np.clip(cos,-1,1)),0)); d=np.abs(m-MPI0); t=ti[i]
    # min por tau: procesar por orden de d
    o=np.argsort(d)[::-1]   # peor primero, mejor sobrescribe
    upd=d[o]<best_d[t[o]]
    # como varias parejas del mismo tau pueden estar en el mismo k, hacer un bucle con np.minimum.at y luego recuperar
    np.minimum.at(best_d,t,d)
    sel=d==best_d[t]; bi[t[sel]]=i[sel]; bj[t[sel]]=j[sel]
has=bi>=0
pairP=np.zeros(nt); pairPx=np.zeros(nt); pairPy=np.zeros(nt); pairPz=np.zeros(nt)
for arr in (bi,bj):
    a=np.maximum(arr,0); pairPx+=np.where(has,R["px"][a],0); pairPy+=np.where(has,R["py"][a],0); pairPz+=np.where(has,R["pz"][a],0)
pairP=np.sqrt(pairPx**2+pairPy**2+pairPz**2)
inpair=np.zeros(n,bool); inpair[bi[has]]=True; inpair[bj[has]]=True
restP=np.bincount(ti,weights=R["P"]*(~inpair),minlength=nt); restN=np.bincount(ti,weights=(~inpair),minlength=nt)
restPmax=np.zeros(nt); np.maximum.at(restPmax,ti[~inpair],R["P"][~inpair])
# resto: m del resto con la pareja? masa total ya en T["mgg_all"]
T["best_d"]=best_d; T["pairP"]=pairP; T["restP"]=restP; T["restPmax"]=restPmax
three=(c=="3g+")
P("## Clase reco 3g+: mejor pareja de fotones y resto, por tipo gen")
P("| gen | N | frac |m_par-m_pi0|<0.05 | <0.03 | P_resto/P_par p10/50/90 | frac P_resto/P_par<0.2 | <0.3 | <0.5 | m(all)<0.25 | m(all)<0.3 |")
P("|---|---|---|---|---|---|---|---|---|---|")
for g,nm in [(1,"1 rho"),(2,"2 pi2pi0"),(3,"3 pi3pi0"),(0,"0 pi"),(-11,"e")]:
    s=three&(gt==g)
    if s.sum()<100: continue
    r=restP[s]/np.maximum(pairP[s],1e-6)
    P("| %s | %d | %.3f | %.3f | %.2f/%.2f/%.2f | %.3f | %.3f | %.3f | %.3f | %.3f |"%(nm,s.sum(),(best_d[s]<0.05).mean(),(best_d[s]<0.03).mean(),*np.percentile(r,[10,50,90]),(r<0.2).mean(),(r<0.3).mean(),(r<0.5).mean(),(T["mgg_all"][s]<0.25).mean(),(T["mgg_all"][s]<0.3).mean()))
# --- criterios: reclasificar 3g+ como rho (2g) ---
crit={}
crit["A: m(all)<0.25"]=T["mgg_all"]<0.25
crit["A: m(all)<0.30"]=T["mgg_all"]<0.30
crit["A: m(all)<0.40"]=T["mgg_all"]<0.40
for w in (0.05,):
    for f in (0.1,0.2,0.3,0.5):
        crit["B: par<%.2f & Prest/Ppar<%.1f"%(w,f)]=(best_d<w)&(restP/np.maximum(pairP,1e-6)<f)
crit["B': par<0.05 & Prest_max/Ppar<0.2"]=(best_d<0.05)&(restPmax/np.maximum(pairP,1e-6)<0.2)
crit["C: par<0.05 & m(all)<0.40"]=(best_d<0.05)&(T["mgg_all"]<0.40)
crit["D: m(all)<0.3 | (par<0.05 & Prest/Ppar<0.3)"]=(T["mgg_all"]<0.3)|((best_d<0.05)&(restP/np.maximum(pairP,1e-6)<0.3))
n2g=(c=="2g").sum(); n1=(gt==1).sum()
P("\n## Reclasificar 3g+ -> 2g con cada criterio (fichero completo)")
P("Ganancia = gen1 3g+ que pasa a 2g / todos los gen1 ; contaminacion = gen2+gen3+otros 3g+ que pasan a 2g, relativo al tamano actual de la clase 2g (%d)"%n2g)
P("| criterio | gen1 recuperados | Δeff(gen1->2g) | gen2 -> 2g | gen3 -> 2g | otros -> 2g | contaminacion añadida / clase 2g | pureza gen1 de la clase 2g resultante |")
P("|---|---|---|---|---|---|---|---|")
base_pure=((c=="2g")&(gt==1)).sum()/n2g
P("| baseline | 0 | 0 | 0 | 0 | 0 | 0 | %.3f |"%base_pure)
for nm,m in crit.items():
    s=three&m; g1=(s&(gt==1)).sum(); g2=(s&(gt==2)).sum(); g3=(s&(gt==3)).sum(); go=(s&~np.isin(gt,[1,2,3])).sum()
    P("| %s | %d | +%.3f | %d | %d | %d | %.3f | %.3f |"%(nm,g1,g1/n1,g2,g3,go,(g2+g3+go)/n2g,(((c=="2g")&(gt==1)).sum()+g1)/(n2g+s.sum())))
# --- calidad del pi0: P(pi0 estimado)/P(pi0 gen) para gen1 3g+ reclasificados: suma de todos vs mejor pareja
P("\n## Calidad del pi0 en gen1 3g+ (criterio D) : P_est/P_gen(pi0)")
s=three&(gt==1)&crit["D: m(all)<0.3 | (par<0.05 & Prest/Ppar<0.3)"]&(gpi0P>0)
for nm,est in [("suma de todos los fotones",T["Pgg_all"]),("solo mejor pareja",pairP)]:
    r=est[s]/gpi0P[s]; P("- %s: p10/50/90 = %.3f/%.3f/%.3f ; frac |r-1|<0.1: %.3f"%(nm,*np.percentile(r,[10,50,90]),(np.abs(r-1)<0.1).mean()))
s2=(c=="2g")&(gt==1)&(gpi0P>0); r=T["Pgg_all"][s2]/gpi0P[s2]; P("- referencia 2g: p10/50/90 = %.3f/%.3f/%.3f ; frac |r-1|<0.1: %.3f"%(*np.percentile(r,[10,50,90]),(np.abs(r-1)<0.1).mean()))
# 1g: el foton conservado, P_reco/P_gen (¿absorbio al otro?)
P("\n## 1g: P_reco/P_gen del foton de pi0 conservado (¿absorbe al perdido?)")
for k in ["1g","2g"]:
    s=(gt[ti]==1)&(c[ti]==k)&(R["cat"]==0); r=R["P"][s]/R["gP"][s]
    P("- %s: p10/50/90 = %.3f/%.3f/%.3f ; frac r>1.15: %.3f ; frac r>1.3: %.3f"%(k,*np.percentile(r,[10,50,90]),(r>1.15).mean(),(r>1.3).mean()))
# en 1g, separar segun dR(gg) gen <0.02 (fusion) o no
two=np.bincount(gi,minlength=nt)==2
idx=np.where(((gt==1)&(c=="1g")&two)[gi])[0]; o=np.lexsort((-G["P"][idx],gi[idx])); idx=idx[o]; a=idx[0::2]; b=idx[1::2]
dgg=ang(G["th"][a],G["phi"][a],G["th"][b],G["phi"][b]); tau_a=gi[a]
# foton reco conservado por tau
kept_r=np.full(nt,-1); s=(R["cat"]==0)&(c[ti]=="1g"); kept_r[ti[s]]=np.where(s)[0]
kk=kept_r[tau_a]; ok=kk>=0
r=R["P"][kk[ok]]/R["gP"][kk[ok]]; dg=dgg[ok]; Pl=np.minimum(G["P"][a],G["P"][b])[ok]
for lo,hi in [(0,0.02),(0.02,0.05),(0.05,0.2),(0.2,2)]:
    m=(dg>=lo)&(dg<hi); P("- 1g con dR(gg) gen en [%g,%g): frac %.3f, P_reco/P_gen del conservado p50 %.3f, frac >1.15: %.3f, P gen del perdido p50 %.2f"%(lo,hi,m.mean(),np.median(r[m]),(r[m]>1.15).mean(),np.median(Pl[m])))
# figuras
fig,ax=plt.subplots(1,3,figsize=(15,4.2))
bins=np.linspace(0,0.3,61)
for g,nm in [(1,"gen 1 (rho)"),(2,"gen 2"),(3,"gen 3")]:
    s=three&(gt==g); ax[0].hist(np.minimum(best_d[s],0.299),bins,histtype="step",density=True,label=nm)
ax[0].set_xlabel("|m(mejor pareja) - m_pi0| [GeV], reco 3g+"); ax[0].legend(); ax[0].grid(alpha=.3)
bins=np.linspace(0,2,51)
for g,nm in [(1,"gen 1 (rho)"),(2,"gen 2"),(3,"gen 3")]:
    s=three&(gt==g)&(best_d<0.05); ax[1].hist(np.minimum(restP[s]/np.maximum(pairP[s],1e-6),1.99),bins,histtype="step",density=True,label=nm)
ax[1].set_xlabel("P(resto)/P(pareja), pareja a masa de pi0"); ax[1].legend(); ax[1].grid(alpha=.3)
s=three&(gt==1)&(gpi0P>0); ax[2].hist(np.clip(T["Pgg_all"][s]/gpi0P[s],0,2.5),np.linspace(0,2.5,61),histtype="step",density=True,label="3g+: suma de todos")
ax[2].hist(np.clip(pairP[s]/gpi0P[s],0,2.5),np.linspace(0,2.5,61),histtype="step",density=True,label="3g+: mejor pareja")
ax[2].hist(np.clip(T["Pgg_all"][s2]/gpi0P[s2],0,2.5),np.linspace(0,2.5,61),histtype="step",density=True,label="2g (referencia)")
ax[2].set_xlabel("P(pi0 estimado)/P(pi0 gen), gen 1"); ax[2].legend(); ax[2].grid(alpha=.3); ax[2].set_yscale("log")
plt.tight_layout(); plt.savefig("../figs/recomb_fig1_criterios.png",dpi=120)
np.savez_compressed("a5_tau.npz",best_d=best_d,pairP=pairP,restP=restP,restPmax=restPmax,mgg_all=T["mgg_all"],Pgg_all=T["Pgg_all"],gpi0P=gpi0P)
out.close()

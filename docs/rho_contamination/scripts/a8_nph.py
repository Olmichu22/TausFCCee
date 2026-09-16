"""Refinar la recombinacion en 3g+: numero de fotones y veto de segunda pareja pi0 (gen 2)."""
import numpy as np
from common import *
T,R,G,C=load(); R=sort_by_tau(R); T,R=photon_pairing(T,R); out=open("a8_nph.md","w")
def P(*a): print(*a); print(*a,file=out)
c=cls(T["rtype"]); ti=R["tau"]; nt=len(T["ev"]); gt=T["gtype"]; n=len(ti)
d5=np.load("a5_tau.npz"); best_d=d5["best_d"]; pairP=d5["pairP"]; restP=d5["restP"]; mall=d5["mgg_all"]
three=(c=="3g+"); ratio=restP/np.maximum(pairP,1e-6)
P("## Numero de fotones en 3g+ por tipo gen")
P("| gen | N | 3 | 4 | 5 | >=6 |"); P("|---|---|---|---|---|---|")
for g in (1,2,3):
    s=three&(gt==g); nph=T["nph"][s]; P("| %d | %d | %.3f | %.3f | %.3f | %.3f |"%(g,s.sum(),(nph==3).mean(),(nph==4).mean(),(nph==5).mean(),(nph>=6).mean()))
# segunda pareja: entre los fotones fuera de la mejor pareja, ¿hay otra pareja con |m-mpi0|<0.05? (solo nph>=4)
# reconstruir la mejor pareja por tau (indices) — recomputar rapido
best=np.full(nt,np.inf); bi=np.full(nt,-1); bj=np.full(nt,-1); kmax=int(T["nph"].max())
pairs=[]
for k in range(1,kmax):
    i=np.arange(0,n-k); j=i+k; same=ti[i]==ti[j]; i=i[same]; j=j[same]
    if len(i)==0: break
    cos=(R["px"][i]*R["px"][j]+R["py"][i]*R["py"][j]+R["pz"][i]*R["pz"][j])/(R["P"][i]*R["P"][j])
    m=np.sqrt(np.maximum(2*R["P"][i]*R["P"][j]*(1-np.clip(cos,-1,1)),0)); d=np.abs(m-MPI0); t=ti[i]
    np.minimum.at(best,t,d); sel=d==best[t]; bi[t[sel]]=i[sel]; bj[t[sel]]=j[sel]
    pairs.append((i,j,d))
I=np.concatenate([p[0] for p in pairs]); J=np.concatenate([p[1] for p in pairs]); D=np.concatenate([p[2] for p in pairs]); Tt=ti[I]
# pareja disjunta de la mejor con d<0.05
disj=(I!=bi[Tt])&(I!=bj[Tt])&(J!=bi[Tt])&(J!=bj[Tt])&(D<0.05)
second=np.zeros(nt,bool); second[Tt[disj]]=True
P("\n## Segunda pareja disjunta a masa de pi0 (|m-m_pi0|<0.05), en 3g+ con nph>=4")
for g in (1,2,3):
    s=three&(gt==g)&(T["nph"]>=4); P("- gen %d: N=%d, frac con segunda pareja: %.3f"%(g,s.sum(),second[s].mean()))
# ROC con refinamientos
n1=(gt==1).sum(); n2g=(c=="2g").sum()
def roc(name,extra):
    P("\n### %s"%name); P("| m(all)< | Prest/Ppar< | Δeff gen1->2g | contaminacion/clase 2g | de gen2 | pureza 2g |"); P("|---|---|---|---|---|---|")
    for X,f in [(0.2,0),(0.25,0),(0.3,0),(0.2,0.1),(0.25,0.1),(0.25,0.2),(0.3,0.2),(0.3,0.3),(0.4,0.3)]:
        rc=three&((mall<X)|((best_d<0.05)&(ratio<f)))&extra
        g1=(rc&(gt==1)).sum(); g2=(rc&(gt==2)).sum(); go=(rc&(gt!=1)).sum()
        P("| %g | %g | +%.4f | %.4f | %.4f | %.3f |"%(X,f,g1/n1,go/n2g,g2/n2g,(((c=="2g")&(gt==1)).sum()+g1)/(n2g+rc.sum())))
roc("sin refinamiento",np.ones(nt,bool))
roc("solo nph==3",T["nph"]==3)
roc("veto segunda pareja pi0",~second)
roc("nph==3 o (nph>=4 sin segunda pareja)",(T["nph"]==3)|~second)
# masa pi+all como discriminante adicional entre gen1 y gen2 dentro de los recombinables
rc=three&((mall<0.25)|((best_d<0.05)&(ratio<0.2)))
P("\n## Dentro de los recombinables (m<0.25 | Prest/Ppar<0.2): m(pi+todos los fotones)")
for g in (1,2):
    s=rc&(gt==g)&np.isfinite(T["mtau_all"]); P("- gen %d: p10/50/90 = %.2f/%.2f/%.2f ; frac m<1.0: %.3f; frac m<1.1: %.3f"%(g,*np.percentile(T["mtau_all"][s],[10,50,90]),(T["mtau_all"][s]<1.0).mean(),(T["mtau_all"][s]<1.1).mean()))
for mcut in (1.0,1.1,1.2):
    rc2=rc&(T["mtau_all"]<mcut); g1=(rc2&(gt==1)).sum(); go=(rc2&(gt!=1)).sum()
    P("- + m(pi+all)<%.1f: Δeff +%.4f, contaminacion %.4f, pureza %.3f"%(mcut,g1/n1,go/n2g,(((c=="2g")&(gt==1)).sum()+g1)/(n2g+rc2.sum())))
out.close()

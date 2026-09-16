import numpy as np, sys
from common import *
T,R,G,C=load()
out=open("a1_baseline.md","w")
def P(*a):
    print(*a); print(*a,file=out)
g1=T["gtype"]==1
P("## Migracion gen 1 (pi pi0), fichero completo. N =",g1.sum())
P(table(g1,T,[0,5,10,15,20,25,30,35,40,46],"gvisP"))
P()
P(table(g1,T,[0,20,30,40,44,46],"gP"))
P("\n## Lo mismo para gen 0 (control):")
P(table(T["gtype"]==0,T,[0,10,20,30,46],"gvisP"))
P("\n## Lo mismo para gen 2 (pi 2pi0):")
P(table(T["gtype"]==2,T,[0,10,20,30,46],"gvisP"))
# ---- destino de los fotones gen del pi0, gen 1, por resultado reco
c=cls(T["rtype"]); gi=G["tau"]
P("\n## Destino de cada foton gen del pi0 (gen 1), por resultado reco")
P("fate: 0 = no hay PFO foton enlazado; 1 = PFO foton en el reco tau emparejado; 2 = PFO foton en otro reco tau; 3 = PFO foton fuera de todo tau")
P("| reco | N fotones | no reco | en tau | otro tau | fuera tau | P gen mediana (GeV) | P mediana no reco | dR(g,pi) mediana no reco |")
P("|---|---|---|---|---|---|---|---|---|")
for k in ["2g","1g","0g","3g+","pi->n","lep","nomatch"]:
    s=g1[gi]&(c[gi]==k); N=s.sum()
    f=G["fate"][s]
    nr=s&(G["fate"]==0)
    P("| %s | %d | %.3f | %.3f | %.3f | %.3f | %.2f | %.2f | %.3f |"%(k,N,(f==0).mean(),(f==1).mean(),(f==2).mean(),(f==3).mean(),np.median(G["P"][s]),np.median(G["P"][nr]) if nr.any() else np.nan,np.median(G["dRpi"][nr]) if nr.any() else np.nan))
# combinaciones por tau: (fate foton1, fate foton2)
P("\n## Combinacion de destinos de los dos fotones por tau (gen 1 con 2 fotones gen)")
ng=np.bincount(gi,minlength=len(T["ev"]))
two=g1&(ng==2)
# ordenar fotones por tau: fila del tau, fate, P
idx=np.where(two[gi])[0]; o=np.lexsort((-G["P"][idx],gi[idx])); idx=idx[o]
tau_i=gi[idx][0::2]; f1=G["fate"][idx][0::2]; f2=G["fate"][idx][1::2]; P1=G["P"][idx][0::2]; P2=G["P"][idx][1::2]
assert (gi[idx][1::2]==tau_i).all()
names={0:"noreco",1:"tau",2:"otro",3:"fuera"}
for k in ["2g","1g","0g","3g+"]:
    s=c[tau_i]==k; N=s.sum(); P("\n### reco %s (N=%d)"%(k,N))
    P("| foton lider | foton blando | fraccion | P lider med | P blando med |"); P("|---|---|---|---|---|")
    combos={}
    for a in range(4):
        for b in range(4):
            m=s&(f1==a)&(f2==b)
            if m.sum()>0.002*N: P("| %s | %s | %.3f | %.2f | %.2f |"%(names[a],names[b],m.sum()/N,np.median(P1[m]),np.median(P2[m])))
# ---- fotones reco en el cono, gen 1, por resultado: categoria de origen
P("\n## Origen de los fotones reco del cono (gen 1), por resultado reco")
P("cat: 0 pi0 del mismo tau, 1 FSR linea tau, 2 rad cargada/ISR, 3 otro, 4 sim, 5 sin match gen (fragmento shower), 6 pi0 del otro tau")
ti=R["tau"]
P("| reco | N fotones | pi0 propio | FSR | rad/ISR | otro | sim | sin match | pi0 otro tau | cluster partido (nshare>=2) |")
P("|---|---|---|---|---|---|---|---|---|---|")
for k in ["2g","1g","3g+"]:
    s=g1[ti]&(c[ti]==k); N=s.sum(); cc=R["cat"][s]
    P("| %s | %d | "%(k,N)+" | ".join("%.3f"%((cc==j).mean()) for j in range(7))+" | %.3f |"%((R["nshare"][s]>=2).mean()))
# 3g+: por tau, cuantos fotones de pi0 propio y de que es el resto
P("\n### reco 3g+: composicion por tau")
n_pi0=np.bincount(ti,weights=(R["cat"]==0),minlength=len(T["ev"])); n_frag=np.bincount(ti,weights=(R["cat"]==5),minlength=len(T["ev"]))
n_fsr=np.bincount(ti,weights=(R["cat"]==1),minlength=len(T["ev"])); n_oth=np.bincount(ti,weights=np.isin(R["cat"],[2,3,4,6]),minlength=len(T["ev"]))
# fotones de pi0 propio distintos (por gen match) -> partidos
distinct=np.zeros(len(T["ev"]))
key=np.where(R["cat"]==0,R["ev"]*1000+R["gm"],-1)
u,ii=np.unique(key,return_index=True); ii=ii[u>=0]; np.add.at(distinct,ti[ii],1)
s=g1&(c=="3g+")
P("| n fotones pi0 propio (PFO) | n gen distintos | fragmentos | FSR | otros | fraccion |"); P("|---|---|---|---|---|---|")
combos={}
for a,b,cf,d,e in zip(n_pi0[s],distinct[s],n_frag[s],n_fsr[s],n_oth[s]):
    combos[(a,b,cf,d,e)]=combos.get((a,b,cf,d,e),0)+1
for kk,v in sorted(combos.items(),key=lambda x:-x[1])[:12]:
    P("| %d | %d | %d | %d | %d | %.3f |"%(kk+(v/s.sum(),)))
P("\n3g+ con algun cluster partido (mas PFO foton que fotones gen distintos): %.3f"%((n_pi0[s]>distinct[s]).mean()))
P("3g+ con algun fragmento sin match: %.3f"%((n_frag[s]>0).mean()))
P("3g+ con algun FSR: %.3f"%((n_fsr[s]>0).mean()))
P("3g+ con foton de pi0 de otro tau/otro: %.3f"%((n_oth[s]>0).mean()))
# 1g: composicion
s=g1&(c=="1g")
P("\n### reco 1g: que es el unico foton")
P("| pi0 propio | fragmento | FSR | otros | "); P("|---|---|---|---|")
P("| %.3f | %.3f | %.3f | %.3f |"%((n_pi0[s]==1).mean(),(n_frag[s]==1).mean(),(n_fsr[s]==1).mean(),(n_oth[s]==1).mean()))
out.close()

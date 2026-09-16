import numpy as np
from common import *
T,R,G,C=load(); out=open("a2_purity.md","w")
def P(*a): print(*a); print(*a,file=out)
c=cls(T["rtype"]); gt=T["gtype"]
gcl=np.full(len(gt),"otro",dtype=object)
for v,nm in [(0,"gen0 pi"),(1,"gen1 rho"),(2,"gen2 pi2pi0"),(3,"gen3"),(10,"gen10 3p"),(11,"gen11"),(-11,"e"),(-13,"mu")]: gcl[gt==v]=nm
names=["gen1 rho","gen0 pi","gen2 pi2pi0","gen3","gen10 3p","gen11","e","mu","otro"]
for k in ["2g","1g","3g+","0g"]:
    P("\n## Composicion de la clase reco %s por gen (filas: P visible reco del tau)"%k)
    P("| rP | N | "+" | ".join(names)+" |"); P("|"+"---|"*(len(names)+2))
    for lo,hi in [(0,5),(5,10),(10,20),(20,30),(30,40),(40,50),(0,50)]:
        s=(c==k)&(T["rP"]>=lo)&(T["rP"]<hi); N=s.sum()
        P("| %g-%g | %d | "%(lo,hi,N)+" | ".join("%.3f"%((gcl[s]==nm).mean()) for nm in names)+" |")
# gen1 en 2g: fraccion con los dos fotones reales del pi0 (cat 0) frente a otras combinaciones
ti=R["tau"]; nt=len(T["ev"])
n0=np.bincount(ti,weights=R["cat"]==0,minlength=nt); nfr=np.bincount(ti,weights=R["cat"]==5,minlength=nt); nfs=np.bincount(ti,weights=R["cat"]==1,minlength=nt)
key=np.where(R["cat"]==0,R["ev"]*1000+R["gm"],-1); u,ii=np.unique(key,return_index=True); ii=ii[u>=0]; dist=np.zeros(nt); np.add.at(dist,ti[ii],1)
s=(c=="2g")&(gt==1)
P("\n## gen1 -> 2g: que son los dos fotones")
combos={}
for a,b,f,e in zip(n0[s],dist[s],nfr[s],nfs[s]): combos[(a,b,f,e)]=combos.get((a,b,f,e),0)+1
P("| PFO pi0 propio | gen distintos | fragmentos | FSR | fraccion |"); P("|---|---|---|---|---|")
for kk,v in sorted(combos.items(),key=lambda x:-x[1])[:8]: P("| %d | %d | %d | %d | %.4f |"%(kk+(v/s.sum(),)))
out.close()

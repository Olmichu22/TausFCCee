import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from common import *
T,R,G,C=load(); out=open("a3_lost.md","w")
def P(*a): print(*a); print(*a,file=out)
c=cls(T["rtype"]); gi=G["tau"]; g1=T["gtype"]==1
# angulo gen gamma-gamma por tau y P del pi0
nt=len(T["ev"])
px,py,pz=cart(G["P"],G["th"],G["phi"])
sx=np.bincount(gi,weights=px,minlength=nt); sy=np.bincount(gi,weights=py,minlength=nt); sz=np.bincount(gi,weights=pz,minlength=nt)
T["gpi0P"]=np.sqrt(sx**2+sy**2+sz**2)
T["gpi0th"]=np.arccos(np.clip(sz/np.maximum(T["gpi0P"],1e-9),-1,1)); T["gpi0phi"]=np.arctan2(sy,sx)
T["dRpipi0"]=ang(T["gpith"],T["gpiphi"],T["gpi0th"],T["gpi0phi"])
G["dRpi0"]=ang(G["th"],G["phi"],T["gpi0th"][gi],T["gpi0phi"][gi])
P("## Fotones gen del pi0 (gen 1) no reconstruidos en el tau: donde estan")
P("dR(gamma, pion gen) de cada foton segun destino, por resultado reco. Cono de reco = 0.4.")
P("| reco | destino | N | dR p10 | p50 | p90 | frac dR>0.4 | P p10 | p50 | p90 | frac P<0.5 | frac P<1 |")
P("|---|---|---|---|---|---|---|---|---|---|---|---|")
names={0:"no reco",1:"en tau",2:"otro tau",3:"fuera tau"}
for k in ["0g","1g","2g","3g+"]:
    for f in (0,1,3):
        s=g1[gi]&(c[gi]==k)&(G["fate"]==f); N=s.sum()
        if N<50: continue
        q=np.percentile(G["dRpi"][s],[10,50,90]); qp=np.percentile(G["P"][s],[10,50,90])
        P("| %s | %s | %d | %.3f | %.3f | %.3f | %.3f | %.2f | %.2f | %.2f | %.3f | %.3f |"%(k,names[f],N,*q,(G["dRpi"][s]>0.4).mean(),*qp,(G["P"][s]<0.5).mean(),(G["P"][s]<1).mean()))
P("\n## Por bin de P visible: fraccion de gen1 en que algun foton del pi0 cae fuera del cono dR>0.4 del pion gen")
ng_out=np.bincount(gi,weights=G["dRpi"]>0.4,minlength=nt); ng_out2=np.bincount(gi,weights=G["dRpi"]>0.3,minlength=nt)
P("| gvisP | N gen1 | >=1 foton dR>0.4 | >=1 foton dR>0.3 | reco 0g | reco 0g si foton fuera | reco 2g si todos dentro |")
P("|---|---|---|---|---|---|---|")
for lo,hi in [(0,5),(5,10),(10,15),(15,20),(20,30),(30,46)]:
    s=g1&(T["gvisP"]>=lo)&(T["gvisP"]<hi); o=s&(ng_out>0); i=s&(ng_out==0)
    P("| %g-%g | %d | %.3f | %.3f | %.3f | %.3f | %.3f |"%(lo,hi,s.sum(),(ng_out[s]>0).mean(),(ng_out2[s]>0).mean(),(c[s]=="0g").mean(),(c[o]=="0g").mean() if o.any() else np.nan,(c[i]=="2g").mean()))
# el foton perdido en 1g: P del foton perdido, P del foton que si esta, y su angulo entre si
P("\n## reco 1g: el foton perdido (no reco) frente al reconstruido")
s1=g1&(c=="1g"); two=np.bincount(gi,minlength=nt)==2
idx=np.where((s1&two)[gi])[0]; o=np.lexsort((-G["P"][idx],gi[idx])); idx=idx[o]
a=idx[0::2]; b=idx[1::2]
lost=np.where(G["fate"][a]==1,b,a); kept=np.where(G["fate"][a]==1,a,b)
m=(G["fate"][kept]==1)&(G["fate"][lost]!=1)
dgg=ang(G["th"][a],G["phi"][a],G["th"][b],G["phi"][b])
P("N taus 1g con exactamente un foton del pi0 en el tau: %d de %d"%(m.sum(),len(a)))
for f,nm in [(0,"no reco"),(3,"fuera tau"),(2,"otro tau")]:
    mm=m&(G["fate"][lost]==f)
    if mm.sum()<50: continue
    P("- perdido = %s (%.3f): P perdido p10/50/90 = %.2f/%.2f/%.2f GeV, P conservado = %.2f/%.2f/%.2f, dR(gg) p10/50/90 = %.3f/%.3f/%.3f, dR(perdido,pi) = %.3f/%.3f/%.3f, frac dR(gg)<0.02: %.3f"%(
        nm,mm.sum()/m.sum(),*np.percentile(G["P"][lost][mm],[10,50,90]),*np.percentile(G["P"][kept][mm],[10,50,90]),*np.percentile(dgg[mm],[10,50,90]),*np.percentile(G["dRpi"][lost][mm],[10,50,90]),(dgg[mm]<0.02).mean()))
# comparar con 2g: dR(gg) y P del foton blando
s2=g1&(c=="2g"); idx=np.where((s2&two)[gi])[0]; o=np.lexsort((-G["P"][idx],gi[idx])); idx=idx[o]; a2=idx[0::2]; b2=idx[1::2]
dgg2=ang(G["th"][a2],G["phi"][a2],G["th"][b2],G["phi"][b2])
P("- control 2g: P blando p10/50/90 = %.2f/%.2f/%.2f, dR(gg) = %.3f/%.3f/%.3f, frac dR(gg)<0.02: %.3f"%(*np.percentile(G["P"][b2],[10,50,90]),*np.percentile(dgg2,[10,50,90]),(dgg2<0.02).mean()))
# eficiencia de reconstruir un foton de pi0 en el tau en funcion de su P y de dR al pion (solo taus emparejados 1-prong)
P("\n## Eficiencia de que un foton gen de pi0 acabe como PFO foton en el tau, vs P y dR(g,pi) (gen1, tau emparejado 1-prong)")
mt=g1&np.isin(c,["0g","1g","2g","3g+"]); s=mt[gi]
Pb=[0,0.2,0.3,0.5,0.75,1,1.5,2,3,5,10,50]; Db=[0,0.02,0.05,0.1,0.2,0.3,0.4,1]
P("| P \\ dR | "+" | ".join("%g-%g"%(a,b) for a,b in zip(Db[:-1],Db[1:]))+" |"); P("|"+"---|"*(len(Db)))
H=np.zeros((len(Pb)-1,len(Db)-1)); Hn=np.zeros_like(H)
ip=np.digitize(G["P"],Pb)-1; idd=np.digitize(G["dRpi"],Db)-1
ok=s&(ip>=0)&(ip<len(Pb)-1)&(idd>=0)&(idd<len(Db)-1)
np.add.at(Hn,(ip[ok],idd[ok]),1); np.add.at(H,(ip[ok],idd[ok]),G["fate"][ok]==1)
for i in range(len(Pb)-1):
    P("| %g-%g | "%(Pb[i],Pb[i+1])+" | ".join("%.2f (%d)"%(H[i,j]/max(Hn[i,j],1),Hn[i,j]) for j in range(len(Db)-1))+" |")
# figura: eficiencia vs P para dR<0.1 y dR>0.1, y distribucion P/dR del foton perdido en 1g
fig,ax=plt.subplots(1,3,figsize=(15,4.2))
for lo,hi,lab in [(0,0.05,"dR<0.05"),(0.05,0.15,"0.05<dR<0.15"),(0.15,0.4,"0.15<dR<0.4")]:
    ss=s&(G["dRpi"]>=lo)&(G["dRpi"]<hi); h,_=np.histogram(G["P"][ss],Pb); hh,_=np.histogram(G["P"][ss&(G["fate"]==1)],Pb)
    ax[0].step(Pb[:-1],hh/np.maximum(h,1),where="post",label=lab)
ax[0].set_xscale("log"); ax[0].set_xlabel("P foton gen del pi0 [GeV]"); ax[0].set_ylabel("frac. reconstruido como PFO foton en el tau"); ax[0].legend(); ax[0].grid(alpha=.3)
bins=np.linspace(0,5,51)
ax[1].hist(G["P"][lost][m&(G["fate"][lost]==0)],bins,histtype="step",label="1g: foton perdido (no reco)",density=True)
ax[1].hist(G["P"][lost][m&(G["fate"][lost]==3)],bins,histtype="step",label="1g: foton fuera del tau",density=True)
ax[1].hist(G["P"][b2],bins,histtype="step",label="2g: foton blando (control)",density=True)
ax[1].set_xlabel("P gen del foton blando [GeV]"); ax[1].legend(); ax[1].grid(alpha=.3)
bins=np.linspace(0,0.6,61)
ax[2].hist(G["dRpi"][lost][m&(G["fate"][lost]==0)],bins,histtype="step",label="1g: no reco",density=True)
ax[2].hist(G["dRpi"][lost][m&(G["fate"][lost]==3)],bins,histtype="step",label="1g: fuera del tau",density=True)
ax[2].hist(G["dRpi"][b2],bins,histtype="step",label="2g: foton blando",density=True)
ax[2].axvline(0.4,color="k",ls="--"); ax[2].set_xlabel("dR(foton blando, pion gen)"); ax[2].legend(); ax[2].grid(alpha=.3)
plt.tight_layout(); plt.savefig("../figs/lost_fig1_foton_perdido.png",dpi=120)
out.close()

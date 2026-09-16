import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from common import *
T,R,G,C=load(); R=sort_by_tau(R); T,R=photon_pairing(T,R); out=open("a4_3g.md","w")
def P(*a): print(*a); print(*a,file=out)
c=cls(T["rtype"]); ti=R["tau"]; nt=len(T["ev"]); gt=T["gtype"]
# 1) fragmentos sin match en gen1->3g+: angulo al foton gen lider frente al pion
gi=G["tau"]; px,py,pz=cart(G["P"],G["th"],G["phi"])
# foton gen lider por tau
lead=np.full(nt,-1); o=np.argsort(G["P"]); lead[gi[o]]=o
s=(gt==1)&(c=="3g+")
fr=s[ti]&(R["cat"]==5); L=lead[ti[fr]]
dRg=ang(R["th"][fr],R["phi"][fr],G["th"][L],G["phi"][L]); dRp=R["dRpi"][fr]
P("## gen1 -> 3g+: fotones reco sin match gen (fragmentos), N=%d"%fr.sum())
P("dR al foton gen lider p10/50/90: %.3f/%.3f/%.3f ; dR al pion reco: %.3f/%.3f/%.3f"%(*np.percentile(dRg,[10,50,90]),*np.percentile(dRp,[10,50,90])))
P("mas cerca del foton lider que del pion: %.3f ; P fragmento p10/50/90: %.2f/%.2f/%.2f ; P/Ppi p50: %.3f ; P/P(foton lider gen) p50: %.3f"%((dRg<dRp).mean(),*np.percentile(R["P"][fr],[10,50,90]),np.median(R["fpi"][fr]),np.median(R["P"][fr]/G["P"][L])))
# clusters partidos: dR entre las dos partes
sp=s[ti]&(R["cat"]==0)&(R["nshare"]>=2)
P("fotones de pi0 partidos (nshare>=2): N=%d, P p50 %.2f, P/P_gen p10/50/90 %.2f/%.2f/%.2f"%(sp.sum(),np.median(R["P"][sp]),*np.percentile(R["P"][sp]/R["gP"][sp],[10,50,90])))
# 2) suma de P de fotones frente al P del pi0 gen
sx=np.bincount(gi,weights=px,minlength=nt); sy=np.bincount(gi,weights=py,minlength=nt); sz=np.bincount(gi,weights=pz,minlength=nt); gpi0P=np.sqrt(sx**2+sy**2+sz**2)
P("\n## Balance: P(todos los fotones reco del cono)/P(pi0 gen), gen 1, por resultado")
P("| reco | N | p10 | p50 | p90 | frac >1.15 |"); P("|---|---|---|---|---|---|")
for k in ["2g","1g","3g+"]:
    ss=(gt==1)&(c==k)&(gpi0P>0); r=T["Pgg_all"][ss]/gpi0P[ss]
    P("| %s | %d | %.3f | %.3f | %.3f | %.3f |"%(k,ss.sum(),*np.percentile(r,[10,50,90]),(r>1.15).mean()))
# 3) masas: m(todos los fotones) y m(pi+fotones) en 3g+ por gen type; y en 2g
P("\n## Masas en la clase reco 3g+ por tipo gen")
P("| gen | N | m(gg..) p10/50/90 | frac m(gg..)<0.25 | m(pi+gg..) p10/50/90 | frac mtau<1.0 | P visible reco p50 |"); P("|---|---|---|---|---|---|---|")
for g,nm in [(1,"1 rho"),(2,"2 pi2pi0"),(3,"3"),(0,"0 pi"),(-11,"e")]:
    ss=(gt==g)&(c=="3g+")
    if ss.sum()<100: continue
    P("| %s | %d | %.3f/%.3f/%.3f | %.3f | %.2f/%.2f/%.2f | %.3f | %.1f |"%(nm,ss.sum(),*np.percentile(T["mgg_all"][ss],[10,50,90]),(T["mgg_all"][ss]<0.25).mean(),*np.nanpercentile(T["mtau_all"][ss],[10,50,90]),(T["mtau_all"][ss]<1.0).mean(),np.median(T["rP"][ss])))
P("\n## Masas en la clase reco 2g por tipo gen (referencia)")
P("| gen | N | m(gg) p10/50/90 | frac |m-m_pi0|<0.05 | m(pi+gg) p10/50/90 |"); P("|---|---|---|---|---|")
for g,nm in [(1,"1 rho"),(2,"2 pi2pi0"),(0,"0 pi"),(-11,"e")]:
    ss=(gt==g)&(c=="2g")
    if ss.sum()<100: continue
    P("| %s | %d | %.3f/%.3f/%.3f | %.3f | %.2f/%.2f/%.2f |"%(nm,ss.sum(),*np.percentile(T["mgg_all"][ss],[10,50,90]),(np.abs(T["mgg_all"][ss]-MPI0)<0.05).mean(),*np.nanpercentile(T["mtau_all"][ss],[10,50,90])))
# 4) emulacion del corte pion_photon_fsr actual sobre 3g+ y el resto
paired=np.isfinite(R["dmpi0"])&(R["dmpi0"]<0.05)
one=np.isin(T["rtype"][ti],np.arange(0,10))&(T["nneu"][ti]==0)
drop_soft=one&~paired&(R["fpi"]<0.05); drop_hard=one&~paired&(R["P"]>2)&(R["mpig"]>1.2)
rt_base=T["rtype"].copy(); rt_fsr=emulate(T,R,~(drop_soft|drop_hard)); rt_soft=emulate(T,R,~drop_soft)
P("\n## Emulacion del corte pion_photon_fsr vigente (sin pareja pi0 & (P/Ppi<0.05 | (P>2 & m(pi+g)>1.2)))")
P("Nota: la emulacion asume que la pareja se evalua con todos los fotones del cono y el pion lider = primer constituyente cargado.")
for g,nm in [(1,"gen 1 (rho)"),(0,"gen 0"),(2,"gen 2"),(3,"gen 3")]:
    ss=gt==g
    P("\n### %s: fraccion en cada clase reco"%nm)
    P("| version | 0g | 1g | 2g | 3g+ |"); P("|---|---|---|---|---|")
    for nm2,rt in [("baseline",rt_base),("solo parte blanda",rt_soft),("blanda+dura (vigente)",rt_fsr)]:
        cc=cls(rt)[ss]; P("| %s | %.3f | %.3f | %.3f | %.3f |"%(nm2,(cc=="0g").mean(),(cc=="1g").mean(),(cc=="2g").mean(),(cc=="3g+").mean()))
np.savez_compressed("a4_rt.npz",rt_base=rt_base,rt_fsr=rt_fsr)
# 5) figuras: m(gg..) 3g+ gen1 vs gen2; balance; dR fragmento
fig,ax=plt.subplots(1,3,figsize=(15,4.2))
bins=np.linspace(0,1.2,61)
for g,nm in [(1,"gen 1 (rho)"),(2,"gen 2 (pi 2pi0)"),(0,"gen 0 (pi)")]:
    ss=(gt==g)&(c=="3g+"); ax[0].hist(T["mgg_all"][ss],bins,histtype="step",density=True,label=nm+" N=%d"%ss.sum())
ax[0].set_xlabel("m(todos los fotones del cono) [GeV], reco 3g+"); ax[0].legend(); ax[0].grid(alpha=.3)
bins=np.linspace(0,2.5,61)
for g,nm in [(1,"gen 1 (rho)"),(2,"gen 2 (pi 2pi0)"),(0,"gen 0 (pi)")]:
    ss=(gt==g)&(c=="3g+"); ax[1].hist(T["mtau_all"][ss],bins,histtype="step",density=True,label=nm)
ax[1].set_xlabel("m(pi + todos los fotones) [GeV], reco 3g+"); ax[1].legend(); ax[1].grid(alpha=.3)
bins=np.linspace(0,0.4,41)
ax[2].hist(dRg,bins,histtype="step",density=True,label="dR(fragmento, foton gen lider)")
ax[2].hist(dRp,bins,histtype="step",density=True,label="dR(fragmento, pion reco)")
ax[2].set_xlabel("gen1 -> 3g+: fotones sin match gen"); ax[2].legend(); ax[2].grid(alpha=.3)
plt.tight_layout(); plt.savefig("../figs/3g_fig1_masas_fragmentos.png",dpi=120)
out.close()

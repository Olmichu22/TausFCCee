"""Firma de conversion en 3g+ (fragmentos colineales en theta con el foton lider) y criterio de fusion de fotones.
Tambien: matriz de transicion del corte pion_photon_fsr vigente en gen 1."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from common import *
T,R,G,C=load(); R=sort_by_tau(R); T,R=photon_pairing(T,R); out=open("a6_conv.md","w")
def P(*a): print(*a); print(*a,file=out)
c=cls(T["rtype"]); ti=R["tau"]; nt=len(T["ev"]); gt=T["gtype"]; n=len(ti)
# foton reco lider por tau
lead=np.full(nt,-1); o=np.argsort(R["P"]); lead[ti[o]]=o
isl=np.zeros(n,bool); isl[lead[lead>=0]]=True
L=lead[ti]
dth=R["th"]-R["th"][L]; dph=(R["phi"]-R["phi"][L]+np.pi)%(2*np.pi)-np.pi
R["dth_lead"]=dth; R["dph_lead"]=dph
P("## Fotones no lider del cono respecto al foton lider: |Δθ| y |Δφ|, por origen (gen 1, reco 3g+)")
P("| origen | N | |Δθ| p10/50/90 | |Δφ| p10/50/90 | frac |Δθ|<0.01 | frac |Δθ|<0.02 | P p50 | P/P_lider p50 |")
P("|---|---|---|---|---|---|---|---|")
s3=(gt[ti]==1)&(c[ti]=="3g+")&~isl
for cat,nm in [(0,"pi0 propio (2o foton real o partido)"),(5,"sin match (fragmento conv/shower)"),(1,"FSR"),(6,"pi0 otro tau"),(2,"rad/ISR")]:
    s=s3&(R["cat"]==cat)
    if s.sum()<50: continue
    P("| %s | %d | %.4f/%.4f/%.4f | %.4f/%.4f/%.4f | %.3f | %.3f | %.2f | %.3f |"%(nm,s.sum(),*np.percentile(np.abs(dth[s]),[10,50,90]),*np.percentile(np.abs(dph[s]),[10,50,90]),(np.abs(dth[s])<0.01).mean(),(np.abs(dth[s])<0.02).mean(),np.median(R["P"][s]),np.median(R["P"][s]/R["P"][L][s])))
# sin match: partidos (nshare>=2 del pi0 propio) por separado
s=s3&(R["cat"]==0)&(R["nshare"]>=2); P("- pi0 propio con cluster partido (nshare>=2): N=%d, |Δθ| p50 %.4f, |Δφ| p50 %.4f, frac |Δθ|<0.01: %.3f"%(s.sum(),np.median(np.abs(dth[s])),np.median(np.abs(dph[s])),(np.abs(dth[s])<0.01).mean()))
s=s3&(R["cat"]==0)&(R["nshare"]<2); P("- pi0 propio segundo foton real: N=%d, |Δθ| p50 %.4f, |Δφ| p50 %.4f, frac |Δθ|<0.01: %.3f"%(s.sum(),np.median(np.abs(dth[s])),np.median(np.abs(dph[s])),(np.abs(dth[s])<0.01).mean()))
P("\n## Referencia: 2o foton en reco 2g y fotones no lider en gen 2 3g+")
s=(gt[ti]==1)&(c[ti]=="2g")&~isl; P("- gen1 2g, foton blando: |Δθ| p10/50/90 %.4f/%.4f/%.4f, frac<0.01: %.3f, frac<0.02: %.3f"%(*np.percentile(np.abs(dth[s]),[10,50,90]),(np.abs(dth[s])<0.01).mean(),(np.abs(dth[s])<0.02).mean()))
s=(gt[ti]==2)&(c[ti]=="3g+")&~isl; P("- gen2 3g+, no lider: |Δθ| p10/50/90 %.4f/%.4f/%.4f, frac<0.01: %.3f, frac<0.02: %.3f"%(*np.percentile(np.abs(dth[s]),[10,50,90]),(np.abs(dth[s])<0.01).mean(),(np.abs(dth[s])<0.02).mean()))
# --- criterio de fusion: fotones con |Δθ|<d respecto a CUALQUIER otro foton mas energetico del cono y |Δφ|<f se fusionan con el.
# implementacion: para cada foton no lider, buscar el foton mas energetico del tau con |Δθ|<d & |Δφ|<f; si existe se absorbe.
def merge_count(d,f):
    absorbed=np.zeros(n,bool)
    kmax=int(T["nph"].max())
    for k in range(1,kmax):
        i=np.arange(0,n-k); j=i+k; same=ti[i]==ti[j]; i=i[same]; j=j[same]
        if len(i)==0: break
        a=np.abs(R["th"][i]-R["th"][j])<d; b=np.abs((R["phi"][i]-R["phi"][j]+np.pi)%(2*np.pi)-np.pi)<f
        m=a&b
        lo=np.where(R["P"][i]<R["P"][j],i,j)  # el menos energetico se absorbe
        absorbed[lo[m]]=True
    return absorbed
P("\n## Criterio de fusion angular: un foton se absorbe en otro mas energetico del cono si |Δθ|<d y |Δφ|<f")
n1=(gt==1).sum(); n2g=(c=="2g").sum()
P("| d | f | gen1 2g | gen1 1g | gen1 3g+ | gen2 2g | gen2 3g+ | gen3 3g+ | gen0 0g | gen0 1g | pureza 2g |"); P("|---|---|---|---|---|---|---|---|---|---|---|")
def row(nm,rt):
    cc=cls(rt); s1=gt==1; s2=gt==2; s3_=gt==3; s0=gt==0
    P("| %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |"%(nm,(cc[s1]=="2g").mean(),(cc[s1]=="1g").mean(),(cc[s1]=="3g+").mean(),(cc[s2]=="2g").mean(),(cc[s2]=="3g+").mean(),(cc[s3_]=="3g+").mean(),(cc[s0]=="0g").mean(),(cc[s0]=="1g").mean(),((cc=="2g")&s1).sum()/max((cc=="2g").sum(),1)))
row("baseline | -",T["rtype"])
res={}
for d in (0.005,0.01,0.02,0.03):
    for f in (0.05,0.1,0.2,0.4):
        ab=merge_count(d,f); rt=emulate(T,R,~ab); res[(d,f)]=ab; row("%g | %g"%(d,f),rt)
# combinado con el corte vigente
paired=np.isfinite(R["dmpi0"])&(R["dmpi0"]<0.05); one=np.isin(T["rtype"][ti],np.arange(0,10))&(T["nneu"][ti]==0)
drop_soft=one&~paired&(R["fpi"]<0.05); drop_hard=one&~paired&(R["P"]>2)&(R["mpig"]>1.2)
P("\n### combinaciones con el corte vigente (pion_photon_fsr)")
P("| version | gen1 2g | gen1 1g | gen1 3g+ | gen2 2g | gen2 3g+ | gen3 3g+ | gen0 0g | gen0 1g | pureza 2g |"); P("|---|---|---|---|---|---|---|---|---|---|")
row("vigente (blanda+dura)",emulate(T,R,~(drop_soft|drop_hard)))
row("solo dura",emulate(T,R,~drop_hard))
for d,f in [(0.01,0.2),(0.02,0.2),(0.02,0.4)]:
    ab=res[(d,f)]
    row("fusion d=%g f=%g + dura"%(d,f),emulate(T,R,~(ab|drop_hard)))
    row("fusion d=%g f=%g + vigente"%(d,f),emulate(T,R,~(ab|drop_soft|drop_hard)))
# --- matriz de transicion del corte vigente en gen 1
P("\n## Matriz de transicion del corte vigente, gen 1 (filas baseline, columnas tras el corte; fraccion de todos los gen 1)")
rt=emulate(T,R,~(drop_soft|drop_hard)); cb=cls(T["rtype"]); ca=cls(rt); s1=gt==1
P("| base \\ tras | 0g | 1g | 2g | 3g+ |"); P("|---|---|---|---|---|")
for a in ["0g","1g","2g","3g+"]:
    P("| %s | "%a+" | ".join("%.4f"%((s1&(cb==a)&(ca==b)).sum()/n1) for b in ["0g","1g","2g","3g+"])+" |")
P("\nSolo parte blanda:"); rt=emulate(T,R,~drop_soft); ca=cls(rt)
P("| base \\ tras | 0g | 1g | 2g | 3g+ |"); P("|---|---|---|---|---|")
for a in ["0g","1g","2g","3g+"]:
    P("| %s | "%a+" | ".join("%.4f"%((s1&(cb==a)&(ca==b)).sum()/n1) for b in ["0g","1g","2g","3g+"])+" |")
# los 1g que la parte blanda manda a 0g: que foton era y m(pi+g)
s=s1&(cb=="1g")&(ca=="0g"); ph=(R["cat"])[lead[s]]
P("\n1g -> 0g por la parte blanda: N=%d, el foton es pi0 propio en %.3f, fragmento %.3f, FSR %.3f; m(pi+g) p10/50/90 = %.2f/%.2f/%.2f; P foton p50 %.2f; P/Ppi p50 %.3f"%(s.sum(),(ph==0).mean(),(ph==5).mean(),(ph==1).mean(),*np.percentile(R["mpig"][lead[s]],[10,50,90]),np.median(R["P"][lead[s]]),np.median(R["fpi"][lead[s]])))
# figura: Δθ vs Δφ por origen
fig,ax=plt.subplots(1,3,figsize=(15,4.2))
for cat,nm in [(0,"pi0 propio"),(5,"sin match (conv./shower)"),(1,"FSR")]:
    s=s3&(R["cat"]==cat); ax[0].hist(np.clip(np.abs(dth[s]),0,0.1),np.linspace(0,0.1,51),histtype="step",density=True,label=nm)
s=(gt[ti]==2)&(c[ti]=="3g+")&~isl; ax[0].hist(np.clip(np.abs(dth[s]),0,0.1),np.linspace(0,0.1,51),histtype="step",density=True,label="gen 2 3g+, todos",ls="--")
ax[0].set_xlabel("|Δθ| respecto al foton lider [rad]"); ax[0].legend(); ax[0].grid(alpha=.3); ax[0].set_yscale("log")
for cat,nm in [(0,"pi0 propio"),(5,"sin match (conv./shower)"),(1,"FSR")]:
    s=s3&(R["cat"]==cat); ax[1].hist(np.clip(np.abs(dph[s]),0,0.4),np.linspace(0,0.4,51),histtype="step",density=True,label=nm)
s=(gt[ti]==2)&(c[ti]=="3g+")&~isl; ax[1].hist(np.clip(np.abs(dph[s]),0,0.4),np.linspace(0,0.4,51),histtype="step",density=True,label="gen 2 3g+, todos",ls="--")
ax[1].set_xlabel("|Δφ| respecto al foton lider [rad]"); ax[1].legend(); ax[1].grid(alpha=.3)
s=s3&(R["cat"]==5); ax[2].scatter(dph[s][:20000],dth[s][:20000],s=1,alpha=.3,label="sin match (gen1 3g+)")
s=(gt[ti]==2)&(c[ti]=="3g+")&~isl&(R["cat"]==0); ax[2].scatter(dph[s][:20000],dth[s][:20000],s=1,alpha=.3,label="gen2 3g+, fotones de pi0",color="C3")
ax[2].set_xlim(-0.4,0.4); ax[2].set_ylim(-0.2,0.2); ax[2].set_xlabel("Δφ"); ax[2].set_ylabel("Δθ"); ax[2].legend(markerscale=8); ax[2].grid(alpha=.3)
plt.tight_layout(); plt.savefig("../figs/conv_fig1_dtheta_dphi.png",dpi=120)
out.close()

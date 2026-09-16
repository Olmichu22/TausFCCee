import sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d=np.load(sys.argv[1]); OUT=open("RESULTS.md","w")
def P(*a):
    print(*a); print(*a,file=OUT)
plt.rcParams.update({"figure.facecolor":"white","axes.facecolor":"white","savefig.dpi":130,"axes.grid":True,"grid.alpha":0.3,"font.size":10})
COL={"FSR tau (shower)":"#0072B2","ISR (haz)":"#E69F00","sin match gen":"#009E73","pi0 otro tau":"#CC79A7","FSR desint. (propio)":"#D55E00","otro":"#999999","pi0 propio (gen pi pi0)":"#000000"}
NM_EN={"FSR tau (shower)":"tau FSR (shower)","ISR (haz)":"ISR (beam)","sin match gen":"no gen match",
       "pi0 otro tau":"other tau's pi0","FSR desint. (propio)":"decay FSR (own)","otro":"other",
       "pi0 propio (gen pi pi0)":"own pi0 (gen pi pi0)"}

ph={k[3:]:d[k] for k in d.files if k.startswith("ph_")}
tau={k[4:]:d[k] for k in d.files if k.startswith("tau_")}
NEV=ph["ievt"].max()+1
P(f"# Origen del foton extra en gen tau->pi nu reconstruido como pi+gamma\n\nMuestra ztt_2M_smearing, {NEV} eventos procesados (fichero completo). Tau reco: cono dR<0.4, P(foton)>0.5 GeV (genminP).\n")

# ---------- clasificacion de cada foton reco asignado
def classify(ph):
    o=ph["gorig"]; tk=ph["gtk"]; own=ph["isown"]; par=np.abs(ph["gpar"]); gmi=ph["gmi"]
    cat=np.full(len(o),"otro",dtype=object)
    cat[gmi<0]="sin match gen"
    cat[(o==1)&own]="FSR desint. (propio)"
    cat[(o==1)&(tk==-1)]="FSR tau (shower)"
    cat[(o==1)&(tk>=0)&~own]="FSR tau desint. otro tau"
    cat[(o==2)&(par==11)]="ISR (haz)"
    cat[(o==2)&(par==211)]="rad. pion"
    cat[(o==0)&own]="pi0 propio"
    cat[(o==0)&~own]="pi0 otro tau"
    cat[(o==3)]="otro (eta...)"
    cat[(o==4)]="simulacion"
    return cat
cat=classify(ph)
CATS=["FSR tau (shower)","ISR (haz)","sin match gen","pi0 otro tau","FSR desint. (propio)","FSR tau desint. otro tau","rad. pion","pi0 propio","otro (eta...)","simulacion","otro"]

def table(mask,label,binvar=None,bins=None,binlabel=""):
    P(f"\n### {label}\n")
    cats=[c for c in CATS if ((cat==c)&mask).sum()>0]
    if binvar is None:
        n=mask.sum(); P("| N | "+" | ".join(cats)+" |"); P("|---|"+"---|"*len(cats))
        P(f"| {n} | "+" | ".join(f"{((cat==c)&mask).sum()/n:.3f}" for c in cats)+" |")
        return
    P(f"| {binlabel} | N | "+" | ".join(cats)+" |"); P("|---|---|"+"---|"*len(cats))
    for lo,hi in zip(bins[:-1],bins[1:]):
        m=mask&(binvar>=lo)&(binvar<hi); n=m.sum()
        if n==0: continue
        P(f"| {lo:g}-{hi:g} | {n} | "+" | ".join(f"{((cat==c)&m).sum()/n:.3f}" for c in cats)+" |")

P("## 1-2. Clasificacion del foton asignado (gen tau->pi nu, RecoTauType==1)\n")
P("Categorias: `FSR tau (shower)` = GenPhotonOrigin 1 con GenPhotonTauKey -1 (foton emitido por una copia Pythia del tau antes de la desintegracion; se asigna al tau por angulo, dR<0.4); "
  "`FSR desint. (propio)` = origen 1 y TauKey == tau propio (foton hijo directo del tau en la desintegracion, tau->pi nu gamma); "
  "`ISR (haz)` = origen 2 con parent |PDG|==11; `rad. pion` = origen 2 con parent 211; `pi0 otro tau` = origen 0 de otro tau; "
  "`sin match gen` = RecoPhotonGenMatchIdx==-1 (el PFO no enlaza a ningun foton del generador: fragmento de shower del pion, cluster partido o link al propio pion).\n")
m0=(ph["gtype"]==0)&(ph["rtype"]==1)
table(m0,"Global, GenTauType==0, RecoTauType==1")
for mode,name in [(10,"pi nu puro (TrueMode 10)"),(20,"K nu (TrueMode 20)"),(23,"K0 pi (TrueMode 23)")]:
    mm=m0&(ph["gmode"]==mode)
    if mm.sum()>0: table(mm,f"GenTauType==0, RecoTauType==1, {name}")
mm=m0&~np.isin(ph["gmode"],[10,20,23])
if mm.sum()>0:
    import collections; P("\nOtros TrueMode en type0/reco1: "+str(collections.Counter(ph["gmode"][mm]).most_common(6)))
    table(mm,"GenTauType==0, RecoTauType==1, resto de TrueMode")
table(m0&(ph["gmode"]==10),"TrueMode 10, RecoTauType==1, por GenVisTauP",ph["gvisP"],np.arange(0,55,5),"GenVisTauP (GeV)")
table(m0&(ph["gmode"]==10),"TrueMode 10, RecoTauType==1, por |cos theta| del tau",np.abs(ph["gcos"]),np.array([0,0.2,0.4,0.6,0.7,0.75,0.8,0.85,0.9,0.95,1.0]),"|cos theta|")
# fraccion de migracion por cos theta
P("\n### Fraccion de migracion 0->1 (gen pi nu puro) por |cos theta| del tau (respecto a los reconstruidos con tipo 0 o 1)\n")
t10=tau["gmode"]==10
P("| bin cos | N(reco0+reco1) | reco1/(reco0+reco1) | reco1 sin match gen | reco1 FSR shower | reco1 ISR |"); P("|---|---|---|---|---|---|")
cb=np.array([0,0.2,0.4,0.6,0.7,0.75,0.8,0.85,0.9,0.95,1.0])
mig_cos=[]
for lo,hi in zip(cb[:-1],cb[1:]):
    mt=t10&(np.abs(tau["gcos"])>=lo)&(np.abs(tau["gcos"])<hi)&np.isin(tau["rtype"],[0,1]); n=mt.sum()
    n1=(mt&(tau["rtype"]==1)).sum()
    mp=m0&(ph["gmode"]==10)&(np.abs(ph["gcos"])>=lo)&(np.abs(ph["gcos"])<hi)
    f=lambda c: ((cat==c)&mp).sum()/n
    P(f"| {lo:g}-{hi:g} | {n} | {n1/n:.4f} | {f('sin match gen'):.4f} | {f('FSR tau (shower)'):.4f} | {f('ISR (haz)'):.4f} |")
    mig_cos.append((0.5*(lo+hi),n1/n,f('sin match gen'),f('FSR tau (shower)'),f('ISR (haz)')))
P("\n### Fraccion de migracion 0->1 por GenVisTauP (gen pi nu puro)\n")
P("| P (GeV) | N(reco0+reco1) | reco1/(reco0+reco1) | sin match gen | FSR shower | ISR | pi0 otro tau |"); P("|---|---|---|---|---|---|---|")
pb=np.arange(0,55,5); mig_p=[]
for lo,hi in zip(pb[:-1],pb[1:]):
    mt=t10&(tau["gvisP"]>=lo)&(tau["gvisP"]<hi)&np.isin(tau["rtype"],[0,1]); n=mt.sum()
    if n==0: continue
    n1=(mt&(tau["rtype"]==1)).sum()
    mp=m0&(ph["gmode"]==10)&(ph["gvisP"]>=lo)&(ph["gvisP"]<hi)
    f=lambda c: ((cat==c)&mp).sum()/n
    P(f"| {lo:g}-{hi:g} | {n} | {n1/n:.4f} | {f('sin match gen'):.4f} | {f('FSR tau (shower)'):.4f} | {f('ISR (haz)'):.4f} | {f('pi0 otro tau'):.4f} |")
    mig_p.append((0.5*(lo+hi),n1/n,f('sin match gen'),f('FSR tau (shower)'),f('ISR (haz)'),f('pi0 otro tau')))
# figura de migracion
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
mp_=np.array(mig_p); mc_=np.array(mig_cos)
for a_,arr,xl in [(ax[0],mp_,"gen visible tau P (GeV)"),(ax[1],mc_,"tau |cos theta|")]:
    a_.plot(arr[:,0],arr[:,1],"k-o",lw=2,ms=5,label="total reco1/(reco0+reco1)")
    for i,(nm,c) in enumerate([("sin match gen","#009E73"),("FSR tau (shower)","#0072B2"),("ISR (haz)","#E69F00")]):
        a_.plot(arr[:,0],arr[:,2+i],"-s",color=c,lw=1.6,ms=4,label=NM_EN[nm])
    a_.set_xlabel(xl); a_.set_ylabel("fraction of gen tau->pi nu"); a_.legend(fontsize=8)
ax[0].set_title("pi -> pi+gamma migration vs P"); ax[1].set_title("pi -> pi+gamma migration vs |cos theta|")
fig.tight_layout(); fig.savefig("fig_migracion_P_cos.png"); plt.close(fig)

# ---------- 3. cinematica
P("\n## 3. Cinematica del foton extra por categoria, comparada con fotones legitimos de pi0 (gen pi pi0, reco 1 o 2)\n")
m1=(ph["gtype"]==1)&np.isin(ph["rtype"],[1,2])&(cat=="pi0 propio")
sets=[("FSR tau (shower)",m0&(cat=="FSR tau (shower)")),("ISR (haz)",m0&(cat=="ISR (haz)")),("sin match gen",m0&(cat=="sin match gen")),("pi0 otro tau",m0&(cat=="pi0 otro tau")),("pi0 propio (gen pi pi0)",m1)]
ratio=ph["phP"]/np.maximum(ph["piP"],1e-3)
P("| categoria | N | P_gamma med (GeV) | P_gamma p10/p90 | P_g/P_pi med | dR(g,pi) med | dR p10/p90 | m(pi+g) med (GeV) | frac dR<0.05 | frac P_g<1 GeV | frac P_g<2 GeV |"); P("|---|---|---|---|---|---|---|---|---|---|---|")
for nm,m in sets:
    if m.sum()==0: continue
    pg=ph["phP"][m]; dr=ph["dR"][m]; r=ratio[m]; ms=ph["rmass"][m]
    P(f"| {nm} | {m.sum()} | {np.median(pg):.2f} | {np.percentile(pg,10):.2f}/{np.percentile(pg,90):.2f} | {np.median(r):.3f} | {np.median(dr):.3f} | {np.percentile(dr,10):.3f}/{np.percentile(dr,90):.3f} | {np.median(ms):.3f} | {(dr<0.05).mean():.3f} | {(pg<1).mean():.3f} | {(pg<2).mean():.3f} |")
# mass only for reco1 for pi0 set
m1r1=(ph["gtype"]==1)&(ph["rtype"]==1)&(cat=="pi0 propio")
P(f"\nMasa reco pi+gamma (solo RecoTauType==1) para pi0 legitimo: med {np.median(ph['rmass'][m1r1]):.3f} GeV (N={m1r1.sum()}); para extra en pi nu: med {np.median(ph['rmass'][m0]):.3f} GeV")

vars_=[("phP","photon P (GeV)",np.linspace(0,30,61),False),("ratio","P_gamma / P_pion",np.linspace(0,1.5,61),False),("dR","dR(photon, pion) [rad]",np.linspace(0,0.4,41),False),("rmass","m(pi+gamma) (GeV)",np.linspace(0,2.0,61),False),("phPlog","photon P (GeV), log axis",np.logspace(-0.3,1.7,41),True)]
fig,axs=plt.subplots(2,3,figsize=(14,8)); axs=axs.ravel()
for i,(v,xl,b,lg) in enumerate(vars_):
    ax=axs[i]
    for nm,m in sets:
        if m.sum()<20: continue
        x={"ratio":ratio,"phPlog":ph["phP"]}.get(v,ph.get(v))[m]
        if v=="rmass": x=x[ph["rtype"][m]==1] if nm.startswith("pi0 propio") else x
        ax.hist(x,bins=b,density=True,histtype="step",lw=1.8,color=COL[nm],label=f"{NM_EN[nm]} (N={len(x)})")
    ax.set_xlabel(xl); ax.set_ylabel("density (norm. to 1)"); ax.legend(fontsize=7)
    if lg: ax.set_xscale("log")
axs[0].set_title("Extra photon in pi nu vs pi0 photon"); axs[5].axis("off")
fig.tight_layout(); fig.savefig("fig_cinematica_1d.png"); plt.close(fig)

# 2D P vs dR
fig,axs=plt.subplots(1,4,figsize=(17,4.2))
for ax,(nm,m) in zip(axs,[s for s in sets if s[0]!="pi0 otro tau"]):
    ax.hist2d(ph["dR"][m],ph["phP"][m],bins=[np.linspace(0,0.4,40),np.logspace(-0.3,1.7,40)],cmap="Blues",cmin=1)
    ax.set_yscale("log"); ax.set_xlabel("dR(photon, pion) [rad]"); ax.set_ylabel("photon P (GeV)"); ax.set_title(f"{NM_EN[nm]} (N={m.sum()})",fontsize=10)
fig.tight_layout(); fig.savefig("fig_P_vs_dR_2d.png"); plt.close(fig)

# cortes sencillos: fraccion de cada categoria que sobrevive
P("\n### Cortes sencillos sobre el foton: fraccion que SOBREVIVE en cada categoria (extra en pi nu vs pi0 legitimo en pi pi0)\n")
mass=ph["rmass"]
cuts=[("P_g > 1 GeV",ph["phP"]>1),("P_g > 2 GeV",ph["phP"]>2),("dR < 0.2",ph["dR"]<0.2),("dR < 0.15",ph["dR"]<0.15),("dR > 0.02",ph["dR"]>0.02),
      ("P_g/P_pi > 0.05",ratio>0.05),("P_g/P_pi > 0.1",ratio>0.1),
      ("P_g>1 y dR<0.2",(ph["phP"]>1)&(ph["dR"]<0.2)),("P_g>1 y dR<0.2 y dR>0.02",(ph["phP"]>1)&(ph["dR"]<0.2)&(ph["dR"]>0.02)),
      ("P_g/P_pi>0.05 y dR<0.2",(ratio>0.05)&(ph["dR"]<0.2))]
names=[s_[0] for s_ in sets if s_[0]!="pi0 otro tau"]
P("| corte | "+" | ".join(names)+" | migracion 0->1 residual (pi nu puro) |"); P("|---|"+"---|"*(len(names)+1))
nt10=(t10&np.isin(tau["rtype"],[0,1])).sum()
for cn,cm in cuts:
    vals=[]
    for nm,m in sets:
        if nm=="pi0 otro tau": continue
        vals.append(f"{(cm&m).sum()/m.sum():.3f}")
    resid=(cm&m0&(ph["gmode"]==10)).sum()/nt10
    P(f"| {cn} | "+" | ".join(vals)+f" | {resid:.4f} |")
P(f"\n(migracion 0->1 sin cortes, pi nu puro: {(m0&(ph['gmode']==10)).sum()/nt10:.4f}; la columna 'pi0 propio' es la fraccion de fotones legitimos de pi0 conservados, evaluada sobre gen pi pi0 con reco 1 o 2)")

# ---------- 4. dependencia con P
P("\n## 4. Dependencia con P: el foton extra escala con el pion?\n")
P("| categoria | bin P_tau | N | P_gamma med | P_g/P_pi med | dR med | corr(P_g,P_pi) |"); P("|---|---|---|---|---|---|---|")
fig,axs=plt.subplots(1,3,figsize=(14,4.2))
for ax,(nm,m) in zip(axs,sets[:3]):
    for lo,hi,c in [(5,15,"#56B4E9"),(20,30,"#0072B2"),(35,50,"#D55E00")]:
        mm=m&(ph["gvisP"]>=lo)&(ph["gvisP"]<hi)
        if mm.sum()<10: continue
        cc=np.corrcoef(ph["phP"][mm],ph["piP"][mm])[0,1] if mm.sum()>2 else np.nan
        P(f"| {nm} | {lo}-{hi} | {mm.sum()} | {np.median(ph['phP'][mm]):.2f} | {np.median(ratio[mm]):.3f} | {np.median(ph['dR'][mm]):.3f} | {cc:.2f} |")
        ax.hist(ph["phP"][mm],bins=np.logspace(-0.3,1.7,31),density=True,histtype="step",lw=1.8,color=c,label=f"P_tau {lo}-{hi} GeV (N={mm.sum()})")
    ax.set_xscale("log"); ax.set_xlabel("photon P (GeV)"); ax.set_ylabel("density"); ax.set_title(NM_EN[nm]); ax.legend(fontsize=8)
fig.suptitle("Extra photon P in three tau-P ranges"); fig.tight_layout(); fig.savefig("fig_Pgamma_por_Ptau.png"); plt.close(fig)
# ratio distributions
fig,axs=plt.subplots(1,3,figsize=(14,4.2))
for ax,(nm,m) in zip(axs,sets[:3]):
    for lo,hi,c in [(5,15,"#56B4E9"),(20,30,"#0072B2"),(35,50,"#D55E00")]:
        mm=m&(ph["gvisP"]>=lo)&(ph["gvisP"]<hi)
        if mm.sum()<10: continue
        ax.hist(ratio[mm],bins=np.linspace(0,1,41),density=True,histtype="step",lw=1.8,color=c,label=f"P_tau {lo}-{hi} GeV")
    ax.set_xlabel("P_gamma / P_pion"); ax.set_ylabel("density"); ax.set_title(NM_EN[nm]); ax.legend(fontsize=8)
fig.suptitle("P_gamma/P_pion ratio by tau-P range"); fig.tight_layout(); fig.savefig("fig_ratio_por_Ptau.png"); plt.close(fig)

# balance de momento: es el foton un trozo del pion?
P("\n### Balance de momento: P_pion/P_vis_gen y (P_pion+P_gamma)/P_vis_gen (mediana y p10/p90)\n")
P("| categoria | bin P_tau | N | P_pi/P_gen med (p10/p90) | (P_pi+P_g)/P_gen med (p10/p90) | frac P_pi/P_gen<0.8 |"); P("|---|---|---|---|---|---|")
fig,axs=plt.subplots(1,2,figsize=(11,4.2))
for nm,m in sets[:3]+[sets[4]]:
    for lo,hi in [(5,15),(20,30),(35,50)]:
        mm=m&(ph["gvisP"]>=lo)&(ph["gvisP"]<hi)
        if mm.sum()<10: continue
        a1=ph["piP"][mm]/ph["gvisP"][mm]; a2=(ph["piP"][mm]+ph["phP"][mm])/ph["gvisP"][mm]
        P(f"| {nm} | {lo}-{hi} | {mm.sum()} | {np.median(a1):.3f} ({np.percentile(a1,10):.2f}/{np.percentile(a1,90):.2f}) | {np.median(a2):.3f} ({np.percentile(a2,10):.2f}/{np.percentile(a2,90):.2f}) | {(a1<0.8).mean():.3f} |")
    mm=m&(ph["gvisP"]>=20)
    axs[0].hist(ph["piP"][mm]/ph["gvisP"][mm],bins=np.linspace(0,1.3,53),density=True,histtype="step",lw=1.8,color=COL[nm],label=NM_EN[nm])
    axs[1].hist((ph["piP"][mm]+ph["phP"][mm])/ph["gvisP"][mm],bins=np.linspace(0,1.3,53),density=True,histtype="step",lw=1.8,color=COL[nm],label=NM_EN[nm])
axs[0].set_xlabel("reco P_pion / gen visible tau P"); axs[1].set_xlabel("reco (P_pion + P_gamma) / gen visible tau P")
for a_ in axs: a_.set_ylabel("density"); a_.legend(fontsize=8); a_.set_yscale("log")
fig.suptitle("Momentum balance (gen P_tau > 20 GeV)"); fig.tight_layout(); fig.savefig("fig_balance_momento.png"); plt.close(fig)

# ---------- 5. reco2 y fisica irreducible
P("\n## 5. RecoTauType==2 y fotones gen legitimos en gen tau->pi nu\n")
m2=(ph["gtype"]==0)&(ph["rtype"]==2)
table(m2,"GenTauType==0, RecoTauType==2 (dos fotones, cada foton una fila)")
# pares en reco2: combinaciones de categorias
import collections
key=list(zip(ph["ievt"][m2],ph["itau"][m2])); cc=collections.defaultdict(list)
for k,c in zip(key,cat[m2]): cc[k].append(c)
comb=collections.Counter(tuple(sorted(v)) for v in cc.values())
P("\nCombinaciones de origen en reco2 (N taus = %d):"%len(cc)); P("| combinacion | N | frac |"); P("|---|---|---|")
for k,v in comb.most_common(8): P(f"| {' + '.join(k)} | {v} | {v/len(cc):.3f} |")
table(m2&(cat=="pi0 otro tau"),"reco2, fotones pi0 de otro tau (comprobacion)")
P("\n### Fotones gen fisicos en el cono del gen tau->pi nu (fichero completo, todos los reco)\n")
t0=tau["gtype"]==0
def frac(m,lab):
    P(f"| {lab} | {m.sum()} | {m.sum()/t0.sum():.5f} |")
P(f"N gen taus type 0: {t0.sum()}; TrueMode 10: {(t0&(tau['gmode']==10)).sum()}\n")
P("| clase | N taus | frac de gen type0 |"); P("|---|---|---|")
frac(t0&(tau["nfsrdec"]>0),"con FSR de desintegracion (origen 1, TauKey propio) - cualquier P")
frac(t0&(tau["nfsrsh"]>0),"con FSR de shower del tau en dR<0.4, P>0.5 GeV")
frac(t0&(tau["nfsrsh_reco_tau"]>0),"  ... y ese foton reconstruido y asignado al reco tau")
frac(t0&(tau["nisr"]>0),"con ISR (parent e) en dR<0.4, P>0.5 GeV")
frac(t0&(tau["nisr_reco_tau"]>0),"  ... y reconstruido y asignado al reco tau")
frac(t0&(tau["npirad"]>0),"con radiacion del pion (origen 2, parent 211) en dR<0.4, P>0.5")
frac(t0&(tau["nfsrsh"]>0)&(tau["rtype"]==0),"con FSR shower P>0.5 en cono pero reco tipo 0 (foton perdido/no asignado)")
frac(t0&(tau["nfsrsh"]>0)&(tau["rtype"]==1),"con FSR shower P>0.5 en cono y reco tipo 1")
frac(t0&(tau["rtype"]==1),"reco tipo 1 (total)")
frac(t0&(tau["rtype"]==1)&(tau["nfsrsh"]==0)&(tau["nisr"]==0),"reco tipo 1 sin ningun foton gen P>0.5 en el cono")
# eficiencia de reconstruir el FSR shower en funcion de su P
P("\nEficiencia de que el FSR de shower (en cono, P>0.5) acabe asignado al reco tau, por P del foton gen:\n")
P("| P_gamma gen (GeV) | N | frac asignada al reco tau | frac reco tau tipo 1 |"); P("|---|---|---|---|")
for lo,hi in [(0.5,1),(1,2),(2,5),(5,10),(10,60)]:
    mm=t0&(tau["fsrPmax"]>=lo)&(tau["fsrPmax"]<hi)&(tau["nfsrsh"]>0)
    if mm.sum()==0: continue
    P(f"| {lo}-{hi} | {mm.sum()} | {(tau['nfsrsh_reco_tau'][mm]>0).mean():.3f} | {(tau['rtype'][mm]==1).mean():.3f} |")
P("\nSpectro de P del FSR de shower (gen, en cono): percentiles 10/50/90 = "+"/".join(f"{x:.2f}" for x in np.percentile(tau["fsrPmax"][t0&(tau["nfsrsh"]>0)],[10,50,90]))+" GeV")
# also type 1 with TrueMode 10 (tau->pi nu gamma classified as pi pi0)
mt1=(ph["gtype"]==1)&(ph["gmode"]==10)
P(f"\nGen tau->pi nu (TrueMode 10) etiquetados GenTauType==1 por un foton FSR duro de la desintegracion: {len(set(zip(ph['ievt'][mt1],ph['itau'][mt1])))} taus con reco 1/2 en las filas (compara con {len(set(zip(ph['ievt'][(ph['gtype']==1)],ph['itau'][(ph['gtype']==1)])))} gen type1 reco1/2).")
# FSR shower photon: P distribution figure gen vs reco eff
fig,ax=plt.subplots(figsize=(6,4.2))
x=tau["fsrPmax"][t0&(tau["nfsrsh"]>0)]; xr=tau["fsrPmax"][t0&(tau["nfsrsh_reco_tau"]>0)]
b=np.logspace(-0.3,1.7,31)
ax.hist(x,bins=b,histtype="step",lw=1.8,color="#0072B2",label="gen shower FSR in cone (P>0.5)")
ax.hist(xr,bins=b,histtype="step",lw=1.8,color="#D55E00",label="... assigned to reco tau")
ax.set_xscale("log"); ax.set_xlabel("gen photon P (GeV)"); ax.set_ylabel("N taus"); ax.legend(fontsize=8); ax.set_title("Tau FSR (shower) in gen tau->pi nu")
fig.tight_layout(); fig.savefig("fig_fsr_gen_vs_reco.png"); plt.close(fig)
OUT.close()

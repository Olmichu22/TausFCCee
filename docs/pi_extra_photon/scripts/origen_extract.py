"""Extrae, por gen tau, el foton asignado al reco tau y su origen.
Salida: npz con tablas planas (una fila por foton reco asignado a taus de interes,
y una fila por gen tau para las cuentas globales)."""
import sys, numpy as np, uproot, awkward as ak
F="/nfs/cms/arqolmo/TausFCCee/Results/TauReco/New2MSample_results0.4_tph0.0_tpi0.0_n0.0_g0.0/Tree_resultsdecayAll_0.4_tph0.0_tpi0.0_n0.0_g0.0.root"
OUT=sys.argv[1]; NMAX=int(sys.argv[2]) if len(sys.argv)>2 else None
STEP=100_000

BR=["GenTauType","GenTauTrueMode","GenVisTauP","GenVisTauTheta","GenVisTauPhi","GenTauTheta","RecoMatchedKey",
    "GenTauHasExtraNeutrals",
    "RecoTauType","RecoTauMass","RecoTauP",
    "RecoTauConstKey","RecoConstP","RecoConstTheta","RecoConstPhi","RecoConstPDG",
    "GenPhotonP","GenPhotonTheta","GenPhotonPhi","GenPhotonTauKey","GenPhotonOrigin","GenPhotonParentPDG",
    "RecoPhotonP","RecoPhotonTheta","RecoPhotonPhi","RecoPhotonTauKey","RecoPhotonGenMatchIdx"]

def dR(th1,ph1,th2,ph2):
    dphi=ph1-ph2
    dphi=ak.where(dphi>np.pi,dphi-2*np.pi,dphi); dphi=ak.where(dphi<-np.pi,dphi+2*np.pi,dphi)
    return np.sqrt((th1-th2)**2+dphi**2)

rows={k:[] for k in ["gtype","gmode","gvisP","gcos","rtype","rmass","rP","rtauKey_evt","ievt","itau",
                     "phP","phTheta","phPhi","gmi","gorig","gtk","gpar","gphP","isown","piP","piTheta","piPhi","dR","piDRgen"]}
taus={k:[] for k in ["gtype","gmode","gvisP","gcos","rtype","nfsr","nfsr_reco","nfsr_reco_tau","npirad","npirad_reco","npirad_reco_tau","ievt","itau","fsrPmax","piradPmax","hasextra","nfsrdec","nfsrdec_reco_tau","nfsrsh","nfsrsh_reco_tau","nisr","nisr_reco_tau","fsrdecPmax"]}

t=uproot.open(F)["Tau_tree"]
nev=0
for a in t.iterate(BR,step_size=STEP,entry_stop=NMAX):
    n=len(a); ievt=np.arange(nev,nev+n); nev+=n
    # ---- gen taus
    itau=ak.local_index(a.GenTauType)
    rk=a.RecoMatchedKey; has=rk>=0; rks=ak.where(has,rk,0)
    rtype=ak.where(has,a.RecoTauType[rks],-99)
    rmass=ak.where(has,a.RecoTauMass[rks],-1.); rP=ak.where(has,a.RecoTauP[rks],-1.)
    gt=ak.zip({"gtype":a.GenTauType,"gmode":a.GenTauTrueMode,"gvisP":a.GenVisTauP,"gcos":np.cos(a.GenTauTheta),
               "gth":a.GenVisTauTheta,"gph":a.GenVisTauPhi,"rk":rk,"rtype":rtype,"rmass":rmass,"rP":rP,"itau":itau,
               "ievt":ak.broadcast_arrays(ievt,itau)[0],"hasextra":a.GenTauHasExtraNeutrals})
    # ---- reco photons con info gen enlazada
    gmi=a.RecoPhotonGenMatchIdx; hg=gmi>=0; gmis=ak.where(hg,gmi,0)
    ph=ak.zip({"P":a.RecoPhotonP,"th":a.RecoPhotonTheta,"phi":a.RecoPhotonPhi,"tk":a.RecoPhotonTauKey,"gmi":gmi,
               "gorig":ak.where(hg,a.GenPhotonOrigin[gmis],-1),"gtk":ak.where(hg,a.GenPhotonTauKey[gmis],-2),
               "gpar":ak.where(hg,a.GenPhotonParentPDG[gmis],0),"gphP":ak.where(hg,a.GenPhotonP[gmis],-1.)})
    # ---- reco constituyentes: pion lider por reco tau
    rc=ak.zip({"P":a.RecoConstP,"th":a.RecoConstTheta,"phi":a.RecoConstPhi,"pdg":a.RecoConstPDG,"tk":a.RecoTauConstKey})
    # ---- gen fotones del evento, con flag "reconstruido" (algun reco photon enlaza a el) y tauKey reco
    gp=ak.zip({"P":a.GenPhotonP,"th":a.GenPhotonTheta,"phi":a.GenPhotonPhi,"tk":a.GenPhotonTauKey,"orig":a.GenPhotonOrigin,
               "par":a.GenPhotonParentPDG,"j":ak.local_index(a.GenPhotonP)})
    pr=ak.cartesian({"g":gp,"r":ph},nested=True)
    m=pr.r.gmi==pr.g.j
    gp_reco=ak.any(m,axis=2)
    gp_rtk=ak.fill_none(ak.firsts(pr.r.tk[m],axis=2),-1)
    gp=ak.with_field(gp,gp_reco,"reco"); gp=ak.with_field(gp,gp_rtk,"rtk")

    # ==== seleccion de taus de interes: tipo 0 (cualquier reco) y tipo 1 con reco 1/2
    sel=(gt.gtype==0)|((gt.gtype==1)&((gt.rtype==1)|(gt.rtype==2)))
    gts=gt[sel]
    # pion lider del reco tau
    pc=ak.cartesian({"t":gts,"c":rc},nested=True)
    mp=(pc.c.tk==pc.t.rk)&(pc.t.rk>=0)&(abs(pc.c.pdg)==211)
    pion=ak.firsts(pc.c[mp],axis=2)
    # fotones asignados
    pp=ak.cartesian({"t":gts,"p":ph},nested=True)
    mph=(pp.p.tk==pp.t.rk)&(pp.t.rk>=0)
    phs=pp.p[mph]
    # broadcast tau y pion a fotones
    tb=ak.broadcast_arrays(gts,phs.P)[0]; pib=ak.broadcast_arrays(pion,phs.P)[0]
    fl=lambda x: ak.to_numpy(ak.flatten(ak.flatten(x,axis=2)))
    keep=fl((tb.gtype==0)|(tb.gtype==1))  # todos
    for k,v in [("gtype",tb.gtype),("gmode",tb.gmode),("gvisP",tb.gvisP),("gcos",tb.gcos),("rtype",tb.rtype),("rmass",tb.rmass),("rP",tb.rP),
                ("rtauKey_evt",tb.rk),("ievt",tb.ievt),("itau",tb.itau),
                ("phP",phs.P),("phTheta",phs.th),("phPhi",phs.phi),("gmi",phs.gmi),("gorig",phs.gorig),("gtk",phs.gtk),("gpar",phs.gpar),("gphP",phs.gphP)]:
        rows[k].append(fl(v))
    rows["isown"].append(fl(phs.gtk==tb.itau))
    piP=ak.fill_none(pib.P,-1.); pith=ak.fill_none(pib.th,0.); piph=ak.fill_none(pib.phi,0.)
    rows["piP"].append(fl(piP)); rows["piTheta"].append(fl(pith)); rows["piPhi"].append(fl(piph))
    rows["dR"].append(fl(dR(phs.th,phs.phi,pith,piph)))
    rows["piDRgen"].append(fl(dR(phs.th,phs.phi,tb.gth,tb.gph)))

    # ==== por gen tau tipo 0: fotones gen FSR del tau y radiacion del pion cerca del tau
    g0=gt[gt.gtype==0]
    pg=ak.cartesian({"t":g0,"g":gp},nested=True)
    drg=dR(pg.g.th,pg.g.phi,pg.t.gth,pg.t.gph)
    # FSR "de desintegracion" (constituyente del tau) y FSR "de shower" (TauKey -1, ancestro = copia del tau; se asigna por angulo)
    fsrdec=(pg.g.tk==pg.t.itau)&(pg.g.orig==1)
    fsrsh=(pg.g.tk==-1)&(pg.g.orig==1)&(drg<0.4)&(pg.g.P>0.5)
    fsr=fsrdec|fsrsh
    # radiacion del pion: origen 2, padre +-211, TauKey -1, dentro de dR<0.4 del tau visible, P>0.5 (corte genminP del reco)
    pirad=(pg.g.orig==2)&(abs(pg.g.par)==211)&(drg<0.4)&(pg.g.P>0.5)
    isr=(pg.g.orig==2)&(abs(pg.g.par)==11)&(drg<0.4)&(pg.g.P>0.5)
    f1=lambda x: ak.to_numpy(ak.flatten(x))
    taus["gtype"].append(f1(g0.gtype)); taus["gmode"].append(f1(g0.gmode)); taus["gvisP"].append(f1(g0.gvisP)); taus["gcos"].append(f1(g0.gcos))
    taus["rtype"].append(f1(g0.rtype)); taus["ievt"].append(f1(g0.ievt)); taus["itau"].append(f1(g0.itau)); taus["hasextra"].append(f1(g0.hasextra))
    taus["nfsr"].append(f1(ak.sum(fsr,axis=2))); taus["nfsr_reco"].append(f1(ak.sum(fsr&pg.g.reco,axis=2)))
    taus["nfsr_reco_tau"].append(f1(ak.sum(fsr&(pg.g.rtk==pg.t.rk)&(pg.t.rk>=0),axis=2)))
    taus["npirad"].append(f1(ak.sum(pirad,axis=2))); taus["npirad_reco"].append(f1(ak.sum(pirad&pg.g.reco,axis=2)))
    taus["npirad_reco_tau"].append(f1(ak.sum(pirad&(pg.g.rtk==pg.t.rk)&(pg.t.rk>=0),axis=2)))
    taus["nfsrdec"].append(f1(ak.sum(fsrdec,axis=2))); taus["nfsrdec_reco_tau"].append(f1(ak.sum(fsrdec&(pg.g.rtk==pg.t.rk)&(pg.t.rk>=0),axis=2)))
    taus["nfsrsh"].append(f1(ak.sum(fsrsh,axis=2))); taus["nfsrsh_reco_tau"].append(f1(ak.sum(fsrsh&(pg.g.rtk==pg.t.rk)&(pg.t.rk>=0),axis=2)))
    taus["nisr"].append(f1(ak.sum(isr,axis=2))); taus["nisr_reco_tau"].append(f1(ak.sum(isr&(pg.g.rtk==pg.t.rk)&(pg.t.rk>=0),axis=2)))
    taus["fsrdecPmax"].append(f1(ak.fill_none(ak.max(ak.where(fsrdec,pg.g.P,-1.),axis=2),-1.)))
    taus["fsrPmax"].append(f1(ak.fill_none(ak.max(ak.where(fsr,pg.g.P,-1.),axis=2),-1.)))
    taus["piradPmax"].append(f1(ak.fill_none(ak.max(ak.where(pirad,pg.g.P,-1.),axis=2),-1.)))
    print(nev,flush=True)

np.savez(OUT,**{("ph_"+k):np.concatenate(v) for k,v in rows.items()},**{("tau_"+k):np.concatenate(v) for k,v in taus.items()})
print("done",nev)

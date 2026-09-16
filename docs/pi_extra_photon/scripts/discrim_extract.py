"""Extrae tablas planas (por tau reco emparejado y por foton constituyente) del Tau_tree.
Salida: tables.npz con dos bloques: tau_* (una fila por gen tau) y ph_* (una fila por foton
constituyente de un reco tau emparejado, con indice a la fila del tau)."""
import sys, numpy as np, uproot, awkward as ak
F="/nfs/cms/arqolmo/TausFCCee/Results/TauReco/New2MSample_results0.4_tph0.0_tpi0.0_n0.0_g0.0/Tree_resultsdecayAll_0.4_tph0.0_tpi0.0_n0.0_g0.0.root"
NMAX=(int(sys.argv[1]) or None) if len(sys.argv)>1 else None
OUT=sys.argv[2] if len(sys.argv)>2 else "tables.npz"
BR=["GenTauType","GenVisTauP","GenTauTrueMode","RecoMatchedKey",
    "RecoTauType","RecoTauDM","RecoTauMass","RecoTauP","RecoTauNConsts",
    "RecoConstP","RecoConstTheta","RecoConstPhi","RecoConstPDG",
    "RecoPhotonP","RecoPhotonTauKey","RecoPhotonGenMatchIdx",
    "GenPhotonOrigin","GenPhotonTauKey","GenPhotonP"]
MPI=0.13957
T={k:[] for k in ["gtype","gvisP","gtrue","rtype","rdm","rmass","rP","nph","ev","gkey",
                  "pipx","pipy","pipz","piP","npi","chpx","chpy","chpz","chE"]}
A={k:[] for k in ["gtype","gvisP","matched","gtrue"]}
P={k:[] for k in ["tau","P","theta","phi","dR","ang","px","py","pz","origin","gtaukey","genP","matched"]}
ev0=0; ntau=0
for a in uproot.open(F)["Tau_tree"].iterate(BR,step_size=100000,entry_stop=NMAX):
    n=len(a); ev=np.arange(ev0,ev0+n)
    rk=a.RecoMatchedKey; has=rk>=0
    rks=ak.where(has,rk,0)
    # constituyentes por (evento, reco tau)
    cnt=ak.flatten(a.RecoTauNConsts)
    ntaus=ak.num(a.RecoTauType)
    def per_tau(x): return ak.unflatten(ak.unflatten(ak.flatten(x),cnt),ntaus)
    cP,cTh,cPh,cPDG=[per_tau(a[k]) for k in ["RecoConstP","RecoConstTheta","RecoConstPhi","RecoConstPDG"]]
    # seleccion: gen taus con reco emparejado
    gtype=a.GenTauType[has]; gvis=a.GenVisTauP[has]; gtrue=a.GenTauTrueMode[has]
    gkey=ak.local_index(a.GenTauType)[has]
    rk2=rk[has]
    rtype=a.RecoTauType[rks][has]; rdm=a.RecoTauDM[rks][has]; rmass=a.RecoTauMass[rks][has]; rP=a.RecoTauP[rks][has]
    tP=cP[rks][has]; tTh=cTh[rks][has]; tPh=cPh[rks][has]; tPDG=cPDG[rks][has]
    evt=ak.broadcast_arrays(ev,gtype)[0]
    # pion lider = primer constituyente (const[0]=lead en buildTauFromPion); si no hay consts -> nan
    nc=ak.num(tP,axis=2)
    lead_ok=nc>0
    def first(x,fill): return ak.fill_none(ak.firsts(x,axis=2),fill)
    lP=first(tP,np.nan); lTh=first(tTh,np.nan); lPh=first(tPh,np.nan); lPDG=first(tPDG,0)
    npi=ak.sum(np.abs(tPDG)==211,axis=2)
    ischg=(np.abs(tPDG)==211)
    cpx=ak.sum((tP*np.sin(tTh)*np.cos(tPh))[ischg],axis=2); cpy=ak.sum((tP*np.sin(tTh)*np.sin(tPh))[ischg],axis=2)
    cpz=ak.sum((tP*np.cos(tTh))[ischg],axis=2); cE=ak.sum(np.sqrt(tP**2+MPI**2)[ischg],axis=2)
    A["gtype"].append(ak.to_numpy(ak.flatten(a.GenTauType))); A["gvisP"].append(ak.to_numpy(ak.flatten(a.GenVisTauP)))
    A["matched"].append(ak.to_numpy(ak.flatten(has))); A["gtrue"].append(ak.to_numpy(ak.flatten(a.GenTauTrueMode)))
    isph=(tPDG==22)
    nph=ak.sum(isph,axis=2)
    fl=lambda x: ak.to_numpy(ak.flatten(x,axis=None))
    lPn=fl(lP); lThn=fl(lTh); lPhn=fl(lPh)
    T["gtype"].append(fl(gtype)); T["gvisP"].append(fl(gvis)); T["gtrue"].append(fl(gtrue))
    T["rtype"].append(fl(rtype)); T["rdm"].append(fl(rdm)); T["rmass"].append(fl(rmass)); T["rP"].append(fl(rP))
    T["nph"].append(fl(nph)); T["ev"].append(fl(evt)); T["gkey"].append(fl(gkey)); T["npi"].append(fl(npi))
    T["chpx"].append(fl(cpx)); T["chpy"].append(fl(cpy)); T["chpz"].append(fl(cpz)); T["chE"].append(fl(cE))
    st=np.sin(lThn); T["pipx"].append(lPn*st*np.cos(lPhn)); T["pipy"].append(lPn*st*np.sin(lPhn)); T["pipz"].append(lPn*np.cos(lThn)); T["piP"].append(lPn)
    # tabla de fotones
    ntau_local=len(lPn)
    # indice de fila del tau, propagado a cada constituyente
    tau_row=ak.unflatten(np.arange(ntau,ntau+ntau_local),ak.num(gtype))
    tau_row_c=ak.broadcast_arrays(tau_row,tP)[0]
    phP=tP[isph]; phTh=tTh[isph]; phPh=tPh[isph]; phRow=tau_row_c[isph]
    evc=ak.broadcast_arrays(evt,tP)[0][isph]; rkc=ak.broadcast_arrays(rk2,tP)[0][isph]
    # lead pion propagado a cada foton
    lPc=ak.broadcast_arrays(lP,tP)[0][isph]; lThc=ak.broadcast_arrays(lTh,tP)[0][isph]; lPhc=ak.broadcast_arrays(lPh,tP)[0][isph]
    pP=fl(phP); pTh=fl(phTh); pPh=fl(phPh)
    lPf=fl(lPc); lThf=fl(lThc); lPhf=fl(lPhc)
    dphi=pPh-lPhf; dphi=np.where(dphi>np.pi,2*np.pi-dphi,dphi); dphi=np.where(dphi<-np.pi,2*np.pi+dphi,dphi)
    dR=np.sqrt((pTh-lThf)**2+dphi**2)
    cosang=np.sin(pTh)*np.sin(lThf)*np.cos(pPh-lPhf)+np.cos(pTh)*np.cos(lThf)
    ang=np.arccos(np.clip(cosang,-1,1))
    P["tau"].append(fl(phRow)); P["P"].append(pP); P["theta"].append(pTh); P["phi"].append(pPh); P["dR"].append(dR); P["ang"].append(ang)
    P["px"].append(pP*np.sin(pTh)*np.cos(pPh)); P["py"].append(pP*np.sin(pTh)*np.sin(pPh)); P["pz"].append(pP*np.cos(pTh))
    # procedencia: emparejar con el bloque RecoPhoton por (evento, reco tau key, P bits)
    evf=fl(evc).astype(np.int64); rkf=fl(rkc).astype(np.int64)
    key_c=(evf<<36)|(rkf<<32)|pP.astype(np.float32).view(np.uint32).astype(np.int64)
    bev=ak.broadcast_arrays(ev,a.RecoPhotonP)[0]
    bP=fl(a.RecoPhotonP); bk=fl(a.RecoPhotonTauKey).astype(np.int64); bev=fl(bev).astype(np.int64); bg=fl(a.RecoPhotonGenMatchIdx)
    selb=bk>=0
    key_b=(bev[selb]<<36)|(bk[selb]<<32)|bP[selb].astype(np.float32).view(np.uint32).astype(np.int64)
    # info gen del foton emparejado
    gidx=bg[selb]; gev=bev[selb]
    gO=a.GenPhotonOrigin; gT=a.GenPhotonTauKey; gPP=a.GenPhotonP
    goff=np.concatenate([[0],np.cumsum(ak.to_numpy(ak.num(gO)))])[:-1]
    gOf=fl(gO); gTf=fl(gT); gPf=fl(gPP)
    ok=gidx>=0
    flat_g=goff[gev-ev0]+np.where(ok,gidx,0)
    origin=np.where(ok,gOf[flat_g],-9); gtk=np.where(ok,gTf[flat_g],-9); genP=np.where(ok,gPf[flat_g],np.nan)
    o=np.argsort(key_b); kb=key_b[o]
    pos=np.searchsorted(kb,key_c); pos=np.clip(pos,0,len(kb)-1)
    m=(kb[pos]==key_c) if len(kb) else np.zeros(len(key_c),bool)
    P["matched"].append(m)
    P["origin"].append(np.where(m,origin[o][pos],-99)); P["gtaukey"].append(np.where(m,gtk[o][pos],-99)); P["genP"].append(np.where(m,genP[o][pos],np.nan))
    ntau+=ntau_local; ev0+=n
    print(ev0,ntau,len(np.concatenate(P["P"])),flush=True)
out={("tau_"+k):np.concatenate(v) for k,v in T.items()}
out.update({("ph_"+k):np.concatenate(v) for k,v in P.items()})
out.update({("all_"+k):np.concatenate(v) for k,v in A.items()})
np.savez(OUT,**out)
print({k:v.shape for k,v in out.items()})

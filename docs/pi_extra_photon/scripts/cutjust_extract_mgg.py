"""Extrae, por foton constituyente sin pareja evidente vs con pareja, la masa invariante
m_gg del par de fotones del mismo tau que MINIMIZA |m_gg - m_pi0| (mismo criterio de
'pareja' que discrim_common.py), pero guardando la masa m_gg en si (sin restar m_pi0).
Reutiliza la logica de discrim_extract.py (incluido el matching a GenPhoton para el flag
truepi0), con el mismo parcheo de padding que cutjust_extract.py para que funcione con el
awkward instalado en este entorno. Salida: cutjust_mgg_tables.npz"""
import sys, numpy as np, uproot, awkward as ak
F="/nfs/cms/arqolmo/TausFCCee/Results/TauReco/New2MSample_results0.4_tph0.0_tpi0.0_n0.0_g0.0/Tree_resultsdecayAll_0.4_tph0.0_tpi0.0_n0.0_g0.0.root"
NMAX=(int(sys.argv[1]) or None) if len(sys.argv)>1 else 600000
MPI0=0.13498
BR=["GenTauType","RecoMatchedKey","RecoTauType","RecoTauNConsts",
    "RecoConstP","RecoConstTheta","RecoConstPhi","RecoConstPDG",
    "RecoPhotonP","RecoPhotonTauKey","RecoPhotonGenMatchIdx",
    "GenPhotonOrigin","GenPhotonTauKey","GenPhotonP"]

ti_l=[]; gpx_l=[]; gpy_l=[]; gpz_l=[]; phP_l=[]
tauGtype_l=[]; tauRtype_l=[]; tauGkey_l=[]
origin_l=[]; gtaukey_l=[]
nrows=0; ev0=0
for a in uproot.open(F)["Tau_tree"].iterate(BR,step_size=100000,entry_stop=NMAX):
    n=len(a); ev=np.arange(ev0,ev0+n)
    rk=a.RecoMatchedKey; has=rk>=0; rks=ak.where(has,rk,0)
    RTN=ak.fill_none(ak.pad_none(a.RecoTauNConsts,1,axis=1),0)
    RT=ak.fill_none(ak.pad_none(a.RecoTauType,1,axis=1),-999)
    cnt=ak.flatten(RTN); ntaus=ak.num(RT)
    def per_tau(x): return ak.unflatten(ak.unflatten(ak.flatten(x),cnt),ntaus)
    cP,cTh,cPh,cPDG=[per_tau(a[k]) for k in ["RecoConstP","RecoConstTheta","RecoConstPhi","RecoConstPDG"]]
    gtype=a.GenTauType[has]; gkey=ak.local_index(a.GenTauType)[has]
    rk2=rk[has]
    tP=cP[rks][has]; tTh=cTh[rks][has]; tPh=cPh[rks][has]; tPDG=cPDG[rks][has]
    rtype=RT[rks][has]
    fl=lambda x: ak.to_numpy(ak.flatten(x,axis=None))
    gtn=fl(gtype); rtn=fl(rtype); gkn=fl(gkey)
    ntau_local=len(gtn)
    tau_row=ak.unflatten(np.arange(nrows,nrows+ntau_local),ak.num(gtype))
    isph=(tPDG==22)
    tau_row_c=ak.broadcast_arrays(tau_row,tP)[0]
    evt=ak.broadcast_arrays(ev,gtype)[0]
    evc=ak.broadcast_arrays(evt,tP)[0][isph]; rkc=ak.broadcast_arrays(rk2,tP)[0][isph]
    phP=fl(tP[isph]); phTh=fl(tTh[isph]); phPh=fl(tPh[isph]); phRow=fl(tau_row_c[isph])
    ti_l.append(phRow)
    phP_l.append(phP)
    gpx_l.append(phP*np.sin(phTh)*np.cos(phPh)); gpy_l.append(phP*np.sin(phTh)*np.sin(phPh)); gpz_l.append(phP*np.cos(phTh))
    tauGtype_l.append(gtn); tauRtype_l.append(rtn); tauGkey_l.append(gkn)

    # --- matching a GenPhoton (identico a discrim_extract.py) para origin/gtaukey
    evf=fl(evc).astype(np.int64); rkf=fl(rkc).astype(np.int64)
    key_c=(evf<<36)|(rkf<<32)|phP.astype(np.float32).view(np.uint32).astype(np.int64)
    bev=ak.broadcast_arrays(ev,a.RecoPhotonP)[0]
    bP=fl(a.RecoPhotonP); bk=fl(a.RecoPhotonTauKey).astype(np.int64); bev=fl(bev).astype(np.int64); bg=fl(a.RecoPhotonGenMatchIdx)
    selb=bk>=0
    key_b=(bev[selb]<<36)|(bk[selb]<<32)|bP[selb].astype(np.float32).view(np.uint32).astype(np.int64)
    gidx=bg[selb]; gev=bev[selb]
    gO=a.GenPhotonOrigin; gT=a.GenPhotonTauKey
    goff=np.concatenate([[0],np.cumsum(ak.to_numpy(ak.num(gO)))])[:-1]
    gOf=fl(gO); gTf=fl(gT)
    ok=gidx>=0
    flat_g=goff[gev-ev0]+np.where(ok,gidx,0)
    origin=np.where(ok,gOf[flat_g],-9); gtk=np.where(ok,gTf[flat_g],-9)
    o=np.argsort(key_b); kb=key_b[o]
    pos=np.searchsorted(kb,key_c); pos=np.clip(pos,0,len(kb)-1)
    m=(kb[pos]==key_c) if len(kb) else np.zeros(len(key_c),bool)
    origin_l.append(np.where(m,origin[o][pos],-99)); gtaukey_l.append(np.where(m,gtk[o][pos],-99))

    nrows+=ntau_local; ev0+=n
    print(ev0,nrows,len(np.concatenate(phP_l)),flush=True)

ti=np.concatenate(ti_l); phP=np.concatenate(phP_l)
gpx=np.concatenate(gpx_l); gpy=np.concatenate(gpy_l); gpz=np.concatenate(gpz_l)
tauGtype=np.concatenate(tauGtype_l); tauRtype=np.concatenate(tauRtype_l); tauGkey=np.concatenate(tauGkey_l)
origin=np.concatenate(origin_l); gtaukey=np.concatenate(gtaukey_l)
gtype=tauGtype[ti]; rtype=tauRtype[ti]
truepi0=(origin==0)&(gtaukey==tauGkey[ti])

# pareja pi0: partner que MINIMIZA |m_gg - m_pi0|; guardamos la masa m_gg real de ese par (no la resta)
n=len(ti); best_d=np.full(n,np.inf); best_m=np.full(n,np.nan)
nph_per_tau=np.bincount(ti,minlength=nrows)
kmax=int(nph_per_tau.max()) if n else 1
for k in range(1,kmax):
    i=np.arange(0,n-k); j=i+k
    same=ti[i]==ti[j]
    i=i[same]; j=j[same]
    if len(i)==0: break
    cos=(gpx[i]*gpx[j]+gpy[i]*gpy[j]+gpz[i]*gpz[j])/(phP[i]*phP[j])
    m=np.sqrt(np.maximum(2*phP[i]*phP[j]*(1-np.clip(cos,-1,1)),0))
    d=np.abs(m-MPI0)
    upd_i=d<best_d[i]
    best_d[i[upd_i]]=d[upd_i]; best_m[i[upd_i]]=m[upd_i]
    upd_j=d<best_d[j]
    best_d[j[upd_j]]=d[upd_j]; best_m[j[upd_j]]=m[upd_j]

out=dict(gtype=gtype,rtype=rtype,origin=origin,gtaukey=gtaukey,truepi0=truepi0,mgg_best=best_m,dmpi0_best=best_d)
np.savez("cutjust_mgg_tables.npz",**out)
print({k:v.shape for k,v in out.items()})

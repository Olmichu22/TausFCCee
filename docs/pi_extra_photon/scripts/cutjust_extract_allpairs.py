"""Complemento de cutjust_extract_mgg.py: para cada tau (gen1/gen2 sobre todo), enumera
TODAS las combinaciones gamma-gamma del cono (no solo la que minimiza |m_gg-m_pi0|), y por
separado guarda, por tau, la masa de la pareja seleccionada (el minimo). Sirve para comparar
'todas las combinaciones' frente a 'la pareja elegida' y ver si el minimo realmente saca un
pico limpio de un fondo combinatorio. Salida: cutjust_allpairs_tables.npz"""
import sys, numpy as np, uproot, awkward as ak
F="/nfs/cms/arqolmo/TausFCCee/Results/TauReco/New2MSample_results0.4_tph0.0_tpi0.0_n0.0_g0.0/Tree_resultsdecayAll_0.4_tph0.0_tpi0.0_n0.0_g0.0.root"
NMAX=(int(sys.argv[1]) or None) if len(sys.argv)>1 else 600000
MPI0=0.13498
BR=["GenTauType","RecoMatchedKey","RecoTauType","RecoTauNConsts",
    "RecoConstP","RecoConstTheta","RecoConstPhi","RecoConstPDG"]

ti_l=[]; gpx_l=[]; gpy_l=[]; gpz_l=[]; phP_l=[]
tauGtype_l=[]; tauRtype_l=[]
nrows=0
for a in uproot.open(F)["Tau_tree"].iterate(BR,step_size=100000,entry_stop=NMAX):
    rk=a.RecoMatchedKey; has=rk>=0; rks=ak.where(has,rk,0)
    RTN=ak.fill_none(ak.pad_none(a.RecoTauNConsts,1,axis=1),0)
    RT=ak.fill_none(ak.pad_none(a.RecoTauType,1,axis=1),-999)
    cnt=ak.flatten(RTN); ntaus=ak.num(RT)
    def per_tau(x): return ak.unflatten(ak.unflatten(ak.flatten(x),cnt),ntaus)
    cP,cTh,cPh,cPDG=[per_tau(a[k]) for k in ["RecoConstP","RecoConstTheta","RecoConstPhi","RecoConstPDG"]]
    gtype=a.GenTauType[has]
    tP=cP[rks][has]; tTh=cTh[rks][has]; tPh=cPh[rks][has]; tPDG=cPDG[rks][has]
    rtype=RT[rks][has]
    fl=lambda x: ak.to_numpy(ak.flatten(x,axis=None))
    gtn=fl(gtype); rtn=fl(rtype)
    ntau_local=len(gtn)
    tau_row=ak.unflatten(np.arange(nrows,nrows+ntau_local),ak.num(gtype))
    isph=(tPDG==22)
    tau_row_c=ak.broadcast_arrays(tau_row,tP)[0]
    phP=fl(tP[isph]); phTh=fl(tTh[isph]); phPh=fl(tPh[isph]); phRow=fl(tau_row_c[isph])
    ti_l.append(phRow); phP_l.append(phP)
    gpx_l.append(phP*np.sin(phTh)*np.cos(phPh)); gpy_l.append(phP*np.sin(phTh)*np.sin(phPh)); gpz_l.append(phP*np.cos(phTh))
    tauGtype_l.append(gtn); tauRtype_l.append(rtn)
    nrows+=ntau_local
    print(nrows,len(np.concatenate(phP_l)),flush=True)

ti=np.concatenate(ti_l); phP=np.concatenate(phP_l)
gpx=np.concatenate(gpx_l); gpy=np.concatenate(gpy_l); gpz=np.concatenate(gpz_l)
tauGtype=np.concatenate(tauGtype_l); tauRtype=np.concatenate(tauRtype_l)

n=len(ti)
nph_per_tau=np.bincount(ti,minlength=nrows)
kmax=int(nph_per_tau.max()) if n else 1
pair_ti=[]; pair_m=[]
best_d=np.full(nrows,np.inf); best_m=np.full(nrows,np.nan)
for k in range(1,kmax):
    i=np.arange(0,n-k); j=i+k
    same=ti[i]==ti[j]
    i=i[same]; j=j[same]
    if len(i)==0: break
    cos=(gpx[i]*gpx[j]+gpy[i]*gpy[j]+gpz[i]*gpz[j])/(phP[i]*phP[j])
    m=np.sqrt(np.maximum(2*phP[i]*phP[j]*(1-np.clip(cos,-1,1)),0))
    d=np.abs(m-MPI0)
    tid=ti[i]  # ti[i]==ti[j], indice de tau para este par
    pair_ti.append(tid); pair_m.append(m)
    upd=d<best_d[tid]
    best_d[tid[upd]]=d[upd]; best_m[tid[upd]]=m[upd]

pair_ti=np.concatenate(pair_ti); pair_m=np.concatenate(pair_m)
pair_gtype=tauGtype[pair_ti]; pair_rtype=tauRtype[pair_ti]

out=dict(pair_gtype=pair_gtype,pair_rtype=pair_rtype,pair_mgg=pair_m,
          tau_gtype=tauGtype,tau_rtype=tauRtype,tau_best_mgg=best_m)
np.savez("cutjust_allpairs_tables.npz",**out)
print({k:v.shape for k,v in out.items()})

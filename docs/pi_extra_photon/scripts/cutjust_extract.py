"""Extrae variables por foton constituyente para justificar el corte final recomendado
(Seccion 4 del REPORT): descartar el foton sin pareja pi0 si P_gamma/P_pion < 0.05, o si
P_gamma > 2 GeV y m(pi+gamma) > 1.2 GeV. Subconjunto de discrim_extract.py + discrim_common.py.
Salida: cutjust_tables.npz (una fila por foton constituyente de un reco tau emparejado)."""
import sys, numpy as np, uproot, awkward as ak
F="/nfs/cms/arqolmo/TausFCCee/Results/TauReco/New2MSample_results0.4_tph0.0_tpi0.0_n0.0_g0.0/Tree_resultsdecayAll_0.4_tph0.0_tpi0.0_n0.0_g0.0.root"
NMAX=(int(sys.argv[1]) or None) if len(sys.argv)>1 else 600000
MPI=0.13957; MPI0=0.13498
BR=["GenTauType","RecoMatchedKey","RecoTauType","RecoTauNConsts",
    "RecoConstP","RecoConstTheta","RecoConstPhi","RecoConstPDG"]
ti_l=[]; phP_l=[]; Ppi_l=[]; mpig_l=[]; gpx_l=[]; gpy_l=[]; gpz_l=[]
tauGtype_l=[]; tauRtype_l=[]
nrows=0
for a in uproot.open(F)["Tau_tree"].iterate(BR,step_size=100000,entry_stop=NMAX):
    rk=a.RecoMatchedKey; has=rk>=0; rks=ak.where(has,rk,0)
    # pad a >=1 reco tau (con 0 constituyentes) por evento para que rks=0 siempre sea un indice valido
    RTN=ak.fill_none(ak.pad_none(a.RecoTauNConsts,1,axis=1),0)
    RT=ak.fill_none(ak.pad_none(a.RecoTauType,1,axis=1),-999)
    cnt=ak.flatten(RTN); ntaus=ak.num(RT)
    def per_tau(x): return ak.unflatten(ak.unflatten(ak.flatten(x),cnt),ntaus)
    cP,cTh,cPh,cPDG=[per_tau(a[k]) for k in ["RecoConstP","RecoConstTheta","RecoConstPhi","RecoConstPDG"]]
    gtype=a.GenTauType[has]
    tP=cP[rks][has]; tTh=cTh[rks][has]; tPh=cPh[rks][has]; tPDG=cPDG[rks][has]
    rtype=RT[rks][has]
    def first(x,fill): return ak.fill_none(ak.firsts(x,axis=2),fill)
    lP=first(tP,np.nan); lTh=first(tTh,np.nan); lPh=first(tPh,np.nan)
    fl=lambda x: ak.to_numpy(ak.flatten(x,axis=None))
    lPn=fl(lP); lThn=fl(lTh); lPhn=fl(lPh); gtn=fl(gtype); rtn=fl(rtype)
    ntau_local=len(lPn)
    tau_row=ak.unflatten(np.arange(nrows,nrows+ntau_local),ak.num(gtype))
    isph=(tPDG==22)
    tau_row_c=ak.broadcast_arrays(tau_row,tP)[0]
    lPc=ak.broadcast_arrays(lP,tP)[0][isph]; lThc=ak.broadcast_arrays(lTh,tP)[0][isph]; lPhc=ak.broadcast_arrays(lPh,tP)[0][isph]
    phP=fl(tP[isph]); phTh=fl(tTh[isph]); phPh=fl(tPh[isph]); phRow=fl(tau_row_c[isph])
    lPf=fl(lPc); lThf=fl(lThc); lPhf=fl(lPhc)
    Epi=np.sqrt(lPf**2+MPI**2)
    E=Epi+phP
    px=lPf*np.sin(lThf)*np.cos(lPhf)+phP*np.sin(phTh)*np.cos(phPh)
    py=lPf*np.sin(lThf)*np.sin(lPhf)+phP*np.sin(phTh)*np.sin(phPh)
    pz=lPf*np.cos(lThf)+phP*np.cos(phTh)
    mpig=np.sqrt(np.maximum(E**2-px**2-py**2-pz**2,0))
    ti_l.append(phRow); phP_l.append(phP); Ppi_l.append(lPf); mpig_l.append(mpig)
    gpx_l.append(phP*np.sin(phTh)*np.cos(phPh)); gpy_l.append(phP*np.sin(phTh)*np.sin(phPh)); gpz_l.append(phP*np.cos(phTh))
    tauGtype_l.append(gtn); tauRtype_l.append(rtn)
    nrows+=ntau_local
    print(nrows,len(np.concatenate(phP_l)),flush=True)

ti=np.concatenate(ti_l); phP=np.concatenate(phP_l); Ppi=np.concatenate(Ppi_l); mpig=np.concatenate(mpig_l)
gpx=np.concatenate(gpx_l); gpy=np.concatenate(gpy_l); gpz=np.concatenate(gpz_l)
tauGtype=np.concatenate(tauGtype_l); tauRtype=np.concatenate(tauRtype_l)
gtype=tauGtype[ti]; rtype=tauRtype[ti]

# pareja pi0 (vectorizado, identico a discrim_common.load: ti esta ordenado por construccion)
n=len(ti); best=np.full(n,np.inf)
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
    np.minimum.at(best,i,d); np.minimum.at(best,j,d)

out=dict(gtype=gtype,rtype=rtype,phP=phP,Ppi=Ppi,mpig=mpig,dmpi0=best)
np.savez("cutjust_tables.npz",**out)
print({k:v.shape for k,v in out.items()})

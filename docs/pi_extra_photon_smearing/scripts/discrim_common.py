"""Carga de tables.npz, variables derivadas por foton y emulacion del re-conteo de fotones."""
import numpy as np
MPI=0.13957; MPI0=0.13498
def load(path="tables.npz"):
    d=np.load(path)
    t={k[4:]:d[k] for k in d if k.startswith("tau_")}
    p={k[3:]:d[k] for k in d if k.startswith("ph_")}
    A={k[4:]:d[k] for k in d if k.startswith("all_")}
    # --- variables por foton ---
    ti=p["tau"]
    p["Ppi"]=t["piP"][ti]; p["Ptau"]=t["rP"][ti]
    p["fpi"]=p["P"]/p["Ppi"]; p["ftau"]=p["P"]/p["Ptau"]
    # masa pi+gamma
    Epi=np.sqrt(t["piP"]**2+MPI**2)[ti]
    E=Epi+p["P"]; px=t["pipx"][ti]+p["px"]; py=t["pipy"][ti]+p["py"]; pz=t["pipz"][ti]+p["pz"]
    p["mpig"]=np.sqrt(np.maximum(E**2-px**2-py**2-pz**2,0))
    # numero de fotones en el tau y masa gamma-gamma con los demas fotones del mismo tau
    p["nph"]=t["nph"][ti]
    n=len(ti); best=np.full(n,np.inf); mgg_max=np.full(n,np.nan)   # min |m_gg-m_pi0|; m_gg con el foton mas energetico
    kmax=int(t["nph"].max())
    partnerP=np.zeros(n)
    for k in range(1,kmax):
        i=np.arange(0,n-k); j=i+k
        same=ti[i]==ti[j]
        i=i[same]; j=j[same]
        if len(i)==0: break
        cos=(p["px"][i]*p["px"][j]+p["py"][i]*p["py"][j]+p["pz"][i]*p["pz"][j])/(p["P"][i]*p["P"][j])
        m=np.sqrt(np.maximum(2*p["P"][i]*p["P"][j]*(1-np.clip(cos,-1,1)),0))
        d=np.abs(m-MPI0)
        np.minimum.at(best,i,d); np.minimum.at(best,j,d)
        # m_gg con el companero mas energetico
        for a,b in ((i,j),(j,i)):
            upd=p["P"][b]>partnerP[a]
            partnerP[a[upd]]=p["P"][b[upd]]; mgg_max[a[upd]]=m[upd]
    p["dmpi0"]=best            # inf si el foton esta solo en el cono
    p["mgg_lead"]=mgg_max
    p["haspartner"]=np.isfinite(best)
    # etiquetas de verdad
    p["gtype"]=t["gtype"][ti]; p["rtype"]=t["rtype"][ti]
    p["truepi0"]=(p["origin"]==0)&(p["gtaukey"]==t["gkey"][ti])
    return t,p,A

def emulate(t,p,keep):
    """Recuenta fotones supervivientes (keep: mascara por foton) y devuelve (rtype, rdm, mass, P) nuevos."""
    nph=np.bincount(p["tau"],weights=keep.astype(float),minlength=len(t["gtype"])).astype(int)
    rt=t["rtype"].copy()
    one=(rt>=0)&(rt<10); three=rt>=10
    rt[one]=np.minimum(nph[one],9); rt[three]=10+nph[three]
    dm=dm_of(rt)
    E=t["chE"].copy(); px=t["chpx"].copy(); py=t["chpy"].copy(); pz=t["chpz"].copy()
    w=keep.astype(float)
    E+=np.bincount(p["tau"],weights=p["P"]*w,minlength=len(E)); px+=np.bincount(p["tau"],weights=p["px"]*w,minlength=len(E))
    py+=np.bincount(p["tau"],weights=p["py"]*w,minlength=len(E)); pz+=np.bincount(p["tau"],weights=p["pz"]*w,minlength=len(E))
    P=np.sqrt(px**2+py**2+pz**2); M=np.sqrt(np.maximum(E**2-P**2,0))
    chg=one|three
    M=np.where(chg,M,t["rmass"]); P=np.where(chg,P,t["rP"])
    return rt,dm,M,P,nph

def dm_of(rt):
    dm=rt.copy(); one=(rt>=0)&(rt<10); three=rt>=10
    dm[one]=np.ceil(rt[one]/2).astype(int); dm[three]=10+np.ceil((rt[three]-10)/2).astype(int)
    return dm

"""Carga de tables_*.npz y variables derivadas."""
import numpy as np, glob
MPI=0.13957; MPI0=0.13498
def ang(th1,ph1,th2,ph2):
    c=np.sin(th1)*np.sin(th2)*np.cos(ph1-ph2)+np.cos(th1)*np.cos(th2)
    return np.arccos(np.clip(c,-1,1))
def load(pattern="tables_*.npz"):
    files=sorted(glob.glob(pattern)); T={};R={};G={};C={}; off=0
    for f in files:
        d=np.load(f); n=len(d["tau_ev"])
        for k in d:
            pre,key=k[:4],k[4:]
            D={"tau_":T,"rph_":R,"gph_":G,"cst_":C}[pre]
            v=d[k]
            if pre!="tau_" and key in ("tau","gtau"): v=np.where(v>=0,v+off,-1)
            D.setdefault(key,[]).append(v)
        off+=n
    T={k:np.concatenate(v) for k,v in T.items()}; R={k:np.concatenate(v) for k,v in R.items()}
    G={k:np.concatenate(v) for k,v in G.items()}; C={k:np.concatenate(v) for k,v in C.items()}
    n=len(T["ev"])
    # pion lider reco por tau (primer constituyente cargado)
    T["piP"]=np.full(n,np.nan); T["pith"]=np.full(n,np.nan); T["piphi"]=np.full(n,np.nan)
    ch=np.abs(C["pdg"])==211
    oo=np.where(ch)[0]; oo=oo[np.argsort(-C["pos"][oo],kind="stable")]  # pos 0 el ultimo -> gana
    T["piP"][C["tau"][oo]]=C["P"][oo]; T["pith"][C["tau"][oo]]=C["th"][oo]; T["piphi"][C["tau"][oo]]=C["phi"][oo]
    T["npi"]=np.bincount(C["tau"][ch],minlength=n); T["nneu"]=np.bincount(C["tau"][np.abs(C["pdg"])==2112],minlength=n)
    T["nph"]=np.bincount(R["tau"],minlength=n)
    # fotones reco: variables por foton
    ti=R["tau"]
    R["dRpi"]=ang(R["th"],R["phi"],T["pith"][ti],T["piphi"][ti]); R["fpi"]=R["P"]/T["piP"][ti]
    R["sametau"]=(R["gtau"]==ti)
    # categoria de origen del foton reco: 0 pi0 mismo tau, 1 FSR linea tau, 2 rad cargada/ISR, 3 otro, 4 sim, 5 sin match(fragmento), 6 pi0 de otro tau
    cat=np.full(len(ti),5); has=R["grow"]>=0
    cat[has&(R["gorigin"]==0)&R["sametau"]]=0
    cat[has&(R["gorigin"]==0)&~R["sametau"]]=6
    for o in (1,2,3,4): cat[has&(R["gorigin"]==o)]=o
    R["cat"]=cat
    # fotones gen pi0: destino
    gi=G["tau"]
    G["dRpi"]=ang(G["th"],G["phi"],T["gpith"][gi],T["gpiphi"][gi])
    fate=np.full(len(gi),0)   # 0 no reco como foton; 1 reco en el tau emparejado; 2 reco en otro tau; 3 reco fuera de todo tau
    rec=G["rtk"]>=-1
    fate[rec&(G["rtk"]==T["rk"][gi])&(T["rk"][gi]>=0)]=1
    fate[rec&(G["rtk"]>=0)&(G["rtk"]!=T["rk"][gi])]=2
    fate[rec&(G["rtk"]==-1)]=3
    G["fate"]=fate
    # emparejar los dos fotones del pi0: angulo gamma-gamma y P del hermano (por tau: primer/segundo)
    return T,R,G,C
def cls(rt):
    """clase de resultado reco para tablas: '-99' sin match, '0','1','2','3+','-20','lep','3p','-1'"""
    out=np.full(len(rt),"3p",dtype=object)
    out[rt==-99]="nomatch"; out[rt==-20]="pi->n"; out[np.isin(rt,[-11,-13])]="lep"; out[rt==-1]="noID"
    out[rt==0]="0g"; out[rt==1]="1g"; out[rt==2]="2g"; out[(rt>=3)&(rt<10)]="3g+"
    return out
CLS=["2g","1g","3g+","0g","pi->n","3p","lep","noID","nomatch"]
def table(sel,T,bins,var="gvisP",classes=CLS,label=None):
    x=T[var]; lines=[]
    hdr="| %s | N | "%var+" | ".join(classes)+" |"; lines.append(hdr); lines.append("|"+"---|"*(len(classes)+2))
    c=cls(T["rtype"])
    for lo,hi in zip(bins[:-1],bins[1:]):
        s=sel&(x>=lo)&(x<hi); N=s.sum()
        if N==0: continue
        lines.append("| %g-%g | %d | "%(lo,hi,N)+" | ".join("%.3f"%((c[s]==k).mean()) for k in classes)+" |")
    s=sel; N=s.sum(); lines.append("| todo | %d | "%N+" | ".join("%.3f"%((c[s]==k).mean()) for k in classes)+" |")
    return "\n".join(lines)

def sort_by_tau(R):
    o=np.argsort(R["tau"],kind="stable")
    for k in R: R[k]=R[k][o]
    return R
def cart(P,th,phi):
    return P*np.sin(th)*np.cos(phi),P*np.sin(th)*np.sin(phi),P*np.cos(th)
def photon_pairing(T,R):
    """Por foton (R ordenado por tau): min|m_gg-m_pi0| con otros fotones del mismo tau, m(pi+g), P/Ppi.
    Ademas por tau: masa de todos los fotones, masa pi+fotones, n fotones."""
    n=len(R["tau"]); ti=R["tau"]
    px,py,pz=cart(R["P"],R["th"],R["phi"]); R["px"],R["py"],R["pz"]=px,py,pz
    ppx,ppy,ppz=cart(T["piP"],T["pith"],T["piphi"]); Epi=np.sqrt(T["piP"]**2+MPI**2)
    E=Epi[ti]+R["P"]; sx=ppx[ti]+px; sy=ppy[ti]+py; sz=ppz[ti]+pz
    R["mpig"]=np.sqrt(np.maximum(E**2-sx**2-sy**2-sz**2,0))
    best=np.full(n,np.inf); partner=np.full(n,-1)
    kmax=int(T["nph"].max())
    for k in range(1,kmax):
        i=np.arange(0,n-k); j=i+k; same=ti[i]==ti[j]; i=i[same]; j=j[same]
        if len(i)==0: break
        cos=(px[i]*px[j]+py[i]*py[j]+pz[i]*pz[j])/(R["P"][i]*R["P"][j])
        m=np.sqrt(np.maximum(2*R["P"][i]*R["P"][j]*(1-np.clip(cos,-1,1)),0)); d=np.abs(m-MPI0)
        for a,b in ((i,j),(j,i)):
            upd=d<best[a]; best[a[upd]]=d[upd]; partner[a[upd]]=b[upd]
    R["dmpi0"]=best; R["partner"]=partner
    nt=len(T["ev"])
    gE=np.bincount(ti,weights=R["P"],minlength=nt); gx=np.bincount(ti,weights=px,minlength=nt); gy=np.bincount(ti,weights=py,minlength=nt); gz=np.bincount(ti,weights=pz,minlength=nt)
    T["mgg_all"]=np.sqrt(np.maximum(gE**2-gx**2-gy**2-gz**2,0)); T["Pgg_all"]=np.sqrt(gx**2+gy**2+gz**2)
    E=Epi+gE; T["mtau_all"]=np.sqrt(np.maximum(E**2-(ppx+gx)**2-(ppy+gy)**2-(ppz+gz)**2,0))
    return T,R
def emulate(T,R,keep):
    """recuenta fotones supervivientes por tau -> nuevo rtype (solo cambia 1-prong 0-9)."""
    nt=len(T["ev"]); nph=np.bincount(R["tau"],weights=keep.astype(float),minlength=nt).astype(int)
    rt=T["rtype"].copy(); one=(rt>=0)&(rt<10); rt[one]=np.minimum(nph[one],9)
    return rt
def eff_table(T,sel_rows,rts,names,bins,var="gvisP",target=2):
    lines=["| %s | N | "%var+" | ".join(names)+" |","|"+"---|"*(len(names)+2)]
    for lo,hi in zip(bins[:-1],bins[1:]):
        s=sel_rows&(T[var]>=lo)&(T[var]<hi); N=s.sum()
        lines.append("| %g-%g | %d | "%(lo,hi,N)+" | ".join("%.3f"%((rt[s]==target).mean()) for rt in rts)+" |")
    s=sel_rows; lines.append("| todo | %d | "%s.sum()+" | ".join("%.3f"%((rt[s]==target).mean()) for rt in rts)+" |")
    return "\n".join(lines)

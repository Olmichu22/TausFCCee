"""Tablas planas del Tau_tree para estudiar la migracion del gen rho (tipo 1).

Salida: tables_<i>.npz con cuatro bloques ligados por indices a la fila del gen tau:
  tau_*   : una fila por gen tau (todos, emparejados o no)
  rph_*   : una fila por foton reco (PandoraPFO) asignado al reco tau emparejado a un gen tau
  gph_*   : una fila por foton gen de pi0 que es constituyente de un gen tau (tipos 1-9)
  cst_*   : una fila por constituyente cargado/neutron del reco tau emparejado
Uso: python extract.py <entry_start> <entry_stop> <out.npz>
"""
import sys, numpy as np, uproot, awkward as ak
F="/nfs/cms/arqolmo/TausFCCee/Results/TauReco/New2MSample_results0.4_tph0.0_tpi0.0_n0.0_g0.0/Tree_resultsdecayAll_0.4_tph0.0_tpi0.0_n0.0_g0.0.root"
E0,E1,OUT=int(sys.argv[1]),int(sys.argv[2]),sys.argv[3]
K=1000  # multiplicador de claves compuestas (ev*K+key)
BR=["GenTauType","GenTauTrueMode","GenVisTauP","GenTauP","GenVisTauTheta","GenVisTauPhi","GenTauHasExtraNeutrals",
    "RecoMatchedKey","RecoTauType","RecoTauDM","RecoTauMass","RecoTauP","RecoTauTheta","RecoTauPhi","RecoTauNConsts",
    "RecoTauConstKey","RecoConstP","RecoConstTheta","RecoConstPhi","RecoConstPDG",
    "GenTauConstKey","GenConstP","GenConstTheta","GenConstPhi","GenConstPDG","GenConstPi0Key","GenConstMCIdx",
    "GenPhotonP","GenPhotonTheta","GenPhotonPhi","GenPhotonTauKey","GenPhotonOrigin","GenPhotonParentPDG","GenPhotonAncestorMCIdx","GenPhotonMCIdx",
    "RecoPhotonP","RecoPhotonTheta","RecoPhotonPhi","RecoPhotonTauKey","RecoPhotonGenMatchIdx"]
def flat(a,x): return ak.to_numpy(ak.flatten(a[x]))
def evof(a,x,ev): return ak.to_numpy(ak.flatten(ak.broadcast_arrays(ev,a[x])[0]))
def posof(a,x): return ak.to_numpy(ak.flatten(ak.local_index(a[x])))
def lookup(keys_sorted, rows_sorted, q):
    """indice de fila para cada clave q (o -1)."""
    i=np.searchsorted(keys_sorted,q); i=np.clip(i,0,len(keys_sorted)-1)
    ok=keys_sorted[i]==q if len(keys_sorted) else np.zeros(len(q),bool)
    return np.where(ok,rows_sorted[np.minimum(i,len(rows_sorted)-1)] if len(rows_sorted) else -1,-1)
T={};R={};G={};C={}
ev0=E0; toff=0
for a in uproot.open(F)["Tau_tree"].iterate(BR,step_size=50000,entry_start=E0,entry_stop=E1):
    n=len(a); ev=ak.Array(np.arange(ev0,ev0+n)); ev0+=n
    # ---------- gen taus ----------
    t={}
    t["ev"]=evof(a,"GenTauType",ev); t["key"]=posof(a,"GenTauType")
    for k,b in [("gtype","GenTauType"),("gtrue","GenTauTrueMode"),("gvisP","GenVisTauP"),("gP","GenTauP"),
                ("gth","GenVisTauTheta"),("gphi","GenVisTauPhi"),("gextra","GenTauHasExtraNeutrals"),("rk","RecoMatchedKey")]:
        t[k]=flat(a,b)
    ntau=len(t["ev"])
    # reco tau emparejado
    rkey=np.where(t["rk"]>=0,t["ev"]*K+t["rk"],-1)
    rt={}
    rt["ev"]=evof(a,"RecoTauType",ev); rt["key"]=posof(a,"RecoTauType")
    for k,b in [("rtype","RecoTauType"),("rdm","RecoTauDM"),("rmass","RecoTauMass"),("rP","RecoTauP"),("rth","RecoTauTheta"),("rphi","RecoTauPhi"),("rnc","RecoTauNConsts")]:
        rt[k]=flat(a,b)
    rkeys=rt["ev"]*K+rt["key"]; o=np.argsort(rkeys); rkeys_s=rkeys[o]; rrow_s=np.arange(len(rkeys))[o]
    rrow=lookup(rkeys_s,rrow_s,rkey); rrow=np.where(t["rk"]>=0,rrow,-1)
    for k in ["rtype","rdm","rmass","rP","rth","rphi","rnc"]:
        t[k]=np.where(rrow>=0,rt[k][np.maximum(rrow,0)],-99)
    # cuantos gen taus apuntan al mismo reco tau (ambiguedad)
    u,inv,cnt=np.unique(rkey,return_inverse=True,return_counts=True)
    t["rshared"]=np.where(rkey>=0,cnt[inv]-1,0)
    # ---------- mapa (ev, reco tau key) -> fila gen tau (primera) ----------
    sel=t["rk"]>=0
    gk=rkey[sel]; gr=np.arange(ntau)[sel]
    o=np.argsort(gk,kind="stable"); gk_s=gk[o]; gr_s=gr[o]
    first=np.concatenate([[True],gk_s[1:]!=gk_s[:-1]]); gk_s=gk_s[first]; gr_s=gr_s[first]
    # ---------- mapa (ev, gen tau key) -> fila gen tau ----------
    tkeys=t["ev"]*K+t["key"]; o=np.argsort(tkeys); tk_s=tkeys[o]; tr_s=np.arange(ntau)[o]
    # ---------- constituyentes reco ----------
    c={}
    c["ev"]=evof(a,"RecoTauConstKey",ev); c["tk"]=flat(a,"RecoTauConstKey")
    for k,b in [("P","RecoConstP"),("th","RecoConstTheta"),("phi","RecoConstPhi"),("pdg","RecoConstPDG")]: c[k]=flat(a,b)
    c["pos"]=posof(a,"RecoTauConstKey")
    c["tau"]=lookup(gk_s,gr_s,c["ev"]*K+c["tk"])
    c["tau"]=np.where(c["tau"]>=0,c["tau"]+toff,-1)
    keep=(c["tau"]>=0)&(np.abs(c["pdg"])!=22)
    for k in c: C.setdefault(k,[]).append(c[k][keep])
    # ---------- fotones gen ----------
    g={}
    g["ev"]=evof(a,"GenPhotonP",ev); g["pos"]=posof(a,"GenPhotonP")
    for k,b in [("P","GenPhotonP"),("th","GenPhotonTheta"),("phi","GenPhotonPhi"),("tk","GenPhotonTauKey"),("origin","GenPhotonOrigin"),("parent","GenPhotonParentPDG"),("anc","GenPhotonAncestorMCIdx"),("mc","GenPhotonMCIdx")]:
        g[k]=flat(a,b)
    g["tau"]=np.where(g["tk"]>=0,lookup(tk_s,tr_s,g["ev"]*K+np.maximum(g["tk"],0)),-1)
    # ---------- fotones reco ----------
    r={}
    r["ev"]=evof(a,"RecoPhotonP",ev); r["pos"]=posof(a,"RecoPhotonP")
    for k,b in [("P","RecoPhotonP"),("th","RecoPhotonTheta"),("phi","RecoPhotonPhi"),("tk","RecoPhotonTauKey"),("gm","RecoPhotonGenMatchIdx")]:
        r[k]=flat(a,b)
    # fila del foton gen emparejado
    gkeys=g["ev"]*K+g["pos"]; o=np.argsort(gkeys); gk2_s=gkeys[o]; gr2_s=np.arange(len(gkeys))[o]
    r["grow"]=np.where(r["gm"]>=0,lookup(gk2_s,gr2_s,r["ev"]*K+np.maximum(r["gm"],0)),-1)
    gr_=np.maximum(r["grow"],0); has=r["grow"]>=0
    r["gorigin"]=np.where(has,g["origin"][gr_],-2); r["gtau"]=np.where(has,g["tau"][gr_],-1)
    r["gP"]=np.where(has,g["P"][gr_],np.nan); r["ganc"]=np.where(has,g["anc"][gr_],-1); r["gparent"]=np.where(has,g["parent"][gr_],0)
    # fila del gen tau cuyo reco tau contiene el foton
    r["tau"]=np.where(r["tk"]>=0,lookup(gk_s,gr_s,r["ev"]*K+np.maximum(r["tk"],0)),-1)
    # cuantos fotones reco apuntan al mismo foton gen (cluster partido)
    key=np.where(has,r["ev"]*K+r["gm"],-1); u,inv,cnt=np.unique(key,return_inverse=True,return_counts=True)
    r["nshare"]=np.where(has,cnt[inv],0)
    # para cada foton gen: primer foton reco emparejado (fila) y su reco tau key
    gm_first=np.full(len(gkeys),-1); 
    idx=np.where(has)[0]
    # ultimo gana; usamos el de mayor P: ordenar por P asc y sobrescribir
    oo=idx[np.argsort(r["P"][idx])]
    gm_first[r["grow"][oo]]=oo
    g["rrow"]=gm_first
    g["rP"]=np.where(gm_first>=0,r["P"][np.maximum(gm_first,0)],np.nan)
    g["rtk"]=np.where(gm_first>=0,r["tk"][np.maximum(gm_first,0)],-2)  # -2 sin reco, -1 reco fuera de tau
    g["rnshare"]=np.where(gm_first>=0,r["nshare"][np.maximum(gm_first,0)],0)
    # gen pi0 photons de constituyentes de taus
    g["tau"]=np.where(g["tau"]>=0,g["tau"]+toff,-1); r["gtau"]=np.where(r["gtau"]>=0,r["gtau"]+toff,-1); r["tau"]=np.where(r["tau"]>=0,r["tau"]+toff,-1)
    keepg=(g["tau"]>=0)&(g["origin"]==0)
    for k in g: G.setdefault(k,[]).append(g[k][keepg])
    keepr=r["tau"]>=0
    for k in r: R.setdefault(k,[]).append(r[k][keepr])
    # ---------- gen pion (cargado) de cada tau ----------
    gc={}
    gc["ev"]=evof(a,"GenTauConstKey",ev); gc["tk"]=flat(a,"GenTauConstKey")
    for k,b in [("P","GenConstP"),("th","GenConstTheta"),("phi","GenConstPhi"),("pdg","GenConstPDG")]: gc[k]=flat(a,b)
    ch=np.isin(np.abs(gc["pdg"]),[211,321])
    row=lookup(tk_s,tr_s,gc["ev"]*K+gc["tk"])
    t["gpiP"]=np.full(ntau,np.nan); t["gpith"]=np.full(ntau,np.nan); t["gpiphi"]=np.full(ntau,np.nan); t["gnch"]=np.zeros(ntau,int)
    np.add.at(t["gnch"],row[ch],1)
    # el de mayor P
    oo=np.where(ch)[0]; oo=oo[np.argsort(gc["P"][oo])]
    t["gpiP"][row[oo]]=gc["P"][oo]; t["gpith"][row[oo]]=gc["th"][oo]; t["gpiphi"][row[oo]]=gc["phi"][oo]
    for k in t: T.setdefault(k,[]).append(t[k])
    toff+=ntau
    print(ev0,ntau,flush=True)
out={}
for pre,D in [("tau_",T),("rph_",R),("gph_",G),("cst_",C)]:
    for k in D: out[pre+k]=np.concatenate(D[k])
np.savez_compressed(OUT,**out)
print("saved",OUT,{k:v.shape for k,v in out.items() if k.endswith("ev")})

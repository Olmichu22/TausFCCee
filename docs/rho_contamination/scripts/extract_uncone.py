"""Fotones reco NO asignados a ningun tau: dR a cada reco tau del evento (para estudiar conos mas anchos)."""
import sys, numpy as np, uproot, awkward as ak
F="/nfs/cms/arqolmo/TausFCCee/Results/TauReco/New2MSample_results0.4_tph0.0_tpi0.0_n0.0_g0.0/Tree_resultsdecayAll_0.4_tph0.0_tpi0.0_n0.0_g0.0.root"
E0,E1,OUT=int(sys.argv[1]),int(sys.argv[2]),sys.argv[3]
BR=["GenTauType","RecoMatchedKey","RecoTauType","RecoTauTheta","RecoTauPhi","RecoTauP",
    "RecoPhotonP","RecoPhotonTheta","RecoPhotonPhi","RecoPhotonTauKey","RecoPhotonGenMatchIdx","GenPhotonOrigin","GenPhotonTauKey"]
O={k:[] for k in ["ev","gkey","gtype","rtype","rP","dR","P","gorigin","gtaukey"]}
ev0=E0
for a in uproot.open(F)["Tau_tree"].iterate(BR,step_size=50000,entry_start=E0,entry_stop=E1):
    n=len(a); ev=ak.Array(np.arange(ev0,ev0+n)); ev0+=n
    free=a.RecoPhotonTauKey<0
    ph=ak.zip({"P":a.RecoPhotonP[free],"th":a.RecoPhotonTheta[free],"phi":a.RecoPhotonPhi[free],"gm":a.RecoPhotonGenMatchIdx[free]})
    has=a.RecoMatchedKey>=0; rk=ak.where(has,a.RecoMatchedKey,0)
    tau=ak.zip({"gkey":ak.local_index(a.GenTauType),"gtype":a.GenTauType,"rtype":a.RecoTauType[rk],"th":a.RecoTauTheta[rk],"phi":a.RecoTauPhi[rk],"rP":a.RecoTauP[rk]})[has]
    pr=ak.cartesian({"t":tau,"p":ph},axis=1)
    t=pr.t; p=pr.p
    cos=np.sin(t.th)*np.sin(p.th)*np.cos(t.phi-p.phi)+np.cos(t.th)*np.cos(p.th)
    dR=np.arccos(np.minimum(np.maximum(cos,-1.0),1.0))
    keep=dR<1.0
    evb=ak.to_numpy(ak.flatten((dR*0+ev)[keep])).astype(np.int64)
    gm=ak.to_numpy(ak.flatten(p.gm[keep]))
    # tabla plana de fotones gen del chunk
    gev=ak.to_numpy(ak.flatten(ak.broadcast_arrays(ev,a.GenPhotonOrigin)[0])).astype(np.int64); gpos=ak.to_numpy(ak.flatten(ak.local_index(a.GenPhotonOrigin)))
    gor_all=ak.to_numpy(ak.flatten(a.GenPhotonOrigin)); gtk_all=ak.to_numpy(ak.flatten(a.GenPhotonTauKey))
    gk=gev*1000+gpos; o=np.argsort(gk); gk=gk[o]; gor_all=gor_all[o]; gtk_all=gtk_all[o]
    q=evb*1000+np.maximum(gm,0); i=np.clip(np.searchsorted(gk,q),0,max(len(gk)-1,0)); ok=(gm>=0)&(gk[i]==q)
    gor=np.where(ok,gor_all[i],-2); gtk=np.where(ok,gtk_all[i],-9)
    for k,v in [("gkey",t.gkey[keep]),("gtype",t.gtype[keep]),("rtype",t.rtype[keep]),("rP",t.rP[keep]),("dR",dR[keep]),("P",p.P[keep])]:
        O[k].append(ak.to_numpy(ak.flatten(v)))
    O["ev"].append(evb); O["gorigin"].append(gor); O["gtaukey"].append(gtk)
    print(ev0,flush=True)
np.savez_compressed(OUT,**{k:np.concatenate(v) for k,v in O.items()})

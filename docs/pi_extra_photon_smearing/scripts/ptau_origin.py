# Migracion 0->1 en bins de GenTauP (momento TOTAL del tau) por origen del foton, y P_gamma vs E_beam-P_tau
import uproot, awkward as ak, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

# Guard (smearing): evento sin reco taus al final de un chunk -> indice 0 fuera del buffer.
# Rellenar a >=1 elemento no cambia nada: esos valores se enmascaran con has/hg.
pad1=lambda x: ak.fill_none(ak.pad_none(x,1,axis=1),0,axis=1)
padl=lambda x: ak.fill_none(ak.pad_none(x,1,axis=1),[],axis=1)
F="/nfs/cms/arqolmo/TausFCCee/Results/TauReco/New2MSample_smearing_results0.4_tph0.0_tpi0.0_n0.0_g0.0/Tree_resultsdecayAll_0.4_tph0.0_tpi0.0_n0.0_g0.0.root"
BR=["GenTauType","GenTauTrueMode","GenTauP","GenVisTauP","RecoMatchedKey","RecoTauType","beamE",
    "RecoPhotonP","RecoPhotonTauKey","RecoPhotonGenMatchIdx","GenPhotonOrigin","GenPhotonTauKey","GenPhotonParentPDG","GenPhotonP"]
out={k:[] for k in ["gP","gvisP","rtype","cat","phP","gphP","beamE"]}
for a in uproot.iterate(F+":Tau_tree",BR,step_size=100000):
    sel=(a.GenTauType==0)&(a.GenTauTrueMode==10)
    rk=a.RecoMatchedKey; has=rk>=0; rks=ak.where(has,rk,0)
    rtype=ak.where(has,pad1(a.RecoTauType)[rks],-99)
    gt=ak.zip({"gP":a.GenTauP,"gvisP":a.GenVisTauP,"rk":rk,"rtype":rtype})[sel]
    gmi=a.RecoPhotonGenMatchIdx; hg=gmi>=0; gmis=ak.where(hg,gmi,0)
    orig=ak.where(hg,pad1(a.GenPhotonOrigin)[gmis],-1); gtk=ak.where(hg,pad1(a.GenPhotonTauKey)[gmis],-2); gpar=ak.where(hg,pad1(a.GenPhotonParentPDG)[gmis],0)
    cat=ak.where(~hg,0, ak.where(orig==1,1, ak.where((orig==2)&(abs(gpar)==11),2, ak.where(orig==0,3,4))))  # 0 sin match,1 FSR,2 ISR,3 pi0,4 otro
    ph=ak.zip({"P":a.RecoPhotonP,"tk":a.RecoPhotonTauKey,"cat":cat,"gP":ak.where(hg,pad1(a.GenPhotonP)[gmis],-1.)})
    pr=ak.cartesian({"t":gt,"p":ph},nested=True)
    m=(pr.p.tk==pr.t.rk)&(pr.t.rk>=0)&(pr.t.rtype==1)
    first=ak.firsts(pr.p[m],axis=2)
    be=ak.broadcast_arrays(a.beamE,gt.gP)[0]
    out["gP"].append(ak.to_numpy(ak.flatten(gt.gP))); out["gvisP"].append(ak.to_numpy(ak.flatten(gt.gvisP)))
    out["rtype"].append(ak.to_numpy(ak.flatten(gt.rtype))); out["beamE"].append(ak.to_numpy(ak.flatten(be)))
    out["cat"].append(ak.to_numpy(ak.flatten(ak.fill_none(first.cat,-1)))); out["phP"].append(ak.to_numpy(ak.flatten(ak.fill_none(first.P,-1.))))
    out["gphP"].append(ak.to_numpy(ak.flatten(ak.fill_none(first.gP,-1.))))
d={k:np.concatenate(v) for k,v in out.items()}; np.savez("ptau_origin.npz",**d)
gP,rt,cat,phP,gphP,be=d["gP"],d["rtype"],d["cat"],d["phP"],d["gphP"],d["beamE"]
print("beamE unique:",np.unique(np.round(be,2))[:5])
bins=np.array([0,10,15,20,25,30,35,40,42,44,45,46,47]); c=0.5*(bins[1:]+bins[:-1])
names={0:"no gen match",1:"tau-line FSR",2:"ISR",3:"pi0 other tau",4:"other"}
print("GenTauP  N   reco0  reco1  reco-20 unm  | reco1: sinmatch FSR ISR  | E_beam-P_tau med | P_gamma med FSR")
rows={k:[] for k in ["tot"]+list(names)}
for lo,hi in zip(bins[:-1],bins[1:]):
    m=(gP>=lo)&(gP<hi); n=m.sum(); den=((rt[m]==0)|(rt[m]==1)).sum()
    r1=m&(rt==1)
    fr={k:((cat[r1]==k).sum()/max(den,1)) for k in names}
    rows["tot"].append(r1.sum()/max(den,1)); [rows[k].append(fr[k]) for k in names]
    fs=r1&(cat==1)
    print(f"{lo:2d}-{hi:2d} {n:7d} {np.mean(rt[m]==0):.3f} {np.mean(rt[m]==1):.3f} {np.mean(rt[m]==-20):.3f} {np.mean(rt[m]==-99):.3f} | {fr[0]:.3f} {fr[1]:.3f} {fr[2]:.3f} | {np.median(be[m]-gP[m]):6.2f} | {np.median(phP[fs]) if fs.sum() else 0:6.2f}")
fig,ax=plt.subplots(1,2,figsize=(13,4.6))
ax[0].plot(c,rows["tot"],"k-o",label="total reco1/(reco0+reco1)")
for k,col in zip([0,1,2],["#009E73","#0072B2","#E69F00"]): ax[0].plot(c,rows[k],"-s",color=col,label=names[k])
ax[0].set_xlabel("total gen tau P (GeV)"); ax[0].set_ylabel("fraction of gen τ→πν (TrueMode 10)"); ax[0].set_title("π → π+γ migration vs the tau's TOTAL P"); ax[0].legend(); ax[0].grid(alpha=.3)
fs=(rt==1)&(cat==1)
ax[1].hist2d(be[fs]-gP[fs],gphP[fs],bins=[np.linspace(0,45,60),np.linspace(0,45,60)],cmap="Blues")
ax[1].plot([0,45],[0,45],"r--",lw=1,label="P_γ = E_beam − P_τ")
ax[1].set_xlabel("E_beam − P_τ gen (GeV)"); ax[1].set_ylabel("gen P of the assigned FSR photon (GeV)"); ax[1].set_title("Tau-line FSR: the photon carries the energy missing from the tau"); ax[1].legend()
plt.tight_layout(); plt.savefig("fig_migracion_vs_GenTauP.png",dpi=130)

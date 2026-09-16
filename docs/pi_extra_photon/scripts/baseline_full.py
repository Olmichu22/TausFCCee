import uproot, awkward as ak, numpy as np, json
f="/nfs/cms/arqolmo/TausFCCee/Results/TauReco/New2MSample_results0.4_tph0.0_tpi0.0_n0.0_g0.0/Tree_resultsdecayAll_0.4_tph0.0_tpi0.0_n0.0_g0.0.root"
bins=np.arange(0,55,5); cats={"reco0":0,"reco1":1,"reco2":2,"reco-20":-20,"unmatched":-99}
H={k:{c:np.zeros(len(bins)-1) for c in list(cats)+["other","N"]} for k in ["GenVisTauP","GenTauP"]}
for a in uproot.iterate(f+":Tau_tree",["GenTauType","GenVisTauP","GenTauP","RecoMatchedKey","RecoTauType","GenTauTrueMode"],step_size=200000):
    sel=(a.GenTauType==0)
    rk=a.RecoMatchedKey
    rt=ak.where(rk>=0, a.RecoTauType[ak.where(rk>=0,rk,0)], -99)
    rtp=ak.to_numpy(ak.flatten(rt[sel]))
    for k in H:
        x=ak.to_numpy(ak.flatten(a[k][sel]))
        H[k]["N"]+=np.histogram(x,bins)[0]
        oth=np.ones_like(rtp,bool)
        for c,v in cats.items():
            H[k][c]+=np.histogram(x[rtp==v],bins)[0]; oth&=(rtp!=v)
        H[k]["other"]+=np.histogram(x[oth],bins)[0]
for k in H:
    print(k); print("bin      N   "+"  ".join(f"{c:>9s}" for c in list(cats)+["other"]))
    for i in range(len(bins)-1):
        n=H[k]["N"][i]
        print(f"{bins[i]:2d}-{bins[i+1]:2d} {int(n):7d} "+"  ".join(f"{H[k][c][i]/max(n,1):9.4f}" for c in list(cats)+["other"]))
json.dump({k:{c:v.tolist() for c,v in d.items()} for k,d in H.items()},open("baseline_full.json","w"))

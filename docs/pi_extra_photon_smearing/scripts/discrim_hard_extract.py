"""Complemento a tables.npz: GenTauP/Theta/Phi alineados con las filas tau_* (gen taus con reco emparejado)
y con las filas all_* (todos los gen taus). Mismo orden de iteracion que extract.py."""
import numpy as np, uproot, awkward as ak
F="/nfs/cms/arqolmo/TausFCCee/Results/TauReco/New2MSample_smearing_results0.4_tph0.0_tpi0.0_n0.0_g0.0/Tree_resultsdecayAll_0.4_tph0.0_tpi0.0_n0.0_g0.0.root"
T={k:[] for k in ["gP","gtheta","gphi"]}; A={"gP":[]}
for a in uproot.open(F)["Tau_tree"].iterate(["GenTauP","GenTauTheta","GenTauPhi","RecoMatchedKey"],step_size=200000):
    has=a.RecoMatchedKey>=0
    fl=lambda x: ak.to_numpy(ak.flatten(x))
    T["gP"].append(fl(a.GenTauP[has])); T["gtheta"].append(fl(a.GenTauTheta[has])); T["gphi"].append(fl(a.GenTauPhi[has]))
    A["gP"].append(fl(a.GenTauP))
np.savez("hard_tables.npz",**{"tau_"+k:np.concatenate(v) for k,v in T.items()},**{"all_"+k:np.concatenate(v) for k,v in A.items()})
print({k:len(np.concatenate(v)) for k,v in T.items()},len(np.concatenate(A["gP"])))

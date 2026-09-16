"""Base vs corregido en datos reales: migracion gen->reco por tipo gen (seccion 5).
Recuperado del codigo inline original; solo cambia BASE/nombre de fichero a ztt_2M_smearing.
Ejecutar desde /nfs/cms/arqolmo/TausFCCee (usa modules/)."""
import sys, os, math
sys.path.insert(0, "/nfs/cms/arqolmo/TausFCCee")
import ROOT
from podio import root_io
from modules import tauReco, myutils
from multiprocessing import Pool

BASE = "/pnfs/ciemat.es/data/calice/local/arqolmo/CLDFCCFullSim/Ztt_SM_2M_SigmaZ"
CFG = {"enable": True, "modes": ["pion_photon_fsr"]}
dRMax = 0.4

def run(fname):
    from collections import Counter
    cb, cc = Counter(), Counter()   # (genType, recoID)
    reader = root_io.Reader([fname])
    for event in reader.get("events"):
        mc = event.get("MCParticles"); pfos = event.get("PandoraPFOs")
        gen = tauReco.findAllGenTaus(mc)
        base = tauReco.findAllTaus(pfos, dRMax, 0., 0., 0., 0.)
        corr = tauReco.findAllTaus(pfos, dRMax, 0., 0., 0., 0., extra_correction=CFG)
        for g in gen.values():
            gid = g.getID()
            if gid < 0 or gid > 12:
                continue
            gv = g.getvisMomentum()
            for coll, cnt in ((base, cb), (corr, cc)):
                best, bdr = -1, 1.0
                for j in range(len(coll)):
                    dr = myutils.dRAngle(coll[j].getMomentum(), gv)
                    if dr < bdr:
                        bdr, best = dr, j
                cnt[(gid, coll[best].getID() if best >= 0 else -99)] += 1
    return cb, cc

if __name__ == "__main__":
    files = [f"{BASE}/out_reco_edm4hep_{i}.root" for i in range(1000, 1008)]
    from collections import Counter
    CB, CC = Counter(), Counter()
    with Pool(8) as p:
        for cb, cc in p.map(run, files):
            CB.update(cb); CC.update(cc)
    for gid in sorted(set(k[0] for k in CB)):
        tot = sum(v for k, v in CB.items() if k[0] == gid)
        if tot < 50: continue
        def frac(C, rid): return sum(v for k, v in C.items() if k[0]==gid and k[1]==rid)/tot
        print(f"gen {gid} (N={tot}):")
        for rid in sorted(set(k[1] for k in list(CB)+list(CC) if k[0]==gid)):
            b, c = frac(CB, rid), frac(CC, rid)
            if b < 0.005 and c < 0.005: continue
            print(f"   reco {rid:>4}: base {b:.4f} -> corr {c:.4f}  ({c-b:+.4f})")

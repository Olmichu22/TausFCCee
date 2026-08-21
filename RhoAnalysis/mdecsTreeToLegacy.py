#!/usr/bin/env python3
"""mdecsTreeToLegacy.py — MDecs tree → flat legacy tree for the cut optimizer.

`RhoAnalysis/optimize_cuts.py` (and `modules/optimize_pso/`) expect the legacy
one-hemisphere-per-branch schema (`recoMesonP`, `lepP`, `genTauID`, …), while the
MDecs pipeline writes a two-hemisphere tree (`tau1_*` / `tau2_*`). This adapter
projects an MDecs tree onto the legacy schema for one hadron+lepton channel, so
the PSO optimizer can be reused unchanged.

Selection (reco ids, rho = 2): one hemisphere with `recoTauID == --target`, the
other with `recoTauID` in `--other` (default: -11 -13). `--other` must not
contain `--target`, so that exactly one hemisphere fills each slot. For a
leptonic target (e.g. `--target -11 --other 0 2 10 -13`) the "meson" slot holds
the lepton, and `--meson-cut` is then a cut on the lepton momentum.
The target hemisphere becomes the "meson" slot (dec0) and the other one the
"lepton" slot (dec1) — the same slot convention `RhoHistFromTree_MDecs_parallel.py` uses, so
the optimized cuts can be fed back as `--meson-cut` / `--lepton-cut` / `--ang` /
`--zmass-cut`.

`genTauID` is the GEN decay id of the target hemisphere (rho = 1), which is what
the optimizer compares against `--selectGEN` to split signal from migrations.

Usage:
    python RhoAnalysis/mdecsTreeToLegacy.py \
        --tree-file Results/RhoAnalysis/<dir>/TTree_MDecs_*.root \
        --target 2 -o Results/RhoAnalysis/<dir>/legacy_rho_lep.root
"""
import argparse
import glob

import numpy as np
import uproot


def build(tree_file, target, others, tree_name="outtree_original"):
    """Return a dict of flat legacy arrays for the target+lepton channel."""
    branches = [
        "beamE", "ZMass",
        "tau1_recoTauID", "tau2_recoTauID",
        "tau1_decayID", "tau2_decayID",
        "tau1_recoVisP", "tau1_recoVisE", "tau1_recoVisTheta", "tau1_recoVisPhi",
        "tau2_recoVisP", "tau2_recoVisE", "tau2_recoVisTheta", "tau2_recoVisPhi",
    ]
    with uproot.open(f"{tree_file}:{tree_name}") as t:
        a = t.arrays(branches, library="np")

    id1 = a["tau1_recoTauID"].astype(int)
    id2 = a["tau2_recoTauID"].astype(int)
    other_set = np.asarray(sorted(others))

    # slot assignment: tau1 = objetivo, o tau2 = objetivo (nunca ambos, porque
    # --other excluye a --target)
    m1 = (id1 == target) & np.isin(id2, other_set)
    m2 = (id2 == target) & np.isin(id1, other_set)

    out = {}
    keep = m1 | m2
    meson_is_tau1 = m1[keep]

    def _sel(field):
        v1, v2 = a[f"tau1_{field}"][keep], a[f"tau2_{field}"][keep]
        return np.where(meson_is_tau1, v1, v2), np.where(meson_is_tau1, v2, v1)

    for src, dst_meson, dst_lep in [
        ("recoVisP",     "recoMesonP",     "lepP"),
        ("recoVisE",     "recoMesonE",     "lepE"),
        ("recoVisTheta", "recoMesonTheta", "lepTheta"),
        ("recoVisPhi",   "recoMesonPhi",   "lepPhi"),
    ]:
        mes, lep = _sel(src)
        out[dst_meson] = mes.astype(np.float32)
        out[dst_lep] = lep.astype(np.float32)

    gen_mes, _ = _sel("decayID")
    out["genTauID"] = gen_mes.astype(np.int32)
    out["beamE"] = a["beamE"][keep].astype(np.float32)
    out["ZMass"] = a["ZMass"][keep].astype(np.float32)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tree-file", required=True,
                   help="MDecs tree (glob accepted)")
    p.add_argument("--tree-name", default="outtree_original")
    p.add_argument("--target", type=int, required=True,
                   help="Reco decay id of the target hemisphere "
                        "(0=pi, 2=rho, 10=a1, -11=e, -13=mu)")
    p.add_argument("--other", type=int, nargs="+", default=[-11, -13],
                   help="Reco decay ids accepted in the other hemisphere")
    p.add_argument("-o", "--outfile", required=True)
    args = p.parse_args()
    if args.target in args.other:
        raise SystemExit("--other must not contain --target (slot assignment "
                         "would be ambiguous)")

    paths = sorted(glob.glob(args.tree_file))
    if not paths:
        raise SystemExit(f"No input matched {args.tree_file}")

    chunks = [build(pth, args.target, args.other, args.tree_name) for pth in paths]
    data = {k: np.concatenate([c[k] for c in chunks]) for k in chunks[0]}

    with uproot.recreate(args.outfile) as f:
        f[args.tree_name] = data
    print(f"[OK] {args.outfile}: {len(data['recoMesonP'])} entries "
          f"(target={args.target}, other={args.other})")


if __name__ == "__main__":
    main()

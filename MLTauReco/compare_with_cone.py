"""Compare ParTauDETR taus with the cone taus of modules/tauReco.

Both algorithms are run through TauAnalysis/TTreesTausLong.py on the same
events, so the gen taus (tauReco.findAllGenTaus) and the gen-reco matching
(tauReco.MatchRecoGenTau) are this repo's for both. From the Z->tautau trees
the hadronic gen taus that a reco tau was matched to are evaluated with
ml-tau-model's tools (mltau/tools/evaluation):

  - decay mode: per-class F1 comparison and confusion matrices
    (decay_mode.HardLabelDecayModeEvaluator, reduced to ml-tau-model's six
    classes h, h+pi0, h+>=2pi0, 3h, 3h+>=pi0, rare),
  - kinematics: visible pT response / resolution vs gen visible pT
    (kinematics.RegressionEvaluator) and the median dR to the gen visible tau
    (kinematics.DeltaREvaluator).

The Z->qq trees give the rate of reco hadronic taus per event (fakes).

  MLTauReco/run_mltau.sh python3 MLTauReco/compare_with_cone.py \
      --trees Cone:z=<cone Z->tautau tree dir> Cone:qq=<cone Z->qq tree dir> \
              ParTauDETR:z=<ML Z->tautau tree dir> ParTauDETR:qq=<ML Z->qq tree dir> \
      --out-dir <plots dir>
"""

import argparse
import glob
import json
import os

import awkward as ak
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot
from omegaconf import OmegaConf

from mltau.tools import features as f
from mltau.tools.evaluation import decay_mode as d
from mltau.tools.evaluation import kinematics as k

ML_TAU_MODEL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "external", "ml-tau-model")
STYLES = {
    "Cone": {"name": "Cone", "marker": "o", "hatch": "/", "color": "tab:orange", "ls": "dashed",
             "label": "Cone", "lw": 3, "marker_size": 15},
    "ParTauDETR": {"name": "ParTauDETR", "marker": "X", "hatch": ".", "color": "tab:blue", "ls": "solid",
                   "label": "ParTauDETR", "lw": 3, "marker_size": 15},
}
TREE_BRANCHES = [
    "GenEventId", "numRecoTaus", "GenTauType", "GenTauQ", "GenVisTauPt", "GenVisTauTheta", "GenVisTauPhi",
    "GenVisTauEta", "RecoMatchedKey", "RecoTauType", "RecoTauDM", "RecoTauQ", "RecoTauPt", "RecoTauEta",
    "RecoTauTheta", "RecoTauPhi",
]


def load_cfg():
    """ml-tau-model's metrics config, plus plot styles for the two algorithms."""
    cdir = os.path.join(ML_TAU_MODEL, "mltau", "config", "metrics")
    metrics = OmegaConf.load(os.path.join(cdir, "metrics.yaml"))
    for sub in ("kinematics", "decay_mode"):
        metrics = OmegaConf.merge(metrics, OmegaConf.load(os.path.join(cdir, f"{sub}.yaml")))
    metrics = OmegaConf.merge(metrics, {"ALGORITHM_PLOT_STYLES": STYLES})
    return OmegaConf.create({"metrics": metrics})


def reduce_dm(dm):
    """Gen ID / RecoTauDM (n_pi0 for 1 prong, 10 + n_pi0 for 3 prongs) -> the
    six classes of ml-tau-model (ntupelizer get_reduced_decaymodes): 0, 1, 2
    (h + >=2 pi0), 10, 11 (3h + >=1 pi0), 15 (rare / not a hadronic tau)."""
    dm = np.asarray(dm)
    out = np.full(dm.shape, 15)
    out[dm == 0] = 0
    out[dm == 1] = 1
    out[(dm >= 2) & (dm < 10)] = 2
    out[dm == 10] = 10
    out[(dm >= 11) & (dm < 20)] = 11
    return out


def read_tree(path):
    files = sorted(glob.glob(os.path.join(path, "Tree_*.root")))
    if len(files) != 1:
        raise FileNotFoundError(f"expected one Tree_*.root in {path}, found {files}")
    return uproot.open(files[0])["Tau_tree"].arrays(TREE_BRANCHES)


def matched_gen_taus(t):
    """One row per hadronic gen tau: gen properties and, if matched, the reco tau."""
    gen_id = t.GenTauType
    hadronic = (gen_id != -11) & (gen_id != -13)
    key = t.RecoMatchedKey
    idx = ak.mask(key, key >= 0)  # None where unmatched
    cols = {
        "gen_id": gen_id, "gen_q": t.GenTauQ, "gen_pt": t.GenVisTauPt, "gen_theta": t.GenVisTauTheta,
        "gen_phi": t.GenVisTauPhi, "gen_eta": t.GenVisTauEta, "matched": key >= 0,
    }
    for name in ("RecoTauType", "RecoTauDM", "RecoTauQ", "RecoTauPt", "RecoTauEta", "RecoTauTheta", "RecoTauPhi"):
        cols[name] = ak.fill_none(t[name][idx], -999)
    out = {c: ak.to_numpy(ak.flatten(v[hadronic])) for c, v in cols.items()}
    out["event"] = ak.to_numpy(ak.flatten(ak.broadcast_arrays(t.GenEventId, gen_id)[0][hadronic]))
    # The tree's event order depends on the number of workers; fix it, so the
    # gen taus line up between the trees of the two algorithms.
    order = np.lexsort((out["gen_pt"], out["event"]))
    return {c: v[order] for c, v in out.items()}


def reco_hadronic_taus_per_event(t):
    """Reco taus with a hadronic ID (>= 0: excludes e / mu and cone non-taus)."""
    return ak.to_numpy(ak.sum(t.RecoTauType >= 0, axis=1))


def dm_metrics(sel, gen):
    truth, pred = reduce_dm(gen["gen_id"][sel]), reduce_dm(gen["RecoTauDM"][sel])
    pt_ratio = gen["RecoTauPt"][sel] / gen["gen_pt"][sel]
    q25, q50, q75 = np.percentile(pt_ratio, [25, 50, 75])
    hadronic_reco = gen["RecoTauType"][sel] >= 0
    return {
        "n": int(sel.sum()),
        "dm_accuracy": float(np.mean(truth == pred)),
        # only taus matched to a reco tau with a hadronic ID (not e / mu, not a
        # cone the cone algorithm itself rejects)
        "dm_accuracy_hadronic_reco": float(np.mean((truth == pred)[hadronic_reco])),
        "charge_accuracy": float(np.mean(gen["RecoTauQ"][sel] == gen["gen_q"][sel])),
        "pt_response": float(q50),
        "pt_resolution": float((q75 - q25) / q50),
    }


def save(fig, out_dir, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out_dir, f"{name}.{ext}"), bbox_inches="tight", dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--trees", nargs="+", required=True,
                    help="<algorithm>:<sample>=<TTreesTausLong output dir>, algorithm Cone | ParTauDETR, sample z | qq")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    cfg = load_cfg()

    trees = {}
    for spec in args.trees:
        key, path = spec.split("=", 1)
        trees[tuple(key.split(":"))] = read_tree(path)
    algorithms = [a for a in ("Cone", "ParTauDETR") if (a, "z") in trees]

    summary = {}
    dm_evaluators, pt_evaluators, dr_evaluators = [], [], []
    gens = {}
    for algo in algorithms:
        gen = matched_gen_taus(trees[(algo, "z")])
        gens[algo] = gen
        sel = gen["matched"]
        n_evt = {s: len(trees[(algo, s)]) for s in ("z", "qq") if (algo, s) in trees}
        summary[algo] = {
            "events": n_evt,
            "hadronic_gen_taus": int(len(sel)),
            "matched_fraction": float(sel.mean()),
            "matched_to_hadronic_reco_fraction": float((sel & (gen["RecoTauType"] >= 0)).mean()),
            "reco_hadronic_taus_per_event": {
                s: float(reco_hadronic_taus_per_event(trees[(algo, s)]).mean()) for s in n_evt
            },
            **dm_metrics(sel, gen),
        }
        dm_evaluators.append(d.HardLabelDecayModeEvaluator(
            predicted=reduce_dm(gen["RecoTauDM"][sel]), truth=reduce_dm(gen["gen_id"][sel]),
            sample="z", algorithm=algo,
        ))
        pt_evaluators.append(k.RegressionEvaluator(
            prediction=gen["RecoTauPt"][sel], truth=gen["gen_pt"][sel],
            bin_edges=cfg.metrics.kinematics.pt.bin_edges.z, algorithm=algo, sample_name="z",
            mode="ratio", variable="pt",
        ))
        dr = f.deltaR_thetaPhi(theta1=gen["RecoTauTheta"][sel], phi1=gen["RecoTauPhi"][sel],
                               theta2=gen["gen_theta"][sel], phi2=gen["gen_phi"][sel])
        dr_evaluators.append(k.DeltaREvaluator(
            deltaR=np.asarray(dr), pt_truth=gen["gen_pt"][sel],
            bin_edges=cfg.metrics.kinematics.pt.bin_edges.z, algorithm=algo,
        ))
        summary[algo]["dm_class_F1"] = dict(zip(
            ["h", "h+pi0", "h+>=2pi0", "3h", "3h+>=pi0", "rare"],
            map(float, dm_evaluators[-1].class_performances["F1"]),
        ))
        summary[algo]["median_dR"] = float(dr_evaluators[-1].median)

    # Gen taus matched by both algorithms (same events in both trees)
    if len(algorithms) == 2:
        a, b = (gens[x] for x in algorithms)
        if np.array_equal(a["event"], b["event"]) and np.allclose(a["gen_pt"], b["gen_pt"]):
            both = a["matched"] & b["matched"]
            summary["common_matched"] = {x: dm_metrics(both, gens[x]) for x in algorithms}

    # --- decay mode: F1 comparison and confusion matrices
    dm_dir = os.path.join(args.out_dir, "decay_mode")
    os.makedirs(dm_dir, exist_ok=True)
    dmcp = d.DecayModeComparisonPlot(cfg=cfg, metric="F1")
    for i, (ev, off) in enumerate(zip(dm_evaluators, d.get_offsets(len(dm_evaluators)))):
        dmcp.add_line(ev, offset=off, annotation_offset=d.annotation_offsets[i % len(d.annotation_offsets)])
    save(dmcp.fig, dm_dir, "decay_mode_F1")
    for ev in dm_evaluators:
        save(d.ConfusionMatrix(ev).fig, dm_dir, f"decay_mode_cm_{ev.algorithm}")

    # --- kinematics: pT response / resolution, median dR
    kin_dir = os.path.join(args.out_dir, "kinematics")
    os.makedirs(kin_dir, exist_ok=True)
    pt_cfg = cfg.metrics.kinematics.pt
    rme = k.RegressionMultiEvaluator(os.path.join(kin_dir, "pt"), cfg, "z", pt_cfg, axhline_loc=1.0)
    rme.combine_results(pt_evaluators)
    save(rme.response_lineplot.fig, os.path.join(kin_dir, "pt"), "responses")
    rme.resolution_lineplot._autoscale_y()
    save(rme.resolution_lineplot.fig, os.path.join(kin_dir, "pt"), "resolutions")
    for algo in algorithms:
        save(rme.bin_distributions_plots[algo].fig, os.path.join(kin_dir, "pt"), f"{algo}_z_bin_contents")
        save(rme.resolution_2d_plots[algo].fig, os.path.join(kin_dir, "pt"), f"{algo}_z_2D_resolution")
    dr_cfg = cfg.metrics.kinematics.deltaR.median_plot
    lp = k.LinePlot(cfg=cfg, xlabel=dr_cfg.xlabel, ylabel=dr_cfg.ylabel, xscale=dr_cfg.xscale,
                    yscale=dr_cfg.yscale, ymin=dr_cfg.ylim[0], ymax=dr_cfg.ylim[1], nticks=dr_cfg.nticks)
    for ev in dr_evaluators:
        lp.add_line(ev.bin_centers, ev.medians, ev.algorithm, label=ev.algorithm)
    lp._autoscale_y()
    os.makedirs(os.path.join(kin_dir, "deltaR"), exist_ok=True)
    save(lp.fig, os.path.join(kin_dir, "deltaR"), "median_plot")

    with open(os.path.join(args.out_dir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

"""Generator origin of the photons that are not reconstructed in the 45 GeV dip.

Companion to ``plot_polar_dip.py``.  That figure shows *where* the photons of the
43.3-45.0 GeV bin are lost (both beam-pipe spikes); this one shows *what they are*.
Every generated photon of the bin is classified with the generator ancestry stored
in the tau tree (``GenPhotonOrigin`` / ``GenPhotonParentPDG``) and, for the ones
with no reconstructed counterpart, drawn as a donut chart split by origin and by
whether the photon was inside the calorimeter acceptance at all.

A photon counts as reconstructed when some photon PFO of the event points back to
it through ``RecoPhotonGenMatchIdx`` (the ``RecoMCTruthLink`` association).

Input is the ``Tau_tree`` produced by TauAnalysis/TTreesTausLong.py.
"""

import argparse
import os

import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot

BRANCHES = [
    "GenPhotonP", "GenPhotonTheta", "GenPhotonOrigin", "GenPhotonParentPDG",
    "RecoPhotonGenMatchIdx",
]

# Borde de aceptancia: |cos(theta)| = 0.99. Por debajo de este angulo respecto al
# haz la eficiencia medida es exactamente 0, por encima ~0.87 (ver REPORT.md).
THETA_ACCEPTANCE = 0.141

# Origen -> (etiqueta en dos lineas, color).  Rojos = ISR, azules = FSR del tau,
# gris = resto.  La etiqueta se parte en la coma para estrechar la leyenda.
CATEGORIES = [
    ("isr_out",  (r"ISR ($e^{\pm}$ radiation),", "outside acceptance"), "#B2182B"),
    ("isr_in",   (r"ISR ($e^{\pm}$ radiation),", "inside acceptance"),  "#E8887A"),
    ("fsr_in",   (r"$\tau$ FSR,",                "inside acceptance"),  "#2166AC"),
    ("fsr_out",  (r"$\tau$ FSR,",                "outside acceptance"), "#92C5DE"),
    ("other",    (r"$\pi^{0}$ / other ancestry", ""),                   "#9E9E9E"),
]


def load_photons(path, tree, p_min, p_max, step_size="200 MB"):
    """Return a record array of generated photons with ``p_min <= P < p_max``."""
    chunks = []
    for a in uproot.open(path)[tree].iterate(BRANCHES, step_size=step_size):
        gen_idx = ak.local_index(a.GenPhotonP)
        matched = a.RecoPhotonGenMatchIdx[a.RecoPhotonGenMatchIdx >= 0]
        pairs = ak.cartesian([gen_idx, matched], axis=1, nested=True)
        reco = ak.any(pairs["0"] == pairs["1"], axis=-1)

        sel = (a.GenPhotonP >= p_min) & (a.GenPhotonP < p_max)
        chunks.append(ak.zip({
            "P": a.GenPhotonP,
            "theta": a.GenPhotonTheta,
            "origin": a.GenPhotonOrigin,
            "parent": a.GenPhotonParentPDG,
            "reco": reco,
        })[sel])
    return ak.flatten(ak.concatenate(chunks))


def classify(df, theta_acc=THETA_ACCEPTANCE):
    """Split the photons into the wedges of the donut."""
    theta = ak.to_numpy(df.theta)
    origin = ak.to_numpy(df.origin)
    parent = np.abs(ak.to_numpy(df.parent))
    # Angulo al haz mas cercano: los dos picos ISR (+z y -z) se pliegan en uno
    theta_beam = np.minimum(theta, np.pi - theta)
    inside = theta_beam >= theta_acc

    is_isr = (origin == 2) & (parent == 11)
    is_fsr = (origin == 1) & (parent == 15)
    return {
        "isr_out": is_isr & ~inside,
        "isr_in": is_isr & inside,
        "fsr_in": is_fsr & inside,
        "fsr_out": is_fsr & ~inside,
        "other": ~(is_isr | is_fsr),
    }, theta_beam


def make_figure(df, detector, p_min, p_max, outdir, theta_acc=THETA_ACCEPTANCE, dpi=200):
    reco = ak.to_numpy(df.reco)
    masks, theta_beam = classify(df, theta_acc)

    lost = ~reco
    n_lost = int(lost.sum())
    counts = np.array([int((masks[key] & lost).sum()) for key, _, _ in CATEGORIES])
    fracs = counts / n_lost

    keep = counts > 0
    colors = [c for (_, _, c), k in zip(CATEGORIES, keep) if k]
    explode = [0.02 if k.startswith("isr") else 0.0
               for (k, _, _), ok in zip(CATEGORIES, keep) if ok]

    fig, ax = plt.subplots(figsize=(5.3, 2.9))
    fig.subplots_adjust(left=0.005, right=0.435, top=0.78, bottom=0.02)

    wedges, _, _ = ax.pie(
        counts[keep], colors=colors, explode=explode, startangle=90,
        counterclock=False, radius=1.0,
        wedgeprops=dict(edgecolor="white", linewidth=1.2),
        autopct=lambda p: f"{p:.1f}%" if p >= 20.0 else "",
        pctdistance=0.62, textprops=dict(fontsize=9, fontweight="bold", color="white"),
    )
    ax.set_aspect("equal")

    labels = []
    for (_, (head, tail), _), n, f, ok in zip(CATEGORIES, counts, fracs, keep):
        if not ok:
            continue
        tally = f"{n:,} ({100 * f:.1f}%)"
        labels.append(f"{head}\n{tail}  {tally}" if tail else f"{head}\n{tally}")
    ax.legend(wedges, labels, loc="center left", bbox_to_anchor=(1.03, 0.5),
              fontsize=7.5, frameon=False, labelspacing=0.7,
              handlelength=0.9, handleheight=0.9, borderpad=0)

    fig.text(0.02, 0.905, f"{detector}", fontsize=11, fontweight="bold", va="bottom")
    fig.text(0.115, 0.908, r"$Z \to \tau^{+}\tau^{-}$, full simulation",
             fontsize=8.5, va="bottom", color="0.3")
    fig.text(0.02, 0.885,
             "Origin of the photons not reconstructed in the efficiency dip,\n"
             rf"$p_{{\gamma}}^{{\mathrm{{gen}}}} \in [{p_min:.1f},\ {p_max:.1f}]$ GeV",
             fontsize=8.5, va="top", color="0.15", linespacing=1.45)

    os.makedirs(outdir, exist_ok=True)
    stem = os.path.join(outdir, f"photon_dip_lost_origin_{detector}")
    for ext in ("png", "pdf"):
        fig.savefig(f"{stem}.{ext}", dpi=dpi)
        print(f"Saved -> {stem}.{ext}")
    plt.close(fig)

    n_tot = len(reco)
    eff = float(reco.mean())
    out_frac = float((theta_beam[lost] < theta_acc).mean())
    print(f"\nBin {p_min}-{p_max} GeV: N = {n_tot}, eff = {eff:.3f}, lost = {n_lost}, "
          f"outside acceptance = {out_frac:.3f} of the losses")
    for (key, _, _), n, f in zip(CATEGORIES, counts, fracs):
        print(f"  {key:8s} {n:5d}  {100 * f:5.1f}%")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-i", "--input", required=True, help="Tau_tree ROOT file")
    ap.add_argument("-t", "--tree", default="Tau_tree")
    ap.add_argument("-d", "--detector", default="CLD")
    ap.add_argument("-o", "--outdir", default="docs/photon_isr_45GeV")
    ap.add_argument("--pmin", type=float, default=43.333)
    ap.add_argument("--pmax", type=float, default=45.0)
    ap.add_argument("--theta-acceptance", type=float, default=THETA_ACCEPTANCE)
    args = ap.parse_args()

    df = load_photons(args.input, args.tree, args.pmin, args.pmax)
    make_figure(df, args.detector, args.pmin, args.pmax, args.outdir,
                theta_acc=args.theta_acceptance)


if __name__ == "__main__":
    main()

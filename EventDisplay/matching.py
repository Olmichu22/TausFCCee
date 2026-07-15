"""matching.py — Two Gen↔PFO matching strategies.

Both return list[dict] with unified column schema suitable for Dash DataTable.
"""

from __future__ import annotations

import math
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from config import PDG_NAMES

try:
    import ROOT as _ROOT
    from modules.myutils import dRAngle as _dRAngle
    _HAS_ROOT = True
except Exception:
    _HAS_ROOT = False


def _pdg_label(pdg: int) -> str:
    return PDG_NAMES.get(int(pdg), str(pdg)) if not (isinstance(pdg, float) and math.isnan(pdg)) else "—"


def _fmt(v, fmt=".3f"):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "—"
    return format(v, fmt)


# ── Strategy 1: dR matching ───────────────────────────────────────────────────

def match_by_dr(gen_particles, pfos, max_dr: float = 0.1) -> list[dict]:
    """
    Nearest-dR matching between gen particles and PFOs.
    One best PFO per gen particle; unmatched gen get pfo_idx=-999.
    """
    rows = []
    for gen in gen_particles:
        best     = None
        best_dr  = float("inf")

        if _HAS_ROOT:
            gen_p4 = _ROOT.TLorentzVector()
            pt = gen.p * math.sin(gen.theta)
            gen_p4.SetXYZT(
                pt * math.cos(gen.phi),
                pt * math.sin(gen.phi),
                gen.p * math.cos(gen.theta),
                gen.energy,
            )
            for pfo in pfos:
                pfo_p4 = _ROOT.TLorentzVector()
                rpt = pfo.p * math.sin(pfo.theta)
                pfo_p4.SetXYZM(
                    rpt * math.cos(pfo.phi),
                    rpt * math.sin(pfo.phi),
                    pfo.p * math.cos(pfo.theta),
                    0.0,
                )
                dr = float(_dRAngle(gen_p4, pfo_p4))
                if dr < best_dr:
                    best_dr  = dr
                    best     = pfo
        else:
            # Fallback: plain eta-phi dR
            gen_eta = -math.log(math.tan(gen.theta / 2.0)) if 0 < gen.theta < math.pi else 0.0
            for pfo in pfos:
                pfo_eta = -math.log(math.tan(pfo.theta / 2.0)) if 0 < pfo.theta < math.pi else 0.0
                dphi = abs(gen.phi - pfo.phi)
                if dphi > math.pi:
                    dphi = 2 * math.pi - dphi
                dr = math.sqrt((gen_eta - pfo_eta) ** 2 + dphi ** 2)
                if dr < best_dr:
                    best_dr = dr
                    best    = pfo

        if best is not None and best_dr <= max_dr:
            rows.append({
                "gen_idx":    gen.idx,
                "gen_pdg":    gen.pdg,
                "gen_name":   _pdg_label(gen.pdg),
                "gen_E":      round(gen.energy, 4),
                "gen_p":      round(gen.p, 4),
                "pfo_idx":    best.idx,
                "pfo_pdg":    best.pdg,
                "pfo_name":   _pdg_label(best.pdg),
                "pfo_E":      round(best.energy, 4),
                "pfo_p":      round(best.p, 4),
                "match_val":  round(best_dr, 5),
                "match_type": "dR",
                "matched":    True,
            })
        else:
            rows.append({
                "gen_idx":    gen.idx,
                "gen_pdg":    gen.pdg,
                "gen_name":   _pdg_label(gen.pdg),
                "gen_E":      round(gen.energy, 4),
                "gen_p":      round(gen.p, 4),
                "pfo_idx":    -999,
                "pfo_pdg":    -999,
                "pfo_name":   "—",
                "pfo_E":      float("nan"),
                "pfo_p":      float("nan"),
                "match_val":  float("nan"),
                "match_type": "dR",
                "matched":    False,
            })

    # Add unmatched PFOs
    matched_pfo_idxs = {r["pfo_idx"] for r in rows if r["matched"]}
    for pfo in pfos:
        if pfo.idx not in matched_pfo_idxs:
            rows.append({
                "gen_idx":    -999,
                "gen_pdg":    -999,
                "gen_name":   "—",
                "gen_E":      float("nan"),
                "gen_p":      float("nan"),
                "pfo_idx":    pfo.idx,
                "pfo_pdg":    pfo.pdg,
                "pfo_name":   _pdg_label(pfo.pdg),
                "pfo_E":      round(pfo.energy, 4),
                "pfo_p":      round(pfo.p, 4),
                "match_val":  float("nan"),
                "match_type": "dR",
                "matched":    False,
            })

    return rows


# ── Strategy 2: RecoMCTruthLink matching ─────────────────────────────────────

def match_by_truth_link(
    truth_links,
    gen_particles,
    pfos,
    weight_mode: str = "decoded",
    dedup: str = "reco",
) -> list[dict]:
    """
    Matching via RecoMCTruthLink collection weights (Bohdan Dudar encoding).

    weight_mode: 'decoded' → prefer track weight, fallback cluster weight.
                 'raw'     → use encoded weight directly.
    dedup:       'reco'    → one gen per PFO (resolves photon fusion).
                 'gen'     → one PFO per gen (legacy).
    """
    if not truth_links:
        return _build_unmatched_table(gen_particles, pfos, "TruthLink")

    # Build arrays
    gen_idxs    = np.array([l.gen_idx    for l in truth_links], dtype=int)
    pfo_idxs    = np.array([l.pfo_idx    for l in truth_links], dtype=int)
    weights_raw = np.array([l.weight     for l in truth_links], dtype=float)
    gen_pdgs    = np.array([l.gen_pdg    for l in truth_links], dtype=int)
    pfo_pdgs    = np.array([l.pfo_pdg    for l in truth_links], dtype=int)
    gen_Es      = np.array([l.gen_energy for l in truth_links], dtype=float)
    pfo_Es      = np.array([l.pfo_energy for l in truth_links], dtype=float)

    if weight_mode == "decoded":
        encoded   = weights_raw.astype(int)
        track_w   = (encoded % 10000) / 1000.0
        cluster_w = (encoded // 10000) / 1000.0
        eff_w     = np.where(track_w > 0.0, track_w, cluster_w)
    else:
        eff_w = weights_raw

    # Deduplication
    if dedup == "reco":
        order = np.argsort(pfo_idxs)
        gi, pi, wi, gpdg, rpdg, ge, re = (
            gen_idxs[order], pfo_idxs[order], eff_w[order],
            gen_pdgs[order], pfo_pdgs[order], gen_Es[order], pfo_Es[order],
        )
        keep = []
        prev = None
        for i in range(len(pi)):
            if pi[i] != prev:
                keep.append(i)
                prev = pi[i]
            else:
                if wi[i] > wi[keep[-1]]:
                    keep[-1] = i
        idx_keep = keep
    else:  # dedup == "gen"
        order = np.argsort(gen_idxs)
        gi, pi, wi, gpdg, rpdg, ge, re = (
            gen_idxs[order], pfo_idxs[order], eff_w[order],
            gen_pdgs[order], pfo_pdgs[order], gen_Es[order], pfo_Es[order],
        )
        keep = []
        prev = None
        for i in range(len(gi)):
            if gi[i] != prev:
                keep.append(i)
                prev = gi[i]
            else:
                if wi[i] > wi[keep[-1]]:
                    keep[-1] = i
        idx_keep = keep

    matched_gen_idxs = set()
    matched_pfo_idxs = set()
    rows = []
    for i in idx_keep:
        matched_gen_idxs.add(int(gi[i]))
        matched_pfo_idxs.add(int(pi[i]))
        rows.append({
            "gen_idx":    int(gi[i]),
            "gen_pdg":    int(gpdg[i]),
            "gen_name":   _pdg_label(int(gpdg[i])),
            "gen_E":      round(float(ge[i]), 4),
            "gen_p":      float("nan"),
            "pfo_idx":    int(pi[i]),
            "pfo_pdg":    int(rpdg[i]),
            "pfo_name":   _pdg_label(int(rpdg[i])),
            "pfo_E":      round(float(re[i]), 4),
            "pfo_p":      float("nan"),
            "match_val":  round(float(wi[i]), 5),
            "match_type": "TruthLink",
            "matched":    True,
        })

    # Unmatched gen
    gen_map = {g.idx: g for g in gen_particles}
    for gen in gen_particles:
        if gen.idx not in matched_gen_idxs:
            rows.append({
                "gen_idx":    gen.idx, "gen_pdg": gen.pdg,
                "gen_name":   _pdg_label(gen.pdg),
                "gen_E":      round(gen.energy, 4), "gen_p": round(gen.p, 4),
                "pfo_idx":    -999, "pfo_pdg": -999, "pfo_name": "—",
                "pfo_E":      float("nan"), "pfo_p": float("nan"),
                "match_val":  float("nan"), "match_type": "TruthLink", "matched": False,
            })

    # Unmatched PFO
    for pfo in pfos:
        if pfo.idx not in matched_pfo_idxs:
            rows.append({
                "gen_idx":    -999, "gen_pdg": -999, "gen_name": "—",
                "gen_E":      float("nan"), "gen_p": float("nan"),
                "pfo_idx":    pfo.idx, "pfo_pdg": pfo.pdg,
                "pfo_name":   _pdg_label(pfo.pdg),
                "pfo_E":      round(pfo.energy, 4), "pfo_p": round(pfo.p, 4),
                "match_val":  float("nan"), "match_type": "TruthLink", "matched": False,
            })

    return rows


def _build_unmatched_table(gen_particles, pfos, match_type: str) -> list[dict]:
    rows = []
    for gen in gen_particles:
        rows.append({
            "gen_idx": gen.idx, "gen_pdg": gen.pdg,
            "gen_name": _pdg_label(gen.pdg),
            "gen_E": round(gen.energy, 4), "gen_p": round(gen.p, 4),
            "pfo_idx": -999, "pfo_pdg": -999, "pfo_name": "—",
            "pfo_E": float("nan"), "pfo_p": float("nan"),
            "match_val": float("nan"), "match_type": match_type, "matched": False,
        })
    for pfo in pfos:
        rows.append({
            "gen_idx": -999, "gen_pdg": -999, "gen_name": "—",
            "gen_E": float("nan"), "gen_p": float("nan"),
            "pfo_idx": pfo.idx, "pfo_pdg": pfo.pdg,
            "pfo_name": _pdg_label(pfo.pdg),
            "pfo_E": round(pfo.energy, 4), "pfo_p": round(pfo.p, 4),
            "match_val": float("nan"), "match_type": match_type, "matched": False,
        })
    return rows

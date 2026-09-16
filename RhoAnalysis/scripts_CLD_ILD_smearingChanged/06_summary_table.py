#!/usr/bin/env python3
"""06_summary_table.py — tabla CLD vs ILD (smearing corregido) de los fits de polarización.

Para cada sin²θ_eff de entrada, modo y detector lee asymmetry_results_def.txt y da:
  * los 4 resultados del fit (A_tau, A_e, g_V/g_A, sin²θ_eff) con su error, junto
    con los valores de entrada (teóricos) usados para generar los pesos;
  * la desviación en sigmas respecto al valor de entrada, (fit − teo)/σ_fit, con el
    error del fit a la luminosidad MC real (fit_MCstat);
  * el error relativo del error del fit de ILD respecto a CLD, σ_ILD/σ_CLD − 1, con
    ambos normalizados a 2M taus generados (fit_2Mtau) para comparar a igual estadística.

Redondeo (tabla markdown): error a 2 cifras significativas hacia arriba; valor
redondeado (half-up) a la última cifra significativa del error. El CSV guarda los
números sin redondear.
"""
import argparse
import csv
import math
import os
from decimal import Decimal, ROUND_CEILING, ROUND_HALF_UP

REPO = "/nfs/cms/arqolmo/TausFCCee"
MODES = ["rho_lep_legacycuts", "ele_nocuts", "ele_lep10", "muon_nocuts", "muon_lep10", "pion_optcuts"]
PARAMS = [("Atau", "A_τ"), ("Ae", "A_e"), ("gv/ga", "g_V/g_A"), ("sin2theta_eff", "sin²θ_eff")]
SCOPES = {
    "MCstat":  "fit a la estadística MC real",
    "2Mtau":   "fit extrapolado a 2M τ generados (misma estadística en ambos detectores)",
    "FCCyear": "fit extrapolado a 7 ab⁻¹ (1 año FCC-ee en el Z)",
}


def theory(sin_eff):
    x = 1.0 - 4.0 * sin_eff
    a = 2.0 * x / (1.0 + x * x)
    return {"Atau": a, "Ae": a, "gv/ga": x, "sin2theta_eff": sin_eff}


def read_asym(path):
    out = {}
    if not os.path.isfile(path):
        return out
    for line in open(path):
        if line.startswith("#"):
            continue
        p = line.split()
        if len(p) == 3:
            out[p[0]] = (p[1], p[2])  # strings: exactos para Decimal
    return out


def round_pair(val, err):
    """(valor, error) -> strings. Error a 2 c.s. hacia arriba; valor al mismo decimal."""
    v, e = Decimal(val), Decimal(err)
    if e <= 0:
        return f"{v}", f"{e}"
    exp = math.floor(math.log10(e))
    q = Decimal(1).scaleb(exp - 1)
    e_r = (e / q).to_integral_value(rounding=ROUND_CEILING) * q
    if e_r >= Decimal(1).scaleb(exp + 1):  # p.ej. 0.0991 -> 0.10
        exp += 1
        q = Decimal(1).scaleb(exp - 1)
        e_r = (e / q).to_integral_value(rounding=ROUND_CEILING) * q
    v_r = v.quantize(q, rounding=ROUND_HALF_UP)
    return f"{v_r:f}", f"{e_r.quantize(q):f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sin-eff", nargs="+", default=["0.2312", "0.2315"])
    ap.add_argument("--binbase", default="Binned_histograms_MDecs/CLD_ILD_smearingChanged_sineff{sin_eff}")
    ap.add_argument("-o", "--outdir", default="Results/CLD_ILD_smearingChanged_Summary")
    ap.add_argument("--scope", choices=list(SCOPES), default="MCstat",
                    help="Fit del que salen valor ± error y nσ: MCstat (estadística MC real), "
                         "2Mtau (extrapolado a 2M τ generados) o FCCyear (7 ab^-1)")
    args = ap.parse_args()
    outdir = os.path.join(REPO, args.outdir)
    os.makedirs(outdir, exist_ok=True)
    scope_sub = f"fit_{args.scope}"

    rows = []
    md = ["# CLD vs ILD (vertex smearing corregido) — fits de polarización (reco)", "",
          "Muestras: CLD = `ztt_2M_smearing`, ILD = `ild_done`. "
          f"Valor ± error: {SCOPES[args.scope]} (`{scope_sub}`), error a 2 c.s. hacia arriba "
          "y valor redondeado a la última cifra del error. "
          f"**nσ** = (fit − entrada)/σ_{args.scope}. "
          "**Δσ/σ_CLD** = σ_ILD/σ_CLD − 1 con ambos detectores normalizados a 2M τ generados "
          "(`fit_2Mtau`); negativo ⇒ ILD mejora.", ""]

    for s in args.sin_eff:
        th = theory(float(s))
        base = os.path.join(REPO, args.binbase.format(sin_eff=s))
        md += [f"## sin²θ_eff de entrada = {s}", "",
               "Valores de entrada: " + ", ".join(
                   f"{lab} = {th[k]:.5f}" if k != "sin2theta_eff" else f"{lab} = {s}"
                   for k, lab in PARAMS), ""]
        head = ["Modo", "Det"]
        for _, lab in PARAMS:
            head += [lab, "nσ"]
        head += [f"Δσ/σ_CLD ({lab})" for _, lab in PARAMS]
        md += ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]

        for mode in MODES:
            fits = {}
            for det in ("CLD", "ILD"):
                fits[det] = {sub: read_asym(os.path.join(base, f"{det}_{mode}", sub, "asymmetry_results_def.txt"))
                             for sub in ("fit_MCstat", "fit_2Mtau", "fit_FCCyear")}
            for det in ("CLD", "ILD"):
                mc, m2 = fits[det]["fit_MCstat"], fits[det]["fit_2Mtau"]
                sel = fits[det][scope_sub]
                if not sel:
                    md.append(f"| {mode} | {det} | " + " | ".join(["—"] * (len(head) - 2)) + " |")
                    continue
                cells = [mode, det]
                rel_cells = []
                for key, lab in PARAMS:
                    val, err = sel[key]
                    pull = (float(val) - th[key]) / float(err)
                    v_r, e_r = round_pair(val, err)
                    cells += [f"{v_r} ± {e_r}", f"{pull:+.2f}"]
                    rel = None
                    if det == "ILD" and m2 and fits["CLD"]["fit_2Mtau"]:
                        rel = float(m2[key][1]) / float(fits["CLD"]["fit_2Mtau"][key][1]) - 1.0
                    rel_cells.append("—" if rel is None else f"{100 * rel:+.1f} %")
                    row = {"sin_eff_input": s, "mode": mode, "detector": det, "param": key,
                           "input": th[key], "fit": float(mc[key][0]) if mc else float(val)}
                    for sc in SCOPES:
                        f_sc = fits[det][f"fit_{sc}"]
                        row[f"err_{sc}"] = float(f_sc[key][1]) if f_sc else ""
                        row[f"pull_{sc}"] = ((float(f_sc[key][0]) - th[key]) / float(f_sc[key][1])
                                             if f_sc else "")
                    row["rel_err_vs_CLD_2Mtau"] = "" if rel is None else rel
                    rows.append(row)
                md.append("| " + " | ".join(cells + rel_cells) + " |")
        md.append("")

    keys = ["sin_eff_input", "mode", "detector", "param", "input", "fit",
            "err_MCstat", "err_2Mtau", "err_FCCyear",
            "pull_MCstat", "pull_2Mtau", "pull_FCCyear", "rel_err_vs_CLD_2Mtau"]
    csv_path = os.path.join(outdir, "polarization_fit_summary.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    md_name = "summary.md" if args.scope == "MCstat" else f"summary_{args.scope}.md"
    md_path = os.path.join(outdir, md_name)
    open(md_path, "w").write("\n".join(md) + "\n")
    print("\n".join(md))
    print(f"\n[OK] {csv_path}\n[OK] {md_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""makeCosBins_MDecs.py — adaptador de salidas MDecs reco → templates BINED legacy.

Toma los TH2 `OptimalReco_vs_CosThetaVis_dec{T}_{cat0}_{cat1}[_variante]` de los
ficheros `HistosMDecs_{id0}_{id1}_<stem>.root` (X = variable óptima reco, Y = cosθ
del visible reco, 100×20), selecciona el canal OBJETIVO + el(los) canal(es) del OTRO
hemisferio, combina por luminosidad (lumi·σ/N_gen) y escribe en formato legacy
(`histo_SIGNAL_{b}`, `histo_SIGNAL_P1/M1_{b}`, `histo_BG_migrations_{b}`,
`histo_BG_{b}`, `..._full`) para reutilizar `fitPolAssym.py` sin cambios.

Esquema (2026-06-09): selección por id reco, señal/fondo por verdad gen.
  - Señal      = dec{slot}_SIGNAL_SIGNAL  (objetivo gen-correcto ∧ otro gen-correcto)
  - BG migrac. = dec{slot}_ALL_ALL − dec{slot}_SIGNAL_SIGNAL  (todo lo no-señal del canal)
  - P1/M1      = variantes de peso reco: per-tau (reco_P1/M1) o corr (reco_corr_P1/M1)

Todas las variables provienen de información RECO (X = var óptima reco, Y = cosθ_vis
reco, pesos reco_*). Reco-puro.

USO:
  python RhoAnalysis/makeCosBins_MDecs.py \
      --sample-dir Results/RhoAnalysis/PolAnalysis_RECO_SM_ALL_FULLSAMPLE_tau_trained0.4_tph0.35_tpi0_n3_g0.0 \
      --target-decay rho --other-decay lep --weights per-tau \
      -o Binned_histograms_MDecs/ -v
"""
import os
import glob
import argparse
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# =====================================================
# Config de eventos (lumi·σ/N_gen) — igual que makeCosBins.py legacy
# =====================================================
EVENT_CONFIG = {
    "Ztt":    {"xsec_pb": 1476.58, "lumi_pb": 6972, "ngen": 10296000},
    "Zqq":    {"xsec_pb": 30170,   "lumi_pb": 6972, "ngen": 1000000},
    "Bhabha": {"xsec_pb": 273500,  "lumi_pb": 6972, "ngen": 5390000},
}

# Mapeo tipo → nombre del histograma legacy que lee fitPolAssym.py
BG_HIST_NAME = {
    "Zqq":    "BG_Zqq",
    "Bhabha": "BG_Bhabha",
}


def compute_weight(event_type, verbose=False):
    cfg = EVENT_CONFIG.get(event_type)
    if cfg is None:
        raise ValueError(f"Event type '{event_type}' no definido en EVENT_CONFIG")
    w = cfg["lumi_pb"] * cfg["xsec_pb"] / cfg["ngen"]
    if verbose:
        print(f"[INFO] {event_type}: weight = {w:.6g} "
              f"(lumi={cfg['lumi_pb']} pb⁻¹, σ={cfg['xsec_pb']} pb, ngen={cfg['ngen']})")
    return w


# id reco por canal (ρ reco = 2; π=0; a1=10; e=-11; µ=-13)
DECAY_ID = {"rho": 2, "pion": 0, "a1": 10, "ele": -11, "muon": -13}
# atajos para --other-decay
OTHER_GROUPS = {"lep": ["ele", "muon"]}  # 'all' se resuelve con lo que exista

BASE2D = "OptimalReco_vs_CosThetaVis"


def _variant_suffix(kind, mode):
    """kind ∈ {nominal,P1,M1}; mode ∈ {per-tau, corr}."""
    if kind == "nominal":
        return ""
    # stem = "reco_corr" if mode == "corr" else "reco"
    # return f"_{stem}_{'P1' if kind == 'P1' else 'M1'}"
    stem = "_corr" if mode =="corr" else ""
    return f"{stem}_{'P1' if kind == 'P1' else 'M1'}"

def _parse_pair_ids(basename, stem):
    """Extrae (id0, id1) de 'HistosMDecs_{id0}_{id1}_{stem}.root'."""
    prefix = "HistosMDecs_"
    if not basename.startswith(prefix):
        return None
    rest = basename[len(prefix):]
    toks = rest.split("_")
    try:
        id0 = int(toks[0])
        id1 = int(toks[1])
    except (ValueError, IndexError):
        return None
    return id0, id1


def find_pair_files(sample_dir, stem, target_id, other_ids, verbose=False):
    """Devuelve [(path, slot)] donde slot es el índice dec del canal objetivo.

    other_ids = None  → 'all': cualquier pareja que contenga el objetivo.
    """
    out = []
    pattern = os.path.join(sample_dir, f"HistosMDecs_*_{stem}.root")
    for path in sorted(glob.glob(pattern)):
        ids = _parse_pair_ids(os.path.basename(path), stem)
        if ids is None:
            continue
        id0, id1 = ids
        if target_id not in (id0, id1):
            continue
        # determinar el id del otro hemisferio
        other = id1 if id0 == target_id else id0
        if other_ids is not None and other not in other_ids:
            continue
        # slot del objetivo (si ambos == objetivo, par simétrico → tomar dec0 y avisar)
        if id0 == target_id and id1 == target_id:
            if verbose:
                print(f"[WARN] par simétrico {os.path.basename(path)}: uso slot dec0 "
                      f"(posible doble-conteo del objetivo).")
            slot = 0
        else:
            slot = 0 if id0 == target_id else 1
        out.append((path, slot))
        if verbose:
            print(f"[INFO] incluyo {os.path.basename(path)} (otro={other}, slot=dec{slot})")
    return out


def load_sum_2d(files_slots, cat, suffix):
    """Suma el TH2 `{BASE2D}_dec{slot}_{cat}{suffix}` sobre los ficheros dados."""
    h = None
    for path, slot in files_slots:
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            raise RuntimeError(f"No se pudo abrir {path}")
        name = f"{BASE2D}_dec{slot}_{cat}{suffix}"
        h2 = f.Get(name)
        if not h2:
            f.Close()
            raise RuntimeError(f"Histograma '{name}' no encontrado en {path}")
        h2 = h2.Clone(f"_tmp_{cat}{suffix}_{slot}")
        h2.SetDirectory(0)
        f.Close()
        if h is None:
            h = h2
        else:
            h.Add(h2)
    if h is None:
        raise RuntimeError(f"Ningún fichero aportó '{cat}{suffix}'")
    return h


def main(sample_dir, stem, target_decay, other_decays, weight_modes,
         signal_type="Ztt", background_dirs=None, background_types=None,
         nBins=None, rebin=1, bg_def="ss", outdir=".", verbose=False):
    # bg_def="ss":       BG = ALL_ALL − SIGNAL_SIGNAL  (ambos hemisferios gen-correctos = señal)
    # bg_def="sa":       BG = ALL_ALL − SIGNAL_ALL      (solo el hemisferio objetivo gen-correcto)
    #                    SIGNAL_ALL = SIGNAL_SIGNAL + SIGNAL_BG; esos eventos desaparecen del total
    # bg_def="rho_only": Replica el comportamiento legacy: señal = cualquier evento donde
    #                    el hemisferio objetivo (ρ) es gen-correcto, sin importar el leptón.
    #                    Signal = SIGNAL_SIGNAL + SIGNAL_BG
    #                    BG     = ALL_ALL − SIGNAL_SIGNAL − SIGNAL_BG  (ídem _sa)
    # background_dirs/background_types: lista de (dir, tipo) para fondos físicos externos
    #   (Bhabha, Zqq). Se toma ALL_ALL del dir externo, se escala con su EVENT_CONFIG,
    #   y se escribe como histo_BG_{tipo}_{b} + se suma al histo_BG total.

    os.makedirs(outdir, exist_ok=True)
    target_id = DECAY_ID[target_decay]

    # resolver other_ids
    if "all" in other_decays:
        other_ids = None
        other_tag = "all"
    else:
        names = []
        for od in other_decays:
            names.extend(OTHER_GROUPS.get(od, [od]))
        other_ids = sorted({DECAY_ID[n] for n in names})
        other_tag = "_".join(other_decays)

    files_slots = find_pair_files(sample_dir, stem, target_id, other_ids, verbose=verbose)
    if not files_slots:
        raise RuntimeError(f"No se encontraron ficheros para objetivo={target_decay} "
                           f"otro={other_decays} en {sample_dir}")

    w_signal = compute_weight(signal_type, verbose)

    # nBins por defecto = nº de bins en Y del 2D de señal
    sig_nom = load_sum_2d(files_slots, "SIGNAL_SIGNAL", "")
    all_nom = load_sum_2d(files_slots, "ALL_ALL", "")
    nY = sig_nom.GetNbinsY()
    if nBins is None:
        nBins = nY
    if nY % nBins != 0:
        raise ValueError(f"nBins={nBins} no divide nY={nY}")
    bin_length = nY // nBins
    if verbose:
        print(f"[INFO] nY={nY}, nBins={nBins}, bin_length={bin_length}, "
              f"señal total={sig_nom.Integral():.1f}, ALL total={all_nom.Integral():.1f}")

    # Fondos externos (Bhabha, Zqq): cargar ALL_ALL y escalar a luminosidad
    bg_inputs = []   # lista de (tag, h2d_escalado)
    if background_dirs:
        bg_types = background_types or []
        if len(bg_types) != len(background_dirs):
            raise ValueError("--bg-dirs y --bg-types deben tener el mismo número de elementos")
        for bg_dir, bg_type in zip(background_dirs, bg_types):
            if bg_type not in EVENT_CONFIG:
                raise ValueError(f"Tipo '{bg_type}' no definido en EVENT_CONFIG: {list(EVENT_CONFIG)}")
            if bg_type not in BG_HIST_NAME:
                raise ValueError(f"Tipo '{bg_type}' no tiene nombre de histograma legacy en BG_HIST_NAME")
            bg_files = find_pair_files(bg_dir, stem, target_id, other_ids, verbose=verbose)
            if not bg_files:
                raise RuntimeError(f"No se encontraron ficheros de fondo '{bg_type}' en {bg_dir}")
            w_bg = compute_weight(bg_type, verbose)
            h2_bg = load_sum_2d(bg_files, "ALL_ALL", "")
            h2_bg.Scale(w_bg)
            bg_inputs.append((BG_HIST_NAME[bg_type], h2_bg))
            if verbose:
                print(f"[INFO] BG externo '{bg_type}': {h2_bg.Integral():.1f} eventos escalados")

    for mode in weight_modes:
        out_path = os.path.join(
            outdir, f"BINED_MDecs_{target_decay}_{other_tag}_{mode}_{bg_def}.root")
        outfile = ROOT.TFile(out_path, "RECREATE")
        if verbose:
            print(f"[INFO] === modo pesos '{mode}' → {out_path} ===")

        # cargar las variantes de señal (nominal + P1 + M1) una sola vez.
        # CONVENCIÓN DE PESOS (Alcaraz/weightsPol, confirmada por closure): el sufijo
        # del peso indica el escenario de ORIGEN, no el destino. El peso reco_P1
        # (newAtau New_Atau=+1) aplicado al SM IMITA el escenario A_τ=-1, y reco_M1
        # imita A_τ=+1. fitPolAssym define hist_p1 como el template P_τ=+1, así que
        # SE INVIERTEN: histo_SIGNAL_P1 ← peso reco_M1, histo_SIGNAL_M1 ← peso reco_P1.
        sig = {
            "SIGNAL":    load_sum_2d(files_slots, "SIGNAL_SIGNAL", _variant_suffix("nominal", mode)),
            "SIGNAL_P1": load_sum_2d(files_slots, "SIGNAL_SIGNAL", _variant_suffix("M1", mode)),
            "SIGNAL_M1": load_sum_2d(files_slots, "SIGNAL_SIGNAL", _variant_suffix("P1", mode)),
        }
        # BG migraciones (nominal)
        bg_mig_2d = all_nom.Clone("_bg_mig_2d")
        bg_mig_2d.SetDirectory(0)
        bg_mig_2d.Add(sig["SIGNAL"], -1.0)
        if bg_def in ("sa", "rho_only"):
            # Restar SIGNAL_BG → BG = solo eventos donde el hemisferio objetivo es BG gen
            sig_bg_nom = load_sum_2d(files_slots, "SIGNAL_BG", "")
            bg_mig_2d.Add(sig_bg_nom, -1.0)
            if bg_def == "rho_only":
                # Añadir SIGNAL_BG a los templates de señal: replica el comportamiento
                # legacy (señal = ρ gen correcto, sin importar verdad gen del leptón).
                sig["SIGNAL"].Add(sig_bg_nom)
                if verbose:
                    print(f"[INFO] rho_only: SIGNAL_BG integral = {sig_bg_nom.Integral():.1f} "
                          f"(añadido a señal)")
                # P1/M1: cargar variante de peso de SIGNAL_BG si existe; si no, usar nominal
                for sig_key, w_suffix in [("SIGNAL_P1", _variant_suffix("M1", mode)),
                                          ("SIGNAL_M1", _variant_suffix("P1", mode))]:
                    try:
                        sig_bg_w = load_sum_2d(files_slots, "SIGNAL_BG", w_suffix)
                        sig[sig_key].Add(sig_bg_w)
                    except RuntimeError:
                        sig[sig_key].Add(sig_bg_nom)
        if verbose:
            print(f"[INFO] bg_def='{bg_def}': Signal={sig['SIGNAL'].Integral():.1f}, "
                  f"BG integral = {bg_mig_2d.Integral():.1f}")

        def _proj_write(h2, name, b_lo, b_hi, scale):
            h1 = h2.ProjectionX(name, b_lo, b_hi)
            if rebin > 1:
                h1.Rebin(rebin)
            h1.Scale(scale)
            h1.SetXTitle("optimal reco var")
            outfile.cd()
            h1.Write(name)
            return h1

        # ---- por bin de cosθ_vis ----
        for b in range(nBins):
            b_lo = b * bin_length + 1
            b_hi = (b + 1) * bin_length
            for key, h2 in sig.items():
                _proj_write(h2, f"histo_{key}_{b}", b_lo, b_hi, w_signal)
            h_mig = _proj_write(bg_mig_2d, f"histo_BG_migrations_{b}", b_lo, b_hi, w_signal)
            # BG externos: escribir individualmente y acumular en total
            h_tot = h_mig.Clone(f"histo_BG_{b}")
            for bg_tag, h2_ext in bg_inputs:
                h_ext = _proj_write(h2_ext, f"histo_{bg_tag}_{b}", b_lo, b_hi, 1.0)
                h_tot.Add(h_ext)
            outfile.cd()
            h_tot.Write(f"histo_BG_{b}")

        # ---- full range ----
        for key, h2 in sig.items():
            _proj_write(h2, f"histo_{key}_full", 1, nY, w_signal)
        h_mig_full = _proj_write(bg_mig_2d, "histo_BG_migrations_full", 1, nY, w_signal)
        h_tot_full = h_mig_full.Clone("histo_BG_full")
        for bg_tag, h2_ext in bg_inputs:
            h_ext_full = _proj_write(h2_ext, f"histo_{bg_tag}_full", 1, nY, 1.0)
            h_tot_full.Add(h_ext_full)
        outfile.cd()
        h_tot_full.Write("histo_BG_full")

        outfile.Close()
        if verbose:
            print(f"[OK] escrito {out_path}")

    if verbose:
        print("[OK] makeCosBins_MDecs completado.")


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        "Binning angular (cosθ_vis) de salidas MDecs reco → templates BINED legacy")
    p.add_argument("--sample-dir", required=True,
                   help="Dir con los HistosMDecs_{id0}_{id1}_<stem>.root de la señal")
    p.add_argument("--stem", default="tau_traineddecayAll_0.4_tph0.35_tpi0_n3_g0.0",
                   help="Sufijo del nombre de fichero tras los ids de par")
    p.add_argument("--target-decay", default="rho", choices=list(DECAY_ID.keys()),
                   help="Canal objetivo (hemisferio que se mide)")
    p.add_argument("--other-decay", nargs="+", default=["lep"],
                   help="Canal(es) del otro hemisferio: rho/pion/a1/ele/muon, "
                        "atajo 'lep' (ele+muon) o 'all' (cualquiera)")
    p.add_argument("--weights", nargs="+", default=["per-tau"],
                   choices=["per-tau", "corr"],
                   help="Modo(s) de peso para los templates P1/M1: per-tau "
                        "(reco_P1/M1) y/o corr (reco_corr_P1/M1, término cruzado)")
    p.add_argument("--signal-type", default="Ztt")
    p.add_argument("--bg-dirs", nargs="+", default=None,
                   help="Dir(s) con HistosMDecs de fondos físicos externos (Bhabha, Zqq). "
                        "Mismo número de elementos que --bg-types.")
    p.add_argument("--bg-types", nargs="+", default=None,
                   help="Tipo(s) de evento para cada --bg-dirs, en el mismo orden. "
                        f"Opciones: {list(BG_HIST_NAME.keys())}. "
                        "Ejemplo: --bg-dirs /ruta/Bhabha --bg-types Bhabha")
    p.add_argument("--bg-def", default="ss", choices=["ss", "sa", "rho_only"],
                   help="Definición del BG: "
                        "'ss' (default) = ALL_ALL − SIGNAL_SIGNAL; "
                        "'sa' = ALL_ALL − SIGNAL_SIGNAL − SIGNAL_BG "
                        "(SIGNAL_BG excluido del total); "
                        "'rho_only' = replica legacy: señal = SIGNAL_SIGNAL + SIGNAL_BG "
                        "(ρ gen correcto sin importar el leptón), "
                        "BG = ALL_ALL − SIGNAL_SIGNAL − SIGNAL_BG")
    p.add_argument("--nBins", type=int, default=None,
                   help="Nº bins en cosθ_vis (def = nº bins Y del 2D, normalmente 20)")
    p.add_argument("--rebin", type=int, default=1)
    p.add_argument("--base2d", default=BASE2D,
                   help="Nombre base del TH2 observable×cosθ_vis a binar. "
                        f"Default: {BASE2D}. Usa 'OptimalNN_vs_CosThetaVis' para "
                        "la variable óptima del ρ vía red neuronal (modelo v2).")
    p.add_argument("-o", "--outdir", default="./Binned_histograms_MDecs/")
    p.add_argument("-v", "--verbose", action="store_true")

    args = p.parse_args()
    # Permite cambiar el observable 2D (p.ej. OptimalNN_vs_CosThetaVis) sin tocar
    # load_sum_2d, que lee el global BASE2D.
    BASE2D = args.base2d
    main(
        sample_dir=args.sample_dir,
        stem=args.stem,
        target_decay=args.target_decay,
        other_decays=args.other_decay,
        weight_modes=args.weights,
        signal_type=args.signal_type,
        background_dirs=args.bg_dirs,
        background_types=args.bg_types,
        nBins=args.nBins,
        rebin=args.rebin,
        bg_def=args.bg_def,
        outdir=args.outdir,
        verbose=args.verbose,
    )

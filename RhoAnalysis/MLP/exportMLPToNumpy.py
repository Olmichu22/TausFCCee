"""
Exporta el MLP de polarización (PyTorch) + el StandardScaler a un .npz plano,
para poder hacer la inferencia en numpy puro dentro del entorno de análisis
(key4hep/ROOT), que NO tiene torch.

Se ejecuta UNA sola vez en un entorno con torch + joblib + sklearn (p.ej. el
conda env `SDHCAL`):

    python RhoAnalysis/MLP/exportMLPToNumpy.py \
        --model-dir MLPolResults/train_pol_results_optimal_v2

Genera <model-dir>/mlp_numpy.npz con:
  - meta:        input_dim, n_layers, feature_names
  - scaler_mean, scaler_scale
  - por capa Linear  i:  layer{i}_type='linear', layer{i}_W, layer{i}_b
  - por capa BatchNorm i: layer{i}_type='batchnorm', layer{i}_gamma, _beta,
                          _running_mean, _running_var, _eps
  - por activación   i:  layer{i}_type='relu' | 'sigmoid'

El forward equivalente (numpy) lo implementa modules/mlpPolInference.py, que
debe consumir exactamente este formato.
"""
import os
import sys
import argparse

import numpy as np
import torch
import joblib

# Reutiliza la MISMA definición de arquitectura del entrenamiento.
from RhoAnalysis.MLP.MLOptimalObservable import TauPolarizationMLP  # noqa: E402

# Orden EXACTO de las 14 features con el que se entrenó el modelo
# (ver RhoAnalysis/MLP/createPolDatasets.py:106-119,173-178). Solo es metadato
# (las features son posicionales); el modelo por defecto (optimal_lepton) usa el
# PIÓN en las 4 primeras, el v2 usaba el meson visible.
FEATURE_NAMES = [
    "PionE", "PionTheta", "PionPhi", "PionP",
    "photon1_E", "photon1_theta", "photon1_phi",
    "photon2_E", "photon2_theta", "photon2_phi",
    "lepE", "lepTheta", "lepPhi", "lepP",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-dir", default="MLPolResults/train_pol_results_optimal_v2",
                    help="Directorio con tau_polarization_mlp.pt y scaler.pkl")
    ap.add_argument("--checkpoint", default="tau_polarization_mlp.pt")
    ap.add_argument("--scaler", default="scaler.pkl")
    ap.add_argument("--out", default="mlp_numpy.npz")
    args = ap.parse_args()

    ckpt_path   = os.path.join(args.model_dir, args.checkpoint)
    scaler_path = os.path.join(args.model_dir, args.scaler)
    out_path    = os.path.join(args.model_dir, args.out)

    for p in (ckpt_path, scaler_path):
        if not os.path.isfile(p):
            print(f"[ERROR] No existe: {p}")
            sys.exit(1)

    # ── Modelo ──────────────────────────────────────────────────────────────
    # weights_only=False: el checkpoint es nuestro y contiene escalares numpy
    # (test_loss, etc.) además del state_dict; torch>=2.6 lo bloquea por defecto.
    ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    hp = ckpt["hyperparams"]
    model = TauPolarizationMLP(
        input_dim=hp["input_dim"],
        hidden_dims=hp["hidden_dims"],
        dropout=hp.get("dropout", 0.0),
    )
    model.load_state_dict(ckpt["model_state_dict"])
    model.eval()

    if hp["input_dim"] != len(FEATURE_NAMES):
        print(f"[ERROR] input_dim={hp['input_dim']} != {len(FEATURE_NAMES)} features esperadas")
        sys.exit(1)

    # ── Scaler ──────────────────────────────────────────────────────────────
    scaler = joblib.load(scaler_path)
    scaler_mean  = np.asarray(scaler.mean_,  dtype=np.float64)
    scaler_scale = np.asarray(scaler.scale_, dtype=np.float64)

    # ── Vuelco capa a capa (model.net es un nn.Sequential) ──────────────────
    out = {
        "input_dim": np.array(hp["input_dim"]),
        "feature_names": np.array(FEATURE_NAMES),
        "scaler_mean": scaler_mean,
        "scaler_scale": scaler_scale,
    }
    types = []
    for i, layer in enumerate(model.net):
        if isinstance(layer, torch.nn.Linear):
            types.append("linear")
            out[f"layer{i}_W"] = layer.weight.detach().numpy().astype(np.float64)
            out[f"layer{i}_b"] = layer.bias.detach().numpy().astype(np.float64)
        elif isinstance(layer, torch.nn.BatchNorm1d):
            types.append("batchnorm")
            out[f"layer{i}_gamma"]        = layer.weight.detach().numpy().astype(np.float64)
            out[f"layer{i}_beta"]         = layer.bias.detach().numpy().astype(np.float64)
            out[f"layer{i}_running_mean"] = layer.running_mean.detach().numpy().astype(np.float64)
            out[f"layer{i}_running_var"]  = layer.running_var.detach().numpy().astype(np.float64)
            out[f"layer{i}_eps"]          = np.array(layer.eps, dtype=np.float64)
        elif isinstance(layer, torch.nn.ReLU):
            types.append("relu")
        elif isinstance(layer, torch.nn.Sigmoid):
            types.append("sigmoid")
        elif isinstance(layer, torch.nn.Dropout):
            # En eval() es identidad: se omite del grafo numpy.
            types.append("identity")
        else:
            print(f"[ERROR] Capa no soportada en posición {i}: {type(layer).__name__}")
            sys.exit(1)

    out["layer_types"] = np.array(types)
    out["n_layers"] = np.array(len(types))

    np.savez(out_path, **out)
    print(f"[OK] Exportado modelo numpy a: {out_path}")
    print(f"     capas: {types}")
    print(f"     features ({len(FEATURE_NAMES)}): {FEATURE_NAMES}")


if __name__ == "__main__":
    main()

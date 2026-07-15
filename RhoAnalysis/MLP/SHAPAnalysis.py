"""
SHAPAnalysis.py
Carga el modelo óptimo entrenado por OptimizeMLObservable.py y genera
análisis de importancia y dependencia de variables usando SHAP.

Requiere: shap >= 0.44, torch, numpy, matplotlib, joblib
"""

import os
import argparse
import numpy as np
import torch
import joblib
import shap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from RhoAnalysis.MLP.MLOptimalObservable import TauPolarizationMLP, load_data, prepare_data

# ── Configuración ──────────────────────────────────────────────────────────────

# FEATURE_NAMES = [
    # 'MesonE',      'MesonTheta',  'MesonPhi',   'MesonP',
    # 'gamma1_E',    'gamma1_theta','gamma1_phi',
    # 'gamma2_E',    'gamma2_theta','gamma2_phi',
    # 'lepE',        'lepTheta',    'lepPhi',      'lepP',
# ]

FEATURE_NAMES = [
    'PionE',      'PionTheta',  'PionPhi',   'PionP',
    'gamma1_E',    'gamma1_theta','gamma1_phi',
    'gamma2_E',    'gamma2_theta','gamma2_phi',
    'lepE',        'lepTheta',    'lepPhi',      'lepP',
]

DEFAULT_CHECKPOINT = "MLPolResults/train_pol_results_optimal_lepton/tau_polarization_mlp.pt"
DEFAULT_SCALER     = "MLPolResults/train_pol_results_optimal_lepton/scaler.pkl"
DEFAULT_OUTPUT_DIR = "MLPolResults/train_pol_results_optimal_lepton/shap"
DEFAULT_P1_PATH    = "MLPolResults/datasets/train_data_pol1_reco_lepton.npz"
DEFAULT_M1_PATH    = "MLPolResults/datasets/train_data_pol-1_reco_lepton.npz"

N_BACKGROUND = 500   # muestras para el background del explainer
N_EXPLAIN    = 1000  # muestras sobre las que calcular SHAP values
TOP_FEATURES = 5     # dependence plots de las N variables más importantes


# ── Wrapper para SHAP (CPU, modo eval, output escalar) ─────────────────────────

class ModelWrapper(torch.nn.Module):
    """Envuelve el modelo; GradientExplainer necesita output 2-D (batch, n_outputs)."""

    def __init__(self, model):
        super().__init__()
        self.model = model

    def forward(self, x):
        return self.model(x)  # shape: (batch, 1)


# ── Carga del modelo ───────────────────────────────────────────────────────────

def load_model(checkpoint_path: str, device: str):
    """Devuelve (model, ckpt) para poder acceder al test set guardado."""
    ckpt = torch.load(checkpoint_path, map_location=device)
    hp   = ckpt['hyperparams']
    model = TauPolarizationMLP(
        input_dim=hp['input_dim'],
        hidden_dims=hp['hidden_dims'],
        dropout=hp['dropout'],
    )
    model.load_state_dict(ckpt['model_state_dict'])
    model.eval()
    return model, ckpt


# ── Plots ──────────────────────────────────────────────────────────────────────

def plot_summary(shap_values, X_scaled, output_dir):
    """Beeswarm + bar de importancia media."""
    shap_exp = shap.Explanation(
        values=shap_values,
        data=X_scaled,
        feature_names=FEATURE_NAMES,
    )

    # Beeswarm
    fig, ax = plt.subplots(figsize=(9, 6))
    shap.plots.beeswarm(shap_exp, max_display=len(FEATURE_NAMES), show=False)
    plt.title('SHAP — impacto de cada variable en el output', pad=12)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'shap_beeswarm.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)

    # Bar (importancia media)
    fig, ax = plt.subplots(figsize=(8, 5))
    shap.plots.bar(shap_exp, max_display=len(FEATURE_NAMES), show=False)
    plt.title('SHAP — importancia media |SHAP|', pad=12)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'shap_bar.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)


def plot_dependence(shap_values, X_scaled, scaler, output_dir, top_n):
    """Dependence plots para las top_n variables más importantes."""
    mean_abs = np.abs(shap_values).mean(axis=0)
    top_idx  = np.argsort(mean_abs)[::-1][:top_n]

    X_orig = scaler.inverse_transform(X_scaled)

    for rank, idx in enumerate(top_idx):
        feat_name = FEATURE_NAMES[idx]
        fig, ax = plt.subplots(figsize=(7, 5))
        sc = ax.scatter(
            X_orig[:, idx], shap_values[:, idx],
            c=X_orig[:, idx], cmap='coolwarm', alpha=0.4, s=6, linewidths=0,
        )
        plt.colorbar(sc, ax=ax, label=feat_name)
        ax.axhline(0, color='black', linewidth=0.8, linestyle='--')
        ax.set_xlabel(feat_name)
        ax.set_ylabel(f'SHAP value ({feat_name})')
        ax.set_title(f'SHAP dependence — {feat_name} (importancia #{rank+1})')
        ax.grid(True, alpha=0.25)
        plt.tight_layout()
        fig.savefig(
            os.path.join(output_dir, f'shap_dependence_{rank+1:02d}_{feat_name}.png'),
            dpi=150, bbox_inches='tight',
        )
        plt.close(fig)


def plot_heatmap(shap_values, output_dir):
    """Heatmap de correlaciones entre SHAP values de distintas variables."""
    corr = np.corrcoef(shap_values.T)
    fig, ax = plt.subplots(figsize=(9, 7))
    im = ax.imshow(corr, cmap='RdBu_r', vmin=-1, vmax=1)
    plt.colorbar(im, ax=ax, label='Correlación')
    ax.set_xticks(range(len(FEATURE_NAMES))); ax.set_xticklabels(FEATURE_NAMES, rotation=45, ha='right')
    ax.set_yticks(range(len(FEATURE_NAMES))); ax.set_yticklabels(FEATURE_NAMES)
    ax.set_title('Correlación entre SHAP values')
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, 'shap_correlation.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)


# ── Main ───────────────────────────────────────────────────────────────────────

def main(args):
    os.makedirs(args.output_dir, exist_ok=True)
    device = 'cpu'  # GradientExplainer requiere tensores en CPU

    # ── Datos — mismo seed=42 que en OptimizeMLObservable ─────────────────────
    print("Cargando datos…")
    X_p1, w_p1, X_m1, w_m1 = load_data(args.p1_path, args.m1_path)
    (X_tr, y_tr, w_tr), _, (X_te, y_te, w_te), _ = prepare_data(X_p1, X_m1, w_p1, w_m1)

    # ── Scaler ─────────────────────────────────────────────────────────────────
    if not os.path.isfile(args.scaler):
        raise FileNotFoundError(
            f"No se encontró el scaler en '{args.scaler}'. "
            "Ejecuta primero OptimizeMLObservable.py para generarlo."
        )
    scaler = joblib.load(args.scaler)
    print(f"Scaler cargado desde '{args.scaler}'")

    print(f"Cargando modelo desde '{args.checkpoint}'…")
    model, ckpt = load_model(args.checkpoint, device)
    wrapper     = ModelWrapper(model).to(device)

    rng = np.random.default_rng(42)

    # Background: subconjunto del train set (ya escalado)
    bg_idx = rng.choice(len(X_tr), size=min(args.n_background, len(X_tr)), replace=False)
    X_bg   = torch.tensor(X_tr[bg_idx], dtype=torch.float32)

    # Subconjunto del test a explicar
    exp_idx = rng.choice(len(X_te), size=min(args.n_explain, len(X_te)), replace=False)
    X_exp   = torch.tensor(X_te[exp_idx], dtype=torch.float32)

    print(f"Calculando SHAP values — background: {len(X_bg)}, explain: {len(X_exp)}…")
    explainer   = shap.GradientExplainer(wrapper, X_bg)
    raw = explainer.shap_values(X_exp)
    # Versiones antiguas devuelven lista [array(n,f)]; nuevas devuelven array(n,f) directamente
    shap_values = np.array(raw[0] if isinstance(raw, list) else raw)

    print(f"SHAP values calculados — shape: {shap_values.shape}")

    # Guardar para reutilización
    np.savez_compressed(
        os.path.join(args.output_dir, 'shap_values.npz'),
        shap_values=shap_values,
        X_exp=X_exp.numpy(),
        exp_idx=exp_idx,
    )

    X_exp_np = X_exp.numpy()
    print("Generando plots…")
    plot_summary(shap_values, X_exp_np, args.output_dir)
    plot_dependence(shap_values, X_exp_np, scaler, args.output_dir, args.top_features)
    plot_heatmap(shap_values, args.output_dir)

    # Ranking de importancia en consola
    mean_abs = np.abs(shap_values).mean(axis=0)
    ranking  = sorted(zip(FEATURE_NAMES, mean_abs), key=lambda x: x[1], reverse=True)
    print("\nImportancia media |SHAP|:")
    for rank, (name, val) in enumerate(ranking, 1):
        bar = '█' * int(val / ranking[0][1] * 30)
        print(f"  {rank:2d}. {name:<15s}  {val:.4f}  {bar}")

    print(f"\nTodo guardado en '{args.output_dir}/'")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Análisis SHAP del modelo de polarización tau")
    parser.add_argument('--checkpoint',   default=DEFAULT_CHECKPOINT)
    parser.add_argument('--scaler',       default=DEFAULT_SCALER)
    parser.add_argument('--output-dir',   default=DEFAULT_OUTPUT_DIR)
    parser.add_argument('--p1-path',      default=DEFAULT_P1_PATH)
    parser.add_argument('--m1-path',      default=DEFAULT_M1_PATH)
    parser.add_argument('--n-background', type=int, default=N_BACKGROUND,
                        help='Nº muestras de background para GradientExplainer')
    parser.add_argument('--n-explain',    type=int, default=N_EXPLAIN,
                        help='Nº muestras del test set a explicar')
    parser.add_argument('--top-features', type=int, default=TOP_FEATURES,
                        help='Nº de dependence plots a generar (top variables)')
    args = parser.parse_args()
    main(args)

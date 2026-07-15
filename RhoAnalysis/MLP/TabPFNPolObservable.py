"""
TabPFN-based tau polarization classifier.
Mirrors MLOptimalObservable.py for direct comparison.

Usage (GPU node):
    export TABPFN_TOKEN=<your_token>
    python RhoAnalysis/TabPFNPolObservable.py \
        [--p1 datasets/train_data_pol1.npz] \
        [--m1 datasets/train_data_pol-1.npz] \
        [--output train_tabpfn_results/] \
        [--max-samples 50000] \
        [--mlp-roc train_pol_results/roc_curve_data.npz]
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc

FEATURE_NAMES = [
    'MesonE', 'MesonTheta', 'MesonPhi', 'MesonP',
    'gamma1_E', 'gamma1_theta', 'gamma1_phi',
    'gamma2_E', 'gamma2_theta', 'gamma2_phi',
    'lepE', 'lepTheta', 'lepPhi', 'lepP',
]


# ── Data ───────────────────────────────────────────────────────────────────────

def load_data(p1_path, m1_path):
    p1 = np.load(p1_path)["data"]
    m1 = np.load(m1_path)["data"]
    print(f"Cargados — P1: {len(p1):,}  M1: {len(m1):,}  Total: {len(p1)+len(m1):,}")
    return p1, m1


def subsample(p1, m1, max_total, seed=42):
    """Stratified subsample to stay within TabPFN row limit."""
    n_each = max_total // 2
    rng = np.random.default_rng(seed)
    if len(p1) > n_each:
        p1 = p1[rng.choice(len(p1), n_each, replace=False)]
    if len(m1) > n_each:
        m1 = m1[rng.choice(len(m1), n_each, replace=False)]
    print(f"Submuestreo → P1: {len(p1):,}  M1: {len(m1):,}")
    return p1, m1


def prepare_splits(p1, m1, test_size=0.30, val_size=0.50, seed=42):
    X = np.concatenate([p1, m1], axis=0)
    y = np.concatenate([np.ones(len(p1)), np.zeros(len(m1))])

    # Class-balance weights
    n_p1, n_m1 = len(p1), len(m1)
    w = np.concatenate([np.full(n_p1, n_m1 / n_p1), np.full(n_m1, n_p1 / n_m1)])
    w[y == 1] *= w[y == 0].sum() / w[y == 1].sum()

    X_tr, X_tmp, y_tr, y_tmp, w_tr, w_tmp = train_test_split(
        X, y, w, test_size=test_size, stratify=y, random_state=seed
    )
    X_val, X_te, y_val, y_te, w_val, w_te = train_test_split(
        X_tmp, y_tmp, w_tmp, test_size=val_size, stratify=y_tmp, random_state=seed
    )

    print(f"Splits — train: {len(y_tr):,}  val: {len(y_val):,}  test: {len(y_te):,}")
    return (X_tr, y_tr, w_tr), (X_val, y_val, w_val), (X_te, y_te, w_te)


# ── Plots ──────────────────────────────────────────────────────────────────────

def plot_output_distribution(preds, labels, weights, output_dir):
    fig, ax = plt.subplots(figsize=(7, 5))
    bins = np.linspace(0, 1, 50)
    ax.hist(preds[labels == 1], bins=bins, weights=weights[labels == 1],
            alpha=0.6, label=r'P1 ($\mathcal{P}_\tau=+1$)', density=True)
    ax.hist(preds[labels == 0], bins=bins, weights=weights[labels == 0],
            alpha=0.6, label=r'M1 ($\mathcal{P}_\tau=-1$)', density=True)
    ax.set_xlabel('Score TabPFN')
    ax.set_ylabel('Densidad')
    ax.set_title('Distribución del discriminante TabPFN (test set)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'output_distribution.png'), dpi=150)
    plt.close()


def plot_roc(fpr, tpr, roc_auc, output_dir):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(fpr, tpr, label=f'TabPFN  (AUC = {roc_auc:.4f})')
    ax.plot([0, 1], [0, 1], 'k--', label='Azar')
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('ROC curve TabPFN (test set)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'roc_curve.png'), dpi=150)
    plt.close()


def plot_output_vs_features(X_te, preds, labels, output_dir):
    fig, axes = plt.subplots(4, 4, figsize=(16, 14))
    axes = axes.flatten()
    for i, name in enumerate(FEATURE_NAMES):
        axes[i].scatter(X_te[:, i], preds, alpha=0.05, s=1, c=labels, cmap='bwr')
        axes[i].set_xlabel(name)
        axes[i].set_ylabel('Score TabPFN')
        axes[i].set_title(f'Score vs {name}')
        axes[i].grid(True, alpha=0.2)
    axes[-1].set_visible(False)
    plt.suptitle('Correlación score TabPFN vs variables de entrada (test set)', y=1.01)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'output_vs_features.png'), dpi=150, bbox_inches='tight')
    plt.close()


def plot_comparison_roc(fpr_tab, tpr_tab, auc_tab, mlp_roc_path, output_dir):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(fpr_tab, tpr_tab, label=f'TabPFN  (AUC = {auc_tab:.4f})')

    if mlp_roc_path and os.path.exists(mlp_roc_path):
        d = np.load(mlp_roc_path)
        ax.plot(d['fpr'], d['tpr'], label=f'MLP     (AUC = {d["auc"]:.4f})', linestyle='--')
        print(f"MLP AUC cargado desde {mlp_roc_path}: {d['auc']:.4f}")
    else:
        print("No se encontró ROC del MLP — solo se dibuja TabPFN")

    ax.plot([0, 1], [0, 1], 'k--', label='Azar')
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('Comparación ROC: TabPFN vs MLP')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'comparison_roc.png'), dpi=150)
    plt.close()


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--p1',          default='datasets/train_data_pol1.npz')
    parser.add_argument('--m1',          default='datasets/train_data_pol-1.npz')
    parser.add_argument('--output',      default='train_tabpfn_results')
    parser.add_argument('--max-samples', type=int, default=50000,
                        help='Max total events (TabPFN limit ~50k)')
    parser.add_argument('--mlp-roc',     default='train_pol_results/roc_curve_data.npz',
                        help='Path to MLP ROC data (npz with fpr/tpr/auc) for comparison plot')
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    # ── Check token ────────────────────────────────────────────────────────────
    if not os.environ.get('TABPFN_TOKEN'):
        print("AVISO: TABPFN_TOKEN no definido — TabPFN pedirá autenticación interactiva")

    # ── Load TabPFN ────────────────────────────────────────────────────────────
    print("Importando TabPFN...")
    from tabpfn import TabPFNClassifier

    # ── Data ───────────────────────────────────────────────────────────────────
    p1, m1 = load_data(args.p1, args.m1)

    if len(p1) + len(m1) > args.max_samples:
        p1, m1 = subsample(p1, m1, args.max_samples)

    (X_tr, y_tr, w_tr), (X_val, y_val, w_val), (X_te, y_te, w_te) = prepare_splits(p1, m1)

    # ── Fit ────────────────────────────────────────────────────────────────────
    import torch
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Device: {device}")

    clf = TabPFNClassifier(device=device, n_estimators=4)
    print("Entrenando TabPFN (fit en train set)...")
    clf.fit(X_tr, y_tr)
    print("Fit completado.")

    # ── Predict on test ────────────────────────────────────────────────────────
    print("Inferencia en test set...")
    preds = clf.predict_proba(X_te)[:, 1]

    # ── Metrics ────────────────────────────────────────────────────────────────
    fpr, tpr, _ = roc_curve(y_te, preds)
    roc_auc = auc(fpr, tpr)
    acc = ((preds > 0.5) == y_te).mean()
    print(f"\nResultados en test set:")
    print(f"  AUC     : {roc_auc:.4f}")
    print(f"  Accuracy: {acc:.4f}")

    # Guardar ROC data para comparación futura con otros modelos
    np.savez(os.path.join(args.output, 'roc_curve_data.npz'),
             fpr=fpr, tpr=tpr, auc=np.float64(roc_auc))

    # ── Plots ──────────────────────────────────────────────────────────────────
    print("Generando plots...")
    plot_output_distribution(preds, y_te, w_te, args.output)
    plot_roc(fpr, tpr, roc_auc, args.output)
    plot_output_vs_features(X_te, preds, y_te, args.output)
    plot_comparison_roc(fpr, tpr, roc_auc, args.mlp_roc, args.output)

    print(f"\nTodo guardado en '{args.output}/'")
    print(f"  output_distribution.png")
    print(f"  roc_curve.png")
    print(f"  output_vs_features.png")
    print(f"  comparison_roc.png")
    print(f"  roc_curve_data.npz  (para comparaciones futuras)")


if __name__ == "__main__":
    main()

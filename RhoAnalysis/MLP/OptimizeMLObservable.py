# train_optuna.py
# Importa todo lo necesario del script base
from RhoAnalysis.MLP.MLOptimalObservable import (
    TauDataset, TauPolarizationMLP,
    load_data, prepare_data,
    weighted_bce, evaluate, train,
    save_evaluation_plots,
)

import os
import torch
import numpy as np
import joblib
import optuna
from optuna.samplers import TPESampler
from torch.utils.data import DataLoader

# ── Configuración ──────────────────────────────────────────────────────────────

OUTPUT_DIR = "MLPolResults/train_pol_results_optimal_lepton"
P1_PATH    = "MLPolResults/datasets/train_data_pol1_reco_lepton.npz"
M1_PATH    = "MLPolResults/datasets/train_data_pol-1_reco_lepton.npz"
N_TRIALS   = 50

# ── Optuna objective ───────────────────────────────────────────────────────────

def make_objective(X_tr, y_tr, w_tr, X_val, y_val, w_val, device):

    def objective(trial):
        n_layers    = trial.suggest_int('n_layers', 2, 5)
        hidden_dims = [
            trial.suggest_categorical(f'h_{i}', [32, 64, 128, 256])
            for i in range(n_layers)
        ]
        dropout      = trial.suggest_float('dropout',      0.0,  0.5,  step=0.1)
        weight_decay = trial.suggest_float('weight_decay', 1e-5, 1e-2, log=True)
        lr           = trial.suggest_float('lr',           1e-4, 1e-2, log=True)
        batch_size   = trial.suggest_categorical('batch_size', [256, 512, 1024])

        train_loader_t = DataLoader(
            TauDataset(X_tr, y_tr, w_tr), batch_size=batch_size, shuffle=True
        )
        val_loader_t = DataLoader(
            TauDataset(X_val, y_val, w_val), batch_size=batch_size, shuffle=False
        )

        model = TauPolarizationMLP(
            input_dim=14, hidden_dims=hidden_dims, dropout=dropout
        ).to(device)

        optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, patience=5, factor=0.5, min_lr=1e-6
        )

        best_val_loss = np.inf
        epochs_no_imp = 0

        for epoch in range(100):
            model.train()
            for X_batch, y_batch, w_batch in train_loader_t:
                X_batch, y_batch, w_batch = (
                    X_batch.to(device), y_batch.to(device), w_batch.to(device)
                )
                optimizer.zero_grad()
                loss = weighted_bce(model(X_batch), y_batch, w_batch)
                loss.backward()
                optimizer.step()

            val_loss, _ = evaluate(model, val_loader_t, device)
            scheduler.step(val_loss)

            trial.report(val_loss, epoch)
            if trial.should_prune():
                raise optuna.exceptions.TrialPruned()

            if val_loss < best_val_loss:
                best_val_loss = val_loss
                epochs_no_imp = 0
            else:
                epochs_no_imp += 1
                if epochs_no_imp >= 15:
                    break

        return best_val_loss

    return objective


# ── Main ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Usando device: {device}")

    X_p1, w_p1, X_m1, w_m1 = load_data(P1_PATH, M1_PATH)

    (X_tr, y_tr, w_tr), (X_val, y_val, w_val), (X_te, y_te, w_te), scaler = \
        prepare_data(X_p1, X_m1, w_p1, w_m1)

    joblib.dump(scaler, os.path.join(OUTPUT_DIR, 'scaler.pkl'))

    # ── Búsqueda ───────────────────────────────────────────────────────────────

    study = optuna.create_study(
        direction='minimize',
        sampler=TPESampler(seed=42),
        pruner=optuna.pruners.MedianPruner(n_startup_trials=5, n_warmup_steps=10),
        study_name='tau_pol_mlp',
        storage=f'sqlite:///{OUTPUT_DIR}/optuna_tau.db',
        load_if_exists=True,
    )
    study.optimize(
        make_objective(X_tr, y_tr, w_tr, X_val, y_val, w_val, device),
        n_trials=N_TRIALS,
        show_progress_bar=True,
    )

    print("\nMejores hiperparámetros:")
    for k, v in study.best_params.items():
        print(f"  {k}: {v}")

    # ── Reentrenar con mejores params ──────────────────────────────────────────

    best        = study.best_params
    hidden_dims = [best[f'h_{i}'] for i in range(best['n_layers'])]

    model = TauPolarizationMLP(input_dim=14, hidden_dims=hidden_dims, dropout=best['dropout'])
    print(f"Arquitectura final: {hidden_dims}  dropout={best['dropout']}")
    print(f"Parámetros: {sum(p.numel() for p in model.parameters()):,}")

    train_loader = DataLoader(TauDataset(X_tr, y_tr, w_tr), batch_size=best['batch_size'], shuffle=True)
    val_loader   = DataLoader(TauDataset(X_val, y_val, w_val), batch_size=best['batch_size'], shuffle=False)
    test_loader  = DataLoader(TauDataset(X_te, y_te, w_te),  batch_size=best['batch_size'], shuffle=False)

    model, history = train(
        model, train_loader, val_loader,
        lr=best['lr'], weight_decay=best['weight_decay'],
        device=device, patience=50, n_epochs=300,
    )

    test_loss, test_acc = evaluate(model, test_loader, device)
    print(f"Test — loss: {test_loss:.4f}  acc: {test_acc:.3f}")

    roc_auc = save_evaluation_plots(model, test_loader, history, scaler, X_te, device, OUTPUT_DIR)

    # ── Guardar checkpoint ─────────────────────────────────────────────────────

    torch.save({
        'model_state_dict' : model.state_dict(),
        'hyperparams'      : {'input_dim': 14, 'hidden_dims': hidden_dims, 'dropout': best['dropout']},
        'best_params'      : best,
        'test_loss'        : test_loss,
        'test_acc'         : test_acc,
        'roc_auc'          : roc_auc,
        'history'          : history,
    }, os.path.join(OUTPUT_DIR, 'tau_polarization_mlp.pt'))

    print(f"\nTodo guardado en '{OUTPUT_DIR}/'")

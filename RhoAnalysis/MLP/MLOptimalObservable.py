import os
import torch
import torch.nn as nn
import numpy as np
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm

# ── Dataset ────────────────────────────────────────────────────────────────────

class TauDataset(Dataset):
    def __init__(self, X, y, weights):
        self.X       = torch.tensor(X,       dtype=torch.float32)
        self.y       = torch.tensor(y,       dtype=torch.float32)
        self.weights = torch.tensor(weights, dtype=torch.float32)

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx], self.weights[idx]


# ── Modelo ─────────────────────────────────────────────────────────────────────

class TauPolarizationMLP(nn.Module):
    def __init__(self, input_dim=14, hidden_dims=[128, 128, 64, 32], dropout=0.0):
        super().__init__()

        layers = []
        in_dim = input_dim
        for h in hidden_dims:
            layers += [
                nn.Linear(in_dim, h),
                nn.BatchNorm1d(h),
                nn.ReLU(),
            ]
            if dropout > 0:
                layers.append(nn.Dropout(dropout))
            in_dim = h

        layers.append(nn.Linear(in_dim, 1))
        layers.append(nn.Sigmoid())

        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


# ── Preparación de datos ───────────────────────────────────────────────────────

def prepare_data(X_p1, X_m1, w_p1, w_m1, test_size=0.30, val_size=0.50, seed=42):
    X = np.concatenate([X_p1, X_m1], axis=0)
    y = np.concatenate([np.ones(len(X_p1)), np.zeros(len(X_m1))])
    w = np.concatenate([w_p1, w_m1])

    # Balanceo de clases
    w[y == 1] *= w[y == 0].sum() / w[y == 1].sum()

    # Splits estratificados: 70 / 15 / 15
    X_tr, X_tmp, y_tr, y_tmp, w_tr, w_tmp = train_test_split(
        X, y, w, test_size=test_size, stratify=y, random_state=seed
    )
    X_val, X_te, y_val, y_te, w_val, w_te = train_test_split(
        X_tmp, y_tmp, w_tmp, test_size=val_size, stratify=y_tmp, random_state=seed
    )

    print(f"Eventos — train: {len(y_tr):,}  val: {len(y_val):,}  test: {len(y_te):,}")
    print(f"  Train  → P1: {int(y_tr.sum()):,}  M1: {int((1-y_tr).sum()):,}")
    print(f"  Val    → P1: {int(y_val.sum()):,}  M1: {int((1-y_val).sum()):,}")
    print(f"  Test   → P1: {int(y_te.sum()):,}  M1: {int((1-y_te).sum()):,}")

    # Normalización — fit SOLO en train
    scaler = StandardScaler()
    X_tr  = scaler.fit_transform(X_tr)
    X_val = scaler.transform(X_val)
    X_te  = scaler.transform(X_te)

    return (X_tr, y_tr, w_tr), (X_val, y_val, w_val), (X_te, y_te, w_te), scaler


# ── Loss ───────────────────────────────────────────────────────────────────────

def weighted_bce(output, target, weights):
    bce = nn.functional.binary_cross_entropy(output.squeeze(), target, reduction='none')
    return (bce * weights).mean()


# ── Accuracy con pesos ─────────────────────────────────────────────────────────

def weighted_accuracy(output, target, weights):
    pred = (output.squeeze() > 0.5).float()
    correct = (pred == target).float()
    return (correct * weights).sum() / weights.sum()


# ── Evaluate ───────────────────────────────────────────────────────────────────

def evaluate(model, loader, device):
    model.eval()
    total_loss, total_acc, total_w = 0.0, 0.0, 0.0
    with torch.no_grad():
        for X_batch, y_batch, w_batch in loader:
            X_batch = X_batch.to(device)
            y_batch = y_batch.to(device)
            w_batch = w_batch.to(device)
            out = model(X_batch)
            loss = weighted_bce(out, y_batch, w_batch)
            acc  = weighted_accuracy(out, y_batch, w_batch)
            n    = w_batch.sum().item()
            total_loss += loss.item() * n
            total_acc  += acc.item()  * n
            total_w    += n
    return total_loss / total_w, total_acc / total_w


# ── Training loop ──────────────────────────────────────────────────────────────

def train(model, train_loader, val_loader, n_epochs=200, lr=1e-3, weight_decay=1e-4,
          patience=20, device='cpu'):

    model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, patience=10, factor=0.5, min_lr=1e-6
    )

    best_val_loss = np.inf
    best_state    = None
    epochs_no_imp = 0
    history       = {'train_loss': [], 'val_loss': [], 'train_acc': [], 'val_acc': []}

    epoch_bar = tqdm(range(n_epochs), desc="Entrenamiento", unit="época")

    for epoch in epoch_bar:

        # ── Train ──
        model.train()
        train_loss, train_acc, train_w = 0.0, 0.0, 0.0

        batch_bar = tqdm(train_loader, desc=f"  Época {epoch+1:3d}", leave=False,
                         unit="batch", dynamic_ncols=True)

        for X_batch, y_batch, w_batch in batch_bar:
            X_batch = X_batch.to(device)
            y_batch = y_batch.to(device)
            w_batch = w_batch.to(device)

            optimizer.zero_grad()
            out  = model(X_batch)
            loss = weighted_bce(out, y_batch, w_batch)
            loss.backward()
            optimizer.step()

            acc = weighted_accuracy(out, y_batch, w_batch)
            n   = w_batch.sum().item()
            train_loss += loss.item() * n
            train_acc  += acc.item()  * n
            train_w    += n

            batch_bar.set_postfix(loss=f"{loss.item():.4f}", acc=f"{acc.item():.3f}")

        train_loss /= train_w
        train_acc  /= train_w

        # ── Validation ──
        val_loss, val_acc = evaluate(model, val_loader, device)

        scheduler.step(val_loss)
        history['train_loss'].append(train_loss)
        history['val_loss'].append(val_loss)
        history['train_acc'].append(train_acc)
        history['val_acc'].append(val_acc)

        # Actualizar barra exterior con métricas de la época
        epoch_bar.set_postfix(
            tr_loss=f"{train_loss:.4f}",
            val_loss=f"{val_loss:.4f}",
            tr_acc=f"{train_acc:.3f}",
            val_acc=f"{val_acc:.3f}",
            lr=f"{optimizer.param_groups[0]['lr']:.2e}",
            best=f"{best_val_loss:.4f}",
            no_imp=epochs_no_imp,
        )

        # ── Early stopping ──
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state    = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            epochs_no_imp = 0
        else:
            epochs_no_imp += 1
            if epochs_no_imp >= patience:
                tqdm.write(f"Early stopping en época {epoch+1} — mejor val loss: {best_val_loss:.4f}")
                break

    model.load_state_dict(best_state)
    tqdm.write(f"\nEntrenamiento finalizado — mejor val loss: {best_val_loss:.4f}")
    return model, history


# ── Carga de datos ─────────────────────────────────────────────────────────────

def load_data(p1_path, m1_path):
    """Lee ficheros npz creados por createPolDatasets.py"""
    p1_data = np.load(p1_path)["data"]
    m1_data = np.load(m1_path)["data"]

    n_p1, n_m1 = len(p1_data), len(m1_data)
    print(f"Cargados — P1: {n_p1:,} eventos  |  M1: {n_m1:,} eventos")

    # Pesos para balanceo de clases
    w_p1 = np.full(n_p1, n_m1 / n_p1)
    w_m1 = np.full(n_m1, n_p1 / n_m1)

    return p1_data, w_p1, m1_data, w_m1
# ── Plots de evaluación ────────────────────────────────────────────────────────

def save_evaluation_plots(model, test_loader, history, scaler, X_te, device, output_dir):
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve, auc

    model.eval()
    all_preds, all_labels, all_weights = [], [], []
    with torch.no_grad():
        for X_batch, y_batch, w_batch in test_loader:
            out = model(X_batch.to(device)).cpu().squeeze().numpy()
            all_preds.append(out)
            all_labels.append(y_batch.numpy())
            all_weights.append(w_batch.numpy())

    preds   = np.concatenate(all_preds)
    labels  = np.concatenate(all_labels)
    weights = np.concatenate(all_weights)

    # ── Plot 1: Curvas de aprendizaje ──────────────────────────────────────────

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    axes[0].plot(history['train_loss'], label='Train')
    axes[0].plot(history['val_loss'],   label='Val')
    axes[0].set_xlabel('Época'); axes[0].set_ylabel('BCE Loss')
    axes[0].set_title('Curvas de loss'); axes[0].legend(); axes[0].grid(True, alpha=0.3)
    axes[1].plot(history['train_acc'], label='Train')
    axes[1].plot(history['val_acc'],   label='Val')
    axes[1].set_xlabel('Época'); axes[1].set_ylabel('Accuracy')
    axes[1].set_title('Curvas de accuracy'); axes[1].legend(); axes[1].grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'learning_curves.png'), dpi=150)
    plt.close()

    # ── Plot 2: Distribución del output ───────────────────────────────────────

    fig, ax = plt.subplots(figsize=(7, 5))
    bins = np.linspace(0, 1, 50)
    ax.hist(preds[labels == 1], bins=bins, weights=weights[labels == 1],
            alpha=0.6, label=r'P1 ($\mathcal{P}_\tau=+1$)', density=True)
    ax.hist(preds[labels == 0], bins=bins, weights=weights[labels == 0],
            alpha=0.6, label=r'M1 ($\mathcal{P}_\tau=-1$)', density=True)
    ax.set_xlabel('Output de la red'); ax.set_ylabel('Densidad')
    ax.set_title('Distribución del discriminante (test set)')
    ax.legend(); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'output_distribution.png'), dpi=150)
    plt.close()

    # ── Plot 3: ROC curve ─────────────────────────────────────────────────────

    fpr, tpr, _ = roc_curve(labels, preds, sample_weight=weights)
    roc_auc = auc(fpr, tpr)
    np.savez(os.path.join(output_dir, 'roc_curve_data.npz'),
             fpr=fpr, tpr=tpr, auc=np.float64(roc_auc))
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(fpr, tpr, label=f'MLP  (AUC = {roc_auc:.4f})')
    ax.plot([0, 1], [0, 1], 'k--', label='Azar')
    ax.set_xlabel('False Positive Rate'); ax.set_ylabel('True Positive Rate')
    ax.set_title('ROC curve (test set)'); ax.legend(); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'roc_curve.png'), dpi=150)
    plt.close()

    # ── Plot 4: Output vs features ────────────────────────────────────────────

    feature_names = [
        'MesonE', 'MesonTheta', 'MesonPhi', 'MesonP',
        'gamma1_E', 'gamma1_theta', 'gamma1_phi',
        'gamma2_E', 'gamma2_theta', 'gamma2_phi',
        'lepE', 'lepTheta', 'lepPhi', 'lepP',
    ]
    fig, axes = plt.subplots(4, 4, figsize=(16, 14))
    axes = axes.flatten()
    X_te_orig = scaler.inverse_transform(X_te)
    for i, name in enumerate(feature_names):
        axes[i].scatter(X_te_orig[:, i], preds, alpha=0.05, s=1, c=labels, cmap='bwr')
        axes[i].set_xlabel(name); axes[i].set_ylabel('Output red')
        axes[i].set_title(f'Output vs {name}'); axes[i].grid(True, alpha=0.2)
    axes[-1].set_visible(False)
    plt.suptitle('Correlación output vs variables de entrada (test set)', y=1.01)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'output_vs_features.png'), dpi=150, bbox_inches='tight')
    plt.close()

    print(f"  AUC test: {roc_auc:.4f}")
    return roc_auc

# ── Main ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":

    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Usando device: {device}")

    X_p1, w_p1, X_m1, w_m1 = load_data(
        p1_path="datasets/train_data_pol1.npz",
        m1_path="datasets/train_data_pol-1.npz",
    )

    (X_tr, y_tr, w_tr), (X_val, y_val, w_val), (X_te, y_te, w_te), scaler = \
        prepare_data(X_p1, X_m1, w_p1, w_m1)

    train_loader = DataLoader(TauDataset(X_tr,  y_tr,  w_tr),  batch_size=512, shuffle=True)
    val_loader   = DataLoader(TauDataset(X_val, y_val, w_val), batch_size=512, shuffle=False)

    model = TauPolarizationMLP(input_dim=14)
    print(f"Parámetros del modelo: {sum(p.numel() for p in model.parameters()):,}")

    model, history = train(model, train_loader, val_loader, device=device, patience=50)

    # Evaluar en test set
    test_loader = DataLoader(TauDataset(X_te, y_te, w_te), batch_size=512, shuffle=False)
    test_loss, test_acc = evaluate(model, test_loader, device)
    print(f"Test — loss: {test_loss:.4f}  acc: {test_acc:.3f}")

    import os
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve, auc

    # ── Configuración de salida ────────────────────────────────────────────────────

    output_dir = "train_pol_results"  # ← cambia aquí si quieres otra carpeta
    os.makedirs(output_dir, exist_ok=True)

    # ── Guardar modelo y scaler ────────────────────────────────────────────────────

    roc_auc = save_evaluation_plots(model, test_loader, history, scaler, X_te, device, output_dir)

    print(f"\nTodo guardado en '{output_dir}/'")
    print(f"  AUC test: {roc_auc:.4f}")
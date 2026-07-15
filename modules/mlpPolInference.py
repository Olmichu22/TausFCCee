"""
Inferencia en numpy puro del MLP de polarización del tau (modelo v2), sin
dependencia de torch. Consume el .npz generado por
RhoAnalysis/MLP/exportMLPToNumpy.py.

Uso típico (en RhoHistFromTree):

    from modules import mlpPolInference
    model = mlpPolInference.load_mlp("MLPolResults/train_pol_results_optimal_v2/mlp_numpy.npz")
    x = mlpPolInference.build_features(meson, gammas, lep)   # (14,) o (N,14)
    score = model.predict(x)                                 # escalar o (N,) en [0,1]

Orden EXACTO de las 14 features (modelo por defecto = optimal_lepton, usa PIÓN):
    [PionE, PionTheta, PionPhi, PionP,
     photon1_E, photon1_theta, photon1_phi,
     photon2_E, photon2_theta, photon2_phi,
     lepE, lepTheta, lepPhi, lepP]
(El modelo v2 usaba Meson=visible; las features son posicionales, lo que cambia
es qué ramas del árbol las alimentan en _compute_nn_optimal de RhoHistFromTree.)
"""
import numpy as np

FEATURE_NAMES = [
    "PionE", "PionTheta", "PionPhi", "PionP",
    "photon1_E", "photon1_theta", "photon1_phi",
    "photon2_E", "photon2_theta", "photon2_phi",
    "lepE", "lepTheta", "lepPhi", "lepP",
]


class _MLP:
    """MLP precargado: aplica scaler + capas (Linear/BatchNorm/ReLU/Sigmoid)."""

    def __init__(self, npz_path):
        d = np.load(npz_path, allow_pickle=True)
        self.input_dim = int(d["input_dim"])
        self.feature_names = [str(s) for s in d["feature_names"]]
        self.scaler_mean = d["scaler_mean"].astype(np.float64)
        self.scaler_scale = d["scaler_scale"].astype(np.float64)
        self.layer_types = [str(s) for s in d["layer_types"]]

        # Pre-extrae los arrays por capa para no tocar el NpzFile en cada llamada.
        self.layers = []
        for i, t in enumerate(self.layer_types):
            if t == "linear":
                self.layers.append(("linear", d[f"layer{i}_W"].astype(np.float64),
                                    d[f"layer{i}_b"].astype(np.float64)))
            elif t == "batchnorm":
                self.layers.append((
                    "batchnorm",
                    d[f"layer{i}_gamma"].astype(np.float64),
                    d[f"layer{i}_beta"].astype(np.float64),
                    d[f"layer{i}_running_mean"].astype(np.float64),
                    d[f"layer{i}_running_var"].astype(np.float64),
                    float(d[f"layer{i}_eps"]),
                ))
            elif t in ("relu", "sigmoid", "identity"):
                self.layers.append((t,))
            else:
                raise ValueError(f"Tipo de capa no soportado: {t}")

        if self.input_dim != len(FEATURE_NAMES):
            raise ValueError(
                f"input_dim={self.input_dim} != {len(FEATURE_NAMES)} features esperadas")

    def predict(self, x):
        """x: array (14,) o (N,14). Devuelve escalar o (N,) en [0,1]."""
        x = np.asarray(x, dtype=np.float64)
        single = (x.ndim == 1)
        if single:
            x = x[None, :]

        # StandardScaler
        x = (x - self.scaler_mean) / self.scaler_scale

        for layer in self.layers:
            kind = layer[0]
            if kind == "linear":
                _, W, b = layer
                x = x @ W.T + b
            elif kind == "batchnorm":
                _, gamma, beta, rmean, rvar, eps = layer
                x = (x - rmean) / np.sqrt(rvar + eps) * gamma + beta
            elif kind == "relu":
                x = np.maximum(x, 0.0)
            elif kind == "sigmoid":
                x = 1.0 / (1.0 + np.exp(-x))
            elif kind == "identity":
                pass

        out = x.reshape(-1)
        return float(out[0]) if single else out


_CACHE = {}


def load_mlp(npz_path):
    """Carga (cacheada por ruta) el modelo numpy."""
    if npz_path not in _CACHE:
        _CACHE[npz_path] = _MLP(npz_path)
    return _CACHE[npz_path]


def build_features(hadron, gammas, lep):
    """Construye el vector de 14 features.

    hadron: (E, Theta, Phi, P) del hadrón del rho. Para el modelo por defecto
            (optimal_lepton) es el PIÓN cargado; para el v2 era el visible (meson).
    gammas: lista/iterable de fotones, cada uno (E, theta, phi). Se toman los
            dos de mayor E (ordenados desc.).
    lep:    (E, Theta, Phi, P) del leptón del hemisferio opuesto.
    """
    g = sorted(gammas, key=lambda p: p[0], reverse=True)
    g1, g2 = g[0], g[1]
    return np.array([
        hadron[0], hadron[1], hadron[2], hadron[3],
        g1[0], g1[1], g1[2],
        g2[0], g2[1], g2[2],
        lep[0], lep[1], lep[2], lep[3],
    ], dtype=np.float64)

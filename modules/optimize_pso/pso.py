# cutopt/pso.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Tuple, Optional
import numpy as np
from joblib import Parallel, delayed

@dataclass(frozen=True)
class PSOConfig:
    n_particles: int = 40
    n_iters: int = 80
    inertia: float = 0.7      # w
    cognitive: float = 1.4    # c1
    social: float = 1.4       # c2
    seed: int = 123
    velocity_clamp: Optional[float] = None  # ej 0.2*(range) si quieres
    patience: int = 20
    tol: float = 1e-6
    njobs: int = 8
    verbose: bool = False


@dataclass
class PSOResult:
    best_x: np.ndarray
    best_f: float
    history_best: List[float]


def pso_minimize(
    f: Callable[[np.ndarray], float],
    bounds: List[Tuple[float, float]],
    cfg: PSOConfig,
) -> PSOResult:
    rng = np.random.default_rng(cfg.seed)
    dim = len(bounds)
    lo = np.array([b[0] for b in bounds], dtype=np.float64)
    hi = np.array([b[1] for b in bounds], dtype=np.float64)
    span = hi - lo

    # init posiciones
    X = lo + rng.random((cfg.n_particles, dim)) * span
    V = (rng.random((cfg.n_particles, dim)) - 0.5) * span * 0.1

    pbest_X = X.copy()
    # pbest_f = np.array([f(x) for x in X], dtype=np.float64)
    pbest_f = np.array(
      Parallel(n_jobs=cfg.njobs, prefer="processes")(
        delayed(f)(x) for x in X
        ),
      dtype=np.float64,
      )
    gbest_idx = int(np.argmin(pbest_f))
    gbest_X = pbest_X[gbest_idx].copy()
    gbest_f = float(pbest_f[gbest_idx])

    hist = [gbest_f]
    no_improve = 0

    for _ in range(cfg.n_iters):
        r1 = rng.random((cfg.n_particles, dim))
        r2 = rng.random((cfg.n_particles, dim))

        V = (
            cfg.inertia * V
            + cfg.cognitive * r1 * (pbest_X - X)
            + cfg.social * r2 * (gbest_X - X)
        )

        if cfg.velocity_clamp is not None:
            vmax = cfg.velocity_clamp * span
            V = np.clip(V, -vmax, vmax)

        X = X + V
        X = np.clip(X, lo, hi)

        # fx = np.array([f(x) for x in X], dtype=np.float64)
        fx = np.array(
            Parallel(n_jobs=cfg.njobs, prefer="processes")(
                delayed(f)(x) for x in X
            ),
            dtype=np.float64,
        )
        improved = fx < pbest_f
        pbest_X[improved] = X[improved]
        pbest_f[improved] = fx[improved]

        new_g_idx = int(np.argmin(pbest_f))
        new_g_f = float(pbest_f[new_g_idx])

        if new_g_f + cfg.tol < gbest_f:
            gbest_f = new_g_f
            gbest_X = pbest_X[new_g_idx].copy()
            no_improve = 0
        else:
            no_improve += 1

        hist.append(gbest_f)
        if no_improve >= cfg.patience:
            break
        if cfg.verbose:
            print(
                f"[PSO] iter={len(hist):3d} "
                f"best_loss={gbest_f:.6e}"
            )


    return PSOResult(best_x=gbest_X, best_f=gbest_f, history_best=hist)

# cutopt/metrics.py
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class LossConfig:
    eps: float = 1e-9
    eff_target: float = 0.80   # “quiero conservar al menos 80%”
    eff_power: float = 2.0     # penalización suave-cuadrática
    eff_lambda: float = 5.0    # peso de la penalización
    # elegir métrica principal:
    use_s_over_b: bool = True  # si False -> S/sqrt(S+B)


def compute_score(S: float, B: float, cfg: LossConfig) -> float:
    if cfg.use_s_over_b:
        return S / (B + cfg.eps)
    return S / np.sqrt(S + B + cfg.eps)


def compute_loss(S: float, B: float, effS: float, cfg: LossConfig) -> float:
    """
    Minimizar loss.
    - Queremos score alto -> loss incluye -score
    - Penalizamos si effS < eff_target
    """
    score = compute_score(S, B, cfg)
    deficit = max(0.0, cfg.eff_target - effS)
    penalty = cfg.eff_lambda * (deficit ** cfg.eff_power)
    return -score + penalty

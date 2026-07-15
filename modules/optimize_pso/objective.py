# cutopt/objective.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple
import numpy as np
import ROOT

from .data import LoadedSample, BranchMap
from .metrics import LossConfig, compute_loss


@dataclass(frozen=True)
class CutParams:
    dR_min: float
    dR_max: float
    mesonP_min: float
    mesonP_max: float
    lepP_min: float
    lepP_max: float
    zmass_min: float
    zmass_max: float


@dataclass
class ObjectiveResult:
    loss: float
    S: float
    B: float
    effS: float
    details: Dict[str, float]


class CutObjective:
    """
    Objetivo que evalúa una configuración de cortes sobre:
      - 1 sample signal (S vs migraciones)
      - N samples background (todo B)
    """

    def __init__(
        self,
        loaded: List[LoadedSample],
        branches: BranchMap,
        selectGEN: int,
        loss_cfg: LossConfig,
    ) -> None:
        self.loaded = loaded
        self.branches = branches
        self.selectGEN = int(selectGEN)
        # if selectGEN == 2:
        #     self.selectGEN = 1
        self.loss_cfg = loss_cfg

        # Cache: total señal baseline en el fichero signal (para eficiencia)
        sig = self._get_signal_sample()
        sig_gen_mask = (sig.genTauID == self.selectGEN) & sig.base_mask
        # print(np.sum(sig_gen_mask))
        # sig_no_gen_mask = (sig.genTauID != self.selectGEN) & sig.base_mask
        # print(np.sum(sig_no_gen_mask))
        
        self._S_baseline = float(np.sum(sig.weights[sig_gen_mask]))

        if self._S_baseline <= 0:
            raise ValueError(
                "S_baseline <= 0. Revisa selectGEN o los branches/cortes fijos."
            )

    def _get_signal_sample(self) -> LoadedSample:
        sigs = [x for x in self.loaded if x.is_signal_file]
        if len(sigs) != 1:
            raise ValueError(f"Se esperaba 1 sample signal, encontrados: {len(sigs)}")
        return sigs[0]

    def _mask_for_params(self, s: LoadedSample, p: CutParams, modify_mask_inplace: bool = True) -> np.ndarray:
        arr = s.arrays

        recoMesonP = arr[self.branches.recoMesonP]
        lepP = arr[self.branches.lepP]
        dR = s.dR
        ZMass = s.ZMass

        # baseline fijo (tauPCut)
        if modify_mask_inplace:
            m = s.base_mask
        else:
            m = s.base_mask.copy()

        # cortes variables
        m &= (recoMesonP >= p.mesonP_min) & (recoMesonP <= p.mesonP_max)
        m &= (lepP >= p.lepP_min) & (lepP <= p.lepP_max)
        m &= (dR >= p.dR_min) & (dR <= p.dR_max)
        m &= (ZMass >= p.zmass_min) & (ZMass <= p.zmass_max)
        return m

    def evaluate(self, p: CutParams,  modify_mask_inplace: bool = True) -> ObjectiveResult:
        # sanity: evitar rangos invertidos
        if p.dR_min > p.dR_max or p.mesonP_min > p.mesonP_max or p.lepP_min > p.lepP_max or p.zmass_min > p.zmass_max:
            return ObjectiveResult(
                loss=1e9, S=0.0, B=1e9, effS=0.0,
                details={"reason": 1.0}
            )

        S = 0.0
        B = 0.0
        B_dict = {}
        sig = self._get_signal_sample()
        m_sig = self._mask_for_params(sig, p, modify_mask_inplace=modify_mask_inplace)

        # señal vs migración dentro del fichero signal
        m_S = m_sig & (sig.genTauID == self.selectGEN)
        m_Bmig = m_sig & (sig.genTauID != self.selectGEN)

        S += float(np.sum(sig.weights[m_S]))
        Bmig = float(np.sum(sig.weights[m_Bmig]))
        B_dict["Bmig"] = Bmig
        B += Bmig

        # fondos externos (3 ficheros bg): todo cuenta como B
        Bext = 0.0
        for bg in self.loaded:
            if bg.is_signal_file:
                continue
            m_bg = self._mask_for_params(bg, p, modify_mask_inplace=modify_mask_inplace)
            B_dict[bg.name] = float(np.sum(bg.weights[m_bg]))
            Bext += B_dict[bg.name]
        B += Bext

        effS = float(S / self._S_baseline)

        loss = float(compute_loss(S, B, effS, self.loss_cfg))
        return ObjectiveResult(
            loss=loss,
            S=S,
            B=B,
            effS=effS,
            details={
                **B_dict,
                "S_baseline": self._S_baseline,
            }
        )

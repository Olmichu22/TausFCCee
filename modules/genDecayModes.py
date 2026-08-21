"""Generator-level tau decay modes built from the tau's direct daughters.

``GenParticle.getID()`` (the ``GenTauType`` branch) is a *visible topology*
label: how many charged prongs and how many pi0 the detector could see. That is
the right quantity to compare against reco, but it merges channels that are
physically different, because an intermediate resonance is invisible to it.
The clearest example is ``tau -> K0 pi``: depending on whether the K0 shows up
as a K0_S (two extra prongs), as a K0_L (nothing) or via pi0, the very same
decay lands on GenTauType 0, 2 or 10.

This module builds the complementary label: the *true* decay mode, taken from
the PDG codes of the tau's direct daughters.

Canonicalisation rules (``canonical_daughters``):

* only direct daughters of the status-2 tau are used;
* neutrinos are dropped (the nu_tau is always there and carries no
  information, and dropping nu_e / nu_mu makes the leptonic modes collapse to
  a single ``11`` / ``13`` label);
* photons hanging directly off the tau are dropped as tau FSR — otherwise
  every radiative event would spawn its own label. That information is not
  lost: it is already exposed through ``GenTauRadPhotonMCIdx``;
* daughters of a tau+ are charge-conjugated so that tau+ and tau- decays share
  one label. Self-conjugate states (pi0, eta, omega, K0_S, K0_L…) are left
  untouched — ``-223`` is not a valid PDG code;
* the list is sorted by ``(|PDG|, PDG)`` descending, which makes the label a
  canonical form of the daughter multiset.

``MODE_TABLE`` then maps each label to a small integer, ``MODE_*``, suitable
as a histogram axis. It was seeded from a scan of 120 002 gen taus of the CLD
``ZTauTau_PolSM_March24_2M_4`` sample, which produced only 37 distinct labels
— all of them listed below. Anything unseen maps to ``MODE_UNKNOWN``.

Note that in this sample the generator writes ``tau -> pi pi0`` directly,
without an intermediate rho or a1: only 2.04 % of the taus have a genuine
resonance daughter. That 2 % is, however, exactly the fraction whose
GenTauType is ambiguous, so this is where the label pays off.
"""

# ── Canonicalisation ──────────────────────────────────────────────────────────

NEUTRINO_PDGS = {12, 14, 16}

# Partículas que son su propia antipartícula: el plegado de carga las deja
# igual, porque -223 (o -111) no son códigos PDG válidos.
SELF_CONJUGATE_PDGS = {
    21, 22, 23, 25,           # gluon, gamma, Z, H
    111, 113, 130, 221, 223,  # pi0, rho0, K0_L, eta, omega
    310, 331, 333, 335,       # K0_S, eta', phi, f2'
    441, 443, 445,            # etac, J/psi, chic2
}


def conjugate_pdg(pdg):
    """Charge-conjugate a PDG code, leaving self-conjugate states untouched."""
    pdg = int(pdg)
    return pdg if pdg in SELF_CONJUGATE_PDGS else -pdg


def canonical_daughters(tau):
    """Canonical direct-daughter PDG list of a gen tau.

    Args:
        tau (MCParticle): generator-status-2 tau.

    Returns:
        list[int]: PDG codes after dropping neutrinos and tau FSR photons,
        folding tau+ onto the tau- convention and sorting canonically.
    """
    pdgs = []
    fold = tau.getCharge() > 0
    for dau in tau.getDaughters():
        # status 0 => secundario de la simulación, no es producto real
        if dau.getGeneratorStatus() == 0:
            continue
        pdg = int(dau.getPDG())
        if abs(pdg) in NEUTRINO_PDGS:
            continue
        # Un fotón colgando directamente del tau solo puede ser FSR del tau
        if abs(pdg) == 22:
            continue
        pdgs.append(conjugate_pdg(pdg) if fold else pdg)
    pdgs.sort(key=lambda p: (abs(p), p), reverse=True)
    return pdgs


def decay_label(pdgs):
    """Comma-separated canonical label of a daughter PDG list."""
    return ",".join(str(int(p)) for p in pdgs)


# ── Mode codes ────────────────────────────────────────────────────────────────

MODE_UNKNOWN = -1

MODE_E              = 0
MODE_MU             = 1

MODE_PI             = 10
MODE_PI_PI0         = 11
MODE_PI_2PI0        = 12
MODE_PI_3PI0        = 13   # 3 o más pi0
MODE_3PI            = 14
MODE_3PI_PI0        = 15
MODE_3PI_2PI0       = 16
MODE_5PI            = 17

MODE_K              = 20
MODE_K_PI0          = 21
MODE_K_2PI0         = 22   # 2 o más pi0
MODE_K0_PI          = 23
MODE_K0_PI_PI0      = 24
MODE_K_PI_PI        = 25
MODE_K_K_PI         = 26
MODE_K_K0           = 27
MODE_K0_K0_PI       = 28
MODE_K0_3PI         = 29

MODE_OMEGA          = 40
MODE_ETA            = 41
MODE_KSTAR          = 42

MODE_NAMES = {
    MODE_UNKNOWN:   "unknown",
    MODE_E:         "e",
    MODE_MU:        "mu",
    MODE_PI:        "pi",
    MODE_PI_PI0:    "pi pi0",
    MODE_PI_2PI0:   "pi 2pi0",
    MODE_PI_3PI0:   "pi >=3pi0",
    MODE_3PI:       "3pi",
    MODE_3PI_PI0:   "3pi pi0",
    MODE_3PI_2PI0:  "3pi >=2pi0",
    MODE_5PI:       "5pi",
    MODE_K:         "K",
    MODE_K_PI0:     "K pi0",
    MODE_K_2PI0:    "K >=2pi0",
    MODE_K0_PI:     "K0 pi",
    MODE_K0_PI_PI0: "K0 pi pi0",
    MODE_K_PI_PI:   "K pi pi",
    MODE_K_K_PI:    "K K pi",
    MODE_K_K0:      "K K0",
    MODE_K0_K0_PI:  "K0 K0 pi",
    MODE_K0_3PI:    "K0 3pi",
    MODE_OMEGA:     "omega X",
    MODE_ETA:       "eta X",
    MODE_KSTAR:     "K* X",
}

# Las 37 etiquetas observadas en 120 002 taus gen de CLD ZTauTau_PolSM_March24_2M_4.
# El comentario de cada línea es la frecuencia medida en esa muestra.
MODE_TABLE = {
    # Leptónicos
    "11":                       MODE_E,             # 17.66 %
    "13":                       MODE_MU,            # 17.22 %

    # Hadrónicos sin extrañeza
    "-211":                     MODE_PI,            # 10.78 %
    "-211,111":                 MODE_PI_PI0,        # 25.55 %
    "-211,111,111":             MODE_PI_2PI0,       #  9.26 %
    "-211,111,111,111":         MODE_PI_3PI0,       #  1.03 %
    "-211,111,111,111,111":     MODE_PI_3PI0,       #  0.10 %
    "211,-211,-211":            MODE_3PI,           #  9.22 %
    "211,-211,-211,111":        MODE_3PI_PI0,       #  4.73 %
    "211,-211,-211,111,111":    MODE_3PI_2PI0,      #  0.47 %
    "211,211,-211,-211,-211":   MODE_5PI,           #  0.08 %

    # Con kaones
    "-321":                     MODE_K,             #  0.75 %
    "-321,111":                 MODE_K_PI0,         #  0.46 %
    "-321,111,111":             MODE_K_2PI0,        #  0.06 %
    "-321,111,111,111":         MODE_K_2PI0,        #  0.05 %
    "-311,-211":                MODE_K0_PI,         #  0.79 %
    "-311,-211,111":            MODE_K0_PI_PI0,     #  0.40 %
    "-311,-211,111,111":        MODE_K0_PI_PI0,     #  0.03 %
    "-321,211,-211":            MODE_K_PI_PI,       #  0.33 %
    "-321,211,-211,111":        MODE_K_PI_PI,       #  0.03 %
    "321,-321,-211":            MODE_K_K_PI,        #  0.14 %
    "321,-321,-211,111":        MODE_K_K_PI,        #  0.00 %
    "-321,311":                 MODE_K_K0,          #  0.17 %
    "-321,311,111":             MODE_K_K0,          #  0.17 %
    "310,-211,130":             MODE_K0_K0_PI,      #  0.13 %
    "-211,130,130":             MODE_K0_K0_PI,      #  0.02 %
    "311,-311,-211,111":        MODE_K0_K0_PI,      #  0.02 %
    "310,310,-211":             MODE_K0_K0_PI,      #  0.02 %
    "-311,211,-211,-211":       MODE_K0_3PI,        #  0.02 %

    # Con resonancia explícita en el registro del generador
    "-321,223":                 MODE_OMEGA,         #  0.04 %
    "223,-211,111,111":         MODE_OMEGA,         #  0.01 %
    "223,211,-211,-211":        MODE_OMEGA,         #  0.01 %
    "221,-211,111":             MODE_ETA,           #  0.14 %
    "221,211,-211,-211":        MODE_ETA,           #  0.02 %
    "-321,221":                 MODE_ETA,           #  0.01 %
    "221,-211,111,111":         MODE_ETA,           #  0.01 %
    "-323,221":                 MODE_KSTAR,         #  0.03 %
}


def true_mode(pdgs_or_label):
    """Mode code of a canonical daughter list (or of its label).

    Args:
        pdgs_or_label (list[int] | str): output of :func:`canonical_daughters`
            or of :func:`decay_label`.

    Returns:
        int: one of the ``MODE_*`` constants, ``MODE_UNKNOWN`` if unlisted.
    """
    label = (pdgs_or_label if isinstance(pdgs_or_label, str)
             else decay_label(pdgs_or_label))
    return MODE_TABLE.get(label, MODE_UNKNOWN)


def mode_name(code):
    """Human-readable name of a mode code."""
    return MODE_NAMES.get(int(code), "unknown")

import sys
import math
import ROOT
from array import array
from podio import root_io
import edm4hep

# ── Physical constants (GeV) ─────────────────────────────────────────────────
_M_TAU = 1.7769    # tau lepton mass
_M_RHO = 0.77545   # rho(770) pole mass
_M_A1  = 1.2300    # a1(1260) PDG pole mass — gives alpha_a1 ≈ 0.021
                   # (old hardcoded 0.12 corresponds to m_a1 ≈ 1.11 GeV)


# ── Internal helpers ─────────────────────────────────────────────────────────

def _alpha_V(mV, mTau=_M_TAU):
    """Alcaraz (2026) eq. (5): alpha_V = (m_tau^2 - 2*mV^2) / (m_tau^2 + 2*mV^2)."""
    return (mTau**2 - 2*mV**2) / (mTau**2 + 2*mV**2)


def _compute_ae_sm(sin2eff=0.2312):
    """Ae_SM from sin^2(theta_eff)."""
    gv_ga = 1 - 4*sin2eff
    return 2*gv_ga / (1 + gv_ga**2)


def _compute_Ptau(costheta, Atau, Ae):
    """Alcaraz (2026) eq. (3): P(z)_tau."""
    c2 = costheta**2
    denom = (1 + c2) + 2*Ae*Atau*costheta
    if abs(denom) < 1e-12:
        return 0.0
    return -(Atau*(1 + c2) + 2*Ae*costheta) / denom


def _compute_H(meson_p4, tau_p4, decay_type):
    """
    H_R(z_R) per Alcaraz (2026) eqs. (4)-(6). Returns None for leptonic channel.

    decay_type 0  (pion/kaon, pseudoscalar): H = z_R                [eq. 4]
    decay_type 1  (rho, vector):             H = alpha_V(m_rho_evt)*z_R  [eq. 5, event mass]
    decay_type 10 (a1, vector):              H = alpha_V(_M_A1)*z_R      [eq. 5, pole mass]

    z_R = (2x - 1 - xi) / (1 - xi),  x = E_R/E_tau,  xi = m_R^2/m_tau^2  [eq. 6]

    Note: for rho (type 1) the event-by-event invariant mass is used, consistent with
    wVariab in optimalVariabRho. For a1 (type 10) the broad resonance makes the pole
    mass more stable numerically.
    """
    if decay_type not in (0, 1, 10):
        return None
    xi = meson_p4.M()**2 / _M_TAU**2
    if tau_p4.E() < 1e-9 or abs(1 - xi) < 1e-9:
        return 0.0
    x = meson_p4.E() / tau_p4.E()
    z_R = max(-1.0, min(1.0, (2*x - 1 - xi) / (1 - xi)))
    if decay_type == 0:
        return z_R
    elif decay_type == 1:
        return _alpha_V(meson_p4.M()) * z_R
    else:  # 10 (a1)
        return _alpha_V(_M_A1) * z_R


def _compute_H_lep(lep_p4, beam_E):
    """Alcaraz (2026) eq. (12): H_ell(x_ell), x_ell = E_lep / E_beam.

    Factored form H_ell = (1 + x - 8x^2) / (5 + 5x - 4x^2), obtained by
    cancelling the common (1-x) factor from numerator and denominator of the
    original g(x)/f(x). Numerically stable at x->1 where the unfactored form
    is 0/0; the limit is H_ell(1) = -1.
    """
    if beam_E <= 0:
        return 0.0
    x = max(0.0, min(1.0, lep_p4.E() / beam_E))
    denom = 5 + 5*x - 4*x**2
    if abs(denom) < 1e-10:
        return 0.0
    return (1 + x - 8*x**2) / denom


def cosThetaStar(tauP4, visP4):
    """Geometric cos(theta*): polar angle of the visible system in the tau rest
    frame, relative to the tau lab flight direction.

    Exact spin-analyzing angle for the single-pion channel (alpha=1), obtained by
    boosting (no energy-fraction approximation). Equals beta_tau * z_R, where z_R
    is the analytic reconstruction used at reco (_compute_H / newAtau). Bounded to
    [-1, 1] by construction; intended for gen-level use where the true tau rest
    frame is available.
    """
    boost = tauP4.BoostVector()
    vis_rest = ROOT.TLorentzVector(visP4)
    vis_rest.Boost(-boost)
    v_vis = vis_rest.Vect()
    v_tau = tauP4.Vect()
    if v_vis.Mag() < 1e-12 or v_tau.Mag() < 1e-12:
        return 0.0
    return math.cos(v_vis.Angle(v_tau))


# ── Public single-tau weight functions ───────────────────────────────────────

def newAtau(TauP4, MesonP4, Type, New_Atau, sin_eff=None):
    """
    Single-tau hadronic reweighting per Alcaraz (2026) eqs. (4)-(6).
    Supports Type=0 (pion/kaon), Type=1 (rho), Type=10 (a1).
    Returns 1.0 for unsupported types.

    W = (1 + P_new * H_R) / (1 + P_SM * H_R)
    """
    if Type not in (0, 1, 10):
        return 1.0
    sin2eff = sin_eff if sin_eff is not None else 0.2312
    Ae = _compute_ae_sm(sin2eff)
    costheta = math.cos(TauP4.Theta())
    P_sm  = _compute_Ptau(costheta, Ae, Ae)
    P_new = _compute_Ptau(costheta, New_Atau, Ae)
    H = _compute_H(MesonP4, TauP4, Type)
    if H is None:
        return 1.0
    denom = 1 + P_sm  * H
    numer = 1 + P_new * H
    return numer / denom if abs(denom) > 1e-12 else 1.0


def newAtauLep(lepP4, lepTauP4, beamE, New_Atau, sin_eff=None):
    """
    Single-tau leptonic reweighting per Alcaraz (2026) eq. (12).

    W = (1 + P_new * H_ell) / (1 + P_SM * H_ell)

    lepP4    : 4-vector del leptón visible (e/μ) — para x_ell = E_lep/E_beam
    lepTauP4 : 4-vector del tau leptónico completo — para la dirección de polarización
    """
    if beamE <= 0:
        return 1.0
    sin2eff = sin_eff if sin_eff is not None else 0.2312
    Ae = _compute_ae_sm(sin2eff)
    costheta = math.cos(lepTauP4.Theta())
    P_sm  = _compute_Ptau(costheta, Ae, Ae)
    P_new = _compute_Ptau(costheta, New_Atau, Ae)
    H_ell = _compute_H_lep(lepP4, beamE)
    denom = 1 + P_sm  * H_ell
    numer = 1 + P_new * H_ell
    return numer / denom if abs(denom) > 1e-12 else 1.0


def newAtauFromH(TauP4, H, New_Atau, sin_eff=None):
    """
    Single-tau weight from a pre-computed spin-analyzing observable H.

    W = (1 + P_new * H) / (1 + P_SM * H)

    H is the full analyzing variable for the channel (already includes any alpha
    factor): omega for rho (wVariab/wVariabRECO), cos(theta*) for the single pion
    (cosThetaStar, alpha=1). More accurate than newAtau(Type=1) for rho because
    omega encodes the full rho -> pi pi0 structure (cosBeta, cosPsi). Suggested by
    J. Alcaraz as a cross-check.

    For rho: at gen level H=omega from wVariab (exact); at reco from wVariabRECO.
    For pion: at gen level H=cos(theta*) from cosThetaStar (exact boost); at reco
    the z_R energy-fraction form is used instead (via newAtau / _compute_H).
    """
    sin2eff = sin_eff if sin_eff is not None else 0.2312
    Ae = _compute_ae_sm(sin2eff)
    costheta = math.cos(TauP4.Theta())
    P_sm  = _compute_Ptau(costheta, Ae, Ae)
    P_new = _compute_Ptau(costheta, New_Atau, Ae)
    denom = 1 + P_sm  * H
    numer = 1 + P_new * H
    return numer / denom if abs(denom) > 1e-12 else 1.0


def newAtauRhoOmega(TauP4, omega, New_Atau, sin_eff=None):
    """Deprecated alias for newAtauFromH (H=omega). Kept for backward compatibility."""
    return newAtauFromH(TauP4, omega, New_Atau, sin_eff=sin_eff)


# ── Public two-tau joint weight functions ────────────────────────────────────

def newAtauJoint(TauP4, H, Hp, New_Atau, sin_eff=None):
    """
    Generic joint two-tau weight given pre-computed spin-analyzing values H and H'.

    W = [1 + P_new*(H + H') + H*H'] / [1 + P_sm*(H + H') + H*H']

    H  : spin-analyzing value for tau- decay (any observable: H_V, omega, H_ell...)
    Hp : spin-analyzing value for tau+ decay (same convention as newAtauJoint_had_had:
         the tau+ sign flip is already absorbed into the sumHH = H + Hp form)
    TauP4: tau- 4-vector (defines costheta for polarization).
    """
    sin2eff = sin_eff if sin_eff is not None else 0.2312
    Ae = _compute_ae_sm(sin2eff)
    costheta = math.cos(TauP4.Theta())
    P_sm  = _compute_Ptau(costheta, Ae, Ae)
    P_new = _compute_Ptau(costheta, New_Atau, Ae)
    sumHH = H + Hp
    cross = H * Hp
    denom = 1 + P_sm  * sumHH + cross
    numer = 1 + P_new * sumHH + cross
    return numer / denom if abs(denom) > 1e-12 else 1.0


def newAtauJoint_had_had(TauP4, MesonP4, OtherTauP4, OtherMesonP4,
                         Type, OtherType, New_Atau, sin_eff=None):
    """
    Joint two-tau hadronic+hadronic weight per Alcaraz (2026) eq. (9).

    W = [1 + P'*(H + H') + H*H'] / [1 + P*(H + H') + H*H']

    TauP4 (tau-) defines z=cosTheta for P(z)_tau. OtherTauP4 is tau+.
    Sign convention: H' for tau+ enters with opposite sign relative to the
    original eq. (9) formula (H - H') because the tau+ decay distribution,
    due to the antineutrino handedness, maps to -H' in the joint weight.
    Net effect: (H - H') -> (H + H') and -H*H' -> +H*H'.
    """
    if Type not in (0, 1, 10) or OtherType not in (0, 1, 10):
        return 1.0
    sin2eff = sin_eff if sin_eff is not None else 0.2312
    Ae = _compute_ae_sm(sin2eff)
    costheta = math.cos(TauP4.Theta())
    P_sm  = _compute_Ptau(costheta, Ae, Ae)
    P_new = _compute_Ptau(costheta, New_Atau, Ae)
    H  = _compute_H(MesonP4,      TauP4,      Type)
    Hp = _compute_H(OtherMesonP4, OtherTauP4, OtherType)
    if H is None or Hp is None:
        return 1.0
    # Hp enters with flipped sign for tau+ (antineutrino handedness correction)
    sumHH = H + Hp
    cross = H * Hp
    denom = 1 + P_sm  * sumHH + cross
    numer = 1 + P_new * sumHH + cross
    return numer / denom if abs(denom) > 1e-12 else 1.0


def newAtauJoint_had_lep(TauHadP4, MesonP4, TauLepP4, LepP4,
                         TypeHad, New_Atau, beamE, sin_eff=None):
    """
    Joint two-tau hadronic+leptonic weight per Alcaraz (2026) eq. (13).

    W = [1 + P'*(H_R + H_ell) + H_R*H_ell] / [1 + P*(H_R + H_ell) + H_R*H_ell]

    TauHadP4 (tau-) defines z=cosTheta for P(z)_tau. LepP4 is from tau+ side.
    H_ell for tau+ enters with flipped sign (same antineutrino handedness
    correction as in newAtauJoint_had_had): (H_R - H_ell) -> (H_R + H_ell)
    and -H_R*H_ell -> +H_R*H_ell.
    """
    if TypeHad not in (0, 1, 10):
        return 1.0
    sin2eff = sin_eff if sin_eff is not None else 0.2312
    Ae = _compute_ae_sm(sin2eff)
    costheta = math.cos(TauHadP4.Theta())
    P_sm  = _compute_Ptau(costheta, Ae, Ae)
    P_new = _compute_Ptau(costheta, New_Atau, Ae)
    H_R   = _compute_H(MesonP4, TauHadP4, TypeHad)
    H_ell = _compute_H_lep(LepP4, beamE)
    if H_R is None:
        return 1.0
    # H_ell enters with flipped sign for tau+ (antineutrino handedness correction)
    sumHH = H_R + H_ell
    cross = H_R * H_ell
    denom = 1 + P_sm  * sumHH + cross
    numer = 1 + P_new * sumHH + cross
    return numer / denom if abs(denom) > 1e-12 else 1.0


# ── Backward-compatible aliases (deprecated) ─────────────────────────────────

def newAtauRHO(TauP4, RhoP4, beamE, TauConst, Type, New_Atau, sin2theta_effective=0.2312):
    """Deprecated alias for newAtau. The beamE/TauConst args are no longer needed."""
    return newAtau(TauP4, RhoP4, Type, New_Atau, sin_eff=sin2theta_effective)


def newAtauRHO2(TauP4, RhoP4, pionP4, beamE, Type, New_Atau, sin2theta_effective=0.2312):
    """Deprecated alias for newAtau. The pionP4/beamE args are no longer needed."""
    return newAtau(TauP4, RhoP4, Type, New_Atau, sin_eff=sin2theta_effective)


# ── Deprecated originals (kept for reference) ────────────────────────────────

def newAtau_depc(TauP4, MesonP4, Type, New_Atau, sin_eff=None):

    if (Type!=0 and Type!=1 and Type!=10): # muons/electrons and others not implemented yet
         return 1
    if sin_eff is not None:
        sin2theta_effective= sin_eff
    else:
        sin2theta_effective= 0.2312
    gv_ga=  1 - 4 *sin2theta_effective
    Ae_sm=  2* gv_ga / (1+gv_ga*gv_ga)
    Atau_sm= Ae_sm

    Ae= Ae_sm

    # Polarization depends on cos(Theta):
    costheta_tau=math.cos(TauP4.Theta()) # this is the theta of the Tau, not the meson
    Ptau_sm= - (Atau_sm * (1+  costheta_tau*costheta_tau) + 2*Ae_sm*costheta_tau) / (1+costheta_tau*costheta_tau + 2*Ae_sm*Atau_sm*costheta_tau)
    Pnew = - ( New_Atau   * (1+  costheta_tau*costheta_tau) + 2*Ae*costheta_tau) / (1+costheta_tau*costheta_tau + 2*Ae* New_Atau *costheta_tau)

    alpha= 1
    if (Type==0):
        alpha= 1
    elif (Type==1):
        alpha= 0.46
    elif (Type==10):
        alpha= 0.12

    # rest frame of the tau?

    boost=ROOT.TVector3()
    boost=TauP4.BoostVector()

    meson_TauRes=ROOT.TLorentzVector()
    meson_TauRes.SetPxPyPzE(MesonP4.Px(),MesonP4.Py(),MesonP4.Pz(),MesonP4.E())
    meson_TauRes.Boost(-boost)

    # theta angle between meson and tau
    # θ = cos-1 [ (a · b) / (|a| |b|) ].
    v1 = meson_TauRes.Vect()
    v2 = TauP4.Vect()
    theta_meson= v1.Angle(v2)
    z = math.cos( theta_meson)

    weight_Pnew = (1+ alpha*Pnew*z) /(1+alpha*Ptau_sm*z)

    return weight_Pnew


def newAtauRHO_depc(TauP4, RhoP4, beamE, TauConst, Type, New_Atau, sin2theta_effective=0.2312):

    if (Type!=1): # this is for RHOs
         weight=newAtau_depc(TauP4, RhoP4,Type,New_Atau, sin_eff=sin2theta_effective)
         return weight

    gv_ga=  1 - 4 *sin2theta_effective
    Ae_sm=  2* gv_ga / (1+gv_ga*gv_ga)
    Atau_sm= Ae_sm

    Ae= Ae_sm

    # Polarization depends on cos(Theta):
    costheta_tau=math.cos(TauP4.Theta()) # this is the theta of the Tau, not the meson
    Ptau_sm= - (Atau_sm * (1+  costheta_tau*costheta_tau) + 2*Ae_sm*costheta_tau) / (1+costheta_tau*costheta_tau + 2*Ae_sm*Atau_sm*costheta_tau)
    Pnew = - ( New_Atau   * (1+  costheta_tau*costheta_tau) + 2*Ae*costheta_tau) / (1+costheta_tau*costheta_tau + 2*Ae* New_Atau *costheta_tau)

    alpha= 0.46

    pion=TauConst[0]
    pi0=TauConst[1]
    if abs(pion.getPDG())!=211:
      pion=TauConst[1]
      pi0=TauConst[0]

    pionP4=ROOT.TLorentzVector()
    pionP4.SetXYZM(pion.getMomentum().x,pion.getMomentum().y,pion.getMomentum().z,pion.getMass())

    mtau=TauP4.M()
    mRho=RhoP4.M()

    boost=ROOT.TVector3()
    boost=TauP4.BoostVector()

    Rho_TauRes=ROOT.TLorentzVector()
    Rho_TauRes.SetPxPyPzE(RhoP4.Px(),RhoP4.Py(),RhoP4.Pz(),RhoP4.E())
    Rho_TauRes.Boost(-boost)

    boostRho=ROOT.TVector3()
    boostRho=RhoP4.BoostVector()

    Pion_RhoRes=ROOT.TLorentzVector()
    Pion_RhoRes.SetPxPyPzE(pionP4.Px(),pionP4.Py(),pionP4.Pz(),pionP4.E())
    Pion_RhoRes.Boost(-boostRho)

    v1 = Rho_TauRes.Vect()
    v2 = TauP4.Vect()
    theta_Rho= v1.Angle(v2)
    z = math.cos( theta_Rho)

    v1 = Pion_RhoRes.Vect()
    v2 = RhoP4.Vect()
    beta= v1.Angle(v2)
    cosBeta = math.cos( beta)

    x=RhoP4.E()/TauP4.E()

    sqrts=beamE*2

    cosPsi= (x * (mtau*mtau + mRho*mRho) - 2*mRho*mRho)/((mtau*mtau-mRho*mRho)*math.sqrt(x*x-4*mRho*mRho/sqrts/sqrts))

    if cosPsi>1:
       cosPsi=1
    if cosPsi<-1:
       cosPsi=-1

    anglePsi=math.acos(cosPsi)

    ratioMass2= mtau*mtau/mRho/mRho

    term_a= 2/3*((1-Ptau_sm*z)- ratioMass2*(1+Ptau_sm*z)) + ratioMass2* (1+Ptau_sm*z)
    term_b= -2/3*((1-Ptau_sm*z-ratioMass2*(1+Ptau_sm*z))*(3*cosPsi*cosPsi-1)/2-3/2*math.sqrt(ratioMass2)*Ptau_sm*math.sin(2*anglePsi) *math.sin(theta_Rho))*(3*cosBeta*cosBeta-1)/2
    den = term_a+term_b

    term_a= 2/3*((1-Pnew*z)- ratioMass2*(1+Pnew*z)) + ratioMass2* (1+Pnew*z)
    term_b= -2/3*((1-Pnew*z-ratioMass2*(1+Pnew*z))*(3*cosPsi*cosPsi-1)/2-3/2*math.sqrt(ratioMass2)*Pnew*math.sin(2*anglePsi) *math.sin(theta_Rho))*(3*cosBeta*cosBeta-1)/2

    num=term_a+term_b

    weight_Pnew =num/den

    return weight_Pnew


def newAtauRHO2_depc(TauP4, RhoP4, pionP4, beamE, Type, New_Atau, sin2theta_effective=0.2312):

    if (Type!=1): # this is for RHOs
         weight=newAtau_depc(TauP4, RhoP4,Type,New_Atau, sin_eff=sin2theta_effective)
         return weight

    gv_ga=  1 - 4 *sin2theta_effective
    Ae_sm=  2* gv_ga / (1+gv_ga*gv_ga)
    Atau_sm= Ae_sm

    Ae=  Ae_sm

    # Polarization depends on cos(Theta):
    costheta_tau=math.cos(TauP4.Theta()) # this is the theta of the Tau, not the meson
    Ptau_sm= - (Atau_sm * (1+  costheta_tau*costheta_tau) + 2*Ae_sm*costheta_tau) / (1+costheta_tau*costheta_tau + 2*Ae_sm*Atau_sm*costheta_tau)
    Pnew = - ( New_Atau   * (1+  costheta_tau*costheta_tau) + 2*Ae*costheta_tau) / (1+costheta_tau*costheta_tau + 2*Ae* New_Atau *costheta_tau)

    mtau=1.7769
    mRho=RhoP4.M()

    boost=ROOT.TVector3()
    boost=TauP4.BoostVector()

    Rho_TauRes=ROOT.TLorentzVector()
    Rho_TauRes.SetPxPyPzE(RhoP4.Px(),RhoP4.Py(),RhoP4.Pz(),RhoP4.E())
    Rho_TauRes.Boost(-boost)

    boostRho=ROOT.TVector3()
    boostRho=RhoP4.BoostVector()

    Pion_RhoRes=ROOT.TLorentzVector()
    Pion_RhoRes.SetPxPyPzE(pionP4.Px(),pionP4.Py(),pionP4.Pz(),pionP4.E())
    Pion_RhoRes.Boost(-boostRho)

    v1 = Rho_TauRes.Vect()
    v2 = TauP4.Vect()
    theta_Rho= v1.Angle(v2)
    z = math.cos( theta_Rho)

    v1 = Pion_RhoRes.Vect()
    v2 = RhoP4.Vect()
    beta= v1.Angle(v2)
    cosBeta = math.cos( beta)

    x=RhoP4.E()/TauP4.E()

    sqrts=beamE*2

    cosPsi= (x * (mtau*mtau + mRho*mRho) - 2*mRho*mRho)/((mtau*mtau-mRho*mRho)*math.sqrt(x*x-4*mRho*mRho/sqrts/sqrts))

    if cosPsi>1:
       print ('What happened?', cosPsi)
       cosPsi=1
    if cosPsi<-1:
       print ('What happened?', cosPsi)
       cosPsi=-1

    anglePsi=math.acos(cosPsi)

    ratioMass2= mtau*mtau/mRho/mRho

    term_a= 2/3*((1-Ptau_sm*z)- ratioMass2*(1+Ptau_sm*z)) + ratioMass2* (1+Ptau_sm*z)
    term_b= -2/3*((1-Ptau_sm*z-ratioMass2*(1+Ptau_sm*z))*(3*cosPsi*cosPsi-1)/2-3/2*math.sqrt(ratioMass2)*Ptau_sm*math.sin(2*anglePsi) *math.sin(theta_Rho))*(3*cosBeta*cosBeta-1)/2
    den = term_a+term_b

    term_a= 2/3*((1-Pnew*z)- ratioMass2*(1+Pnew*z)) + ratioMass2* (1+Pnew*z)
    term_b= -2/3*((1-Pnew*z-ratioMass2*(1+Pnew*z))*(3*cosPsi*cosPsi-1)/2-3/2*math.sqrt(ratioMass2)*Pnew*math.sin(2*anglePsi) *math.sin(theta_Rho))*(3*cosBeta*cosBeta-1)/2

    num=term_a+term_b

    weight_Pnew =num/den

    return weight_Pnew


def newAtauLep_depc(lepP4, lepTauP4, beamE, New_Atau, sin_eff=None):
    """
    lepP4    : 4-vector del leptón visible (e/μ) — para x = E_lep/E_beam
    lepTauP4 : 4-vector del tau leptónico completo — para la dirección de polarización
    """
    if sin_eff is not None:
        sin2theta_effective = sin_eff
    else:
        sin2theta_effective = 0.2312

    if beamE <= 0:
        return 1.0

    x = lepP4.E() / beamE
    x = max(0.0, min(1.0, x))

    a = 1.0/3.0 * (5.0 - 9.0*x*x + 4.0*x*x*x)
    b = 1.0/3.0 * (1.0 - 9.0*x*x + 8.0*x*x*x)

    if abs(a) < 1e-10:
        return 1.0

    gv_ga   = 1 - 4*sin2theta_effective
    Ae_sm   = 2*gv_ga / (1 + gv_ga*gv_ga)
    Atau_sm = Ae_sm
    Ae      = Ae_sm

    costheta = math.cos(lepTauP4.Theta())  # dirección del tau leptónico completo

    Ptau_sm = -(Atau_sm*(1+costheta*costheta) + 2*Ae_sm*costheta) / \
               (1+costheta*costheta + 2*Ae_sm*Atau_sm*costheta)
    Pnew    = -(New_Atau*(1+costheta*costheta) + 2*Ae*costheta) / \
               (1+costheta*costheta + 2*Ae*New_Atau*costheta)

    return (1 + Pnew * b/a) / (1 + Ptau_sm * b/a)

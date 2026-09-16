import numpy as np
import ROOT
from modules.ParticleObjects import GenParticle, RecoParticle
import warnings
warnings.filterwarnings("once", category=UserWarning)
from modules import myutils
from modules import genDecayModes

import logging
try:
   logger = logging.getLogger("processing")
except:
   logger = None

import ROOT

 

def MatchRecoGenTau(genTau, recoTaus, nTausType, maxDRMatch=1, selectDecay=-777):
   """ Find the reconstructed tau that is closest to the generator level tau using the angle between the momenta.
      Args:
            genTau: generator level tau.
            recoTaus: list of reconstructed taus.
            maxDRMatch: maximum angle between the momenta.
      Returns:
            Tuple: Tuple with the index of the closest reconstructed tau and the number of taus of the same type.
      """
   findMatch=-1
   genVisTauP4 = genTau.getvisMomentum()
   nRecoTaus = len(recoTaus)
   
   for j in range(0,nRecoTaus):
      recoTauP4=recoTaus[j].getMomentum()
      recoTauId=recoTaus[j].getID()

      # we want to study migrations: keep all the decays but count how many are good 
      # careful, at reco level we count photons and at gen level pi0s: difference in the
      # decay mode (1 gen can be 1,2 reco)

      recoDM=recoTauId
      if recoTauId==2:
         recoDM=1
      elif (recoTauId>=11 and recoTauId<15):
         recoDM=11
      elif recoTauId>=3 and recoTauId<10:
         recoDM=3

      if selectDecay!=-777 and selectDecay==recoDM:
            nTausType+=1

      # but remove at least the leptonic ones / failed ID
      # if recoTauId<0:
      #    continue

      angleMatch=myutils.dRAngle(recoTauP4, genVisTauP4)

      # find closest
      if angleMatch<maxDRMatch:
         maxDRMatch=angleMatch
         findMatch=j
      # if logger:
      #    logger.debug(
      #       f"genTau: {genVisTauP4.P()}, recoTau: {recoTauP4.P()} id {recoDM} idx {j}, angleMatch: {angleMatch}, maxDRMatch: {maxDRMatch}, match: {findMatch}"
      #    )
         
   return findMatch, nTausType

def get_visible_final_state(particle, exclude_neutrinos=True):
    daughters = list(particle.getDaughters())
    pdg = abs(particle.getPDG())
    is_final = particle.getGeneratorStatus() == 1

    # Trata el pi0 como hoja, aunque decaiga en el generador
    if pdg == 111:
        return [particle]

    if is_final or len(daughters) == 0:
        if exclude_neutrinos and pdg in (12, 14, 16):
            return []
        return [particle]

    result = []
    for d in daughters:
        if d.getGeneratorStatus() == 0:
            continue
        result.extend(get_visible_final_state(d, exclude_neutrinos))
    return result

# ── Gen-level provenance helpers ──────────────────────────────────────────────

# Origin codes for a gen photon. Shared by GenParticle.constOrigin and by the
# GenConstOrigin / GenPhotonOrigin branches of the tau tree.
PHOTON_ORIGIN_NOT_A_PHOTON = -1
PHOTON_ORIGIN_PI0          = 0   # gamma from a pi0 decay
PHOTON_ORIGIN_TAU_FSR      = 1   # radiated by the tau itself
PHOTON_ORIGIN_CHARGED_RAD  = 2   # radiated by a charged decay product (pi, K, e, mu)
PHOTON_ORIGIN_OTHER        = 3   # any other generator-level ancestor
PHOTON_ORIGIN_SIMULATION   = 4   # secondary created by the detector simulation

# PDGs que la clasificación de visTauGen sí registra, más neutrinos y fotones
# (que llevan su propio marcado de origen). Todo neutro fuera de este conjunto
# queda sin contar en el ID del decay y por eso se etiqueta como "neutro extra".
_TAGGED_OR_COUNTED_PDGS = {22, 111, 12, 14, 16}

# Partículas cargadas de las que un fotón puede radiar dentro del decay.
_RADIATING_CHARGED_PDGS = {11, 13, 211, 321, 323}

_MAX_ANCESTOR_DEPTH = 200


def _mc_index(mcp):
   """Return the per-event object index of an MCParticle (-1 if unavailable)."""
   try:
      return int(mcp.getObjectID().index)
   except Exception:
      return -1


def _walk_ancestors(mcp, max_depth=_MAX_ANCESTOR_DEPTH):
   """Yield the ancestors of *mcp*, following the first parent at each step.

   Guards against cycles and against runaway chains, so it is safe to call on
   arbitrary MCParticle collections.
   """
   seen = set()
   cur = mcp
   for _ in range(max_depth):
      try:
         parents = list(cur.getParents())
      except Exception:
         return
      if not parents:
         return
      cur = parents[0]
      idx = _mc_index(cur)
      if idx in seen:
         return
      seen.add(idx)
      yield cur


def _is_simulation_secondary(mcp):
   """True when the particle was created by the simulation, not the generator."""
   try:
      status = int(mcp.getGeneratorStatus())
   except Exception:
      status = -1
   if status == 0:
      return True
   if status == 1:
      return False
   try:
      return bool(mcp.isCreatedInSimulation())
   except Exception:
      return False


def photon_origin_info(mcp):
   """Classify where a gen photon comes from.

   Returns:
      tuple: ``(origin_code, parent_pdg, ancestor_mc_idx)`` where *origin_code*
      is one of the PHOTON_ORIGIN_* constants, *parent_pdg* the signed PDG of
      the first non-photon ancestor and *ancestor_mc_idx* its object index
      (which groups together the two photons of the same pi0).
   """
   try:
      if abs(int(mcp.getPDG())) != 22:
         return (PHOTON_ORIGIN_NOT_A_PHOTON, 0, -1)
   except Exception:
      return (PHOTON_ORIGIN_NOT_A_PHOTON, 0, -1)

   if _is_simulation_secondary(mcp):
      return (PHOTON_ORIGIN_SIMULATION, 0, -1)

   for anc in _walk_ancestors(mcp):
      try:
         anc_pdg = int(anc.getPDG())
      except Exception:
         break
      abs_pdg = abs(anc_pdg)
      if abs_pdg == 22:
         continue          # copia del propio fotón, sigue subiendo
      idx = _mc_index(anc)
      if abs_pdg == 111:
         return (PHOTON_ORIGIN_PI0, anc_pdg, idx)
      if abs_pdg == 15:
         return (PHOTON_ORIGIN_TAU_FSR, anc_pdg, idx)
      if abs_pdg in _RADIATING_CHARGED_PDGS:
         return (PHOTON_ORIGIN_CHARGED_RAD, anc_pdg, idx)
      return (PHOTON_ORIGIN_OTHER, anc_pdg, idx)

   return (PHOTON_ORIGIN_OTHER, 0, -1)


def classify_photon_origin(mcp):
   """Return only the origin code of :func:`photon_origin_info`."""
   return photon_origin_info(mcp)[0]


def classify_tau_origin(mcp):
   """Trace where a gen tau comes from, separating primary from radiative taus.

   A tau produced through ``e+e- -> gamma*/Z -> tau tau`` also hangs from a
   photon, so the parent PDG alone cannot flag a secondary tau. The chain is
   therefore followed past the first non-tau ancestor: only if another tau shows
   up above it (``tau -> gamma -> tau tau``) is the candidate secondary.

   Returns:
      tuple: ``(origin_pdg, is_secondary, mother_tau_mc_idx, rad_photon_mc_idx)``.
      The two indices are -1 for a primary tau.
   """
   origin_pdg = 0
   is_secondary = False
   mother_idx = -1
   photon_idx = -1
   passed_non_tau = False

   for anc in _walk_ancestors(mcp):
      try:
         abs_pdg = abs(int(anc.getPDG()))
      except Exception:
         break

      if not passed_non_tau:
         if abs_pdg == 15:
            continue       # copia del propio tau
         passed_non_tau = True
         try:
            origin_pdg = int(anc.getPDG())
         except Exception:
            origin_pdg = 0
         if abs_pdg == 22:
            photon_idx = _mc_index(anc)
         continue

      if abs_pdg == 22 and photon_idx < 0:
         photon_idx = _mc_index(anc)
      if abs_pdg == 15:
         # Hay un tau por encima del primer ancestro no-tau: radiación secundaria.
         is_secondary = True
         mother_idx = _mc_index(anc)
         break

   if not is_secondary:
      photon_idx = -1

   return origin_pdg, is_secondary, mother_idx, photon_idx


def is_extra_neutral(mcp):
   """True for neutral decay products that no visTauGen counter registers.

   Covers K0_L, neutrons, Lambdas… — particles that end up in ``const`` and in
   the visible 4-momentum but leave the decay-mode ID untouched. Photons and
   pi0s are excluded on purpose: they carry their own origin tagging. K0_S is
   never seen here because it decays in the generator and its charged pions are
   counted as prongs.
   """
   try:
      if float(mcp.getCharge()) != 0.0:
         return False
      return abs(int(mcp.getPDG())) not in _TAGGED_OR_COUNTED_PDGS
   except Exception:
      return False


# Check a generator level tau candidate, find the decay,
# and compute visible (meson) variables
def visTauGen(candTau, getHelicity=False):
   """ Check a generator level tau candidate, find the decay, and compute visible (meson) variables.

   Args:
       candTau (Particle Object): The generator level tau candidate.
       getHelicity (bool): Whether to compute the helicity of the tau.
   Returns:
       Tuple: Tuple with the visible 4-momentum, the tau ID, the charge, the true 4-momentum, the maximum angle between constituents, the number of constituents, and the constituents.

   Besides the visible-topology ``ID``, the returned dict also carries the
   *true* decay mode taken from the tau's direct daughters
   (``decayDaughterPDG`` / ``trueMode``, see :mod:`modules.genDecayModes`).
   That label is purely additive: it never feeds back into ``ID``.
   """
   countPionsTauGen=0
   countPi0TauGen=0
   countMuonDecay=0
   countElectronDecay=0
   countOther=0

   genTauP4=ROOT.TLorentzVector()
   
   genTauP4.SetXYZM(candTau.getMomentum().x,candTau.getMomentum().y,candTau.getMomentum().z,candTau.getMass())

   # visible 4-momentum
   visTauP4=ROOT.TLorentzVector()
   visTauP4.SetXYZM(0,0,0,0)
   chargeTau=0
   daughters=candTau.getDaughters()
   tauID=-1

   maxAngleConsts=0
   nConsts=0
   const={}
   # Etiquetado adicional: no altera el ID del decay ni el 4-momento visible.
   constOrigin={}
   extraNeutrals={}
   nExtraNeutrals=0
   # Neutrinos: se guardan aparte, no entran en const ni en visTauP4.
   neutrinos={}
   nNeutrinos=0
   if getHelicity:
      try:
         helicity = candTau.getHelicity()
      except AttributeError:
         warnings.warn(
               "getHelicity is True but the tau candidate does not have helicity information. Setting helicity to None.",
               UserWarning,
               stacklevel=2
         )
         helicity = None
   else: 
      helicity = None

   # loop over daughter particles of the tau
   
   # IMPORTANT CHANGE -> First identify the final state products, then clasify
   # exclude_neutrinos=False: los neutrinos se recogen aparte más abajo, el
   # bucle los sigue saltando antes de tocar visTauP4 / const / maxAngle.
   final_daughters = get_visible_final_state(candTau, exclude_neutrinos=False)
   for dTau in final_daughters:
         if dTau.getGeneratorStatus() == 0:
         # Secondary, not a real product
            continue
         dauP4=ROOT.TLorentzVector()
         dauP4.SetXYZM(dTau.getMomentum().x,dTau.getMomentum().y,dTau.getMomentum().z,dTau.getMass())
         dauPDG=abs(dTau.getPDG())

         #print ('...dau',dauP4.P(),dauP4.Theta(),dauP4.Phi(),dau.getMass())

         # we want to compare the reco P4 to the 'visible' gen P4: skip neutrinos
         # PDG ID of Neutrinos

         if (dauPDG==12 or dauPDG==14 or dauPDG==16):
            # En un decay leptónico hay dos: el nu_tau y el nu_l. Guardarlos
            # es la única forma de separarlos, porque genP4-visP4 solo da la suma.
            neutrinos[nNeutrinos]=dTau
            nNeutrinos+=1
            continue 

         # lepton decays 
         if dauPDG==13 :
            countMuonDecay+=1
            #continue # either filter here or at the analysis level
         if dauPDG==11 :
            countElectronDecay+=1               
            #continue # either filter here or at the analysis level
 
         # in this Pythia sample the tau decay directly goes to pi0/pi, without the rho/a1
         # to be checked in KKMC and Whizard...
         # if there was a rho, we would need an additional step
         
         # 211 -> Pions 321 -> Kaons 323 -> K*(892) 111 -> Pi0
         if dauPDG==211 or dauPDG==321 or dauPDG==323:   # kaons and pions paired together 
            countPionsTauGen+=1
         elif dauPDG==111 :
            countPi0TauGen+=1
         # Charged particles that are not electrons or muons
         elif dTau.getCharge()!=0 and (dauPDG!=11 and dauPDG!=13):
            logger.warning(
               f"Found a charged particle with PDG {dauPDG} and charge {dTau.getCharge()} in the tau decay. "
               f"This is not expected and may indicate an issue with the tau decay reconstruction."
            )
            countOther+=1

         # compute the angle of the constituents (cone size) for further studies 
         dR=myutils.dRAngle(genTauP4,dauP4)
         if maxAngleConsts<dR:
               maxAngleConsts=dR

         # Neutros que ningún contador de arriba registra (K0_L, n, Lambda...).
         # Siguen entrando en const y en el 4-momento visible: solo se etiquetan,
         # para poder saber a posteriori que un ID=0 "en realidad" traía un K0.
         if is_extra_neutral(dTau):
            extraNeutrals[nExtraNeutrals]=dTau
            nExtraNeutrals+=1
            logger.debug(
               f"Extra neutral in the tau decay: PDG {dTau.getPDG()}, "
               f"P {dauP4.P():.3f} GeV (not counted in the decay-mode ID)."
            )

         constOrigin[nConsts]=classify_photon_origin(dTau)
         const[nConsts]=dTau
         nConsts+=1

         # Sum the visible 4-momentum and charge
         chargeTau+=dTau.getCharge()
         visTauP4+=dauP4

   # now encode the ID in a int
   # this could be much more elegant, simple for now
   if countMuonDecay>0:
      tauID=-13
   elif countElectronDecay>0:
      tauID=-11
   elif countOther>0: # refinement: check what these are 
      tauID=-2
   elif abs(chargeTau)==1: 
      if (countPionsTauGen==1):
               tauID=countPi0TauGen
      elif (countPionsTauGen==3):
               tauID=countPi0TauGen+10

   # return an object with the visible pt, ID, charge, and the true Pt 
   # a future step would be to define a class for the tau
   if tauID == 0:
      logger.debug(f"Tau Visible Momentum: {visTauP4.P()}")
      cum_momentum = ROOT.TLorentzVector()
      cum_momentum.SetXYZM(0, 0, 0, 0)
      for const_key in const:
         daup4 = ROOT.TLorentzVector()
         daup4.SetXYZM(const[const_key].getMomentum().x, const[const_key].getMomentum().y, const[const_key].getMomentum().z, const[const_key].getMass())
         cum_momentum += daup4
         logger.debug(f"Constituent {const_key} PDG {const[const_key].getPDG()}: {daup4.P()}")
         logger.debug(f"Total Visible Momentum: {cum_momentum.P()}")
   # Modo real a partir de los hijos directos: complementa al ID de topología
   # visible, no lo modifica.
   decayDaughterPDG = genDecayModes.canonical_daughters(candTau)
   trueMode         = genDecayModes.true_mode(decayDaughterPDG)
   if trueMode == genDecayModes.MODE_UNKNOWN and decayDaughterPDG:
      logger.debug(
         "Gen decay label not in MODE_TABLE: %s",
         genDecayModes.decay_label(decayDaughterPDG)
      )

   return {"visP4": visTauP4, "ID": tauID, "charge": chargeTau, "genP4": genTauP4, "maxAngleConsts": maxAngleConsts, "nConsts": nConsts, "const": const, "helicity": helicity,
           "constOrigin": constOrigin, "extraNeutrals": extraNeutrals,
           "hasExtraNeutrals": nExtraNeutrals > 0, "neutrinos": neutrinos,
           "decayDaughterPDG": decayDaughterPDG, "trueMode": trueMode}

# Reversed procedure for reconstructed pfos
# Starting from a pion, find particles in a cone around it, and 
# build the tau 

def track_momentum_error_p4_extremes(p4):
    """
    Devuelve dos TLorentzVector:
      - p4_max : con p + sigma_p
      - p4_min : con p - sigma_p
    usando la resolución FCC-ee del track momentum.
    """

    p = p4.P()
    pt = p4.Pt()

    # --- Momentum resolution ---
    # sigma_p/p = A * pt ⊕ B
    A = 0.02e-3
    B = 1e-3
    frac_sigma = np.sqrt((A * pt)**2 + B**2)

    sigma_p = frac_sigma * p

    # --- new momenta ---
    pmax = p + sigma_p
    pmin = max(p - sigma_p, 1e-8)

    # Reconstruct 4-vectors keeping same direction & mass
    eta = p4.Eta()
    phi = p4.Phi()
    m     = p4.M()

    def make_p4(newp):

        out = ROOT.TLorentzVector()
        out.SetPtEtaPhiM(newp, eta, phi, m)
        return out

    return make_p4(pmax), make_p4(pmin)
 
def electromagnetic_energy_error_p4_extremes(p4, cfg):
   """
   Devuelve dos TLorentzVector construidos explícitamente con SetPxPyPzE:
     - p4_max : con E + sigma_E
     - p4_min : con E - sigma_E
   usando la resolución FCC-ee del ECAL para fotones.
   """
   E = p4.E()

   # --- EM energy resolution ---
   # sigma_E/E = (factor / sqrt(E)) ⊕ const
   factor = cfg.get("factor", 0.)
   const = cfg.get("const", 0.01)
   stoch = factor / np.sqrt(max(E, 1e-12))

   frac_sigma = np.sqrt(stoch**2 + const**2)
   sigma_E = frac_sigma * E

   Emax = E + sigma_E
   Emin = max(E - sigma_E, 1e-8)

   # Mantener la dirección; recalcular px,py,pz a partir de theta,phi y la nueva E (p=E para fotón)
   theta = p4.Theta()
   phi   = p4.Phi()

   # pT = p / cosh(eta) where p = E (masa nula)
   def make_p4_from_E(newE):
      px = newE * np.sin(theta) * np.cos(phi)
      py = newE * np.sin(theta) * np.sin(phi)
      pz = newE * np.cos(theta)
      out = ROOT.TLorentzVector()
      out.SetPxPyPzE(px, py, pz, newE)
      return out

   return make_p4_from_E(Emax), make_p4_from_E(Emin)

 
def electromagnetic_direction_error_p4_extremes(p4, cfg):
    """
    Return one TLorentzVector constructed with SetPxPyPzE:
      - p4_smeeared :  theta smeared (along the direction) and phi random (0, 2pi)
    """
    # Get unitary vector of the original momentum
    x = p4.x
    y = p4.y
    z = p4.z
    p = np.sqrt(x**2 + y**2 + z**2)
    v0 = np.array([x/(p + 1e-8), y/(p + 1e-8), z/(p + 1e-8)])

    # random theta with std from cfg
    theta = np.random.normal(0, cfg.get("sigma_theta", 0.01))
    phi   = np.random.uniform(0, 2 * np.pi)

    # --- Rotate v0 inside a cone ---
    def get_orthogonal(v):
        if abs(v[0]) <= abs(v[1]) and abs(v[0]) <= abs(v[2]):
            ref = np.array([1., 0., 0.])
        elif abs(v[1]) <= abs(v[2]):
            ref = np.array([0., 1., 0.])
        else:
            ref = np.array([0., 0., 1.])
        orth = np.cross(v, ref)
        return orth / np.linalg.norm(orth)

    def rotate_around_axis(v, axis, angle):
        """Rodrigues' rotation formula"""
        axis = axis / np.linalg.norm(axis)
        return (v * np.cos(angle)
                + np.cross(axis, v) * np.sin(angle)
                + axis * np.dot(axis, v) * (1 - np.cos(angle)))

    orth      = get_orthogonal(v0)
    perp_axis = rotate_around_axis(orth, v0, phi)   # eje perp aleatorio
    new_dir   = rotate_around_axis(v0, perp_axis, theta)  # desviación theta

    # Rebuild p4 keeping original magnitude
    px = p * new_dir[0]
    py = p * new_dir[1]
    pz = p * new_dir[2]
    E  = p4.E()

    out = ROOT.TLorentzVector()
    out.SetPxPyPzE(px, py, pz, E)
    return out


# ---------------------------------------------------------------------------
# Correcciones extra sobre el tau ya construido (extraTauRecoCorrection)
# ---------------------------------------------------------------------------
# La reconstruccion base cuenta todo lo que cae en el cono: un foton extra
# convierte un tau->pi nu (DM 0) en un pi+gamma (DM 1). El origen de ese foton
# esta estudiado en docs/pi_extra_photon/REPORT.md: ~55 % es FSR real de la
# linea del tau y ~43 % un fragmento del shower del pion. Un corte fijo en P
# del foton no vale porque destroza pi 2pi0 / pi 3pi0, asi que la correccion
# se aplica solo a los fotones que no tienen pareja de pi0 en el cono.
#
# El registro permite anadir modos nuevos sin tocar buildTauFromPion: basta
# decorar una funcion con @registerExtraCorrection("nombre") y pedirla desde
# el YAML (seccion extra_reco_correction) o por linea de comandos.

EXTRA_CORRECTIONS = {}


def registerExtraCorrection(name):
   """Decorador para registrar un modo de correccion extra bajo `name`."""
   def _wrap(fn):
      EXTRA_CORRECTIONS[name] = fn
      return fn
   return _wrap


def availableExtraCorrections():
   """Nombres de los modos de correccion registrados."""
   return sorted(EXTRA_CORRECTIONS.keys())


def particleP4(part):
   """4-momento de un PFO/GenParticle, tolerando los dos accesores (edm4hep vs objetos propios)."""
   p4 = ROOT.TLorentzVector()
   mom = part.getMomentum()
   try:
      p4.SetXYZM(mom.x, mom.y, mom.z, part.getMass())
   except AttributeError:
      p4.SetXYZM(mom.X(), mom.Y(), mom.Z(), part.getMass())
   return p4


def assignTauID(countPions, countPhotons, countNeutrons):
   """ID del tau a partir del recuento de constituyentes (misma convencion que el tree).

   Args:
      countPions (int): numero de pi+- (kaones incluidos) en el cono.
      countPhotons (int): numero de fotones (ojo: fotones, no pi0s).
      countNeutrons (int): numero de neutrones sobre el corte.

   Returns:
      int: 0-9 para 1 prong, 10+N para 3 prong, -20 / -21 para el misID
         pion->neutron de Pandora, -1 si la combinacion no es un tau.
   """
   if countPions == 1 and countNeutrons == 0:
      # careful: here counting photons and not pi0s. Account for merged/lost photons.
      return countPhotons if countPhotons < 10 else 9

   if countPions == 3 and countNeutrons == 0:
      return countPhotons + 10  # 3 pions + photons

   if countPions == 3 and countNeutrons > 0:  # Pandora FIXME: pion -> neutron misID
      # Mismo caso que el -20 pero con 3 prongs: el neutron se suma al P4 del
      # tau, asi que estos eventos sesgaban la resolucion de la DM10 al colarse
      # como ID 10. Id propio para poder verlos como categoria aparte.
      return -21

   if countPions == 1 and countNeutrons > 0:  # Future FIXME: Pandora pion->neutron misID issue
      return -20  # To not interact with the other IDs

   return -1


def _countConstituents(const, PNeutron=0, minP_photon=0, minP_pion=0):
   """Recuenta piones / fotones / neutrones de un diccionario de constituyentes."""
   countPions = countPhotons = countNeutrons = 0
   for part in const.values():
      pdg = abs(part.getPDG())
      if pdg == 211:  # mismo criterio que buildTauFromPion (Pandora da 211 a todo hadron cargado)
         countPions += 1
      elif pdg == 22:
         countPhotons += 1
      elif pdg == 2112:
         countNeutrons += 1
   return countPions, countPhotons, countNeutrons


def _rebuildTauState(state):
   """Recalcula P4, carga, cono maximo, recuentos e ID desde `state["const"]`.

   Se usa despues de que un modo de correccion quite constituyentes. El primer
   constituyente (clave 0) sigue siendo el pion lider, que nunca se elimina.
   """
   const = {i: part for i, part in enumerate(state["const"].values())}
   lead_p4 = particleP4(const[0])

   tauP4 = ROOT.TLorentzVector(0, 0, 0, 0)
   chargeTau = 0
   maxCone = 0.0
   for part in const.values():
      p4 = particleP4(part)
      tauP4 += p4
      chargeTau += part.getCharge()
      dR = myutils.dRAngle(p4, lead_p4)
      if dR > maxCone:
         maxCone = dR

   countPions, countPhotons, countNeutrons = _countConstituents(const)

   state["const"] = const
   state["nConst"] = len(const)
   state["p4"] = tauP4
   state["charge"] = chargeTau
   state["maxCone"] = maxCone
   state["counts"] = {"pions": countPions, "photons": countPhotons, "neutrons": countNeutrons}
   state["id"] = assignTauID(countPions, countPhotons, countNeutrons)
   return state


_PI0_MASS = 0.1349768

_PION_PHOTON_FSR_DEFAULTS = {
   # Un fotón con pareja a masa de pi0 se protege: es un pi0 de verdad.
   "pi0_mass": _PI0_MASS,
   "pi0_mass_window": 0.05,
   # Fotón blando pegado al pión: fragmento del shower hadronico.
   "soft_frac": 0.05,
   # Fotón duro: FSR de la linea del tau. Un rho no pasa de m_rho, asi que
   # m(pi+gamma) alta senala que el fotón no viene de un pi0 perdido.
   "hard_p_min": 2.0,
   "mass_min": 1.2,
   # Solo 1 prong sin neutrones (DM 0/1/...): es donde vive la migracion.
   "max_photons": 0,   # 0 = sin limite
}


@registerExtraCorrection("pion_photon_fsr")
def recoverPionFromExtraPhoton(state, params=None):
   """Recupera tau->pi nu quitando fotones de FSR / shower del cono.

   Un fotón se elimina del tau si **no** tiene pareja de pi0 en el cono y ademas
   cumple una de las dos condiciones:

   - ``P_gamma / P_pion < soft_frac``: fragmento del shower hadronico del pion.
   - ``P_gamma > hard_p_min`` y ``m(pi+gamma) > mass_min``: FSR duro de la linea
     del tau (un rho no puede superar su masa, un pi0 huerfano da m(pi+gamma)
     por debajo de 1 GeV en el 90 % de los casos).

   Args:
      state (dict): estado del tau (ver :func:`extraTauRecoCorrection`).
      params (dict, optional): sobreescribe :data:`_PION_PHOTON_FSR_DEFAULTS`.

   Returns:
      dict: el estado, recalculado si se quito algun fotón.
   """
   cfg = dict(_PION_PHOTON_FSR_DEFAULTS)
   cfg.update(params or {})

   counts = state["counts"]
   if counts["pions"] != 1 or counts["neutrons"] != 0 or counts["photons"] < 1:
      return state
   if cfg["max_photons"] and counts["photons"] > cfg["max_photons"]:
      return state

   const = state["const"]
   pion_p4 = particleP4(const[0])
   photons = [(k, particleP4(p)) for k, p in const.items()
              if k != 0 and abs(p.getPDG()) == 22]
   if not photons:
      return state

   # Fotones con pareja a masa de pi0: intocables.
   paired = set()
   for i, (ki, p4i) in enumerate(photons):
      for kj, p4j in photons[i + 1:]:
         if abs((p4i + p4j).M() - cfg["pi0_mass"]) < cfg["pi0_mass_window"]:
            paired.add(ki)
            paired.add(kj)

   dropped = []
   for k, p4g in photons:
      if k in paired:
         continue
      frac = p4g.P() / pion_p4.P() if pion_p4.P() > 0 else 0.
      if frac < cfg["soft_frac"]:
         dropped.append(k)
      elif p4g.P() > cfg["hard_p_min"] and (pion_p4 + p4g).M() > cfg["mass_min"]:
         dropped.append(k)

   if not dropped:
      return state

   for k in dropped:
      del const[k]
   state["const"] = const
   state.setdefault("corrections", []).append(("pion_photon_fsr", len(dropped)))
   return _rebuildTauState(state)


def extraTauRecoCorrection(tauP4, tauID, chargeTau, maxConeTau, nConsts, const, cfg):
   """Aplica las correcciones extra configuradas sobre un tau ya construido.

   Cada modo recibe un estado con las claves ``p4``, ``id``, ``charge``,
   ``maxCone``, ``nConst``, ``const`` (dict indexado desde 0, con el pion lider
   en el 0) y ``counts`` (``pions`` / ``photons`` / ``neutrons``); devuelve el
   estado, recalculado con :func:`_rebuildTauState` si ha tocado constituyentes.

   Args:
      tauP4 (TLorentzVector), tauID (int), chargeTau (float), maxConeTau (float),
      nConsts (int), const (dict): salida de :func:`buildTauFromPion`.
      cfg (dict | list | str | None): configuracion. Un dict con ``enable``,
         ``modes`` (lista de nombres registrados) y ``params`` (dict por modo);
         una lista o un string se interpretan como la lista de modos con los
         parametros por defecto. ``None`` o vacio: no se hace nada.

   Returns:
      Tuple: ``(tauP4, tauID, chargeTau, maxConeTau, nConsts, const)`` corregidos.
   """
   if not cfg:
      return tauP4, tauID, chargeTau, maxConeTau, nConsts, const

   if isinstance(cfg, str):
      cfg = {"modes": [cfg]}
   elif isinstance(cfg, (list, tuple)):
      cfg = {"modes": list(cfg)}
   if not cfg.get("enable", True):
      return tauP4, tauID, chargeTau, maxConeTau, nConsts, const

   modes = cfg.get("modes") or []
   if isinstance(modes, str):
      modes = [modes]
   if not modes:
      return tauP4, tauID, chargeTau, maxConeTau, nConsts, const

   all_params = cfg.get("params") or {}
   countPions, countPhotons, countNeutrons = _countConstituents(const)
   state = {
      "p4": tauP4,
      "id": tauID,
      "charge": chargeTau,
      "maxCone": maxConeTau,
      "nConst": nConsts,
      "const": const,
      "counts": {"pions": countPions, "photons": countPhotons, "neutrons": countNeutrons},
      "corrections": [],
   }

   for mode in modes:
      fn = EXTRA_CORRECTIONS.get(mode)
      if fn is None:
         raise KeyError(
            f"extraTauRecoCorrection: modo '{mode}' desconocido. "
            f"Disponibles: {availableExtraCorrections()}"
         )
      state = fn(state, all_params.get(mode, {}))

   if state["corrections"] and logger is not None:
      logger.debug("extraTauRecoCorrection: %s -> ID %d", state["corrections"], state["id"])

   return (state["p4"], state["id"], state["charge"], state["maxCone"],
           state["nConst"], state["const"])


def buildTauFromPion(lead, allPfs, DRCone=1, minP_photon=0, minP_pion=0, PNeutron=1, genminP = 0.5, charge_condition=True, extra_correction=None):
   """ Starting from a pion, find particles in a cone around it, and build the tau.

   Args:
      lead (Particle Object): Pion candidate.
      allPfs (Particle Collection): All particles in the event.
      DRCone (int, optional): Radius of the cone. Defaults to 1.
      minP_photon (int, optional): Minimum photon momentum. Defaults to 0.
      minP_pion (int, optional): Minimum pion momentum. Defaults to 0.
      PNeutron (int, optional): Minimum neutron momentum. Defaults to 10.
      genminP (int, optional): Minimum general level momentum. Defaults to 0.5.
      extra_correction (dict, optional): Configuracion de las correcciones extra
         aplicadas al tau ya construido (ver extraTauRecoCorrection). None = ninguna.

   Returns:
      Tuple: Tuple with the 4-momentum of the tau, the tau ID, the charge, the maximum angle between constituents, the number of constituents, and the constituents.
   """
   countPions=1
   countPhotons=0

   # Initialize charge with the pion charge
   chargeTau=lead.getCharge()
   tauID=-1

   # Initialize the 4-momentum of the tau with the pion 4-momentum
   leadP4=ROOT.TLorentzVector()
   try:
      leadP4.SetXYZM(lead.getMomentum().x,lead.getMomentum().y,lead.getMomentum().z,lead.getMass())
   except AttributeError:
      leadP4.SetXYZM(lead.getMomentum().X(),lead.getMomentum().Y(),lead.getMomentum().Z(),lead.getMass())
   tauP4=ROOT.TLorentzVector()
   tauP4=leadP4

   maxConeTau=0
   # Constituents of the tau
   const={}

   nConsts=1
   const[0]=lead      
   countNeutrons=0

#        print ('...lead',leadP4.P(),leadP4.Theta(),math.cos(leadP4.Theta())) # leadP4.Phi(),lead.getMass())
   # Set to avoid duplicates
   # found_pions_id = set()
   for cand in allPfs:
      if (cand==lead):
         continue 

      candP4=ROOT.TLorentzVector()
      try:
         candP4.SetXYZM(cand.getMomentum().x,cand.getMomentum().y,cand.getMomentum().z,cand.getMass())
      except AttributeError:
         candP4.SetXYZM(cand.getMomentum().X(),cand.getMomentum().Y(),cand.getMomentum().Z(),cand.getMass())

      candPDG=abs(cand.getPDG())

#             print ('...cand',candP4.P(),candP4.Theta(),candP4.Phi(),cand.getMass())
      # SYSTEMATICS ERRORS
      # if candPDG == 22 and get_extremes:
      #    # Apply momentum smearing to photons if specified
      #    candP4_max, candP4_min = electromagnetic_energy_error_p4_extremes(candP4)
      # # elif candPDG == 211:
      #    # Apply momentum smearing to pions if specified
      #    candP4 = modifyMomentumWithErrors(candP4, err_d_pion, err_p_pion)
      # max angle? check backgrounds as well
      dR=myutils.dRAngle(candP4,leadP4) 
      if (dR>DRCone):  # cut to be tuned
         continue 
      # how low should we go in P?
      if (candP4.P()<genminP):
         continue 

      # now check ID and clean
      # Ignore events with electrons and muons
      if abs(candPDG)==11 or abs(candPDG)==13: 
         continue
      # Counting neutrons
      elif (candPDG==2112 and candP4.P()>PNeutron): # Pandora FIXME: pion -> neutron misID 
         countNeutrons+=1
      # Counting pions and kaons
      elif candPDG==211 and candP4.P()>minP_pion:   # ignoring the difference between kaons and pions for now 
         countPions+=1
         # found_pions_id.add(key)
      # Counting photons (should be 2 x pi0s)
      elif candPDG==22 and candP4.P()>minP_photon:  # careful: here counting photons and not pi0s. Account for merged/lost photons.  
         countPhotons+=1
      else: 
         continue

      if maxConeTau<dR:
         maxConeTau=dR

      # Sum the charge and 4-momentum of the tau
      chargeTau+=cand.getCharge()
      tauP4+=candP4
      const[nConsts]=cand
      nConsts+=1

   #print (tauP4.Pt(),chargeTau,countPions,countPhotons)

   # set the ID: only valid combinations can be a tau (charge, constituents compatible with
   # # tau decay). can be refined in the future. 
   # print("DEBUG")
   # for cn in const:
   #    print(const[cn].getPDG())
   # print("\n")
   
   if abs(chargeTau)==1 or not charge_condition:
      tauID = assignTauID(countPions, countPhotons, countNeutrons)

      # Correcciones extra configurables (p.ej. quitar el foton de FSR que
      # convierte tau->pi nu en pi+gamma). No hace nada si no se configura.
      if extra_correction:
         tauP4, tauID, chargeTau, maxConeTau, nConsts, const = extraTauRecoCorrection(
            tauP4, tauID, chargeTau, maxConeTau, nConsts, const, extra_correction)


      # return an object with P4, ID, Charge, AngleMax, nConsts, constIdx 
      # should be a class in the future
      # logger.debug(
      #    f"Tau con carga absoluta 1: "
      #    f"chargeTau: {chargeTau}, countPhotons: {countPhotons}, countPions: {countPions}, countNeutrons: {countNeutrons}, tauP4.P(): {tauP4.P()}, math.cos(tauP4.Theta()): {math.cos(tauP4.Theta())}, tauID: {tauID}"
      # )
      # print (tauP4.P(),math.cos(tauP4.Theta()),chargeTau,countPions,countPhotons,tauID)
      # if tauID == -1:
      #    logger.warning(
      #       f"Tau ID is -1, which is unexpected. "
      #       f"chargeTau: {chargeTau}, countPhotons: {countPhotons}, countPions: {countPions}, countNeutrons: {countNeutrons}"
      #    )
      return (tauP4, tauID, chargeTau, maxConeTau, nConsts, const), None

   else:
      # logger.debug(
      #    f"Tau con carga absoluta distinta de 1: "
      #    f"chargeTau: {chargeTau}, countPhotons: {countPhotons}, countPions: {countPions}, countNeutrons: {countNeutrons}"
      # )
      tauP4.SetXYZM(0,0,0,0) # safety, always return an object 
      return (tauP4,-1,0,0,nConsts,const), None

def normal_sample(mean, stddev, size=1):
   """Generate samples from a normal distribution.

   Args:
      mean (float): Mean of the distribution.
      stddev (float): Standard deviation of the distribution.
      size (int, optional): Number of samples to generate. Defaults to 1.
   Returns:
      samples (list): List of generated samples.
   """
   if not isinstance(mean, (list, np.ndarray)):
         samples = np.random.normal(mean, stddev, size)
   else:
      samP = []
      for i, m in enumerate(mean):
         samP.append(np.random.normal(m, stddev[i], size))
      samples = np.array(samP).T
   #
         
   return samples                 

# loop over all gen taus 
def findAllGenTaus(mc_particles, getHelicity=False):
   """ Find all generator level taus.
   
   Args:
       mc_particles (Particle Collection): All particles in the event.
       
   Returns:
      genTaus (dict): Dictionary with the generator level taus containing tuples with the visible 4-momentum, the tau ID, and the charge.

   Each tau also carries its provenance (``originPDG``, ``isSecondary``,
   ``motherTauKey``…) so radiative ``tau -> gamma -> tau tau`` chains can be
   told apart from the primary pair and linked back to their mother.
   """
   genTaus={}
   nGenTaus=0
   # mc_idx de cada tau y de todas sus copias -> clave en genTaus. Necesario
   # porque el ancestro que se encuentra subiendo suele ser una copia (status
   # 3/1) de la madre, mientras que en genTaus solo están las de status 2.
   tau_copy_key={}
   for particle in mc_particles:
      # only taus
      if abs(particle.getPDG()) != 15:
         continue
      # in the pythia sample we need to check the genStatus:
      # (in some events we have several copies of the tau)
      
      # 2: final state tau (to not double count)
      if particle.getGeneratorStatus()!=2:
         continue
#        print ("genTau!",particle.getGeneratorStatus())

      # tauP4=ROOT.TLorentzVector()
      # tauP4.SetXYZM(particle.getMomentum().x,particle.getMomentum().y,particle.getMomentum().z,particle.getMass())

      genTau_data=visTauGen(particle, getHelicity=getHelicity)

      origin_pdg, is_secondary, mother_mc_idx, rad_photon_mc_idx = classify_tau_origin(particle)
      genTau_data["mcIdx"]          = _mc_index(particle)
      genTau_data["originPDG"]      = origin_pdg
      genTau_data["isSecondary"]    = is_secondary
      genTau_data["motherTauMCIdx"] = mother_mc_idx
      genTau_data["radPhotonMCIdx"] = rad_photon_mc_idx

      genTau = GenParticle(**genTau_data)

      if genTau.getCharge()<0:
         genTau.setPDG(15)
      else:
         genTau.setPDG(-15)
      # visTauP4=genTau[0]
      # genTauId=genTau[1]

      # Registra el tau y su cadena de copias para poder resolver la madre.
      tau_copy_key[_mc_index(particle)] = nGenTaus
      for anc in _walk_ancestors(particle):
         if abs(int(anc.getPDG())) != 15:
            break
         tau_copy_key[_mc_index(anc)] = nGenTaus

      genTaus[nGenTaus]=genTau
      nGenTaus+=1

   # Se resuelve al final: la madre puede aparecer después que la hija.
   for genTau in genTaus.values():
      mother_mc_idx = genTau.getMotherTauMCIdx()
      if mother_mc_idx >= 0:
         genTau.setMotherTauKey(tau_copy_key.get(mother_mc_idx, -1))

   return genTaus

def modifyMomentumWithErrors(p4, err_d: float = 0., err_p: float = 0.):
      if err_d>0.:
         theta = p4.Theta()
         phi = p4.Phi()
         p = p4.P()
         m = p4.M()
         new_values =  normal_sample([theta, phi], [err_d, err_d], 1)[0]
         new_theta = new_values[0]
         new_phi = new_values[1]
         if new_theta < 0:
            new_theta = abs(new_theta)
         if new_theta > np.pi:
            new_theta = 2 * np.pi - new_theta
         if new_phi < 0:
            new_phi += 2 * np.pi
         if new_phi > 2 * np.pi:
            new_phi -= 2 * np.pi
         
         new_direction = np.array([np.sin(new_theta) * np.cos(new_phi),
                                    np.sin(new_theta) * np.sin(new_phi),
                                    np.cos(new_theta)])
         new_Pxyz = new_direction * p
         p4.SetPxPyPzE(new_Pxyz[0], new_Pxyz[1], new_Pxyz[2], np.sqrt(p**2 + m**2))
      if err_p>0.:
         theta = p4.Theta()
         phi = p4.Phi()
         p = p4.P()
         m = p4.M()
         new_p = normal_sample([p], [err_p], 1)[0]
         new_Pxyz = np.array([np.sin(theta) * np.cos(phi),
                              np.sin(theta) * np.sin(phi),
                              np.cos(theta)]) * new_p
         p4.SetPxPyPzE(new_Pxyz[0], new_Pxyz[1], new_Pxyz[2], np.sqrt(new_p**2 + m**2))
      return p4
      
# function to find all reco taus starting from PFO collection 
def findAllTaus(pfos,
                dRMax,
                minP_photon,
                minP_pion,
                PNeutron,
                genminP,
                charge_condition=True,
                extra_correction=None,):
   """ Find all tau candidates starting from PFO collection by recognizing the decay products.

   Args:
      pfos (PandoraPFOs): PandoraPFOs collection
      dRMax (float): Maximum cone radius.
      minP_photon (float): Minimum photon momentum.
      minP_pion (float): Minimum pion momentum.
      PNeutron (float): Minimum neutron momentum.
      genminP (float): Minimum general level momentum.
      extra_correction (dict, optional): Configuracion de las correcciones extra
         (ver extraTauRecoCorrection). None = reconstruccion base sin tocar.

   Returns:
       taus (dict): Dictionary with the tau candidates containing tuples with the visible 4-momentum, the tau ID, and the charge.
   """
   taus={}
   nTaus=0
   # Dict to avoid duplicates
   # id_pfos = {i: pfos[i] for i in range(len(pfos))}
   found_pions_id = set()
   # for key, pf in id_pfos.items():
   for pf in pfos:
      if (abs(pf.getPDG())!=211): 
         continue 
      # if key in found_pions_id:
      #    continue
      pionP4 = ROOT.TLorentzVector()
      try:
         pionP4.SetXYZM(pf.getMomentum().x,pf.getMomentum().y,pf.getMomentum().z,pf.getMass())
      except AttributeError as e:
         pionP4.SetXYZM(pf.getMomentum().X(),pf.getMomentum().Y(),pf.getMomentum().Z(),pf.getMass())
      
      # SYSTEMATICS ERRORS
      # pionP4 = modifyMomentumWithErrors(pionP4, err_d_pion, err_p_pion)
      
      if pionP4.P() < minP_pion or  pionP4.P() < genminP:
         continue

      recoTau_data, pions_id = buildTauFromPion(pf, pfos, dRMax, minP_photon, minP_pion, PNeutron, genminP, charge_condition, extra_correction)
      recoTau = RecoParticle(recoTau_data[0], recoTau_data[1], recoTau_data[2], recoTau_data[3], recoTau_data[4], recoTau_data[5])
      # logger.debug(
      #    f"Id del RecoTau {recoTau.getID()}"
      # )
      if recoTau.getCharge()<0:
         recoTau.setPDG(15)
      else:
         recoTau.setPDG(-15)
      
      # # Add the main pion to pions_id
      # if pions_id:
      #    pions_id.add(key)
      #    found_pions_id.update(pions_id)
      # else:
      #    found_pions_id.add(key)
      
      
      candTauP4=recoTau.getMomentum()
      # candTauId=recoTau[1]
      # candTauCharge=recoTau[2]

      # # FIXME: this is very ugly, angular separation between taus to avoid duplicates
      # # in these samples most of the events have 1 tau (and 1 prong)
      # # in a real scenario this could be slow and we could veto events 
      duplicate=False
      for i in range(0,nTaus):
         if (myutils.dRAngle(candTauP4,taus[i].getMomentum())<0.05): duplicate=True
      if (duplicate==True): continue

      taus[nTaus]=recoTau
      nTaus+=1
      #print ("...",pf.getObjectID().index,candTauP4.Pt(),candTauP4.Phi(),candTauP4.Theta(),candTauId,candTauCharge)

   return taus


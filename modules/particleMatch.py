import math
import ROOT
from modules.ParticleObjects import GenParticle, RecoParticle, GenRecoMatched
from typing import Dict, List
from modules import myutils
from dataclasses import dataclass, field
import logging

try:
    logger = logging.getLogger("processing")
except:
    logger = None

@dataclass
class ParticleMatchResults:
    """Container for the results of particle matching.
    Attributes:
          unmatched_gen (Dict[int, GenParticle]): Dictionary of unmatched generator particles.
          unmatched_reco (Dict[int, RecoParticle]): Dictionary of unmatched reconstructed particles.
          gen_matched_with_reco (List[GenRecoMatched]): List of generator particles matched with reconstructed particles.
          reco_matched_with_other (List[GenRecoMatched]): List of reconstructed particles matched with other types of particles.
          gen_matched_with_other (List[GenRecoMatched]): List of generator particles matched with other types of particles.
          all_gen (Dict[int, GenParticle]): Dictionary of all generator particles considered in the matching.
          all_reco (Dict[int, RecoParticle]): Dictionary of all reconstructed particles considered in the matching.
          all_matched (List[GenRecoMatched]): List of all matched particles regardless of type.
    """
    unmatched_gen: Dict[int, GenParticle] = field(default_factory=dict)
    unmatched_reco: Dict[int, RecoParticle] = field(default_factory=dict)
    gen_matched_with_reco: List[GenRecoMatched] = field(default_factory=list)
    reco_matched_with_other: List[GenRecoMatched] = field(default_factory=list)
    gen_matched_with_other: List[GenRecoMatched] = field(default_factory=list)
    all_gen: Dict[int, GenParticle] = field(default_factory=dict)
    all_reco: Dict[int, RecoParticle] = field(default_factory=dict)
    all_matched: List[GenRecoMatched] = field(default_factory=list)
# ----------------------------------------------------------------------

def GetRecoPions(pfos):
    """Find all reco level pion.

    Args:
        mc_particles (Particle Collection): All particles in the event.

    Returns:
       genPions (dict): Dictionary with the generator level pions containing tuples with the visible 4-momentum, the tau ID, and the charge.
    """
    recoPions = {}
    nRecoPions = 0
    for particle in pfos:
          # only pions
        if abs(particle.getPDG()) != 211:
            continue

        try:
            dauP4 = ROOT.TLorentzVector()
            dauP4.SetXYZM(
                particle.getMomentum().x,
                particle.getMomentum().y,
                particle.getMomentum().z,
                particle.getMass(),
            )
            recoPion = RecoParticle(
                dauP4, 211, particle.getCharge(), 0, 0, None, 211, nRecoPions)
        except AttributeError:
            recoPion = particle
            recoPion.setIdx(nRecoPions)

        recoPions[nRecoPions] = recoPion
        nRecoPions += 1

    return recoPions


# ----------------------------------------------------------------------
# Utilidad: recorrer recursivamente las hijas de un tau y recoger productos finales
# ----------------------------------------------------------------------
def GetGenTauDecayProducts(mc_particles, only_final_state=True):
    """
    Recorre todos los taus generadores (status==2) y obtiene sus productos de decaimiento.
    Si only_final_state=True, devuelve solo partículas finales (status==1).

    Returns:
        genProds (dict): idx -> GenParticle (de tu clase), conteniendo 4-momento y PDG de cada producto.
    """
    genProds = {}
    n = 0
    for particle in mc_particles:
        # Solo taus generadores finales (para no contar duplicados)
        if abs(particle.getPDG()) != 15:
            continue
        if particle.getGeneratorStatus() != 2:
            continue

        stack = list(particle.getDaughters())
        while stack:
            d = stack.pop()
            d_status = getattr(d, "getGeneratorStatus", lambda: None)()
            d_pdg = abs(d.getPDG())

            # Si queremos solo productos finales (status==1), seguimos bajando si no lo son
            if only_final_state and d_status is not None and d_status != 1:
                # seguir explorando descendencia
                for dd in d.getDaughters():
                    stack.append(dd)
                continue

            dauP4 = ROOT.TLorentzVector()
            dauP4.SetXYZM(
                d.getMomentum().x, d.getMomentum().y, d.getMomentum().z, d.getMass()
            )
            genPart = GenParticle(
                dauP4, d_pdg, d.getCharge(), dauP4, 0, 0, None, d_pdg, d, n
            )
            genProds[n] = genPart
            n += 1

    return genProds


def GetGenPions(mc_particles, generator_status=True):
    """Find all generator level pion.

    Args:
        mc_particles (Particle Collection): All particles in the event.

    Returns:
       genPions (dict): Dictionary with the generator level pions containing tuples with the visible 4-momentum, the tau ID, and the charge.
    """
    genPions = {}
    nGenPion = 0

    for mc in mc_particles:
        mcPDG = abs(mc.getPDG())
        if mcPDG != 211:
            continue  # solo piones cargados
        if mc.getGeneratorStatus() != 1 and generator_status:
            continue  # solo estado final
        mcP4 = ROOT.TLorentzVector()
        mcP4.SetXYZM(
            mc.getMomentum().x, mc.getMomentum().y, mc.getMomentum().z, mc.getMass()
        )

        # # Corte angular: aceptancia del detector
        # costheta = math.cos(mcP4.Theta())
        if abs(math.cos(mcP4.Theta())) > 0.95:
            continue  # Cambio para código María

        genPion = GenParticle(mcP4, mcPDG, mc.getCharge(),
                              mcP4, 0, 0, None, mcPDG, mc, nGenPion)
        genPions[nGenPion] = genPion
        nGenPion += 1

    return genPions


def GetRecoParticles(pfos):
    recoParticles = {}
    nRecoParticles = 0
    for particle in pfos:
        try:
            dauP4 = ROOT.TLorentzVector()
            dauP4.SetXYZM(
                particle.getMomentum().x,
                particle.getMomentum().y,
                particle.getMomentum().z,
                particle.getMass(),
            )
            recoParticle = RecoParticle(
                dauP4,
                nRecoParticles,
                particle.getCharge(),
                0,
                0,
                None,
                particle.getPDG(),
                nRecoParticles
            )
        except AttributeError:
            recoParticle = particle
            recoParticle.setIdx(nRecoParticles)
        recoParticles[nRecoParticles] = recoParticle
        nRecoParticles += 1
    return recoParticles


def MatchGenWithRecoParticles(
    genDaus, recoParticles, maxDRMatch=1, non_considered_particles=[], force=False
):
    """
    Matches generator-level (truth) particles with reconstructed particles 
    based on angular separation (ΔR).  

    For each generated daughter particle, the function finds the closest 
    reconstructed particle within a given ΔR cone. Certain PDG types can 
    be excluded from the matching. Returns a list of matched Gen-Reco pairs.

    Args:
        genDaus (dict): Dictionary of generated (truth-level) daughter particles, 
            where each value provides access to methods like `getMomentum()`, 
            `getPDG()`, and `getIdx()`.
        recoParticles (dict): Dictionary of reconstructed particles, with similar 
            interface to `genDaus`. Each entry must support `getMomentum()` and `getPDG()`.
        maxDRMatch (float, optional): Maximum ΔR distance allowed for a valid match. 
            Defaults to 1.
        non_considered_particles (list[int], optional): List of PDG IDs to ignore 
            during the matching. Defaults to an empty list.
        force (bool, optional): If True, may override certain matching constraints 
            (currently unused). Defaults to False.

    Returns:
        list[GenRecoMatched]: List of `GenRecoMatched` objects, each containing 
            the generator particle, its index, and the corresponding reconstructed 
            particles within the ΔR cone.

    Notes:
        - Skips generator particles whose PDG ID is in `non_considered_particles`.
        - Avoids matching photons (PDG 22) with charged pions (PDG 211) and vice versa.
        - Uses `myutils.dRAngle()` to compute ΔR between generated and reconstructed momenta.
        - The matching stores not only the best match but also all reconstructed 
          particles found within the ΔR cone.
    """
    gen_reco_matches = []

    for genDau in genDaus:
        findMatch = -1
        genP4 = genDaus[genDau].getMomentum()

        minDR = maxDRMatch
        gendauPDG = abs(genDaus[genDau].getPDG())
        minDRs = []
        particles_in_cone = []
        id_particles_in_cone = []
        # if gendauPDG != 211:
        #    continue
        # ignorar partículas no consideradas
        if gendauPDG in non_considered_particles:
            continue

        for recoParticle in recoParticles:
            #  if recoParticle in reco_gen_match.values():  # quitar condición
            #        # ya está enlazado, no lo volvemos a usar
            #        continue

            recoPartP4 = recoParticles[recoParticle].getMomentum()
            angleMatch = myutils.dRAngle(recoPartP4, genP4)
            recoPDG = abs(recoParticles[recoParticle].getPDG())

            if gendauPDG == 211 and recoPDG == 22:
                continue
            if gendauPDG == 22 and recoPDG == 211:
                continue

            if angleMatch < minDR:
                minDR = angleMatch
                findMatch = recoParticle
                minDRs.append(angleMatch)
                particles_in_cone.append(recoParticles[recoParticle])
                id_particles_in_cone.append(recoParticle)

        # if minDR > 0.1: # CAMBIO CODIGO MARIA
        #     continue
        if findMatch != -1:
            # Si no hay match con pión, usa la partícula más cercana
            matched_obj = GenRecoMatched(
                genDaus[genDau].getIdx(), genDaus[genDau], id_particles_in_cone, particles_in_cone, minDRs)
            gen_reco_matches.append(matched_obj)

    return gen_reco_matches


def MatchedUnmatchedParticles(
    mc_particles, pfos, maxDRMatch=1, non_considered_particles=[], force=False
):
    """Get matched and unmatched particles from generator and reconstructed particles.
    Args:
        mc_particles (Particle Collection): All particles in the event.
        pfos (Particle Collection): All particles in the event.
        maxDRMatch (float, optional): Maximum angle between the momenta. Defaults to 1.
    Returns:
        Tuple: Tuple with dictionaries of matched and unmatched pions and photons.
    """
    # GetGenDaughters(mc_particles) GetGenPions
    genDaus = GetGenPions(mc_particles)
    # genDaus = GetGenTauDecayProducts(mc_particles, only_final_state=True)
    # recoParticles = GetRecoPions(pfos)
    recoParticles = GetRecoParticles(pfos)

    reco_gen_match = MatchGenWithRecoParticles(
        genDaus, recoParticles, maxDRMatch, non_considered_particles, force
    )
    matching_results: dict[str, ParticleMatchResults] = {}
    for pdg_id in [211, 22]:
        matching_result = get_particle_matching_results(
            genDaus, recoParticles, reco_gen_match, pdg_id)
        matching_results["all_matched"] = reco_gen_match
        matching_results[str(pdg_id)] = ParticleMatchResults(**matching_result)

    return matching_results


def get_particle_matching_results(genDaus, recoParticles, reco_gen_match, pdg_code):
    """
    Generalizes the matching process between generated (gen) and reconstructed (reco) particles.

    Args:
        genDaus (dict): Dictionary containing the generated particles.
        recoParticles (dict): Dictionary containing the reconstructed particles.
        reco_gen_match (list): List of matching objects between gen and reco.
        pdg_code (int): PDG ID of the particle to filter (e.g., 211 for pions, 22 for photons).

    Returns:
        dict: Dictionary containing all matching results, organized by particle type.
    """

    # Diccionarios de partículas no emparejadas y totales
    unmatched_gen = dict()
    unmatched_reco = dict()
    all_gen = dict()
    all_reco = dict()

    # --- Partículas generadas sin match ---
    for genDau in genDaus:
        if genDaus[genDau] not in reco_gen_match and abs(genDaus[genDau].getPDG()) == abs(pdg_code):
            unmatched_gen[genDau] = genDaus[genDau]
            all_gen[genDau] = genDaus[genDau]

    # --- Partículas reconstruidas sin match ---
    for recoParticle in recoParticles:
        if (
            recoParticles[recoParticle] not in reco_gen_match
            and abs(recoParticles[recoParticle].getPDG()) == abs(pdg_code)
        ):
            unmatched_reco[recoParticle] = recoParticles[recoParticle]
            all_reco[recoParticle] = recoParticles[recoParticle]

    # --- Listas de coincidencias ---
    gen_matched_with_reco = []
    gen_matched_with_other = []
    reco_matched_with_other = []

    # --- Procesar coincidencias ---
    for matched_obj in reco_gen_match:
        recoPDG = abs(matched_obj.getMatchedRecoPDG())
        genPDG = abs(matched_obj.getGenPDG())
        genID = matched_obj.getGenID()
        genParticle = matched_obj.getGenParticle()
        recoID = matched_obj.getMatchedRecoID()
        recoParticle = matched_obj.getMatchedRecoParticle()

        if recoPDG == abs(pdg_code) and genPDG == abs(pdg_code):
            gen_matched_with_reco.append(matched_obj)
            all_gen[genID] = genParticle
            all_reco[recoID] = recoParticle
        elif recoPDG == abs(pdg_code):
            reco_matched_with_other.append(matched_obj)
            all_reco[recoID] = recoParticle
        elif genPDG == abs(pdg_code):
            gen_matched_with_other.append(matched_obj)
            all_gen[genID] = genParticle

    # --- Resultado final ---
    results = {
        f"unmatched_gen": unmatched_gen,
        f"unmatched_reco": unmatched_reco,
        f"gen_matched_with_reco": gen_matched_with_reco,
        f"reco_matched_with_other": reco_matched_with_other,
        f"gen_matched_with_other": gen_matched_with_other,
        f"all_gen": all_gen,
        f"all_reco": all_reco,
    }

    return results

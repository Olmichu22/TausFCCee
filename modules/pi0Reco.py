import ROOT
from modules.ParticleObjects import RecoParticle
from modules import myutils
import itertools

import logging
logger = logging.getLogger("pi0mass")

PI0INVARIANTMASS = 0.1349768 # GeV

def getPi0Mass(photons: dict[int,RecoParticle], strategy):
  """
  Calculates the invariant mass of a neutral pion candidate from photon pairs.

  Depending on the number of input photons and the selected strategy, the function
  finds the best photon pair based on either their invariant mass proximity to the
  π⁰ nominal mass or their angular distance (ΔR).  

  If exactly two photons are provided, the function directly computes their invariant mass.
  For more than two photons, it searches over all possible photon pairs and selects the best
  candidate according to the specified strategy.

  Args:
      photons (dict[int, RecoParticle]): Collection of photon objects.  
      strategy (dict): Dictionary defining the pairing strategy.  
          Supported keys:
            - `"mass"` (*float*): Selects the pair whose invariant mass is closest
              to the π⁰ invariant mass (`PI0INVARIANTMASS`). Optional threshold can
              be set with the `"mass"` value.
            - `"distance"` (*float*): Selects the pair with the smallest angular
              distance (ΔR) between photons. A maximum allowed distance can be set
              with the `"distance"` value.

  Returns:
      tuple:
          - **mass** (*float or None*): Calculated invariant mass of the selected photon pair,  
            or `None` if no valid pair is found.
          - **non_matched_photons** (*dict or list*): Remaining photons not used in the best pair.  
            If no valid pair is found, returns the original photon collection.

  Raises:
      AttributeError: If a photon object does not provide compatible `getMomentum()` attributes
          (`x`, `y`, `z` vs. `X()`, `Y()`, `Z()`).

  Notes:
      - For fewer than two photons, the function returns `(None, photons)`.
      - If `strategy["mass"]` is provided, the best pair is chosen as the one
        whose invariant mass is closest to `PI0INVARIANTMASS`.
      - If `strategy["distance"]` is provided, the best pair is chosen based on
        minimal ΔR distance between photon momenta (computed via `myutils.dRAngle()`).
      - Uses `ROOT.TLorentzVector` to construct 4-momenta for invariant mass calculation.
      - Logs detailed information and warnings using the global `logger` instance.
  """
  # Conditions to not calculate the mass
  if len(photons) < 2:
    # Return non matched photon
    return None, photons
  elif len(photons)==2:
    # Only two photons, return mass
    P1 = ROOT.TLorentzVector()
    P2 = ROOT.TLorentzVector()
    try:
      P1.SetXYZM(
                  photons[0].getMomentum().x,
                  photons[0].getMomentum().y,
                  photons[0].getMomentum().z,
                  photons[0].getMass(),
              )
      P2.SetXYZM(
                  photons[1].getMomentum().x,
                  photons[1].getMomentum().y,
                  photons[1].getMomentum().z,
                  photons[1].getMass(),
              )
    except AttributeError:
      P1.SetXYZM(
                photons[0].getMomentum().X(),
                photons[0].getMomentum().Y(),
                photons[0].getMomentum().Z(),
                photons[0].getMass(),
              )
      P2.SetXYZM(
                photons[1].getMomentum().X(),
                photons[1].getMomentum().Y(),
                photons[1].getMomentum().Z(),
                photons[1].getMass(),
              )
    return cumulatedPhotonsMass(P1, P2), None
  
  # Posible pairs combinations
  combinations = itertools.combinations(photons.keys(), 2)

  photon_momentums = {}
  for i in photons.keys():
    P = ROOT.TLorentzVector()
    try:
      P.SetXYZM(
                  photons[i].getMomentum().x,
                  photons[i].getMomentum().y,
                  photons[i].getMomentum().z,
                  photons[i].getMass(),
              )
    except AttributeError:
      P.SetXYZM(
                photons[i].getMomentum().X(),
                photons[i].getMomentum().Y(),
                photons[i].getMomentum().Z(),
                photons[i].getMass(),
              )
    photon_momentums[i] = P
  logger.debug(f"Trying to find best pair of photons with strategy {strategy} {list(strategy.keys())[0]}")
  # Two strategies: mass and distance
  if list(strategy.keys())[0] == "mass":
    # Minimum acceptable mass
    min_mass = strategy.get("mass", -1)
    
    best_mass = 99999999
    best_pair = None
    for i,j in list(combinations):
      mass = cumulatedPhotonsMass(photon_momentums[i], photon_momentums[j])
      
      logger.debug(f"Mass of pair {i} and {j} is {mass}")
      
      # Evaluate new pair mass
      if abs(best_mass-PI0INVARIANTMASS) > abs(mass-PI0INVARIANTMASS):
        
        logger.debug(f"Changing best mass {best_mass} to {mass}")
        
        best_mass = mass
        best_pair = [i,j]
        
    
    # Evaluate result    
    if best_pair is not None and best_mass > min_mass:
      # Get non matched photons
      non_matched_photons = [k for k in photons.keys() if k not in best_pair]
      non_matched_photons = {k:photons[k] for k in non_matched_photons}
      logger.debug(f"Non matched photons are {non_matched_photons.keys()}")
      return best_mass, non_matched_photons
    
    else:
      logger.warning(f"No best pair found. Best pair is {best_pair} with mass {best_mass}")
      return None, photons
        
  elif list(strategy.keys())[0] == "distance":
    if strategy["distance"] != -1:
      max_distance = strategy["distance"]
    else:
      max_distance = 99999999
    
    min_distance_pair = None
    min_distance_value = 99999999
    for i,j in list(combinations):
      distance = myutils.dRAngle(photon_momentums[i], photon_momentums[j])
      
      logger.debug(f"Distance of pair {i} and {j} is {distance}")
      
      # Evaluate new pair mass
      if abs(min_distance_value) > abs(distance):
        
        logger.debug(f"Changing best distance {min_distance_value} to {distance}")
        
        min_distance_value = distance
        min_distance_pair = [i,j]
      
    if min_distance_pair is not None and min_distance_value < max_distance:
      logger.debug(f"Best distance pair is {min_distance_pair}")
      
      non_matched_photons = [k for k in photons.keys() if k not in min_distance_pair]
      non_matched_photons = {k:photons[k] for k in non_matched_photons}

      logger.debug(f"Non matched photons are {non_matched_photons}")
      
      mass = cumulatedPhotonsMass(photon_momentums[min_distance_pair[0]], photon_momentums[min_distance_pair[1]])
      
      return mass, non_matched_photons
    else:
      logger.warning(f"No best pair found. Best pair is {min_distance_pair} with distance {min_distance_value}")
      return None, photons
    
  else:
    logger.warning("Unknown strategy for pi0 mass calculation")
    return None, None
  

def cumulatedPhotonsMass(P1, P2):
  
  P1 += P2
  return P1.M()
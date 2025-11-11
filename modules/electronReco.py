import ROOT
from modules.ParticleObjects import RecoParticle

def findAllElectrons(pfos, minPt):
  """ Find all tau candidates starting from PFO collection by recognizing the decay products.

  Args:
      pfos (PandoraPFOs): PandoraPFOs collection
      minPt (float): Minimum particle momentum.
  Returns:
      electrons (dict): Dictionary with the electron candidates containing RecoParticle Objects.
  """
  electrons={}
  nElectrons=0
  for pf in pfos:
    if (abs(pf.getPDG())!=11): 
        continue 

    
    electronP4 = ROOT.TLorentzVector()
    try:
        electronP4.SetXYZM(pf.getMomentum().x,pf.getMomentum().y,pf.getMomentum().z,pf.getMass())
    except AttributeError:
        # Handle the case where getMass() is not available
        electronP4.SetXYZM(pf.getMomentum().X(), pf.getMomentum().Y(), pf.getMomentum().Z(),  pf.getMass())
    if electronP4.P()<minPt:
        continue
    electronpdg = pf.getPDG()
    electroncharge = pf.getCharge()
    
    electron = RecoParticle(p4 = electronP4, ID = -11, charge = electroncharge, PDGID=electronpdg)      

    electrons[nElectrons]=electron
    nElectrons+=1
    #print ("...",pf.getObjectID().index,candTauP4.Pt(),candTauP4.Phi(),candTauP4.Theta(),candTauId,candTauCharge)

  return electrons
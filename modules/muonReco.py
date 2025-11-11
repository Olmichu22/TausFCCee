import ROOT
from modules.ParticleObjects import RecoParticle

def findAllMuons(pfos, minPt):
  """ Find all tau candidates starting from PFO collection by recognizing the decay products.

  Args:
      pfos (PandoraPFOs): PandoraPFOs collection
      minPt (float): Minimum particle momentum.
  Returns:
      muons (dict): Dictionary with the muon candidates containing RecoParticle Objects.
  """
  muons={}
  nMuons=0
  for pf in pfos:
    if (abs(pf.getPDG())!=13): 
        continue 

    
    muonP4 = ROOT.TLorentzVector()
    try:
        muonP4.SetXYZM(pf.getMomentum().x,pf.getMomentum().y,pf.getMomentum().z,pf.getMass())
    except AttributeError:
        # Handle the case where getMass() is not available
        muonP4.SetXYZM(pf.getMomentum().X(), pf.getMomentum().Y(), pf.getMomentum().Z(),  pf.getMass())
    if muonP4.P()<minPt:
        continue
    muonpdg = pf.getPDG()
    muoncharge = pf.getCharge()
    
    muon = RecoParticle(p4 = muonP4, ID = -13, charge = muoncharge, PDGID=muonpdg)      

    muons[nMuons]=muon
    nMuons+=1
    #print ("...",pf.getObjectID().index,candTauP4.Pt(),candTauP4.Phi(),candTauP4.Theta(),candTauId,candTauCharge)

  return muons
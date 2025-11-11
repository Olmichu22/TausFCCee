import ROOT
from modules.ParticleObjects import GenParticle, RecoParticle

from modules import myutils

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

# Check a generator level tau candidate, find the decay, 
# and compute visible (meson) variables 
def visTauGen(candTau):
   """ Check a generator level tau candidate, find the decay, and compute visible (meson) variables.

   Args:
       candTau (Particle ObjefindAllGenZs
   Returns:
       Tuple: Tuple with the visible 4-momentum, the tau ID, the charge, the true 4-momentum, the maximum angle between constituents, the number of constituents, and the constituents.
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

   # loop over daughter particles of the tau 
   for dTau in daughters:
         dauP4=ROOT.TLorentzVector()
         dauP4.SetXYZM(dTau.getMomentum().x,dTau.getMomentum().y,dTau.getMomentum().z,dTau.getMass())
         dauPDG=abs(dTau.getPDG())

         #print ('...dau',dauP4.P(),dauP4.Theta(),dauP4.Phi(),dau.getMass())

         # we want to compare the reco P4 to the 'visible' gen P4: skip neutrinos
         # PDG ID of Neutrinos

         if (dauPDG==12 or dauPDG==14 or dauPDG==16):
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
   return (visTauP4,tauID,chargeTau,genTauP4,maxAngleConsts,nConsts,const)                 

# Reversed procedure for reconstructed pfos
# Starting from a pion, find particles in a cone around it, and 
# build the tau 
def buildTauFromPion(lead, allPfs, DRCone=1, minP_photon=0, minP_pion=0, PNeutron=1, genminP = 0.5, charge_condition=True):
   """ Starting from a pion, find particles in a cone around it, and build the tau.

   Args:
      lead (Particle Object): Pion candidate.
      allPfs (Particle Collection): All particles in the event.
      DRCone (int, optional): Radius of the cone. Defaults to 1.
      minP_photon (int, optional): Minimum photon momentum. Defaults to 0.
      minP_pion (int, optional): Minimum pion momentum. Defaults to 0.
      PNeutron (int, optional): Minimum neutron momentum. Defaults to 10.
      genminP (int, optional): Minimum general level momentum. Defaults to 0.5.
   
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


   # set the ID: only valid combinations can be a tau (charge, constituents compatible with
   # tau decay). can be refined in the future. 
   if abs(chargeTau)==1 or not charge_condition:
      if (countPions==1 and countNeutrons==0):
         if countPhotons<10:
            tauID=countPhotons
         else:
            tauID=9

      elif (countPions==3): 
         tauID = countPhotons+10 

      elif (countPions==1 and countNeutrons>0): # Future FIXME: Pandora pion->neutron misID issue 
         # tauID=15
         tauID = -20 # To not interact with the other IDs 


      return (tauP4,tauID,chargeTau,maxConeTau,nConsts,const), None

   else:

      tauP4.SetXYZM(0,0,0,0) # safety, always return an object 
      return (tauP4,-1,0,0,0,dict()), None
                 

# loop over all gen taus 
def findAllGenTaus(mc_particles):
   """ Find all generator level taus.
   
   Args:
       mc_particles (Particle Collection): All particles in the event.
       
   Returns:
      genTaus (dict): Dictionary with the generator level taus containing tuples with the visible 4-momentum, the tau ID, and the charge.
   """
   genTaus={}
   nGenTaus=0
   for particle in mc_particles:
      # only taus
      if abs(particle.getPDG()) != 15:
         continue
      # in the pythia sample we need to check the genStatus:
      # (in some events we have several copies of the tau)
      
      # 2: final state tau (to not double count)
      if particle.getGeneratorStatus()!=2:
         continue

      genTau_data=visTauGen(particle)
      genTau = GenParticle(genTau_data[0], genTau_data[1], genTau_data[2], genTau_data[3], genTau_data[4], genTau_data[5], genTau_data[6])
      if genTau.getCharge()<0:
         genTau.setPDG(15)
      else:
         genTau.setPDG(-15)


      genTaus[nGenTaus]=genTau
      nGenTaus+=1

   return genTaus

# function to find all reco taus starting from PFO collection 
def findAllTaus(pfos, dRMax, minP_photon, minP_pion, PNeutron, genminP, charge_condition=True):
   """ Find all tau candidates starting from PFO collection by recognizing the decay products.

   Args:
      pfos (PandoraPFOs): PandoraPFOs collection
      dRMax (float): Maximum cone radius.
      minP_photon (float): Minimum photon momentum.
      minP_pion (float): Minimum pion momentum.
      PNeutron (float): Minimum neutron momentum.
      genminP (float): Minimum general level momentum.

   Returns:
       taus (dict): Dictionary with the tau candidates containing tuples with the visible 4-momentum, the tau ID, and the charge.
   """
   taus={}
   nTaus=0

   for pf in pfos:
      if (abs(pf.getPDG())!=211): 
         continue 

      pionP4 = ROOT.TLorentzVector()
      try:
         pionP4.SetXYZM(pf.getMomentum().x,pf.getMomentum().y,pf.getMomentum().z,pf.getMass())
      except AttributeError as e:
         pionP4.SetXYZM(pf.getMomentum().X(),pf.getMomentum().Y(),pf.getMomentum().Z(),pf.getMass())
         
      if pionP4.P() < minP_pion or  pionP4.P() < genminP:
         continue

      recoTau_data, pions_id = buildTauFromPion(pf, pfos, dRMax, minP_photon, minP_pion, PNeutron, genminP, charge_condition)
      recoTau = RecoParticle(recoTau_data[0], recoTau_data[1], recoTau_data[2], recoTau_data[3], recoTau_data[4], recoTau_data[5])

      if recoTau.getCharge()<0:
         recoTau.setPDG(15)
      else:
         recoTau.setPDG(-15)

      
      
      candTauP4=recoTau.getMomentum()


      duplicate=False
      for i in range(0,nTaus):
         if (myutils.dRAngle(candTauP4,taus[i].getMomentum())<0.05): duplicate=True
      if (duplicate==True): continue

      taus[nTaus]=recoTau
      nTaus+=1

   return taus


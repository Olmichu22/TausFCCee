import os, math

import ROOT
from ROOT import TFile, TTree, TH1F, TH2F, std
import numpy as np
from podio import root_io

import pprint
import yaml

from modules import tauReco, electronReco, muonReco
from modules import myutils

# ----------------------------------------------------------------------------
# Load config (necessary for set up the logger)
default_config = "config/default/taurecolong.yaml"
# Output Configuration
outputbasepath = "Results/TauReco/"

general_configs = myutils.setup_analysis_config(default_config,
                                                outputbasepath)



loggers = general_configs["loggers"]

run_config = general_configs["config"]
# Cut Configuration
dRMax = 0.4 # Max distance to look for components in the pion cone (rads)

logger_config = loggers["config"]
logger_io = loggers["io"]
logger_process = loggers["processing"]
logger_pi0mass = loggers["pi0mass"]

minPTauPhoton =run_config["cuts"]["TauPhotonPCut"] # Min Photon P Cut
minPTauPion = run_config["cuts"]["TauPionPCut"] # Min Pion P Cut
PNeutron = run_config["cuts"]["NeutronCut"] # Min Neutron P Cut
dRMatch = run_config["cuts"]["MatchedGenMaxDR"] # Min DR for start matching
generalPCut = run_config["cuts"]["generalPCut"]

selectDecay = general_configs["decay"] # Not Used

outputpath = general_configs["outputpath"] # Complete outputpath

filename = "Tree_"+ general_configs["fileOutName"] # Root file output name

fileOutName=os.path.join(general_configs["outputpath"], filename) # Complete Root file output path

# Continue with the rest of configs


# General Configuration
# ------------------------------------------------------------------------

sample=run_config["general"]["sample"] # Prefix of files to look for
test_arg = general_configs["flags"]["test"]


logger_config.info("Configuration loaded!")
logger_config.info("Configuration:\n%s", pprint.pformat(general_configs, indent=4))

# If there are results provied from neural network
gatr_results_path = general_configs["args"].gatr_result

filenames, mlpf_results = myutils.get_root_trees_path(sample, gatr_results_path, loggers, test_arg)
# mlpf results is an empty dict if no neural network results provided

if test_arg:
    logger_io.info("Running in test mode, limiting to 10 files.")
    filenames = filenames[:10]  # Limit to 10 files for testing


reader = root_io.Reader(filenames)
logger_io.info("Read %d files", len(filenames))
logger_io.info("First %s files.", filenames[:10]) 

# String to identify the performed cuts
cut_string = general_configs["decay_str"]

# Configs and reading finished
# ----------------------------------------------------------------------

# collections to use
genparts = "MCParticles"
pfobjects = "PandoraPFOs"

outfile = ROOT.TFile(fileOutName, "RECREATE")
Tau_tree = TTree("Tau_tree", f"Tree {cut_string}_Sample_{sample}")

# Defining many histogram
numGenTaus = np.array([0], dtype=int)
GenEventId = ROOT.std.vector('int')()
GenTauPt = ROOT.std.vector('float')()
GenVisTauPt = ROOT.std.vector('float')()
GenTauP = ROOT.std.vector('float')()
GenVisTauP = ROOT.std.vector('float')()
GenTauType = ROOT.std.vector('int')()
GenVisTauMass = ROOT.std.vector('float')()
GenTauQ = ROOT.std.vector('float')()
GenTauEta = ROOT.std.vector('float')()
GenTauTheta = ROOT.std.vector('float')()
GenTauDR = ROOT.std.vector('float')()
GenTauNConsts = ROOT.std.vector('int')()
GenTauNConstKey = ROOT.std.vector('int')()
GenTauConstKey = ROOT.std.vector('int')() # Permite identificar a que Tau pertenece el constituyente
GenMatchedKey = ROOT.std.vector('int')()  # Key for matched gen tau
GenConstP = ROOT.std.vector('float')()
GenConstTheta = ROOT.std.vector('float')()
GenConstEta = ROOT.std.vector('float')()
GenConstPDG = ROOT.std.vector('int')()
Tau_tree.Branch("numGenTaus", numGenTaus, "numGenTaus/I")
Tau_tree.Branch("GenEventId", GenEventId)
Tau_tree.Branch("GenTauPt", GenTauPt)
Tau_tree.Branch("GenVisTauPt", GenVisTauPt)
Tau_tree.Branch("GenTauP", GenTauP)
Tau_tree.Branch("GenVisTauP", GenVisTauP)
Tau_tree.Branch("GenTauType", GenTauType)
Tau_tree.Branch("GenVisTauMass", GenVisTauMass)
Tau_tree.Branch("GenTauQ", GenTauQ)
Tau_tree.Branch("GenTauEta", GenTauEta)
Tau_tree.Branch("GenTauTheta", GenTauTheta)
Tau_tree.Branch("GenTauDR", GenTauDR)
Tau_tree.Branch("GenTauNConsts", GenTauNConsts)
Tau_tree.Branch("GenTauNConstKey", GenTauNConstKey)
Tau_tree.Branch("GenTauConstKey", GenTauConstKey)
Tau_tree.Branch("GenMatchedKey", GenMatchedKey)
Tau_tree.Branch("GenConstP", GenConstP)
Tau_tree.Branch("GenConstTheta", GenConstTheta)
Tau_tree.Branch("GenConstEta", GenConstEta)
Tau_tree.Branch("GenConstPDG", GenConstPDG)

numRecoTaus = np.array([0], dtype=int)
RecoTauPt = ROOT.std.vector('float')()
RecoTauP = ROOT.std.vector('float')()
RecoTauMass = ROOT.std.vector('float')()
RecoTauType = ROOT.std.vector('int')()
RecoTauDM = ROOT.std.vector('int')()
RecoTauQ = ROOT.std.vector('float')()
RecoTauEta = ROOT.std.vector('float')()
RecoTauTheta = ROOT.std.vector('float')()
RecoTauDR = ROOT.std.vector('float')()
RecoTauNConsts = ROOT.std.vector('int')()
RecoTauNConstKey = ROOT.std.vector('int')()
RecoTauConstKey = ROOT.std.vector('int')()  # For identifying the parent tau
RecoMatchedKey = ROOT.std.vector('int')()  # Key for matched reco tau
RecoConstP = ROOT.std.vector('float')()
RecoConstTheta = ROOT.std.vector('float')()
RecoConstEta = ROOT.std.vector('float')()
RecoConstPDG = ROOT.std.vector('int')()
Tau_tree.Branch("RecoTauPt", RecoTauPt)
Tau_tree.Branch("RecoTauP", RecoTauP)
Tau_tree.Branch("RecoTauMass", RecoTauMass)
Tau_tree.Branch("RecoTauType", RecoTauType)
Tau_tree.Branch("RecoTauDM", RecoTauDM)
Tau_tree.Branch("RecoTauQ", RecoTauQ)
Tau_tree.Branch("RecoTauEta", RecoTauEta)
Tau_tree.Branch("RecoTauTheta", RecoTauTheta)
Tau_tree.Branch("RecoTauDR", RecoTauDR)
Tau_tree.Branch("RecoTauNConsts", RecoTauNConsts)
Tau_tree.Branch("RecoTauNConstKey", RecoTauNConstKey)
Tau_tree.Branch("RecoTauConstKey", RecoTauConstKey)
Tau_tree.Branch("RecoMatchedKey", RecoMatchedKey)
Tau_tree.Branch("RecoConstP", RecoConstP)
Tau_tree.Branch("RecoConstTheta", RecoConstTheta)
Tau_tree.Branch("RecoConstEta", RecoConstEta)
Tau_tree.Branch("RecoConstPDG", RecoConstPDG)


vectors = [
GenEventId,
GenTauPt,
GenVisTauPt,
GenTauP,
GenVisTauP,
GenTauType,
GenVisTauMass,
GenTauQ,
GenTauEta,
GenTauTheta,
GenTauDR,
GenTauNConsts,
GenTauNConstKey,
GenTauConstKey,
GenMatchedKey,
GenConstP,
GenConstTheta,
GenConstEta,
GenConstPDG,
RecoTauPt,
RecoTauP,
RecoTauMass,
RecoTauType,
RecoTauDM,
RecoTauQ,
RecoTauEta,
RecoTauTheta,
RecoTauDR,
RecoTauNConsts,
RecoTauNConstKey,
RecoTauConstKey,  
RecoMatchedKey,  
RecoConstP,
RecoConstTheta,
RecoConstEta,
RecoConstPDG,
]


countEvents = 0
# run over all events
for eventid, event in enumerate(reader.get("events")):
    logger_process.debug("Processing event %d", eventid)
    if countEvents % 1000 == 0:
        logger_process.info("Processing event %d", countEvents)
    countEvents += 1

    mc_particles = event.get(genparts)
    pfos = event.get(pfobjects)

    genTaus = tauReco.findAllGenTaus(mc_particles)
    nGenTaus = len(genTaus)

    logger_process.debug(
        "Found %d gen taus. Details:\n%s",
        nGenTaus,
        "\n".join("GenTau %d: %s" % (i, tau) for i, tau in genTaus.items()),
    )

    if gatr_results_path is not None and not general_configs["args"].test_pfo:
        particles = mlpf_results.get(eventid, {})
        recoTau = tauReco.findAllTaus(
            particles, dRMax, minPTauPhoton, minPTauPion, PNeutron, generalPCut, charge_condition=False
        )
        recoElectrons = electronReco.findAllElectrons(particles, generalPCut)
        recoMuons = muonReco.findAllMuons(particles, generalPCut)  
    else:
        recoTau = tauReco.findAllTaus(
            pfos, dRMax, minPTauPhoton, minPTauPion, PNeutron, generalPCut
        )
        recoElectrons = electronReco.findAllElectrons(pfos, generalPCut)
        recoMuons = muonReco.findAllMuons(pfos, generalPCut)
    
    nRecoTaus = len(recoTau)
    nRecoElectrons = len(recoElectrons)
    nRecoMuons = len(recoMuons)

    recoTaus = {}
    pidx = 0
    for taui in range(nRecoTaus):
        recoTaus[pidx] = recoTau[taui]
        pidx += 1
    for elei in range(nRecoElectrons):
        recoTaus[pidx] = recoElectrons[elei]
        pidx += 1
    for mui in range(nRecoMuons):
        recoTaus[pidx] = recoMuons[mui]
        pidx += 1
    nRecoTaus = len(recoTaus)
    
    logger_process.debug(
        "Found %d reconstructed taus. Details:\n%s",
        nRecoTaus,
        "\n".join("RecoTau %d: %s" % (i, tau) for i, tau in recoTaus.items()),
    )

    foundGen = False

    nGenTausType = 0
    nTausType = 0
    nGenTausHad = 0

    for v in vectors:
            v.clear()
    
    numGenTaus[0] = nGenTaus
    numRecoTaus[0] = nRecoTaus
    
    # Fill with all the information
    for i in range(0, nGenTaus):
        
        genVisTauP4 = genTaus[
            i
        ].getvisMomentum()  
        genTauId = genTaus[i].getID()
        genTauQ = genTaus[i].getCharge()
        genTauP4 = genTaus[i].getMomentum()
        genTauDR = genTaus[
            i
        ].getMaxAngle()  # Maximum angle between the tau and its constituents
        genTauNConsts = genTaus[i].getnConst()
        genTauConsts = genTaus[i].getDaughters()

        # Fill histograms
        GenEventId.push_back(eventid)  # Event ID
        GenTauPt.push_back(genTauP4.Pt())  # Transverse momentum
        GenVisTauPt.push_back(genVisTauP4.Pt())  # Visible transverse momentum
        GenTauP.push_back(genTauP4.P())  # Momentum
        GenVisTauP.push_back(genVisTauP4.P())  # Visible momentum
        GenVisTauMass.push_back(genVisTauP4.M())  # Visible mass
        GenTauType.push_back(genTauId)  # Tau decay type
        GenTauQ.push_back(genTauQ)  # Tau charge

        GenTauEta.push_back(genTauP4.Eta())  # Pseudo-rapidity
        GenTauTheta.push_back(genTauP4.Theta())  # Theta angle

        GenTauDR.push_back(genTauDR)  # Angle of Tau Constituents
        GenTauNConstKey.push_back(i)  # Key for Tau Constituents
        GenTauNConsts.push_back(genTauNConsts)  # Number of Tau Constituents

        countPionsRun = 0


        # For each generator level tau, check the constituents and fill histograms
        for c in range(0, genTauNConsts):
            const = genTauConsts[c]
            constP4 = ROOT.TLorentzVector()
            constP4.SetXYZM(
                const.getMomentum().x,
                const.getMomentum().y,
                const.getMomentum().z,
                const.getMass(),
            )
            GenTauConstKey.push_back(i)  # Add the key to identify the tau
            GenConstPDG.push_back(const.getPDG())  # PDG ID of the constituent
            GenConstP.push_back(constP4.P())  # Momentum of the constituent
            GenConstTheta.push_back(constP4.Theta())  # Theta angle of the constituent
            GenConstEta.push_back(constP4.Eta())  # Pseudo-rapidity of the constituent
        
     # Fill with all the information
    for i in range(0, nRecoTaus):
        
        recoTauP4 = recoTaus[
            i
        ].getMomentum()
        recoTauId = recoTaus[i].getID()
        recoTauQ = recoTaus[i].getCharge()
        recoTauP4 = recoTaus[i].getMomentum()
        recoTauDR = recoTaus[
            i
        ].getMaxCone()  # Maximum angle between the tau and its constituents
        recoTauNConsts = recoTaus[i].getnConst()
        recoTauConsts = recoTaus[i].getDaughters()
        
        recoDM = recoTauId
        n_pi0s = 0
        
        if recoTauId < 10 and recoTauId >= 0:
            nPhotons = recoTauId
            n_pi0s = math.ceil(nPhotons / 2)
            recoDM = n_pi0s
        elif recoTauId >= 10:
            nPhotons = recoTauId - 10
            n_pi0s = math.ceil(nPhotons / 2)
            recoDM = 10 + n_pi0s
        
        RecoTauPt.push_back(recoTauP4.Pt())  # Transverse momentum
        RecoTauP.push_back(recoTauP4.P())  # Momentum
        RecoTauMass.push_back(recoTauP4.M())  # Mass
        RecoTauType.push_back(recoTauId)  # Tau decay type
        RecoTauDM.push_back(recoDM)  # Decay mode
        RecoTauQ.push_back(recoTauQ)  # Tau charge
        RecoTauEta.push_back(recoTauP4.Eta())  # Pseudo-rapidity
        RecoTauTheta.push_back(recoTauP4.Theta())  # Theta angle
        RecoTauDR.push_back(recoTauDR)  # Angle of Tau Constituents

        RecoTauNConstKey.push_back(i)  # Key for Tau Constituents
        RecoTauNConsts.push_back(recoTauNConsts)  # Number of Tau Constituents

        countPionsRun = 0

        # print ("all GEN")
        # Look inside the generator level tau: check the constituents (decay products)

        # For each generator level tau, check the constituents and fill histograms
        for c in range(0, recoTauNConsts):
            const = recoTauConsts[c]
            constP4 = ROOT.TLorentzVector()
            try:
                constP4.SetXYZM(
                    const.getMomentum().x,
                    const.getMomentum().y,
                    const.getMomentum().z,
                    const.getMass(),
                )
            except AttributeError: # Ocurre si cargamos un objeto ParticleObject en vez de Particle
                constP4.SetXYZM(
                    const.getMomentum().X(),
                    const.getMomentum().Y(),
                    const.getMomentum().Z(),
                    const.getMass(),
                )
            RecoTauConstKey.push_back(i)  # Add the key to identify the tau
            RecoConstPDG.push_back(const.getPDG())  # PDG ID of the constituent
            RecoConstP.push_back(constP4.P())  # Momentum of the constituent
            RecoConstTheta.push_back(constP4.Theta())  # Theta angle of the constituent
            RecoConstEta.push_back(constP4.Eta())  # Pseudo-rapidity of the constituent

    
    for i in range(0, nGenTaus):
        
        genVisTauP4 = genTaus[
            i
        ].getvisMomentum()  
        genTauId = genTaus[i].getID()
        genTauQ = genTaus[i].getCharge()
        genTauP4 = genTaus[i].getMomentum()
        genTauDR = genTaus[
            i
        ].getMaxAngle()  # Maximum angle between the tau and its constituents
        genTauNConsts = genTaus[i].getnConst()
        genTauConsts = genTaus[i].getDaughters()
        
        # Compare with reconstructed taus using angle matching
        findMatch, nTausType = tauReco.MatchRecoGenTau(
            genTaus[i], recoTaus, nTausType, maxDRMatch=dRMatch, selectDecay=selectDecay
        )
        GenMatchedKey.push_back(i)  # Key for matched gen tau
        RecoMatchedKey.push_back(findMatch)  # Key for matched reco tau
        
        
    Tau_tree.Fill()
        
logger_process.info("Found %d events", countEvents)


output_config_file = outputpath + "config.yaml"
with open(output_config_file, "w") as file:
    yaml.dump(run_config, file)
    logger_io.info("Configuration file saved to %s", output_config_file)

Tau_tree.SetDirectory(outfile)

Tau_tree.Write()

logger_io.info("Output file %s", outputpath + fileOutName)
logger_io.info("End of job")
outfile.Close()

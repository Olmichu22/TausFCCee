import os, math, sys

import ROOT
from ROOT import TH1F
from podio import root_io

from modules import tauReco, electronReco, muonReco
from modules import myutils
import logging
import argparse

# ----------------------------------------------------------------------------
# Simple arg parse

parser = argparse.ArgumentParser(description="Configure the analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

parser.add_argument(
        "-v", "--verbose",
        action="count",
        default=0,
        help="Increase verbosity level: -v for INFO, -vv for DEBUG",
    )
parser.add_argument("--log-dir", default="./logs", help="Path to save the logs.")
args = parser.parse_args()

log_dir = args.log_dir
# Create path if path does not exist
os.makedirs(log_dir, exist_ok=True)

# ----------------------------------------------------------------------------
# Config logging
lvl = logging.WARNING if args.verbose == 0 else logging.INFO if args.verbose == 1 else logging.DEBUG
handlers = []
if args.verbose < 2:
    handlers = [logging.StreamHandler(sys.stdout), logging.FileHandler(os.path.join(log_dir, "app.log"), mode="w")]
elif args.verbose == 2:
    sh = logging.StreamHandler(sys.stdout); sh.setLevel(logging.DEBUG)
    fh = logging.FileHandler(os.path.join(log_dir, "app.log"), mode="w"); fh.setLevel(logging.DEBUG)
    handlers = [sh, fh]
else:
    handlers = [logging.FileHandler(os.path.join(log_dir, "app.log"), mode="w")]

logging.basicConfig(
    level=lvl,
    format="%(asctime)s, %(levelname)s, [%(name)s] - %(message)s",
    handlers=handlers
)

loggers = {
    "config": logging.getLogger("config"),
    "io": logging.getLogger("io"),
    "processing": logging.getLogger("processing"),
}
loggers["config"].info("Logging set up!")

logger_config = loggers["config"]
logger_io = loggers["io"]
logger_process = loggers["processing"]


# Cut Configuration
dRMax=0.4 # Max distance to look for components in the pion cone (rads)


minPTauPhoton = 0. # Min Photon P Cut
minPTauPion = 0. # Min Pion P Cut
PNeutron = 0. # Min Neutron P Cut
dRMatch = 1 # Max DR for start gen-reco matching
generalPCut = 0. # General Momentum Cut

selectDecay = -777 # Not Used (all decays)

outputpath = "Results/" # Complete outputpath
os.makedirs(outputpath, exist_ok=True)

filename = "TauRecoTestHist.root" # Root file output name

fileOutName = os.path.join(outputpath, filename) # Complete Root file output path

# Continue with the rest of configs


# General Configuration
# ------------------------------------------------------------------------

# ------------------------------------- INSTRUCTIONS -----------------------------
# If more than one file provided, use the combination of path + file_prefix to read all the files (if test is false)
# path -> main path of the root files
# file_prefix -> files must be named as {file_prefix}_{id}.root with id = [1....N]

# If you want to read one file, use the function "get_single_root_file".

test_arg = False # True to read less files

# USE THIS FLAG TO SELECT THE CASE

read_one_file = True

if not read_one_file:
    # More than one file
    path = "path/to/root/files"
    file_prefix = "root_file_prefix"

    filenames, mlpf_results = myutils.get_root_trees_path(sample="",
                                                        gatr_results_path=None,
                                                        loggers=loggers,
                                                        test_arg=test_arg,
                                                        path=path,
                                                        file_prefix=file_prefix
                                                        )
    # mlpf results is an empty dict as no neural network results provided
else:
    file_path = "/pnfs/ciemat.es/data/cms/store/user/cepeda/FCC/FullSim/ZTauTau_SMPol_25Sept_MuonFix/out_reco_edm4hep_edm4hep_1.root"
    filenames = myutils.get_single_root_file(file_path=file_path, loggers=loggers)

if test_arg:
    logger_io.info("Running in test mode, limiting to 10 files.")
    filenames = filenames[:10]  # Limit to 10 files for testing


reader = root_io.Reader(filenames)
logger_io.info("Read %d files", len(filenames))
logger_io.info("First %s files.", filenames[:10]) 


# Configs and reading finished
# ----------------------------------------------------------------------

# Test Histograms
# Gen Level
hGenTauPt = TH1F("histoGenTauPt", "", 250, 0, 50)
hGenTauVisPt = TH1F("histoGenTauVisPt", "", 250, 0, 50)
hGenTauP = TH1F("histoGenTauP", "", 250, 0, 50) 
hGenTauVisP = TH1F("histoGenTauVisP", "", 250, 0, 50)

# Decay Type
hGenTauType = TH1F("histoGenTauType", "", 41, -21, 20)

# Mass
hGenTauMass = TH1F("histoGenTauMass", "", 500, 0, 10)
hGenTauVisMass = TH1F("histoGenTauVisMass", "", 500, 0, 10)

# Reco Level
hRecoTauPt = TH1F("histoRecoTauPt", "", 250, 0, 50)
hRecoTauP = TH1F("histoRecoTauP", "", 250, 0, 50) 

# Decay Type
hRecoTauType = TH1F("histoRecoTauType", "", 41, -21, 20)

# Mass
hRecoTauMass = TH1F("histoRecoTauMass", "", 500, 0, 10)


# Resolution
hResTauP = TH1F("histoResTauP", "", 500, -1, 1)



# collections to use
genparts = "MCParticles"
pfobjects = "PandoraPFOs"


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


    recoTau = tauReco.findAllTaus( # Hadronic decays
        pfos, dRMax, minPTauPhoton, minPTauPion, PNeutron, generalPCut
    )
    recoElectrons = electronReco.findAllElectrons(pfos, generalPCut) # Leptonic decay -> e
    recoMuons = muonReco.findAllMuons(pfos, generalPCut) # Leptonic decay -> mu
    
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

    
    # Loop over gen taus
    for i in range(0, nGenTaus):
        
        genVisTauP4 = genTaus[
            i
        ].getvisMomentum()  
        genTauId = genTaus[i].getID() # This is the identificator of the decay
        genTauQ = genTaus[i].getCharge()
        genTauP4 = genTaus[i].getMomentum()
        genTauDR = genTaus[
            i
        ].getMaxAngle()  # Maximum angle between the tau and its constituents
        genTauNConsts = genTaus[i].getnConst()
        genTauConsts = genTaus[i].getDaughters()
        genTauMass = genTaus[i].getMass()
        genTauVisMass = genTaus[i].getVisMass()
        
        if genTauVisMass < 0.2:
            print(genTaus[i])
        # Fill histograms
        hGenTauPt.Fill(genTauP4.Pt())
        hGenTauVisPt.Fill(genVisTauP4.Pt())
        hGenTauP.Fill(genTauP4.P()) 
        hGenTauVisP.Fill(genVisTauP4.P())

        # Decay Type
        hGenTauType.Fill(genTauId)

        # Mass
        hGenTauMass.Fill(genTauMass)
        hGenTauVisMass.Fill(genTauVisMass)


        # NOTE Objects at genTaus and recoTaus are GenParticle and RecoParticle (modules/ParticleObjects.py)
        # Constituents are PandoraPFOs, which have slightly different methods
        
        # For each generator level tau, we can loop over constituents
        for c in range(0, genTauNConsts):
            const = genTauConsts[c]
            constP4 = ROOT.TLorentzVector()
            constP4.SetXYZM(
                const.getMomentum().x, # For GenParticle or RecoParticle is getMomentum().X()
                const.getMomentum().y,
                const.getMomentum().z,
                const.getMass(), # This is the same for GenParticle or RecoParticle 
            )
            
        
    # Loop over reco taus
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
        recoTauMass = recoTaus[i].getMass()
        recoDM = recoTauId
        n_pi0s = 0
        
        # Infer number of pi0s based on number of photons
        if recoTauId < 10 and recoTauId >= 0:
            nPhotons = recoTauId
            n_pi0s = math.ceil(nPhotons / 2)
            recoDM = n_pi0s
        elif recoTauId >= 10:
            nPhotons = recoTauId - 10
            n_pi0s = math.ceil(nPhotons / 2)
            recoDM = 10 + n_pi0s
        
        # Fill the histograms
        # Reco Level
        hRecoTauPt.Fill(recoTauP4.Pt())
        hRecoTauP.Fill(recoTauP4.P())

        # Decay Type
        hRecoTauType.Fill(recoDM) #recoTauId has a different codification
                                    #recoDM is the equivalent.

        # Mass
        hRecoTauMass.Fill(recoTauMass)


        # For each reco level tau, we can loop over constituents
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
            except AttributeError: # IGNORE THIS
                constP4.SetXYZM(
                    const.getMomentum().X(),
                    const.getMomentum().Y(),
                    const.getMomentum().Z(),
                    const.getMass(),
                )
    # Try to match gen a reco taus
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
        if findMatch == -1:
            continue
            # Gen Tau with key i will be matched with reco tau with key findMatch. (-1) means unmatched.
        recoTau = recoTaus[findMatch]
        recotauP4 = recoTau.getMomentum()
        hResTauP.Fill((recoTauP4.P() - genVisTauP4.P()) / genVisTauP4.P())
        
        
logger_process.info("Found %d events", countEvents)


outfile = ROOT.TFile(fileOutName, "RECREATE")

hGenTauPt.Write()
hGenTauVisPt.Write()
hGenTauP.Write()
hGenTauVisP.Write()

# Decay Type
hGenTauType.Write()

# Mass
hGenTauMass.Write()
hGenTauVisMass.Write()

# Reco Level
hRecoTauPt.Write()
hRecoTauP.Write()

# Decay Type
hRecoTauType.Write()

# Mass
hRecoTauMass.Write()


# Resolution
hResTauP.Write()


logger_io.info("Output file %s", outputpath + fileOutName)
logger_io.info("End of job")
outfile.Close()

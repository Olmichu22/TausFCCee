import ROOT
from modules import myutils


class Track:
    """Reconstructed track class to store momentum and track object.
    Attributes:
      p4 (TLorentzVector): 4-momentum of the particle.
      associations (list): List of associated objects (e.g., hits, clusters).
    """

    def __init__(self, p4=None, associations: list = []):

        self.p4 = p4 if p4 is not None else ROOT.TLorentzVector(0, 0, 0, 0)
        self.associations = associations

    def getMomentum(self):
        return self.p4

    def getAssociations(self):
        return self.associations

    def setMomentum(self, p4):
        self.p4 = p4

    def setAssociations(self, associations: list):
        self.associations = associations

    def __str__(self):
        return "Track with momentum: (%f, %f, %f, %f) and %d associations" % (
            self.p4.Px(),
            self.p4.Py(),
            self.p4.Pz(),
            self.p4.E(),
            len(self.associations),
        )

    def __repr__(self):
        return self.__str__()


class Particle:
    """Standard particle class to store the information of a particle.
    Attributes:
        p4 (TLorentzVector): 4-momentum of the particle.
        PDGID (int): PDG ID of the particle.
        ID (int): Decay ID of the particle.
        charge (int): Charge of the particle.
        idx (int): Index of the particle in a collection.
    """

    def __init__(self, p4=None, PDGID=-1, ID=-1, charge=0, idx=-1):
        """Constructor of the Particle class."""
        self.p4 = p4 if p4 is not None else ROOT.TLorentzVector(0, 0, 0, 0)
        self.PDGID = PDGID
        self.ID = ID
        self.charge = charge
        self.idx = idx

    def getIdx(self):
        return self.idx

    def getPDG(self):
        return self.PDGID

    def getID(self):
        return self.ID

    def getCharge(self):
        return self.charge

    def getMomentum(self):
        return self.p4

    def getMass(self):
        return self.p4.M()

    def setIdx(self, idx):
        self.idx = idx

    def setPDG(self, PDGID):
        self.PDGID = PDGID

    def setID(self, ID):
        self.ID = ID

    def setCharge(self, charge):
        self.charge = charge

    def setMomentum(self, p4):
        self.p4 = p4

    def ShowInfo(self):
        """Print the information of the particle."""
        print(
            "Particle PDG: %d, ID: %d, Charge: %d, Mass: %f"
            % (self.PDGID, self.ID, self.charge, self.p4.M())
        )
        print("Particle Index: %d" % self.idx)

    def __str__(self):
        return "PDG: %d, ID: %d, Charge: %d, Mass: %f, Index: %d" % (
            self.PDGID,
            self.ID,
            self.charge,
            self.p4.M(),
            self.idx,
        )


class RecoParticle(Particle):
    """Particle class to store the information of a reco level particle.
    Attributes:
        p4 (TLorentzVector): 4-momentum of the particle.
        ID (int): Decay ID of the particle.
        charge (int): Charge of the particle.
        maxCone (float): Maximum angle between the guide particle and the rest of the constituents.
        nConst (int): Number of constituents of the particle.
        const (dict): Dictionary with the constituents of the particle.
        PDGID (int):  PDG ID of the particle.
        idx (int): Unique identificator for the particle (per event).
    """

    def __init__(
        self,
        p4=None,
        ID=-1,
        charge=0,
        maxCone=0,
        nConst=0,
        const=None,
        PDGID=-1,
        idx=-1,
    ):
        """Constructor of the Particle class."""
        # Initialize the parent class
        super(RecoParticle, self).__init__(p4, PDGID, ID, charge, idx)

        self.maxCone = maxCone
        self.nConst = nConst
        self.const = const if const is not None else {}

    def getMomentum(self):
        return self.p4

    def getMass(self):
        return self.p4.M()

    def getID(self):
        return self.ID

    def getPDG(self):
        return self.PDGID

    def getCharge(self):
        return self.charge

    def getMaxCone(self):
        return self.maxCone

    def getnConst(self):
        return self.nConst

    def getDaughters(self):
        return self.const

    def setMomentum(self, p4):
        self.p4 = p4

    def setID(self, ID):
        self.ID = ID

    def setPDG(self, PDGID):
        self.PDGID = PDGID

    def setCharge(self, charge):
        self.charge = charge

    def setMaxCone(self, maxCone):
        self.maxCone = maxCone

    def setDaughters(self, const):
        self.const = const
        self.nConst = len(const)

    def addDaughter(self, daughter):
        self.const[self.nConst] = daughter
        self.nConst += 1

    def ShowInfo(self):
        """Print the information of the particle and its constituents."""
        print(
            "RecoParticle ID: %d, Charge: %d, Mass: %f, nConst: %d, Angle: %f, Index: %d"
            % (self.ID, self.charge, self.p4.M(), self.nConst, self.maxCone, self.idx)
        )
        for i in range(self.nConst):
            print("  Constituent %d:" % i)
            print(self.const[i])

    def __str__(self):
        return "ID: %d, Charge: %d, Mass: %f, nConst: %d, Angle: %f" % (
            self.ID,
            self.charge,
            self.p4.M(),
            self.nConst,
            self.maxCone,
        )

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return self is other
        if isinstance(other, GenRecoMatched):
            return other.matchedRecoID == self.getIdx()
        return False


class GenParticle(Particle):
    """Particle class to store the information of a gen level particle.
    Attributes:
        visp4 (TLorentzVector): 4-momentum of the visible part of the particle.
        ID (int): For daughters, it will be the PID. For taus, it is the decay identificator.
        charge (int): Charge of the particle.
        visCharge (int): Charge of the visible part of the particle.
        genP4 (TLorentzVector): 4-momentum of the particle.
        maxAngleConsts (float): Maximum angle between the constituents of the particle.
        nConst (int): Number of constituents of the particle.
        const (dict): Dictionary with the constituents of the particle.
        PDGID (int):  PDG ID of the particle.
        mcp (MCParticle): Original MCParticle object from the simulation.
        idx (int): Unique identificator for the particle (per event).
    """

    def __init__(
        self,
        visP4=None,
        ID=-1,
        charge=0,
        genP4=None,
        maxAngleConsts=0,
        nConsts=0,
        const=None,
        PDGID=-1,
        mcp=None,
        idx=-1,
    ):
        """Constructor of the Particle class."""
        # Initialize the parent class
        super(GenParticle, self).__init__(genP4, PDGID, ID, charge, idx)

        self.visp4 = visP4 if visP4 is not None else ROOT.TLorentzVector(0, 0, 0, 0)
        self.maxAngle = maxAngleConsts
        self.nConst = nConsts
        self.const = const if const is not None else {}
        self.mcp = mcp

    def getPDG(self):
        return self.PDGID

    def getID(self):
        return self.ID

    def getCharge(self):
        return self.charge

    def getMomentum(self):
        return self.p4

    def getvisMomentum(self):
        return self.visp4

    def getMass(self):
        return self.p4.M()

    def getVisMass(self):
        return self.visp4.M()

    def getDaughters(self):
        return self.const

    def getnConst(self):
        return self.nConst

    def getMaxAngle(self):
        return self.maxAngle

    def setID(self, ID):
        self.ID = ID

    def setPDG(self, PDGID):
        self.PDGID = PDGID

    def setCharge(self, charge):
        self.charge = charge

    def setMomentum(self, candP):
        self.p4.SetXYZM(
            candP.getMomentum().x,
            candP.getMomentum().y,
            candP.getMomentum().z,
            candP.getMass(),
        )

    def setvisMomentum(self, visP4):
        self.visp4 = visP4

    def setDaughters(self, const):
        self.const = const
        self.nConst = len(const)

    def addDaughter(self, daughter):
        self.const[self.nConst] = daughter
        self.nConst += 1

    def setnConst(self, nConst):
        self.nConst = nConst

    def setMaxAngle(self, maxAngle):
        self.maxAngle = maxAngle

    def ShowInfo(self):
        """Print the information of the particle and its constituents."""
        print(
            "GenParticle ID: %d, Charge: %d, Mass: %f, VisMass: %f, nConst: %d, Angle: %f, Index: %d"
            % (
                self.ID,
                self.charge,
                self.p4.M(),
                self.visp4.M(),
                self.nConst,
                self.maxAngle,
                self.idx,
            )
        )
        for i in range(self.nConst):
            print("  Constituent %d:" % i)
            print(self.const[i])

    def __str__(self):
        return (
            "ID: %d, Charge: %d, Mass: %f, VisMass: %f, nConst: %d, Angle: %f, Index: %d"
            % (
                self.ID,
                self.charge,
                self.p4.M(),
                self.visp4.M(),
                self.nConst,
                self.maxAngle,
                self.idx,
            )
        )

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return self is other
        if isinstance(other, GenRecoMatched):
            return other.genID == self.getIdx()
        return False


class GenRecoMatched:
    """Container for the results of particle matching.
    Attributes:
        genID (int): ID of the generator particle.
        genParticle (GenParticle): The generator particle object.
        genPDG (int): PDG code of the generator particle.
        recoIDs (List[int]): List of IDs of the reconstructed particles matched to the generator particle.
        recoParticles (List[RecoParticle]): List of reconstructed particle objects matched to the generator particle.
        recoPDGs (List[int]): List of PDG codes of the reconstructed particles matched to the generator particle.
        dRs (List[float]): List of angular distances (dR) between the generator particle and each reconstructed particle.
        minDR (float): Minimum angular distance (dR) between the generator particle and the nearest reconstructed particle.
        nearestRecoParticle (RecoParticle): The nearest reconstructed particle to the generator particle.
        nearestRecoID (int): ID of the nearest reconstructed particle.
        nearestRecoPDG (int): PDG code of the nearest reconstructed particle.
        matchedRecoParticle (RecoParticle): The reconstructed particle matched to the generator particle (initially the nearest).
        matchedRecoID (int): ID of the matched reconstructed particle.
        matchedRecoPDG (int): PDG code of the matched reconstructed particle.
        matchedDR (float): Angular distance (dR) between the generator particle and the matched reconstructed particle.
    Methods:
        getGenID(): Returns the ID of the generator particle.
        getGenParticle(): Returns the generator particle object.
        getGenPDG(): Returns the PDG code of the generator particle.
        getRecoIDs(): Returns the list of IDs of the reconstructed particles.
        getRecoParticles(): Returns the list of reconstructed particle objects.
        getdRs(): Returns the list of angular distances (dR) between the generator particle and each reconstructed particle.
        getMinDR(): Returns the minimum angular distance (dR) between the generator particle and the nearest reconstructed particle.
        getNearestRecoParticle(): Returns the nearest reconstructed particle to the generator particle.
        getNearestRecoID(): Returns the ID of the nearest reconstructed particle.
        getNearestRecoPDG(): Returns the PDG code of the nearest reconstructed particle.
        getMatchedRecoParticle(): Returns the reconstructed particle matched to the generator particle.
        getMatchedRecoID(): Returns the ID of the matched reconstructed particle.
        getMatchedRecoPDG(): Returns the PDG code of the matched reconstructed particle.
        getMatchedDR(): Returns the angular distance (dR) between the generator particle and the matched reconstructed particle.
        setMatchedRecoParticle(idx, recoParticle=None): Sets the matched reconstructed particle by ID or by providing a RecoParticle object.
        __eq__(other): Compares the GenRecoMatched object with another GenParticle or RecoParticle based on their IDs.
    Note:
        The class assumes that the input lists (recoIDs, recoParticles, dRs) are of the same length and contain at least one element.
    """

    def __init__(
        self, genID, genParticle, recoIDs: list, recoParticles: list, dRs: list = None
    ):

        assert len(recoIDs) == len(
            recoParticles
        ), "recoIDs and recoParticles must have the same length"
        assert (
            len(recoIDs) > 0
        ), "recoIDs and recoParticles must have at least one element"
        if dRs is not None:
            assert len(dRs) == len(
                recoParticles
            ), "dRs must have the same length as recoParticles"

        self.genID = genID
        self.genParticle: GenParticle = genParticle
        self.genPDG = genParticle.getPDG()

        self.recoIDs = recoIDs
        self.recoParticles = recoParticles

        self.dRs = (
            dRs
            if dRs is not None
            else [
                myutils.dRAngle(genParticle.getvisMomentum(), rp.getMomentum())
                for rp in recoParticles
            ]
        )

        # sort recoIDs, recoParticles and dRs by dRs
        sorted_indices = sorted(range(len(self.dRs)), key=lambda i: self.dRs[i])
        self.dRs = [self.dRs[i] for i in sorted_indices]
        self.recoIDs = [self.recoIDs[i] for i in sorted_indices]
        self.recoParticles = [self.recoParticles[i] for i in sorted_indices]
        self.recoPDGs = [rp.getPDG() for rp in recoParticles]

        self.minDR = self.dRs[0]
        self.nearestRecoParticle: RecoParticle = self.recoParticles[0]
        self.nearestRecoID = self.recoIDs[0]
        self.nearestRecoPDG = self.recoPDGs[0]

        self.matchedRecoParticle: RecoParticle = self.nearestRecoParticle
        self.matchedRecoID = self.nearestRecoID
        self.matchedRecoPDG = self.nearestRecoPDG
        self.matchedDR = self.minDR

    def getGenID(self) -> int:
        return self.genID

    def getGenParticle(self) -> GenParticle:
        return self.genParticle

    def getGenPDG(self) -> int:
        return self.genPDG

    def getRecoIDs(self) -> list[int]:
        return self.recoIDs

    def getRecoParticles(self) -> list[RecoParticle]:
        return self.recoParticles

    def getdRs(self) -> list[float]:
        return self.dRs

    def getMinDR(self) -> float:
        return self.minDR

    def getNearestRecoParticle(self) -> RecoParticle:
        return self.nearestRecoParticle

    def getNearestRecoID(self) -> int:
        return self.nearestRecoID

    def getNearestRecoPDG(self) -> int:
        return self.nearestRecoPDG

    def getMatchedRecoParticle(self) -> RecoParticle:
        return self.matchedRecoParticle

    def getMatchedRecoID(self) -> int:
        return self.matchedRecoID

    def getMatchedRecoPDG(self) -> int:
        return self.matchedRecoPDG

    def getMatchedDR(self) -> float:
        return self.matchedDR

    def setMatchedRecoParticle(self, idx: int, recoParticle: RecoParticle = None):
        self.matchedRecoID = idx
        if idx in self.recoIDs:
            i = self.recoIDs.index(idx)
            self.matchedRecoParticle = self.recoParticles[i]
            self.matchedRecoPDG = self.recoPDGs[i]
            self.matchedDR = self.dRs[i]

        elif recoParticle is not None:
            self.matchedRecoParticle = recoParticle
            self.matchedRecoPDG = recoParticle.getPDG()
            self.matchedDR = myutils.dRAngle(
                self.genParticle.getMomentum(), recoParticle.getMomentum()
            )
            self.recoIDs.append(idx)
            self.recoParticles.append(recoParticle)
            self.recoPDGs.append(recoParticle.getPDG())
            self.dRs.append(self.matchedDR)
        else:
            raise ValueError(
                "Either idx must be in recoIDs or recoParticle must be provided"
            )

    def __eq__(self, other):
        if isinstance(other, GenParticle):
            return self.genID == other.getIdx()
        elif isinstance(other, RecoParticle):
            return self.matchedRecoID == other.getIdx()
        return False

    def __str__(self):
        return f"GenRecoMatched(GenID: {self.genID}, GenPDG: {self.genPDG}, MatchedRecoID: {self.matchedRecoID}, MatchedRecoPDG: {self.matchedRecoPDG}, MatchedDR: {self.matchedDR})"

    def __repr__(self):
        return self.__str__()

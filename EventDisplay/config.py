"""config.py — All configurable constants for the FCC Event Display.

The display supports several detector concepts (CLD, ILD, …). Each is described
by a *profile* below (collection names + geometry). Call ``select_detector()``
to switch the active profile; that reassigns the module-level constants that the
rest of the package reads (``DETECTOR_GROUPS``, ``MC_COLLECTION``, …). The
active profile defaults to CLD (see the ``select_detector("CLD")`` call at the
bottom of this file), preserving the original behaviour.
"""

# ── Detector profiles ─────────────────────────────────────────────────────────
# Each profile maps the display's abstract concepts to the concrete EDM4HEP
# collection names + geometry of one detector.
#
#   detector_groups        : abstract group → {collections: [...], color}
#                            The "tracker" group is drawn from SimTrackerHits /
#                            tracks (see below), not read as calo hits.
#   mc / pfo / track        : names of the MCParticle / PFO / Track collections.
#   calo_link_collection    : CalorimeterHit → MCParticle link (for gen-coloring
#                            of calo hits). "" if the detector has no such link.
#   sim_tracker_collections : SimTrackerHitCollection names (getEDep + MCParticle
#                            link) drawn as tracker hits.
#   truth_link_variants     : RecoMCTruthLink-like collections to try, in order.
#   geometry_versions       : version tag → compact-XML path relative to k4geo_DIR.
#   geometry_fallback       : version tag → approximate barrel/endcap dims (mm)
#                            used when the XML can't be parsed. Dims that differ
#                            between versions (e.g. ILD TPC inner radius) live here.

_CLD_PROFILE = {
    "detector_groups": {
        "tracker": {"collections": ["SiTracks_Refitted"], "color": "#9b59b6"},
        "ecal":    {"collections": ["ECALBarrel", "ECALEndcap", "ECALOther"], "color": "#f39c12"},
        "hcal":    {"collections": ["HCALBarrel", "HCALEndcap", "HCALOther"], "color": "#2ecc71"},
        "muon":    {"collections": ["MUON"], "color": "#3498db"},
    },
    "mc_collection":       "MCParticles",
    "pfo_collection":      "PandoraPFOs",
    "track_collection":    "SiTracks_Refitted",
    "calo_link_collection": "CalohitMCTruthLink",
    "sim_tracker_collections": [
        "VertexBarrelHits", "VertexEndcapHits",
        "InnerTrackerBarrelHits", "InnerTrackerEndcapHits",
        "OuterTrackerBarrelHits", "OuterTrackerEndcapHits",
    ],
    "truth_link_variants": ["RecoMCTruthLink", "PandoraPFOsToMCParticles"],
    "geometry_versions": {
        "v06": "FCCee/CLD/compact/CLD_o2_v06/CLD_o2_v06.xml",
    },
    "default_geometry_version": "v06",
    # Approximate CLD geometry (mm) used when XML is not available
    "geometry_fallback": {
        "v06": {
            "ecal_barrel": {"rmin": 1470.0, "rmax": 1788.0, "zhalf": 2350.0},
            "ecal_endcap": {"rmin": 200.0,  "rmax": 2088.0, "zpos":  2450.0, "zthick": 400.0},
            "hcal_barrel": {"rmin": 1808.0, "rmax": 3188.0, "zhalf": 2850.0},
            "hcal_endcap": {"rmin": 300.0,  "rmax": 3188.0, "zpos":  2900.0, "zthick": 1100.0},
            "muon_barrel": {"rmin": 4200.0, "rmax": 4600.0, "zhalf": 5800.0},
        },
    },
}

# Calo/yoke dims are identical in ILD_FCCee v01 and v02 (k4geo DD4hep constants);
# only the TPC inner radius differs (v01: 365 mm, v02: 701 mm).
_ILD_CALO_DIMS = {
    "ecal_barrel": {"rmin": 1804.8, "rmax": 2028.0, "zhalf": 2350.0},
    "ecal_endcap": {"rmin": 400.0,  "rmax": 2095.8, "zpos":  2411.8, "zthick": 223.2},
    "hcal_barrel": {"rmin": 2058.0, "rmax": 3395.5, "zhalf": 2350.0},
    "hcal_endcap": {"rmin": 350.0,  "rmax": 3225.5, "zpos":  2650.0, "zthick": 1287.0},
    "muon_barrel": {"rmin": 4475.0, "rmax": 7776.0, "zhalf": 4047.0},
}

_ILD_PROFILE = {
    "detector_groups": {
        "tracker": {"collections": ["MarlinTrkTracks"], "color": "#9b59b6"},
        "ecal":    {"collections": ["EcalBarrelCollectionRec", "EcalEndcapsCollectionRec",
                                    "EcalEndcapRingCollectionRec"], "color": "#f39c12"},
        "hcal":    {"collections": ["HcalBarrelCollectionRec", "HcalEndcapsCollectionRec",
                                    "HcalEndcapRingCollectionRec"], "color": "#2ecc71"},
        "muon":    {"collections": ["MUON"], "color": "#3498db"},
    },
    "mc_collection":       "MCParticles",
    "pfo_collection":      "PandoraPFOs",
    "track_collection":    "MarlinTrkTracks",
    # ILD has no direct CalorimeterHit→MCParticle link (only Cluster↔MCParticle
    # and CaloHit↔SimCaloHit); leave empty → calo hits are not gen-colored.
    "calo_link_collection": "",
    "sim_tracker_collections": [
        "VertexBarrelCollection", "VertexEndcapCollection",
        "InnerTrackerBarrelCollection", "InnerTrackerEndcapCollection",
        "TPCCollection", "SETCollection",
    ],
    "truth_link_variants": ["RecoMCTruthLink", "MCTruthRecoLink"],
    "geometry_versions": {
        "v01": "FCCee/ILD_FCCee/compact/ILD_FCCee_v01/ILD_FCCee_v01.xml",
        "v02": "FCCee/ILD_FCCee/compact/ILD_FCCee_v02/ILD_FCCee_v02.xml",
    },
    "default_geometry_version": "v02",
    # Approximate ILD_FCCee geometry (mm), resolved from k4geo DD4hep constants
    "geometry_fallback": {
        "v01": {**_ILD_CALO_DIMS,
                "tpc_barrel": {"rmin": 365.0, "rmax": 1769.8, "zhalf": 2350.0}},
        "v02": {**_ILD_CALO_DIMS,
                "tpc_barrel": {"rmin": 701.0, "rmax": 1769.8, "zhalf": 2350.0}},
    },
}

DETECTOR_PROFILES: dict = {
    "CLD": _CLD_PROFILE,
    "ILD": _ILD_PROFILE,
}

# ── Active detector state (set by select_detector) ────────────────────────────
# These names are rebound by select_detector() and read by data_io/geometry/app.
ACTIVE_DETECTOR:         str  = "CLD"
GEOMETRY_VERSION:        str  = "v06"
DETECTOR_GROUPS:         dict = {}
COLLECTION_TO_GROUP:     dict = {}
MC_COLLECTION:           str  = ""
PFO_COLLECTION:          str  = ""
TRACK_COLLECTION:        str  = ""
CALO_LINK_COLLECTION:    str  = ""
SIM_TRACKER_COLLECTIONS: list = []
TRUTH_LINK_VARIANTS:     list = []
GEOMETRY_RELPATH:        str  = ""
GEOMETRY_FALLBACK:       dict = {}


def select_detector(name: str, geometry_version: str | None = None) -> None:
    """Activate a detector profile, rebinding the module-level constants.

    Args:
        name: profile key, case-insensitive ("CLD" or "ILD").
        geometry_version: geometry tag (e.g. "v01"/"v02" for ILD); if None or
            unknown, the profile's default version is used.
    """
    global ACTIVE_DETECTOR, GEOMETRY_VERSION
    global DETECTOR_GROUPS, COLLECTION_TO_GROUP
    global MC_COLLECTION, PFO_COLLECTION, TRACK_COLLECTION, CALO_LINK_COLLECTION
    global SIM_TRACKER_COLLECTIONS, TRUTH_LINK_VARIANTS
    global GEOMETRY_RELPATH, GEOMETRY_FALLBACK

    key = name.upper()
    if key not in DETECTOR_PROFILES:
        raise ValueError(
            f"Unknown detector '{name}'. Choose one of {sorted(DETECTOR_PROFILES)}."
        )
    p = DETECTOR_PROFILES[key]

    ACTIVE_DETECTOR         = key
    DETECTOR_GROUPS         = p["detector_groups"]
    COLLECTION_TO_GROUP     = {
        col: grp
        for grp, info in DETECTOR_GROUPS.items()
        for col in info["collections"]
    }
    MC_COLLECTION           = p["mc_collection"]
    PFO_COLLECTION          = p["pfo_collection"]
    TRACK_COLLECTION        = p["track_collection"]
    CALO_LINK_COLLECTION    = p["calo_link_collection"]
    SIM_TRACKER_COLLECTIONS = p["sim_tracker_collections"]
    TRUTH_LINK_VARIANTS     = p["truth_link_variants"]

    versions = p["geometry_versions"]
    ver = geometry_version if geometry_version in versions else p["default_geometry_version"]
    GEOMETRY_VERSION        = ver
    GEOMETRY_RELPATH        = versions[ver]
    GEOMETRY_FALLBACK       = p["geometry_fallback"][ver]


# ── Particle filters ──────────────────────────────────────────────────────────
NEUTRINO_PDGS    = {12, 14, 16}
VALID_GEN_STATUS = {1}

# PDG → human-readable name
PDG_NAMES: dict = {
    11:    "e⁻",    -11:   "e⁺",
    13:    "μ⁻",    -13:   "μ⁺",
    15:    "τ⁻",    -15:   "τ⁺",
    12:    "νe",    -12:   "ν̄e",
    14:    "νμ",    -14:   "ν̄μ",
    16:    "ντ",    -16:   "ν̄τ",
    22:    "γ",
    111:   "π⁰",
    211:   "π⁺",   -211:  "π⁻",
    130:   "K⁰L",
    310:   "K⁰S",
    321:   "K⁺",   -321:  "K⁻",
    2212:  "p",    -2212: "p̄",
    2112:  "n",    -2112: "n̄",
    3122:  "Λ",    -3122: "Λ̄",
    3222:  "Σ⁺",   -3222: "Σ̄⁻",
    3112:  "Σ⁻",   -3112: "Σ̄⁺",
    411:   "D⁺",   -411:  "D⁻",
    421:   "D⁰",   -421:  "D̄⁰",
    431:   "Ds⁺",  -431:  "Ds⁻",
    521:   "B⁺",   -521:  "B⁻",
    511:   "B⁰",   -511:  "B̄⁰",
}


def pdg_label(pdg: int) -> str:
    """Return 'name (pdg)' string, e.g. 'e⁻ (11)'."""
    name = PDG_NAMES.get(pdg, f"PDG {pdg}")
    return f"{name} ({pdg})"


# ── Themes ────────────────────────────────────────────────────────────────────
THEMES: dict = {
    "dark": {
        "app_bg":    "#1a1a2e",
        "card_bg":   "#16213e",
        "plot_bg":   "#0f3460",
        "sidebar_bg":"#0d1b35",
        "text":      "#e0e0e0",
        "muted":     "#888888",
        "accent":    "#e94560",
        "grid":      "#334455",
        "border":    "#334",
        "tab_bg":    "#16213e",
        "tab_sel":   "#0f3460",
        "input_bg":  "#0f1a30",
        "btn_bg":    "#e94560",
        "btn_text":  "#ffffff",
        "table_header_bg":  "#0f3460",
        "table_header_txt": "#e0e0e0",
        "table_cell_bg":    "#16213e",
        "table_cell_txt":   "#e0e0e0",
        "table_border":     "#334",
    },
    "light": {
        "app_bg":    "#f0f4f8",
        "card_bg":   "#ffffff",
        "plot_bg":   "#eef2ff",
        "sidebar_bg":"#e8edf5",
        "text":      "#1a1a2e",
        "muted":     "#555555",
        "accent":    "#c0392b",
        "grid":      "#cccccc",
        "border":    "#d0d0d0",
        "tab_bg":    "#ffffff",
        "tab_sel":   "#dce6ff",
        "input_bg":  "#f9f9f9",
        "btn_bg":    "#c0392b",
        "btn_text":  "#ffffff",
        "table_header_bg":  "#dce6ff",
        "table_header_txt": "#1a1a2e",
        "table_cell_bg":    "#ffffff",
        "table_cell_txt":   "#1a1a2e",
        "table_border":     "#d0d0d0",
    },
}
DEFAULT_THEME = "dark"

# ── Initialize the default active detector (CLD, original behaviour) ──────────
select_detector("CLD")

"""data_io.py — EDM4HEP ROOT file reader and event data extraction.

Two loading modes (set in UI):
  'lazy'  — re-opens file and iterates to requested event index on demand.
  'eager' — pre-loads all events into memory at file open time.
"""

from __future__ import annotations

import math
import sys
import os
from dataclasses import dataclass, field, asdict
from typing import Optional

# Allow running from EventDisplay/ or project root
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

try:
    from podio import root_io as _podio_root_io
    _HAS_PODIO = True
except ImportError:
    _HAS_PODIO = False

# Detector-dependent collection names are read live from `config` (via the
# `config.X` attribute) so that config.select_detector() takes effect even
# though this module is imported once. Only truly invariant filters are bound.
import config
from config import NEUTRINO_PDGS, VALID_GEN_STATUS

try:
    from modules.NeutralRecover import omega_to_pt as _omega_to_pt
    _HAS_OMEGA = True
except Exception:
    _HAS_OMEGA = False


# ── Data classes ──────────────────────────────────────────────────────────────

@dataclass
class HitData:
    group:      str
    collection: str
    x: float; y: float; z: float
    energy: float
    gen_idx:    Optional[int]
    pfo_idx:    Optional[int]
    # True depositing MCParticle PDG (may be a Geant4 secondary), plus its
    # primary generator ancestor (idx + PDG). Filled for hits with gen_idx set.
    gen_pdg:      Optional[int] = None
    primary_idx:  Optional[int] = None
    primary_pdg:  Optional[int] = None


@dataclass
class TrackData:
    pfo_idx: int
    ox: float; oy: float; oz: float    # origin
    px: float; py: float; pz: float    # momentum (GeV)


@dataclass
class GenData:
    idx: int; pdg: int; status: int
    energy: float; p: float; theta: float; phi: float


@dataclass
class PFOData:
    idx: int; pdg: int; charge: float
    energy: float; p: float; theta: float; phi: float


@dataclass
class TruthLinkRaw:
    gen_idx:    int
    pfo_idx:    int
    weight:     float
    gen_pdg:    int
    pfo_pdg:    int
    gen_energy: float
    pfo_energy: float
    gen_px: float; gen_py: float; gen_pz: float
    pfo_px: float; pfo_py: float; pfo_pz: float


@dataclass
class EventData:
    global_idx:  int
    file_path:   str
    local_idx:   int
    gen_particles:  list = field(default_factory=list)   # list[GenData]
    pfos:           list = field(default_factory=list)   # list[PFOData]
    hits:           list = field(default_factory=list)   # list[HitData]
    tracks:         list = field(default_factory=list)   # list[TrackData]
    truth_links:    list = field(default_factory=list)   # list[TruthLinkRaw]


# ── Helpers ───────────────────────────────────────────────────────────────────

def _momentum_to_angles(px: float, py: float, pz: float):
    p = math.sqrt(px * px + py * py + pz * pz)
    if p < 1e-12:
        return 0.0, 0.0, 0.0
    theta = math.acos(max(-1.0, min(1.0, pz / p)))
    phi   = math.atan2(py, px)
    return p, theta, phi


def _safe(obj, method, default=None):
    try:
        return getattr(obj, method)()
    except Exception:
        return default


def _safe_index(obj, default=-999):
    try:
        return obj.getObjectID().index
    except Exception:
        return default


# ── Per-event extraction ──────────────────────────────────────────────────────

def _extract_gen(event, collections_cfg: dict) -> list:
    out = []
    try:
        mc_parts = event.get(config.MC_COLLECTION)
    except Exception:
        return out
    for part in mc_parts:
        try:
            status = int(part.getGeneratorStatus())
            if status not in VALID_GEN_STATUS:
                continue
            pdg = int(part.getPDG())
            if abs(pdg) in NEUTRINO_PDGS:
                continue
            m = part.getMomentum()
            p, theta, phi = _momentum_to_angles(float(m.x), float(m.y), float(m.z))
            if p <= 0:
                continue
            out.append(GenData(
                idx=_safe_index(part), pdg=pdg, status=status,
                energy=float(part.getEnergy()), p=p, theta=theta, phi=phi,
            ))
        except Exception:
            continue
    return out


def _extract_pfos(event) -> list:
    out = []
    try:
        pfos = event.get(config.PFO_COLLECTION)
    except Exception:
        return out
    for pfo in pfos:
        try:
            m = pfo.getMomentum()
            p, theta, phi = _momentum_to_angles(float(m.x), float(m.y), float(m.z))
            out.append(PFOData(
                idx=_safe_index(pfo),
                pdg=int(pfo.getPDG()),
                charge=float(pfo.getCharge()),
                energy=float(pfo.getEnergy()),
                p=p, theta=theta, phi=phi,
            ))
        except Exception:
            continue
    return out


def _build_hit_maps(event, collections_cfg: dict):
    """Returns (type_map, energy_map, pos_map) keyed by (collID, index)."""
    type_map   = {}
    energy_map = {}
    pos_map    = {}
    for grp, info in collections_cfg.items():
        if grp == "tracker":
            continue  # tracker hits handled separately via tracks
        for coll_name in info["collections"]:
            try:
                coll = event.get(coll_name)
            except Exception:
                continue
            for hit in coll:
                try:
                    oid = hit.getObjectID()
                    key = (oid.collectionID, oid.index)
                    type_map[key]   = (grp, coll_name)
                    energy_map[key] = float(hit.getEnergy())
                    pos = hit.getPosition()
                    pos_map[key]    = (float(pos.x), float(pos.y), float(pos.z))
                except Exception:
                    continue
    return type_map, energy_map, pos_map


def _build_hit_to_gen(event, type_map: dict) -> dict:
    out = {}
    if not config.CALO_LINK_COLLECTION:
        return out  # detector has no CaloHit→MCParticle link (e.g. ILD)
    try:
        links = list(event.get(config.CALO_LINK_COLLECTION))
    except Exception:
        return out
    for link in links:
        try:
            hit_obj = link.getRec()
            gen_obj = link.getSim()
            oid = hit_obj.getObjectID()
            key = (oid.collectionID, oid.index)
            if key in type_map:
                out[key] = _safe_index(gen_obj)
        except Exception:
            continue
    return out


def _build_hit_to_pfo(event, type_map: dict) -> dict:
    out = {}
    try:
        pfos = event.get(config.PFO_COLLECTION)
    except Exception:
        return out
    for pfo in pfos:
        pfo_idx = _safe_index(pfo)
        try:
            for cluster in pfo.getClusters():
                for hit in cluster.getHits():
                    try:
                        oid = hit.getObjectID()
                        key = (oid.collectionID, oid.index)
                        if key in type_map:
                            out[key] = pfo_idx
                    except Exception:
                        continue
        except Exception:
            continue
    return out


def _extract_hits(event, collections_cfg: dict) -> list:
    type_map, energy_map, pos_map = _build_hit_maps(event, collections_cfg)
    hit_to_gen = _build_hit_to_gen(event, type_map)
    hit_to_pfo = _build_hit_to_pfo(event, type_map)
    out = []
    for key, (grp, coll) in type_map.items():
        x, y, z = pos_map[key]
        out.append(HitData(
            group=grp, collection=coll,
            x=x, y=y, z=z,
            energy=energy_map.get(key, 0.0),
            gen_idx=hit_to_gen.get(key),
            pfo_idx=hit_to_pfo.get(key),
        ))
    return out


def _extract_tracks(event) -> list:
    out = []
    try:
        pfos = event.get(config.PFO_COLLECTION)
    except Exception:
        return out
    for pfo in pfos:
        pfo_idx = _safe_index(pfo)
        try:
            for track in pfo.getTracks():
                states = track.getTrackStates()
                if len(states) == 0:
                    continue
                st = states[0]
                if _HAS_OMEGA:
                    pt = _omega_to_pt(float(st.omega), isclic=False)
                else:
                    pt = abs(1.0 / float(st.omega)) * 0.3 * 2.0 if abs(st.omega) > 1e-9 else 0.0
                px = pt * math.cos(float(st.phi))
                py = pt * math.sin(float(st.phi))
                pz = float(st.tanLambda) * pt
                out.append(TrackData(
                    pfo_idx=pfo_idx,
                    ox=0.0, oy=0.0, oz=0.0,
                    px=px, py=py, pz=pz,
                ))
        except Exception:
            continue
    return out


def _extract_sim_tracker_hits(event) -> list:
    """Extract SimTrackerHit collections; each hit carries a direct MCParticle link."""
    out = []
    for coll_name in config.SIM_TRACKER_COLLECTIONS:
        try:
            coll = event.get(coll_name)
        except Exception:
            continue
        for hit in coll:
            try:
                pos = hit.getPosition()
                edep = float(hit.getEDep())
                # Try EDM4HEP ≥ 2 API first, fall back to older name
                try:
                    mc_part = hit.getMCParticle()
                except AttributeError:
                    mc_part = hit.getParticle()
                gen_idx = _safe_index(mc_part) if mc_part is not None else None
                out.append(HitData(
                    group="tracker", collection=coll_name,
                    x=float(pos.x), y=float(pos.y), z=float(pos.z),
                    energy=edep,
                    gen_idx=gen_idx,
                    pfo_idx=None,
                ))
            except Exception:
                continue
    return out


def _extract_truth_links(event) -> list:
    """Try known RecoMCTruthLink collection variants; return raw link records."""
    out = []
    candidates = list(config.TRUTH_LINK_VARIANTS)
    try:
        available = list(event.getAvailableCollections())
        for name in available:
            if "RecoMCTruthLink" in name and name not in candidates:
                candidates.append(name)
    except Exception:
        pass

    links_raw = []
    for coll_name in candidates:
        try:
            links_raw = list(event.get(coll_name))
            if links_raw:
                break
        except Exception:
            continue

    for link in links_raw:
        try:
            if hasattr(link, "getSim") and hasattr(link, "getRec"):
                sim_obj = link.getSim()
                rec_obj = link.getRec()
            else:
                sim_obj = link.getTo()
                rec_obj = link.getFrom()

            gen_status = _safe(sim_obj, "getGeneratorStatus", -1)
            if gen_status not in VALID_GEN_STATUS:
                continue
            gen_pdg = int(_safe(sim_obj, "getPDG", 0))
            if abs(gen_pdg) in NEUTRINO_PDGS:
                continue

            gen_m = _safe(sim_obj, "getMomentum")
            pfo_m = _safe(rec_obj, "getMomentum")
            out.append(TruthLinkRaw(
                gen_idx    =_safe_index(sim_obj),
                pfo_idx    =_safe_index(rec_obj),
                weight     =float(_safe(link, "getWeight", 0.0)),
                gen_pdg    =gen_pdg,
                pfo_pdg    =int(_safe(rec_obj, "getPDG", -999)),
                gen_energy =float(_safe(sim_obj, "getEnergy", float("nan"))),
                pfo_energy =float(_safe(rec_obj, "getEnergy", float("nan"))),
                gen_px     =float(gen_m.x) if gen_m is not None else float("nan"),
                gen_py     =float(gen_m.y) if gen_m is not None else float("nan"),
                gen_pz     =float(gen_m.z) if gen_m is not None else float("nan"),
                pfo_px     =float(pfo_m.x) if pfo_m is not None else float("nan"),
                pfo_py     =float(pfo_m.y) if pfo_m is not None else float("nan"),
                pfo_pz     =float(pfo_m.z) if pfo_m is not None else float("nan"),
            ))
        except Exception:
            continue
    return out


def _build_mc_meta(event) -> dict:
    """Return {mc_idx: (pdg, status, primary_idx, primary_pdg)} for all MCParticles.

    The primary ancestor is found by walking getParents() until a generator
    particle (generatorStatus==1 or not createdInSimulation) is reached; this
    lets the display attribute a Geant4 secondary (delta ray, decay-in-flight)
    back to the physics particle it came from.
    """
    meta: dict = {}
    try:
        mc = list(event.get(config.MC_COLLECTION))
    except Exception:
        return meta

    def _primary_of(part):
        seen = set()
        cur = part
        while True:
            idx = _safe_index(cur)
            if idx in seen:      # cycle guard
                return cur
            seen.add(idx)
            try:
                is_primary = (int(cur.getGeneratorStatus()) == 1
                              or not bool(cur.isCreatedInSimulation()))
            except Exception:
                is_primary = True
            if is_primary:
                return cur
            try:
                parents = cur.getParents()
            except Exception:
                return cur
            if len(parents) == 0:
                return cur
            cur = parents[0]     # follow the primary lineage

    for part in mc:
        idx = _safe_index(part)
        try:
            pdg    = int(part.getPDG())
            status = int(part.getGeneratorStatus())
        except Exception:
            pdg, status = 0, -1
        anc = _primary_of(part)
        meta[idx] = (pdg, status, _safe_index(anc), int(_safe(anc, "getPDG", 0)))
    return meta


def _extract_event(event, global_idx: int, file_path: str, local_idx: int,
                   collections_cfg: dict) -> EventData:
    hits = _extract_hits(event, collections_cfg)
    hits += _extract_sim_tracker_hits(event)

    # Annotate each hit's true MCParticle with PDG + primary ancestor
    mc_meta = _build_mc_meta(event)
    for h in hits:
        info = mc_meta.get(h.gen_idx) if h.gen_idx is not None else None
        if info is not None:
            h.gen_pdg, _status, h.primary_idx, h.primary_pdg = info

    return EventData(
        global_idx    =global_idx,
        file_path     =file_path,
        local_idx     =local_idx,
        gen_particles =_extract_gen(event, collections_cfg),
        pfos          =_extract_pfos(event),
        hits          =hits,
        tracks        =_extract_tracks(event),
        truth_links   =_extract_truth_links(event),
    )


# ── EventReader ───────────────────────────────────────────────────────────────

class EventReader:
    """Reads EDM4HEP ROOT files. Supports eager and lazy loading modes."""

    def __init__(self):
        self.filepath:    str  = ""
        self.collections: dict = config.DETECTOR_GROUPS   # current config
        self._cache:      list = []                 # list[EventData] for eager
        self._n_events:   int  = 0
        self._mode:       str  = "lazy"             # "eager"|"lazy"

    # ── Public API ────────────────────────────────────────────────────────────

    def load(self, filepath: str, mode: str = "lazy",
             collections_override: dict | None = None) -> int:
        """Open file and pre-load (eager) or index (lazy). Returns event count."""
        if not _HAS_PODIO:
            raise RuntimeError("podio not available — run inside Key4HEP environment.")
        self.filepath = filepath
        self._mode = mode
        self._cache = []
        self.collections = collections_override or config.DETECTOR_GROUPS

        # Count events (fast scan)
        reader = _podio_root_io.Reader(filepath)
        events = list(reader.get("events"))  # materialise to know length
        self._n_events = len(events)

        if mode == "eager":
            for local_idx, ev in enumerate(events):
                self._cache.append(
                    _extract_event(ev, local_idx, filepath, local_idx, self.collections)
                )

        return self._n_events

    def get_event(self, idx: int) -> EventData:
        """Return EventData for event at index idx."""
        if not self.filepath:
            raise RuntimeError("No file loaded. Call load() first.")
        idx = max(0, min(idx, self._n_events - 1))

        if self._mode == "eager" and self._cache:
            return self._cache[idx]

        # Lazy: iterate to idx
        reader = _podio_root_io.Reader(self.filepath)
        for i, ev in enumerate(reader.get("events")):
            if i == idx:
                return _extract_event(ev, idx, self.filepath, idx, self.collections)
        raise IndexError(f"Event {idx} not found in file.")

    @property
    def n_events(self) -> int:
        return self._n_events


# ── Serialization for dcc.Store ───────────────────────────────────────────────

def event_to_dict(ev: EventData) -> dict:
    return asdict(ev)


def event_from_dict(d: dict) -> EventData:
    ev = EventData(
        global_idx=d["global_idx"],
        file_path=d["file_path"],
        local_idx=d["local_idx"],
    )
    ev.gen_particles = [GenData(**g) for g in d.get("gen_particles", [])]
    ev.pfos          = [PFOData(**p) for p in d.get("pfos", [])]
    ev.hits          = [HitData(**h) for h in d.get("hits", [])]
    ev.tracks        = [TrackData(**t) for t in d.get("tracks", [])]
    ev.truth_links   = [TruthLinkRaw(**l) for l in d.get("truth_links", [])]
    return ev

"""geometry.py — Load detector geometry and build Plotly wireframe traces.

The active detector (CLD, ILD, …) is selected via config.select_detector();
this module reads the active geometry path/fallback from `config`.

Sources (in priority order):
  1. Explicit xml_path argument (can be set via UI)
  2. k4geo_DIR env var + config.GEOMETRY_RELPATH
  3. config.GEOMETRY_FALLBACK hardcoded values
"""

from __future__ import annotations

import math
import os
import sys
import xml.etree.ElementTree as ET
from typing import Optional

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np

try:
    import plotly.graph_objects as go
    _HAS_PLOTLY = True
except ImportError:
    _HAS_PLOTLY = False

import config


# ── Path resolution ───────────────────────────────────────────────────────────

def default_geometry_path() -> Optional[str]:
    k4geo = os.environ.get("k4geo_DIR", "")
    if not k4geo or not config.GEOMETRY_RELPATH:
        return None
    path = os.path.join(k4geo, config.GEOMETRY_RELPATH)
    return path if os.path.isfile(path) else None


# ── XML parsing ───────────────────────────────────────────────────────────────

def _float_attr(elem, key, default=None):
    val = elem.get(key)
    if val is None:
        return default
    val = val.strip().rstrip("*mm").strip()
    try:
        return float(val.split("*")[0])
    except ValueError:
        return default


def _find_constant(root, name: str) -> Optional[float]:
    for c in root.iter("constant"):
        if c.get("name") == name:
            val = c.get("value", "")
            try:
                return float(val.split("*")[0])
            except ValueError:
                pass
    return None


def load_detector_geometry(xml_path: str) -> dict:
    """
    Parse a DD4HEP compact XML file and extract key barrel/endcap dimensions.

    Missing volumes are filled from the active detector's geometry fallback.
    Returns a dict like config.GEOMETRY_FALLBACK, or empty dict on failure.

    Note: many detectors (e.g. ILD) define envelopes via included files and
    DD4hep constants rather than inline <dimensions>, so this parser usually
    yields nothing for them and the fallback dominates.
    """
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
    except Exception:
        return {}

    geo = {}

    # Strategy: look for detector elements by name patterns
    for det in root.iter("detector"):
        name = det.get("name", "")
        dtype = det.get("type", "")

        # ECAL barrel
        if "ECalBarrel" in name or ("ecal" in name.lower() and "barrel" in name.lower()):
            dims = det.find(".//dimensions") or det.find(".//barrel_envelope")
            if dims is not None:
                rmin   = _float_attr(dims, "inner_r") or _float_attr(dims, "rmin")
                rmax   = _float_attr(dims, "outer_r") or _float_attr(dims, "rmax")
                zhalf  = _float_attr(dims, "zhalf") or _float_attr(dims, "z_max")
                if rmin and rmax and zhalf:
                    geo["ecal_barrel"] = {"rmin": rmin, "rmax": rmax, "zhalf": zhalf}

        # ECAL endcap
        elif "ECalEndcap" in name or ("ecal" in name.lower() and "endcap" in name.lower()):
            dims = det.find(".//dimensions") or det.find(".//endcap_envelope")
            if dims is not None:
                rmin  = _float_attr(dims, "inner_r") or _float_attr(dims, "rmin")
                rmax  = _float_attr(dims, "outer_r") or _float_attr(dims, "rmax")
                zpos  = _float_attr(dims, "zmin") or _float_attr(dims, "z_offset")
                zthk  = _float_attr(dims, "zhalf") or _float_attr(dims, "dz")
                if rmin and rmax and zpos:
                    geo["ecal_endcap"] = {"rmin": rmin, "rmax": rmax,
                                          "zpos": zpos, "zthick": zthk or 0}

        # HCAL barrel
        elif "HCalBarrel" in name or ("hcal" in name.lower() and "barrel" in name.lower()):
            dims = det.find(".//dimensions") or det.find(".//barrel_envelope")
            if dims is not None:
                rmin  = _float_attr(dims, "inner_r") or _float_attr(dims, "rmin")
                rmax  = _float_attr(dims, "outer_r") or _float_attr(dims, "rmax")
                zhalf = _float_attr(dims, "zhalf") or _float_attr(dims, "z_max")
                if rmin and rmax and zhalf:
                    geo["hcal_barrel"] = {"rmin": rmin, "rmax": rmax, "zhalf": zhalf}

        # HCAL endcap
        elif "HCalEndcap" in name or ("hcal" in name.lower() and "endcap" in name.lower()):
            dims = det.find(".//dimensions") or det.find(".//endcap_envelope")
            if dims is not None:
                rmin  = _float_attr(dims, "inner_r") or _float_attr(dims, "rmin")
                rmax  = _float_attr(dims, "outer_r") or _float_attr(dims, "rmax")
                zpos  = _float_attr(dims, "zmin") or _float_attr(dims, "z_offset")
                zthk  = _float_attr(dims, "zhalf") or _float_attr(dims, "dz")
                if rmin and rmax and zpos:
                    geo["hcal_endcap"] = {"rmin": rmin, "rmax": rmax,
                                          "zpos": zpos, "zthick": zthk or 0}

        # Muon barrel
        elif "Muon" in name and "barrel" in name.lower():
            dims = det.find(".//dimensions")
            if dims is not None:
                rmin  = _float_attr(dims, "inner_r") or _float_attr(dims, "rmin")
                rmax  = _float_attr(dims, "outer_r") or _float_attr(dims, "rmax")
                zhalf = _float_attr(dims, "zhalf") or _float_attr(dims, "z_max")
                if rmin and rmax and zhalf:
                    geo["muon_barrel"] = {"rmin": rmin, "rmax": rmax, "zhalf": zhalf}

    # Fill missing keys from the active detector's fallback
    for key, val in config.GEOMETRY_FALLBACK.items():
        if key not in geo:
            geo[key] = val

    return geo


# Backwards-compatible alias
load_cld_geometry = load_detector_geometry


def get_geometry(xml_path: Optional[str] = None) -> dict:
    """Return geometry dict, trying xml_path, then env, then fallback."""
    if xml_path and os.path.isfile(xml_path):
        geo = load_detector_geometry(xml_path)
        if geo:
            return geo
    default = default_geometry_path()
    if default:
        geo = load_detector_geometry(default)
        if geo:
            return geo
    return dict(config.GEOMETRY_FALLBACK)


# ── Plotly wireframe builders ─────────────────────────────────────────────────

def _circle_xy(r: float, z: float, n: int = 64):
    """Return (x, y, z_arr, None) arrays for a circle at height z."""
    angles = np.linspace(0, 2 * math.pi, n + 1)
    x = r * np.cos(angles)
    y = r * np.sin(angles)
    z_arr = np.full(n + 1, z)
    return x.tolist(), y.tolist(), z_arr.tolist()


def _cylinder_traces(rmin: float, rmax: float, zhalf: float,
                     color: str, name: str, n: int = 64) -> list:
    """Barrel wireframe: inner & outer circles at ±z, plus 4 connecting lines."""
    traces = []
    for r, label in [(rmin, "inner"), (rmax, "outer")]:
        for z in [zhalf, -zhalf]:
            x, y, z_arr = _circle_xy(r, z, n)
            traces.append(go.Scatter3d(
                x=x, y=y, z=z_arr,
                mode="lines",
                line=dict(color=color, width=1),
                name=f"{name} {label} z={'+'if z>0 else '-'}",
                legendgroup=name,
                showlegend=(label == "outer" and z > 0),
                hoverinfo="skip",
                opacity=0.5,
            ))
    # 4 vertical lines connecting ±z at rmin and rmax
    for r in [rmin, rmax]:
        for angle in [0, math.pi / 2, math.pi, 3 * math.pi / 2]:
            xv = r * math.cos(angle)
            yv = r * math.sin(angle)
            traces.append(go.Scatter3d(
                x=[xv, xv], y=[yv, yv], z=[zhalf, -zhalf],
                mode="lines",
                line=dict(color=color, width=1),
                name=name, legendgroup=name, showlegend=False,
                hoverinfo="skip", opacity=0.4,
            ))
    return traces


def _annulus_traces(rmin: float, rmax: float, zpos: float,
                    color: str, name: str, n: int = 64) -> list:
    """Endcap wireframe: two circles at ±zpos."""
    traces = []
    for z_sign in [1, -1]:
        z = z_sign * zpos
        for r in [rmin, rmax]:
            x, y, z_arr = _circle_xy(r, z, n)
            traces.append(go.Scatter3d(
                x=x, y=y, z=z_arr,
                mode="lines",
                line=dict(color=color, width=1),
                name=f"{name} endcap z={'+'if z>0 else '-'}",
                legendgroup=name,
                showlegend=(r == rmax and z_sign == 1),
                hoverinfo="skip",
                opacity=0.4,
            ))
        # 4 radial spokes
        for angle in [0, math.pi / 2, math.pi, 3 * math.pi / 2]:
            traces.append(go.Scatter3d(
                x=[rmin * math.cos(angle), rmax * math.cos(angle)],
                y=[rmin * math.sin(angle), rmax * math.sin(angle)],
                z=[z, z],
                mode="lines",
                line=dict(color=color, width=1),
                name=name, legendgroup=name, showlegend=False,
                hoverinfo="skip", opacity=0.3,
            ))
    return traces


def build_geometry_traces(geometry: dict) -> list:
    """Return list of Plotly Scatter3d traces for all detector volumes."""
    if not _HAS_PLOTLY or not geometry:
        return []

    group_colors = {grp: info["color"] for grp, info in config.DETECTOR_GROUPS.items()}
    traces = []

    mapping = [
        ("tpc_barrel",   "tracker", "TPC",        "barrel"),
        ("ecal_barrel",  "ecal",  "ECAL Barrel",  "barrel"),
        ("ecal_endcap",  "ecal",  "ECAL Endcap",  "endcap"),
        ("hcal_barrel",  "hcal",  "HCAL Barrel",  "barrel"),
        ("hcal_endcap",  "hcal",  "HCAL Endcap",  "endcap"),
        ("muon_barrel",  "muon",  "Muon Barrel",  "barrel"),
    ]
    for key, grp, label, shape in mapping:
        g = geometry.get(key)
        if not g:
            continue
        color = group_colors.get(grp, "#aaaaaa")
        if shape == "barrel":
            traces += _cylinder_traces(g["rmin"], g["rmax"], g["zhalf"], color, label)
        else:
            traces += _annulus_traces(g["rmin"], g["rmax"], g["zpos"], color, label)

    return traces

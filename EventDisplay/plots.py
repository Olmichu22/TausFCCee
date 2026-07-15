"""plots.py — Plotly figure builders for the FCC Event Display."""

from __future__ import annotations

import math
import sys
import os
from collections import defaultdict
from typing import Optional

import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from config import DETECTOR_GROUPS, THEMES, PDG_NAMES

_TAB20 = px.colors.qualitative.Dark24 + px.colors.qualitative.Light24

# Consistent color per PDG name
_PDG_COLOR_MAP: dict = {}
def _pdg_color(name: str) -> str:
    if name not in _PDG_COLOR_MAP:
        _PDG_COLOR_MAP[name] = _TAB20[len(_PDG_COLOR_MAP) % len(_TAB20)]
    return _PDG_COLOR_MAP[name]


def _theme(theme: str) -> dict:
    return THEMES.get(theme, THEMES["dark"])


def _axis3d(theme: str) -> dict:
    t = _theme(theme)
    return dict(
        backgroundcolor=t["plot_bg"],
        gridcolor=t["grid"],
        showbackground=True,
        tickfont=dict(size=9, color=t["text"]),
        title_font=dict(size=10, color=t["text"]),
    )


def _layout_base(theme: str, height: int = 420) -> dict:
    t = _theme(theme)
    return dict(
        paper_bgcolor=t["app_bg"],
        plot_bgcolor=t["plot_bg"],
        font=dict(color=t["text"], size=11),
        margin=dict(l=50, r=20, t=50, b=50),
        height=height,
        legend=dict(
            font=dict(size=10, color=t["text"]),
            bgcolor=t["card_bg"],
            bordercolor=t["border"],
            borderwidth=1,
        ),
    )


def _pdg_name(pdg) -> str:
    if pdg is None:
        return "?"
    return PDG_NAMES.get(int(pdg), str(pdg))


# ── 3D Event Display ──────────────────────────────────────────────────────────

def make_3d_figure(
    hits,
    tracks,
    gen_particles,
    pfos,
    geometry_traces: list,
    *,
    active_groups: set,
    color_mode: str = "detector",
    energy_range: tuple = (0.0, 1e9),
    gen_filter: Optional[list] = None,
    pfo_filter: Optional[list] = None,
    show_geometry: bool = False,
    show_tracks: bool = True,
    hide_secondary_hits: bool = False,
    momentum_scale: float = 60.0,
    pfo_to_gen: dict = {},
    theme: str = "dark",
) -> go.Figure:
    """
    Build the main 3D event display figure.

    color_mode: 'detector' → color by detector group
                'gen'      → color by associated gen particle index (tab20)
                'pfo'      → color by associated PFO index (tab20)
    pfo_to_gen: mapping pfo_idx → gen_idx (for track filtering by gen)
    """
    t = _theme(theme)
    e_min, e_max = float(energy_range[0]), float(energy_range[1])

    gen_idx_set = set(gen_filter) if gen_filter else None
    pfo_idx_set = set(pfo_filter) if pfo_filter else None

    gen_colors = {g.idx: _TAB20[i % len(_TAB20)] for i, g in enumerate(gen_particles)}
    pfo_colors = {p.idx: _TAB20[i % len(_TAB20)] for i, p in enumerate(pfos)}
    pfo_pdg = {p.idx: p.pdg for p in pfos}   # for hover PDG labels

    buckets: dict = defaultdict(lambda: {"x": [], "y": [], "z": [], "text": []})

    def _owner_key(hit):
        if color_mode == "gen":
            return ("gen", hit.gen_idx)
        elif color_mode == "pfo":
            return ("pfo", hit.pfo_idx)
        else:
            return ("det", hit.group)

    for hit in hits:
        if hit.group not in active_groups:
            continue
        if hit.energy < e_min or hit.energy > e_max:
            continue
        if gen_idx_set is not None and hit.gen_idx not in gen_idx_set:
            continue
        # Geant4 secondaries: depositing MCParticle differs from its primary
        # ancestor (delta rays, decays in flight). Hits without truth info
        # (e.g. ILD calo hits) are kept.
        if (hide_secondary_hits and hit.gen_idx is not None
                and hit.primary_idx is not None
                and hit.primary_idx != hit.gen_idx):
            continue
        # Hits with pfo_idx=None (e.g. SimTrackerHits) bypass the PFO filter
        if pfo_idx_set is not None and hit.pfo_idx is not None and hit.pfo_idx not in pfo_idx_set:
            continue

        key = (hit.group, _owner_key(hit))
        buckets[key]["x"].append(hit.x)
        buckets[key]["y"].append(hit.y)
        buckets[key]["z"].append(hit.z)
        # Sim (true MCParticle) line: id + PDG, and primary ancestor if secondary
        if hit.gen_idx is not None:
            sim_line = f"sim {hit.gen_idx} ({_pdg_name(hit.gen_pdg)})"
            if hit.primary_idx is not None and hit.primary_idx != hit.gen_idx:
                sim_line += f" ← primaria {hit.primary_idx} ({_pdg_name(hit.primary_pdg)})"
        else:
            sim_line = "sim None"
        # PFO line: id + PDG
        if hit.pfo_idx is not None:
            pfo_line = f"pfo {hit.pfo_idx} ({_pdg_name(pfo_pdg.get(hit.pfo_idx))})"
        else:
            pfo_line = "pfo None"
        buckets[key]["text"].append(
            f"{hit.group}|{hit.collection}<br>E={hit.energy:.4f} GeV"
            f"<br>{sim_line}<br>{pfo_line}"
        )

    traces = []
    for (grp, owner_info), pts in buckets.items():
        if not pts["x"]:
            continue

        mode_type, owner_val = owner_info
        if mode_type == "gen":
            color = gen_colors.get(owner_val, "#aaaaaa") if owner_val is not None else "#555555"
            name_label = _pdg_name(next((g.pdg for g in gen_particles if g.idx == owner_val), None))
            trace_name = f"Gen {owner_val} ({name_label}) [{grp}]" if owner_val is not None else f"{grp} unassoc."
        elif mode_type == "pfo":
            color = pfo_colors.get(owner_val, "#aaaaaa") if owner_val is not None else "#555555"
            name_label = _pdg_name(next((p.pdg for p in pfos if p.idx == owner_val), None))
            trace_name = f"PFO {owner_val} ({name_label}) [{grp}]" if owner_val is not None else f"{grp} unassoc."
        else:
            color = DETECTOR_GROUPS.get(grp, {}).get("color", "#aaaaaa")
            trace_name = grp.upper()

        traces.append(go.Scatter3d(
            x=pts["x"], y=pts["y"], z=pts["z"],
            mode="markers",
            marker=dict(size=3, color=color, opacity=0.85, line=dict(width=0)),
            name=trace_name,
            text=pts["text"],
            hovertemplate="%{text}<extra></extra>",
        ))

    # ── Track vectors ────────────────────────────────────────────────────────
    # Tracks always come from PFOs (pfo.getTracks()). A gen particle can only
    # "own" a track indirectly: via the pfo_to_gen map.
    # If a PFO has no gen match (pfo_to_gen has no entry for it), its tracks
    # are shown as "unmatched" (gray) rather than hidden — so they remain
    # visible even when a gen filter is active.
    if show_tracks and tracks:
        tracker_color = DETECTOR_GROUPS.get("tracker", {}).get("color", "#9b59b6")

        if color_mode == "detector":
            # Single "TRACKER" trace — all track segments combined
            tx, ty, tz = [], [], []
            for tv in tracks:
                if gen_idx_set is not None:
                    tv_gen = pfo_to_gen.get(tv.pfo_idx)
                    if tv_gen is not None and tv_gen not in gen_idx_set:
                        continue
                if pfo_idx_set is not None and tv.pfo_idx not in pfo_idx_set:
                    continue
                ex = tv.ox + tv.px * momentum_scale
                ey = tv.oy + tv.py * momentum_scale
                ez = tv.oz + tv.pz * momentum_scale
                tx += [tv.ox, ex, None]
                ty += [tv.oy, ey, None]
                tz += [tv.oz, ez, None]
            if tx:
                traces.append(go.Scatter3d(
                    x=tx, y=ty, z=tz,
                    mode="lines",
                    line=dict(color=tracker_color, width=3),
                    name="TRACKER",
                    legendgroup="TRACKER",
                    hoverinfo="skip",
                    opacity=0.9,
                ))

        else:
            # Group tracks by owner (gen or pfo) — same color as the hits
            track_buckets: dict = defaultdict(lambda: {"x": [], "y": [], "z": []})
            for tv in tracks:
                tv_gen = pfo_to_gen.get(tv.pfo_idx)  # None if PFO unmatched

                if color_mode == "gen":
                    owner = tv_gen  # None = unmatched PFO
                    if gen_idx_set is not None:
                        # Hide only if gen is known and NOT in the filter
                        if owner is not None and owner not in gen_idx_set:
                            continue
                    color_key = ("gen", owner)
                else:  # pfo
                    owner = tv.pfo_idx
                    if pfo_idx_set is not None and owner not in pfo_idx_set:
                        continue
                    color_key = ("pfo", owner)

                ex = tv.ox + tv.px * momentum_scale
                ey = tv.oy + tv.py * momentum_scale
                ez = tv.oz + tv.pz * momentum_scale
                track_buckets[color_key]["x"] += [tv.ox, ex, None]
                track_buckets[color_key]["y"] += [tv.oy, ey, None]
                track_buckets[color_key]["z"] += [tv.oz, ez, None]

            for (mode_t, owner_val), pts in track_buckets.items():
                if not pts["x"]:
                    continue
                if mode_t == "gen":
                    if owner_val is None:
                        color = "#777777"
                        tname = "Track (unmatched gen)"
                    else:
                        color = gen_colors.get(owner_val, tracker_color)
                        plabel = _pdg_name(next(
                            (g.pdg for g in gen_particles if g.idx == owner_val), None))
                        tname = f"Track gen {owner_val} ({plabel})"
                else:  # pfo
                    color = pfo_colors.get(owner_val, tracker_color)
                    plabel = _pdg_name(next(
                        (p.pdg for p in pfos if p.idx == owner_val), None))
                    gen_info = f"→gen {pfo_to_gen[owner_val]}" if owner_val in pfo_to_gen else ""
                    tname = f"Track PFO {owner_val} ({plabel}) {gen_info}".strip()

                traces.append(go.Scatter3d(
                    x=pts["x"], y=pts["y"], z=pts["z"],
                    mode="lines",
                    line=dict(color=color, width=3),
                    name=tname,
                    hoverinfo="skip",
                    opacity=0.9,
                ))

    if show_geometry:
        for tr in geometry_traces:
            traces.append(tr)

    layout = go.Layout(
        paper_bgcolor=t["app_bg"],
        scene=dict(
            xaxis=dict(title="x [mm]", **_axis3d(theme)),
            yaxis=dict(title="y [mm]", **_axis3d(theme)),
            zaxis=dict(title="z [mm]", **_axis3d(theme)),
            bgcolor=t["plot_bg"],
            aspectmode="data",
            # Orient the beam/z-axis horizontally (detector lying on its side):
            # "up" = y, and the eye sits mostly in the x-y plane so z runs across
            # the screen instead of pointing up. The user can still rotate freely.
            camera=dict(
                up=dict(x=0, y=1, z=0),
                eye=dict(x=1.5, y=1.0, z=0.4),
                center=dict(x=0, y=0, z=0),
            ),
        ),
        margin=dict(l=0, r=0, t=30, b=0),
        height=580,
        legend=dict(
            font=dict(size=10, color=t["text"]),
            bgcolor=t["card_bg"],
            bordercolor=t["border"],
            borderwidth=1,
        ),
        uirevision="3d-camera",
    )
    return go.Figure(
        data=traces if traces else [go.Scatter3d(x=[], y=[], z=[], mode="markers")],
        layout=layout,
    )


# ── Histogram figures ─────────────────────────────────────────────────────────

_HIST_OPTIONS = [
    {"label": "Hit energy by detector group",    "value": "hit_energy_by_detector"},
    {"label": "Gen particle energy by PDG",      "value": "gen_energy"},
    {"label": "PFO energy by PDG",               "value": "pfo_energy"},
    {"label": "Gen PDG distribution",            "value": "pdg_gen"},
    {"label": "PFO PDG distribution",            "value": "pdg_pfo"},
    {"label": "Hits per PFO",                    "value": "hits_per_pfo"},
    {"label": "Hits per PFO by detector group",  "value": "hits_per_pfo_by_group"},
    {"label": "Hit multiplicity by group",       "value": "hit_mult_by_group"},
    {"label": "Gen vs PFO energy scatter",       "value": "gen_vs_pfo_energy"},
]


def make_histogram_figure(event_data, hist_type: str, theme: str = "dark") -> go.Figure:
    t = _theme(theme)
    lb = _layout_base(theme)

    # ── Hit energy by detector ────────────────────────────────────────────────
    if hist_type == "hit_energy_by_detector":
        e_all = [h.energy for h in event_data.hits if h.energy > 1e-9]
        if not e_all:
            return _empty_fig("No hit energy data > 0", theme)

        # Pre-compute log10 values — Plotly histogram bins work correctly in
        # log-space this way; using type="log" on a regular histogram causes
        # binning artifacts.
        log_all = [math.log10(e) for e in e_all]
        x_min = min(log_all) - 0.2
        x_max = max(log_all) + 0.2

        # Build custom tick labels at integer/half-integer decades
        tick_vals, tick_text = [], []
        for exp in range(math.floor(x_min), math.ceil(x_max) + 1):
            for sub in [0, 0.5]:
                v = exp + sub
                if x_min <= v <= x_max:
                    tick_vals.append(v)
                    tick_text.append(f"10^{v:.1g}" if sub else f"10^{exp}")

        htraces = []
        for grp, info in DETECTOR_GROUPS.items():
            log_vals = [math.log10(h.energy)
                        for h in event_data.hits
                        if h.group == grp and h.energy > 1e-9]
            if log_vals:
                htraces.append(go.Histogram(
                    x=log_vals, name=grp.upper(),
                    marker_color=info["color"], opacity=0.75,
                    nbinsx=60,
                    xbins=dict(start=x_min, end=x_max,
                               size=(x_max - x_min) / 60),
                ))
        fig = go.Figure(data=htraces, layout=go.Layout(
            barmode="overlay",
            xaxis=dict(title="Hit Energy [GeV]",
                       tickmode="array",
                       tickvals=tick_vals,
                       ticktext=tick_text,
                       range=[x_min, x_max]),
            yaxis=dict(title="N hits"),
            title="Hit Energy by Detector Group",
            **lb,
        ))

    # ── Gen energy subplots by PDG ────────────────────────────────────────────
    elif hist_type == "gen_energy":
        pdg_groups: dict = defaultdict(list)
        for g in event_data.gen_particles:
            pdg_groups[_pdg_name(g.pdg)].append(g.energy)

        if not pdg_groups:
            return _empty_fig("No gen particles", theme)

        n_pdg = len(pdg_groups)
        n_cols = min(n_pdg, 3)
        n_rows = math.ceil(n_pdg / n_cols)
        fig = make_subplots(
            rows=n_rows, cols=n_cols,
            subplot_titles=list(pdg_groups.keys()),
            shared_xaxes=False, shared_yaxes=False,
            horizontal_spacing=0.08, vertical_spacing=0.12,
        )
        for i, (name, energies) in enumerate(pdg_groups.items()):
            r, c = divmod(i, n_cols)
            fig.add_trace(
                go.Histogram(x=energies, name=name,
                             marker_color=_pdg_color(name), opacity=0.8,
                             nbinsx=max(5, len(energies)),
                             showlegend=True),
                row=r + 1, col=c + 1,
            )
        fig.update_layout(
            title="Gen Particle Energy by PDG",
            height=max(350, n_rows * 220),
            **{k: v for k, v in lb.items() if k != "height"},
        )
        fig.update_annotations(font_color=t["text"])
        fig.update_xaxes(
            title_text="Energy [GeV]",
            gridcolor=t["grid"],
            tickfont=dict(color=t["text"]),
            title_font=dict(color=t["text"]),
        )
        fig.update_yaxes(
            title_text="N particles",
            gridcolor=t["grid"],
            tickfont=dict(color=t["text"]),
            title_font=dict(color=t["text"]),
        )

    # ── PFO energy subplots by PDG ────────────────────────────────────────────
    elif hist_type == "pfo_energy":
        pdg_groups: dict = defaultdict(list)
        for p in event_data.pfos:
            pdg_groups[_pdg_name(p.pdg)].append(p.energy)

        if not pdg_groups:
            return _empty_fig("No PFOs", theme)

        n_pdg = len(pdg_groups)
        n_cols = min(n_pdg, 3)
        n_rows = math.ceil(n_pdg / n_cols)
        fig = make_subplots(
            rows=n_rows, cols=n_cols,
            subplot_titles=list(pdg_groups.keys()),
            shared_xaxes=False, shared_yaxes=False,
            horizontal_spacing=0.08, vertical_spacing=0.12,
        )
        for i, (name, energies) in enumerate(pdg_groups.items()):
            r, c = divmod(i, n_cols)
            fig.add_trace(
                go.Histogram(x=energies, name=name,
                             marker_color=_pdg_color(name), opacity=0.8,
                             nbinsx=max(5, len(energies)),
                             showlegend=True),
                row=r + 1, col=c + 1,
            )
        fig.update_layout(
            title="PFO Energy by PDG",
            height=max(350, n_rows * 220),
            **{k: v for k, v in lb.items() if k != "height"},
        )
        fig.update_annotations(font_color=t["text"])
        fig.update_xaxes(
            title_text="Energy [GeV]",
            gridcolor=t["grid"],
            tickfont=dict(color=t["text"]),
            title_font=dict(color=t["text"]),
        )
        fig.update_yaxes(
            title_text="N PFOs",
            gridcolor=t["grid"],
            tickfont=dict(color=t["text"]),
            title_font=dict(color=t["text"]),
        )

    # ── PDG distribution (gen) ────────────────────────────────────────────────
    elif hist_type == "pdg_gen":
        from collections import Counter
        cnt = Counter(_pdg_name(g.pdg) for g in event_data.gen_particles)
        keys, vals = zip(*sorted(cnt.items(), key=lambda x: -x[1])) if cnt else ([], [])
        colors = [_pdg_color(k) for k in keys]
        fig = go.Figure(data=[go.Bar(x=list(keys), y=list(vals), marker_color=colors)],
                        layout=go.Layout(
                            xaxis=dict(title="PDG", tickangle=-30),
                            yaxis=dict(title="Count"),
                            title="Gen Particle PDG Distribution",
                            **lb,
                        ))

    # ── PDG distribution (pfo) ────────────────────────────────────────────────
    elif hist_type == "pdg_pfo":
        from collections import Counter
        cnt = Counter(_pdg_name(p.pdg) for p in event_data.pfos)
        keys, vals = zip(*sorted(cnt.items(), key=lambda x: -x[1])) if cnt else ([], [])
        colors = [_pdg_color(k) for k in keys]
        fig = go.Figure(data=[go.Bar(x=list(keys), y=list(vals), marker_color=colors)],
                        layout=go.Layout(
                            xaxis=dict(title="PDG", tickangle=-30),
                            yaxis=dict(title="Count"),
                            title="PFO PDG Distribution",
                            **lb,
                        ))

    # ── Hits per PFO (simple) ─────────────────────────────────────────────────
    elif hist_type == "hits_per_pfo":
        from collections import Counter
        cnt = Counter(h.pfo_idx for h in event_data.hits if h.pfo_idx is not None)
        if not cnt:
            return _empty_fig("No PFO-associated hits", theme)
        pfo_ids = sorted(cnt.keys())
        fig = go.Figure(
            data=[go.Bar(
                x=[f"PFO {i}" for i in pfo_ids],
                y=[cnt[i] for i in pfo_ids],
                marker_color=t["accent"],
            )],
            layout=go.Layout(
                xaxis=dict(title="PFO", tickangle=-45),
                yaxis=dict(title="N hits"),
                title="Hits per PFO",
                **lb,
            ),
        )

    # ── Hits per PFO × detector group (stacked) ───────────────────────────────
    elif hist_type == "hits_per_pfo_by_group":
        pfo_grp: dict = defaultdict(lambda: defaultdict(int))
        for h in event_data.hits:
            if h.pfo_idx is not None:
                pfo_grp[h.pfo_idx][h.group] += 1
        if not pfo_grp:
            return _empty_fig("No PFO-associated hits", theme)
        pfo_ids = sorted(pfo_grp.keys())
        htraces = []
        for grp, info in DETECTOR_GROUPS.items():
            htraces.append(go.Bar(
                name=grp.upper(),
                x=[f"PFO {i}" for i in pfo_ids],
                y=[pfo_grp[i].get(grp, 0) for i in pfo_ids],
                marker_color=info["color"],
            ))
        fig = go.Figure(data=htraces, layout=go.Layout(
            barmode="stack",
            xaxis=dict(title="PFO", tickangle=-45),
            yaxis=dict(title="N hits"),
            title="Hits per PFO by Detector Group",
            **lb,
        ))

    # ── Hit multiplicity by group ─────────────────────────────────────────────
    elif hist_type == "hit_mult_by_group":
        from collections import Counter
        cnt = Counter(h.group for h in event_data.hits)
        labels = list(cnt.keys())
        colors = [DETECTOR_GROUPS.get(g, {}).get("color", "#aaa") for g in labels]
        fig = go.Figure(
            data=[go.Bar(x=[l.upper() for l in labels], y=list(cnt.values()),
                         marker_color=colors)],
            layout=go.Layout(
                xaxis=dict(title="Detector Group"),
                yaxis=dict(title="N hits"),
                title="Hit Multiplicity by Group",
                **lb,
            ),
        )

    # ── Gen vs PFO energy scatter ─────────────────────────────────────────────
    elif hist_type == "gen_vs_pfo_energy":
        gen_map = {g.idx: g for g in event_data.gen_particles}
        pfo_map = {p.idx: p for p in event_data.pfos}
        matched_gen_e, matched_pfo_e, hover = [], [], []
        seen: set = set()
        for h in event_data.hits:
            pair = (h.gen_idx, h.pfo_idx)
            if None in pair or pair in seen:
                continue
            seen.add(pair)
            if h.gen_idx in gen_map and h.pfo_idx in pfo_map:
                ge = gen_map[h.gen_idx]
                pe = pfo_map[h.pfo_idx]
                matched_gen_e.append(ge.energy)
                matched_pfo_e.append(pe.energy)
                hover.append(f"gen {h.gen_idx} ({_pdg_name(ge.pdg)}) ↔ pfo {h.pfo_idx} ({_pdg_name(pe.pdg)})")

        fig = go.Figure(
            data=[go.Scatter(
                x=matched_gen_e, y=matched_pfo_e, mode="markers",
                marker=dict(color=t["accent"], size=7, opacity=0.8),
                text=hover, hovertemplate="%{text}<extra></extra>",
            )],
            layout=go.Layout(
                xaxis=dict(title="Gen Energy [GeV]"),
                yaxis=dict(title="PFO Energy [GeV]"),
                title="Gen vs PFO Energy",
                **lb,
            ),
        )
        if matched_gen_e:
            mx = max(max(matched_gen_e), max(matched_pfo_e))
            fig.add_trace(go.Scatter(
                x=[0, mx], y=[0, mx], mode="lines",
                line=dict(dash="dash", color=t["muted"], width=1),
                showlegend=False,
            ))

    else:
        fig = go.Figure(layout=go.Layout(title=f"Unknown: {hist_type}", **lb))

    return fig


def _empty_fig(msg: str, theme: str) -> go.Figure:
    t = _theme(theme)
    return go.Figure(layout=go.Layout(
        paper_bgcolor=t["app_bg"],
        plot_bgcolor=t["plot_bg"],
        margin=dict(l=50, r=20, t=50, b=50),
        height=420,
        title=dict(text=msg, font=dict(color=t["muted"], size=13)),
    ))


def make_empty_3d(theme: str = "dark") -> go.Figure:
    t = _theme(theme)
    return go.Figure(layout=go.Layout(
        paper_bgcolor=t["app_bg"],
        scene=dict(bgcolor=t["plot_bg"]),
        margin=dict(l=0, r=0, t=30, b=0),
        height=580,
        title=dict(text="Load a file to begin", font=dict(color=t["muted"], size=14)),
    ))


def make_empty_hist(theme: str = "dark") -> go.Figure:
    return _empty_fig("Load a file to begin", theme)

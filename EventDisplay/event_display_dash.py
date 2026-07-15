#!/usr/bin/env python3
"""FCC Simulation Event Display — Dash + Plotly

Usage:
    source /nfs/cms/arqolmo/ExamplesFCCFullSim/setupKey4Hep.sh
    source ~/.venv/fcc-display/bin/activate
    python EventDisplay/event_display_dash.py -i sim.root --port 8050

Then open http://localhost:8050 in your browser (or SSH-tunnel the port).
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))

import dash
from dash import dcc, html, Input, Output, State, ctx, ClientsideFunction
try:
    import dash_bootstrap_components as dbc
    _HAS_DBC = True
except ImportError:
    _HAS_DBC = False

from dash import dash_table

import config
from config import (
    THEMES, DEFAULT_THEME, PDG_NAMES, pdg_label,
)
from data_io import EventReader, EventData
from matching import match_by_dr, match_by_truth_link
from geometry import get_geometry, build_geometry_traces, default_geometry_path
from plots import (
    make_3d_figure, make_histogram_figure, make_empty_3d, make_empty_hist,
    _HIST_OPTIONS,
)

_reader = EventReader()

# ── CSS variables injected into the page ─────────────────────────────────────

_CSS_VARS_DARK = """
  --cc-app-bg:     #1a1a2e;
  --cc-card-bg:    #16213e;
  --cc-sidebar-bg: #0d1b35;
  --cc-plot-bg:    #0f3460;
  --cc-text:       #e0e0e0;
  --cc-muted:      #888888;
  --cc-accent:     #e94560;
  --cc-border:     #334455;
  --cc-input-bg:   #0f1a30;
  --cc-btn-bg:     #e94560;
  --cc-btn-text:   #ffffff;
  --cc-tab-bg:     #16213e;
  --cc-tab-sel:    #0f3460;
"""

_INDEX_STRING = """<!DOCTYPE html>
<html>
  <head>
    {%metas%}
    <title>{%title%}</title>
    {%favicon%}
    {%css%}
    <style>
      :root { """ + _CSS_VARS_DARK + """ }
      * { box-sizing: border-box; }
      body { margin: 0; padding: 0; background: var(--cc-app-bg); color: var(--cc-text); }
      .drag-handle {
        width: 5px; cursor: col-resize; flex-shrink: 0;
        background: var(--cc-border);
        transition: background 0.15s;
      }
      .drag-handle:hover { background: var(--cc-accent); }
      /* Make dcc.Graph container resizable vertically */
      .resizable-graph { resize: vertical; overflow: auto; min-height: 300px; }
    </style>
  </head>
  <body>
    {%app_entry%}
    <footer>
      {%config%}
      {%scripts%}
      {%renderer%}
    </footer>
  </body>
</html>"""

# ── Table column definitions ──────────────────────────────────────────────────

MATCH_TABLE_COLUMNS = [
    {"name": "Gen idx",  "id": "gen_idx"},
    {"name": "Gen PDG",  "id": "gen_pdg"},
    {"name": "Gen name", "id": "gen_name"},
    {"name": "Gen E",    "id": "gen_E"},
    {"name": "PFO idx",  "id": "pfo_idx"},
    {"name": "PFO PDG",  "id": "pfo_pdg"},
    {"name": "PFO name", "id": "pfo_name"},
    {"name": "PFO E",    "id": "pfo_E"},
    {"name": "Value",    "id": "match_val"},
    {"name": "Type",     "id": "match_type"},
]
GEN_TABLE_COLUMNS = [
    {"name": "idx",    "id": "idx"},   {"name": "PDG",    "id": "pdg"},
    {"name": "Name",   "id": "name"},  {"name": "Status", "id": "status"},
    {"name": "E [GeV]","id": "energy"},{"name": "p [GeV]","id": "p"},
    {"name": "θ",      "id": "theta"}, {"name": "φ",      "id": "phi"},
]
PFO_TABLE_COLUMNS = [
    {"name": "idx",    "id": "idx"},   {"name": "PDG",    "id": "pdg"},
    {"name": "Name",   "id": "name"},  {"name": "Charge", "id": "charge"},
    {"name": "E [GeV]","id": "energy"},{"name": "p [GeV]","id": "p"},
    {"name": "θ",      "id": "theta"}, {"name": "φ",      "id": "phi"},
]

# ── Shared table style (always uses CSS vars for theme) ───────────────────────

_TABLE_STYLE = {
    "style_header": {
        "backgroundColor": "var(--cc-sidebar-bg)",
        "color": "var(--cc-text)",
        "fontWeight": "600",
        "border": "1px solid var(--cc-border)",
        "fontSize": "11px",
    },
    "style_data": {
        "backgroundColor": "var(--cc-card-bg)",
        "color": "var(--cc-text)",
        "border": "1px solid var(--cc-border)",
        "fontSize": "11px",
    },
    "style_data_conditional": [
        {"if": {"filter_query": "{matched} = false"},
         "opacity": "0.6", "fontStyle": "italic"},
    ],
}

# ── Layout helpers (all colors via CSS variables) ─────────────────────────────

def _sec(text: str) -> html.Div:
    return html.Div(text, style={
        "color": "var(--cc-accent)", "fontSize": "12px",
        "fontWeight": "600", "letterSpacing": "1px",
        "marginBottom": "5px", "marginTop": "6px",
    })


def _lbl(text: str) -> html.Div:
    return html.Div(text, style={
        "color": "var(--cc-muted)", "fontSize": "11px", "marginBottom": "3px",
    })


def _hr():
    return html.Hr(style={"borderColor": "var(--cc-border)", "margin": "8px 0"})


def _btn(bid: str, label: str, **extra):
    style = {
        "backgroundColor": "var(--cc-btn-bg)",
        "color": "var(--cc-btn-text)",
        "border": "none", "borderRadius": "4px",
        "padding": "5px 12px", "cursor": "pointer", "fontSize": "12px",
        **extra.pop("style", {}),
    }
    return html.Button(label, id=bid, style=style, **extra)


def _input_style(width="100%"):
    return {
        "width": width, "padding": "5px 8px",
        "borderRadius": "4px", "fontSize": "11px",
        "backgroundColor": "var(--cc-input-bg)",
        "color": "var(--cc-text)",
        "border": "1px solid var(--cc-border)",
        "boxSizing": "border-box",
    }


def _checklist_style():
    return {
        "inputStyle": {"marginRight": "4px"},
        "labelStyle": {
            "color": "var(--cc-text)", "fontSize": "11px",
            "display": "block", "marginBottom": "2px",
        },
    }


# ── Sidebar ───────────────────────────────────────────────────────────────────

def _build_sidebar() -> html.Div:
    return html.Div(
        id="sidebar",
        style={
            "backgroundColor": "var(--cc-sidebar-bg)",
            "borderRadius": "8px",
            "padding": "12px",
            "overflowY": "auto",
            "width": "280px",
            "minWidth": "160px",
            "maxWidth": "600px",
            "flexShrink": "0",
            "height": "100%",
        },
        children=[
            # ── File ──────────────────────────────────────────────────────────
            _sec("File"),
            dcc.Input(id="file-path-input", type="text",
                      placeholder="Path to .root file…",
                      style={**_input_style(), "marginBottom": "6px"}),
            html.Div(style={"display": "flex", "gap": "6px",
                            "alignItems": "center", "marginBottom": "4px"},
                     children=[
                         _btn("load-file-btn", "Load"),
                         dcc.RadioItems(
                             id="load-mode-radio",
                             options=[{"label": " Lazy", "value": "lazy"},
                                      {"label": " Eager", "value": "eager"}],
                             value="lazy", inline=True,
                             inputStyle={"marginRight": "3px"},
                             labelStyle={"color": "var(--cc-text)", "fontSize": "11px",
                                         "marginRight": "8px"},
                         ),
                     ]),
            html.Div(id="load-status-label",
                     style={"color": "var(--cc-muted)", "fontSize": "10px",
                            "marginBottom": "6px"}),
            _hr(),

            # ── Event navigation ──────────────────────────────────────────────
            _sec("Event Navigation"),
            html.Div(style={"display": "flex", "alignItems": "center",
                            "gap": "6px", "marginBottom": "6px"},
                     children=[
                         _btn("prev-event-btn", "◄",
                              style={"padding": "4px 10px"}),
                         dcc.Input(id="event-index-input", type="number",
                                   value=0, min=0, step=1,
                                   style={**_input_style("60px"),
                                          "textAlign": "center"}),
                         _btn("next-event-btn", "►",
                              style={"padding": "4px 10px"}),
                         html.Span(id="total-events-label",
                                   style={"color": "var(--cc-muted)",
                                          "fontSize": "11px"}),
                     ]),
            _hr(),

            # ── 3D filters ────────────────────────────────────────────────────
            _sec("3D Filters"),
            _lbl("Detector groups:"),
            dcc.Checklist(
                id="group-checklist",
                options=[{"label": f"  {k.upper()}", "value": k}
                         for k in config.DETECTOR_GROUPS],
                value=list(config.DETECTOR_GROUPS.keys()),
                **_checklist_style(),
            ),
            html.Div(style={"marginTop": "8px"}, children=[
                _lbl("Hit energy range [GeV]:"),
                dcc.RangeSlider(
                    id="energy-slider", min=0, max=10, step=0.001,
                    value=[0, 10], marks={0: "0", 10: "10"},
                    tooltip={"placement": "bottom", "always_visible": False},
                ),
            ]),
            html.Div(style={"marginTop": "8px"}, children=[
                _lbl("Color mode:"),
                dcc.RadioItems(
                    id="color-mode-radio",
                    options=[
                        {"label": " Detector group", "value": "detector"},
                        {"label": " Gen particle",   "value": "gen"},
                        {"label": " PFO",            "value": "pfo"},
                    ],
                    value="detector",
                    **_checklist_style(),
                ),
            ]),
            html.Div(style={"marginTop": "8px"}, children=[
                dcc.Checklist(
                    id="view-options-check",
                    options=[
                        {"label": "  Show geometry", "value": "geometry"},
                        {"label": "  Show tracks",   "value": "tracks"},
                        {"label": "  Hide secondary hits", "value": "hide_secondaries"},
                    ],
                    value=["tracks"],
                    **_checklist_style(),
                ),
            ]),
            html.Div(style={"marginTop": "8px"}, children=[
                _lbl("Filter by Gen PDG:"),
                dcc.Dropdown(id="gen-pdg-filter", options=[], value=None,
                             multi=True, placeholder="All gen particles",
                             style={"fontSize": "11px",
                                    "backgroundColor": "var(--cc-input-bg)"}),
            ]),
            html.Div(style={"marginTop": "6px"}, children=[
                _lbl("Filter by PFO PDG:"),
                dcc.Dropdown(id="pfo-pdg-filter", options=[], value=None,
                             multi=True, placeholder="All PFOs",
                             style={"fontSize": "11px",
                                    "backgroundColor": "var(--cc-input-bg)"}),
            ]),
            _hr(),

            # ── Matching ──────────────────────────────────────────────────────
            _sec("Matching"),
            _lbl("dR threshold:"),
            dcc.Slider(id="dr-max-slider", min=0.0, max=0.5, step=0.01,
                       value=0.1,
                       marks={0: "0", 0.1: "0.1", 0.3: "0.3", 0.5: "0.5"},
                       tooltip={"placement": "bottom", "always_visible": False}),
            _hr(),

            # ── Collections override ──────────────────────────────────────────
            _sec("Collections Override"),
            _lbl("Edit JSON + Apply to reload:"),
            dcc.Textarea(
                id="collections-textarea",
                value=json.dumps(
                    {g: i["collections"] for g, i in config.DETECTOR_GROUPS.items()},
                    indent=2,
                ),
                style={
                    "width": "100%", "height": "110px", "fontSize": "10px",
                    "backgroundColor": "var(--cc-input-bg)",
                    "color": "var(--cc-text)",
                    "border": "1px solid var(--cc-border)",
                    "borderRadius": "4px", "resize": "vertical",
                    "boxSizing": "border-box",
                },
            ),
            _btn("apply-collections-btn", "Apply",
                 style={"marginTop": "4px", "fontSize": "11px"}),
            _hr(),

            # ── Geometry ──────────────────────────────────────────────────────
            _sec("Geometry XML"),
            dcc.Input(
                id="geometry-path-input", type="text",
                value=default_geometry_path() or "",
                placeholder="Path to compact XML…",
                style={**_input_style(), "fontSize": "10px", "marginBottom": "4px"},
            ),
            _btn("load-geometry-btn", "Load Geometry",
                 style={"fontSize": "11px"}),
            html.Div(id="geometry-status-label",
                     style={"color": "var(--cc-muted)", "fontSize": "10px",
                            "marginTop": "4px"}),
        ],
    )


# ── Main content ──────────────────────────────────────────────────────────────

def _tab_style():
    return {
        "backgroundColor": "var(--cc-tab-bg)",
        "color": "var(--cc-text)",
        "border": "1px solid var(--cc-border)",
    }


def _tab_sel_style():
    return {
        "backgroundColor": "var(--cc-tab-sel)",
        "color": "var(--cc-accent)",
        "border": "1px solid var(--cc-border)",
        "fontWeight": "600",
    }


def _build_main_content() -> html.Div:
    return html.Div(
        id="main-content",
        style={"flex": "1 1 auto", "minWidth": "400px",
               "overflow": "auto", "height": "100%"},
        children=[
            dcc.Tabs(
                id="main-tabs", value="tab-3d",
                style={"fontSize": "12px"},
                children=[
                    # ── 3D View ────────────────────────────────────────────
                    dcc.Tab(label="3D View", value="tab-3d",
                            style=_tab_style(), selected_style=_tab_sel_style(),
                            children=[
                                html.Div(
                                    className="resizable-graph",
                                    children=[
                                        dcc.Graph(
                                            id="3d-graph",
                                            figure=make_empty_3d(DEFAULT_THEME),
                                            config={"displayModeBar": True,
                                                    "modeBarButtonsToRemove":
                                                        ["resetCameraLastSave3d"]},
                                            style={"height": "100%"},
                                        ),
                                    ],
                                    style={"height": "600px"},
                                ),
                            ]),

                    # ── Particles ──────────────────────────────────────────
                    dcc.Tab(label="Particles", value="tab-particles",
                            style=_tab_style(), selected_style=_tab_sel_style(),
                            children=[
                                dcc.Tabs(
                                    id="match-tabs", value="tab-dr",
                                    style={"fontSize": "11px", "marginTop": "8px"},
                                    children=[
                                        dcc.Tab(label="Gen Particles", value="tab-gen",
                                                style=_tab_style(),
                                                selected_style=_tab_sel_style(),
                                                children=[
                                                    dash_table.DataTable(
                                                        id="gen-table",
                                                        columns=GEN_TABLE_COLUMNS,
                                                        data=[], page_size=25,
                                                        sort_action="native",
                                                        filter_action="native",
                                                        **_TABLE_STYLE,
                                                    ),
                                                ]),
                                        dcc.Tab(label="PFOs", value="tab-pfo",
                                                style=_tab_style(),
                                                selected_style=_tab_sel_style(),
                                                children=[
                                                    dash_table.DataTable(
                                                        id="pfo-table",
                                                        columns=PFO_TABLE_COLUMNS,
                                                        data=[], page_size=25,
                                                        sort_action="native",
                                                        filter_action="native",
                                                        **_TABLE_STYLE,
                                                    ),
                                                ]),
                                        dcc.Tab(label="Gen↔PFO (dR)", value="tab-dr",
                                                style=_tab_style(),
                                                selected_style=_tab_sel_style(),
                                                children=[
                                                    dash_table.DataTable(
                                                        id="match-table-dr",
                                                        columns=MATCH_TABLE_COLUMNS,
                                                        data=[], page_size=25,
                                                        sort_action="native",
                                                        filter_action="native",
                                                        **_TABLE_STYLE,
                                                    ),
                                                ]),
                                        dcc.Tab(label="Gen↔PFO (TruthLink)",
                                                value="tab-truth",
                                                style=_tab_style(),
                                                selected_style=_tab_sel_style(),
                                                children=[
                                                    dash_table.DataTable(
                                                        id="match-table-truth",
                                                        columns=MATCH_TABLE_COLUMNS,
                                                        data=[], page_size=25,
                                                        sort_action="native",
                                                        filter_action="native",
                                                        **_TABLE_STYLE,
                                                    ),
                                                ]),
                                    ],
                                ),
                            ]),

                    # ── Histograms ──────────────────────────────────────────
                    dcc.Tab(label="Histograms", value="tab-histograms",
                            style=_tab_style(), selected_style=_tab_sel_style(),
                            children=[
                                html.Div(style={"padding": "8px"}, children=[
                                    dcc.Dropdown(
                                        id="hist-type-dropdown",
                                        options=_HIST_OPTIONS,
                                        value="hit_energy_by_detector",
                                        clearable=False,
                                        style={"fontSize": "12px",
                                               "marginBottom": "8px",
                                               "maxWidth": "420px"},
                                    ),
                                    dcc.Graph(id="histogram-graph",
                                              figure=make_empty_hist(DEFAULT_THEME)),
                                ]),
                            ]),
                ],
            ),
        ],
    )


# ── App factory ───────────────────────────────────────────────────────────────

def build_app() -> dash.Dash:
    external = [dbc.themes.DARKLY] if _HAS_DBC else []
    app = dash.Dash(__name__, external_stylesheets=external,
                    suppress_callback_exceptions=True)
    app.title = "FCC Event Display"
    app.index_string = _INDEX_STRING

    app.layout = html.Div(
        id="app-root",
        style={"backgroundColor": "var(--cc-app-bg)", "minHeight": "100vh",
               "fontFamily": "monospace", "padding": "12px",
               "color": "var(--cc-text)"},
        children=[
            # Stores
            dcc.Store(id="theme-store",          data=DEFAULT_THEME),
            dcc.Store(id="theme-dummy",           data=None),
            dcc.Store(id="events-store",          data=None),
            dcc.Store(id="current-event-store",   data=0),
            dcc.Store(id="geometry-store",        data=None),
            dcc.Store(id="collections-store",     data=None),

            # Header
            html.Div(
                style={"display": "flex", "alignItems": "center",
                       "marginBottom": "10px", "gap": "14px"},
                children=[
                    html.H4("FCC · Event Display",
                            style={"color": "var(--cc-accent)",
                                   "letterSpacing": "3px",
                                   "margin": "0", "flexGrow": "1"}),
                    html.Span(id="event-meta-label",
                              style={"color": "var(--cc-muted)",
                                     "fontSize": "11px"}),
                    html.Button(
                        "🌙", id="theme-toggle-btn",
                        style={
                            "backgroundColor": "transparent",
                            "border": "1px solid var(--cc-border)",
                            "borderRadius": "20px",
                            "color": "var(--cc-text)",
                            "cursor": "pointer",
                            "fontSize": "16px",
                            "padding": "2px 10px",
                        },
                    ),
                ],
            ),

            # Main flex layout
            html.Div(
                id="layout-container",
                style={"display": "flex", "gap": "0",
                       "height": "calc(100vh - 70px)"},
                children=[
                    _build_sidebar(),
                    html.Div(id="drag-handle-h", className="drag-handle"),
                    _build_main_content(),
                ],
            ),
        ],
    )

    _register_callbacks(app)
    return app


# ── Callbacks ─────────────────────────────────────────────────────────────────

def _register_callbacks(app: dash.Dash):

    # ── Theme: CSS variables via clientside callback ───────────────────────────
    app.clientside_callback(
        """
        function(theme) {
            var dark = {
                '--cc-app-bg':     '#1a1a2e',
                '--cc-card-bg':    '#16213e',
                '--cc-sidebar-bg': '#0d1b35',
                '--cc-plot-bg':    '#0f3460',
                '--cc-text':       '#e0e0e0',
                '--cc-muted':      '#888888',
                '--cc-accent':     '#e94560',
                '--cc-border':     '#334455',
                '--cc-input-bg':   '#0f1a30',
                '--cc-btn-bg':     '#e94560',
                '--cc-btn-text':   '#ffffff',
                '--cc-tab-bg':     '#16213e',
                '--cc-tab-sel':    '#0f3460',
            };
            var light = {
                '--cc-app-bg':     '#f0f4f8',
                '--cc-card-bg':    '#ffffff',
                '--cc-sidebar-bg': '#e8edf5',
                '--cc-plot-bg':    '#eef2ff',
                '--cc-text':       '#1a1a2e',
                '--cc-muted':      '#555555',
                '--cc-accent':     '#c0392b',
                '--cc-border':     '#d0d0d0',
                '--cc-input-bg':   '#f9f9f9',
                '--cc-btn-bg':     '#c0392b',
                '--cc-btn-text':   '#ffffff',
                '--cc-tab-bg':     '#ffffff',
                '--cc-tab-sel':    '#dce6ff',
            };
            var vars = (theme === 'light') ? light : dark;
            var root = document.documentElement;
            for (var k in vars) { root.style.setProperty(k, vars[k]); }
            return null;
        }
        """,
        Output("theme-dummy", "data"),
        Input("theme-store", "data"),
    )

    # ── Theme toggle button ────────────────────────────────────────────────────
    @app.callback(
        Output("theme-store",      "data"),
        Output("theme-toggle-btn", "children"),
        Input("theme-toggle-btn",  "n_clicks"),
        State("theme-store",       "data"),
        prevent_initial_call=True,
    )
    def toggle_theme(n_clicks, current):
        new = "light" if current == "dark" else "dark"
        return new, ("☀" if new == "light" else "🌙")

    # ── Drag handle: horizontal resize ────────────────────────────────────────
    app.clientside_callback(
        """
        function(_id) {
            var handle  = document.getElementById('drag-handle-h');
            var sidebar = document.getElementById('sidebar');
            if (!handle || !sidebar) return null;
            var dragging = false, startX = 0, startW = 0;
            handle.addEventListener('mousedown', function(e) {
                dragging = true;
                startX   = e.clientX;
                startW   = sidebar.getBoundingClientRect().width;
                document.body.style.userSelect = 'none';
                document.body.style.cursor     = 'col-resize';
                e.preventDefault();
            });
            document.addEventListener('mousemove', function(e) {
                if (!dragging) return;
                var newW = Math.max(160, Math.min(600, startW + e.clientX - startX));
                sidebar.style.width = newW + 'px';
            });
            document.addEventListener('mouseup', function() {
                dragging = false;
                document.body.style.userSelect = '';
                document.body.style.cursor     = '';
            });
            return null;
        }
        """,
        Output("drag-handle-h", "title"),   # harmless dummy output
        Input("drag-handle-h",  "id"),
    )

    # ── Geometry load ──────────────────────────────────────────────────────────
    @app.callback(
        Output("geometry-store",        "data"),
        Output("geometry-status-label", "children"),
        Input("load-geometry-btn",      "n_clicks"),
        State("geometry-path-input",    "value"),
        prevent_initial_call=False,
    )
    def load_geometry(n_clicks, xml_path):
        geo = get_geometry(xml_path or None)
        traces = build_geometry_traces(geo)
        serialized = [
            {"x": list(tr.x), "y": list(tr.y), "z": list(tr.z),
             "name": tr.name, "color": tr.line.color if tr.line else "#aaa",
             "legendgroup": tr.legendgroup,
             "showlegend": bool(tr.showlegend),
             "opacity": float(tr.opacity) if tr.opacity else 0.4}
            for tr in traces
        ]
        n_vol = len({tr.legendgroup for tr in traces})
        status = f"Loaded {n_vol} volumes from {xml_path or 'fallback defaults'}"
        return serialized, status

    # ── File load ──────────────────────────────────────────────────────────────
    @app.callback(
        Output("events-store",         "data"),
        Output("load-status-label",    "children"),
        Output("total-events-label",   "children"),
        Output("energy-slider",        "min"),
        Output("energy-slider",        "max"),
        Output("energy-slider",        "value"),
        Output("gen-pdg-filter",       "options"),
        Output("pfo-pdg-filter",       "options"),
        Output("current-event-store",  "data",   allow_duplicate=True),
        Input("load-file-btn",         "n_clicks"),
        State("file-path-input",       "value"),
        State("load-mode-radio",       "value"),
        State("collections-textarea",  "value"),
        prevent_initial_call=True,
    )
    def load_file(n_clicks, filepath, load_mode, collections_json):
        if not filepath:
            return (dash.no_update,) * 9
        try:
            coll_override = _parse_collections(collections_json)
            n = _reader.load(filepath, mode=load_mode,
                             collections_override=coll_override)
        except Exception as exc:
            return (None, f"Error: {exc}", "—", 0, 10, [0, 10], [], [], 0)

        ev0 = _reader.get_event(0)
        energies = [h.energy for h in ev0.hits if h.energy > 1e-9]
        e_max = round(max(energies) * 1.1, 6) if energies else 10.0

        gen_pdgs = sorted({g.pdg for g in ev0.gen_particles})
        pfo_pdgs = sorted({p.pdg for p in ev0.pfos})
        gen_opts = [{"label": pdg_label(p), "value": p} for p in gen_pdgs]
        pfo_opts = [{"label": pdg_label(p), "value": p} for p in pfo_pdgs]

        store = {"loaded": True, "filepath": filepath, "n_events": n}
        return (store,
                f"Loaded: {os.path.basename(filepath)} ({n} events, mode={load_mode})",
                f"/ {n}", 0.0, e_max, [0.0, e_max],
                gen_opts, pfo_opts, 0)

    # ── Event navigation ───────────────────────────────────────────────────────
    @app.callback(
        Output("current-event-store", "data"),
        Output("event-index-input",   "value"),
        Input("prev-event-btn",       "n_clicks"),
        Input("next-event-btn",       "n_clicks"),
        Input("event-index-input",    "value"),
        State("current-event-store",  "data"),
        State("events-store",         "data"),
        prevent_initial_call=True,
    )
    def navigate_event(prev, nxt, manual_val, current, events_store):
        if not events_store or not events_store.get("loaded"):
            return 0, 0
        n = events_store.get("n_events", 1)
        tid = ctx.triggered_id
        if tid == "prev-event-btn":
            idx = max(0, (current or 0) - 1)
        elif tid == "next-event-btn":
            idx = min(n - 1, (current or 0) + 1)
        else:
            idx = max(0, min(n - 1, int(manual_val or 0)))
        return idx, idx

    # ── 3D view ────────────────────────────────────────────────────────────────
    @app.callback(
        Output("3d-graph",           "figure"),
        Output("event-meta-label",   "children"),
        Input("current-event-store", "data"),
        Input("group-checklist",     "value"),
        Input("energy-slider",       "value"),
        Input("color-mode-radio",    "value"),
        Input("view-options-check",  "value"),
        Input("gen-pdg-filter",      "value"),
        Input("pfo-pdg-filter",      "value"),
        Input("theme-store",         "data"),
        State("events-store",        "data"),
        State("geometry-store",      "data"),
        prevent_initial_call=True,
    )
    def update_3d(event_idx, groups, energy_range, color_mode,
                  view_opts, gen_pdg_filter, pfo_pdg_filter,
                  theme, events_store, geo_store):
        if not events_store or not events_store.get("loaded"):
            return make_empty_3d(theme or DEFAULT_THEME), "No file loaded"
        try:
            ev = _reader.get_event(int(event_idx or 0))
        except Exception as exc:
            return make_empty_3d(theme or DEFAULT_THEME), f"Error: {exc}"

        view_opts   = view_opts or []
        show_geo    = "geometry" in view_opts
        show_tracks = "tracks" in view_opts
        hide_sec    = "hide_secondaries" in view_opts

        # Rebuild Plotly geometry traces from stored serialized data
        geo_traces = []
        if show_geo and geo_store:
            import plotly.graph_objects as _go
            for s in geo_store:
                geo_traces.append(_go.Scatter3d(
                    x=s["x"], y=s["y"], z=s["z"],
                    mode="lines",
                    line=dict(color=s["color"], width=1),
                    name=s["name"],
                    legendgroup=s.get("legendgroup", s["name"]),
                    showlegend=s.get("showlegend", False),
                    hoverinfo="skip",
                    opacity=s.get("opacity", 0.4),
                ))

        # pfo → gen map for track coloring/filtering.
        # RecoMCTruthLink is many-to-many: one PFO can be linked to several MCParticles
        # with different weights (track contribution vs. cluster contribution).
        # We decode the Bohdan Dudar encoding and prefer the MCParticle with the highest
        # track_w > 0 (must be charged). If all track_w == 0 (neutral PFO), use cluster_w.
        pfo_to_gen: dict = {}
        if ev.truth_links:
            _candidates: dict = {}   # pfo_idx → (best_track_w, best_cluster_w, gen_idx)
            for lnk in ev.truth_links:
                if lnk.pfo_idx < 0 or lnk.gen_idx < 0:
                    continue
                w = int(lnk.weight)
                track_w   = (w % 10000) / 1000.0
                cluster_w = (w // 10000) / 1000.0
                prev = _candidates.get(lnk.pfo_idx)
                if prev is None or (track_w, cluster_w) > (prev[0], prev[1]):
                    _candidates[lnk.pfo_idx] = (track_w, cluster_w, lnk.gen_idx)
            for pfo_idx, (_, _, gen_idx) in _candidates.items():
                pfo_to_gen[pfo_idx] = gen_idx
        else:
            for m in match_by_dr(ev.gen_particles, ev.pfos, max_dr=0.4):
                if m["matched"] and m["pfo_idx"] >= 0:
                    pfo_to_gen[m["pfo_idx"]] = m["gen_idx"]

        gen_filter = None
        if gen_pdg_filter:
            pdg_set = set(gen_pdg_filter)
            gen_filter = [g.idx for g in ev.gen_particles if g.pdg in pdg_set]

        pfo_filter = None
        if pfo_pdg_filter:
            pdg_set = set(pfo_pdg_filter)
            pfo_filter = [p.idx for p in ev.pfos if p.pdg in pdg_set]

        fig = make_3d_figure(
            ev.hits, ev.tracks, ev.gen_particles, ev.pfos, geo_traces,
            active_groups   =set(groups or []),
            color_mode      =color_mode or "detector",
            energy_range    =tuple(energy_range) if energy_range else (0, 1e9),
            gen_filter      =gen_filter,
            pfo_filter      =pfo_filter,
            show_geometry   =show_geo,
            show_tracks     =show_tracks,
            hide_secondary_hits=hide_sec,
            pfo_to_gen      =pfo_to_gen,
            theme           =theme or DEFAULT_THEME,
        )
        meta = (
            f"Event {ev.global_idx} | {os.path.basename(ev.file_path)} "
            f"| local #{ev.local_idx} | "
            f"{len(ev.gen_particles)} gen · {len(ev.pfos)} PFOs · {len(ev.hits)} hits"
        )
        return fig, meta

    # ── Particle tables ────────────────────────────────────────────────────────
    @app.callback(
        Output("gen-table",         "data"),
        Output("pfo-table",         "data"),
        Output("match-table-dr",    "data"),
        Output("match-table-truth", "data"),
        Input("current-event-store","data"),
        Input("dr-max-slider",      "value"),
        State("events-store",       "data"),
        prevent_initial_call=True,
    )
    def update_particle_tables(event_idx, dr_max, events_store):
        if not events_store or not events_store.get("loaded"):
            return [], [], [], []
        try:
            ev = _reader.get_event(int(event_idx or 0))
        except Exception:
            return [], [], [], []

        gen_rows = [
            {"idx": g.idx, "pdg": g.pdg,
             "name": PDG_NAMES.get(int(g.pdg), str(g.pdg)),
             "status": g.status, "energy": round(g.energy, 4),
             "p": round(g.p, 4), "theta": round(g.theta, 4),
             "phi": round(g.phi, 4)}
            for g in ev.gen_particles
        ]
        pfo_rows = [
            {"idx": p.idx, "pdg": p.pdg,
             "name": PDG_NAMES.get(int(p.pdg), str(p.pdg)),
             "charge": p.charge, "energy": round(p.energy, 4),
             "p": round(p.p, 4), "theta": round(p.theta, 4),
             "phi": round(p.phi, 4)}
            for p in ev.pfos
        ]
        dr_rows    = match_by_dr(ev.gen_particles, ev.pfos,
                                 max_dr=float(dr_max or 0.1))
        truth_rows = match_by_truth_link(ev.truth_links,
                                         ev.gen_particles, ev.pfos)
        return gen_rows, pfo_rows, _clean(dr_rows), _clean(truth_rows)

    # ── Histograms ─────────────────────────────────────────────────────────────
    @app.callback(
        Output("histogram-graph",    "figure"),
        Input("current-event-store", "data"),
        Input("hist-type-dropdown",  "value"),
        Input("theme-store",         "data"),
        State("events-store",        "data"),
        prevent_initial_call=True,
    )
    def update_histogram(event_idx, hist_type, theme, events_store):
        if not events_store or not events_store.get("loaded"):
            return make_empty_hist(theme or DEFAULT_THEME)
        try:
            ev = _reader.get_event(int(event_idx or 0))
        except Exception:
            return make_empty_hist(theme or DEFAULT_THEME)
        return make_histogram_figure(
            ev, hist_type or "hit_energy_by_detector",
            theme=theme or DEFAULT_THEME,
        )

    # ── Collections apply ──────────────────────────────────────────────────────
    @app.callback(
        Output("events-store",      "data",    allow_duplicate=True),
        Output("load-status-label", "children", allow_duplicate=True),
        Input("apply-collections-btn", "n_clicks"),
        State("collections-textarea",  "value"),
        State("load-mode-radio",       "value"),
        State("events-store",          "data"),
        prevent_initial_call=True,
    )
    def apply_collections(n_clicks, collections_json, load_mode, events_store):
        if not events_store or not events_store.get("loaded"):
            return dash.no_update, "Load a file first."
        filepath = events_store.get("filepath", "")
        try:
            coll_override = _parse_collections(collections_json)
            _reader.load(filepath, mode=load_mode,
                         collections_override=coll_override)
            return dict(events_store), f"Collections reloaded."
        except Exception as exc:
            return dash.no_update, f"Collections error: {exc}"


# ── Helpers ───────────────────────────────────────────────────────────────────

def _parse_collections(json_str):
    if not json_str:
        return None
    raw = json.loads(json_str)
    return {
        grp: {
            "collections": cols,
            "color": config.DETECTOR_GROUPS.get(grp, {}).get("color", "#aaaaaa"),
        }
        for grp, cols in raw.items()
    }


def _clean(rows: list) -> list:
    out = []
    for r in rows:
        row = {}
        for k, v in r.items():
            row[k] = "—" if isinstance(v, float) and math.isnan(v) else v
        out.append(row)
    return out


# ── Entry point ───────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="FCC Event Display — Dash + Plotly")
    p.add_argument("-i", "--input", default="", help="Initial ROOT file")
    p.add_argument("--detector", default="CLD",
                   choices=sorted(config.DETECTOR_PROFILES),
                   help="Detector concept: collections + geometry (default: CLD)")
    p.add_argument("--geometry-version", default=None,
                   help="Geometry version tag (e.g. v01/v02 for ILD, v06 for CLD); "
                        "default is the detector's default version")
    p.add_argument("--port",  type=int, default=8050)
    p.add_argument("--host",  default="0.0.0.0")
    p.add_argument("--debug", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    config.select_detector(args.detector, args.geometry_version)
    print(f"Detector: {config.ACTIVE_DETECTOR}  geometry: {config.GEOMETRY_VERSION} "
          f"({config.GEOMETRY_RELPATH})")
    app = build_app()

    if args.input and os.path.isfile(args.input):
        print(f"Pre-loading {args.input}…")
        try:
            n = _reader.load(args.input, mode="lazy")
            print(f"  {n} events available.")
        except Exception as exc:
            print(f"  Warning: {exc}")

    print(f"\nServer at  http://{args.host}:{args.port}")
    print("Ctrl+C to stop.\n")
    app.run(debug=args.debug, port=args.port, host=args.host)


if __name__ == "__main__":
    main()

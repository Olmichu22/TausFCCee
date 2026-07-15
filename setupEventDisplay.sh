# Environment setup for the Dash Event Display (EventDisplay/).
#
# Usage (must be *sourced*, not executed):
#   source setupEventDisplay.sh
#
# What it does:
#   1. Loads the Key4hep stack (release with podio >= 1.7, required to read
#      the reco EDM4hep files) unless one is already loaded.
#   2. Creates the dash/plotly virtualenv on first use (installs it if missing).
#   3. Activates the virtualenv on top of the Key4hep python.
#
# Overridable via environment variables:
#   KEY4HEP_RELEASE   Key4hep release tag        (default: 2026-04-08)
#   FCC_DISPLAY_VENV  virtualenv location        (default: ~/.venv/fcc-display-latest)

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
    echo "ERROR: this script must be sourced:  source $0" >&2
    exit 1
fi

KEY4HEP_RELEASE="${KEY4HEP_RELEASE:-2026-04-08}"
FCC_DISPLAY_VENV="${FCC_DISPLAY_VENV:-$HOME/.venv/fcc-display-latest}"

# 1) Key4hep stack (skip only if the *same* release is already sourced;
#    a different release cannot be unloaded — needs a fresh shell)
if [ -z "$KEY4HEP_STACK" ]; then
    echo ">> Loading Key4hep release $KEY4HEP_RELEASE ..."
    source /cvmfs/sw.hsf.org/key4hep/setup.sh -r "$KEY4HEP_RELEASE" || return 1
elif [[ "$KEY4HEP_STACK" == *"/releases/$KEY4HEP_RELEASE/"* ]]; then
    echo ">> Key4hep release $KEY4HEP_RELEASE already loaded."
else
    echo "ERROR: a different Key4hep release is already loaded in this shell:" >&2
    echo "       $KEY4HEP_STACK" >&2
    echo "       Start a fresh shell (without sourcing any other Key4hep setup) and retry." >&2
    return 1
fi

# 2) Create the venv if it does not exist yet.
#    --system-site-packages so podio/ROOT from Key4hep remain importable.
#    NOTE: the venv is tied to the python version of the Key4hep release used
#    to create it; if a new release changes python, delete the venv and re-source.
if [ ! -f "$FCC_DISPLAY_VENV/bin/activate" ]; then
    echo ">> Virtualenv not found — creating it at $FCC_DISPLAY_VENV ..."
    python3 -m venv --system-site-packages "$FCC_DISPLAY_VENV" || return 1
    "$FCC_DISPLAY_VENV/bin/pip" install --upgrade pip
    "$FCC_DISPLAY_VENV/bin/pip" install "dash>=4" plotly dash-bootstrap-components || return 1
fi

# 3) Activate it
source "$FCC_DISPLAY_VENV/bin/activate"
echo ">> Event Display environment ready (venv: $FCC_DISPLAY_VENV)"
echo ">> Run:  python EventDisplay/event_display_dash.py -i <file.root> [--detector CLD|ILD]"

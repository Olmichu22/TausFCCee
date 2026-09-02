# FCC analysis dependencies

The frozen FCC study was validated with Key4hep `2026-08-21`; the exact IFIC
setup reference is recorded in `configs/environments/key4hep_2026-08-21.yaml`.
ROOT, podio, edm4hep and the Key4hep C++ stack are external environment
dependencies and are not installed by pip.

Because the Key4hep CVMFS Python prefix is read-only, use the
`--system-site-packages` virtual environment documented in `QUICKSTART.md` and
install this repository with `pip --no-build-isolation -e .` inside it.

Python-only helpers use NumPy, pandas, matplotlib, PyArrow and PyYAML. Minimal
package metadata is provided by `pyproject.toml` without converting or moving
the established upstream module layout.

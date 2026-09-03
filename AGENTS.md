# Repository guardrails

`TausFCCee` owns the maintained FCC analysis layer: G geometric association,
HitAnalysis, frozen truth/PID/tau/photon-origin definitions, G/L comparisons,
migration studies, and scientific documentation. It also retains broader
upstream/legacy tau-analysis code.

For maintained FCC work, start with `readme.md`, `docs/QUICKSTART.md`,
`docs/ANALYSIS_RECIPES.md`, `docs/SCIENTIFIC_DEFINITIONS.md`, and
`docs/TROUBLESHOOTING.md`. Later legacy sections of `readme.md` may use other
environment, path, and output conventions. Treat current code, configs, tests,
and maintained documentation as authoritative.

Preserve these boundaries:

- the 2026-09-01 scientific freeze;
- G belongs here; L_direct and L_ancestor belong to `FCC-tau-workflow` and
  must not be reimplemented here;
- cross-repository exchange is data-only, with no runtime imports or absolute
  clone dependencies;
- `(sample, source_file_id, event_in_file)` is the cross-repository event key;
- raw HitAnalysis file ordinals must be normalized through the input manifest
  before joining workflow products;
- P8C and P8O remain distinct reconstruction provenances.

Use the narrowest relevant validators, distinguish maintained FCC behavior
from legacy behavior, and never claim an unrun test or production step passed.
Do not commit large generated data or silently overwrite outputs.

Final pushes are manual unless the project owner explicitly authorizes one.
Avoid history rewriting and do not alter remotes without explicit approval.

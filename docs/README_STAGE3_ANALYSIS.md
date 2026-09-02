# Stage 3 maintained FCC analysis layer

Stage 3 preserves the upstream TausFCCee layout and adds the minimum portable
FCC analysis structure. The authoritative geometric association remains in
`modules/NeutralRecover.py`, called by the hardened parallel HitAnalysis. Its
selected-truth, PID, tau-origin and photon-origin definitions are centralized
in `modules/fcc_truth_definitions.py`.

Workflow L_direct/L_ancestor products enter only as data through a versioned
manifest and Parquet schema checks in `modules/fcc_workflow_interface.py`.
There is no import or runtime path dependency on FCC-tau-workflow, and no
second L_direct/L_ancestor assignment implementation.

Maintained executable analysis lives in `scripts/analysis/`: a portable
HitAnalysis wrapper and parameterized migration-matrix postprocessing.
Freeze-specific equal-statistics diagnostics were classified as scientific
support rather than copied into the active runtime path. Their final reports
are under `docs/scientific/supporting_reports/`.

Stage 4 must provide final collaborator-facing operational documentation. In
particular, `RECONSTRUCTION_CHAIN.md` must document the validated REC-only
truth-link production mode and a future integrated reconstruction/linker mode
labelled `SUPPORTED DESIGN / NOT YET REGRESSION-VALIDATED`, including the
required A/B validation and future skip-if-valid-relation behavior. Stage 3
does not implement that design.

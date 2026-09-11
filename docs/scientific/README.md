# FCC tau scientific documentation

**SCIENTIFIC ANALYSIS FROZEN AT 2026-09-01**

This directory contains scientific reports and freeze metadata. Operational
instructions are intentionally separate: use the repository
[quick start](../QUICKSTART.md), [analysis recipes](../ANALYSIS_RECIPES.md),
[scientific definitions](../SCIENTIFIC_DEFINITIONS.md), and
[troubleshooting guide](../TROUBLESHOOTING.md) to run maintained software.

## Primary scientific freeze

- [FCC_TAU_W_P8O_RECONSTRUCTION_CONSOLIDATED_REPORT_20260901.md](consolidated/FCC_TAU_W_P8O_RECONSTRUCTION_CONSOLIDATED_REPORT_20260901.md)
- [PDF](consolidated/FCC_TAU_W_P8O_RECONSTRUCTION_CONSOLIDATED_REPORT_20260901.pdf)
- [TeX](consolidated/FCC_TAU_W_P8O_RECONSTRUCTION_CONSOLIDATED_REPORT_20260901.tex)

## Supporting current reports

- [Migration matrices](supporting_reports/MIGRATION_MATRICES_W_P8_G_LDIRECT_LANCESTOR_20260828.md)
- [Tau topology](supporting_reports/TAU_TOPOLOGY_W_P8O_EQUALSTATS_AUDIT_20260831.md)
- [Reconstructed PID](supporting_reports/PID_AUDIT_W_P8O_EQUALSTATS_20260901.md)
- [PFO reconstruction efficiency](supporting_reports/PFO_RECONSTRUCTION_EFFICIENCY_W_P8O_EQUALSTATS_20260901.md)
- [Photon-origin discovery](supporting_reports/PHOTON_ORIGIN_DISCOVERY_W_P8O_EQUALSTATS_20260901.md)
- [Photon-origin reconstruction](supporting_reports/PHOTON_ORIGIN_RECO_AUDIT_W_P8O_EQUALSTATS_20260901.md)

## Interpretation limits

W `parentless_unresolved` photons are not proven ISR, and
`tau_direct_daughter_unresolved` photons are not proven FSR. Only the frozen
P8O `explicit_ISR` category carries its defined positive stored-genealogy
evidence. No operational page changes these interpretations.

Small freeze manifests and validation metadata are under
[`freezes/20260901/`](freezes/20260901/). They preserve hashes and references
to machine-readable artifacts kept outside Git. Large Parquet/CSV collections,
plot trees, and production files are not duplicated here.


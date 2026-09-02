# FCC analysis consolidation provenance — Stage 3

## Baseline and sources

- Target base commit: `a2fd560254dcd2e747478273608da4f31959844b`.
- Primary analysis snapshot:
  `/lhome/ific/a/airqui/FCC/archive/consolidation_20260901/snapshots/TausFCCee_pre_consolidation_20260901/`.
- FCC workflow-analysis snapshot:
  `/lhome/ific/a/airqui/FCC/archive/consolidation_20260901/snapshots/ild-tau_pre_consolidation_20260901/`.
- Workflow contract consulted read-only: `fcc_tau_association_v1`.
- Scientific freeze date: `2026-09-01`.
- Complete candidate classification and source hashes:
  `provenance/source_port_manifest_stage3.csv`.

## Audited tracked modifications

- `.gitignore` — **REQUIRED_FCC_MAINTAINED**: runtime/product exclusions;
  scientific manifests, JSON and selected PDFs remain versionable.
- `HitAnalysis/particle_level_analisis_parallel.py` —
  **REQUIRED_FCC_MAINTAINED**: tuple/composite event identity, loud worker
  failures and hardened association behavior; Stage 3 adds explicit output-root
  and G-threshold options without changing defaults.
- `modules/NeutralRecover.py` — **REQUIRED_FCC_MAINTAINED**: detector-signal
  gate correction and authoritative G; selected-truth predicates now use the
  centralized byte-equivalent definition.
- `modules/myutils.py` — **REQUIRED_FCC_MAINTAINED**: prediction keys no longer
  assume 1000 events/file; explicit output root added with legacy default intact.
- `modules/ConfusionMatrixParticleLevel.py` — **GENERAL_UPSTREAM_WORTHY**:
  figures are closed after saving; scientific content is unchanged.

No tracked modification was classified as a required local diagnostic or
obsolete behavior.

## Critical source hashes

- Hardened HitAnalysis: `850afb0c8e0949fd60a2c70a868d4cdf0bbbb9452c9816b390520d2aa6c07262`.
- Authoritative G source: `ec43f463efe497ee48c9404858590de95250deff5970db1befc2ebd2cd20abf0`.
- Event hardening source: `c098ac97129292cdbd320feadd000fc150674e6bee10f55003cec2e58c4ea45`.
- Migration analysis source: `10b294045d26dd77f3a6555a2eab9bd19100cba5dbc2687520378453df261d20`.
- Consolidated report MD: `46df78cc41b467a36502aeac955cd741304fef35819cf74527d8442fb3e3665f`.
- Consolidated report PDF: `d1cdfc773a99aa2840c2915d8d6afd056d16ea7a152fd723388b5f8926bc8449`.

## Boundary

G remains implemented once, in TausFCCee. L_direct and L_ancestor are not
implemented here: exported workflow data are consumed through manifests and
the `fcc_tau_association_v1` physical schema. Freeze-specific diagnostic
builders and large artifacts were not copied into the active code path.

No scientific definitions were intentionally changed during Stage 3. No
simulation, reconstruction, linker, HitAnalysis production or scientific
sample analysis was run during this port.

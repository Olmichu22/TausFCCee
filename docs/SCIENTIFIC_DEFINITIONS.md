# Frozen FCC scientific definitions

This page is the human-readable reference for the maintained FCC analysis
definitions. Machine-readable authorities are
`modules/fcc_truth_definitions.py`, `modules/NeutralRecover.py`,
`configs/analysis/geometric_association_v1.yaml`, and the workflow association
contract. Changing these definitions requires a separately reviewed scientific
change; operational recipes must not reinterpret them.

## Selected truth (`selected_truth_v1`)

An MCParticle is selected truth when all conditions hold:

- `generatorStatus == 1`;
- `abs(PDG)` is not 12, 14, or 16—electron, muon, and tau neutrinos are
  explicitly excluded;
- momentum magnitude is at least `1e-10 GeV`.

## G: geometric selected-truth association

G compares selected truth with non-zero-momentum `PandoraPFOs` using

```text
d = sqrt((theta_gen-theta_reco)^2 + wrap(phi_gen-phi_reco)^2)
```

and the strict condition `d < 0.1`. Default deduplication is `reco`: if several
truth particles claim the same PFO, the smallest-distance truth keeps it and
the losing truths become unmatched. Selected truths without a retained PFO are
unmatched; unused PFOs are fakes. The detector-signal gate in the authoritative
implementation is preserved.

`AssocMaxDR=0.1` belongs to G. It is not the tau/neutral-recovery cone
`dRMax=0.4`.

## Event identity

The authoritative cross-file event key is
`(sample, source_file_id, event_in_file)`. Within G products the essential
pair is `(source_file_id, event_in_file)`; add sample when joining datasets.
PFO and MCParticle identities add `pfo_index` and `mc_index` respectively.

Raw HitAnalysis Parquet currently stores a deterministic zero-based input-file
ordinal in its `source_file_id` column. Before a cross-repository join, map
that ordinal to the stable `source_file_id` declared by the corresponding
input-manifest entry. The declared stable ID—not the local ordinal—is the
authoritative cross-repository field. The maintained smoke validator performs
and checks this mapping explicitly.

Cantor or packed integer event IDs are derived conveniences. Never assume
`file_index*1000+event` or `file_index*2000+event`, and never use a convenience
integer in place of the composite key across repositories.

## Reconstructed PID mapping (`fcc_reco_pid_v1`)

| Absolute PDG | Category |
|---:|---|
| 11 | electron |
| 13 | muon |
| 22 | photon |
| 211 | charged pion |
| 310 | K0S |
| 2112 | neutron |
| 3122 | Lambda |

`999` is an association sentinel for unmatched outcomes, not a reconstructed
PID. Unsupported reconstructed PDGs are errors in maintained matrix inputs;
they must not be silently folded into a listed category.

## `tau_origin`

A selected truth particle has tau origin when recursive traversal of its stored
MCParticle parents contains any particle with `abs(PDG)==15`. Traversal uses
the stored graph, rejects invalid indices, and treats a cycle as fatal. It does
not infer missing parents.

## `photon_origin_v1`

The frozen categories derive from stored ancestry and the already validated
beam-chain predicate:

- `tau_decay_pi0`;
- `tau_direct_daughter_unresolved`;
- `tau_other_descendant`;
- `explicit_ISR` when the frozen beam-chain evidence is positive;
- `parentless_unresolved`;
- `electron_parent_non_tau_unresolved`;
- `positron_parent_non_tau_unresolved`;
- `photon_chain_non_tau_unresolved`;
- `resolved_non_tau_other`.

Interpretation limits are mandatory:

- P8O `explicit_ISR` has the defined positive stored-genealogy evidence.
- W `parentless_unresolved` is not proven ISR and must not be called ISR.
- `tau_direct_daughter_unresolved` is not proven FSR and must not be called
  FSR.

No additional ISR/FSR interpretation is part of `photon_origin_v1`.

## Fiducial reconstructed outcome

Maintained migration matrices first choose the representative PFO, then apply
the reconstructed-level fiducial selection:

- PFO energy `> 1.0 GeV`;
- `1 degree < theta_reco < 179 degrees`.

An assigned representative that fails becomes
`assigned_but_fails_reco_selection`. Association-unmatched and ambiguous
outcomes remain distinct. This is not a truth-level selection.

## G, L_direct, and L_ancestor

- G answers which selected truth and PFO are geometrically compatible.
- L_direct answers which immediate detector-level MC contributor wins the
  frozen exact packed-weight rule. The relation weight is decoded as
  `W = 10000*C + T`, `T = int(W) % 10000`, and `C = int(W) // 10000`. A
  relation with packed zero remains a relation; no relation is a distinct
  state. The immediate linked MC may have any `generatorStatus`.
- L_ancestor starts only from a unique L_direct MC. If it satisfies
  `selected_truth_v1`, it is retained at depth zero; otherwise a parent-only
  breadth-first traversal promotes to exactly one qualifying selected truth at
  the nearest depth. Several qualifying nearest ancestors are ambiguous, none
  gives `no_selected_ancestor`, and a stored cycle is fatal. It never reruns or
  reweights T/C and cannot repair missing or ambiguous L_direct.

When an MC-to-PFO inversion needs one reconstructed object, the authoritative
`truthlink_representative` helper ranks the underlying direct TruthLink weights:
maximum T then C on the track branch, or maximum C on the cluster branch. Exact
terminal ties remain `ambiguous_multiple_pfo`; index, order, PID, energy,
residual, and PDG are not tie-breakers. This selects one PFO for one-truth/
one-PFO observables and does not redefine the persisted PFO-to-MC assignments.

They answer different questions, are complementary diagnostics, and must not
be substituted for one another. L_direct/L_ancestor are produced by
FCC-tau-workflow and consumed here as versioned data; they are not implemented
again in TausFCCee.

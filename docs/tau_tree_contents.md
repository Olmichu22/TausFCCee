# `Tau_tree` — branch reference

Contents of the ROOT trees produced by [TauAnalysis/TTreesTausLong.py](TauAnalysis/TTreesTausLong.py).

* Tree name: `Tau_tree`, 99 branches, **one entry per event**.
* Every physics object lives in a `std::vector` branch, so each entry holds a variable
  number of taus, constituents and photons.
* The `config.yaml` written next to the ROOT file is a snapshot of the exact cuts,
  sample and flags of that production. Read it before comparing two files.

## How the file is produced

`TTreesTausLong.py` splits the input EDM4hep files over `--n-workers` forked workers,
each writing `tmp_chunk_<id>.root`, and merges them with `hadd`. If any worker fails the
job aborts without merging, so a file that exists is always complete. Per event it:

1. reads `MCParticles` and `PandoraPFOs`
2. builds gen taus (`tauReco.findAllGenTaus`, status-2 taus) and reco tau / electron /
   muon candidates (`tauReco.findAllTaus`, `electronReco`, `muonReco`);
3. walks each gen tau's decay tree down to final-state, non-neutrino leaves;
4. matches gen to reco taus by the angle between the gen *visible* momentum and the reco
   momentum, and gen to reco photons through `RecoMCTruthLink`.

## Quick start

```python
import ROOT
r = ROOT.RDataFrame("Tau_tree", "Tree_effisdecayAll_0.4_tph0.0_tpi0.0_n0.0_g0.0.root")
r.Filter("Sum(GenTauType == 1) > 0").Histo1D("GenVisTauP")
```

With `uproot` (no key4hep environment needed):

```python
import uproot, awkward as ak
t = uproot.open("Tree_....root")["Tau_tree"]
a = t.arrays(["GenTauType", "GenTauTrueMode", "GenVisTauP"])
ak.sum(a.GenTauTrueMode == 11)          # true tau -> pi pi0
```

## Index conventions

Branches come in blocks. Inside a block all vectors have the same length and are aligned
element by element; blocks are linked through explicit key branches.

| Block | Length per entry | Prefix |
| --- | --- | --- |
| Gen taus | `numGenTaus` | `GenTau*`, `GenVisTau*`, `GenEventId`, `GenMatchedKey` |
| Gen tau constituents | Σ `GenTauNConsts` | `GenConst*`, `GenTauConstKey` |
| Gen decay daughters | Σ `GenTauNDecayDaughters` | `GenDecayDaughter*` |
| Gen extra neutrals | `numGenExtraNeutrals` = Σ `GenTauNExtraNeutrals` | `GenExtraNeutral*` |
| Gen neutrinos | `numGenNus` = Σ `GenTauNNus` | `GenNu*` |
| Reco taus | `numRecoTaus` | `RecoTau*` |
| Reco tau constituents | Σ `RecoTauNConsts` | `RecoConst*`, `RecoTauConstKey` |
| Gen photons (whole event) | `numGenPhotons` | `GenPhoton*` |
| Reco photons (`PandoraPFOs`) | `numRecoPhotons` | `RecoPhoton*` |

Two kinds of index appear throughout:

* **Key** (`*Key`) — position inside a block of *this* entry. `GenTauConstKey[j] == 2`
  means constituent `j` belongs to the third gen tau of the event. `-1` means "no link".
* **MC index** (`*MCIdx`) — `getObjectID().index` in the `MCParticles` collection of the
  original EDM4hep event. Stable within an event, so it is the safe way to join across
  blocks; `-1` means unavailable.

## Event-level scalars

| Branch | Type | Meaning |
| --- | --- | --- |
| `numGenTaus` | `int` | Gen taus in the entry (`generatorStatus == 2`) |
| `numRecoTaus` | `int` | Reco candidates = hadronic taus + electrons + muons, in that order |
| `numGenPhotons` | `int` | Generator photons with `generatorStatus == 1`, whole event |
| `numRecoPhotons` | `int` | Photon PFOs in `PandoraPFOs` |
| `numGenExtraNeutrals` | `int` | Untagged neutrals summed over all gen taus of the entry |
| `numGenNus` | `int` | Neutrinos summed over all gen taus of the entry |
| `beamE` | `double` | Energy of `MCParticles[0]` — the beam energy (45.6 GeV here) |

## Gen taus

| Branch | Type | Meaning |
| --- | --- | --- |
| `GenEventId` | `int` | Global event id, repeated once per gen tau; unique across the file |
| `GenTauP`, `GenTauPt`, `GenTauEta`, `GenTauTheta`, `GenTauPhi`, `GenTauMass` | `float` | Full tau 4-momentum, neutrinos included |
| `GenVisTauP`, `GenVisTauPt`, `GenVisTauEta`, `GenVisTauTheta`, `GenVisTauPhi`, `GenVisTauMass` | `float` | Visible 4-momentum, neutrinos removed |
| `GenTauType` | `int` | Visible-topology decay code, see below |
| `GenTauTrueMode` | `int` | True decay mode from the direct daughters, see below |
| `GenTauQ` | `float` | Visible charge |
| `GenTauDR` | `float` | Largest opening angle between constituents |
| `GenTauNConsts` | `int` | Number of final-state constituents of this tau |
| `GenTauNConstKey` | `int` | Own index of the tau (`0, 1, …`); redundant with the position |
| `GenMatchedKey` | `int` | Gen tau index of the gen↔reco pair |
| `RecoMatchedKey` | `int` | Reco tau matched to `GenMatchedKey[i]`, `-1` if none within `MatchedGenMinDR` |

`GenEventId` encodes provenance as `(file_index * 1000) + event_in_file`.

Matching is nearest-neighbour in angle between `GenVisTau` and `RecoTau` momenta, with no
uniqueness requirement: two gen taus can in principle point at the same reco tau.

### `GenTauType` — visible topology

Charged prongs and pi0s as a detector could count them (`tauReco.findAllGenTaus`).


| Code | Decay |
| --- | --- | 
| `0` | 1 prong, 0 pi0 | 
| `N` (1–9) | 1 prong, N pi0 |
| `10 + N` | 3 prongs, N pi0 |
| `-11` | tau → e ν ν | 
| `-13` | tau → μ ν ν | 
| `-2` | other charged content (unclassified) |
| `-1` | no valid assignment (e.g. \|charge\| ≠ 1) |


### `GenTauTrueMode` — true decay mode

Built from the PDG codes of the
tau's **direct daughters**: neutrinos and tau-FSR photons dropped, tau+ charge-conjugated
onto the tau- convention, sorted canonically, then looked up in `MODE_TABLE`.

This is complementary to `GenTauType`, not a refinement of it. `GenTauType` is the right
quantity to compare against reco, but it merges physically different channels whenever an
intermediate resonance is invisible: `tau → K0 π` lands on `GenTauType` 0, 2 or 10
depending on whether the K0 shows up as K0_L, via pi0, or as a K0_S with two extra prongs.
`GenTauTrueMode` separates them.

| Code | Mode | Code | Mode |
| --- | --- | --- | --- |
| `-1` | unknown (label not in the table) | `20` | K |
| `0` | e | `21` | K pi0 |
| `1` | mu | `22` | K ≥2pi0 |
| `10` | pi | `23` | K0 pi |
| `11` | pi pi0 | `24` | K0 pi pi0 |
| `12` | pi 2pi0 | `25` | K pi pi |
| `13` | pi ≥3pi0 | `26` | K K pi |
| `14` | 3pi | `27` | K K0 |
| `15` | 3pi pi0 | `28` | K0 K0 pi |
| `16` | 3pi ≥2pi0 | `29` | K0 3pi |
| `17` | 5pi | `40` / `41` / `42` | omega X / eta X / K* X |

`modules.genDecayModes.mode_name(code)` returns these names. Anything unseen maps to `-1`.

The daughter list itself is kept, flattened over taus:

| Branch | Type | Meaning |
| --- | --- | --- |
| `GenTauNDecayDaughters` | `int` | Per gen tau: how many canonical daughters |
| `GenDecayDaughterTauKey` | `int` | Gen tau this daughter belongs to |
| `GenDecayDaughterPDG` | `int` | Canonical (charge-folded) PDG code |

So an uncatalogued mode can always be inspected without reprocessing:
group `GenDecayDaughterPDG` by `GenDecayDaughterTauKey` where `GenTauTrueMode == -1`.

## Gen tau constituents

One row per final-state, non-neutrino leaf of the gen tau decay tree. Pi0s are expanded
into their two photons.

| Branch | Type | Meaning |
| --- | --- | --- |
| `GenTauConstKey` | `int` | Gen tau this constituent belongs to |
| `GenConstPDG` | `int` | PDG code |
| `GenConstP`, `GenConstTheta`, `GenConstEta`, `GenConstPhi` | `float` | Momentum |
| `GenConstPi0Key` | `int` | Index of the parent pi0 within the tau, `-1` if not from a pi0. Both photons of a pi0 share the value |
| `GenConstMCIdx` | `int` | MC index — join key to `GenPhotonMCIdx` |
| `GenConstOrigin` | `int` | Photon origin, table below (`-1` for non-photons) |

### Neutrinos

Neutrinos are **not** constituents: they are excluded from `GenConst*`, from
`GenTauNConsts` and from the visible 4-momentum, exactly as before. They are stored in a
block of their own, so the visible quantities keep their meaning while the neutrino
kinematics stay available.

| Branch | Type | Meaning |
| --- | --- | --- |
| `GenTauNNus` | `int` | Per gen tau: how many neutrinos (1 hadronic, 2 leptonic) |
| `GenNuTauKey` | `int` | Gen tau this neutrino belongs to |
| `GenNuMCIdx` | `int` | MC index |
| `GenNuPDG` | `int` | Signed PDG code: `+16` for a tau-, `-16` for a tau+, and `-12`/`-14` (resp. `+12`/`+14`) for the lepton antineutrino of a leptonic decay. Not charge-folded, unlike `GenDecayDaughterPDG` — use `abs()` to select by flavour alone |
| `GenNuP`, `GenNuTheta`, `GenNuEta`, `GenNuPhi` | `float` | Momentum |

The point of the block is the leptonic channels: `GenTau* - GenVisTau*` gives only the
**sum** of the two neutrinos, whereas `GenNuPDG` separates the nu_tau from the nu_l.
For the hadronic ones it is just the nu_tau, and the subtraction still works:

```python
tau = ROOT.TLorentzVector(); tau.SetPtEtaPhiM(GenTauPt[i], GenTauEta[i], GenTauPhi[i], GenTauMass[i])
vis = ROOT.TLorentzVector(); vis.SetPtEtaPhiM(GenVisTauPt[i], GenVisTauEta[i], GenVisTauPhi[i], GenVisTauMass[i])
nu  = tau - vis          # equals the sum over GenNu* with GenNuTauKey == i
```

Same collection path as the constituents (`get_visible_final_state`), so simulation
secondaries (`generatorStatus == 0`) are excluded here too. In memory the equivalent is
`GenParticle.getNeutrinos()`, a dict keyed like `getDaughters()` but disjoint from it.

### Photon origin (`GenConstOrigin`, `GenPhotonOrigin`)

| Code | Name | Meaning |
| --- | --- | --- |
| `-1` | not-a-photon | The particle is not a photon |
| `0` | pi0 | From a pi0 decay |
| `1` | tau-FSR | Radiated directly by a tau |
| `2` | charged-rad | Radiated by another charged particle (e, μ, π, K…) |
| `3` | other | Any other generator ancestry (η → γγ lands here) |
| `4` | simulation | Produced by the detector simulation, not the generator |

Event-wide counts in the reference file: charged-rad 2 719, pi0 2 306, tau-FSR 445,
other 2. Restricted to tau constituents the picture flips — a tau's constituent photons
are almost always pi0 daughters.

## Gen extra neutrals

Neutral final-state particles inside the tau decay that the decay-mode encoding does
**not** count (K0_L, neutrons, Λ…). They are already included in the constituents and in
the visible 4-momentum; only the ID ignores them. A `tau → π K0_L ν` therefore still gets
`GenTauType == 0`, and this block is how you find it. Only K0_L (PDG 130) appears in the
reference file, 21 of them.

| Branch | Type | Meaning |
| --- | --- | --- |
| `GenTauHasExtraNeutrals` | `int` | Per gen tau: 1 if it has any |
| `GenTauNExtraNeutrals` | `int` | Per gen tau: how many |
| `GenExtraNeutralTauKey` | `int` | Gen tau this neutral belongs to |
| `GenExtraNeutralMCIdx` | `int` | MC index |
| `GenExtraNeutralPDG` | `int` | PDG code (130 = K0_L) |
| `GenExtraNeutralP`, `…Theta`, `…Eta`, `…Phi` | `float` | Momentum |

## Tau provenance

Most taus come straight from the hard process. A few come from `τ → γ → τ τ`
(radiative splitting) and are flagged here; the reference file contains none.

| Branch | Type | Meaning |
| --- | --- | --- |
| `GenTauMCIdx` | `int` | MC index of the tau |
| `GenTauOriginPDG` | `int` | PDG of the first non-tau ancestor (23 for Z, 22 for a photon) |
| `GenTauIsSecondary` | `int` | 1 when another tau appears above that ancestor |
| `GenTauMotherTauKey` | `int` | Gen tau **key** of the mother tau in this entry, `-1` if primary/unresolved |
| `GenTauMotherTauMCIdx` | `int` | MC index of the mother tau |
| `GenTauRadPhotonMCIdx` | `int` | MC index of the intermediate radiated photon |

`e+e- → γ* → τ τ` is *not* secondary: the ancestry is walked past the first non-tau
ancestor and the flag is raised only if a further tau is found. Pythia tau copy chains are
resolved, so `GenTauMotherTauKey` points at the status-2 tau kept in the tree.

## Reco taus

Hadronic tau candidates first, then electrons, then muons — all in the same block.

| Branch | Type | Meaning |
| --- | --- | --- |
| `RecoTauP`, `RecoTauPt`, `RecoTauEta`, `RecoTauTheta`, `RecoTauPhi`, `RecoTauMass` | `float` | Candidate 4-momentum |
| `RecoTauType` | `int` | Raw reco code, counting **photons** not pi0s |
| `RecoTauDM` | `int` | Photon count folded into a pi0 count, comparable to `GenTauType` |
| `RecoTauQ` | `float` | Charge |
| `RecoTauDR` | `float` | Cone size of the candidate (`getMaxCone`) |
| `RecoTauNConsts` | `int` | Number of constituents |
| `RecoTauNConstKey` | `int` | Own index of the candidate |

`RecoTauType`: `N` = 1 prong + N photons (capped at `9`), `10 + N` = 3 prongs + N photons,
`-11` / `-13` for electron and muon candidates, `-20` for the 1 pion + neutron Pandora
misidentification, `-1` when no valid assignment was possible. `RecoTauDM` maps
`N → ceil(N/2)` and `10 + N → 10 + ceil(N/2)`, i.e. two photons make one pi0, so it can be
histogrammed directly against `GenTauType`.

Constituents, one row per PFO of the candidate:

| Branch | Type | Meaning |
| --- | --- | --- |
| `RecoTauConstKey` | `int` | Reco candidate this constituent belongs to |
| `RecoConstPDG` | `int` | PDG code assigned by the reconstruction |
| `RecoConstP`, `RecoConstTheta`, `RecoConstEta`, `RecoConstPhi` | `float` | Momentum |

## Gen and reco photons

`GenPhoton*` covers **all** generator photons of the event, not only those inside taus.

| Branch | Type | Meaning |
| --- | --- | --- |
| `GenPhotonP`, `…Pt`, `…Eta`, `…Theta`, `…Phi` | `float` | Momentum |
| `GenPhotonMCIdx` | `int` | MC index |
| `GenPhotonTauKey` | `int` | Gen tau it belongs to, `-1` if not a tau constituent |
| `GenPhotonOrigin` | `int` | Origin code (table above) |
| `GenPhotonParentPDG` | `int` | PDG of the direct parent |
| `GenPhotonAncestorMCIdx` | `int` | MC index of the parent — the two photons of a pi0 share it |
| `RecoPhotonP`, `…Pt`, `…Eta`, `…Theta`, `…Phi` | `float` | Reco momentum |
| `RecoPhotonPFOIdx` | `int` | Index in `PandoraPFOs` |
| `RecoPhotonTauKey` | `int` | Reco tau it belongs to, `-1` otherwise |
| `RecoPhotonGenMatchIdx` | `int` | **Position** in the `GenPhoton*` block matched via `RecoMCTruthLink`, `-1` if unmatched |

Note the asymmetry: `RecoPhotonGenMatchIdx` is a position in the gen photon arrays, not an
MC index, so it indexes `GenPhoton*` directly. `RecoPhotonTauKey` is resolved by matching
the rounded momentum vector against the reco tau constituents, so in MLPF/GATr mode — where
reco taus come from the prediction and reco photons from `PandoraPFOs` — it is `-1`
everywhere, which is the intended behaviour.

The truth link is built by `build_truth_links`, controlled by the CLI flags
`--weight-mode` (`decoded` splits the packed track/calo weights, `raw` uses them as is),
`--dedup-mode` (`reco` keeps the best gen per reco, `gen` the other way round),
`--skip-gen-status-filter` and `--max-gen-pdg`.

## Recipes

**Which tau does a gen photon belong to, and where does it come from:**

```cpp
bool from_tau = GenPhotonTauKey[j] >= 0;   // and that same value is the tau key
bool from_pi0 = GenPhotonOrigin[j] == 0;
```

**From a reco photon back to the gen tau:**

```cpp
int g   = RecoPhotonGenMatchIdx[k];        // -1 if no truth link
int tau = (g >= 0) ? GenPhotonTauKey[g] : -1;
```

**Gen–reco migration matrix:** pair `GenTauType[GenMatchedKey[i]]` with
`RecoTauDM[RecoMatchedKey[i]]`, skipping `RecoMatchedKey[i] < 0`.

**Where a visible topology hides something:**

```python
# 1-prong taus that actually carry a K0
r.Filter("Sum(GenTauType == 0 && GenTauTrueMode == 23) > 0")
# or, model-independently, the untagged neutrals
r.Filter("Sum(GenTauType == 0 && GenTauHasExtraNeutrals == 1) > 0")
```

**Constituents of one gen tau** (block scan, since only the key is stored):

```python
mask = ak.Array(a.GenTauConstKey) == i     # i = gen tau key
p    = a.GenConstP[mask]
```

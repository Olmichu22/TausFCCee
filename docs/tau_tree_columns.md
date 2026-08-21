# `Tau_tree` — contents and usage

Output of `TauAnalysis/TTreesTausLong.py`. Reference file (sample `ztt_2M`):

```
Results/TauReco/results0.4_tph0.0_tpi0.0_n0.0_g0.0/Tree_resultsdecayAll_0.4_tph0.0_tpi0.0_n0.0_g0.0.root
```

* Tree name: `Tau_tree` — 2 000 000 entries, 86 branches, ~2.0 GB.
* **One entry per event.** Every physics object lives in a `std::vector` branch, so
  each entry holds a variable number of taus, constituents and photons.
* The `config.yaml` written next to the ROOT file records the exact cuts and sample
  used for that production.

## Quick start

```python
import ROOT
ROOT.EnableImplicitMT(8)
r = ROOT.RDataFrame("Tau_tree", "Tree_resultsdecayAll_0.4_tph0.0_tpi0.0_n0.0_g0.0.root")

# 1-prong + 1 pi0 gen taus, visible momentum
r.Filter("Sum(GenTauType == 1) > 0").Histo1D("GenVisTauP")
```

With `uproot` (no key4hep needed):

```python
import uproot, awkward as ak
t = uproot.open("Tree_....root")["Tau_tree"]
a = t.arrays(["GenTauType", "GenVisTauP", "GenTauHasExtraNeutrals"])
ak.sum(a.GenTauHasExtraNeutrals)
```

## Index conventions

Branches come in blocks; within a block all vectors have the same length and are
aligned element by element. Blocks are linked through explicit key branches.

| Block | Length per entry | Prefix |
| --- | --- | --- |
| Gen taus | `numGenTaus` | `GenTau*`, `GenVisTau*`, `GenEventId` |
| Gen tau constituents | Σ `GenTauNConsts` | `GenConst*`, `GenTauConstKey` |
| Gen extra neutrals | `numGenExtraNeutrals` | `GenExtraNeutral*` |
| Reco taus | `numRecoTaus` | `RecoTau*` |
| Reco tau constituents | Σ `RecoTauNConsts` | `RecoConst*`, `RecoTauConstKey` |
| Gen photons (whole event) | `numGenPhotons` | `GenPhoton*` |
| Reco photons (PandoraPFOs) | `numRecoPhotons` | `RecoPhoton*` |

Two distinct kinds of index appear throughout:

* **Key** (`*Key`) — position inside a block of *this* entry (e.g. `GenTauConstKey[j] = 2`
  means constituent `j` belongs to the third gen tau of the event). `-1` means "no link".
* **MC index** (`*MCIdx`) — `getObjectID().index` of the particle in the `MCParticles`
  collection of the original EDM4hep event. Stable within an event, so it is the safe
  way to join across blocks; `-1` means unavailable.

## Event-level scalars

| Branch | Type | Meaning |
| --- | --- | --- |
| `numGenTaus` | `int` | Gen taus in the entry (`generatorStatus == 2`) |
| `numRecoTaus` | `int` | Reconstructed tau candidates |
| `numGenPhotons` | `int` | Generator photons with `generatorStatus == 1`, whole event |
| `numRecoPhotons` | `int` | Photon PFOs in `PandoraPFOs` |
| `numGenExtraNeutrals` | `int` | Untagged neutrals summed over all gen taus of the entry |
| `beamE` | `double` | Energy of `MCParticles[0]`, i.e. the beam energy |

## Gen taus

| Branch | Type | Meaning |
| --- | --- | --- |
| `GenEventId` | `int` | Global event id, repeated once per gen tau. Unique across the whole file |
| `GenTauP`, `GenTauPt`, `GenTauEta`, `GenTauTheta`, `GenTauPhi`, `GenTauMass` | `float` | Full tau 4-momentum (including neutrinos) |
| `GenVisTauP`, `GenVisTauPt`, `GenVisTauEta`, `GenVisTauTheta`, `GenVisTauPhi`, `GenVisTauMass` | `float` | Visible 4-momentum (neutrinos removed) |
| `GenTauType` | `int` | Decay-mode code, see below |
| `GenTauQ` | `float` | Visible charge |
| `GenTauDR` | `float` | Largest opening angle between constituents |
| `GenTauNConsts` | `int` | Number of final-state constituents of this tau |
| `GenTauNConstKey` | `int` | Own index of the tau (`0, 1, …`); redundant with the position |
| `GenMatchedKey` | `int` | Gen tau index in the gen↔reco matching pair |
| `RecoMatchedKey` | `int` | Reco tau matched to `GenMatchedKey[i]`, or `-1` if none within `MatchedGenMaxDR` |

`GenEventId` encodes provenance: `(file_index * 1000) + event_in_file`, the same key used
by `myutils.get_root_trees_path` for MLPF/GATr predictions. `1000` is the assumed number
of events per input file, and the job aborts if a file exceeds it.

### Decay-mode encoding (`GenTauType`)

| Code | Decay | Count in `ztt_2M` |
| --- | --- | --- |
| `0` | 1 prong, 0 pi0 | 479 425 |
| `1` | 1 prong, 1 pi0 | 1 050 478 |
| `2`, `3`, … | 1 prong, N pi0 | 380 647 / 48 003 / … |
| `10 + N` | 3 prongs, N pi0 | 408 736 (`10`), 195 238 (`11`), … |
| `-11` | tau → e ν ν | 709 365 |
| `-13` | tau → μ ν ν | 694 206 |
| `-2` | other charged content (unclassified) | — |
| `-1` | no valid assignment | — |

`RecoTauType` follows the same idea but counts **photons** instead of pi0s (`N` photons
→ code `N`, capped at `9`; `10 + N` for 3 prongs; `-20` for the 1 pion + neutron Pandora
misidentification). `RecoTauDM` is the photon count folded back into a pi0 count
(`ceil(N/2)`), so it can be compared directly against `GenTauType`.

## Gen tau constituents

One row per final-state particle of a gen tau. Pi0s are expanded into their two photons.

| Branch | Type | Meaning |
| --- | --- | --- |
| `GenTauConstKey` | `int` | Gen tau this constituent belongs to |
| `GenConstPDG` | `int` | PDG code |
| `GenConstP`, `GenConstTheta`, `GenConstEta`, `GenConstPhi` | `float` | Momentum |
| `GenConstPi0Key` | `int` | Index of the parent pi0 within the tau, `-1` if not from a pi0. Both photons of a pi0 share the value |
| `GenConstMCIdx` | `int` | MC index — join key to `GenPhotonMCIdx` |
| `GenConstOrigin` | `int` | Photon origin, see the table below (`-1` for non-photons) |

## Photon origin (`GenConstOrigin`, `GenPhotonOrigin`)

| Code | Name | Meaning |
| --- | --- | --- |
| `-1` | not-a-photon | The particle is not a photon |
| `0` | pi0 | Comes from a pi0 decay |
| `1` | tau-FSR | Radiated directly by a tau |
| `2` | charged-rad | Radiated by another charged particle (e, μ, π, K…) |
| `3` | other | Any other generator ancestry (η → γγ shows up here) |
| `4` | simulation | Produced by the detector simulation, not the generator |

Event-wide counts in `ztt_2M`: charged-rad 5 467 842, pi0 4 422 870, tau-FSR 855 913,
other 7 409. Restricted to tau constituents the picture flips — pi0 4 420 793,
tau-FSR 6 981, charged-rad 1 — since a tau's constituent photons are almost always
pi0 daughters.

## Extra neutrals

Neutral final-state particles that sit inside the tau decay but are **not** counted by
the decay-mode encoding (currently only K(L)0 in this sample: 40 232 of them). They are
already included in `const` and in the visible 4-momentum; only the ID ignores them.
So a `tau → π K0_L ν` still gets `GenTauType == 0`, and this block is how you find it.

| Branch | Type | Meaning |
| --- | --- | --- |
| `GenTauHasExtraNeutrals` | `int` | Per gen tau: 1 if it has any |
| `GenTauNExtraNeutrals` | `int` | Per gen tau: how many |
| `GenExtraNeutralTauKey` | `int` | Gen tau this neutral belongs to |
| `GenExtraNeutralMCIdx` | `int` | MC index |
| `GenExtraNeutralPDG` | `int` | PDG code (130 = K(L)0) |
| `GenExtraNeutralP`, `…Theta`, `…Eta`, `…Phi` | `float` | Momentum |

Particles counted elsewhere — photons, pi0, neutrinos — are never tagged here. Neutrons
are also counted by the reco side but not by the gen ID, so they would appear in this
block; none show up in `ztt_2M`. Filter on `GenExtraNeutralPDG` if you need a subset.

## Tau provenance

Most taus come straight from the hard process. A few (97 out of 4 000 097 in `ztt_2M`)
come from `tau → γ → τ τ`, and those are flagged here.

| Branch | Type | Meaning |
| --- | --- | --- |
| `GenTauMCIdx` | `int` | MC index of the tau |
| `GenTauOriginPDG` | `int` | PDG of the first non-tau ancestor (23 for Z, 22 for a photon) |
| `GenTauIsSecondary` | `int` | 1 when another tau appears above that ancestor |
| `GenTauMotherTauKey` | `int` | Gen tau **key** of the mother tau in this entry, `-1` if unresolved |
| `GenTauMotherTauMCIdx` | `int` | MC index of the mother tau |
| `GenTauRadPhotonMCIdx` | `int` | MC index of the intermediate radiated photon |

`e+e- → γ* → τ τ` is *not* secondary: the ancestry is walked past the first non-tau
ancestor and the flag is only raised if a further tau is found. Pythia tau copy chains
are resolved, so `GenTauMotherTauKey` always points at the status-2 tau kept in the tree
(zero unresolved mothers in `ztt_2M`).

## Gen and reco photons

`GenPhoton*` covers **all** generator photons of the event, not just those inside taus.

| Branch | Type | Meaning |
| --- | --- | --- |
| `GenPhotonP`, `…Pt`, `…Eta`, `…Theta`, `…Phi` | `float` | Momentum |
| `GenPhotonMCIdx` | `int` | MC index |
| `GenPhotonTauKey` | `int` | Gen tau it belongs to, `-1` if it is not a tau constituent |
| `GenPhotonOrigin` | `int` | Origin code (table above) |
| `GenPhotonParentPDG` | `int` | PDG of the direct parent |
| `GenPhotonAncestorMCIdx` | `int` | MC index of the parent — the two photons of a pi0 share it |
| `RecoPhotonP`, `…Pt`, `…Eta`, `…Theta`, `…Phi` | `float` | Reco momentum |
| `RecoPhotonPFOIdx` | `int` | Index in `PandoraPFOs` |
| `RecoPhotonTauKey` | `int` | Reco tau it belongs to, `-1` otherwise |
| `RecoPhotonGenMatchIdx` | `int` | **Position** in the `GenPhoton*` block matched via `RecoMCTruthLink`, `-1` if unmatched |

Note the asymmetry: `RecoPhotonGenMatchIdx` is a position in the gen photon arrays, not
an MC index, so it indexes `GenPhoton*` directly.

## Recipes

**Is this photon from a tau, and from which one?**

```cpp
// per gen photon j
bool from_tau = GenPhotonTauKey[j] >= 0;              // and which tau: that same value
bool from_pi0 = GenPhotonOrigin[j] == 0;
```

**From a reco photon back to the gen tau:**

```cpp
int g = RecoPhotonGenMatchIdx[k];                     // -1 if no truth link
int tau = (g >= 0) ? GenPhotonTauKey[g] : -1;
```

**Which tau does a secondary tau come from:**

```cpp
// per gen tau i
if (GenTauIsSecondary[i]) {
   int mother = GenTauMotherTauKey[i];                // key into the same gen tau block
   // GenTauP[mother], GenTauType[mother], …
}
```

**Taus whose 1-prong decay hides a K(L)0:**

```python
r.Filter("Sum(GenTauType == 0 && GenTauHasExtraNeutrals == 1) > 0")
```

**Gen–reco migration matrix:** pair `GenTauType[GenMatchedKey[i]]` with
`RecoTauDM[RecoMatchedKey[i]]`, skipping `RecoMatchedKey[i] < 0`.

## Consistency checks

Validated on the full 2M-event production: all block lengths coherent,
`GenEventId` unique over the 2 000 000 entries, every tau-constituent photon
resolvable in `GenPhotonMCIdx`, and no secondary tau with an unresolved mother key.
`test/check_gen_flags.py` reruns the provenance statistics on raw EDM4hep files, and
`test/test_gen_provenance.py` covers the classification logic without needing podio.

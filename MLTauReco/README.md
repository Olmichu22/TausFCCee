# ParTauDETR taus

Replaces the cone taus of `modules/tauReco.findAllTaus` with taus reconstructed
by ParTauDETR, a transformer trained on CLD Z->tautau / Z->qq jets at 91 GeV.
The model code comes from the
[ml-tau-model](https://github.com/HEP-KBFI/ml-tau-model) submodule in
`external/ml-tau-model` (branch `TauTr_upgrades`). Its own `ml-tau-data`
submodule supplies the ntupelizer that built the training data, so fetch both
levels:

```bash
git submodule update --init --recursive
```

The model needs torch, lightning, omegaconf and fastjet's Python bindings, which
the key4hep stack does not provide. It therefore runs as a separate first step,
in the ml-tau container (`run_mltau.sh`), and writes its taus to parquet files.
The analysis then reads them inside its usual key4hep event loop.

## 1. Run the model

```bash
export MLTAU_CONTAINER=<apptainer image with torch, lightning, omegaconf, fastjet>
MLTauReco/run_mltau.sh python3 MLTauReco/process_edm4hep.py <EDM4hep .root files> --out-dir <predictions dir>
```

- **What it does:** each event is clustered into jets (ee_genkt, R = 0.4, all
  PandoraPFOs) and every jet goes through the model. Per jet the model
  predicts the tau-ID score p(tau) and the visible tau daughters (pi+- / pi0,
  with 4-momentum and charge). The daughters give the visible tau p4, charge
  and decay mode.
- **Inputs:** only reco collections are read (`PandoraPFOs`,
  `SiTracks_Refitted` track states, `PrimaryVertices`). No MC truth is needed.
- **Speed:** CPU only, about 10 s per 100 events. Spread large samples over
  parallel jobs with `--threads 1`. Existing outputs are skipped unless
  `--overwrite` is given.
- **Model:** taken from `--model-dir`, by default the bundled `MLTauReco/model/`:
  `model.ckpt` (30 MB), the input scaler `scaler/cand_feature_scaler.npz` and
  `bundle.json` (provenance, test performance, tau-ID working points). Another
  directory with the same layout can be passed instead.
- **Container:** `run_mltau.sh` takes the image from `MLTAU_CONTAINER`
  (required). If the input files are outside the directories apptainer mounts
  by default, list their directories in `MLTAU_BINDS` (comma-separated).

The output is one `<input name>_mltaus.parquet` per input file, one row per jet:

| column | content |
|---|---|
| `event_idx` | entry in the input file = position in the `root_io.Reader` loop |
| `event_number`, `jet_idx` | `EventHeader.eventNumber`; jet within the event |
| `jet_pt/eta/phi/energy`, `jet_n_cands` | the reco jet |
| `cand_pfo_idx` | the jet's constituents, as indices into the event's PandoraPFOs |
| `n_reco_leptons` | reco e/mu among the constituents (such jets were not in training) |
| `tau_score` | p(tau): cut on it to select taus |
| `tau_px/py/pz/energy`, `tau_pt/eta/phi/mass` | visible tau = sum of the predicted daughters |
| `tau_charge`, `tau_n_charged`, `tau_n_neutral` | from the daughters |
| `tau_decaymode` | 5 (n_charged - 1) + n_pi0: 0, 1, 2 = h + n pi0, 10, 11 = 3h + n pi0; -1 = no charged daughter |
| `daughter_px/py/pz/energy/charge/meson_class/objectness` | per daughter, pT-ordered; meson_class 0 = pi+-, 1 = pi0 |

Daughters are written for every jet, but they only mean something on jets the
tau ID accepts, because the daughter objectness was trained on tau jets only.

## 2. Use the taus in the analysis

In `TauAnalysis/TTreesTausLong.py`, `--mltau-predictions` replaces
`tauReco.findAllTaus`. Electrons and muons still come from
`electronReco` / `muonReco`.

This step runs in the key4hep stack like the rest of the analysis (ROOT,
podio), not in the ml-tau container. Run it from the repository root:

```bash
source setupKey4Hep.sh
python3 -m TauAnalysis.TTreesTausLong -c config/default/taurecolong_CLD.yaml \
    --input-list <the same EDM4hep files> --prefix MLTau_ \
    --mltau-predictions <predictions dir>
    # optional: --mltau-tau-score-cut 0.9736 (80% eff.)  --mltau-keep-lepton-jets
```

The tree is written to
`Results/TauReco/<prefix><outfile><cuts>/Tree_<outfile>decayAll_<cuts>.root`
(relative to the working directory), e.g.
`Results/TauReco/MLTau_effis0.4_tph0.0_tpi0.0_n0.0_g0.0/`. The prefix only names
the directory: the ParTauDETR taus fill the usual `RecoTau*` branches.

The cone-reconstruction options `--cone-axis` and the photon systematics
(`--sys-err`) do not apply to ParTauDETR taus and are refused together with
`--mltau-predictions`. `--clean-reco-taus` works as for the cone taus.

Only jets passing the tau-score cut (by default the 90%-efficiency working
point, 0.9085) become reco taus in the tree. **To study the effect of the
score on the taus** (ROC curves, other working points, efficiency vs. score),
run with `--mltau-tau-score-cut 0` so that every jet with a predicted daughter
is written, and cut on `RecoTauTag` afterwards. Jets with a reco e/mu are still
dropped unless `--mltau-keep-lepton-jets` is also given.

In your own event loop, use `modules/mlTauReco.py`:

```python
from modules import mlTauReco
mltaus = mlTauReco.MLTauReader("<predictions dir>")
for local_event, event in enumerate(root_io.Reader([filename]).get("events")):
    recoTau_raw = mltaus.findAllTaus(filename, local_event)   # {i: RecoParticle}
```

How these taus differ from the cone taus:

- **ID:** `getID()` uses the cone's photon-counting convention, with each
  predicted pi0 counted as two photons: 2 n_pi0 for 1 prong, 10 + 2 n_pi0 for
  3 prongs. So `RecoTauDM = ceil(ID / 2)` and `MatchRecoGenTau` give the pi0
  count. `tau.n_pi0` holds it directly.
- **Constituents:** `getDaughters()` returns the predicted pi+- / pi0
  (PDG +-211 / 111), not PFOs. That means no `RecoPhotonTauKey` links, and
  `extra_reco_correction` does not apply. The jet's PFOs are in
  `tau.cand_pfo_idx`.
- **Cone size:** `getMaxCone()` (`RecoTauDR`) is the largest 3D angle between
  a predicted daughter and the visible tau, not the theta-phi distance
  (`myutils.dRAngle`) to the seed pion that the cone taus use. Do not compare
  the two directly.
- **Default selection:**
  - `tau_score` >= the 90%-efficiency working point (0.9085).
  - Jets with a reco e/mu are dropped (`veto_lepton_jets`).
- **Extra attributes:** `tau.tau_score`, `tau.jet_idx`, `tau.jet_p4`.

Extra branches in the `TTreesTausLong` tree, written only with
`--mltau-predictions` (-1, or -999 for eta / phi, for the electron / muon
entries):

| branch | content |
|---|---|
| `RecoTauTag` | tau-ID score p(tau) |
| `RecoTauJetIdx`, `RecoTauJetPt/Eta/Phi/E` | the jet the tau was predicted in |
| `RecoTauJetNConsts` | PFOs in that jet |
| `RecoConstObjectness` | per predicted daughter (`RecoConst*`, keyed by `RecoTauConstKey`) |
| `RecoConstCharge` | per predicted daughter |

The predicted daughters themselves are the `RecoConst*` entries: P, theta,
eta, phi, PDG (+-211 / 111), charge. Their mass is the pi+- / pi0 mass.

## Model and performance

The bundled model is the U1-log run of ml-tau-model's TauTr_upgrades study:
20k steps, `significance_transform=log`. On the training test jets (hadronic
taus, gen-matched jets without reco e/mu):

| metric | value |
|---|---|
| visible pT resolution (IQR/median) | 3.43% |
| median dR to the gen visible tau | 2.2e-3 |
| decay-mode accuracy | 0.915 |
| charge accuracy | 0.953 |
| tau-ID AUC vs Z->qq jets | 0.991 |

Tau-ID working points:

| tau efficiency | tau_score cut | qq-jet misID |
|---|---|---|
| 50% | 0.9965 | 2.8e-4 |
| 70% | 0.9888 | 7.9e-4 |
| 80% | 0.9736 | 1.8e-3 |
| 90% | 0.9085 | 6.2e-3 |
| 95% | 0.7257 | 2.1e-2 |

### ParTauDETR vs. cone taus in the analysis

Both algorithms were run through `TTreesTausLong` on the same 25k Z->tautau
events (`reco_p8_ee_Z_tautau_ecm91_700000` - `700249`) and 25k Z->qq events
(`reco_p8_ee_Z_qq_ecm91_600000` - `600249`), ParTauDETR at the default 90%
working point. Gen taus, gen-reco matching, decay mode and charge are this
repo's definitions, evaluated with `compare_with_cone.py` over all 32,619
hadronic gen taus. To redo it on other trees (in the ml-tau container, since it
uses ml-tau-model's evaluation tools):

```bash
MLTauReco/run_mltau.sh python3 MLTauReco/compare_with_cone.py \
    --trees Cone:z=<dir> Cone:qq=<dir> ParTauDETR:z=<dir> ParTauDETR:qq=<dir> \
    --out-dir <plots dir>
```

| | ParTauDETR | cone (`tauReco`) |
|---|---|---|
| matched to any reco tau | 0.901 | 0.861 |
| matched to a hadronic reco tau | 0.821 | 0.648 |
| decay-mode accuracy | 0.855 | 0.679 |
| decay-mode accuracy, hadronic reco taus only | 0.937 | 0.900 |
| charge accuracy | 0.946 | 0.960 |
| pT response (median reco / gen) | 1.005 | 1.002 |
| pT resolution (IQR / median) | 3.9% | 3.7% |
| median dR to the gen visible tau | 1.9e-3 | 4.5e-3 |
| hadronic reco taus per Z->qq event (fakes) | 0.056 | 2.37 |

Decay-mode F1 per class:

| | h | h+pi0 | h+>=2pi0 | 3h | 3h+>=pi0 |
|---|---|---|---|---|---|
| ParTauDETR | 0.928 | 0.904 | 0.852 | 0.899 | 0.838 |
| cone | 0.860 | 0.805 | 0.736 | 0.668 | 0.619 |

- **ParTauDETR is better in** decay mode (every class, most for 3 prongs),
  efficiency, direction, and fake rate (about 40x fewer in Z->qq).
- **The cone is slightly better in** charge accuracy and pT resolution, but
  only overall. On the 26,591 gen taus both algorithms match, they are equal:
  charge 0.963 (ParTauDETR) vs. 0.965 (cone), pT resolution 3.46% vs. 3.52%.
  The difference comes from the harder taus that only ParTauDETR finds.

Reco taus per Z->tautau event. Each event has 2 gen taus, of which on average
1.305 decay hadronically:

| | ParTauDETR | cone |
|---|---|---|
| hadronic reco taus / event | 1.103 | 0.910 |
| ... matched to a hadronic gen tau | 1.071 | 0.845 |
| ... not matched (fakes) | 0.032 | 0.064 |

So ParTauDETR finds 82% of the hadronic taus (1.071 / 1.305), with half the
cone's fake rate. Events with 0 / 1 / 2 / 3 ParTauDETR taus: 5046 / 12329 /
7618 / 7.

These numbers are lower than those on the training test jets above because
the selection differs: here every hadronic gen tau counts, with this repo's
angular matching, instead of only gen-matched jets without reco e/mu.

## Validation

`process_edm4hep` was checked jet by jet against the model's training test
ntuples (rebuilding the same raw-file jets and comparing model inputs and
outputs):

- **Z->qq:** identical.
- **Z->tautau:** agrees within the ntuples' float32 rounding of the
  candidates' eta / phi. The tau score moves by at most 1e-3.
- **Undefined IP sign:** in single-track jets the impact-parameter sign is
  undefined (the PCA is perpendicular to the jet axis), so the two can differ
  there in sign only.

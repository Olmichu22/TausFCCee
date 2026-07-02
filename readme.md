# Tau Reconstruction at FCC-ee
This repository contains tools and scripts for reconstructing tau leptons at the Future Circular Collider (FCC-ee) using the EDM4hep data format. The focus is on validating tau reconstruction algorithms (such as PandoraPFO).

The main tools included are:
- Identification of the gennerator-level tau leptons and their decay products (decay type).
- Reconstruction of hadronic and leptonic tau decays from PandoraPFO collections.
- Matching reconstructed taus to generator-level taus.

# Simple Tau Reconstruction Test

This script (`test/simple_taureco_test.py`) provides a **basic example** for reconstructing **tau leptons** and related decay products from ROOT event data using the **EDM4hep** format.

It reads generator-level (`MCParticles`) and reconstructed-level (`PandoraPFOs`) collections, performs tau reconstruction, fills histograms for analysis, and writes the results to a ROOT output file.

The script performs the following tasks:

* Reconstruct **tau leptons** from **hadronic** and **leptonic** (**electrons** and **muons**) decays from PFO collections.
* Match reconstructed taus to their generator-level counterparts.
* Produce validation histograms (momentum, mass, decay type, and P resolution).
* Save the output in a ROOT file for further analysis.

---

# Tree-Based Tau Reconstruction Workflow

An alternative two-step workflow stores all reconstruction results in a **TTree** first and then derives histograms from it. This separates the (slow) EDM4hep loop from the (fast) histogram filling and allows re-running the analysis without reprocessing the raw data.

## Step 1 — Build the TTree (`test/test_tau_reco_tree.py`)

Reads EDM4hep input files, runs the same tau reconstruction and gen–reco matching as `simple_taureco_test.py`, and writes a ROOT file containing a `Tau_tree` with one entry per event.

Run with:

```bash
python test/test_tau_reco_tree.py [options]
```

Output: `Results/TauReco/Tree_<fileOutName>.root`

### Key branches

| Branch            | Type            | Description                                          |
| ----------------- | --------------- | ---------------------------------------------------- |
| `numGenTaus`      | `int`           | Number of gen-level taus in the event                |
| `GenTauPt`        | `vector<float>` | Gen tau transverse momentum                          |
| `GenVisTauPt`     | `vector<float>` | Visible tau transverse momentum                      |
| `GenTauP`         | `vector<float>` | Gen tau total momentum                               |
| `GenVisTauP`      | `vector<float>` | Visible tau total momentum                           |
| `GenTauMass`      | `vector<float>` | Full gen tau invariant mass (~1.777 GeV)             |
| `GenVisTauMass`   | `vector<float>` | Visible tau invariant mass                           |
| `GenTauType`      | `vector<int>`   | Decay mode ID (see decay ID table below)             |
| `RecoMatchedKey`  | `vector<int>`   | Index of matched reco tau for each gen tau (-1 if unmatched) |
| `RecoTauPt`       | `vector<float>` | Reco tau transverse momentum                         |
| `RecoTauP`        | `vector<float>` | Reco tau total momentum                              |
| `RecoTauMass`     | `vector<float>` | Reco tau invariant mass                              |
| `RecoTauDM`       | `vector<int>`   | Reco tau decay mode ID                               |

## Step 2 — Fill histograms (`test/hist_from_tree.py`)

Reads the TTree produced in Step 1 and fills the same histograms as `simple_taureco_test.py`. The momentum resolution is computed using `RecoMatchedKey[i]` directly from the tree (no re-running of the matching algorithm).

Run with:

```bash
python test/hist_from_tree.py [options]
```

Output: `Results/TauReco/Hist_<fileOutName>.root`

---

## Requirements

To set up the environment, ensure you have access to cvmfs. Load the environment with:

```bash
source setupKey4Hep.sh
```

---

## Usage

Run the script with:

```bash
python test/simple_tau_reco_test.py [options]
```

### Optional Arguments

| Argument          | Description                                         | Default  |
| ----------------- | --------------------------------------------------- | -------- |
| `-v`, `--verbose` | Increase verbosity (`-v` for INFO, `-vv` for DEBUG) | 0        |
| `--log-dir`       | Directory to save log files                         | `./logs` |

### Example

```bash
python test/simple_tau_reco_test.py -v --log-dir=./mylogs
```

---

## Output

All results are saved in:

```
Results/TauRecoTestHist.root
```

This ROOT file includes histograms such as:

| Histogram            | Description                                      |
| -------------------- | ------------------------------------------------ |
| `histoGenTauPt`      | Generator-level tau transverse momentum          |
| `histoGenTauVisPt`   | Visible tau transverse momentum                  |
| `histoGenTauP`       | Gen tau total momentum                           |
| `histoGenTauVisP`    | Visible tau total momentum                       |
| `histoGenTauMass`    | True tau invariant mass (~1.777 GeV)             |
| `histoGenTauVisMass` | Visible tau invariant mass                       |
| `histoGenTauType`    | Gen tau decay mode ID                            |
| `histoRecoTauPt`     | Reconstructed tau transverse momentum            |
| `histoRecoTauP`      | Reconstructed tau total momentum                 |
| `histoRecoTauType`   | Reco tau decay mode (hadronic/leptonic)          |
| `histoRecoTauMass`   | Reconstructed tau invariant mass                 |
| `histoResTauP`       | Momentum resolution ((Reco − GenVis) / GenVis)  |

> **Note:** `hist_from_tree.py` produces the same set of histograms.


### Description of decay ID:
Decay types are identified by integers:
- `-11`: tau → e $\nu$
- `-13`: tau → $\mu\nu$
For hadronic decays:
- One-prong decays:
  - `0` + number of $π^0$ (e.g., `1` for $\pi^{\pm}\pi^0\nu$, `2` for $\pi^{\pm}2\pi^0\nu$ ($\rho$ decay), etc.)
- Three-prong decays:
  - `10` + number of $π^0$ (e.g., `10` for $3\pi^{\pm}\nu$, `11` for $3\pi^{\pm}\pi^0\nu$ ($a_1$ decay), etc.)

- `-1` indicates non recognized tau decay (for example, if charge is not conserved).
- `-20`is for reco taus with reco neutral hadrons (known issue with PandoraPFO).

---

## Main Modules and Functions

### **1. `myutils`**

Utilities for I/O and ROOT file management.

* `get_root_trees_path(path, file_prefix, ...)`
  Reads multiple ROOT files with a common prefix.
* `get_single_root_file(file_path, loggers)`
  Reads a single ROOT file for analysis.

### **2. `tauReco`**

Main reconstruction logic for hadronic tau decays.

* `findAllGenTaus(mc_particles)`
  Identifies generator-level tau leptons.
* `findAllTaus(pfos, dRMax, minPTauPhoton, ...)`
  Reconstructs hadronic tau candidates from PFOs.
* `MatchRecoGenTau(gen_tau, reco_taus, ..., maxDRMatch, selectDecay)`
  Matches generator-level and reconstructed taus using angular distance.

### **3. `electronReco`**

Reconstructs tau → e decays.

* `findAllElectrons(pfos, generalPCut)`
  Finds reconstructed electrons in the event.

### **4. `muonReco`**

Reconstructs tau → μ decays.

* `findAllMuons(pfos, generalPCut)`
  Finds reconstructed muons in the event.


## Workflow Summary

1. **Parse configuration and logging options.**
2. **Prepare histograms** for generator and reconstructed variables.
3. **Read input ROOT files** via `podio.root_io.Reader`.
4. **Loop through events**:

   * Extract MCParticles and PandoraPFOs.
   * Reconstruct taus, electrons, and muons.
   * Match reconstructed and generator-level taus.
   * Fill histograms accordingly.
5. **Write all histograms** to the output ROOT file.
6. **Log execution summary.**

---

## Output Example

After execution, expect a summary like:

```
INFO [io] - Read 1 files
INFO [processing] - Processing event 0
INFO [processing] - Found 2 gen taus and 2 reco taus
INFO [io] - Output file Results/TauRecoTestHist.root
INFO [io] - End of job
```

---

## Configuration Parameters

| Parameter                                  | Description                  | Default    |
| ------------------------------------------ | ---------------------------- | ---------- |
| `dRMax`                                    | Max ΔR for tau constituents  | 0.4        |
| `minPTauPhoton`, `minPTauPion`, `PNeutron` | Minimum momentum cuts        | 0.0        |
| `dRMatch`                                  | Max ΔR for gen–reco matching | 1.0        |
| `generalPCut`                              | General minimum momentum cut | 0.0        |
| `outputpath`                               | Directory for ROOT output    | `Results/` |

---

# Particle-Level Hit Analysis (parallel)

`HitAnalysis/particle_level_analisis_parallel.py` measures **particle-level
reconstruction performance**: it matches every generator particle to a
reconstructed one (PandoraPFO or MLPF/GATr), and produces **PID confusion
matrices**, **reconstruction-efficiency curves**, and **energy-resolution
plots** — split by energy bin and aggregated. It is the parallel (multi-worker)
version of `HitAnalysis/particle_level_analisis.py` and is the recommended way
to run over more than a handful of files.

The matching is done **two independent ways**, and every output is produced for
both:

- **`dR`** — geometric matching by angular distance between gen and reco.
- **`truthlink`** — matching via the `RecoMCTruthLink` collection (MC-truth
  energy-deposit weights). This branch is empty if the input files do not carry
  a readable `RecoMCTruthLink` (depends on the podio version that wrote them).

---

## 1. Environment

The script needs the Key4hep stack (ROOT, podio, edm4hep). **For ILD work use:**

```bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08
```

If that release misbehaves, fall back to the nightlies:

```bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
```

All commands below are run from the repository root.

---

## 2. Requirement: enable the truth-link association in the config

The confusion matrices are built from a per-event association table that is only
produced when the config enables the hit-type map. A ready-made config is
provided:

```yaml
# config/default/taurecolong_optimal_neutral_reco.yaml
neutral_recover:
  enable: true
  return_hit_type_map: true   # <-- required, otherwise the output is empty
```

Always pass it with `-c` (see below). Without `return_hit_type_map: true` the
run finishes cleanly but writes **no** matrices.

---

## 3. Choosing the input — two options

### Option A — by sample name (`config/samples/samples.yaml`)

Give a **sample name** (or alias) defined in `config/samples/samples.yaml`; the
script resolves it to the directory on disk and reads every `.root` there.

```bash
python HitAnalysis/particle_level_analisis_parallel.py \
    -c config/default/taurecolong_optimal_neutral_reco.yaml \
    -f ztt \
    --n-workers 8
```

To register a new sample, add an entry to `config/samples/samples.yaml`:

```yaml
default_base: /pnfs/ciemat.es/data/cms/store/user/cepeda/FCC/FullSim/
default_file_prefix: out_reco_edm4hep_edm4hep

samples:
  my_sample:                        # name passed to -f
    folder: My_Sample_Directory     # relative to default_base
    aliases: [mine]                 # optional alternative names
  # other forms: `path:` (absolute dir), `folders:`/`paths:` (several dirs merged)
```

The sample can also be fixed in the config under `general.sample`; the `-f` flag
overrides it.

### Option B — by explicit file list (`--input-list`)

Skip `samples.yaml` entirely and pass one or more absolute ROOT paths. Handy for
a quick test on a single file:

```bash
python HitAnalysis/particle_level_analisis_parallel.py \
    -c config/default/taurecolong_optimal_neutral_reco.yaml \
    --input-list /pnfs/.../out_reco_edm4hep_edm4hep_1.root \
    --n-workers 4
```

`--input-list` accepts several files separated by spaces and takes precedence
over the sample name.

---

## 4. Key command-line options

| Option                     | Default                                         | Description                                                                 |
| -------------------------- | ----------------------------------------------- | --------------------------------------------------------------------------- |
| `-c`, `--config`           | `config/default/taurecolong.yaml`               | Analysis config. Use the `_optimal_neutral_reco` one to get matrices.       |
| `-f`, `--sample`           | from config `general.sample`                    | Sample name/alias in `config/samples/samples.yaml`.                         |
| `--input-list FILE [...]`  | `None`                                           | Explicit ROOT file(s); bypasses `samples.yaml`.                             |
| `--n-workers`              | all CPU cores                                    | Number of parallel workers.                                                 |
| `--gatr-result PATH`       | `None`                                           | Use MLPF/GATr predictions instead of PandoraPFOs.                          |
| `--dedup-mode {gen,reco}`  | `reco`                                           | Deduplication side for `RecoMCTruthLink` matching.                          |
| `--weight-mode {raw,decoded}` | `decoded`                                     | How to interpret the truth-link weight (decoded splits track/cluster).      |
| `--max-gen-pdg N`          | `10000`                                          | Ignore gen particles with `|PDG| > N` (excludes nuclear fragments).         |
| `--skip-gen-status-filter` | off                                              | Keep non-final (`generatorStatus != 1`) gen particles.                      |
| `--all-plot PDG [...]`     | `22`                                             | PDGs for the combined energy-resolution plots.                              |
| `-v`, `-vv`                | warnings only                                    | Increase log verbosity (INFO / DEBUG).                                       |

---

## 5. Output structure

Everything is written under `Results/TauReco/<auto-named-by-cuts>/` (the
sub-directory name encodes the cut values). Each result type has a `dR/` and a
`truthlink/` sub-folder:

```
Results/TauReco/<run>/
├── association_results_full_dR.csv         # one row per gen–reco match (dR)
├── association_results_full_truthlink.csv  # one row per gen–reco match (truth-link)
├── config.yaml                             # snapshot of the config used
├── worker_*.log                            # per-worker logs
│
├── confusion_matrices_particle_level/
│   ├── dR/                                  # and truthlink/
│   │   ├── confusion_matrix_bin_00_absolute.png    # per energy bin:
│   │   ├── confusion_matrix_bin_00_efficiency.png  #   abs. counts / efficiency
│   │   ├── confusion_matrix_bin_00_purity.png      #   / purity
│   │   ├── ...                                      # bins 00..07
│   │   ├── confusion_matrix_general_absolute.png   # ALL bins combined:
│   │   ├── confusion_matrix_general_efficiency.png #   abs / efficiency / purity
│   │   ├── confusion_matrix_general_purity.png
│   │   └── confusion_matrix_all_bins.pdf           # every matrix in one PDF
│   └── truthlink/ ...
│
├── efficiency_plots/                        # efficiency vs |p_gen|
│   ├── dR/ , truthlink/
│   │   ├── efficiency_<genpid>_<recopid>.png       # one gen→reco pair
│   │   └── efficiency_global_<genpid>.png          # all destinations + total
│
├── efficiency_plots_theta/                  # efficiency vs θ_gen (same layout)
│
└── energy_distributions/                    # (reco−true)/true resolution vs E
    ├── dR/ , truthlink/
    │   ├── residual_resolution_std_<mig>.png
    │   ├── residual_resolution_iqr84_16_<mig>.png
    │   ├── residual_resolution_std90_<mig>.png
    │   └── residual_resolution_combined_<mig>.png
```

How to read the matrices:
- **rows = gen particle, columns = reco particle** (PDG labels).
- **Efficiency** = matrix normalised by row (gen) → fraction of each gen species
  reconstructed as each reco species.
- **Purity** = matrix normalised by column (reco) → fraction of each reco species
  that truly came from each gen species.
- `general` = the same three matrices summed over all energy bins; the `bin_NN`
  files break it down by the gen-energy bins `[0,1,5,10,20,30,45,100,∞]` GeV.

The `_dR.csv` / `_truthlink.csv` tables hold the raw matches (`Gen_pid`,
`Reco_pid`, energies, momenta, `event_id`, and `dR`) if you want to re-plot or
cross-check the aggregated figures yourself.

---

## 6. Quick start (copy-paste)

```bash
# 1. Environment (ILD)
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08

# 2. Run on a single file to try it out
python HitAnalysis/particle_level_analisis_parallel.py \
    -c config/default/taurecolong_optimal_neutral_reco.yaml \
    --input-list /pnfs/.../out_reco_edm4hep_edm4hep_1.root \
    --n-workers 4

# 3. Inspect the results
ls Results/TauReco/*/confusion_matrices_particle_level/dR/
```

> The single-worker, non-parallel `HitAnalysis/particle_level_analisis.py`
> produces the same `dR` confusion matrices (written directly under
> `confusion_matrices_particle_level/`, without the `dR`/`truthlink` split) and
> is useful for debugging, but the parallel version is preferred for real runs.

---

## License

This project is intended for research and educational use in particle reconstruction studies.
Please cite appropriately if used in your work.

---

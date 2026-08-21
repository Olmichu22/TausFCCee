# Tau Reconstruction at FCC-ee
This repository contains tools and scripts for reconstructing tau leptons at the Future Circular Collider (FCC-ee) using the EDM4hep data format. The focus is on validating tau reconstruction algorithms (such as PandoraPFO).

The main tools included are:
- Identification of the gennerator-level tau leptons and their decay products (decay type).
- Reconstruction of hadronic and leptonic tau decays from PandoraPFO collections.
- Matching reconstructed taus to generator-level taus.
- **Tau polarization analysis** (`RhoAnalysis/`): extraction of the polarization
  asymmetry `A_τ` from the optimal polarimeter observable, with event reweighting,
  background handling, selection-cut optimization and an MLP observable.
  See *Physics of Tau Polarization at the Z Pole* below for the theory, formulas
  and their mapping to the code.

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
├── association_results_full_dR.*           # one row per gen–reco match (dR)
├── association_results_full_truthlink.*    # one row per gen–reco match (truth-link)
│                                           #   .parquet if pyarrow/fastparquet is
│                                           #   available, else .pkl.gz (read back
│                                           #   with pd.read_parquet / pd.read_pickle)
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

# Physics of Tau Polarization at the Z Pole

Theory reference for the `RhoAnalysis/` pipeline. Everything below follows
J. Alcaraz / FCC-CIEMAT team, *Differential cross section distributions for tau
polarization at the Z pole*, 2 June 2026 (`Reweight_tautau.pdf`); equation
numbers in brackets are the ones in that note. The implementation lives in
`modules/weightsPol.py` and `modules/optimalVariabRho.py` — see the
formula→code map in §7.

## 1. Differential cross section

In the $m_\tau/E_\tau \to 0$ limit, for $e^+e^- \to \tau^+\tau^-$ with both taus
decaying to a pseudoscalar or charged vector resonance plus one neutrino
($\tau^-\to R^-\nu_\tau$, $\tau^+\to R'^+\bar\nu_\tau$):

$$\frac{d\sigma}{dz\,dz_R\,dz_{R'}} = \left[\frac{3}{8}(1+z^2) + A_{FB}\,z\right] F(z, z_R, z_{R'}) \tag{1}$$

$$F(z, z_R, z_{R'}) = 1 + \mathcal{P}(z)_\tau\,\big(H_R(z_R) + H_{R'}(z_{R'})\big) + H_R(z_R)\,H_{R'}(z_{R'}) \tag{2}$$

The angular variables — **the definitions matter, see §6**:

| Variable | Definition |
| --- | --- |
| $z \equiv \cos\theta$ | Polar angle of the **$\tau^-$**, positive Z axis along the direction of the **colliding electron** |
| $z_R \equiv \cos\theta_R^*$ | Polar angle of $R^-$ in the **$\tau^-$ rest frame**, using the $\tau^-$ flight direction as $+Z$ |
| $z_{R'} \equiv \cos\theta_{R'}^*$ | Polar angle of $R'^+$ in the **$\tau^+$ rest frame**, using the $\tau^+$ flight direction as $+Z$ |

with $A_{FB} = \tfrac{3}{4}\mathcal{A}_e\mathcal{A}_\tau$ the forward-backward
charge asymmetry. Each spin analyzer is therefore computed in **its own tau's
rest frame**, which is what makes the formula charge-symmetric in the decay part.

Eq. (2) can be rewritten as a helicity decomposition [7–8]:

$$F = \left(\frac{1+\mathcal{P}(z)_\tau}{2}\right)(1+H_R)(1+H_{R'}) + \left(\frac{1-\mathcal{P}(z)_\tau}{2}\right)(1-H_R)(1-H_{R'})$$

i.e. at a given $z$ the process splits into two subprocesses with probabilities
$(1\pm\mathcal{P})/2$ in which the $\tau^-$ has positive / negative helicity (and
the $\tau^+$ the opposite). Within this approximation the two decays can be
implemented independently.

## 2. The polarization $\mathcal{P}(z)_\tau$

$\mathcal{P}(z)_\tau$ is the polarization **of the $\tau^-$** at $z=\cos\theta$.
At the Z peak, ignoring $\gamma^*$ exchange:

$$\mathcal{P}(z)_\tau \approx -\frac{\mathcal{A}_\tau(1+z^2) + 2\mathcal{A}_e z}{(1+z^2) + 2\mathcal{A}_e\mathcal{A}_\tau z} \tag{3}$$

The lepton asymmetry parameters follow from the effective couplings,
$\mathcal{A}_\ell = 2(g_V^\ell/g_A^\ell)/\big(1+(g_V^\ell/g_A^\ell)^2\big)$ with
$g_V^\ell/g_A^\ell = 1 - 4\sin^2\theta_\text{eff}$. This is the single entry point
for $\sin^2\theta_\text{eff}$ in the whole chain (`--sin-eff`, default `0.2312`).

**The structure of eq. (3) drives the entire measurement strategy:** the
$\mathcal{A}_\tau(1+z^2)$ term is **even** in $z$, while $2\mathcal{A}_e z$ and
$2\mathcal{A}_e\mathcal{A}_\tau z$ are **odd**. So $\mathcal{A}_\tau$ is accessible
from the $z$-integrated polarization, whereas $\mathcal{A}_e$ lives *entirely* in
the odd modulation — which is why the $\cos\theta$ binning
(`makeCosBins_MDecs.py`) exists at all, and why the sign convention of §6 is not
negotiable.

## 3. Spin analyzers $H$ per channel

**Hadronic** — $H_R$ depends on whether $R$ is a pseudoscalar ($\pi$, $K$) or a
vector resonance $V$ of mass $m_V$ ($\rho$, $a_1$):

$$H_\pi(z_\pi) = z_\pi \tag{4}$$

$$H_V(z_V) = \alpha_V\, z_V = \left(\frac{m_\tau^2 - 2m_V^2}{m_\tau^2 + 2m_V^2}\right) z_V \tag{5}$$

The polar angle is measurable from $x = E_R/E_\tau$, the fraction of the tau
energy carried by $R$ in the laboratory system:

$$z_R \equiv \cos\theta_R^* = \frac{2x - 1 - \xi}{1 - \xi}, \qquad \xi = m_R^2/m_\tau^2 \tag{6}$$

The dilution factor $\alpha_V$ is what makes the vector channels much less
sensitive than the pion: $\alpha_\rho \approx 0.46$, and with the PDG $a_1$ pole
mass $\alpha_{a_1} \approx 0.021$ (nearly blind).

**Leptonic** — replacing $R'$ by a purely leptonic decay, $z_{R'}$ is substituted
by $x_\ell$, the tau energy fraction taken by the charged lepton in the lab:

$$F(z, z_R, x_\ell) = f(x_\ell)\left[1 + \mathcal{P}(z)_\tau (H_R(z_R) + H_\ell(x_\ell)) + H_R(z_R) H_\ell(x_\ell)\right] \tag{10}$$

$$f(x_\ell) = \tfrac{1}{3}(5 - 9x_\ell^2 + 4x_\ell^3), \qquad g(x_\ell) = \tfrac{1}{3}(1 - 9x_\ell^2 + 8x_\ell^3) \tag{11,12}$$

$$H_\ell(x_\ell) = \frac{g(x_\ell)}{f(x_\ell)} = \frac{1 + x_\ell - 8x_\ell^2}{5 + 5x_\ell - 4x_\ell^2} \tag{13}$$

The right-hand factored form (common $(1-x)$ cancelled) is the one implemented:
it is numerically stable at $x\to1$, where the unfactored ratio is $0/0$ and the
limit is $H_\ell(1) = -1$.

## 4. The optimal variable $\omega$

The normalized differential decay distribution of the $\tau$ can always be
expressed via an **optimal variable** $\omega$ that absorbs all the information
from the full decay chain (Davier, Duflot, Le Diberder, Rougé, 1993):

$$\frac{1}{N}\frac{dN}{d\omega} = f(\omega)\,(1 + \mathcal{P}\,\omega) \tag{15}$$

For the $\rho$ this is strictly more powerful than $H_V = \alpha_V z_V$, because
$\omega$ encodes the full $\rho \to \pi\pi^0$ substructure ($\cos\beta$,
$\cos\psi$) rather than the energy fraction alone. It is computed by
`optimalVariabRho.wVariab` and is the **default** for the $\rho$ channel;
`--no-omega-weights` falls back to eq. (5).

Note $\omega$ plays the role of $H$: for the single pion the two coincide
($\alpha_\pi = 1$), and gen-level code may use the exact boosted angle
(`cosThetaStar`) instead of the analytic $z_R$ reconstruction.

## 5. Event reweighting

To move from the reference scenario $(\mathcal{A}_e, \mathcal{A}_\tau)$ to an
alternative $(\mathcal{A}'_e, \mathcal{A}'_\tau)$, the event weight is the ratio
of eq. (2) evaluated with the new and old polarizations:

$$\mathcal{W} = \frac{1 + \mathcal{P}'(z)_\tau (H_R + H_{R'}) + H_R H_{R'}}{1 + \mathcal{P}(z)_\tau (H_R + H_{R'}) + H_R H_{R'}} \tag{9}$$

with the had–lep [14] and fully general optimal-variable [16] versions following
the same pattern ($H\to\omega$, $H'\to\omega'$). A **per-tau** (single-hemisphere)
version drops the partner terms, $\mathcal{W} = (1+\mathcal{P}'H)/(1+\mathcal{P}H)$.

Two properties worth internalizing:

- **The weight does not "add" polarization — it divides by the SM density with
  which the events were generated.** With $\mathcal{A}'_\tau = \pm1$ the numerator
  collapses to $(1\mp H_R)(1\mp H_{R'})$ exactly, independent of $z$, so *all* the
  $z$ dependence sits in the denominator. A wrong $z$ means dividing by a density
  that did not generate the sample. This is what the templates `P1` / `M1`
  (`weight_P1`, `weight_M1`) are.
- **The angular factor of eq. (1) is not included** in eq. (9): the
  $[\tfrac{3}{8}(1+z^2) + A_{FB}z]$ terms depend on $A_{FB}$ and hence on
  $\mathcal{A}_e, \mathcal{A}_\tau$ themselves. For a combined analysis that also
  uses the observed charge asymmetry, those factors must be reinstated in both
  numerator and denominator.

## 6. Sign convention for $z$ — mandatory

$z$ is **always the $\cos\theta$ of the $\tau^-$**, never of the hemisphere being
analyzed. Since the taus are back-to-back,
$\cos\theta_{\tau^+} = -\cos\theta_{\tau^-}$, and by §2 the odd terms of eq. (3)
flip: $z \to -z \iff \mathcal{A}_e \to -\mathcal{A}_e$. Using the hemisphere's own
angle therefore inverts $\mathcal{A}_e$ in ~50% of events, destroying exactly the
odd information the measurement is after. $\mathcal{A}_\tau$ (even part) survives.

What **does not** carry a charge sign: the spin analyzers. $\omega$, $z_R$ and
$H_\ell$ are each computed in their own tau's rest frame (eq. 2 definitions
above), which is natural for both charges — the antineutrino handedness of the
$\tau^+$ decay is already absorbed, since $d\Gamma(\tau^+) \propto 1 - h^+\omega^+
= 1 + h^-\omega^+$. Both hemispheres are governed by the *same*
$\mathcal{P}(z)_\tau$ of the $\tau^-$.

In this repository the $+Z$ axis convention is satisfied by construction: the
Pythia cards in `../SimProdScripts/CLDFCC_sim/` set `Beams:idA = 11` (electron)
with no frame overrides, and Pythia8 sends beam A along $+z$. (The CLD sim
applies a 15 mrad `crossingAngleBoost`, a known second-order tilt that does not
affect the sign.)

The charge sign is carried through the tree branches `tauPDG` (gen) and
`recoCharge` (reco). For background on why this convention is enforced explicitly
rather than left implicit, see
[`docs/plan_signo_costheta_taum.md`](docs/plan_signo_costheta_taum.md) and
`../ExamplesFCCFullSim/docs/signo_costheta_carga_tau.md`.

## 7. Formula → code map

| Formula | Implementation |
| --- | --- |
| $\mathcal{A}_\ell(\sin^2\theta_\text{eff})$ | `weightsPol._compute_ae_sm` |
| (3) $\mathcal{P}(z)_\tau$ | `weightsPol._compute_Ptau` (also inlined in `optimalVariabRho.wVariab`) |
| (4)(5)(6) $H_\pi$, $H_V$, $z_R$ | `weightsPol._compute_H`, `_alpha_V` |
| (13) $H_\ell$ | `weightsPol._compute_H_lep` |
| Exact boosted $\cos\theta^*$ (π, gen) | `weightsPol.cosThetaStar` |
| (15) optimal variable $\omega$ | `optimalVariabRho.wVariab` |
| Per-tau weight, hadronic | `weightsPol.newAtau` |
| Per-tau weight, leptonic | `weightsPol.newAtauLep` |
| Per-tau weight from precomputed $H$ | `weightsPol.newAtauFromH` (`newAtauRhoOmega` alias) |
| (16) joint weight, generic | `weightsPol.newAtauJoint` |
| (9) joint weight, had–had | `weightsPol.newAtauJoint_had_had` |
| (14) joint weight, had–lep | `weightsPol.newAtauJoint_had_lep` |

Masses are hardcoded in `weightsPol.py`: $m_\tau = 1.7769$, $m_\rho = 0.77545$,
$m_{a_1} = 1.2300$ GeV. For the $\rho$, $\alpha_V$ uses the **event-by-event**
invariant mass (consistent with `wVariab`); for the broad $a_1$ the pole mass is
used instead, being numerically more stable.

## 8. Corrections not included

The formulas above are the leading approximation. Known, deliberately omitted
effects, in rough order of size:

- **QED radiation (ISR/FSR).** All expressions assume a pure $e^+e^-\to\tau^+\tau^-$
  process. Near the Z pole these are mostly soft-collinear and expected at the
  permille level, and ISR/FSR interference can be neglected *at the peak* (not
  away from it). Mitigations if needed: for FSR, use generator information before
  radiation, or "dress" each charged particle with its nearby hard photons; for
  ISR, boost the $\tau^+\tau^-$ system to the Collins–Soper frame (two pure
  boosts: longitudinal $\vec\beta_z = (0,0,-p_z/E)$, then transverse
  $\vec\beta_T = (-p_x/\sqrt{p_T^2+M^2}, -p_y/\sqrt{p_T^2+M^2}, 0)$). Note that
  with large ISR the $\tau\tau$ invariant mass shifts and the relevant
  polarization is the one at the *updated* mass.
- **$\gamma^*$ exchange.** Shifts both the charge asymmetry and the polarization
  [17–26]. At the Z peak $\mathrm{Re}(\chi)=0$ and the correction reduces to
  eq. (27); for $z=0$ in the SM it is $\Delta\mathcal{P}(0)_\tau \approx -0.0002$,
  a 1.3 permille relative effect. It grows rapidly away from the pole.
- **Transverse spin correlations** between the decay products of the two taus.
  These exist even in the SM because $m_\tau \neq 0$, and are neglected by the
  factorized treatment. Kinematically they show up as aplanarity effects once all
  decay products are considered (Bernabéu, Rius, Pich 1991).

## 9. References

1. J. Alcaraz / FCC-CIEMAT team, *Differential cross section distributions for tau
   polarization at the Z pole*, 2 June 2026 — `../ExamplesFCCFullSim/Reweight_tautau.pdf`.
2. J. Bernabéu, N. Rius, A. Pich, *Tau spin correlations at the Z peak: Aplanarities
   of the decay products*, Phys. Lett. B **257** (1991) 219–226.
3. J. M. John, A. Tapadar, Z. Wąs, *On 'τ' Spin Use with KKMCee*, arXiv:2509.04400.
4. M. Davier, L. Duflot, F. Le Diberder, A. Rougé, *The Optimal method for the
   measurement of tau polarization*, Phys. Lett. B **306** (1993) 411–417.
5. J. C. Collins, D. E. Soper, *Angular distribution of dileptons in high-energy
   hadron collisions*, Phys. Rev. D **16** (1977) 2219.

---

# Tau Polarization Analysis (`RhoAnalysis/` — ρ / MDecs pipeline)

`RhoAnalysis/` contains the pipeline that **extracts the tau polarization
asymmetry `A_τ`** from the *optimal polarimeter observable* of each tau decay,
channel by channel. The theory behind it is summarized in the section above. It
is a three-stage, mostly-parallel workflow:

1. **Stage 1 — build a flat `TTree`** (one entry per event, two taus per entry)
   holding every kinematic quantity, the optimal observable `ω`/`optimalVar`, the
   gen helicity, and the **polarization reweighting weights** `weight_P1`/
   `weight_M1` (reweight the sample to `A_τ = ±1`). Two entry points:
   - `genOnlyRHOTree_MDecs_parallel.py` — **gen-level only** (fast; truth
     kinematics, reco branches mirror gen). Ideal for closure tests.
   - `analysisRHOTree_MDecs_parallel.py` — **full reco** (reads EDM4hep
     PandoraPFO / MLPF-GATr collections, runs tau reco + gen–reco matching,
     then fills the same branches plus their `reco_*` counterparts).
2. **Stage 2 — histograms** (`RhoHistFromTree_MDecs_parallel.py`): reads the tree
   and fills the `Omega_*` / observable distributions per decay category, with the
   nominal, `P1`, `M1` and correlated (`corr_P1`/`corr_M1`) reweighted variants,
   split into signal/background categories.
3. **Extraction** — `makeCosBins_MDecs.py` bins the observable in `cos θ` and
   `fitPolAssym.py` fits `A_τ` (and `A_e`) from the reweighted templates.

The polarization physics lives in the shared modules
`modules/weightsPol.py` (reweighting weights, Alcaraz joint formulas) and
`modules/optimalVariabRho.py` (optimal ρ observable `ω`). Auxiliary helpers moved
into `modules/` during the migration: `modules/rhoTreeUtils.py` (per-entry
variable extraction, `make_p4`, histogram filling) and
`modules/rhoParallelUtils.py` (file-splitting, per-worker logging, ROOT merge).

> **Decay-ID codes** used by `--decay-modes` / `--single-decay` / `--decay-pair`:
> `0` = π/K, `2` = ρ (π±π⁰), `10` = a₁, `-11` = e, `-13` = μ.
> Gen-only trees store ρ as `1` (the pipeline remaps reco `2 → 1`), so always
> pass `--only-gen` to Stage 2 when the tree came from `genOnlyRHOTree`.

---

## 1. Environment

The pipeline runs on the standard Key4hep stack (all commands from the repo root):

```bash
source setupKey4Hep.sh          # = source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-10-03
```

This release already provides ROOT, podio, edm4hep, numpy/pandas, `uproot` and
`joblib` (used by the cut optimizer) and `torch` (used by the MLP observable).
**MLP *training* additionally needs** `optuna` (hyper-parameter search), and
`shap` / `tabpfn` for the analysis variants — these are **not** in the base
stack; install them into a `--system-site-packages` virtualenv on top of Key4hep
(same pattern as `setupEventDisplay.sh`). Plain **inference** of the trained
observable (`modules/mlpPolInference.py`) is numpy-only and needs none of them.

---

## 2. Stage 1a — gen-level tree (`genOnlyRHOTree_MDecs_parallel.py`)

```bash
python RhoAnalysis/genOnlyRHOTree_MDecs_parallel.py \
    --sample ztt \
    --config config/default/taurecolong_optimal.yaml \
    --n-workers 8
# quick single-file test: add  --input-list /pnfs/.../out_reco_edm4hep_edm4hep_1.root
```

Output: `Results/RhoAnalysis/<prefix><cuts>/tau_trained*.root`, TTree
`outtree_original` with per-tau branches `tau{1,2}_{decayID,omega,optimalVar,
weight_P1,weight_M1,genHelicity,...}`.

| Option | Default | Description |
| --- | --- | --- |
| `-f`, `--sample` | from config | Sample name/alias in `config/samples/samples.yaml` (`ztt`, `p1`, `m1`, …). |
| `--input-list FILE [...]` | — | Explicit ROOT file(s); bypasses `samples.yaml`. |
| `-c`, `--config` | `config/default/taurecolong.yaml` | Analysis config (use `taurecolong_optimal.yaml`). |
| `--decay-modes ID [...]` | all | Keep only pairs whose taus are in this list. |
| `--sin-eff VAL` | `0.2312` | sin²θ_eff for the reweighting weights. |
| `--n-workers N` | min(n_files, n_cpus) | Parallel workers. |
| `--prefix STR` | `PolAnalysis_GEN_<SAMPLE>_NewWeights_` | Output-dir prefix. |
| `-v` / `-vv` | warnings | Verbosity (INFO / DEBUG). |

Cut flags `--tauCut --TauPhotonPCut --TauPionPCut --generalPCut --dRMax
--NeutronCut --MatchedGenMinDR` override the values in the config.

## 3. Stage 1b — full-reco tree (`analysisRHOTree_MDecs_parallel.py`)

Same output schema as Stage 1a plus the `reco_*` branches. Reads real detector
collections, so it is heavier — run it over a sample or a file list:

```bash
python RhoAnalysis/analysisRHOTree_MDecs_parallel.py \
    --sample ztt \
    --config config/default/taurecolong_optimal.yaml \
    --n-workers 8
```

Extra options on top of the common ones: `--gatr-result PATH` (use MLPF/GATr
predictions instead of PandoraPFOs), `-e/--electron-cut`, `-u/--muon-cut`,
`--lepton-xor-p`, `--sys-err config/systematics/err_sys.yml`.

## 4. Stage 2 — histograms (`RhoHistFromTree_MDecs_parallel.py`)

```bash
python RhoAnalysis/RhoHistFromTree_MDecs_parallel.py \
    --tree-file "Results/RhoAnalysis/<...>/tau_trained*.root" \
    --only-gen \                                   # required for gen-only trees
    --single-decay 2 \                             # or  --decay-pair 2 -13
    --hist-config-mdecs config/histograms/rho_analysis_config_mdecs.yml \
    --config config/default/taurecolong_optimal.yaml \
    --compute-weights \
    --n-workers 8
```

Output: `Results/RhoAnalysis/<...>/HistosMDecs_*.root` with the `Omega_*` family
(nominal + `_P1` / `_M1` / `_corr_P1` / `_corr_M1` reweighted variants) split by
signal/background category.

| Option | Default | Description |
| --- | --- | --- |
| `--tree-file PATH` | — | Input tree (**required**). |
| `--single-decay D` / `--decay-pair D0 D1` | — | Channel to histogram (**one required**). |
| `--only-gen` | off | Remap ρ `2→1`; **always set for `genOnlyRHOTree` trees**. |
| `--compute-weights` | off | Recompute reweighting from kinematics (needed with `--sin-eff`). |
| `--sin-eff VAL` | tree weights | Recompute weights at this sin²θ_eff. |
| `--no-omega-weights` | ω on | Use the simplified `H_V = α_V·z_R` weight for ρ instead of the full ω. |
| `--use-nn-optimal` + `--nn-model-path` | off | Use the MLP observable (see §8) as the ρ optimal variable. |
| `--meson-cut`/`--lepton-cut`/`--zmass-cut MIN MAX`, `--ang MIN MAX`, `--cut "EXPR"` | open | Event selection cuts. |
| `--n-workers N` | min(n_entries, n_cpus) | Parallel workers. |

## 5. One-shot orchestrator (`runTreeHistPipeline_MDecs.py`)

Runs Stage 1 → Stage 2 for a set of named runs described in a pipeline YAML
(examples in `config/pipeline/`, incl. `config/pipeline/PolPipeline/`). Use
`--dry-run` to print the exact commands without executing:

```bash
python RhoAnalysis/runTreeHistPipeline_MDecs.py \
    --pipeline-config config/pipeline/tree_hist_pipeline_example.yaml --dry-run

# only regenerate histograms from an existing tree:
python RhoAnalysis/runTreeHistPipeline_MDecs.py \
    --pipeline-config <cfg> --hist-only --tree-file <tree.root>
```

## 6. `cos θ` binning + `A_τ` fit

```bash
# 1) bin the observable in cos(theta) per decay pair
python RhoAnalysis/makeCosBins_MDecs.py \
    --sample-dir Results/RhoAnalysis/<...> \
    --target-decay rho --other-decay lep \
    -o ./Binned_histograms_MDecs/

# 2) fit A_tau (and A_e) from the reweighted templates
python RhoAnalysis/fitPolAssym.py \
    -i Binned_histograms_MDecs/BINED_*.root \
    --nBins 20 --bg-mode total -o ./Binned_histograms_MDecs/
```

`fitPolAssym.py` supports `--no-bg`, `--bg-mode {total,mig,ext,split}`, `--chi2`,
`--perfect` (truth closure), and luminosity scaling (`--lumi-base`,
`--lumi-target`, `--extra-legend`).

## 7. Selection-cut optimization (`optimize_cuts.py` / `evaluate_cuts.py`)

PSO-based optimization of the selection cuts (dR, meson-P, lepton-P, Z-mass)
against a signal-vs-background figure of merit, backed by `modules/optimize_pso/`:

```bash
python RhoAnalysis/optimize_cuts.py \
    --signal-root <signal_tree.root> --bg-root <bg1.root> <bg2.root> \
    --selectGEN 2 --eff-target 0.90 --particles 500 --iters 1000
```

`evaluate_cuts.py` takes the same inputs plus `--cut-list <cuts.csv>` to evaluate
a fixed list of working points. Score is `S/√(S+B)` by default, `S/(B+ε)` with
`--use-s-over-b`.

## 8. MLP optimal observable (`RhoAnalysis/MLP/`)

Trains and exports a neural-network polarimeter observable for the ρ channel:
`createPolDatasets.py` (build training datasets), `MLOptimalObservable.py` +
`OptimizeMLObservable.py` (train / hyper-optimize, needs `optuna`),
`exportMLPToNumpy.py` (export weights to the numpy-only inferencer
`modules/mlpPolInference.py`), `SHAPAnalysis.py` (feature importance, needs
`shap`). The exported model is consumed at Stage 2 via `--use-nn-optimal
--nn-model-path <weights.npz>`.

## 9. Quick start (copy-paste, gen-level closure)

```bash
source setupKey4Hep.sh
F=/pnfs/ciemat.es/data/cms/store/user/cepeda/FCC/FullSim/ZTauTau_SMPol_25Sept_MuonFix/out_reco_edm4hep_edm4hep_1.root

# Stage 1 (gen tree, single file)
python RhoAnalysis/genOnlyRHOTree_MDecs_parallel.py \
    --input-list "$F" -c config/default/taurecolong_optimal.yaml \
    -o quicktest --n-workers 1

# Stage 2 (rho histograms) — the tree stem follows the -o value ("quicktest")
TREE=$(find Results/RhoAnalysis/quicktest* -name '*.root' ! -name 'HistosMDecs_*' | head -1)
python RhoAnalysis/RhoHistFromTree_MDecs_parallel.py \
    --tree-file "$TREE" --only-gen --single-decay 2 --compute-weights \
    --hist-config-mdecs config/histograms/rho_analysis_config_mdecs.yml \
    --config config/default/taurecolong_optimal.yaml --n-workers 1
```

## 10. CLD vs ILD detector comparison (`RhoAnalysis/scripts_CLD_ILD/`)

A packaged end-to-end run of the whole pipeline over the two full-simulation samples
(`ztt_2M` = CLD, `ild_fcc` = ILD), for the π, ρ and a₁ channels each paired with a
leptonic hemisphere, using the **correlated** (two-hemisphere) weight templates.
The base configuration is `config/default/taurecolong.yaml` (no kinematic cuts) and
the only background is the internal migrations of the ττ sample.

```bash
bash RhoAnalysis/scripts_CLD_ILD/01_build_trees.sh        # reco + gen trees, 16 workers
bash RhoAnalysis/scripts_CLD_ILD/run_all_after_trees.sh   # stages 2 → 6, unattended
```

| Script | Stage |
| --- | --- |
| `01_build_trees.sh` | Reco (`analysisRHOTree`) and gen (`genOnlyRHOTree`) trees, decay modes `0 2 10 -11 -13`. |
| `02_histograms_nocuts.sh` | Six hadron+lepton channel pairs per tree, via `runTreeHistPipeline_MDecs.py` and the configs in `config/pipeline/CLD_ILD/`. |
| `03_optimize_cuts.sh` | PSO cut optimization per detector and channel. |
| `03b_make_optcut_configs.py` | Writes `<DET>_reco_optcuts.yaml` from the PSO results. |
| `04_binning_and_fit.py` | `makeCosBins_MDecs.py` (`--weights corr --bg-def ss`) + `fitPolAssym.py`, once at the MC-equivalent luminosity and once rescaled to `LUMI_FCC_YEAR_FB`. |
| `05_compare_plots.py` | Merges the e/μ files per channel and runs `CompareAlgs.py`: optimal variable for CLD/ILD/gen, kinematics, and the reweighted fit templates. |
| `06_summary.py` | Collects every fit and cut result into `Results/CLD_ILD_Summary/`. |

Outputs: figures in `TauPolOutputs/CLD_ILD/`, binned templates and fits in
`Binned_histograms_MDecs/CLD_ILD/`, cuts in `Results/CutOptimization_CLD_ILD/`.

Two helpers were added for this study and are reusable on their own:

- `RhoAnalysis/mdecsTreeToLegacy.py` projects a two-hemisphere MDecs tree onto the
  legacy one-hemisphere schema (`recoMesonP`, `lepP`, `genTauID`, …) that
  `optimize_cuts.py` and `modules/optimize_pso/` expect.
- `makeCosBins_MDecs.py` gained `--signal-ngen` / `--signal-lumi-pb`, needed whenever
  the processed sample is not the one hardcoded in `EVENT_CONFIG`; `optimize_cuts.py`
  now accepts any number of background files (including none), any `--selectGEN`, and
  an `--outdir`.

---

# Interactive Event Display (`EventDisplay/`)

A web-based event display built with **Dash + Plotly**: interactive 3D view of an
EDM4hep event (MC particles, PFOs, tracker/calorimeter hits), particle tables and
summary histograms, all in a single browser page.

Main features:

* **Detector profiles** for **CLD** and **ILD** (`--detector`), each with its own
  collection names and DD4hep compact-XML geometry (resolved automatically from
  `$k4geo_DIR`; selectable version with `--geometry-version`, e.g. `v06` for CLD,
  `v01`/`v02` for ILD). The ILD profile includes the TPC envelope.
* 3D geometry envelopes parsed directly from the compact XML.
* Hover info with PDG name, momentum and primary-particle ancestry.
* Filter to hide secondary particles.
* Gen–reco matching by ΔR or by truth links (`matching.py`).
* Lazy podio reading: only the requested event is decoded.

## Environment setup

The display needs a **Key4hep** stack (for `podio` ≥ 1.7, matching the reco files)
plus a small virtualenv with `dash`/`plotly` on top of it. Both are handled by
`setupEventDisplay.sh`, which must be **sourced**:

```bash
source setupEventDisplay.sh
```

On first use it creates the virtualenv (default `~/.venv/fcc-display-latest`,
override with `FCC_DISPLAY_VENV`) with `--system-site-packages` and installs
`dash`, `plotly` and `dash-bootstrap-components`. On later uses it simply loads
Key4hep (release `2026-04-08` by default, override with `KEY4HEP_RELEASE`) and
activates the virtualenv.

> The virtualenv is tied to the python version of the Key4hep release used to
> create it. If a newer release changes the python version, delete the venv
> directory and source the script again.

## Running

```bash
python EventDisplay/event_display_dash.py -i /path/to/reco_edm4hep.root --detector CLD
```

Options:

| Option               | Default   | Description                                        |
| -------------------- | --------- | -------------------------------------------------- |
| `-i, --input`        | —         | ROOT file to pre-load at startup (can also be loaded from the UI) |
| `--detector`         | `CLD`     | Detector concept: `CLD` or `ILD`                   |
| `--geometry-version` | per-detector | Compact-XML version tag (`v06` CLD, `v01`/`v02` ILD) |
| `--port`             | `8050`    | HTTP port                                          |
| `--host`             | `0.0.0.0` | Bind address                                       |
| `--debug`            | off       | Dash debug mode (auto-reload)                      |

Then open `http://<machine>:8050` in a browser. If running on a remote node
(e.g. `gaeui05`), either browse to `http://gaeui05:8050` from inside the network
or tunnel the port: `ssh -L 8050:localhost:8050 <node>` and open
`http://localhost:8050`.

---

## License

This project is intended for research and educational use in particle reconstruction studies.
Please cite appropriately if used in your work.

---

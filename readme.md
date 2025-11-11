# Tau Reconstruction at FCC-ee
This repository contains tools and scripts for reconstructing tau leptons at the Future Circular Collider (FCC-ee) using the EDM4hep data format. The focus is on validating tau reconstruction algorithms (such as PandoraPFO).

The main tools included are:
- Identification of the gennerator-level tau leptons and their decay products (decay type).
- Reconstruction of hadronic and leptonic tau decays from PandoraPFO collections.
- Matching reconstructed taus to generator-level taus.

# Simple Tau Reconstruction Test

This script (`test/simple_tau_reco_test.py`) provides a **basic example** for reconstructing **tau leptons** and related decay products from ROOT event data using the **EDM4hep** format.

It reads generator-level (`MCParticles`) and reconstructed-level (`PandoraPFOs`) collections, performs tau reconstruction, fills histograms for analysis, and writes the results to a ROOT output file.


The script performs the following tasks:

* Reconstruct **tau leptons** from **hadronic** and **leptonic** (**electrons** and **muons**) decays from PFO collections.
* Match reconstructed taus to their generator-level counterparts.
* Produce validation histograms (momentum, mass, decay type, and P resolution).
* Save the output in a ROOT file for further analysis.

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

| Histogram          | Description                             |
| ------------------ | --------------------------------------- |
| `histoGenTauPt`    | Generator-level tau transverse momentum |
| `histoGenTauVisPt` | Visible tau transverse momentum         |
| `histoGenTauP`     | Tau total momentum                      |
| `histoGenTauMass`  | True tau invariant mass                 |
| `histoRecoTauPt`   | Reconstructed tau transverse momentum   |
| `histoRecoTauType` | Tau decay mode (hadronic/leptonic)      |
| `histoResTauP`     | Momentum resolution ((Reco - Gen)/Gen)  |


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

## License

This project is intended for research and educational use in particle reconstruction studies.
Please cite appropriately if used in your work.

---

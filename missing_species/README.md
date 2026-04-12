# Missing Species Robustness Test

Tests whether ancestral karyotype reconstruction remains accurate when one or
more extant species are removed from the input.

## Scripts

| Script | Purpose |
|--------|---------|
| `evolution_simulator.py` | Simulate chromosome evolution (copied from root) |
| `ancestor_reconstruction.py` | Reconstruct ancestral karyotypes (configured for reduced data) |
| `drop_species.py` | Remove species from karyotype and tree files |
| `missing_species_experiment.py` | Batch experiment across multiple random scenarios |

## Quick Start

### Single run

**Step 1: Simulate**

```bash
python evolution_simulator.py
```

**Step 2: Drop species**

```bash
python drop_species.py --drop Sp3 Sp7
```

Writes reduced files to `output_simulator_reduced/`. Also handles `_preWGD`
variants automatically (e.g. dropping `Sp3` also removes `Sp3_preWGD` if present).

**Step 3: Reconstruct**

```bash
python ancestor_reconstruction.py
```

Output goes to `output_reconstruction_reduced/`.

## Batch Experiment

Runs N scenarios, each with a randomly chosen species count (13–23) and a
random number of dropped species (1 or 2). No full-species reconstruction is
performed — only the reduced dataset is reconstructed and compared against the
true root karyotype.

```bash
python missing_species_experiment.py           # 10 scenarios, default seed=42
python missing_species_experiment.py --n 20    # 20 scenarios
python missing_species_experiment.py --skip-existing  # resume interrupted run
```

Results are written to `output_experiments_missing/summary.tsv`.

### Summary columns

| Column | Description |
|--------|-------------|
| `Scenario` | Scenario ID |
| `Seed` | Random seed used |
| `SpeciesFull` | Total species simulated |
| `SpeciesReduced` | Species after dropping |
| `DroppedCount` | Number of species dropped (1 or 2) |
| `DroppedSpecies` | Names of dropped species |
| `TrueChrCount` | True ancestral chromosome count |
| `ReconChrCount` | Reconstructed chromosome count |
| `ExactMatch` | Chromosomes with exact gene-set match |
| `ExactMatchRate` | Exact match / true count |
| `GeneCoverage` | Fraction of true root genes recovered |

## Results (10 scenarios, seed=42)

| Species range | Drop range | Exact match | Gene coverage |
|---|---|---|---|
| 13–23 | 1–2 | **100%** (10/10) | **100%** (10/10) |

Reconstruction is robust to 1–2 missing species across all tested configurations.

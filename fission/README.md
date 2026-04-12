# Chromosome Fission Analysis

Tools for simulating chromosome fission events and detecting them in ancestral
karyotype reconstruction.

## Scripts

| Script | Purpose |
|--------|---------|
| `evolution_simulator.py` | Simulate chromosome evolution along a phylogenetic tree |
| `fission_simulator.py` | Introduce fission events into species karyotypes |
| `ancestor_reconstruction.py` | Reconstruct ancestral karyotypes from extant species |
| `genome_error_detector.py` | Detect fission events and misassembly errors |

## Quick Start

### Step 1: Run Evolution Simulator

```bash
python evolution_simulator.py
```

Generates:
- `output_simulator/karyotypes_species_with_outgroup.txt` — species karyotypes (pre-fission)
- `output_simulator/karyotypes_true_root.txt` — true root ancestor karyotype
- `output_simulator/tree.nwk` — phylogenetic tree (no outgroup)
- `output_simulator/tree_with_outgroup.nwk` — phylogenetic tree (with outgroup)

### Step 2: Run Fission Simulator

```bash
python fission_simulator.py \
    --input output_simulator/karyotypes_species_with_outgroup.txt \
    --output-karyotypes output_simulator/karyotypes_fission.txt \
    --report output_simulator/fission_report.tsv \
    --fission-rate 0.01 \
    --seed 42
```

Generates:
- `output_simulator/karyotypes_fission.txt` — karyotypes after fission events
- `output_simulator/fission_report.tsv` — ground-truth fission event list

> You can also manually split chromosomes in `karyotypes_species_with_outgroup.txt`
> to create `karyotypes_fission.txt`.

### Step 3: Reconstruct Ancestors

```bash
# Set in ancestor_reconstruction.py:
#   input_karyotypes = "output_simulator/karyotypes_fission.txt"
#   output_dir = "output_reconstruction_fission"
python ancestor_reconstruction.py
```

Generates:
- `output_reconstruction_fission/ancestor_gene_sets_by_node.tsv` — reconstructed
  ancestor karyotypes per node

### Step 4: Detect Fission Events

```bash
python genome_error_detector.py \
    --ancestors output_reconstruction_fission/ancestor_gene_sets_by_node.tsv \
    --species   output_simulator/karyotypes_fission.txt \
    --tree      output_simulator/tree.nwk \
    --output    error_detection_report.txt
```

| Argument | Default | Description |
|----------|---------|-------------|
| `--ancestors` | `output_reconstruction_fission/ancestor_gene_sets_by_node.tsv` | Reconstructed ancestor TSV |
| `--species` | `output_simulator/karyotypes_fission.txt` | Species karyotype file |
| `--tree` | `output_simulator/tree.nwk` | Newick tree (used for SIB_DOM misassembly detection) |
| `--output` | `error_detection_report.txt` | Output report path |

## Detection Logic

### Foundation: reconstruction as a reference

Fission detection is grounded in `ancestor_gene_sets_by_node.tsv` — the output
of `ancestor_reconstruction.py`. The reconstruction propagates chromosomes
top-down from the Root through every internal node, inferring each node's
karyotype via a hierarchy of evolutionary events:

| Event | What it produces |
|-------|-----------------|
| **Shared** | Chromosome conserved intact across two child lineages |
| **NCF** (nested chromosome fusion) | One chromosome nested inside another; residue fragments peeled off |
| **EEJ** (end-to-end joining) | Two chromosomes fused end-to-end |
| **RCT** (reciprocal chromosomal translocation) | Segments exchanged between two chromosomes |
| **Inversion** | Internal segment flipped within a chromosome |
| **WGD** (whole-genome duplication) | All chromosomes duplicated at a _preWGD node |

Every chromosome that is a natural product of any of these events — at any
node from Root down to the species' direct ancestor — will appear in the TSV
with its exact gene set. The detector uses this as a reference: a species
chromosome whose gene set is already explained by the reconstruction is a
legitimate evolutionary product.

### Fission Detection

For each root chromosome R and each extant species S:

1. Find all S chromosomes that overlap R's gene set.
2. Require 2+ **pure fragments** — chromosomes whose genes all belong to a
   single gene family (same alphabetic prefix).
3. Require the fragments to be **disjoint** — overlapping fragments indicate
   WGD duplication, not fission.
4. Require combined fragment coverage ≥ 50% of R.

Breakpoints are located at adjacent gene pairs in root gene order where the
assignment switches between fragments.


### Why a Fissioned Chromosome Remains Intact in the Reconstructed Root

Shared-Strict has higher priority than all other rules. Take `Sp9(18a/18b)` —
two fragments from splitting one ancestral chromosome — as an example:

```
--- Node Anc5 Inference ---   (left=Sp8, right=Sp9)
[Shared-Strict] Sp8(18) == Sp9(19) -> Sp8(18)
[Nested-Strict] Sp9(18a) in Sp8(19) -> Sp9(18a) + Residue Sp8(19):115-139
[Residue] Sp8(19) -> Sp8(19):115-139 size=25 telomeres=True
```

`Sp8(18)` is confirmed as the Anc5 ancestor. `Sp9(18a)` and `Sp9(18b)` are also considered potential ancestral chromosomes of Anc5.

```
--- Node Anc2 Inference ---   (left=Anc5, right=Anc6)
[Shared-Strict] Sp8(18) == Sp3(20) -> Sp8(18)
[Nested-Strict-LocalConflict] Sp9(18a) in Sp5(38) rejected  (overlap=114 with inferred Sp8(18))
[Nested-Strict-LocalConflict] Sp9(18b) in Sp5(38) rejected  (overlap=25  with inferred Sp8(18))
```

The intact chromosome propagates up via Shared-Strict at every node; fission
fragments are blocked by LocalConflict. The Root ends up with a single intact
chromosome — identical to the pre-fission ancestor — so the detector bypasses
the ancestor filter and looks for disjoint pure fragments directly in extant
species.

> **Limitation**: When the fissioned chromosome is itself the product of a prior evolutionary
fusion (EEJ or NCF), one of the resulting fragments may contain genes from
two ancestral families. Such a fragment fails the single-gene-family check
and cannot be detected as a pure fission piece by this method. Only fissions
where all resulting fragments are single-family are reliably detected.

## Output Files

| File | Description |
|------|-------------|
| `fission_report.tsv` | Ground-truth fission list from the simulator |
| `output_reconstruction_fission/ancestor_gene_sets_by_node.tsv` | Reconstructed ancestor karyotypes |
| `error_detection_report.txt` | Detection results from `genome_error_detector.py` |

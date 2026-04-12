# Assembly Error Simulation

Tools for simulating chromosome misassembly errors and evaluating their impact
on ancestral karyotype reconstruction.

## Scripts

| Script | Purpose |
|--------|---------|
| `evolution_simulator.py` | Simulate chromosome evolution along a phylogenetic tree |
| `misassembly_simulator.py` | Introduce misassembly errors into species karyotypes |
| `ancestor_reconstruction.py` | Reconstruct ancestral karyotypes from extant species |
| `genome_error_detector.py` | Detect misassembly errors and fission events |

## Quick Start

### Step 1: Run Evolution Simulator

```bash
python evolution_simulator.py
```

Generates:
- `output_simulator/karyotypes_species_with_outgroup.txt` — species karyotypes
- `output_simulator/karyotypes_true_root.txt` — true root ancestor karyotype
- `output_simulator/tree.nwk` — phylogenetic tree (no outgroup)
- `output_simulator/tree_with_outgroup.nwk` — phylogenetic tree (with outgroup)

### Step 2: Simulate Misassembly Errors

```bash
python misassembly_simulator.py \
    --input output_simulator/karyotypes_species_with_outgroup.txt \
    --output-karyotypes output_simulator/karyotypes_misassembled.txt \
    --report output_simulator/misassembly_report.tsv \
    --error-rate 0.05 \
    --seed 42
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--error-rate` | 0.01 | Fraction of chromosomes to misassemble |
| `--max-errors` | None | Optional cap on total misassemblies |
| `--max-fragment-ratio` | 0.5 | Max fragment length as fraction of source chromosome |
| `--seed` | 42 | Random seed |

### Step 3: Reconstruct Ancestors (clean)

```bash
# Set in ancestor_reconstruction.py:
#   input_karyotypes = "output_simulator/karyotypes_species_with_outgroup.txt"
#   output_dir = "output_reconstruction"
python ancestor_reconstruction.py
```

### Step 4: Reconstruct Ancestors (misassembled)

```bash
# Set in ancestor_reconstruction.py:
#   input_karyotypes = "output_simulator/karyotypes_misassembled.txt"
#   output_dir = "output_reconstruction_misassembled"
python ancestor_reconstruction.py
```

### Step 5: Detect Errors

```bash
python genome_error_detector.py \
    --ancestors output_reconstruction_misassembled/ancestor_gene_sets_by_node.tsv \
    --species   output_simulator/karyotypes_misassembled.txt \
    --tree      output_simulator/tree.nwk \
    --output    error_detection_report.txt
```

| Argument | Default | Description |
|----------|---------|-------------|
| `--ancestors` | `output_reconstruction_misassembled/ancestor_gene_sets_by_node.tsv` | Reconstructed ancestor TSV |
| `--species` | `output_simulator/karyotypes_misassembled.txt` | Species karyotype file |
| `--tree` | `output_simulator/tree.nwk` | Newick tree (required for SIB_DOM) |
| `--output` | `error_detection_report.txt` | Output report path |
| `--min-dup-contiguity` | 0.7 | Min contiguity score for DUP detection |

## Detection Methods

### Foundation: reconstruction as a filter

The detector's primary signal comes from `ancestor_gene_sets_by_node.tsv` —
the output of `ancestor_reconstruction.py`. The reconstruction algorithm
propagates chromosomes top-down from the Root through every internal node,
inferring each node's karyotype via a hierarchy of evolutionary events:

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
with its exact gene set. The detector uses this as a whitelist: if a species
chromosome's gene set is already explained by the reconstruction, it is a
legitimate evolutionary product and is not flagged.

Only chromosomes whose gene sets are **absent from the entire reconstruction
hierarchy** are candidates for misassembly detection.

---

`genome_error_detector.py` then applies three methods to the remaining
candidates, plus a fission detector.

### Misassembly: DUP

Detects a duplicated gene block — a contiguous run of second-occurrence copies
of genes already present earlier in the chromosome. Requires contiguity ≥ 0.7
(at least 70% of the span is duplicate genes).

### Misassembly: SANDWICH

Detects a foreign gene block flanked by the chromosome's dominant gene family.
The algorithm scans for every maximal contiguous block that contains no genes
from the main root chromosome, then applies three filters:

| Filter | What it removes |
|--------|----------------|
| WGD check | Main-family count exceeds root chromosome size → WGD product |
| Ancestor check | Full chromosome gene set matches any reconstructed ancestor → natural rearrangement |
| Junction novelty | Either flanking junction pair is present in ancestor adjacencies → natural rearrangement junction |

The junction novelty check is the key discriminator: misassembly breakpoints
create novel gene adjacencies that never appeared in any ancestor chromosome,
while natural rearrangement junctions (NCF, EEJ, RCT) are preserved in the
reconstruction hierarchy. Both junction pairs must be absent from all ancestor
adjacencies for a block to be flagged.

### Misassembly: SIB_DOM

Catches misassemblies where the reconstruction accepted the error and promoted
it to a high-level ancestor node (so the ancestor filter passes). Compares the
count of each root chromosome's genes in the candidate against the maximum
found in any non-WGD sibling branch. A ratio > 1.0 with excess ≥ 5 genes
flags the chromosome as suspicious — the inflated count cannot be explained by
normal inheritance through the reconstruction hierarchy.

### Fission

Detects a root chromosome split into 2+ disjoint pure fragments (single-gene-
family chromosomes) in an extant species. Overlapping fragments are excluded
(they indicate WGD duplication, not fission). Requires ≥ 50% coverage of the
root chromosome across all fragments.

## Output Files

| File | Description |
|------|-------------|
| `misassembly_report.tsv` | Ground-truth misassembly list from the simulator |
| `output_reconstruction/` | Reconstruction from clean karyotypes |
| `output_reconstruction_misassembled/` | Reconstruction from misassembled karyotypes |
| `error_detection_report.txt` | Detection results from `genome_error_detector.py` |

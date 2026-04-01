# Chromosome Evolution Simulation and Ancestral Karyotype Reconstruction

Tools for simulating chromosome evolution and reconstructing ancestral karyotypes.

## Quick Start

### Install Dependencies

```bash
pip install ete3 matplotlib
```

### Basic Workflow

```bash
# 1. Run simulation
python evolution_simulator.py

# 2. Reconstruct ancestors
python ancestor_reconstruction.py

# 3. Run experiments (optional)
python experiment_table.py --param-mode random --num-scenarios 10
```

## Scripts

| Script | Purpose |
|--------|---------|
| `evolution_simulator.py` | Simulate chromosome evolution along phylogenetic tree |
| `ancestor_reconstruction.py` | Reconstruct ancestral karyotypes ([detailed guide](ANCESTOR_RECONSTRUCTION.md)) |
| `experiment_table.py` | Run batch experiments for accuracy evaluation |

## Output Files

### Simulation (`output_simulator/`)

| File | Description |
|------|-------------|
| `tree.nwk` | Phylogenetic tree |
| `karyotypes_species_with_outgroup.txt` | Extant species karyotypes |
| `karyotypes_true_root.txt` | True ancestral karyotype |
| `events.txt` | Branch-level event log |

### Reconstruction (`output_reconstruction/`)

| File | Description |
|------|-------------|
| `ancestors_by_node.tsv` | Per-node ancestral chromosomes |
| `ancestor_gene_sets_by_node.tsv` | Per-node ancestral gene sets |
| `root_rct_outgroup_decisions.tsv` | RCT ancestor decisions |
| `validation_report.txt` | Lineage isolation validation |
| `inference_log.txt` | Full details of ancestral karyotype reconstruction process |

## Model Assumptions

- **No gene loss**: Pipeline assumes no gene deletions/insertions during simulation and reconstruction. Coverage check warns if genes are missing.
- **Unique homology mapping**: Gene IDs are globally unique and represent one-to-one orthology across species.
- **Restricted event types**: Core events are RCT, EEJ, NCF and inversions. Fission is handled separately in `fission/` module.
- **Telomere flag**: `telomeres=True/False` distinguishes complete chromosomes from fragments, affecting strict matching and Root validation.
- **Root child event constraint**: If a chromosome pair (X, Y) undergoes RCT at one root child branch, the same pair cannot undergo NCF or EEJ at the other root child branch. This prevents ambiguous scenarios where RCT + NCF/EEJ combinations on the same chromosome pair would complicate ancestral state inference.

## Reconstruction Algorithm

Multi-stage approach:

| Stage | Method | Implementation |
|-------|--------|----------------|
| 1. Strict Inference | Shared/Nested pattern detection (bottom-up) | `_compare_and_merge()` |
| 2. Residual Discovery | Iterative cross-branch residual matching | `_discover_from_residuals()` |
| 3. RCT Decision | Outgroup-based ancestral state determination | `root_rct_outgroup_decision()` |
| 4. Final Promotion | Promote pending ancestors with sufficient support | `_finalize_root_pending_ancestors()` |

### 1. Strict Inference (postorder, bottom-up)

Two strict patterns are inferred in `_compare_and_merge()`:

- **Shared (Strict)**: Two chromosomes have identical gene sets (100% equality) and both are telomere-complete.
- **Nested (Strict)**: One chromosome's gene set is a true subset of the other, forming a contiguous block in the larger chromosome (inversions allowed; gaps not allowed).

During Nested inference, remaining genes are peeled as residue fragments and marked as `telomeres=False`.

### 2. Residual Discovery (iterative)

Key workflow in `_discover_from_residuals()`:

1. For each internal node, collect residual pools from all descendant leaves
2. Subtract confirmed ancestral gene sets from residuals
3. Compare leftover fragments between sibling branches
4. If leftovers are identical (100% set equality), promote as new ancestral chromosome with `provenance="Residual (Shared)"` and `telomeres=True`

### 3. Root RCT Outgroup Decision

When outgroup is available, the algorithm determines ancestral RCT states by:

1. Scanning candidate chromosome pairs from left/right child subtrees
2. Computing fusion point signatures (original vs translocation pairs)
3. Comparing against outgroup adjacency patterns
4. Selecting configuration with higher outgroup support

### Telomere Flag (`telomeres=True/False`)

- **`telomeres=True`**: Complete chromosome candidates
  - Extant input chromosomes default to True
  - Strict Shared/Nested inferred chromosomes are True
  - Residual-Shared promoted chromosomes are True

- **`telomeres=False`**: Fragment residues
  - Peeled residues during Nested inference
  - Not treated as complete ancestral chromosomes

Root validated chromosomes require both: allowed `provenance` AND `telomeres=True`.

## Key Concepts

- **RCT (Reciprocal Chromosome Translocation)**: Two chromosomes exchange segments
- **EEJ (End-to-End Joining)**: Two chromosomes fuse end-to-end
- **NCF (Nested Chromosome Fragment)**: Fragment inserted into chromosome
- **Fusion Point**: Adjacent gene pair indicating chromosome fusion event
- **Lineage Isolation**: Fusion points should appear in single lineages only

## Configuration

### Simulator

```python
CONFIG = dict(
    num_modern_species=16,
    num_ancestor_chromosomes=16,
    min_genes_per_chr=100,
    max_genes_per_chr=1000,
    rearrangement_counts=dict(
        inversion_prob=0.8,
        translocation_rct=0.4,
        fusion_ncf=0.2,
        fusion_eej=0.3,
    ),
    wgd_probability=0.05,
)
```

**Note**: The phylogenetic tree is built as a strict binary tree using ete3's `populate()` method. 

### Reconstructor

```python
CONFIG = dict(
    input_tree="output_simulator/tree.nwk",
    input_karyotypes="output_simulator/karyotypes_species_with_outgroup.txt",
    input_true_root_karyotype="output_simulator/karyotypes_true_root.txt",
    output_dir="output_reconstruction",
    min_block_size=2,
    enable_cross_branch_shared=True,
    enable_cross_branch_nested=True,
    enable_conflict_check=True,
    enable_root_rct_outgroup_decision=True,
    outgroup_name="Outgroup",
)
```

## Subfolders

| Folder | Description |
|--------|-------------|
| `fission/` | Chromosome fission simulation and conflict analysis |
| `assembly_errors/` | Misassembly error simulation and impact evaluation |

See [fission/README.md](fission/README.md) and [assembly_errors/README.md](assembly_errors/README.md) for module-specific workflows.

## File Structure

```
karyotype-phylogenomics-simulator/
├── evolution_simulator.py
├── ancestor_reconstruction.py
├── experiment_table.py
├── README.md
├── fission/                    # Fission analysis module
│   ├── README.md
│   ├── fission_simulator.py
│   ├── fission_point_conflict_analyzer.py
│   └── ...
├── assembly_errors/            # Misassembly simulation module
│   ├── README.md
│   ├── misassembly_simulator.py
│   └── ...
├── output_simulator/           # Simulation outputs
└── output_reconstruction/      # Reconstruction outputs
```

## License

Research use only.
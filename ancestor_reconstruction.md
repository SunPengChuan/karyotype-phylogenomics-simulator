# Ancestral Karyotype Reconstruction — Complete Guide

## Overview

`ancestor_reconstruction.py` is the core algorithm for reconstructing ancestral karyotypes at all nodes of a phylogenetic tree. It uses a multi-stage bottom-up and top-down approach combining strict pattern matching, residual discovery, outgroup validation, and post-processing logic.

**Key Features**:
- Bottom-up strict inference (Shared/Nested patterns)
- Iterative residual discovery across branches
- Outgroup-based RCT decision making
- Six post-processing logic passes for edge cases
- Top-down intermediate node completion
- Comprehensive validation and logging

---

### Basic Usage

```bash
# Run reconstruction (uses hardcoded CONFIG paths)
python ancestor_reconstruction.py

# Outputs will be in output_reconstruction/
```


## Input Files

All paths are configured in `CONFIG` dict (lines 14-26):

| File | Path | Description |
|------|------|-------------|
| Tree | `output_simulator/tree.nwk` | Phylogenetic tree (Newick format) |
| Karyotypes | `output_simulator/karyotypes_species_with_outgroup.txt` | Extant species + outgroup karyotypes |
| True Root | `output_simulator/karyotypes_true_root.txt` | Ground truth for validation |

### Karyotype File Format

```
Species  Chr_ID  [Length]  Gene1 Gene2 Gene3 ...
Sp1      Chr1    500       A1 A2 A3 B1 B2 ...
Sp1      Chr2    300       C1 C2 C3 ...
Outgroup Chr1    450       A1 A2 A3 B1 B2 ...
```

- Gene names: `<Chr_prefix><number>` (e.g., `A123`, `E45`)
- Chr prefix (A, B, C, ...) indicates ancestral chromosome family
- Genes with `-` prefix indicate reversed orientation (handled as unsigned)

---

## Output Files

All outputs in `output_reconstruction/`:

| File | Description |
|------|-------------|
| `ancestors_by_node.tsv` | Chromosome IDs per node |
| `ancestor_gene_sets_by_node.tsv` | Full gene sets per node (main output) |
| `validation_report.txt` | Lineage isolation validation for fusion points |
| `inference_log.txt` | Complete reconstruction trace (~4000+ lines) |
| `tree_with_ancestors.png` | Phylogenetic tree visualization with ancestor counts |

---

## Algorithm Overview

### Stage 1: Bottom-Up Strict Inference

**Method**: `_compare_and_merge()` (postorder traversal)

Two strict patterns detected between sibling nodes:

#### Shared (Strict)
- **Condition**: Two chromosomes have identical gene sets (100% equality)
- **Requirement**: Both must have `telomeres=True`
- **Result**: Promoted as ancestral chromosome with `provenance="Shared (Strict)"`

#### Nested (Strict)
- **Condition**: One chromosome's genes form a contiguous block in another
- **Requirement**: Container chromosome must have `telomeres=True`
- **Result**: Subset promoted as ancestor; remaining genes peeled as residue with `telomeres=False`

**Key Implementation Details**:
- Only telomeric chromosomes (`telomeres=True`) can be Nested-Strict candidates
- After Nested match, residue fragments inherit telomere status from source via `get_telo_k1/k2(id)`
- Range-based IDs (e.g., `Sp14(3):375-891`) created by `peel_residue()` may have `telomeres=True` if they cover a chromosomal end

---

### Stage 2: Residual Discovery

**Method**: `_discover_from_residuals()` (iterative, multi-pass)

**Algorithm**:
1. Collect residual pools from all descendant leaves for each internal node
2. Subtract already-confirmed ancestral gene sets
3. Compare leftover fragments between sibling branches
4. If leftovers are identical (100% set equality), promote as `provenance="Residual (Shared)"` with `telomeres=True`

**Purpose**: Recovers ancestral chromosomes that were modified in ALL descendants (no intact copy preserved), so bottom-up Shared/Nested couldn't detect them.

---

### Stage 3: Root Post-Processing

After bottom-up inference completes, only chromosomes with Shared-Strict or Nested-Strict provenance are in `confirmed_ancestors[Root]`. Many valid root chromosomes remain unconfirmed because:
- They were modified in all descendants (no intact copy for Shared matching)
- They underwent complex rearrangements (RCT, fusion chains)
- They exist only as residual fragments in `node_karyotypes`

Root post-processing applies specialized logic to recover these missing chromosomes:

#### RCT-Residual-Validated (~line 1480)
Merges residual fragments at root using outgroup-based RCT validation.
- **Input**: Residual fragments in `node_karyotypes[Root]`
- **Method**: Find pairs of residuals that form valid RCT configurations based on outgroup adjacency
- **Cross-Family Guard**: Rejects merges spanning multiple gene families
- **Output**: `provenance="RCT-Residual-Validated"`

#### Root-Anchored-Residual (~line 1570)
Promotes residuals when their anchor chromosome is split across root's children.
- **Input**: Residual fragments in `node_karyotypes[Root]`
- **Method**: If a residual's genes appear split across both root children, promote it
- **Output**: `provenance="Root-Anchored-Residual"`

**Note**: The original "Logic 1-6" framework described in earlier documentation has been simplified. Current implementation focuses on RCT-based validation and residual promotion.

---

### Stage 4: Same-Family Merge (2 passes)

**Problem**: Logic 5 and Logic 6 can independently identify disjoint subsets of the same gene family, creating multiple root chromosomes for one family.

**Solution** (~line 2161 in `run()`):
- Group validated root chromosomes by gene family prefix
- Merge all same-family disjoint groups into one chromosome
- **Provenance**: `Same-Family-Merge`

**Two passes**:
1. After initial Logic 1-6
2. After Missing-Family-Rescue (to merge rescued fragments)

---

### Stage 5: Missing-Family Rescue

**Problem**: A gene family may be completely absent from root due to complex fusion patterns across all branches.

**Solution** (~line 2209):
- **Detection**: Compare all input genes (union across all species) against root gene coverage
- **Search scope**:
  1. First search `confirmed_ancestors` at intermediate nodes
  2. If not found, search `node_karyotypes` for telomeric chromosomes
- **Validation**: Candidate must contain complete missing family genes with no other-family contamination
- **Provenance**: `Missing-Family-Rescue`

**Key improvement** (2024 fix): Extended search to `node_karyotypes` with telomere filter, enabling recovery of chromosomes involved in complex RCT+EEJ fusion chains (e.g., scenario_049: G and Y families).

---

### Stage 6: Top-Down Intermediate Node Completion

**Method**: `_topdown_complete_intermediate_nodes()` (~line 2361)

**Purpose**: Fill intermediate nodes with complete karyotypes using validated root as anchor.

**Algorithm**:
1. Start from validated root karyotype
2. For each internal node (preorder, skip root/leaves):
   - Keep bottom-up confirmed (Shared/Nested) as highest-trust
   - For each parent chromosome not yet covered:
     - **Case 1 (Intact)**: Exact match in `node_karyotypes[N]` with `telomeres=True` → add
     - **Case 2 (Split)**: Distributed across multiple telomeric pieces → add all pieces
     - **Case 3 (Fusion)**: Already in confirmed via Shared/Nested → skip
     - **Case 4 (Mixed)**: Complex rearrangement → skip

**Safety**: Only adds chromosomes already in `node_karyotypes[N]`; overlap check prevents gene duplication.

---

## Key Concepts

### Telomeres Flag

**`telomeres=True`**: Complete chromosome or end-spanning fragment
- Extant input chromosomes (default)
- Shared/Nested confirmed chromosomes
- Residual-Shared promoted chromosomes
- Range-based residues covering chromosomal ends

**`telomeres=False`**: Internal fragment
- Peeled residues from Nested inference (middle sections)
- Cannot participate as Nested-Strict candidates
- Not validated as root ancestors

### Gene Families

Genes named `<prefix><number>` (e.g., `E292`, `C77`):
- Alphabetic prefix = ancestral chromosome family
- By simulator design: each family has EXACTLY ONE chromosome at root
- Used in Logic 5 cross-family guard to reject spurious merges

### Provenance Types

| Provenance | Source | Telomeres |
|------------|--------|-----------|
| `Shared (Strict)` | Bottom-up identical match | True |
| `Nested (Strict)` | Bottom-up subset match | True |
| `Residual (Shared)` | Cross-branch residual match | True |
| `Logic1-RCT-Outgroup-Validated` | Outgroup validation | True |
| `Logic2-Cross-Node-Confirmed` | Multi-node confirmation | True |
| `Logic5-RCT-Residual-Validated` | Residual merge | True |
| `Same-Family-Merge` | Post-processing merge | True |
| `Missing-Family-Rescue` | Family rescue | True |
| `Unmatched` | No pattern found | Varies |

### Data Structures

**`confirmed_ancestors[node]`**: Only Shared/Nested provenances (highest trust)

**`node_karyotypes[node]`**: Complete inferred karyotype including Residue/Unmatched

**`topdown_karyotypes[node]`**: Top-down completed karyotypes for export

**`validated_root_chrs`**: Final root chromosomes after all post-processing

---

## Configuration Reference

```python
CONFIG = dict(
    input_tree="output_simulator/tree.nwk",
    input_karyotypes="output_simulator/karyotypes_species_with_outgroup.txt",
    input_true_root_karyotype="output_simulator/karyotypes_true_root.txt",
    output_dir="output_reconstruction",
    tree_viz_output_dir="output_simulator",
    enable_tree_viz=True,           # Generate tree PNG with ancestor counts
    enable_root_dotplot=True,       # Generate root dotplot
    min_block_size=2,               # Minimum gene block size for Nested matching
    enable_residual_rct_discovery=False,
    outgroup_name="Outgroup",       # Name of outgroup in karyotype file
)
```

---

## Class Architecture

| Class | Role |
|-------|------|
| `DataLoader` | Parse karyotype file and tree |
| `Karyotype` | Hold chromosome→genes mapping with attributes |
| `AncestorReconstructor` | Main reconstruction engine |

### Key Methods (AncestorReconstructor)

| Method | Line | Purpose |
|--------|------|---------|
| `_compare_and_merge()` | ~680 | Shared/Nested bottom-up inference |
| `_discover_from_residuals()` | ~900 | Cross-branch residual matching |
| `root_rct_outgroup_decision()` | ~1100 | RCT outgroup validation |
| `_finalize_root_pending_ancestors()` | ~1300 | Logic 1-6 post-processing |
| `_topdown_complete_intermediate_nodes()` | ~2361 | Top-down node completion |
| `_export_node_ancestor_sets()` | ~2480 | Write TSV outputs |
| `run()` | ~2100 | Main orchestration |


## Running in Batch Mode

```bash
# Run 30 scenarios and collect accuracy results
python experiment_table.py --num-scenarios 30

# Skip simulation, only collect existing results
python experiment_table.py --skip-run --num-scenarios 30

# Disable collinearity visualization (faster)
python experiment_table.py --no-visualize --num-scenarios 30
```

Batch results written to `output_experiments/results.tsv`.

---

## Verifying After Changes

```bash
# Scenario 006 (expect 9)
cp ancestor_reconstruction.py output_experiments/scenario_006/
cd output_experiments/scenario_006 && python ancestor_reconstruction.py


# Scenario 010 (expect 23)
cp ancestor_reconstruction.py output_experiments/scenario_010/
cd output_experiments/scenario_010 && python ancestor_reconstruction.py

# Scenario 016 (expect 5)
cp ancestor_reconstruction.py output_experiments/scenario_016/
cd output_experiments/scenario_016 && python ancestor_reconstruction.py

# Full batch regression (expect 30/30)
python experiment_table.py --no-visualize --skip-existing --num-scenarios 30

```

---

## Reading the Output

### `ancestor_gene_sets_by_node.tsv`

Each row is one ancestral chromosome at a node:

```
node    chr_id              genes
Root    RootChr1            A1,A2,A3,A4,...
Root    RootChr2            B1,B2,B3,...
Anc1    Sp1(1)              A1,A2,A3,...
Anc1    Sp1(1):200-500      B5,B6,B7,...
```

- **Root rows**: validated root karyotype (all Logic passes applied)
- **Intermediate node rows**: top-down completed karyotypes
- Gene order within a chromosome is preserved from source

### `inference_log.txt`

Read this to trace why a specific chromosome was or wasn't promoted:

```
[Shared-Strict] Sp1(1) == Sp3(1) → Root confirmed A1..A450
[Nested-Strict] Sp1(2) ⊂ Sp3(3) → Root confirms B1..B200, peel residue B201..B400
[Logic6] Anchored residue B201..B400 promoted
[Same-Family-Merge] B-family: merge {B1..B200} + {B201..B400} → single chr
```

---

See [README.md](README.md) for project overview and other scripts.

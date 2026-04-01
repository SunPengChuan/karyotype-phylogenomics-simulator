#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations
import os
from collections import Counter, defaultdict
import random
import string
import copy
import statistics
import matplotlib
matplotlib.use('Agg')  # Use non-GUI backend to avoid Qt conflicts
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
# Import Tree last to avoid Qt initialization issues
from ete3 import Tree
from ete3.treeview import TreeStyle, NodeStyle, TextFace
import itertools

# Type aliases
Gene = tuple[str, str]  # (gene_name, orientation)
Karyotype = dict[str, list[Gene]]  # chromosome_id -> list of genes
ConfigDict = dict

# -----------------------------
# 1. CONFIGURATION
# -----------------------------
CONFIG = dict(
    # Tree settings
    num_modern_species=19,
    # num_ancestor_nodes is auto-calculated for strict binary tree: species - 2
    seed=68,
    
    # Karyotype settings
    num_ancestor_chromosomes=21,  # A-P
    min_genes_per_chr=100,
    max_genes_per_chr=1000,
    
    # Evolution settings (Counts per branch)
    # Reduced frequencies for simplified evolution
    rearrangement_counts=dict(
        inversion_prob=0.922,   # Probability of inversion per branch (0-1)
        translocation_rct=0.72,  # Probability of RCT per branch (0-1)
        fusion_ncf=0.489,       # Reduced to 0.2 (20% chance per branch)
        fusion_eej=0.268,       # Reduced to 0.3 (30% chance per branch)
    ),
    rct_pairing_mode="random",
    eej_end_mode="random",
    
    # Whole Genome Duplication (WGD)
    wgd_probability=0.027,     # Reduced from 0.1
    force_wgd_nodes=[],       # List of node names to force WGD (e.g., ["Anc5"])
    enable_pre_wgd_nodes=True,
    wgd_pre_event_fraction=0.5,
    wgd_rng_seed_offset=1000,
    
    # Collinearity Plot Settings
    visualize_collinearity=False,
    collinearity_dir="species_vs_ancestor_plots",
    
    # Validation
    validate_lineage=False,
    validation_report="validation_report.txt",
    validation_strict=True,

    avoid_duplicate_root_fusions=True,

    enable_outgroup_for_rct=True,
    outgroup_name="Outgroup",
    outgroup_seed_offset=1000003,
    outgroup_min_block_len=20,
    outgroup_break_rate=0,
    outgroup_allow_block_reverse=False,

    rct_significance=dict(
        true_root_path="output_simulator/karyotypes_true_root.txt",
        karyotypes_path="output_simulator/karyotypes_species_with_outgroup.txt",
        out_path="output_simulator/rct_significance.tsv",
        outgroup_name="Outgroup",
        seed=68,
        replicates=200,
        min_block_len=20,
        break_rate=0.02,
        allow_block_reverse=False,
        min_chr_len=20,
    ),
    
    # Output
    save_dir="output_simulator",
    output_files=dict(
        species_karyotypes="karyotypes_species_with_outgroup.txt",
        true_root_karyotype="karyotypes_true_root.txt",
        events="events.txt",
        tree_img="tree_with_wgd.png",
        tree_nwk="tree.nwk",
        tree_nwk_with_outgroup="tree_with_outgroup.nwk"
    )
)

# -----------------------------
# 2. HELPER FUNCTIONS
# -----------------------------
def get_chr_labels(n: int):
    base = list(string.ascii_uppercase)
    if n <= 26: return base[:n]
    return base + [a + b for a in base for b in base][: n - 26]

def format_gene(gene_tuple: Gene) -> str:
    return f"-{gene_tuple[0]}" if gene_tuple[1] == "-" else gene_tuple[0]

def format_gene_id(gene_tuple: Gene) -> str:
    return gene_tuple[0]

def dedup_by_gene_id(genes: list[Gene]) -> list[Gene]:
    seen: set[str] = set()
    out: list[Gene] = []
    for g in genes:
        gid = g[0]
        if gid in seen:
            continue
        seen.add(gid)
        out.append(g)
    return out


class _RCTOutgroupSignificance:
    @staticmethod
    def _read_karyotypes_species(path: str) -> dict[str, Karyotype]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"karyotypes file not found: {path}")
        with open(path, "r", encoding="utf-8") as f:
            lines = [line.rstrip("\n") for line in f]
        start_idx = 0
        for i, line in enumerate(lines):
            if line.startswith("-"):
                start_idx = i + 1
                break
        if start_idx == 0 and lines and lines[0].startswith("Species"):
            start_idx = 1

        out: dict[str, Karyotype] = {}
        for raw in lines[start_idx:]:
            parts = raw.strip().split()
            if len(parts) < 3:
                continue
            sp, cid, *genes = parts
            out.setdefault(sp, {})[cid] = genes  # type: ignore
        return out

    @staticmethod
    def _read_true_root(path: str) -> Karyotype:
        data = _RCTOutgroupSignificance._read_karyotypes_species(path)
        if "TrueRoot" not in data:
            raise ValueError("TrueRoot not found in true root karyotypes file.")
        return data["TrueRoot"]

    @staticmethod
    def _segment_into_blocks(genes: list[str], rng: random.Random, min_block_len: int, break_rate: float) -> list[list[str]]:
        n = len(genes)
        if n <= min_block_len:
            return [list(genes)]
        target_blocks = max(1, int(n * break_rate))
        target_blocks = min(target_blocks, max(1, n // min_block_len))
        blocks: list[list[str]] = []
        pos = 0
        remaining_blocks = target_blocks
        while pos < n:
            remaining = n - pos
            if remaining_blocks <= 1:
                blocks.append(list(genes[pos:]))
                break
            max_len = remaining - (remaining_blocks - 1) * min_block_len
            if max_len < min_block_len:
                max_len = min_block_len
            if max_len > remaining:
                max_len = remaining
            blen = rng.randint(min_block_len, max_len)
            blocks.append(list(genes[pos : pos + blen]))
            pos += blen
            remaining_blocks -= 1
        return blocks

    @staticmethod
    def build_outgroup_by_block_shuffle(true_root_karyo: Karyotype, seed: int, min_block_len: int, break_rate: float, allow_block_reverse: bool) -> Karyotype:
        rng = random.Random(seed)
        template: list[tuple[str, int]] = []
        blocks: list[list[str]] = []
        for cid in sorted(true_root_karyo.keys()):
            genes = [g[0] for g in true_root_karyo[cid]]  # extract gene names
            template.append((cid, len(genes)))
            blocks.extend(_RCTOutgroupSignificance._segment_into_blocks(genes, rng, min_block_len, break_rate))
        rng.shuffle(blocks)

        out: Karyotype = {}
        for cid, target_len in template:
            assembled: list[Gene] = []
            while len(assembled) < target_len:
                if not blocks:
                    raise RuntimeError("Outgroup build failed: insufficient blocks.")
                b = blocks.pop()
                if allow_block_reverse and rng.random() < 0.5:
                    b = list(reversed(b))
                need = target_len - len(assembled)
                if len(b) <= need:
                    assembled.extend((g, "+") for g in b)
                else:
                    head = b[:need]
                    tail = b[need:]
                    assembled.extend((g, "+") for g in head)
                    blocks.insert(rng.randint(0, len(blocks)), tail)
            out[cid] = assembled

        # Validation skipped for simplicity
        return out

    @staticmethod
    def _gene_to_true_chr(true_root_karyo: Karyotype) -> dict[str, str]:
        m: dict[str, str] = {}
        for cid, genes in true_root_karyo.items():
            for g in genes:
                m[g[0]] = cid
        return m

    @staticmethod
    def _transitions_for_chromosome(genes: list[Gene], gene_to_chr: dict[str, str]) -> int:
        labels = [gene_to_chr.get(g[0]) for g in genes if g[0] in gene_to_chr]
        if len(labels) < 2:
            return 0
        # Count transitions where consecutive labels differ
        return sum(1 for i in range(len(labels) - 1) if labels[i] != labels[i + 1])

    @staticmethod
    def rct_like_score(karyo: Karyotype, gene_to_chr: dict[str, str], min_chr_len: int) -> int:
        score = 0
        for genes in karyo.values():
            if len(genes) < min_chr_len:
                continue
            if _RCTOutgroupSignificance._transitions_for_chromosome(genes, gene_to_chr) == 1:
                score += 1
        return score

    @staticmethod
    def mean_std(values: list[float]) -> tuple[float, float]:
        # Use statistics module for mean and stdev
        if not values:
            return 0.0, 0.0
        try:
            return statistics.mean(values), statistics.stdev(values)
        except statistics.StatisticsError:
            return 0.0, 0.0

    @staticmethod
    def main() -> None:
        cfg_raw = CONFIG.get("rct_significance", {})
        cfg: ConfigDict = cfg_raw if isinstance(cfg_raw, dict) else {}
        true_root_path = str(cfg.get("true_root_path", ""))
        karyotypes_path = str(cfg.get("karyotypes_path", ""))
        out_path = str(cfg.get("out_path", ""))
        outgroup_name = str(cfg.get("outgroup_name", "Outgroup"))
        seed = int(cfg.get("seed", 42))
        replicates = int(cfg.get("replicates", 200))
        min_block_len = int(cfg.get("min_block_len", 50))
        break_rate = float(cfg.get("break_rate", 0.02))
        allow_block_reverse = bool(cfg.get("allow_block_reverse", False))
        min_chr_len = int(cfg.get("min_chr_len", 50))

        true_root = _RCTOutgroupSignificance._read_true_root(true_root_path)
        gene_to_chr = _RCTOutgroupSignificance._gene_to_true_chr(true_root)
        all_data = _RCTOutgroupSignificance._read_karyotypes_species(karyotypes_path)
        species = sorted([s for s in all_data.keys() if s != outgroup_name])

        null_scores: list[int] = []
        for i in range(int(replicates)):
            og = _RCTOutgroupSignificance.build_outgroup_by_block_shuffle(
                true_root_karyo=true_root,
                seed=int(seed) + 1000003 + i,
                min_block_len=int(min_block_len),
                break_rate=float(break_rate),
                allow_block_reverse=bool(allow_block_reverse),
            )
            null_scores.append(_RCTOutgroupSignificance.rct_like_score(og, gene_to_chr, int(min_chr_len)))

        null_mean, null_sd = _RCTOutgroupSignificance.mean_std([float(x) for x in null_scores])

        with open(out_path, "w", encoding="utf-8") as f:
            f.write("\t".join(["species", "score_obs", "null_mean", "null_sd", "pvalue", "z", "replicates", "seed"]) + "\n")
            for sp in species:
                obs = _RCTOutgroupSignificance.rct_like_score(all_data[sp], gene_to_chr, int(min_chr_len))
                ge = sum(1 for x in null_scores if x >= obs)
                p = (1 + ge) / (1 + len(null_scores))
                z = 0.0 if null_sd == 0 else (obs - null_mean) / null_sd
                f.write(f"{sp}\t{obs}\t{null_mean:.6f}\t{null_sd:.6f}\t{p:.6g}\t{z:.6f}\t{len(null_scores)}\t{seed}\n")


def rct_outgroup_significance_main():
    _RCTOutgroupSignificance.main()


def _genes_multiset(karyo: Karyotype) -> dict[str, int]:
    # Use collections.Counter for efficient counting
    return Counter(g[0] for genes in karyo.values() for g in genes)


def _segment_chromosome_into_blocks(genes: list[Gene], rng: random.Random, min_block_len: int, break_rate: float) -> list[list[Gene]]:
    n = len(genes)
    if n <= min_block_len:
        return [list(genes)]
    target_blocks = max(1, int(n * break_rate))
    target_blocks = min(target_blocks, max(1, n // min_block_len))
    blocks: list[list[Gene]] = []
    pos = 0
    remaining_blocks = target_blocks
    while pos < n:
        remaining = n - pos
        if remaining_blocks <= 1:
            blocks.append(list(genes[pos:]))
            break
        min_len = min_block_len
        max_len = remaining - (remaining_blocks - 1) * min_block_len
        if max_len < min_len:
            max_len = min_len
        if max_len > remaining:
            max_len = remaining
        blen = rng.randint(min_len, max_len)
        blocks.append(list(genes[pos : pos + blen]))
        pos += blen
        remaining_blocks -= 1
    return blocks


def build_outgroup_karyotype_by_block_shuffle(true_root_karyo: Karyotype, seed: int, min_block_len: int, break_rate: float, allow_block_reverse: bool = False) -> Karyotype:
    rng = random.Random(seed)
    blocks = []
    length_template = []
    for cid in sorted(true_root_karyo.keys()):
        genes = list(true_root_karyo[cid])
        length_template.append((cid, len(genes)))
        blocks.extend(_segment_chromosome_into_blocks(genes, rng, int(min_block_len), float(break_rate)))

    rng.shuffle(blocks)
    out_karyo: Karyotype = {}
    for cid, target_len in length_template:
        assembled = []
        while len(assembled) < target_len:
            if not blocks:
                raise RuntimeError("Outgroup assembly failed: insufficient blocks to fill chromosomes.")
            b = blocks.pop()
            if allow_block_reverse and rng.random() < 0.5:
                b = list(reversed(b))
            need = target_len - len(assembled)
            if len(b) <= need:
                assembled.extend(b)
            else:
                head = b[:need]
                tail = b[need:]
                assembled.extend(head)
                insert_at = rng.randint(0, len(blocks))
                blocks.insert(insert_at, tail)
        out_karyo[cid] = assembled

    if any(len(genes) != l for (cid, l), genes in zip(length_template, [out_karyo[cid] for cid, _ in length_template])):
        raise RuntimeError("Outgroup assembly failed: chromosome lengths do not match template.")

    a = _genes_multiset(true_root_karyo)
    b = _genes_multiset(out_karyo)
    if a != b:
        missing = [g for g in a if g not in b]
        extra = [g for g in b if g not in a]
        raise RuntimeError(f"Outgroup gene multiset mismatch. missing={len(missing)} extra={len(extra)}")
    return out_karyo


def wrap_tree_with_outgroup(tree: Tree, outgroup_name: str) -> Tree:
    super_root = Tree(name="SuperRoot")
    super_root.add_child(tree)
    super_root.add_child(name=outgroup_name)
    return super_root

def flip_orient(orient: str) -> str:
    return "+" if orient == "-" else "-"

def reverse_segment(genes: list[Gene]) -> list[Gene]:
    return [(g[0], flip_orient(g[1])) for g in reversed(genes)]

def get_context(genes: list[Gene], index: int, window: int = 2) -> str:
    start, end = max(0, index - window), min(len(genes), index + window)
    return f"{' '.join(format_gene(g) for g in genes[start:index])} | {' '.join(format_gene(g) for g in genes[index:end])}"

def build_gene_origin_map(root_karyo: Karyotype) -> dict[str, str]:
    origin: dict[str, str] = {}
    for cid, genes in root_karyo.items():
        for g_name, _ in genes:
            origin[g_name] = cid
    return origin

def extract_fusion_points(karyo: Karyotype) -> set[tuple[str, str]]:
    points: set[tuple[str, str]] = set()
    for _, genes in karyo.items():
        if not genes or len(genes) < 2:
            continue
        # Use zip to avoid range indexing
        for g1, g2 in zip(genes, genes[1:]):
            a, b = (g1[0], g2[0]) if g1[0] <= g2[0] else (g2[0], g1[0])
            points.add((a, b))
    return points

def build_gene_to_chrids(karyo: Karyotype) -> dict[str, set[str]]:
    # Use defaultdict for efficient set construction
    gene_to_chrids: dict[str, set[str]] = defaultdict(set)
    for cid, genes in karyo.items():
        for g_name, _ in genes:
            gene_to_chrids[g_name].add(cid)
    return gene_to_chrids

def validate_lineage_isolation(cfg: ConfigDict, tree: Tree, engine: Any, root_karyo: Karyotype) -> bool:
    if not cfg.get("validate_lineage", True):
        return True

    report_path = os.path.join(cfg["save_dir"], cfg.get("validation_report", "validation_report.txt"))
    name_to_node = {n.name: n for n in tree.traverse("preorder")}

    # Precompute fusion points and gene-to-chromosome mappings
    fp_by_node = {node_name: extract_fusion_points(karyo) for node_name, karyo in engine.node_karyotypes.items()}
    gene_to_chrids_by_node = {node_name: build_gene_to_chrids(karyo) for node_name, karyo in engine.node_karyotypes.items()}

    # Build introduced fusion points per node
    intro_nodes_by_fp = {}
    for node_name, node in name_to_node.items():
        if node.is_root():
            continue
        parent_name = node.up.name
        introduced_raw = fp_by_node.get(node_name, set()) - fp_by_node.get(parent_name, set())
        parent_gene_to_chrids = gene_to_chrids_by_node.get(parent_name, {})
        # Filter out fusion points from same chromosome
        introduced = set()
        for fp in introduced_raw:
            g1, g2 = fp
            c1 = parent_gene_to_chrids.get(g1, set())
            c2 = parent_gene_to_chrids.get(g2, set())
            if c1 and c2 and (c1 & c2):
                continue
            introduced.add(fp)
        for fp in introduced:
            intro_nodes_by_fp.setdefault(fp, []).append(node_name)

    # Collect fusion points in leaves
    leaves = [n.name for n in tree.iter_leaves()]
    valid_fps = set(intro_nodes_by_fp.keys())
    fp_in_leaves = {}
    for leaf in leaves:
        for fp in fp_by_node.get(leaf, set()):
            if fp not in valid_fps:
                continue
            fp_in_leaves.setdefault(fp, []).append(leaf)

    # Find violations
    violations = []
    for fp, intro_nodes in intro_nodes_by_fp.items():
        if len(intro_nodes) <= 1:
            continue
        involved_leaves = fp_in_leaves.get(fp, [])
        lca_name = ""
        lca_has_fp = ""
        if len(involved_leaves) >= 2:
            nodes = [name_to_node[n] for n in involved_leaves if n in name_to_node]
            if nodes:
                lca = tree.get_common_ancestor(nodes)
                lca_name = lca.name
                lca_has_fp = "yes" if fp in fp_by_node.get(lca.name, set()) else "no"
        violations.append((fp, intro_nodes, involved_leaves, lca_name, lca_has_fp))

    # Write report
    root_name = tree.get_tree_root().name
    leaf_union = set()
    for l in leaves:
        leaf_union |= fp_by_node.get(l, set())
    total_fp_leaf_union = len(leaf_union)
    total_fp_interchrom = len(valid_fps)
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("=== LINEAGE ISOLATION VALIDATION (fusion points; parent-chromosome based; inversions excluded) ===\n\n")
        f.write(f"Nodes: {len(name_to_node)}\n")
        f.write(f"Leaves: {len(leaves)}\n")
        f.write(f"Total fusion points across leaves (raw union): {total_fp_leaf_union}\n")
        f.write(f"Total fusion points introduced by inter-chrom events: {total_fp_interchrom}\n")
        f.write(f"Violations (same fusion point introduced on multiple branches): {len(violations)}\n\n")

        if violations:
            f.write("--- Violations ---\n")
            for (a, b), intro_nodes, involved_leaves, lca_name, lca_has_fp in sorted(violations, key=lambda x: (x[0][0], x[0][1])):
                f.write(f"FusionPoint: {a} -- {b}\n")
                f.write(f"IntroducedOn: {', '.join(sorted(intro_nodes))}\n")
                if involved_leaves:
                    f.write(f"LeavesWithPoint: {', '.join(sorted(involved_leaves))}\n")
                if lca_name:
                    f.write(f"LCA: {lca_name} (contains point: {lca_has_fp})\n")
                f.write("\n")


    if violations and cfg.get("validation_strict", True):
        raise RuntimeError(f"Lineage isolation validation failed. See: {report_path}")
    return not violations

def choose_wgd_targets(tree: Tree, cfg: ConfigDict) -> set[str]:
    rng = random.Random(cfg["seed"] + cfg.get("wgd_rng_seed_offset", 0))
    targets: set[str] = set()
    forced = set(cfg.get("force_wgd_nodes", []))
    for node in tree.traverse("preorder"):
        if node.is_root():
            continue
        if node.name in forced:
            targets.add(node.name)
            continue
        if rng.random() < cfg.get("wgd_probability", 0.0):
            targets.add(node.name)
    return targets

def insert_pre_wgd_nodes(tree: Tree, wgd_targets: set[str], suffix: str = "_preWGD") -> Tree:
    existing = set(n.name for n in tree.traverse("preorder"))
    for name in list(wgd_targets):
        n = tree.search_nodes(name=name)
        if not n:
            continue
        node = n[0]
        if node.is_root() or not node.up:
            continue
        parent = node.up
        mid_name = f"{name}{suffix}"
        if mid_name in existing:
            continue
        dist = getattr(node, "dist", 0.0)
        node.detach()
        mid = parent.add_child(name=mid_name)
        try:
            mid.dist = dist
        except Exception:
            pass
        mid.add_child(node)
        existing.add(mid_name)
    return tree

def build_event_plan(tree: Tree, cfg: ConfigDict, wgd_targets: set[str], suffix: str = "_preWGD") -> dict[str, dict[str, Any]]:
    base = cfg["rearrangement_counts"]
    pre_frac = cfg.get("wgd_pre_event_fraction", 0.5)
    plan: dict[str, dict[str, Any]] = {}

    p_ncf = float(base.get("fusion_ncf", 0.0))
    p_eej = float(base.get("fusion_eej", 0.0))
    p_rct = float(base.get("translocation_rct", 0.0))
    p_ncf_pre = 1.0 - (1.0 - p_ncf) ** pre_frac
    p_ncf_post = 1.0 - (1.0 - p_ncf) ** (1.0 - pre_frac)
    p_eej_pre = 1.0 - (1.0 - p_eej) ** pre_frac
    p_eej_post = 1.0 - (1.0 - p_eej) ** (1.0 - pre_frac)
    p_rct_pre = 1.0 - (1.0 - p_rct) ** pre_frac
    p_rct_post = 1.0 - (1.0 - p_rct) ** (1.0 - pre_frac)

    splits: dict[str, tuple[dict[str, Any], dict[str, Any]]] = {}
    for target in sorted(wgd_targets):
        splits[target] = (
            dict(
                translocation_rct=p_rct_pre,
                inversion_prob=float(base.get("inversion_prob", 0.0)),
                fusion_ncf_p=p_ncf_pre,
                fusion_eej_p=p_eej_pre,
            ),
            dict(
                translocation_rct=p_rct_post,
                inversion_prob=float(base.get("inversion_prob", 0.0)),
                fusion_ncf_p=p_ncf_post,
                fusion_eej_p=p_eej_post,
            ),
        )

    for node in tree.traverse("preorder"):
        if node.is_root():
            continue
        if node.name not in plan:
            plan[node.name] = dict(
                inversion_prob=float(base.get("inversion_prob", 0.0)),
                translocation_rct=float(base.get("translocation_rct", 0.0)),
                fusion_ncf_p=float(base.get("fusion_ncf", 0.0)),
                fusion_eej_p=float(base.get("fusion_eej", 0.0)),
            )

    for target, (pre_plan, post_plan) in splits.items():
        plan[target] = post_plan
        plan[f"{target}{suffix}"] = pre_plan
    return plan

# -----------------------------
# 3. TREE BUILDING
# -----------------------------
def rebuild_parent_children(tree: Tree, info: dict[str, dict[str, Any]]) -> None:
    for n in tree.traverse():
        if n.name in info:
            info[n.name].update({"parent": n.up.name if n.up else None, "children": [c.name for c in n.children]})

def build_tree(cfg: ConfigDict) -> tuple[Tree, dict[str, dict[str, Any]]]:
    """
    Build a strict binary tree using ete3's populate method.
    
    For a strict binary tree:
    - Each internal node has exactly 2 children
    - Leaves = Internal nodes + 1
    - Ancestors = Species - 2 (automatically determined)
    """
    random.seed(cfg["seed"])
    num_species = cfg["num_modern_species"]
    
    info: dict[str, dict[str, Any]] = {}
    
    # Use ete3's populate to create a random binary tree
    tree = Tree(name="Root")
    tree.populate(size=num_species, names_library=[f"Sp{i}" for i in range(1, num_species + 1)])
    info["Root"] = {"type": "root"}
    
    # Rename internal nodes to Anc1, Anc2, ... and set info
    ancestor_counter = 0
    for node in tree.traverse():
        if node.is_leaf():
            info[node.name] = {"type": "sp"}
        elif node.name != "Root":
            ancestor_counter += 1
            node.name = f"Anc{ancestor_counter}"
            info[node.name] = {"type": "anc"}

    print(f"Tree built: {num_species} species, {ancestor_counter} ancestors (strict binary)")

    rebuild_parent_children(tree, info)
    return tree, info

def collapse_unary_ancestors(tree: Tree, info: dict[str, dict[str, Any]]) -> tuple[Tree, dict[str, dict[str, Any]]]:
    """
    Collapses unary nodes (nodes with only 1 child).
    If a node P has only 1 child C:
    - If P is Root: C becomes Root.
    - If P is Anc: Remove P, connect C directly to P's parent.
    """
    changed = True
    while changed:
        changed = False
        # Iterate postorder to handle bottom-up
        for p in list(tree.traverse("postorder")):
            # Skip leaves (SpX)
            if p.is_leaf(): continue
            
            # Check if unary
            if len(p.children) == 1:
                child = p.children[0]
                
                # Logic:
                # If P is Root -> C becomes Root (keep C's name if C is Anc/Sp)
                # If P is Anc -> Remove P. Connect C to P.up
                
                if p.is_root():
                    # Root has 1 child. Promote child to root.
                    child.detach()
                    tree = child
                    # If child was Anc or Sp, it keeps its name.
                    # But we need to update info to know it's now root?
                    # Or just keep it as is. Tree root is just the top node.
                    # info[child.name]["type"] = "root" # Optional, or keep as "anc"
                    info.pop(p.name, None) # Remove old root info
                    changed = True
                    # Rebuild info immediately or at end?
                    # Need to update 'tree' variable reference in outer scope? 
                    # We return 'tree', so it's fine.
                else:
                    # P is an internal node with 1 child.
                    # Connect child to P's parent.
                    parent = p.up
                    if parent:
                        child.detach()
                        parent.add_child(child)
                        p.detach()
                        info.pop(p.name, None)
                        changed = True
    
    rebuild_parent_children(tree, info)
    return tree, info

def prune_non_species_leaves(tree: Tree, info: dict[str, dict[str, Any]]) -> tuple[Tree, dict[str, dict[str, Any]]]:
    for leaf in list(tree.iter_leaves()):
        if leaf.name.startswith("Sp"):
            continue
        info.pop(leaf.name, None)
        leaf.detach()
    rebuild_parent_children(tree, info)
    return tree, info

# -----------------------------
# 4. EVOLUTION ENGINE
# -----------------------------
def init_root_karyotype(cfg: ConfigDict) -> Karyotype:
    random.seed(cfg["seed"])
    return {cid: [(f"{cid}{i+1}", "+") for i in range(random.randint(cfg["min_genes_per_chr"], cfg["max_genes_per_chr"]))] 
            for cid in get_chr_labels(cfg["num_ancestor_chromosomes"])}

class EvolutionEngine:
    def __init__(self, config: ConfigDict):
        self.cfg: ConfigDict = config
        self.events: list[dict[str, Any]] = []
        self.node_karyotypes: dict[str, Karyotype] = {}
        self.wgd_nodes: set[str] = set()
        self.used_root_child_fusion_pairs: set[frozenset[str]] = set()
        self.used_root_child_rct_pairs: set[frozenset[str]] = set()

    def apply_evolution(self, tree: Tree, root_karyo: Karyotype, wgd_targets: set[str] | None = None, event_plan: dict[str, dict[str, Any]] | None = None) -> None:
        wgd_targets = set(wgd_targets or [])
        event_plan = event_plan or {}
        suffix = "_preWGD"
        for node in tree.traverse("preorder"):
            if node.is_root():
                self.node_karyotypes[node.name] = copy.deepcopy(root_karyo)
                continue
                
            curr = copy.deepcopy(self.node_karyotypes[node.up.name])

            def apply_plan(karyo: Karyotype, node_obj: Any, name_for_log: str, plan_for_node: dict[str, Any]) -> None:
                is_root_child = bool(node_obj.up and node_obj.up.is_root())
                if random.random() < float(plan_for_node.get("inversion_prob", 0.0)):
                    self._apply_inversion(karyo, name_for_log)
                
                # RCT: enabled on root children, but record the pair for sibling constraint
                if random.random() < float(plan_for_node.get("translocation_rct", 0.0)):
                    used_pair = self._apply_rct(karyo, name_for_log)
                    if used_pair is not None and is_root_child and self.cfg.get("avoid_duplicate_root_fusions", True):
                        self.used_root_child_rct_pairs.add(used_pair)
                
                # NCF/EEJ: forbidden if RCT or fusion already used this pair in sibling branch
                if random.random() < float(plan_for_node.get("fusion_ncf_p", 0.0)):
                    forbidden = None
                    if is_root_child and self.cfg.get("avoid_duplicate_root_fusions", True):
                        forbidden = self.used_root_child_fusion_pairs | self.used_root_child_rct_pairs
                    used_pair = self._apply_ncf(karyo, name_for_log, forbidden_pairs=forbidden)
                    if used_pair is not None and forbidden is not None:
                        self.used_root_child_fusion_pairs.add(used_pair)
                if random.random() < float(plan_for_node.get("fusion_eej_p", 0.0)):
                    forbidden = None
                    if is_root_child and self.cfg.get("avoid_duplicate_root_fusions", True):
                        forbidden = self.used_root_child_fusion_pairs | self.used_root_child_rct_pairs
                    used_pair = self._apply_eej(karyo, name_for_log, forbidden_pairs=forbidden)
                    if used_pair is not None and forbidden is not None:
                        self.used_root_child_fusion_pairs.add(used_pair)

            if node.name.endswith(suffix) and node.children and node.children[0].name in wgd_targets:
                self._log(node.name, "PreWGD", list(curr.keys()), "Before WGD", f"Split to {node.children[0].name}")

            if node.name in wgd_targets and not (node.up and node.up.name == f"{node.name}{suffix}"):
                pre_name = f"{node.name}{suffix}"
                self._log(pre_name, "PreWGD", list(curr.keys()), "Before WGD", f"Split from {node.up.name}->{node.name}")
                pre_plan = event_plan.get(pre_name, {})
                apply_plan(curr, node, pre_name, pre_plan)
                self.node_karyotypes[pre_name] = copy.deepcopy(curr)

            if node.name in wgd_targets:
                curr = self._apply_wgd(curr, node.name)
                self.wgd_nodes.add(node.name)

            apply_plan(curr, node, node.name, event_plan.get(node.name, {}))
            self.node_karyotypes[node.name] = curr

    def _log(self, node: str, etype: str, chrs: list[str], pos: str, ctx: str) -> None:
        self.events.append({"Node": node, "Type": etype, "Chromosomes": chrs, "Position": pos, "Context": ctx})

    def _apply_wgd(self, karyo, node):
        new_k = {}
        for c, g in karyo.items(): new_k[f"{c}_1"], new_k[f"{c}_2"] = copy.deepcopy(g), copy.deepcopy(g)
        self._log(node, "WGD", list(karyo.keys()), "Doubled", f"{len(karyo)}->{len(new_k)}")
        return new_k

    def _apply_inversion(self, karyo, node):
        if not karyo: return
        cid = random.choice(list(karyo.keys()))
        genes = karyo[cid]
        if len(genes) < 2: return
        p1, p2 = sorted(random.sample(range(len(genes)+1), 2))
        if p1 == p2: return
        karyo[cid] = genes[:p1] + [(g[0], flip_orient(g[1])) for g in reversed(genes[p1:p2])] + genes[p2:]
        self._log(node, "Inversion", [cid], f"{p1}-{p2}", "Inverted")

    def _apply_rct(self, karyo, node):
        if len(karyo) < 2: return None
        c1, c2 = random.sample(list(karyo.keys()), 2)
        g1, g2 = karyo[c1], karyo[c2]
        if len(g1) < 4 or len(g2) < 4: return None
        b1, b2 = random.randint(2, len(g1) - 2), random.randint(2, len(g2) - 2)
        aL, aR = g1[:b1], g1[b1:]
        bL, bR = g2[:b2], g2[b2:]

        mode = self.cfg.get("rct_pairing_mode", "random")
        if mode == "random":
            mode = random.choice(["tailswap_only", "alt_only"])

        if mode == "alt_only":
            karyo[c1], karyo[c2] = aL + reverse_segment(bL), reverse_segment(aR) + bR
            pairing = "AL-BL/AR-BR"
            ctx = "Swapped with alternate end pairing"
        else:
            karyo[c1], karyo[c2] = aL + bR, bL + aR
            pairing = "AL-BR/BL-AR"
            ctx = "Swapped tails"

        self._log(node, "RCT", [c1, c2], f"{c1}:{b1}|{c2}:{b2}", f"{ctx}; pairing={pairing}")
        return frozenset((c1, c2))

    def _apply_ncf(self, karyo, node, forbidden_pairs=None):
        if len(karyo) < 2: return
        forbidden_pairs = forbidden_pairs or set()
        keys = list(karyo.keys())
        candidates = [(a, b) for a, b in itertools.combinations(keys, 2) if frozenset((a, b)) not in forbidden_pairs]
        if not candidates:
            self._log(node, "NCF", [], "Skipped", "Skipped: no available unique chromosome pair")
            return None
        a, b = random.choice(candidates)
        base_pair = frozenset((a, b))
        recipient_cid, donor_cid = (a, b) if random.choice([True, False]) else (b, a)
        recipient = karyo[recipient_cid]
        donor = karyo[donor_cid]

        if len(recipient) < 2:
            return None

        insert_pos = random.randint(1, len(recipient) - 1)
        if random.choice([True, False]):
            donor = reverse_segment(donor)
            donor_orient = "TH"
        else:
            donor_orient = "HT"

        karyo[recipient_cid] = recipient[:insert_pos] + donor + recipient[insert_pos:]
        del karyo[donor_cid]
        self._log(
            node,
            "NCF",
            [recipient_cid, donor_cid],
            f"{recipient_cid}:{insert_pos}",
            f"Inserted donor={donor_cid} into recipient={recipient_cid}; donor_orient={donor_orient}",
        )
        return base_pair

    def _apply_eej(self, karyo, node, forbidden_pairs=None):
        if len(karyo) < 2: return
        forbidden_pairs = forbidden_pairs or set()
        keys = list(karyo.keys())
        candidates = [(a, b) for a, b in itertools.combinations(keys, 2) if frozenset((a, b)) not in forbidden_pairs]
        if not candidates:
            self._log(node, "EEJ", [], "Skipped", "Skipped: no available unique chromosome pair")
            return None
        c1, c2 = random.choice(candidates)
        base_pair = frozenset((c1, c2))
        g1, g2 = karyo[c1], karyo[c2]

        mode = self.cfg.get("eej_end_mode", "random")
        if mode == "tail_only":
            end1 = "T"
            end2 = "T" if random.choice([True, False]) else "H"
        elif mode in {"HH", "HT", "TH", "TT"}:
            end1, end2 = mode[0], mode[1]
        else:
            end1 = random.choice(["H", "T"])
            end2 = random.choice(["H", "T"])

        if end1 == "H":
            g1 = reverse_segment(g1)
        if end2 == "T":
            g2 = reverse_segment(g2)

        karyo[c1] = g1 + g2
        del karyo[c2]
        self._log(node, "EEJ", [c1, c2], "End-End", f"Fused; mode={end1}{end2}")
        return base_pair

# -----------------------------
# 5. OUTPUT & VISUALIZATION
# -----------------------------
def save_results(cfg, engine, tree):
    os.makedirs(cfg["save_dir"], exist_ok=True)
    ingroup_tree = tree
    if cfg.get("enable_outgroup_for_rct", False):
        roots = tree.search_nodes(name="Root")
        if roots:
            ingroup_tree = roots[0].copy(method="deepcopy")
        with open(os.path.join(cfg["save_dir"], cfg["output_files"]["tree_nwk_with_outgroup"]), "w") as f:
            f.write(tree.write(format=1))
        with open(os.path.join(cfg["save_dir"], cfg["output_files"]["tree_nwk"]), "w") as f:
            f.write(ingroup_tree.write(format=1))
    else:
        with open(os.path.join(cfg["save_dir"], cfg["output_files"]["tree_nwk"]), "w") as f:
            f.write(tree.write(format=1))
    
    with open(os.path.join(cfg["save_dir"], cfg["output_files"]["events"]), "w", encoding="utf-8") as f:
        f.write("Node\tType\tChromosomes\tPosition\tDetails\n")
        f.write("-" * 80 + "\n")
        for ev in engine.events:
            # Format the event string for better readability
            # Chromosomes involved
            chrs_str = ", ".join(ev['Chromosomes'])
            
            # Write tab-separated or formatted columns
            f.write(f"{ev['Node']:<10}\t{ev['Type']:<15}\t{chrs_str:<15}\t{ev['Position']:<20}\t{ev['Context']}\n")

    leaf_names = set(n.name for n in tree.iter_leaves())

    def write_karyo(fname, filter_nodes=None):
        with open(os.path.join(cfg["save_dir"], fname), "w", encoding="utf-8") as f:
            f.write(f"{'Species':<15} {'Chr':<15} {'Genes'}\n{'-'*80}\n")
            for node in sorted(engine.node_karyotypes):
                # Filter out Anc nodes strictly.
                if node.startswith("Anc"): continue
                
                # If a filter list is provided, respect it (as long as it's not Anc)
                if filter_nodes and node not in filter_nodes: continue
                
                # Determine if we should use numeric IDs (for Species Sp* or any leaf node)
                is_species_like = node.startswith("Sp") or node in leaf_names
                
                # Get chromosomes and sort them
                sorted_cids = sorted(engine.node_karyotypes[node].items())
                
                for idx, (cid, genes) in enumerate(sorted_cids):
                    display_cid = cid
                    if is_species_like:
                        display_cid = str(idx + 1)

                    f.write(f"{node:<10} {display_cid:<15} {' '.join(format_gene_id(g) for g in genes)}\n")

    ingroup_leaves = [n.name for n in ingroup_tree.iter_leaves() if n.name.startswith("Sp")]
    if cfg.get("enable_outgroup_for_rct", False):
        write_karyo(cfg["output_files"]["species_karyotypes"], [n.name for n in tree.iter_leaves()])
        legacy = os.path.join(cfg["save_dir"], "karyotypes_species.txt")
        keep = os.path.join(cfg["save_dir"], cfg["output_files"]["species_karyotypes"])
        if os.path.exists(legacy) and os.path.abspath(legacy) != os.path.abspath(keep):
            try:
                os.remove(legacy)
            except OSError:
                pass
    else:
        write_karyo(cfg["output_files"]["species_karyotypes"], ingroup_leaves)
    
    # Prune tree visualization:
    # Although we collapsed unaries in 'tree', let's ensure the visualization reflects it.
    # The 'tree' object passed here IS the collapsed one from main().
    # However, ete3 rendering might show internal node names if we don't hide them.
    
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.scale = 60
    # Do NOT show internal node names automatically, we will label WGDs specifically.
    ts.show_branch_length = False 
    ts.show_branch_support = False
    
    for n in tree.traverse():
        # Label WGD nodes
        if n.name in engine.wgd_nodes:
            # Add a face for WGD
            n.add_face(TextFace(" ★ WGD", fsize=10, fgcolor="red", bold=True), column=0, position="branch-top")
            n.set_style(NodeStyle(fgcolor="red", size=10))
            
            # Optionally show the node name if it's WGD? 
            # Or just let it be.
        else:
            # Standard node style
            n.set_style(NodeStyle(fgcolor="blue" if not n.is_leaf() else "green", size=5))
            
            # If you want to verify no redundant nodes, maybe print names for debugging?
            # But user wants clean tree.
            # If we collapsed correctly, there are no unaries.
            
    tree.render(os.path.join(cfg["save_dir"], cfg["output_files"]["tree_img"]), tree_style=ts)

def visualize_collinearity(cfg, engine, root_karyo):
    """
    Generate Dot Plots for each species vs Ancestor (Root).
    Style:
    - Top Axis: Ancestor Chromosomes (Color Bars at Top)
    - Left Axis: Species Chromosomes (Numbered 1, 2, 3...)
    - Dots: Colored by Ancestor Chromosome (Larger size)
    - No Left Color Bar
    """
    if not cfg.get("visualize_collinearity", False):
        return

    out_dir = os.path.join(cfg["save_dir"], cfg.get("collinearity_dir", "species_vs_ancestor_plots"))
    os.makedirs(out_dir, exist_ok=True)
    
    # 1. Prepare Ancestor X-axis Layout
    anc_gene_map = {}
    anc_chr_order = sorted(root_karyo.keys()) # ['A', 'B', 'C'...]
    
    current_x = 0
    anc_chr_boundaries = [] # (Name, StartX, EndX, Color)
    
    # Colors
    cmap = plt.get_cmap('tab20')
    anc_colors = {cid: cmap(i % 20) for i, cid in enumerate(anc_chr_order)}
    
    # Map gene to ancestor ID for coloring
    gene_to_anc_id = {}
    
    for cid in anc_chr_order:
        genes = root_karyo[cid]
        start_x = current_x
        for g_tuple in genes:
            g_name = g_tuple[0]
            anc_gene_map[g_name] = current_x
            gene_to_anc_id[g_name] = cid
            current_x += 1
        end_x = current_x
        anc_chr_boundaries.append((cid, start_x, end_x, anc_colors[cid]))
    
    total_genes_x = current_x

    # 2. Generate Plot for each Species
    species_nodes = [n for n in engine.node_karyotypes.keys() if n.startswith("Sp")]
    
    for sp_name in species_nodes:
        sp_karyo = engine.node_karyotypes[sp_name]
        sp_chr_order = sorted(sp_karyo.keys(), key=lambda x: (0, int(x)) if str(x).isdigit() else (1, str(x)))
        
        fig, ax = plt.subplots(figsize=(15, 15)) # Larger figure for clarity
        
        # Prepare Y-axis
        current_y = 0
        sp_yticks = []
        sp_chr_boundaries = [] # (Name, StartY, EndY)
        
        # Collect points
        points_x = []
        points_y = []
        points_c = [] 
        
        for idx, sp_cid in enumerate(sp_chr_order):
            genes = sp_karyo[sp_cid]
            if not genes: continue
            
            start_y = current_y
            
            # Collect points for this chromosome using list comprehension
            for i, g_tuple in enumerate(genes):
                raw_name = g_tuple[0]
                if raw_name in anc_gene_map:
                    points_x.append(anc_gene_map[raw_name])
                    points_y.append(current_y + i)
                    anc_id = gene_to_anc_id.get(raw_name)
                    points_c.append(anc_colors.get(anc_id, 'gray'))
            
            end_y = current_y + len(genes)
            # Use number for label (1, 2, 3...)
            label_num = str(idx + 1)
            sp_chr_boundaries.append((label_num, start_y, end_y))
            sp_yticks.append((start_y + len(genes)/2, label_num))
            
            current_y += len(genes) # Gap removed (was + 5)
            
        total_genes_y = current_y

        # Plot Points
        # Increased dot size (s=35)
        ax.scatter(points_x, points_y, c=points_c, s=35, marker='.', edgecolors='none')
        
        # Draw Ancestor Boundaries (Vertical Grid & Top Bar)
        bar_height = total_genes_y * 0.02
        for (acid, asx, aex, acol) in anc_chr_boundaries:
            # Grid line
            ax.axvline(x=asx, color='black', linestyle='-', linewidth=0.5)
            ax.axvline(x=aex, color='black', linestyle='-', linewidth=0.5)
            
            # Top Color Bar (At Y=0, extending upwards)
            # Since Y is inverted (0 is top), negative Y is "above" 0.
            # We place rectangle at y = -bar_height
            rect = mpatches.Rectangle((asx, -bar_height), aex-asx, bar_height, 
                                      color=acol, clip_on=False)
            ax.add_patch(rect)
            
            # Label (Above the bar)
            # Push text higher up (negative y is up)
            # Adjusted to be closer to the bar (-1.5 instead of -2.5)
            ax.text((asx+aex)/2, -bar_height * 1.5, acid, 
                    ha='center', va='bottom', fontsize=12, fontweight='bold')

        # Draw Species Boundaries (Horizontal Grid only, No Left Bar)
        for (label_num, ssy, sey) in sp_chr_boundaries:
            # Grid line
            ax.axhline(y=ssy, color='black', linestyle='-', linewidth=0.5)
            ax.axhline(y=sey, color='black', linestyle='-', linewidth=0.5)
            
            # Label (Just number, Left side)
            ax.text(-total_genes_x * 0.01, (ssy+sey)/2, label_num, 
                    ha='right', va='center', fontsize=12, fontweight='bold')

        # Setup Axes Limits
        ax.set_xlim(0, total_genes_x)
        ax.set_ylim(0, total_genes_y)
        
        # Turn off standard ticks
        ax.set_xticks([])
        ax.set_yticks([])
        
        # Invert Y to put first chr at top
        ax.invert_yaxis()
        
        # Labels
        ax.set_ylabel(sp_name, fontsize=16, labelpad=40, fontweight='bold')
        # Increase pad to prevent overlap with chromosome letters
        ax.set_title("Root Ancestor", fontsize=16, pad=60, fontweight='bold') 
        
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"DotPlot_{sp_name}.png"), dpi=200)
        plt.close(fig)

def visualize_pairwise_species(cfg, engine):
    """
    Generate Dot Plots for pairwise species comparisons.
    Style:
    - Top Axis: Species A Chromosomes (Numbered)
    - Left Axis: Species B Chromosomes (Numbered)
    - Dots: Black (No Ancestor Color Assumption)
    - All pairwise combinations
    """
    if not cfg.get("visualize_collinearity", False):
        return

    out_dir = os.path.join(cfg["save_dir"], "pairwise_plots")
    os.makedirs(out_dir, exist_ok=True)
    
    species_nodes = sorted([n for n in engine.node_karyotypes.keys() if n.startswith("Sp")])
    
    # Generate ALL pairwise comparisons
    comparisons = list(itertools.combinations(species_nodes, 2))

    print(f"Generating {len(comparisons)} pairwise plots (All Pairs)...")

    for spA_name, spB_name in comparisons:
        spA_karyo = engine.node_karyotypes[spA_name]
        spB_karyo = engine.node_karyotypes[spB_name]

        # 1. Prepare SpA (Top/X-axis) Layout
        spA_gene_map = collections.defaultdict(list)
        spA_chr_order = sorted(spA_karyo.keys(), key=lambda x: (0, int(x)) if str(x).isdigit() else (1, str(x)))
        
        current_x = 0
        spA_chr_boundaries = [] # (Name, StartX, EndX)
        
        for cid in spA_chr_order:
            genes = spA_karyo[cid]
            start_x = current_x
            for g_tuple in genes:
                g_name = g_tuple[0] # "A1"
                spA_gene_map[g_name].append(current_x)
                current_x += 1
            end_x = current_x
            
            label = str(spA_chr_order.index(cid) + 1)
            spA_chr_boundaries.append((label, start_x, end_x))
            
        total_genes_x = current_x
        
        # 2. Prepare SpB (Left/Y-axis) and Plot
        spB_chr_order = sorted(spB_karyo.keys(), key=lambda x: (0, int(x)) if str(x).isdigit() else (1, str(x)))
        
        fig, ax = plt.subplots(figsize=(15, 15))
        
        current_y = 0
        spB_yticks = []
        spB_chr_boundaries = []
        
        points_x = []
        points_y = []
        
        for idx, spB_cid in enumerate(spB_chr_order):
            genes = spB_karyo[spB_cid]
            if not genes: continue
            
            start_y = current_y
            
            # Collect points for this chromosome
            for i, g_tuple in enumerate(genes):
                raw_name = g_tuple[0]
                if raw_name in spA_gene_map:
                    # Plot against ALL copies in SpA
                    for x_pos in spA_gene_map[raw_name]:
                        points_x.append(x_pos)
                        points_y.append(current_y + i)
            
            end_y = current_y + len(genes)
            label_num = str(idx + 1)
            spB_chr_boundaries.append((label_num, start_y, end_y))
            spB_yticks.append((start_y + len(genes)/2, label_num))
            
            current_y += len(genes) # No gap
            
        total_genes_y = current_y
        
        # Plot Dots (Black)
        ax.scatter(points_x, points_y, c='black', s=15, marker='.', edgecolors='none', alpha=0.5)
        
        # Draw SpA Boundaries (Top) - No Color Bars
        # Just Grid and Labels
        for (label, sx, ex) in spA_chr_boundaries:
            ax.axvline(x=sx, color='black', linestyle='-', linewidth=0.5)
            ax.axvline(x=ex, color='black', linestyle='-', linewidth=0.5)
            # Label
            ax.text((sx+ex)/2, -total_genes_y * 0.01, label, 
                    ha='center', va='bottom', fontsize=12, fontweight='bold')
            
        # Draw SpB Boundaries (Left)
        for (label, sy, ey) in spB_chr_boundaries:
            ax.axhline(y=sy, color='black', linestyle='-', linewidth=0.5)
            ax.axhline(y=ey, color='black', linestyle='-', linewidth=0.5)
            ax.text(-total_genes_x * 0.01, (sy+ey)/2, label, 
                    ha='right', va='center', fontsize=12, fontweight='bold')
            
        # Setup Axes
        ax.set_xlim(0, total_genes_x)
        ax.set_ylim(0, total_genes_y)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.invert_yaxis()
        
        # Labels
        ax.set_ylabel(spB_name, fontsize=16, labelpad=40, fontweight='bold')
        ax.set_title(spA_name, fontsize=16, pad=40, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"Pairwise_{spA_name}_vs_{spB_name}.png"), dpi=150)
        plt.close(fig)

# -----------------------------
# 6. MAIN
# -----------------------------
def main():
    print("Evolution Simulator Starting...")
    
    # 1. Build Tree
    print("Building phylogenetic tree...")
    tree, info = build_tree(CONFIG)
    tree, info = collapse_unary_ancestors(tree, info)
    tree, info = prune_non_species_leaves(tree, info)
    tree, info = collapse_unary_ancestors(tree, info)
    wgd_targets = choose_wgd_targets(tree, CONFIG)
    if CONFIG.get("enable_pre_wgd_nodes", True):
        tree = insert_pre_wgd_nodes(tree, wgd_targets)
    event_plan = build_event_plan(tree, CONFIG, wgd_targets)
    
    # 2. Init Root Karyotype
    print("Initializing root karyotype...")
    root_karyo = init_root_karyotype(CONFIG)
    os.makedirs(CONFIG["save_dir"], exist_ok=True)
    true_root_path = os.path.join(CONFIG["save_dir"], CONFIG["output_files"]["true_root_karyotype"])
    with open(true_root_path, "w", encoding="utf-8") as f:
        f.write(f"{'Species':<15} {'Chr':<15} {'Genes'}\n{'-'*80}\n")
        for cid in sorted(root_karyo.keys()):
            genes = root_karyo[cid]
            dedup_genes = dedup_by_gene_id(genes)
            f.write(f"{'TrueRoot':<10} {cid:<15} {' '.join(format_gene_id(g) for g in dedup_genes)}\n")
    
    # 3. Simulate Evolution
    print("Simulating evolution (WGD & Rearrangements)...")
    if CONFIG.get("enable_outgroup_for_rct", False):
        tree = wrap_tree_with_outgroup(tree, CONFIG.get("outgroup_name", "Outgroup"))
    engine = EvolutionEngine(CONFIG)
    engine.apply_evolution(tree, root_karyo, wgd_targets=wgd_targets, event_plan=event_plan)
    if CONFIG.get("enable_outgroup_for_rct", False):
        og_seed = int(CONFIG.get("seed", 0)) + int(CONFIG.get("outgroup_seed_offset", 0))
        outgroup_karyo = build_outgroup_karyotype_by_block_shuffle(
            true_root_karyo=root_karyo,
            seed=og_seed,
            min_block_len=int(CONFIG.get("outgroup_min_block_len", 50)),
            break_rate=float(CONFIG.get("outgroup_break_rate", 0.02)),
            allow_block_reverse=bool(CONFIG.get("outgroup_allow_block_reverse", False)),
        )
        og_name = CONFIG.get("outgroup_name", "Outgroup")
        engine.node_karyotypes[og_name] = outgroup_karyo
    
    # 4. Save Results
    print("Saving results...")
    save_results(CONFIG, engine, tree)
    
    # 5. Visualize Collinearity
    print("Generating species-vs-ancestor dot plots...")
    visualize_collinearity(CONFIG, engine, root_karyo)
    visualize_pairwise_species(CONFIG, engine)

    validate_lineage_isolation(CONFIG, tree, engine, root_karyo)
    
    print(f"Done! Check {CONFIG['save_dir']}/ for outputs.")

if __name__ == "__main__":
    main()
